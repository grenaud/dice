/*
 * nucular
 * Date: Jan-22-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <gzstream.h>
#include <iomanip>      // std::setprecision

#include "utils.h"
#include "libnuc.h"

using namespace std;



int main (int argc, char *argv[]) {

   string line;
   igzstream myFile;
   string filename = string(argv[1]);
   myFile.open(filename.c_str(), ios::in);
   vector<freqSite> * dataToAdd =  new vector<freqSite>();


   bool firstIteration  = true;
   bool twoPopMode  = false;
   bool threePopMode = false;

   if (myFile.good()){
       getline (myFile,line);//header
       while ( getline (myFile,line)){
	   vector<string> fields = allTokens(line,'\t');
	   if(fields.size() != 4 &&
	      fields.size() != 5){
	       cerr<<"Line "<<line<<" does not contain 4 or 5 fields but "<<fields.size()<<" fields"<<endl;
	       exit(1);
	   }
	   
	   if(firstIteration){
	       twoPopMode   =  ( fields.size() == 4 );
	       threePopMode =  ( fields.size() == 5 );
	       firstIteration=false;
	   }

	   freqSite toaddF;

	   if(twoPopMode){
	       toaddF.ancCount      = destringify<int>         (fields[0]);
	       toaddF.derCount      = destringify<int>         (fields[1]);
	       toaddF.panelFreqCont = destringify<long double >(fields[2]);
	       toaddF.num           = destringify<int>         (fields[3]);
	   }else{
	       if(threePopMode){
		   toaddF.ancCount      = destringify<int>         (fields[0]);
		   toaddF.derCount      = destringify<int>         (fields[1]);
		   toaddF.panelFreqCont = destringify<long double >(fields[2]);
		   toaddF.panelFreqAdmx = destringify<long double >(fields[3]);
		   toaddF.num           = destringify<int>         (fields[4]);
	       }else{
		   cerr<<"Line "<<line<<" does not contain 4 or 5 fields"<<endl;
		   exit(1);
	       }
	   }
	   dataToAdd->push_back(toaddF);

       }
       myFile.close();
   }else{
       cerr << "Unable to open file "<<filename<<endl;
       return 1;
   }
   cerr<<"done"<<endl;

   if( twoPopMode  == threePopMode ){
       cerr << "Internal error, cannot have two and three pops at once"<<endl;
       return 1;
   }

   // Set lower boundaries for optimization algorithm
   long double elower         = 0.00001;
   long double rlower         = 0.00001;
   long double tau_Clower     = 0.000001;
   long double tau_Alower     = 0.000001;
   long double admixratelower = 0.000001;
   long double admixtimelower = 0.05;

   // Set upper boundaries for optimization algorithm
   long double eupper         = 0.1;
   long double rupper         = 0.5;
   long double tau_Cupper     = 1.0;
   long double tau_Aupper     = 1.0;
   long double admixrateupper = 0.5;
   long double admixtimeupper = 0.11;


   //Constants
   long double innerdriftY =   0.16;
   long double innerdriftZ =   0.16;
   long double nC          =  20.0 ;
   long double nB          =  20.0 ;

   //Vector Variables
   long double e_i         = randomLongDouble(elower,         eupper);
   long double r_i         = randomLongDouble(rlower,         rupper);
   long double tau_C_i     = randomLongDouble(tau_Clower,     tau_Cupper);
   long double tau_A_i     = randomLongDouble(tau_Alower,     tau_Aupper);
   long double admixrate_i = randomLongDouble(admixratelower, admixrateupper);
   long double admixtime_i = randomLongDouble(admixtimelower, admixtimeupper);


   long double e_i_1;
   long double r_i_1;
   long double tau_C_i_1;
   long double tau_A_i_1;
   long double admixrate_i_1;
   long double admixtime_i_1;

   //to test

   // long double e          = 0.06959408;
   // long double r          = 0.31454792;
   // long double tau_C      = 0.34873276; 
   // long double tau_A      = 0.02521097; 
   // long double admixrate  = 0.16149144; 
   // long double admixtime  = 0.10068321;

   // cout<<LogFinalThreeP(dataToAdd,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,true)<<endl;
   // cout<<"done "<<endl;
   // return 1;


   long double x_il;
   long double x_i_1l;

   if(twoPopMode){
       x_il = LogFinalTwoP(  dataToAdd,e_i,r_i,tau_C_i,tau_A_i,true);
   }else{
       x_il = LogFinalThreeP(dataToAdd,e_i,r_i,tau_C_i,tau_A_i,admixrate_i,admixtime_i,innerdriftY,innerdriftZ,nC,nB,true);
   }

   for(unsigned int chain=0;chain<100000;chain++){
       cout<<"current "<<std::setprecision(10)<<x_il<<"\t"<<e_i<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<"\t"<<admixrate_i<<"\t"<<admixtime_i<<endl;
       long double partition=1000.0;

       //e
       long double facte = fmod((long double)(randomProb()), (eupper-elower)/partition );
       cout<<facte<<endl;
       if(randomBool()){
            e_i_1=e_i+facte;
        }else{
            e_i_1=e_i-facte;
        }

       //r
       long double factr = fmod((long double)(randomProb()), (rupper-rlower)/partition );
       if(randomBool()){
            r_i_1=r_i+factr;
        }else{
            r_i_1=r_i-factr;
        }

       //tau_C
       long double facttau_C = fmod((long double)(randomProb()), (tau_Cupper-tau_Clower)/partition);
       if(randomBool()){
            tau_C_i_1=tau_C_i+facttau_C;
        }else{
            tau_C_i_1=tau_C_i-facttau_C;
        }

       //tau_A
       long double facttau_A = fmod( (long double)(randomProb()), (tau_Aupper-tau_Alower)/partition);
       if(randomBool()){
            tau_A_i_1=tau_A_i+facttau_A;
        }else{
            tau_A_i_1=tau_A_i-facttau_A;
        }

       if(!twoPopMode){

	   //admix_rate  
	   long double factadmixrate = fmod( (long double)(randomProb()), (admixrateupper-admixratelower)/partition);
	   if(randomBool()){
	       admixrate_i_1=admixrate_i+factadmixrate;
	   }else{
	       admixrate_i_1=admixrate_i-factadmixrate;
	   }

	   //admix_time  
	   long double factadmixtime = fmod( (long double)(randomProb()), (admixtimeupper-admixtimelower)/partition);
	   if(randomBool()){
	       admixtime_i_1=admixtime_i+factadmixtime;
	   }else{
	       admixtime_i_1=admixtime_i-factadmixtime;
	   }
       }

       if(twoPopMode){
	   x_i_1l    = LogFinalTwoP(  dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,                                                          true);
       }else{
	   x_i_1l    = LogFinalThreeP(dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,admixrate_i_1,admixtime_i_1,innerdriftY,innerdriftZ,nC,nB,true);
       }

       long double acceptance = min( (long double)(1.0)  , expl(x_i_1l-x_il) );

       cout<< "new   "<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;
       cout<< "ratio "<<std::setprecision(10)<<expl(x_i_1l-x_il)<<"\t"<<(x_i_1l-x_il)<<endl;

       if( (long double)(randomProb()) < acceptance){
	   e_i           =  e_i_1;
	   r_i           =  r_i_1;
	   tau_C_i       =  tau_C_i_1;
	   tau_A_i       =  tau_A_i_1;	  
	   admixrate_i   =  admixrate_i_1;
	   admixtime_i   =  admixtime_i_1;
	   x_il      = x_i_1l;

	   cout<<"new state"<<endl;
       }else{
	   cout<<"reject"<<endl;
       }

       cout<<endl;
       //break;
   }

   return 0;
}

