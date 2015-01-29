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

   if (myFile.good()){
       getline (myFile,line);//header
       while ( getline (myFile,line)){
	   vector<string> fields = allTokens(line,'\t');
	   if(fields.size() != 4){
	       cerr<<"Line "<<line<<" does not contain 4 fields"<<endl;
	       exit(1);
	   }

	   freqSite toaddF;
	   toaddF.ancCount  = destringify<int>         (fields[0]);
	   toaddF.derCount  = destringify<int>         (fields[1]);
	   toaddF.panelFreq = destringify<long double >(fields[2]);
	   toaddF.num       = destringify<int>         (fields[3]);
	   dataToAdd->push_back(toaddF);

       }
       myFile.close();
   }else{
       cerr << "Unable to open file "<<filename<<endl;
       return 1;
   }

   // Set lower boundaries for optimization algorithm
   long double elower  = 0.00001;
   long double rlower  = 0.00001;
   long double tau_Clower  = 0.000001;
   long double tau_Alower  = 0.000001;

   // Set upper boundaries for optimization algorithm
   long double eupper  = 0.1;
   long double rupper  = 0.5;
   long double tau_Cupper  = 1.0;
   long double tau_Aupper  = 1.0;

   long double e_i     = randomLongDouble(elower,eupper);
   long double r_i     = randomLongDouble(rlower,rupper);
   long double tau_C_i = randomLongDouble(tau_Clower,tau_Cupper);
   long double tau_A_i = randomLongDouble(tau_Alower,tau_Aupper);

   // long double e_i     = 0.002; 
   // long double r_i     = 0.038809028;
   // long double tau_C_i = 0.531137046;
   // long double tau_A_i = 0.560341912;
   // long double e_i     = 0.001105315; 
   // long double r_i     = 0.038809028;
   // long double tau_C_i = 0.531137046;
   // long double tau_A_i = 0.560341912;



   long double e_i_1;
   long double r_i_1;
   long double tau_C_i_1;
   long double tau_A_i_1;


   // long double e     = 0.02228953;
   // long double r     = 0.33323844;
   // long double tau_C = 0.41876066;
   // long double tau_A = 0.71627768;
   long double x_il    = LogFinalTwoP(dataToAdd,e_i,r_i,tau_C_i,tau_A_i,true);
   long double x_i_1l;

   for(unsigned int chain=0;chain<100000;chain++){
       cout<<"current "<<std::setprecision(10)<<x_il<<"\t"<<e_i<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<endl;
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

       x_i_1l    = LogFinalTwoP(dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,true);
       long double acceptance = min( (long double)(1.0)  , expl(x_i_1l-x_il) );
       cout<< "new "<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<acceptance<<endl;


       cout<< "ratio "<<std::setprecision(10)<<expl(x_i_1l-x_il)<<"\t"<<(x_i_1l-x_il)<<endl;

       if((long double)(randomProb()) < acceptance){
	   e_i     =  e_i_1;
	   r_i     =  r_i_1;
	   tau_C_i =  tau_C_i_1;
	   tau_A_i =  tau_A_i_1;	  
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

