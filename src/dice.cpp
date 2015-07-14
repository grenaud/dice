/*
 * dice
 * Date: Jan-22-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <gzstream.h>
#include <iomanip>      // std::setprecision
#include <random>

#include "utils.h"
#include "libnuc.h"

using namespace std;


pair<long double,long double> paramsComma(string tosplit){
    vector<string> s=allTokens(tosplit,',');
    if(s.size() != 2){
	cerr<<"No single comma was found in parameter "<<tosplit<<endl;
	exit(1);
    }

    pair<long double,long double> t =  make_pair<long double,long double>(destringify<long double>(s[0]),destringify<long double>(s[1]));

    if(t.second<t.first){
	cerr<<"The second parameter must be greater than the first "<<tosplit<<endl;
	exit(1);
    }
    return t;
}


int main (int argc, char *argv[]) {
    string   outLog  = "/dev/stdout";
    ofstream outLogFP;
    bool twoPopMode   = false;
    bool threePopMode = false;
    double step = 1000;
    int maxChains = 100000;


    //Constants
    long double innerdriftY =    0.16;
    long double innerdriftZ =    0.16;
    long double nC          =  100.0 ;
    long double nB          =  100.0 ;
    


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

    long double e_i         ; //TS
    long double eTS_i       ; //TV

    long double r_i         ;
    long double tau_C_i     ;
    long double tau_A_i     ;
    long double admixrate_i ;
    long double admixtime_i ;

    bool e_i_0         = false; //TS
    bool eTS_i_0       = false;

    bool r_i_0         = false;
    bool tau_C_i_0     = false;
    bool tau_A_i_0     = false;
    bool admixrate_i_0 = false;
    bool admixtime_i_0 = false;


    const string usage=string("\t"+string(argv[0])+
                              " [options]  [input file]"+"\n\n"+

                              "\t\t"+"-2p" +"\t\t\t\t"+"Use 2pop mode (default: none)"+"\n"+
                              "\t\t"+"-3p" +"\t\t\t\t"+"Use 3pop mode (default: none)"+"\n"+

                              "\t\t"+"-o     [output log]" +"\t\t"+"Output log (default: stdout)"+"\n"+
			      
                              "\n\tComputation options:\n"+
                              "\t\t"+"-s     [step]" +"\t\t\t"+"MCMC interval space step (default: "+stringify(step)+")"+"\n"+
                              "\t\t"+"-c     [#chains]" +"\t\t"+"Max. number of Markov chains (default: "+stringify(maxChains)+")"+"\n"+

                              "\n\tStarting values:\n"+
			      "\t\t"+"-e0     [error]"+"\t\t\t"+"Error rate         (default: random)"+"\n"+
			      "\t\t"+"-ets0   [error]"+"\t\t\t"+"Error rate at TSs  (default: random)"+"\n"+
			      "\t\t"+"-r0     [cont]" +"\t\t\t"+"Contamination rate (default: random)"+"\n"+
			      "\t\t"+"-tA0    [tauA]" +"\t\t\t"+"Tau ancient genome        (default: random)"+"\n"+
			      "\t\t"+"-tC0    [tauC]" +"\t\t\t"+"Tau anchor    (default: random)"+"\n"+
			      "\t\t"+"-aR0    [admR]" +"\t\t\t"+"Admixture time     (default: random)"+"\n"+
			      "\t\t"+"-aT0    [admT]" +"\t\t\t"+"Admixture rate     (default: random)"+"\n"+
                              
			      "\n\tRange for parameter values:\n"+
			      "\t\t"+"-e     el,eh"+"\t\t\t"+"Error rate range          (default: "+stringify(elower)         +","+stringify(eupper)+" )"+"\n"+
			      "\t\t"+"-r     rl,rh" +"\t\t\t"+"Contamination rate range  (default: "+stringify(rlower)         +","+stringify(rupper)+" )"+"\n"+
			      "\t\t"+"-tA    tauAl,tauAh" +"\t\t"+"Tau Archaic range         (default: "+stringify(tau_Alower)     +","+stringify(tau_Aupper)+"   )"+"\n"+
			      "\t\t"+"-tC    tauCl,tauCh" +"\t\t"+"Tau Contaminant range     (default: "+stringify(tau_Clower)     +","+stringify(tau_Cupper)+"   )"+"\n"+
			      "\t\t"+"-aR    admRl,admRh" +"\t\t"+"Admixture time range      (default: "+stringify(admixratelower) +","+stringify(admixrateupper)+" )"+"\n"+
			      "\t\t"+"-aT    admTl,admTh" +"\t\t"+"Admixture rate range      (default: "+stringify(admixtimelower) +","+stringify(admixtimeupper)+" )"+"\n"+
			      
			      "\n\tPopulation specific constants:\n"+
			      "\t\t"+"-idy     [drift]" +"\t\t"+"Inner drift Y (default: "+stringify(innerdriftY)+")"+"\n"+
			      "\t\t"+"-idz     [drift]" +"\t\t"+"Inner drift Z (default: "+stringify(innerdriftZ)+")"+"\n"+
			      "\t\t"+"-nc      [num c]"   +"\t\t"+"Number nC (default: "+stringify(nC)+")"+"\n"+
			      "\t\t"+"-nb      [num b]"   +"\t\t"+"Number nB (default: "+stringify(nB)+")"+"\n"+


                              "");


    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<"Usage:"<<endl;
        cout<<""<<endl;
        cout<<usage<<endl;
        return 1;
    }

    string cwdProg=getCWD(argv[0]);

    int lastOpt=1;
    for(int i=1;i<(argc);i++){ //all but the last 3 args
	if(string(argv[i])[0] != '-'  ){
	    lastOpt=i;
	    break;
	}

        if(string(argv[i]) == "-nc"  ){
	    nC  = destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-nb"  ){
	    nB  = destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-idy"  ){
	    innerdriftY  = destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-idz"  ){
	    innerdriftZ  = destringify<long double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-aT0"  ){
	    admixtime_i  = destringify<double>(argv[i+1]);
	    admixtime_i_0=true;
            i++;
            continue;
        }

	if(string(argv[i]) == "-aT"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     admixtimelower = t.first;
	     admixtimeupper = t.second;	     
             i++;
             continue;
         }


        if(string(argv[i]) == "-aR0"  ){
	    admixrate_i  = destringify<double>(argv[i+1]);
	    admixrate_i_0=true;
            i++;
            continue;
        }

         if(string(argv[i]) == "-aR"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     admixratelower = t.first;
	     admixrateupper = t.second;	     
             i++;
             continue;
         }


        if(string(argv[i]) == "-tC0"  ){
	    tau_C_i = destringify<double>(argv[i+1]);
	    tau_C_i_0=true;
            i++;
            continue;
        }

         if(string(argv[i]) == "-tC"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     tau_Clower = t.first;
	     tau_Cupper = t.second;	     
             i++;
             continue;
         }

        if(string(argv[i]) == "-tA0"  ){
	    tau_A_i = destringify<double>(argv[i+1]);
	    tau_A_i_0=true;
            i++;
            continue;
        }

         if(string(argv[i]) == "-tA"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     tau_Alower = t.first;
	     tau_Aupper = t.second;	     
             i++;
             continue;
         }


        if(string(argv[i]) == "-e0"  ){
	    e_i = destringify<double>(argv[i+1]);
	    e_i_0=true;
            i++;
            continue;
        }

        if(string(argv[i]) == "-ets0"  ){
	    eTS_i   = destringify<double>(argv[i+1]);
	    eTS_i_0 = true;
            i++;
            continue;
        }

         if(string(argv[i]) == "-e"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     elower = t.first;
	     eupper = t.second;	     
             i++;
             continue;
         }

        if(string(argv[i]) == "-r0"  ){
	    r_i = destringify<double>(argv[i+1]);
	    r_i_0=true;
            i++;
            continue;
        }


         if(string(argv[i]) == "-r"  ){
	     pair<long double,long double> t = paramsComma(string(argv[i+1]));
	     rlower = t.first;
	     rupper = t.second;	     
             i++;
             continue;
         }



        if(string(argv[i]) == "-2p" ){
	    twoPopMode   = true;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
            threePopMode = true;
            continue;
        }

	if( string(argv[i]) == "-o"  ){
            outLog=string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-s"  ){
            step = destringify<double>(argv[i+1]);
            i++;
            continue;
        }


        if(string(argv[i]) == "-c"  ){
            maxChains = destringify<int>(argv[i+1]);
            i++;
            continue;
        }


        cerr<<"Wrong option "<<string(argv[i])<<endl;
        return 1;
    }


    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
        cerr << "Unable to write to output file "<<outLog<<endl;
        return 1;
    }


    if(twoPopMode == false && threePopMode==false){
	cerr<<"Either specify -2p or -3p "<<endl;
        return 1;
    }

    // cout<<lastOpt<<endl;

   string line;
   igzstream myFile;
   string filename = string(argv[lastOpt+0]);
   myFile.open(filename.c_str(), ios::in);
   vector<freqSite> * dataToAdd =  new vector<freqSite>();


   bool firstIteration  = true;
   bool has4Cols        = false;
   bool has5Cols        = false;
   bool has6Cols        = false;
   bool hasTSInfo           = false;

   if (myFile.good()){
       getline (myFile,line);//header
       vector<string> fieldsHeader = allTokens(line,'\t');
       // cout<<line<<endl;
       if(fieldsHeader[ fieldsHeader.size() - 1 ] == "Ts"){
	   hasTSInfo=true;
       }


       while ( getline (myFile,line)){
	   vector<string> fields = allTokens(line,'\t');

	   // cout<<line<<"\t"<<hasTSInfo<<endl;

	   // if(fields.size() != 4 &&
	   //    fields.size() != 5){
	   //     cerr<<"Line "<<line<<" does not contain 4 or 5 fields but "<<fields.size()<<" fields"<<endl;
	   //     exit(1);
	   // }
	   // unsigned int fieldSize = fields.size();
	   // if(hasTS){	       
	   //}
	   
	   
	   if(firstIteration){
	       has4Cols =  ( (fields.size() - (unsigned int)hasTSInfo)  == 4 );
	       has5Cols =  ( (fields.size() - (unsigned int)hasTSInfo)  == 5 );
	       has6Cols =  ( (fields.size() - (unsigned int)hasTSInfo)  == 6 );

	       if(has4Cols == false && 
		  has5Cols == false &&
		  has6Cols == false ){
		   cerr<<"Line "<<line<<" does not contain 4, 5 or 6 fields but "<<fields.size()<<" fields"<<endl;
		   exit(1);
	       }
	       firstIteration=false;
	   }

	   freqSite toaddF;

	   if(twoPopMode){

	       if(has5Cols){
		   // cerr<<"Line "<<line<<" does not contain 4  but "<<fields.size()<<" fields"<<endl;
		   // exit(1);
		   toaddF.ancCount          = destringify<int>         (fields[0]);
		   toaddF.derCount          = destringify<int>         (fields[1]);
		   toaddF.panelFreqAnchor   = destringify<long double >(fields[2]);
		   toaddF.panelFreqCont     = destringify<long double >(fields[3]);
		   toaddF.num               = destringify<int>         (fields[4]);

		   if(hasTSInfo)
		       toaddF.isTS          = ( fields[5] == "1" );
		       

	       }else{
		   if(has4Cols){
		       toaddF.ancCount          = destringify<int>         (fields[0]);
		       toaddF.derCount          = destringify<int>         (fields[1]);
		       toaddF.panelFreqCont     = destringify<long double >(fields[2]); //anchor is the contaminant
		       toaddF.num               = destringify<int>         (fields[3]);
		       
		       if(hasTSInfo)
			   toaddF.isTS          = ( fields[4] == "1" );

		   }else{
		       cerr<<"Line "<<line<<" does not contain 4 or 5 fields but "<<fields.size()<<" fields"<<endl;
		       exit(1);
		   }
	       }

	   }else{
	       if(threePopMode){
		   if(has5Cols){
		       toaddF.ancCount        = destringify<int>         (fields[0]);
		       toaddF.derCount        = destringify<int>         (fields[1]);
		       toaddF.panelFreqCont   = destringify<long double >(fields[2]); //w and y Contaminant is the admixing population
		       toaddF.panelFreqAnchor = destringify<long double >(fields[3]); //z anchor 
		       toaddF.num             = destringify<int>         (fields[4]); //count

		       if(hasTSInfo)
			   toaddF.isTS          = ( fields[5] == "1" );

		   }else{ //6 cols
		       if(has6Cols){
			   toaddF.ancCount        = destringify<int>         (fields[0]);
			   toaddF.derCount        = destringify<int>         (fields[1]);
			   toaddF.panelFreqAdmx   = destringify<long double >(fields[2]); //y recipient admixing
			   toaddF.panelFreqAnchor = destringify<long double >(fields[3]); //z anchor
			   toaddF.panelFreqCont   = destringify<long double >(fields[4]); //w contamination freq
			   toaddF.num             = destringify<int>         (fields[5]);

			   if(hasTSInfo)
			       toaddF.isTS          = ( fields[6] == "1" );

			   
		       }else{
			   cerr<<"Line "<<line<<" does not contain 5 or 6 fields but "<<fields.size()<<" fields"<<endl;
			   exit(1);
		       }
		   }
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
   //cout<<"done\t"<<has5Cols<<"\t"<<has6Cols<<"\t"<<twoPopMode<<"\t"<<threePopMode<<endl;
   //return 1;
   if( twoPopMode  == threePopMode ){
       cerr << "Internal error, cannot have two and three pops at once"<<endl;
       return 1;
   }


    //vector Variables
   if(!e_i_0)
       e_i       = randomLongDouble(elower,         eupper);

   if(!eTS_i_0)
       eTS_i     = randomLongDouble(elower,         eupper);

   if(!r_i_0)
     r_i         = randomLongDouble(rlower,         rupper);

   if(!tau_C_i_0)
     tau_C_i     = randomLongDouble(tau_Clower,     tau_Cupper);

   if(!tau_A_i_0)
     tau_A_i     = randomLongDouble(tau_Alower,     tau_Aupper);

   if(!admixrate_i_0)
     admixrate_i = randomLongDouble(admixratelower, admixrateupper);

   if(!admixtime_i_0)
     admixtime_i = randomLongDouble(admixtimelower, admixtimeupper);

    long double e_i_1;
    long double eTS_i_1;

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
       x_il = LogFinalTwoP(  dataToAdd,e_i,r_i,tau_C_i,tau_A_i,                                                              has4Cols,hasTSInfo,eTS_i);
   }else{
       x_il = LogFinalThreeP(dataToAdd,e_i,r_i,tau_C_i,tau_A_i,admixrate_i,admixtime_i,innerdriftY,innerdriftZ,nC,nB,cwdProg,has5Cols,hasTSInfo,eTS_i);
   }
   outLogFP<<"chain"<<"\tllik"<<"\terror"<<"\tContRate"<<"\ttau_C"<<"\ttau_A"<<"\tadmixrate"<<"\tadmixtime\tacceptance";
   if(hasTSInfo){
       outLogFP<<"\terrorTV";
   }
   outLogFP<<endl;

   int accept=0;

   random_device rd;
   default_random_engine dre (rd());

   for(int chain=0;chain<maxChains;chain++){
     
       long double partition= (long double)(step);


        // e_i_1         = randomLongDouble(elower,         eupper);
	// r_i_1         = randomLongDouble(rlower,         rupper);
        // tau_C_i_1     = randomLongDouble(tau_Clower,     tau_Cupper);
        // tau_A_i_1     = randomLongDouble(tau_Alower,     tau_Aupper);
        // admixrate_i_1 = randomLongDouble(admixratelower, admixrateupper);
        // admixtime_i_1 = randomLongDouble(admixtimelower, admixtimeupper);

       //e
       normal_distribution<long double> distribution_e(e_i,             (eupper-elower)/partition  );
       e_i_1      = distribution_e(dre);
       // e_i_1      = e_i;

       if(e_i_1 <= elower     ||  e_i_1 >= eupper     ){
	   e_i_1      = e_i;
	   //chain--;
	   //continue;
       }

       if(hasTSInfo){
	   normal_distribution<long double> distribution_eTS(eTS_i,     (eupper-elower)/partition  );
	   eTS_i_1      = distribution_eTS(dre);
	   // e_i_1      = e_i;
	   
	   if(eTS_i_1 <= elower     ||  eTS_i_1 >= eupper     ){
	       eTS_i_1      = eTS_i;
	       //chain--;
	       //continue;
	   }
       }

       normal_distribution<long double> distribution_r(r_i,     (rupper-rlower)/partition  );
       r_i_1      = distribution_r(dre);
       // r_i_1      = r_i;


       if(r_i_1 <= rlower     ||  r_i_1 >= rupper     ){
	   r_i_1      = r_i;
	   //chain--;
	   //continue;
       }
       // cout<<tau_Aupper<<endl;
       // cout<<tau_Cupper<<endl;

       normal_distribution<long double> distribution_tau_C(tau_C_i, (tau_Cupper-tau_Clower)/partition  );
       tau_C_i_1  = distribution_tau_C(dre);

       if(tau_C_i_1 <= tau_Clower ||  tau_C_i_1 >= tau_Cupper ){
	   tau_C_i_1  = tau_C_i;
	   //chain--;
	   //continue;
       }

       normal_distribution<long double> distribution_tau_A(tau_A_i, (tau_Aupper-tau_Alower)/partition  );
       tau_A_i_1  = distribution_tau_A(dre);
     
       if(tau_A_i_1 <= tau_Alower ||  tau_A_i_1 >= tau_Aupper ){
	   tau_A_i_1  = tau_A_i;
       }


       // cout<<"tC\t"<<tau_C_i<<"\t"<<tau_C_i_1<<"\t"<<(tau_C_i_1-tau_C_i)<< endl;
       // cout<<"tA\t"<<tau_A_i<<"\t"<<tau_A_i_1<<"\t"<<(tau_A_i_1-tau_A_i)<<endl<<endl;
      


	  
       if(chain!=0){
	   outLogFP<<chain<<"\t"<<std::setprecision(10)<<x_il<<"\t"<<e_i<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<"\t"<<admixrate_i<<"\t"<<admixtime_i<<"\t"<<double(accept)/double(chain);

	   if(hasTSInfo){
	       outLogFP<<"\t"<<eTS_i;
	   }

	   outLogFP<<endl;
       }
       // cout<<"e"<<e_i<<"\te_1\t"<<e_i_1<<endl;
       // return 1;


       // long double facte = fmod((long double)(randomProb()), (eupper-elower)/partition );
       // //cout<<facte<<endl;
       // if(randomBool()){
       //      r_i_1=e_i+facte;
       //  }else{
       // 	   e_i_1=e_i-facte;
       //  }

       // //r
       // long double factr = fmod((long double)(randomProb()), (rupper-rlower)/partition );
       // if(randomBool()){
       //      r_i_1=r_i+factr;
       //  }else{
       //      r_i_1=r_i-factr;
       //  }

       // //tau_C
       // long double facttau_C = fmod((long double)(randomProb()), (tau_Cupper-tau_Clower)/partition);
       // if(randomBool()){
       //      tau_C_i_1=tau_C_i+facttau_C;
       //  }else{
       //      tau_C_i_1=tau_C_i-facttau_C;
       //  }

       // //tau_A
       // long double facttau_A = fmod( (long double)(randomProb()), (tau_Aupper-tau_Alower)/partition);
       // if(randomBool()){
       //      tau_A_i_1=tau_A_i+facttau_A;
       //  }else{
       //      tau_A_i_1=tau_A_i-facttau_A;
       //  }

       if(!twoPopMode){

	   //admix_rate  

	   normal_distribution<long double> distribution_admixrate(admixrate_i, (admixrateupper-admixratelower)/partition  );
	   admixrate_i_1  = distribution_admixrate(dre);
     
	   if(admixrate_i_1 <= admixratelower ||  admixrate_i_1 >= admixrateupper ){
	       admixrate_i_1  = admixrate_i;
	   }

	   normal_distribution<long double> distribution_admixtime(admixtime_i, (admixtimeupper-admixtimelower)/partition  );
	   admixtime_i_1  = distribution_admixtime(dre);
     
	   if(admixtime_i_1 <= admixtimelower ||  admixtime_i_1 >= admixtimeupper ){
	       admixtime_i_1  = admixtime_i;
	   }
	   
	   // long double factadmixrate = fmod( (long double)(randomProb()), (admixrateupper-admixratelower)/partition);
	   // if(randomBool()){
	   //     admixrate_i_1=admixrate_i+factadmixrate;
	   // }else{
	   //     admixrate_i_1=admixrate_i-factadmixrate;
	   // }

	   //admix_time  
	   // long double factadmixtime = fmod( (long double)(randomProb()), (admixtimeupper-admixtimelower)/partition);
	   // if(randomBool()){
	   //     admixtime_i_1=admixtime_i+factadmixtime;
	   // }else{
	   //     admixtime_i_1=admixtime_i-factadmixtime;
	   // }
       }
       //cout<<"it\t"<<has5Cols<<"\t"<<has6Cols<<"\t"<<twoPopMode<<"\t"<<threePopMode<<endl;
       if(twoPopMode){
	   x_i_1l    = LogFinalTwoP(  dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,                                                                  has4Cols,hasTSInfo,eTS_i_1 );
       }else{
	   x_i_1l    = LogFinalThreeP(dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,admixrate_i_1,admixtime_i_1,innerdriftY,innerdriftZ,nC,nB,cwdProg,has5Cols,hasTSInfo,eTS_i_1 );
       }

       long double acceptance = min( (long double)(1.0)  , expl(x_i_1l-x_il) );

       //cout<< "new   "<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;
       //outLogFP<< "ratio "<<std::setprecision(10)<<expl(x_i_1l-x_il)<<"\tnew "<<(x_i_1l)<<"\told "<<(x_il)<<"\t"<<(x_i_1l-x_il)<<"\t"<<acceptance<<endl;

       //outLogFP<<chain<<"p\t"<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;

       if( (long double)(randomProb()) < acceptance){
	   e_i           =  e_i_1;
	   eTS_i         =  eTS_i_1;

	   r_i           =  r_i_1;
	   tau_C_i       =  tau_C_i_1;
	   tau_A_i       =  tau_A_i_1;	  
	   admixrate_i   =  admixrate_i_1;
	   admixtime_i   =  admixtime_i_1;
	   x_il      = x_i_1l;
	   accept++;
	   //outLogFP<<"new state"<<endl;
       }else{
	   //outLogFP<<"reject"<<endl;
       }

       // cout<<endl;
       //break;
   }
   outLogFP.close();

   return 0;
}

