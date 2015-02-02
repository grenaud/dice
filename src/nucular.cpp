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
	   toaddF.ancCount      = destringify<int>         (fields[0]);
	   toaddF.derCount      = destringify<int>         (fields[1]);
	   toaddF.panelFreqCont = destringify<long double >(fields[2]);
	   toaddF.num           = destringify<int>         (fields[3]);
	   dataToAdd->push_back(toaddF);

       }
       myFile.close();
   }else{
       cerr << "Unable to open file "<<filename<<endl;
       return 1;
   }

   //long double 
   long double e     = 0.001105315;
   long double r     =  0.038809028;
   long double tau_C = 0.531137046;
   long double tau_A =  0.560341912;
   cout<<LogFinalTwoP(dataToAdd,e,r,tau_C,tau_A,true)<<endl;

   return 0;
}

