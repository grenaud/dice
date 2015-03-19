/*
 * endoCaller
 * Date: Mar-26-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#if defined(__CYGWIN__) 
#define atanl(X) atan(X)
#define logl(X) log(X)
#define sqrtl(X) sqrt(X)
#endif
  

// #define DEBUGMST


//GLOBAL CONSTANTS

#define MAXCOV    5000
// #define INDELERRORPROB 1.0e-5 // http://genomebiology.com/2011/12/11/R112
// #define LOGRATIOFORINDEL 50  //beyond that difference in log, a indel will be called
#define MAXMAPPINGQUAL 257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits

// #define IGNOREINDELBOUND 5  //ignore INDEL if there are within this amount of bp of the. 5 is good since it offsets the cost of a gap in a standard SW scoring scheme

// #ifdef	IGNOREINDELBOUND
// #define IGNOREINDELLENGTH 35  //ignore reads of length less than this fpr INDEL calling, this cutoffs was decided because of presence of noise before 35bp
// #endif

#define MIN(a,b) (((a)<(b))?(a):(b))

   
#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>
#include <gzstream.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <math.h>
#include <limits>

#include "utils.h"
#include "miscfunc.h"
#include "ReconsReferenceBAM.h"
#include "MistarParser.h"
#include "mistarOperations.h"

using namespace BamTools;
using namespace std;
const long double PI  = atanl(1.0L)*4;   







char   offsetQual=33;

double likeMatch[MAXMAPPINGQUAL];
double likeMismatch[MAXMAPPINGQUAL];
double likeMatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
double likeMismatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

double likeMatchProb[MAXMAPPINGQUAL];
double likeMismatchProb[MAXMAPPINGQUAL];
double likeMatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
double likeMismatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

double probMapping[MAXMAPPINGQUAL];
double probMismapping[MAXMAPPINGQUAL];



probSubstition illuminaErrorsProb;
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;
vector<probSubstition> sub5pC;
vector<probSubstition> sub3pC;

probSubstition defaultSubMatch;

string dnaAlphabet="ACGT";
map<int, PHREDgeno> pos2phredgeno;

// map<string, long double > read2endoProb; //map seq id to probability that the read is endogenous using a deamination model
// long double read2endoProbInit=false;



template <typename T>
inline string arrayToStringInt(const T toPrint[] ,const int size,const string separator=","){
    if(size == 0){
    	return "";
    }
    string toReturn="";
    for(int i=0;i<(size-1);i++){
    	toReturn+=(stringify(int(toPrint[i]))+separator);
    }
    toReturn+=(stringify(int(toPrint[ size -1 ])));
    return toReturn;
}


inline void transformRef(char * refeBase,char * readBase){
    if( (*refeBase) == 'M'){
	(*refeBase)=(*readBase);
    }
    
}


// inline bool hasIinfirstOrLastTwoBases(const string & reconstructedReference){
//     if(reconstructedReference.length() <= 4){
// 	cerr<<"ERROR read has length less than 4 bp"<<endl;
// 	exit(1);
//     }

//     for(unsigned int j=0;j<2;j++){
// 	if(reconstructedReference[j] == 'I')
// 	    return true;
//     }


//     for(unsigned int j=(reconstructedReference.length()-2);
// 	j<(reconstructedReference.length());
// 	j++){
// 	if(reconstructedReference[j] == 'I')
// 	    return true;
//     }

//     return false;
// }

// inline bool deletionsNextToTwo(const BamAlignment  * al){
//     vector<int> lengthOfNonDels;
//     vector<CigarOp> cigarData=al->CigarData;
//     bool foundDel=false;
//     for(unsigned int i=0;i<cigarData.size();i++){
//         if(cigarData[i].Type == 'D'){
// 	    foundDel=true;
// 	}else{
// 	    lengthOfNonDels.push_back(cigarData[i].Length);
// 	}
//     }

//     if(foundDel){
// 	if(lengthOfNonDels[0]<=2)
// 	    return true;
// 	if(lengthOfNonDels[ lengthOfNonDels.size() -1 ]<=2)
// 	    return true;
//     }
    
//     return false;
// }

 

// //checks for an 'R' or 'S' for soft clip
// inline bool hasBadCharacter(const string & reconstructedReference){

//     for(unsigned int j=0;j<(reconstructedReference.length());j++){
// 	if(reconstructedReference[j] == 'R'  || 
// 	   reconstructedReference[j] == 'S' ){
// 	    return true;
// 	}
//     }
//     return false;
// }

// // if we skip the alignment and cannot get a deamination for this read
// inline bool skipAlign(const string & reconstructedReference,const BamAlignment  * al,unsigned int * skipped){
//     if(hasBadCharacter(reconstructedReference)){
// 	(*skipped)++;
// 	return true;
//     }
	

//     if(hasIinfirstOrLastTwoBases(reconstructedReference)){
// 	(*skipped)++;
// 	return true;
//     }

//     if(deletionsNextToTwo(al)){
// 	(*skipped)++;
// 	return true;
//     }

//     return false;
// }




// //! A method that calls the best nucleotide given the likelihood for the 4 nucleotides
// /*!
//   This method is called by callSingleNucleotide and  will use the information stored in likeBaseNoindel to find the most likely nucleotide and compute the error in this assignment. It can be used for both the contaminant and endogenous.

//   \param bestNuc : The best nucleotide will be stored here
//   \param sumLogLikeAll : Sum of the log-likelihood for every base will be stored here
//   \param sumLogLikeOnlyBest : Log-likelihood for the best nucleotide will be stored here
//   \param sumLogLikeAllButBest : Sum of the log-likelihood for every base excluding the best will be stored here
//   \param sumLogForNucs[] : The log-likelihood of the sum of the remaining bases (ex: For A, only consider C,G,T) for all 4 bases the will be stored here
//   \param likeBaseNoindel: The pre-computed likelihood for all bases
//   \param infoPPos: The vector of structure populated by the bam reader, needed to get the coverage to break ties in likelihood, unlikely to be used
//   \param i:        The position in the genome

// */
// inline void callBestNucleotideGivenLikelihood( int         & bestNuc,
// 					       long double & sumLogLikeAll,        // sum of all the logs
// 					       long double & sumLogLikeOnlyBest,   // sum of the logs for the best
// 					       long double & sumLogLikeAllButBest, // sum of the logs for all but the best
// 					       long double sumLogForNucs[],
// 					       const long double likeBaseNoindel [],
// 					       const vector<singlePosInfo> & infoPPos,
// 					       const int i			       
// 					       ){

//     //init
//     for(int nuc=0;nuc<4;nuc++){	    
// 	sumLogForNucs[nuc]   = 0.0;
//     }

//     sumLogLikeAll          = 0.0; // sum of all the logs
//     sumLogLikeOnlyBest     = 0.0; // sum of the logs for the best
//     sumLogLikeAllButBest   = 0.0; // sum of the logs for all but the best
//     // bool sumLogLikeAllB           = true;
//     // bool sumLogLikeOnlyBestB      = true;
//     // bool sumLogLikeAllButBestB    = true;
//     // bool  sumLogForNucsB[4]; //sumLogForNucs has to be initialized
//     // for(int nuc=0;nuc<4;nuc++){
//     // 	sumLogForNucsB[nuc]=true;
//     // }
    
//     long double bestLike=-INFINITY;
//     //int    bestNuc=-1;

//     //Determining most likely nucleotide
//     for(unsigned int nuc=0;nuc<4;nuc++){
// 	// if(i==146){
// 	//     cout<<nuc<<"\t"<<likeBaseNoindel[nuc]<<endl;
// 	// }
// 	if(likeBaseNoindel[nuc] > bestLike){
// 	    bestLike=likeBaseNoindel[nuc];
// 	    bestNuc=nuc;
// 	}
//     }

//     //If there are more than one with equal likelihood (this is highly unlikely)
//     //take the one with the greatest coverage 
//     vector<int> bestNucs;

//     for(unsigned int nuc=0;nuc<4;nuc++){
// 	if(likeBaseNoindel[nuc] == bestLike){
// 	    bestNucs.push_back(nuc);
// 	}
//     }
	 
//     if(bestNucs.size() > 1){ // multiple equally likely nuc, use coverage to call best one
// 	// cerr<<"size "<<bestNucs.size()<<endl;
// 	// return 1;
// 	int bestCov=-1;
// 	int bestCovN=bestNucs[0];
	     
// 	for(unsigned int bc=0;bc<bestNucs.size();bc++){
// 	    if(infoPPos[i].covPerBase[ bestNucs[bc] ] > bestCov){
// 		bestCov  =infoPPos[i].covPerBase[ bestNucs[bc] ];
// 		bestCovN =                        bestNucs[bc];		     
// 	    }
// 	}
	     
// 	bestNuc = bestCovN;
//     }
//     //end
	 
    
    

//     //computing the probability of error
//     for(int nuc=0;nuc<4;nuc++){
// 	// cout<<(i+1)<<"\tnuc\t"<<dnaAlphabet[nuc]<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<<infoPPos[i].likeBaseNoindel[nuc]<<"\t"<<likeBaseNoindel[nuc]<<"\t"<<pow(10.0,infoPPos[i].likeBaseNoindel[nuc])<<endl;
	     
// 	for(int nuc2=0;nuc2<4;nuc2++){
// 	    if(nuc!=nuc2){
		
// 		//sumLogForNucs[nuc]         += pow(10.0,likeBaseNoindel[nuc2]);
// 		sumLogForNucs[nuc]  = oplusInit(sumLogForNucs[nuc] , likeBaseNoindel[nuc2]);

// 	    }
// 	}


// 	//sumLogLikeAll              +=  pow(10.0,likeBaseNoindel[nuc]);
// 	sumLogLikeAll =  oplusInit(sumLogLikeAll,likeBaseNoindel[nuc]);
	
// 	if(nuc==bestNuc){	    
// 	    //sumLogLikeOnlyBest    +=  pow(10.0,likeBaseNoindel[nuc]);
// 	    sumLogLikeOnlyBest   = oplusInit(sumLogLikeOnlyBest,likeBaseNoindel[nuc]);
// 	}else{

// 	    //sumLogLikeAllButBest  +=  pow(10.0,likeBaseNoindel[nuc]);
// 	    sumLogLikeAllButBest  = oplusInit( sumLogLikeAllButBest,likeBaseNoindel[nuc]);

// 	}		 	    
//     }//end for nuc

    
// }

// //! A method that calls potential insertion in the sample/deletions in the reference
// /*!
//   This method is called by printLogAndGenome(). 
//   When we assume we have a single contaminant :
//      We use insertion2loglikeEndoCont to call both insertions in the contaminant and endogenous. We marginalize over each possible insertion (and no insertion) for the other and determine the most likely insert (or lack thereof) 
//   When we cannot assume we have a single contaminant:
//       Just use insertion2loglike to find the most likely insert (or lack thereof)

//   \param i : Position on the mitonchondria
//   \param genomeRef : The reference genome 
//   \param infoPPos: The vector of structure populated by the bam reader
//   \param singleCont: Boolean as to we assume that we have a single contaminant or not
//   \param minQual: PHRED quality threshold, beyong this we print N instead of the base
//   \param genomeToPrint: String on which the endogenous genome will be printed
//   \param genomeToPrintC: String on which the contaminant genome will be printed
//   \param logToPrint:  Pointer to the string stream for the endogenous log
//   \param logToPrintC: Pointer to the string stream for the contaminant log
//   \param setFlags:  Boolean to say whether we skip printing to the genome/log or not
//   \param endoIndel:   Boolean to know if the endogenous has a deletion
//   \param contIndel:   Boolean to know if the contaminanthas a deletion

// */
// void insertionInSample(const int i,
// 		       const string & genomeRef,
// 		       const vector<singlePosInfo> & infoPPos, 
// 		       const bool singleCont,
// 		       const int minQual,			  
// 		       string & genomeToPrint,
// 		       string & genomeToPrintC,
// 		       stringstream * logToPrint,
// 		       stringstream * logToPrintC,
// 		       const bool setFlags,
// 		       bool & endoIndel,
// 		       bool & contIndel){


//     if(infoPPos[i].allInserts.size() != 0){  // there are potential insertions


// 	if(singleCont){
// 	    //calling the endogenous
// 	    string       bestInsertEndo        = "";
// 	    long double  bestInsertLogLikeEndo = -1.0*numeric_limits<long double>::infinity();
// 	    bool         initializeBEndo       = false;
// 	    long double  sumLogLikeEndo        = 0.0;

// 	    //iterate for each endogenous insert
// 	    for(set<string>::const_iterator it1 = infoPPos[i].allInserts.begin(); 
//  		it1 != infoPPos[i].allInserts.end(); 
//  		++it1) {

// 		long double sumLogLikeForThisIns=0.0;

// 		//marginalize over each possible contaminant
// 		for(set<string>::const_iterator it2 = infoPPos[i].allInserts.begin(); 
// 		    it2 != infoPPos[i].allInserts.end(); 
// 		    ++it2) {
// 		    pair<string,string> keytouse (*it1,*it2);
// 		    // infoPPos[i].insertion2loglikeEndoCont[ keytouse ] = 0.0;

// #ifdef DEBUGINS
// 		    if(i==DEBUGINS){
// 			//cout<<"i\t"<<i<<"\t"<<"inse="<<keytouse.first<<"#"<<"\tinsc="<<keytouse.second<<"#"<<"\t" << infoPPos[i].insertion2loglikeEndoCont.at( keytouse ) <<endl;
// 			cout<<"i\t"<<i<<"\t"<<"inse=#"<<keytouse.first<<"#"<<infoPPos[i].insertion2count.at(keytouse.first)<<"\tinsc=#"<<keytouse.second<<"#"<<infoPPos[i].insertion2count.at(keytouse.second)<<"\t" << infoPPos[i].insertion2loglikeEndoCont.at( keytouse ) <<endl;
// 		    }
// #endif	
	
// 		    sumLogLikeForThisIns =  oplusInit(sumLogLikeForThisIns,
// 						      infoPPos[i].insertion2loglikeEndoCont.at( keytouse ));
// 		}
		 
// 		if(!initializeBEndo){
// 		    bestInsertEndo            = *it1;
// 		    bestInsertLogLikeEndo     = sumLogLikeForThisIns;
// 		    initializeBEndo           = true;
// 		}else{
// 		    if(bestInsertLogLikeEndo < sumLogLikeForThisIns){
// 			bestInsertEndo        = *it1;
// 			bestInsertLogLikeEndo = sumLogLikeForThisIns;
// 		    }
// 		}

// 		//cout<<i<<"\tins=#"<<*it1<<"#\t"<<sumLogLikeForThisIns<<"\t"<<setFlags<<"\t"<<sumLogLikeEndo<<"\tbest=\t"<<bestInsertEndo<<"#\t"<<bestInsertLogLikeEndo<<endl;

// 		sumLogLikeEndo = oplusInit(sumLogLikeEndo,sumLogLikeForThisIns);
// 		//cout<<i<<"\tins=#"<<*it1<<"#\t"<<sumLogLikeForThisIns<<"\t"<<setFlags<<"\t"<<sumLogLikeEndo<<"\tbest=\t"<<bestInsertEndo<<"#\t"<<bestInsertLogLikeEndo<<endl;

// 	    }//end for each endogenous insert

     
// 	     string       bestInsertCont        = "";
// 	     long double  bestInsertLogLikeCont = -1.0*numeric_limits<long double>::infinity();;
// 	     bool         initializeBCont       = false;
// 	     long double  sumLogLikeCont        = 0.0;

// 	     //iterate for each contaminant insert
// 	     for(set<string>::const_iterator it2 = infoPPos[i].allInserts.begin(); 
// 		 it2 != infoPPos[i].allInserts.end(); 
// 		 ++it2) {

// 		 long double sumLogLikeForThisIns=0.0;

// 		 //marginalize over each possible endegenous		
// 		 for(set<string>::const_iterator it1 = infoPPos[i].allInserts.begin(); 
// 		     it1 != infoPPos[i].allInserts.end(); 
// 		     ++it1) {
// 		     pair<string,string> keytouse (*it1,*it2);
// 		     // infoPPos[i].insertion2loglikeEndoCont[ keytouse ] = 0.0;
// 		     sumLogLikeForThisIns =  oplusInit(sumLogLikeForThisIns,
// 						       infoPPos[i].insertion2loglikeEndoCont.at( keytouse ));
// 		 }
		 
// 		 if(!initializeBCont){
// 		     bestInsertCont            = *it2;
// 		     bestInsertLogLikeCont     = sumLogLikeForThisIns;
// 		     initializeBCont           = true;
// 		 }else{
// 		     if(bestInsertLogLikeCont < sumLogLikeForThisIns){
// 			 bestInsertCont        = *it2;
// 			 bestInsertLogLikeCont = sumLogLikeForThisIns;
// 		     }
// 		 }
// 		 sumLogLikeCont = oplusInit(sumLogLikeCont,sumLogLikeForThisIns);
// 	     }//end for each endogenous insert
	     
// 	     endoIndel=(bestInsertEndo != "");
// 	     contIndel=(bestInsertCont != "");

// 	     if(setFlags){
// 		 return ;
// 	     }


// 	     if(bestInsertEndo != ""){ //more likely there is an insert in the endogenous
// 		long double  sumLogLikeButTheBestEndo=log( pow(10.0,sumLogLikeEndo) - pow(10.0,bestInsertLogLikeEndo) )/log(10);
// 		//cout<<i<<"\t"<<setFlags<<"\t"<<sumLogLikeEndo<<"\tbest=\t"<<bestInsertEndo<<"#\t"<<bestInsertLogLikeEndo<<"\t"<<sumLogLikeButTheBestEndo<<"\t"<<(sumLogLikeEndo==bestInsertLogLikeEndo)<<"\t"<<(sumLogLikeEndo==bestInsertLogLikeEndo)<<endl;
// 		long double qualInsToPrint ;

// 		if(sumLogLikeEndo==bestInsertLogLikeEndo)
// 		    qualInsToPrint = 1.0/INDELERRORPROB; //There is no likely alternative
// 		else
// 		    qualInsToPrint = -10.0*( sumLogLikeButTheBestEndo - sumLogLikeEndo);
		
// 		for(unsigned int k=0;k<(bestInsertEndo.size());k++){
// 		    (*logToPrint)<<(i+1)<<"i\t"<<"-"<<"\t"<<bestInsertEndo[k]<<"\t"<<qualInsToPrint<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].insertion2count.at(bestInsertEndo)<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		}

// 		//cout<<"inse "<<i<<"\t"<<qualInsToPrint<<"\t"<<minQual<<"\t"<<bestInsertEndo<<endl;
// 		if(qualInsToPrint >= minQual){
// 		    //cout<<"inse "<<i<<"\t"<<qualInsToPrint<<"\t"<<minQual<<"\t"<<bestInsertEndo<<endl<<genomeToPrint<<endl<<genomeToPrintC<<endl;
// 		    genomeToPrint+=bestInsertEndo;
// 		    //cout<<"inse "<<i<<"\t"<<qualInsToPrint<<"\t"<<minQual<<"\t"<<bestInsertEndo<<endl<<genomeToPrint<<endl<<genomeToPrintC<<endl;
// 		}
// 	    }


// 	     if(bestInsertCont != ""){ //more likely there is an insert in the endogenous
// 		long double  sumLogLikeButTheBestCont=log( pow(10.0,sumLogLikeCont) - pow(10.0,bestInsertLogLikeCont) )/log(10);
// 		long double qualInsToPrint ;

// 		if(sumLogLikeCont==bestInsertLogLikeCont)
// 		    qualInsToPrint = 1.0/INDELERRORPROB; //There is no likely alternative
// 		else
// 		    qualInsToPrint = -10.0*( sumLogLikeButTheBestCont - sumLogLikeCont);


// 		for(unsigned int k=0;k<(bestInsertCont.size());k++){
// 		    (*logToPrintC)<<(i+1)<<"i\t"<<"-"<<"\t"<<bestInsertCont[k]<<"\t"<<qualInsToPrint<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].insertion2count.at(bestInsertCont)<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		}
// 		//cout<<"insc "<<i<<"\t"<<qualInsToPrint<<"\t"<<minQual<<"\t"<<bestInsertEndo<<endl;		
// 		if(qualInsToPrint >= minQual)
// 		    genomeToPrintC+=bestInsertCont;
// 	    }





	    
// 	}else{ //cannot assume we have a single contaminant
// 	    //for each potential insert
// 	    string       bestInsert        = "";
// 	    long double  bestInsertLogLike = -1.0*numeric_limits<long double>::infinity();
// 	    bool         initializeB       = false;
// 	    long double  sumLogLike        = 0.0;

// 	    for(set<string>::const_iterator it1 = infoPPos[i].allInserts.begin(); 
// 		it1 != infoPPos[i].allInserts.end(); 
// 		++it1) {
// 		if(!initializeB){
// 		    bestInsertLogLike     = infoPPos[i].insertion2loglike.at(*it1);
// 		    bestInsert            = *it1;		    
// 		    initializeB = true;
// 		}else{
// 		    if( bestInsertLogLike < infoPPos[i].insertion2loglike.at(*it1) ){
// 			bestInsertLogLike = infoPPos[i].insertion2loglike.at(*it1);
// 			bestInsert        = *it1;		    
// 		    }
// 		}
// 		sumLogLike = oplusInit(sumLogLike,infoPPos[i].insertion2loglike.at(*it1));
// 	    }
	    
// 	    endoIndel=(bestInsert != "");
// 	    if(setFlags){
// 		return ;
// 	    }
	    
// 	    if(bestInsert != ""){ //more likely there is an insert
// 		long double  sumLogLikeButTheBest=log( pow(10.0,sumLogLike) - pow(10.0,bestInsertLogLike) )/log(10);
// 		long double qualInsToPrint ;

// 		if(sumLogLike==bestInsertLogLike)
// 		    qualInsToPrint = 1.0/INDELERRORPROB; //There is no likely alternative
// 		else
// 		    qualInsToPrint = -10.0*( sumLogLikeButTheBest - sumLogLike);


// 		for(unsigned int k=0;k<(bestInsert.size());k++){
// 		    (*logToPrint)<<(i+1)<<"i\t"<<"-"<<"\t"<<bestInsert[k]<<"\t"<<qualInsToPrint<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].insertion2count.at(bestInsert)<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		}
// 		//adding in fasta file if qual is higher than threshold
// 		if(qualInsToPrint >= minQual)
// 		    genomeToPrint+=bestInsert;
// 	    }
	    
// 	}

	
	

	
//     }




    

// }

// //! A method that calls potential deletion in the sample/insertion in the reference
// /*!
//   This method is called by printLogAndGenome(). 

//   When we assume we have a single contaminant :
//      We use llikDeletionBoth, llikDeletionEndo, llikDeletionEndo, llikDeletionCont, llikDeletionNone 
//      from infoPPos to find the most likely state for both the endogenous and the contaminant.

//   When we cannot assume we have a single contaminant:
//      Use llikDeletion llikNoDeletion to find out wether a deletion is more likely than the

//   \param i : Position on the mitonchondria
//   \param genomeRef : The reference genome 
//   \param infoPPos: The vector of structure populated by the bam reader
//   \param singleCont: Boolean as to we assume that we have a single contaminant or not
//   \param minQual: PHRED quality threshold, beyong this we print N instead of the base
//   \param logToPrint:  Pointer to the string stream for the endogenous log
//   \param logToPrintC: Pointer to the string stream for the contaminant log
//   \param skipEndo:  Boolean set by the method for the endogenous sample, set to 1 if the sample has a deletion hence no need to call a base
//   \param skipCont:  Boolean set by the method for the contaminant, set to 1 if the contaminant has a deletion hence no need to call a base
//   \param setFlags:  Boolean to say whether we skip printing to the genome/log or not
//   \param endoIndel:   Boolean to know if the endogenous has a deletion
//   \param contIndel:   Boolean to know if the contaminanthas a deletion
// */
// void deletionInSample(const int i,
// 		      const string & genomeRef,
// 		      const vector<singlePosInfo> & infoPPos,
// 		      const bool singleCont,
// 		      const int minQual,
// 		      string & genomeToPrint,
// 		      string & genomeToPrintC,			  
// 		      stringstream * logToPrint,
// 		      stringstream * logToPrintC,
// 		      bool & skipEndo,
// 		      bool & skipCont,
// 		      const bool setFlags,
// 		      bool & endoIndel,
// 		      bool & contIndel){
    
//     if(  infoPPos[i].numDel >= 0){

	    

// 	if(singleCont){
// 	    long double logLikeDel[] = { infoPPos[i].llikDeletionBoth,infoPPos[i].llikDeletionEndo,infoPPos[i].llikDeletionCont,infoPPos[i].llikDeletionNone };
	    
// 	    pair<long double,long double> maxAndSecondMax = firstAndSecondHighestArray( logLikeDel,4 );

// 	    long double sumLogLikeAll=0.0;
// 	    for(int n=0;n<4;n++){
// 		sumLogLikeAll =  oplusInit(sumLogLikeAll,logLikeDel[n]);
// 	    }

// 	    // cout<<maxAndSecondMax.first<<"\t"<<maxAndSecondMax.second<<endl;

// 	    // if(i==513){
// 	    // 	cout<<"POSllik\t"<<infoPPos[i].llikDeletionBoth<<"\t"<<infoPPos[i].llikDeletionEndo<<"\t"<<infoPPos[i].llikDeletionCont<<"\t"<<infoPPos[i].llikDeletionNone<<endl;
// 	    // }

// 	    //both the contaminant and endegenous have a deletion wrt the reference, no need to call any base
// 	    if( (infoPPos[i].llikDeletionBoth == maxAndSecondMax.first) &&
// 		((infoPPos[i].llikDeletionBoth - maxAndSecondMax.second) > LOGRATIOFORINDEL)){
	

// 		endoIndel=true;
// 		contIndel=true;
		    
// 		if(setFlags){
// 		    return ;
// 		}

// 		// cout<<"both"<<endl;
// 		long double sumLogLikeCase=0.0;
// 		for(int n=0;n<4;n++){
// 		    if(n!=0)
// 			sumLogLikeCase =  oplusInit(sumLogLikeCase,logLikeDel[n]);
// 		}

// 		long double qualDelToPrint = (-10.0*(sumLogLikeCase-sumLogLikeAll  )/log(10.0));

// 		(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<(qualDelToPrint)<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

// 		if(qualDelToPrint >= minQual){ //deletion is high quality, do nothing

// 		}else{ //deletion is low quality, put an N
// 		    genomeToPrint+="N";
// 		}

// 		(*logToPrintC)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<(qualDelToPrint)<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

// 		if(qualDelToPrint >= minQual){ //deletion is high quality, do nothing

// 		}else{ //deletion is low quality, put an N
// 		    genomeToPrintC+="N";
// 		}

// 		PHREDgeno toadd;
		
// 		toadd.ref       = genomeRef[i];
// 		toadd.consensus = 'D';
		
// 		pos2phredgeno[   (i+1)   ] = toadd;
		
// 		//continue;
// 		skipEndo=true;
// 		skipCont=true;
// 	    }
	    

// 	    //deletion only in the endogenous, need to call the base for the contaminant
// 	    if( (infoPPos[i].llikDeletionEndo == maxAndSecondMax.first) &&
// 		((infoPPos[i].llikDeletionEndo - maxAndSecondMax.second) > LOGRATIOFORINDEL)){

// 		endoIndel=true;
// 		contIndel=false;

// 		if(setFlags){
// 		    //		    cout<<"i="<<i<<"\t"<<endoIndel<<"\t"<<contIndel<<"\tendoDel"<<endl;
// 		    return ;
// 		}
// 		// cout<<"endo"<<endl;

// 		long double sumLogLikeCase=0.0;
// 		for(int n=0;n<4;n++){
// 		    if(n!=1)
// 			sumLogLikeCase =  oplusInit(sumLogLikeCase,logLikeDel[n]);
// 		}

// 		long double qualDelToPrint =(-10.0*(sumLogLikeCase-sumLogLikeAll  )/log(10.0));
// 		(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<< qualDelToPrint<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		//(*logToPrintD)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<(-10.0*(sumLogLikeCase-sumLogLikeAll )/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		if(qualDelToPrint >= minQual){ //deletion is high quality, do nothing

// 		}else{ //deletion is low quality, put an N
// 		    genomeToPrint+="N";
// 		}

// 		PHREDgeno toadd;
		
// 		toadd.ref       = genomeRef[i];
// 		toadd.consensus = 'D';
		
// 		pos2phredgeno[   (i+1)   ] = toadd;
		
// 		//continue;
// 		skipEndo=true;
// 		skipCont=false;
// 	    }



// 	    //deletion only in the contaminant, need to call the base for the endogenous
// 	    if( (infoPPos[i].llikDeletionCont == maxAndSecondMax.first) &&
// 		((infoPPos[i].llikDeletionCont - maxAndSecondMax.second) > LOGRATIOFORINDEL)){
// 		// cout<<"cont"<<endl;
// 		endoIndel=false;
// 		contIndel=true;

// 		if(setFlags){
// 		    return ;
// 		}

// 		long double sumLogLikeCase=0.0;
// 		for(int n=0;n<4;n++){
// 		    if(n!=2)
// 			sumLogLikeCase =  oplusInit(sumLogLikeCase,logLikeDel[n]);
// 		}

// 		//(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<(-10.0*(sumLogLikeCase-sumLogLikeAll  )/log(10.0))<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
// 		long double qualDelToPrint = (-10.0*(sumLogLikeCase-sumLogLikeAll )/log(10.0)) ;

// 		(*logToPrintC)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<qualDelToPrint <<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;

// 		if(qualDelToPrint >= minQual){ //deletion is high quality, do nothing

// 		}else{ //deletion is low quality, put an N
// 		    genomeToPrintC+="N";
// 		}

// 		// PHREDgeno toadd;
		
// 		// toadd.ref       = genomeRef[i];
// 		// toadd.consensus = 'D';
		
// 		// pos2phredgeno[   (i+1)   ] = toadd;

// 		//continue;
// 		skipEndo=false;
// 		skipCont=true;
// 	    }





// 	}else{
	
// 	    // double llikDel  =0;
// 	    // double llikNoDel=0;
		
// 	    if( infoPPos[i].llikDeletion - infoPPos[i].llikNoDeletion > LOGRATIOFORINDEL){
// 		endoIndel=true;
// 		if(setFlags){
// 		    return ;
// 		}
// 		long double qualDelToPrint = -10.0*( (infoPPos[i].llikNoDeletion - oplus(infoPPos[i].llikDeletion,infoPPos[i].llikNoDeletion))/log(10.0));

// 		(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<"D"<<"\t"<<qualDelToPrint<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<infoPPos[i].numDel<<"\t0.0\t0.0\t0.0\t0.0"<<endl;
		
// 		if(qualDelToPrint >= minQual){ //deletion is high quality, do nothing

// 		}else{ //deletion is low quality, put an N
// 		    genomeToPrint+="N";
// 		}

// 		PHREDgeno toadd;
		
// 		toadd.ref       = genomeRef[i];
// 		toadd.consensus = 'D';
		
// 		pos2phredgeno[   (i+1)   ] = toadd;
		
// 		//continue;
// 		skipEndo=true;
// 		skipCont=true;

// 	    }
// 	}
// 	//TODO add contamination detection
	     
//     }

// }






// //! A method that prints the log for bases without any coverage
// /*!
//   This method is called by printLogAndGenome() and just prints a simple line saying there was no coverage.

//   \param i : Position on the mitonchondria
//   \param genomeRef : The reference genome 
//   \param infoPPos: The vector of structure populated by the bam reader
//   \param outLogCflag: Flag to say we print to the contaminant log as well.
//   \param genomeToPrint: String on which the endogenous genome will be printed
//   \param genomeToPrintC: String on which the contaminant genome will be printed
//   \param logToPrint:  Pointer to the string stream for the endogenous log
//   \param logToPrintC: Pointer to the string stream for the contaminant log
//   \param skipEndo:  Boolean if we skip the  endogenous sample, set to 1 if the endogenous sample has no coverage
//   \param skipCont:  Boolean set by the method for the contaminant, set to 1 if the contaminant has no coverage
// */
// void noCoverage(const int i,
// 		const string & genomeRef,
// 		const vector<singlePosInfo> & infoPPos,
// 		const bool outLogCflag,
// 		string & genomeToPrint,
// 		string & genomeToPrintC,
// 		stringstream * logToPrint,
// 		stringstream * logToPrintC,
// 		bool & skipEndo,
// 		bool & skipCont){


//     if(infoPPos[i].cov == 0){
// 	genomeToPrint+="N";
// 	(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\tN\t0\t0\t0\t0\t0.0\t0.0\t0.0\t0.0"<<endl;
// 	if(outLogCflag){
// 	    genomeToPrintC+="N";
// 	    (*logToPrintC)<<(i+1)<<"\t"<<genomeRef[i]<<"\tN\t0\t0\t0\t0\t0.0\t0.0\t0.0\t0.0"<<endl;	    
// 	}
	
// 	skipEndo=true;
// 	skipCont=true;	
//     }

// }


// //! A method that computes the most likely single nucleotide
// /*!
//   This method is called by printLogAndGenome() and either computes:
//   When we assume we have a single contaminant :
//      use likeBaseNoindelCont to compute the most likely endogenous base and contaminant
//      we marginalize over each contaminant base to call the endogenous base and vice-versa
//   When we cannot assume we have a single contaminant:
//      use likeBaseNoindel for all four endogenous nucleotides
  
//   It calls callBestNucleotideGivenLikelihood() for each possibility 
//   \param i : Position on the mitonchondria
//   \param genomeRef : The reference genome 
//   \param infoPPos: The vector of structure populated by the bam reader
//   \param singleCont: Boolean as to we assume that we have a single contaminant or not
//   \param minQual: PHRED quality threshold, beyong this we print N instead of the base
//   \param genomeToPrint: String on which the endogenous genome will be printed
//   \param genomeToPrintC: String on which the contaminant genome will be printed
//   \param logToPrint:  Pointer to the string stream for the endogenous log
//   \param logToPrintC: Pointer to the string stream for the contaminant log
//   \param skipEndo:  If there was a deletion, we do not print to the endogenous sample
//   \param skipCont:  If there was a deletion, we do not print to the contaminant
// */
// void callSingleNucleotide(const int i,
// 			  const string & genomeRef,
// 			  const vector<singlePosInfo> & infoPPos,
// 			  const bool singleCont,
// 			  const int minQual,			  
// 			  string & genomeToPrint,
// 			  string & genomeToPrintC,
// 			  stringstream * logToPrint,
// 			  stringstream * logToPrintC,
// 			  bool & skipEndo,
// 			  bool & skipCont,
// 			  const bool endoIndelCurrent,
// 			  const bool contIndelCurrent){
//     int    bestNuc=-1;
//     int    bestNucC=-1;

//     long double sumLogLikeAll           = 0.0; // sum of all the logs
//     long double sumLogLikeOnlyBest      = 0.0; // sum of the logs for the best
//     long double sumLogLikeAllButBest    = 0.0; // sum of the logs for all but the best
//     long double sumLogLikeAllC          = 0.0; // sum of all the logs
//     long double sumLogLikeOnlyBestC     = 0.0; // sum of the logs for the best
//     long double sumLogLikeAllButBestC   = 0.0; // sum of the logs for all but the best

//     // bool sumLogLikeAllB           = true;
//     // bool sumLogLikeOnlyBestB      = true;
//     // bool sumLogLikeAllButBestB    = true;
	 
//     // int nuc=0;
//     // sumLogLikeAll              =  infoPPos[i].likeBaseNoindel[nuc];  //oplus= log10( pow(10,x)+pow(10,y) )
//     // if(nuc==bestNuc)
//     //     sumLogLikeOnlyBest     =  infoPPos[i].likeBaseNoindel[nuc];
//     // else
//     //     sumLogLikeAllButBest   =  infoPPos[i].likeBaseNoindel[nuc];

//     long double sumLogForNucs[4];  //sum of the logs for all but the nucleotide itself
//     long double sumLogForNucsC[4]; //sum of the logs for all but the nucleotide itself

//     long double likeBaseNoindel[4];  //log likelihood for each possible endogenous base
//     long double likeBaseNoindelC[4]; //log likelihood for each possible endogenous base


//     for(unsigned int nuc=0;nuc<4;nuc++){ 
// 	sumLogForNucs[nuc]     = 0.0;
// 	sumLogForNucsC[nuc]  = 0.0;
		  
// 	likeBaseNoindel[nuc]  = 0.0;
// 	likeBaseNoindelC[nuc] = 0.0;
//     }



//     if(singleCont){//if single contaminant need to marginalized over each contaminant
	     
// 	//Calling the endogenous base
// 	for(unsigned int nuce=0;nuce<4;nuce++){ //iterate over each possible endogenous base

		 
// 	    likeBaseNoindel[nuce]      = 0.0;

		 
// 	    for(unsigned int  nucc=0;nucc<4;nucc++){ //marginalize over each possible contaminant base A,C,G,T
		     		     
// 		//likeBaseNoindel[nuce]  +=  pow(10.0,infoPPos[i].likeBaseNoindelCont[nuce][nucc])*0.25;

// 		// const bool 
// 		// 			  const bool contIndelCurrent){
// 		if(endoIndelCurrent)
// 		    likeBaseNoindel[nuce]  =  oplusInit(likeBaseNoindel[nuce], infoPPos[i].likeBaseNoindelContNoBoundary[nuce][nucc] + log(0.25)/log(10) );		
// 		else
// 		    likeBaseNoindel[nuce]  =  oplusInit(likeBaseNoindel[nuce],           infoPPos[i].likeBaseNoindelCont[nuce][nucc] + log(0.25)/log(10) );
// 		// if(i==145)
// 		//     cout<<"E"<<(i+1)<<"\t"<<endoIndelCurrent<<contIndelCurrent<<"\tllik\t"<<dnaAlphabet[nuce]<<"\t"<<dnaAlphabet[nucc]<<"\t"<<(infoPPos[i].likeBaseNoindelCont[nuce][nucc])<<"\t"<<likeBaseNoindel[nuce]<<endl;
// 	    }
// 	    // if(i==146)
// 	    //     cout<<"E"<<(i+1)<<"\tllik\t"<<dnaAlphabet[nuce]<<"\t"<<likeBaseNoindel[nuce]<<endl;

// 	    //likeBaseNoindel[nuce] = log(likeBaseNoindel[nuce])/log(10);
// 	}

// 	//Calling the contaminant base
// 	for(unsigned int nucc=0;nucc<4;nucc++){ //iterate over each possible contaminant base
		 
// 	    likeBaseNoindelC[nucc]      = 0.0;

// 	    for(unsigned int nuce=0;nuce<4;nuce++){ //iterate over each possible endogenous base A,C,G,T
		     
// 		//likeBaseNoindelC[nucc]  += pow(10.0,infoPPos[i].likeBaseNoindelCont[nuce][nucc])*0.25;
// 		if(contIndelCurrent)
// 		    likeBaseNoindelC[nucc]  = oplusInit(likeBaseNoindelC[nucc], infoPPos[i].likeBaseNoindelContNoBoundary[nuce][nucc] + log(0.25)/log(10) );
// 		else
// 		    likeBaseNoindelC[nucc]  = oplusInit(likeBaseNoindelC[nucc],           infoPPos[i].likeBaseNoindelCont[nuce][nucc] + log(0.25)/log(10) );
		     
// 		// if(i==146)
// 		//cout<<"C"<<(i+1)<<"\t"<<endoIndelCurrent<<contIndelCurrent<<"\tllik\t"<<dnaAlphabet[nuce]<<"\t"<<dnaAlphabet[nucc]<<"\t"<<(infoPPos[i].likeBaseNoindelCont[nuce][nucc])<<"\t"<<likeBaseNoindelC[nucc]<<endl;
// 	    }
// 	    // if(i==146)
// 	    //     cout<<"C"<<(i+1)<<"\tllik\t"<<dnaAlphabet[nucc]<<"\t"<<likeBaseNoindelC[nucc]<<endl;
// 	    //likeBaseNoindelC[nucc] = log ( likeBaseNoindelC[nucc] )/log(10);
// 	}	     


	     
	     

//     }else{

// 	for(unsigned int nuc=0;nuc<4;nuc++){ //iterate over each possible endogenous base
// 	    if(endoIndelCurrent)
// 		likeBaseNoindel[nuc] = infoPPos[i].likeBaseNoindelNoBoundary[nuc];
// 	    else
// 		likeBaseNoindel[nuc] = infoPPos[i].likeBaseNoindel[nuc];
// 	}

//     }


//     //Begin determining best endogenous nucleotide

//     callBestNucleotideGivenLikelihood( bestNuc,
// 				       sumLogLikeAll,          // sum of all the logs
// 				       sumLogLikeOnlyBest,     // sum of the logs for the best
// 				       sumLogLikeAllButBest,   // sum of the logs for all but the best
// 				       sumLogForNucs,
// 				       likeBaseNoindel,
// 				       infoPPos,
// 				       i );

//     // if(i==146){
//     //     for(unsigned int nuc=0;nuc<4;nuc++){ //iterate over each possible endogenous base
//     // 	 cout<<"E\t"<<dnaAlphabet[nuc]<<"\tbest="<<dnaAlphabet[bestNuc]<<"\t"<<likeBaseNoindel[nuc]<<"\t"<<likeBaseNoindelC[nuc]<<endl;
//     // 	 //likeBaseNoindel[nuc] = infoPPos[i].likeBaseNoindel[nuc];
//     //     }		
//     // }

//     if(singleCont){//if single contaminant, call the contaminant

// 	callBestNucleotideGivenLikelihood( bestNucC,
// 					   sumLogLikeAllC,        // sum of all the logs
// 					   sumLogLikeOnlyBestC,   // sum of the logs for the best
// 					   sumLogLikeAllButBestC, // sum of the logs for all but the best
// 					   sumLogForNucsC,
// 					   likeBaseNoindelC,
// 					   infoPPos,
// 					   i );

// 	// if(i==146){
// 	// 	 for(unsigned int nuc=0;nuc<4;nuc++){ //iterate over each possible endogenous base
// 	// 	     cout<<"C\t"<<dnaAlphabet[nuc]<<"\tbest="<<dnaAlphabet[bestNucC]<<"\t"<<likeBaseNoindel[nuc]<<"\t"<<likeBaseNoindelC[nuc]<<endl;
// 	// 	     //likeBaseNoindel[nuc] = infoPPos[i].likeBaseNoindel[PHREDgeno];
// 	// 	 }		
// 	// }

//     }

//     // if(i==513){
//     // 	cout<<"REF1 "<<pos2phredgeno[   (i+1)   ].ref<<"\t"<<pos2phredgeno[   (i+1)   ].consensus<<endl;
//     // }



//     //cout<<"single "<<skipEndo<<"\t"<<skipCont<<"\t"<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<bestNuc<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<<sumLogLikeAllButBest<<"\t"<<sumLogLikeAll<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<endl;//<<endl;    
//     if(!skipEndo  && !skipCont){  //need to define both

// 	PHREDgeno toadd;

// 	if( (-10.0*(sumLogLikeAllButBest-sumLogLikeAll)) >= minQual){
// 	    genomeToPrint+=dnaAlphabet[bestNuc];
// 	}else{
// 	    genomeToPrint+="N";
// 	}

// 	int covEndoToPrint;
// 	if(endoIndelCurrent)
// 	    covEndoToPrint=infoPPos[i].covPerBaseNoBoundary[bestNuc];
// 	else
// 	    covEndoToPrint=infoPPos[i].covPerBase[bestNuc];

// 	long double valToPrintE= (-10*(sumLogLikeAllButBest-sumLogLikeAll)+0.0);
// 	if(valToPrintE==0.0) 
// 	    valToPrintE = 0.0;
	
	
	
// 	(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<< valToPrintE <<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<covEndoToPrint;//<<endl;


// 	for(int nuc=0;nuc<4;nuc++){
// 	    long double valToPrint =  (-10*(sumLogForNucs[nuc]-sumLogLikeAll));
// 	    if(valToPrint==0.0) 
// 		valToPrint = 0.0;

// 	    (*logToPrint)<<"\t"<<      valToPrint;
// 	    toadd.phred[nuc]  =        valToPrint;
// 	    toadd.perror[nuc] =    pow(10.0,(sumLogForNucs[nuc]-sumLogLikeAll));
// 	}
// 	(*logToPrint)<<endl;
// 	toadd.ref       = genomeRef[i];
// 	toadd.consensus = dnaAlphabet[bestNuc];



// 	if(singleCont){//if single contaminant, call the contaminant
		
// 	    if( (-10.0*(sumLogLikeAllButBestC-sumLogLikeAllC)) >= minQual){
// 		genomeToPrintC+=dnaAlphabet[bestNucC];
// 	    }else{
// 		genomeToPrintC+="N";
// 	    }

// 	    int covContToPrint;
// 	    if(contIndelCurrent)
// 		covContToPrint=infoPPos[i].covPerBaseNoBoundary[bestNucC];
// 	    else
// 		covContToPrint=infoPPos[i].covPerBase[bestNucC];
// 	    long double valToPrintC= (-10*(sumLogLikeAllButBestC-sumLogLikeAllC)+0.0);
// 	    if(valToPrintC==0.0) 
// 		valToPrintC = 0.0;
// 	    (*logToPrintC)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNucC]<<"\t"<< valToPrintC <<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<covContToPrint;//<<endl;
	 
// 	    for(int nuc=0;nuc<4;nuc++){
// 		long double valToPrint= (-10*(sumLogForNucsC[nuc]-sumLogLikeAllC)); 
// 		if(valToPrint==0.0) 
// 		    valToPrint = 0.0;

// 		(*logToPrintC)<<"\t"<<      (valToPrint);
// 		toadd.phredC[nuc]  =        (valToPrint);
// 		toadd.perrorC[nuc] = pow(10.0,(sumLogForNucsC[nuc]-sumLogLikeAllC));
// 	    }
// 	    (*logToPrintC)<<endl;
// 	}
	
// 	pos2phredgeno[   (i+1)   ] = toadd;

//     }

    
//     if(skipEndo  && !skipCont){  //need to define contamination

// 	if(singleCont){//if single contaminant, call the contaminant
		
// 	    if( (-10.0*(sumLogLikeAllButBestC-sumLogLikeAllC)) >= minQual){
// 		genomeToPrintC+=dnaAlphabet[bestNucC];
// 	    }else{
// 		genomeToPrintC+="N";
// 	    }

// 	    int covContToPrint;
// 	    if(contIndelCurrent)
// 		covContToPrint=infoPPos[i].covPerBaseNoBoundary[bestNucC];
// 	    else
// 		covContToPrint=infoPPos[i].covPerBase[bestNucC];

// 	    long double valToPrintC= (-10*(sumLogLikeAllButBestC-sumLogLikeAllC)+0.0);
// 	    if(valToPrintC==0.0) 
// 		valToPrintC = 0.0;

// 	    (*logToPrintC)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNucC]<<"\t"<<valToPrintC<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<covContToPrint;//<<endl;
	 
	    
// 	    for(int nuc=0;nuc<4;nuc++){
// 		long double valToPrint= (-10*(sumLogForNucsC[nuc]-sumLogLikeAllC)); 
// 		if(valToPrint==0.0) 
// 		    valToPrint = 0.0;

// 		(*logToPrintC)<<"\t"<<                        (valToPrint);
// 		pos2phredgeno[   (i+1)   ].phredC[nuc]  =     (valToPrint);
// 		pos2phredgeno[   (i+1)   ].perrorC[nuc] = pow(10.0,(sumLogForNucsC[nuc]-sumLogLikeAllC));
// 	    }
// 	    (*logToPrintC)<<endl;
// 	}

//     }

    
//     if(!skipEndo  && skipCont){  //need not to define endogenous
	
// 	PHREDgeno toadd;

// 	if( (-10.0*(sumLogLikeAllButBest-sumLogLikeAll)) >= minQual){
// 	    genomeToPrint+=dnaAlphabet[bestNuc];
// 	}else{
// 	    genomeToPrint+="N";
// 	}


// 	int covEndoToPrint;
// 	if(endoIndelCurrent)
// 	    covEndoToPrint=infoPPos[i].covPerBaseNoBoundary[bestNuc];
// 	else
// 	    covEndoToPrint=infoPPos[i].covPerBase[bestNuc];

	
// 	long double valToPrintE= (-10*(sumLogLikeAllButBest-sumLogLikeAll)+0.0);
// 	if(valToPrintE==0.0) 
// 	    valToPrintE = 0.0;


// 	(*logToPrint)<<(i+1)<<"\t"<<genomeRef[i]<<"\t"<<dnaAlphabet[bestNuc]<<"\t"<<valToPrintE<<"\t"<<infoPPos[i].mapqAvg<<"\t"<<infoPPos[i].cov<<"\t"<<covEndoToPrint;//<<endl;


// 	for(int nuc=0;nuc<4;nuc++){
// 	    long double valToPrint= (-10*(sumLogForNucs[nuc]-sumLogLikeAll)); 
// 	    if(valToPrint==0.0) 
// 		valToPrint = 0.0;

// 	    (*logToPrint)<<"\t"<<      valToPrint;
// 	    toadd.phred[nuc]     =     valToPrint;
// 	    toadd.perror[nuc] = pow(10.0,(sumLogForNucs[nuc]-sumLogLikeAll));
// 	}
// 	(*logToPrint)<<endl;
// 	toadd.ref       = genomeRef[i];
// 	toadd.consensus = dnaAlphabet[bestNuc];

// 	pos2phredgeno[   (i+1)   ] = toadd;

//     }



//     // if(i==513){
//     // 	cout<<"REF2 "<<pos2phredgeno[   (i+1)   ].ref<<"\t"<<pos2phredgeno[   (i+1)   ].consensus<<endl;
//     // }



//     // cout<<(i+1)<<"\t"<<toadd.ref<<"\t"<<toadd.consensus<<endl;

// }


// //! A method that calls each subroutine for call the endogenous (and contaminant) after the BAM file was read
// /*!
//   This method is called by printLogAndGenome() and just prints a simple line saying there was no coverage.

//   \param sizeGenome: the actual (biological) size of the mitochondrial. The reference can be longer if the first base pairs were copied at the end
//   \param infoPPos: The vector of structure populated by the bam reader.
//   \param outSeq:  String containing the filename of the endogenous sequence
//   \param outLog:  String containing the filename of the endogenous log file
//   \param genomeRef : The reference genome 
//   \param minQual: PHRED quality threshold, beyong this we print N instead of the base for single nucleotides
//   \param nameMT: Name of the produced fasta record in the fasta file for the endogenous sample
//   \param singleCont: Boolean as to we assume that we have a single contaminant or not
//   \param outSeqC:  String containing the filename of the contaminant sequence
//   \param outLogC:  String containing the filename of the contaminat log file
//   \param outSeqCflag: Boolean flag to say whether we print to the contaminant fasta sequence or not
//   \param outLogCflag: Boolean flag to say whether we print to the contaminant log file or not
//   \param nameMTC:  Name of the produced fasta record in the fasta file for the endogenous sample
// */				
// void  printLogAndGenome(const int sizeGenome,
// 			const vector<singlePosInfo> & infoPPos,
// 			const string outSeq,
// 			const string outLog, 
// 			const string genomeRef,
// 			const int minQual,
// 			const string nameMT,
// 			const bool singleCont,
// 			const string outSeqC,
// 			const string outLogC,
// 			const bool   outSeqCflag,
// 			const bool   outLogCflag,
// 			const string nameMTC){


//     // cerr<<outSeq<<"\t"<<outLog<<"\t#"<<outSeqC<<"#\t"<<outLogC<<"#\t"<<outSeqCflag<<"#\t"<<outLogCflag<<"#\t"<<endl;
//     // 
//     ofstream outSeqFP ;
//     ofstream outLogFP;

//     ofstream outSeqFPC;
//     ofstream outLogFPC;


//     outSeqFP.open(outSeq.c_str());

//     if (!outSeqFP.is_open()){
// 	cerr << "Unable to write to seq file "<<outSeq<<endl;
// 	exit(1);
//     }

//     outLogFP.open(outLog.c_str());

//     if (!outLogFP.is_open()){
// 	cerr << "Unable to write to qual file "<<outLog<<endl;
// 	exit(1);
//     }


//     if(outSeqCflag){
// 	outSeqFPC.open(outSeqC.c_str());
	
// 	if (!outSeqFPC.is_open()){
// 	    cerr << "Unable to write to seq file "<<outSeqC<<endl;
// 	    exit(1);
// 	}
//     }

//     if(outLogCflag){
// 	outLogFPC.open(outLogC.c_str());

// 	if (!outLogFPC.is_open()){
// 	    cerr << "Unable to write to qual file "<<outLogC<<endl;
// 	    exit(1);
// 	}

//     }

//     string genomeToPrint="";
//     string genomeToPrintC="";

//     stringstream logToPrint;
//     stringstream logToPrintC;


//     logToPrint<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt\n";
//     if(outLogCflag){
// 	logToPrintC<<"pos\trefBase\tbase\tqual\tavgmapq\tcov\tsupp\tpa\tpc\tpg\tpt\n";
//     }

//     //genomeRef
//     bool previousEndoDel=false;
//     bool previousContDel=false;

//     for(int i=0;i<sizeGenome;i++){
// 	bool skipEndo=false;
// 	bool skipCont=false;

// 	//if position ahead or behind had an indel
// 	//need if we define IGNOREINDELBOUND 
// 	//to call nuc. around the indel
// 	bool  endoIndel=false;
// 	bool  contIndel=false;


// 	if( !isResolvedDNA(genomeRef[i]) ){
// 	    cerr<<"Skipping position "<<i<<" due to unknown base, found  = "<<genomeRef[i]<<endl;
// 	}
	

// 	//need to check for indel in the position in front
// #ifdef	IGNOREINDELBOUND
// 	if(i<(sizeGenome-1)){
// 	    bool  endoIndelD=false;
// 	    bool  contIndelD=false;
// 	    bool  endoIndelI=false;
// 	    bool  contIndelI=false;
	
// 	    deletionInSample(i+1,
// 			     genomeRef,
// 			     infoPPos,
// 			     singleCont,
// 			     minQual,
// 			     genomeToPrint,
// 			     genomeToPrintC,
// 			     &logToPrint,
// 			     &logToPrintC,
// 			     skipEndo,
// 			     skipCont,
// 			     true,
// 			     endoIndelD,
// 			     contIndelD);

	    
// 	    insertionInSample(i,
// 	    		      genomeRef,
// 	    		      infoPPos,
// 	    		      singleCont,
// 			      minQual,
// 	    		      genomeToPrint,
// 	    		      genomeToPrintC,
// 	    		      &logToPrint,
// 	    		      &logToPrintC,
// 	    		      true,
// 	    		      endoIndelI,
// 	    		      contIndelI);
	    
// 	    endoIndel = (endoIndelI || endoIndelD);
// 	    contIndel = (contIndelI || contIndelD);

// 	}
// #endif

// 	//There are 4 possibilities:
// 	//1) There is a deletion in the sample (or insertion in the reference)
// 	//2) There is no coverage
// 	//3) There are no bases
// 	//4) There is a insertion in the sample (or deletion in the reference)

	
// 	bool  endoIndelCurrent=false;
// 	bool  contIndelCurrent=false;
// 	bool  endoIndelCurrentD=false;
// 	bool  contIndelCurrentD=false;
// 	bool  endoIndelCurrentI=false;
// 	bool  contIndelCurrentI=false;


// 	// cout<<i<<"\t"<<endoIndel<<"\t"<<previousEndoDel<<"\t"<<endoIndelCurrent<<"\t"<<contIndel<<"\t"<<previousContDel<<"\t"<<contIndelCurrent<<endl;

// 	/////////////////////////////////////////
// 	//                                     //
// 	//        Deletion in the sample       //
// 	//                                     //
// 	/////////////////////////////////////////
// 	deletionInSample(i,
// 			 genomeRef,
// 			 infoPPos,
// 			 singleCont,
// 			 minQual,
// 			 genomeToPrint,
// 			 genomeToPrintC,
// 			 &logToPrint,
// 			 &logToPrintC,
// 			 skipEndo,
// 			 skipCont,
// 			 false,
// 			 endoIndelCurrentD,
// 			 contIndelCurrentD);

// 	// cout<<i<<"\t"<<endoIndel<<"\t"<<previousEndoDel<<"\t"<<endoIndelCurrent<<"\t"<<contIndel<<"\t"<<previousContDel<<"\t"<<contIndelCurrent<<endl;

// 	endoIndelCurrent = ( endoIndelCurrentD || endoIndelCurrentI);
// 	contIndelCurrent = ( contIndelCurrentD || contIndelCurrentI);

// 	if(skipCont && skipEndo){
// 	    previousEndoDel=endoIndelCurrent;
// 	    previousContDel=contIndelCurrent;
// 	    continue;
// 	}

// 	/////////////////////////////////////////
// 	//                                     //
// 	//            No coverage              //
// 	//                                     //
// 	/////////////////////////////////////////
// 	noCoverage(i,
// 		   genomeRef,
// 		   infoPPos,
// 		   outLogCflag,	
// 		   genomeToPrint,
// 		   genomeToPrintC,
// 		   &logToPrint,
// 		   &logToPrintC,
// 		   skipEndo,
// 		   skipCont);
// 	if(skipCont && skipEndo){
// 	    previousEndoDel=endoIndelCurrent;
// 	    previousContDel=contIndelCurrent;
// 	    continue;
// 	}
	     

// 	/////////////////////////////////////////
// 	//                                     //
// 	//      Calling single nucleotides     //
// 	//                                     //
// 	/////////////////////////////////////////
// 	callSingleNucleotide(i,
// 			     genomeRef,
// 			     infoPPos,
// 			     singleCont,
// 			     minQual,
// 			     genomeToPrint,
// 			     genomeToPrintC,
// 			     &logToPrint,
// 			     &logToPrintC,
// 			     skipEndo,
// 			     skipCont,
// 			     endoIndel ||  previousEndoDel,
// 			     contIndel ||  previousContDel
// 			     );





// 	/////////////////////////////////////////
// 	//                                     //
// 	//      Insertions in the sample       //
// 	//                                     //
// 	/////////////////////////////////////////

// 	insertionInSample(i,
// 			  genomeRef,
// 			  infoPPos,
// 			  singleCont,
// 			  minQual,
// 			  genomeToPrint,
// 			  genomeToPrintC,
// 			  &logToPrint,
// 			  &logToPrintC,
// 			  false,
// 			  endoIndelCurrentI,
// 			  contIndelCurrentI);

// 	endoIndelCurrent = ( endoIndelCurrentD || endoIndelCurrentI);
// 	contIndelCurrent = ( contIndelCurrentD || contIndelCurrentI);

// 	previousEndoDel=endoIndelCurrent;
// 	previousContDel=contIndelCurrent;


//     }//end for each position in the genome


//      string genomeToPrintCopy="";
//      for(unsigned int i=1;i<(genomeToPrint.size()+1);i++){
// 	 genomeToPrintCopy+=genomeToPrint[i-1];
// 	 if(i!=0 && (i%80 == 0))
// 	     genomeToPrintCopy+="\n";
//      }

//      string genomeToPrintCopyC="";
//      for(unsigned int i=1;i<(genomeToPrintC.size()+1);i++){
// 	 genomeToPrintCopyC+=genomeToPrintC[i-1];
// 	 if(i!=0 && (i%80 == 0))
// 	     genomeToPrintCopyC+="\n";
//      }



//      outSeqFP<<nameMT+"\n"+genomeToPrintCopy+"\n";
//      outSeqFP.close() ;

//      if(outSeqCflag)
//      if(singleCont){//if single contaminant need to marginalized over each contaminant
// 	 outSeqFPC<<nameMTC+"\n"+genomeToPrintCopyC+"\n";
// 	 outSeqFPC.close() ;
//      }

//      outLogFP<<logToPrint.str();
//      outLogFP.close();

//      if(outLogCflag)
//      if(singleCont){//if single contaminant need to marginalized over each contaminant
// 	 outLogFPC<<logToPrintC.str();
// 	 outLogFPC.close();
//      }
// }





// //! A method that computes by how much a string overlaps another using the prefix
// /*!
  
//   \param s1: first string
//   \param s2: second string
//   \return The fraction that a string overlaps another using the prefix
// */
// inline long double string2sub(const string & s1,const string & s2){
//     if(s1==s2)
// 	return 1.0;
//     // else
//     //  	return 0.0;
//     int minStrSize = int(min(s1.size(),s2.size()));
//     int maxStrSize = int(max(s1.size(),s2.size()));
//     int ident=0;

//     for(int i=0;i<minStrSize;i++){
// 	if(s1[i] == s2[i]){
// 	    ident++;
// 	}else
// 	    break;
//     }

//     return (long double)(ident) / (long double)(maxStrSize);
// }




// //! A method that computes the probability of observing a given insertion 
// /*!
//   This method is called by insertionInSample to compute the probability of 
//   observing a given insertion given an insertion in the model.
  
//   \param modelIns:   The string from the "model"
//   \param obserIns:   The observed string
//   \return The probability of observing obserIns given modelIns
// */
// inline long double insertPair2Prob(const string & modelIns,const string & obserIns){
//     long double overlapFrac = string2sub(modelIns,obserIns);//1 if identical, 0 otherwise
//     //if(modelIns == obserIns)
//     //                     //correct                        + incorrect
//     long double toReturn = overlapFrac*(1.0-INDELERRORPROB) + (1.0-overlapFrac)*INDELERRORPROB;
//     //cout<<"insertPair2Prob\t"<<modelIns<<"\t"<<obserIns<<"\t"<<toReturn<<endl;
//     return toReturn;
//     // else
//     // 	return (    INDELERRORPROB);
// }
// // inline long double insertPair2Prob(const string & modelIns,const string & obserIns){
// //     if(modelIns == obserIns)
// // 	return (1.0-INDELERRORPROB);
// //     else
// // 	return (    INDELERRORPROB);
// // }


//unsigned int posFound;
//vector<positionInfo> * positionsInfoFound;
// map<string, map<unsigned int,contaminationInfo> > contFreq;

//! Object to iterate over read at a single position
/*!
  This object is called by iterateOverReads() and the Visit() method is used for each position

*/				
class MyPileupVisitor : public PileupVisitor {
  
public:
  MyPileupVisitor(const RefVector& references,
		  Fasta * fastaReference,
		  // vector<singlePosInfo> * infoPPos,
		  // const int sizeGenome,
		  vector<MistarParser * > * vectorOfMP,
		  unsigned int coordFirst,
		  const bool ignoreMQ
		  // const long double contaminationPrior,
		  //const bool singleCont
)
    : PileupVisitor()
      // , m_references(references)
    , m_fastaReference(fastaReference)
    // , m_infoPPos(infoPPos)
    // , sizeGenome(sizeGenome)
    , m_vectorOfMP(vectorOfMP)
    , ignoreMQ(ignoreMQ)
    // , contaminationPrior(contaminationPrior)
    // , singleCont(singleCont)
  { 
     
    
    initFiles(*vectorOfMP,
	      // atLeastOneHasData,
	      hasData,
	      popSizePerFile,
	      vecAlleleRecords,
	      chr1,
	      coordFirst);

    hasCoordinate = vector<bool>(m_vectorOfMP->size(),true);//dummy value


  }
  ~MyPileupVisitor(void) { }
  
  // PileupVisitor interface implementation
public:

  void Visit(const PileupPosition& pileupData) {
    //cout<<"visit"<<endl;
    char referenceBase = 'N';
	    
    ///unsigned int posAlign = pileupData.Position+1;
    //    int posVector=int(pileupData.Position)%sizeGenome;
	    

    unsigned int posAlign = pileupData.Position+1;
    // int posVector=int(pileupData.Position)%sizeGenome;
	    
    // if( (posAlign%100) == 0){
    // 	cerr<<"pos  = "<<posAlign<<endl;
    // }
    // cout<<endl<<"pos = "<<posAlign<<"\t"<<posVector;
	    

    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
	cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file at chr "<<pileupData.RefId<<" position "<<(posAlign-1) << endl;
	exit(1);
    }


    //re-iterate over each read to compute the likelihood for each insert (or pair of inserts)
    for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){				


		
	if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
	    pileupData.PileupAlignments[i].IsNextInsertion &&
	    (pileupData.PileupAlignments[i].InsertionLength>0)){ //has insertion		    
	    continue;//ignore
	}

	if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
	    pileupData.PileupAlignments[i].IsNextInsertion &&
	    (pileupData.PileupAlignments[i].InsertionLength == 0)){ //has deletion
	    continue;
	}

	//base that was read
	char b   = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
	//quality score
	char q   = pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual;
	int  m   = int(pileupData.PileupAlignments[i].Alignment.MapQuality);
	
	//skip unresolved
	if(b == 'N')
	    continue;	    


	// BEGIN DEAMINATION COMPUTATION
	//zero base distance to the 5p/3p end
	int dist5p=-1;
	int dist3p=-1;

	if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){
	    dist5p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
	    dist3p = pileupData.PileupAlignments[i].PositionInAlignment;
	}else{
	    dist5p = pileupData.PileupAlignments[i].PositionInAlignment;
	    dist3p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
	}
		    		    
	probSubstition * probSubMatchToUseEndo = &defaultSubMatch ;
	probSubstition * probSubMatchToUseCont = &defaultSubMatch ;


	if(dist5p <= (int(sub5p.size()) -1)){
	    probSubMatchToUseEndo = &sub5p[  dist5p ];			

	}

	if(dist5p <= (int(sub5pC.size()) -1)){
	    probSubMatchToUseCont = &sub5pC[ dist5p ];
	}
		    
	if(dist3p <= (int(sub3p.size()) -1)){
	    probSubMatchToUseEndo = &sub3p[  dist3p ];
	}
		    
	if(dist3p <= (int(sub3pC.size()) -1)){
	    probSubMatchToUseCont = &sub3pC[ dist3p ];
	}

	//we have substitution probabilities for both... take the closest
	if(dist5p <= (int(sub5p.size()) -1) &&
	   dist3p <= (int(sub3p.size()) -1) ){
		    
	    if(dist5p < dist3p){
		probSubMatchToUseEndo = &sub5p[  dist5p ];
	    }else{
		probSubMatchToUseEndo = &sub3p[  dist3p ];
	    }
		    
	}

	if(dist5p <= (int(sub5pC.size()) -1) &&
	   dist3p <= (int(sub3pC.size()) -1) ){
		    
	    if(dist5p < dist3p){
		probSubMatchToUseCont = &sub5pC[ dist5p ];
	    }else{
		probSubMatchToUseCont = &sub3pC[ dist3p ];
	    }
		    
	}

	// END DEAMINATION COMPUTATION


	
	// vector<bool> hasData;
	// vector<AlleleRecords *> vecAlleleRecords;
	// vector<bool>  hasCoordinate (m_vectorOfMP->size(),true);//dummy value
	bool stayLoop=true;


	while(stayLoop){

#ifdef DEBUGMST
	    cout<<"posAlign "<<posAlign<<endl;
#endif

	    bool allHaveCoordinate=true;

	    for(unsigned int i=0;i<m_vectorOfMP->size();i++){ 
		//cout<<"hasData["<<i<<"]"<<hasData[i]<<endl;
		if(hasData[i]){		
		    if(posAlign  != vecAlleleRecords[i]->coordinate){
			allHaveCoordinate=false;
		    }
		}else{
		    stayLoop=false;
		    break;
		}
	    }

	
	    //we print
	    if(allHaveCoordinate){
#ifdef DEBUGMST
		cout<<"same posAlign "<<posAlign<<endl;
#endif

		cout<<"i="<<i<<"\t"<<posAlign<<"\t"<<b<<"\t"<<referenceBase<<"\t"<<int(q)<<endl;	    
		printAllele(*m_vectorOfMP,
			    hasData,
			    hasCoordinate,
			    popSizePerFile,
			    vecAlleleRecords,
			    chr1,
			    posAlign,
			    false);


		// 	seekdata:
		allHaveCoordinate=false;
	    
		for(unsigned int i=0;i<m_vectorOfMP->size();i++){ 
		 
		    if(!hasData[i] ){
			cerr<<"Invalid state"<<endl;
			exit(1);
		    }
		    hasData[i]  =  m_vectorOfMP->at(i)->hasData();
		    if(hasData[i]){
			vecAlleleRecords[i] = m_vectorOfMP->at(i)->getData() ;
		    }else{
			stayLoop=false;
			break;
		    }
		
		}


		// //all have add getData called, we need to reposition to the maximum coord
		// bool needToSetCoord=true;
		// for(unsigned int i=0;i<m_vectorOfMP->size();i++){ 
		
		//     if(needToSetCoord){
		// 	posAlign  = vecAlleleRecords[i]->coordinate;
		// 	needToSetCoord=false;
		//     }else{
		// 	posAlign  = max(posAlign,vecAlleleRecords[i]->coordinate);
		//     }

		// }

		continue;

	    }else{
		// cerr<<"Invalid state"<<endl;
		// return 1;  
	    
	    
		for(unsigned int i=0;i<m_vectorOfMP->size();i++){ 
#ifdef DEBUGMST
		    cout<<"coord["<<i<<"] "<< vecAlleleRecords[i]->coordinate<<endl;
#endif

		     if(posAlign  < vecAlleleRecords[i]->coordinate){ //The BAM file is behind the MST files
			 // posAlign = vecAlleleRecords[i]->coordinate;
			 // continue;
			 stayLoop=false;
			 break;
		     }

		    if(posAlign == vecAlleleRecords[i]->coordinate){ //fine
		    }

		    if(posAlign >  vecAlleleRecords[i]->coordinate){ //running behind
			hasData[i]  =   m_vectorOfMP->at(i)->hasData();
			if(hasData[i]){
			    vecAlleleRecords[i] = m_vectorOfMP->at(i)->getData() ;
			}else{
			    stayLoop=false;
			    break;
			}
		    }
		
		}

	    
	    }//end different coord

	
	}//end stay loop



	//}//end, no insert
    }//end for each read
	




    // cout<<"end visit"<<endl;



    
  }//end visit()
        
private:

    vector<bool> hasData;
    vector<int> popSizePerFile;
    vector<AlleleRecords *> vecAlleleRecords;
    string chr1;
    unsigned int coordCurrent;
    vector<bool>  hasCoordinate;

    RefVector m_references;
    Fasta * m_fastaReference;
    vector<MistarParser * > * m_vectorOfMP;
    // vector<singlePosInfo> * m_infoPPos;
    // int sizeGenome;
    bool ignoreMQ;
    // long double contaminationPrior;
    // bool singleCont;
    
    //        ostream*  m_out;
};

























//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){

    for(int i=0;i<2;i++){
        likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )    /log(10);         
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);

	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
    }


    //Computing for quality scores 2 and up
    for(int i=2;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )     /log(10);          
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10);

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;
    }


    //Adding mismapping probability
    for(int m=0;m<MAXMAPPINGQUAL;m++){

	double incorrectMappingProb   =     pow(10.0,m/-10.0); //m
	double correctMappingProb     = 1.0-pow(10.0,m/-10.0); //1-m
	
	probMapping[m]    = correctMappingProb;    //1-m
	probMismapping[m] = incorrectMappingProb;  //m

#ifdef DEBUG1
	cerr<<"m\t"<<m<<"\t"<<incorrectMappingProb<<"\t"<<correctMappingProb<<endl;
#endif
	
    	for(int i=0;i<2;i++){
    	    likeMatchMQ[m][i]           = log(  correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0   )/log(10);         
    	    likeMismatchMQ[m][i]        = log(  correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0   )/log(10);
    	    likeMatchProbMQ[m][i]       = correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0;
    	    likeMismatchProbMQ[m][i]    = correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0;
    	}


    	//Computing for quality scores 2 and up
    	for(int i=2;i<MAXMAPPINGQUAL;i++){
	    //  (1-m)(1-e) + m/4  = 1-m-e+me +m/4  = 1+3m/4-e+me
    	    likeMatchMQ[m][i]         = log(  correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0    )/log(10);    
	    //  (1-m)(e/3) + m/4  = e/3 -me/3 + m/4
    	    likeMismatchMQ[m][i]      = log(  correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0    )/log(10);    
	    
    	    likeMatchProbMQ[m][i]           = correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0;
    	    likeMismatchProbMQ[m][i]        = correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0;
    	}


#ifdef DEBUG1
    	for(int i=0;i<MAXMAPPINGQUAL;i++){
	    cerr<<"m\t"<<m<<"\t"<<i<<"\t"<<likeMatchMQ[m][i]<<"\t"<<likeMismatchMQ[m][i]<<"\t"<<likeMatchProbMQ[m][i]<<"\t"<<likeMismatchProbMQ[m][i]<<endl;
	}
#endif

    }

}



//! Main method
/*!
  The main:
    calls initScores(), 
    captures the arguments
    reads the deamination and Illumina error profiles
    calls     iterateOverReads() to populated infoPPos and printLogAndGenome to print the information contained therein
    if we use deamination or length priors:
       Call computePriorOnReads() to compute the prob. that each read is endogenous
       Recalls iterateOverReads() to populated infoPPos and printLogAndGenome to print the information contained therein using the new prior on the reads
*/
int main (int argc, char *argv[]) {
    setlocale(LC_ALL, "POSIX");
    int sizeGenome=0;
    string output  = "/dev/stdout";
    string outlog  = "/dev/stderr";
    // string nameMT  = "MT";

    // string outSeqC  = "";
    // string outLogC  = "";
    // bool   outSeqCflag  = false;
    // bool   outLogCflag  = false;


    // string nameMTC  = "MTc";
    // bool userWantsContProduced = false;
  
    int minQual                = 0;
    bool ignoreMQ              = false;


    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
    initScores();
    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////


    long double locatione=0.0;
    long double locationc=0.0;
    long double scalee   =0.0;
    long double scalec   =0.0;



    ////////////////////////////////////
    // BEGIN Parsing arguments        //
    ////////////////////////////////////

    //    return 1;
    string errFile     = getCWD(argv[0])+"illuminaProf/error.prof";
    string deam5pfreqE = getCWD(argv[0])+"deaminationProfile/none.prof";
    string deam3pfreqE = getCWD(argv[0])+"deaminationProfile/none.prof";
    string deam5pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";
    string deam3pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";

    // substitutionRates freqIlluminaError;
    vector<substitutionRates>    deam5PsubE;
    vector<substitutionRates>    deam3PsubE;
    vector<substitutionRates>    deam5PsubC;
    vector<substitutionRates>    deam3PsubC;

    // bool deamread                  = false;
    // bool useLengthPrior            = false;
    // long double contaminationPrior       = 0.5;
    // long double contaminationPriorNoDeam = 0.0;

    // bool singleCont                = false;
    // bool specifiedContPrior        = false;
    // bool specifiedLoce             = false;
    // bool specifiedLocc             = false;
    // bool specifiedScalee           = false;
    // bool specifiedScalec           = false;


    const string usage=string("\nThis program takes an aligned BAM file for a mitonchondria and calls a\nconsensus for the endogenous material\n\n\t"+
			      string(argv[0])+			      
			      " [options] [fasta file] [bam file] [region to use] [freq for pop1] [freq for pop2] ... "+"\n\n"+
			      
			      "The BAM file has to be indexed and the region has to be in chr:start-end format\n\n"+
			      "\n\tOutput options:\n"+	
			      "\t\t"+"-out  [output]" +"\t\t"+"Output file (default: stdout)"+"\n"+
			      "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      // "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      // "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+

			      "\n\tDeamination options:\n"+				      
			      "\t\t"+"-deam5p  [.prof file]" +"\t\t"+"5p deamination frequency for the endogenous  (default: "+deam5pfreqE+")"+"\n"+
			      "\t\t"+"-deam3p  [.prof file]" +"\t\t"+"3p deamination frequency for the endogenous  (default: "+deam3pfreqE+")"+"\n"+
			      "\t\t"+"-deam5pc [.prof file]" +"\t\t"+"5p deamination frequency for the contaminant (default: "+deam5pfreqC+")"+"\n"+
			      "\t\t"+"-deam3pc [.prof file]" +"\t\t"+"3p deamination frequency for the contaminant (default: "+deam3pfreqC+")"+"\n"+


			      // "\t\t"+"-deamread" +"\t\t\t"+"Set a prior on reads according to their deamination pattern (default: "+ booleanAsString(deamread) +")"+"\n"+
			      // "\t\t"+"-cont    [cont prior]"+"\t\t"+"If the -deamread option is specified, this is the contamination prior (default: "+ stringify(contaminationPrior) +")"+"\n"+
			      // "\t\t"+"                  "+"\t\t"+"Otherwise the contamination prior will be  "+ stringify(contaminationPriorNoDeam) +""+"\n\n"+

			      // "\t\t"+"-single"+"\t\t\t\t"+"Try to determine the contaminant under the assumption that there is a single\n\t\t\t\t\t\tone  (default: "+ booleanAsString(singleCont) +")"+"\n\n"+

			      // "\t\tIf the -single option is used, the following are available:\n"+
			      // "\t\t\t"+"-seqc  [fasta file]" +"\t"+"Output contaminant as fasta file (default: none)"+"\n"+
			      // "\t\t\t"+"-logc  [log file]" +"\t"+"Output contaminant as log  (default: none)"+"\n"+
			      // "\t\t\t"+"-namec [name]" +"\t\t"  +"Name of contaminant sequence (default "+nameMTC+") "+"\n"+




			      // "\n\tLength options:\n"+				      
			      // "\t\t"+"--loce"+  "\t\t\t\t"+"Location for lognormal dist for the endogenous sequences (default none)"+"\n"+
			      // "\t\t"+"--scalee"+"\t\t\t"+"Scale for lognormal dist for the endogenous sequences (default none)"+"\n"+
			      // "\t\t"+"--locc"+  "\t\t\t\t"+"Location for lognormal dist for the contaminant sequences (default none)"+"\n"+
			      // "\t\t"+"--scalec"+"\t\t\t"+"Scale for lognormal dist for the contaminant sequences (default none)"+"\n"+

			      // "\n\tComputation options:\n"+	
			      // "\t\t"+"-nomq" +"\t\t\t\t"+"Ignore mapping quality (default: "+booleanAsString(ignoreMQ)+")"+"\n"+
			      // "\t\t"+"-err" +"\t\t\t\t"+"Illumina error profile (default: "+errFile+")"+"\n"+
			      "\t\t"+"--phred64" +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+

			      // "\n\tReference options:\n"+	
			      // "\t\t"+"-l [length]" +"\t\t\t"+"Actual length of the genome used for"+"\n"+
			      // "\t\t"+"  " +"\t\t\t\t"+"the reference as been wrapped around"+"\n"+
			      // "\t\t"+"  " +"\t\t\t\t"+"by default, the length of the genome will be used "+"\n"+
			      
			      "");
			      
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<usage<<endl;
	return 1;
    }


    int lastOpt=1;

    for(int i=1;i<(argc);i++){ 

	if(string(argv[i])[0] != '-'  ){
	    lastOpt=i;
	    break;
	}

	if(string(argv[i]) == "--phred64"  ){
	    offsetQual=64;
	    continue;
	}


	// if(string(argv[i]) == "--loce" ){
	//     locatione =destringify<long double>(argv[i+1]);
	//     i++;
	//     specifiedLoce=true;
	//     continue;
	// }

	// if(string(argv[i]) == "--scalee" ){
	//     scalee =destringify<long double>(argv[i+1]);
	//     i++;
	//     specifiedScalee=true;
	//     continue;
	// }

	// if(string(argv[i]) == "--locc" ){
	//     locationc =destringify<long double>(argv[i+1]);
	//     i++;
	//     specifiedLocc=true;
	//     continue;
	// }

	// if(string(argv[i]) == "--scalec" ){
	//     scalec =destringify<long double>(argv[i+1]);
	//     i++;
	//     specifiedScalec=true;
	//     continue;
	// }

	// if(string(argv[i]) == "-err"  ){
	//     errFile=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

	if(string(argv[i]) == "-deam5p"  ){
	    deam5pfreqE=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-deam3p"  ){
	    deam3pfreqE=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-deam5pc"  ){
	    deam5pfreqC=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-deam3pc"  ){
	    deam3pfreqC=string(argv[i+1]);
	    i++;
	    continue;
	}

	// if(string(argv[i]) == "-deamread"  ){
	//     deamread=true;
	//     continue;
	// }


	// if(string(argv[i]) == "-single"  ){
	//     singleCont=true;
	//     continue;
	// }


	// if(string(argv[i]) ==  "-cont" ){
	//     contaminationPrior=destringify<long double>(argv[i+1]);
	//     specifiedContPrior=true;
	//     i++;
	//     continue;
	// }

	// if(strcmp(argv[i],"-deam5") == 0 ){
	//     deam5File=string(argv[i+1]);
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-nomq" ){
	//     ignoreMQ=true;
	//     continue;
	// }

	if(string(argv[i]) == "-out" ){
	    output=string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-log" ){
	    outlog=string(argv[i+1]);
	    i++;
	    continue;
	}

	// if(string(argv[i]) == "-name" ){
	//     nameMT=string(argv[i+1]);
	//     i++;
	//     continue;
	// }


	// if(string(argv[i]) == "-seqc" ){
	//     outSeqC=string(argv[i+1]);
	//     outSeqCflag = true;
	//     userWantsContProduced=true;
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-logc" ){
	//     outLogC=string(argv[i+1]);
	//     outLogCflag = true;
	//     userWantsContProduced=true;
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-namec" ){
	//     nameMTC=string(argv[i+1]);
	//     userWantsContProduced=true;
	//     i++;
	//     continue;
	// }

	// if(string(argv[i]) == "-qual" ){
	//     minQual=destringify<int>(argv[i+1]);
	//     i++;
	//     continue;
	// }



	// if(string(argv[i]) == "-l" ){
	//     sizeGenome=atoi(argv[i+1]);
	//     i++;
	//     continue;
	// }


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    if(output == outlog){
	cerr<<"Error: The sequence output is the same as the log"<<endl;
	return 1;	
    }


    vector<string>  filesAlFreq;
    string fastaFile    = string(argv[lastOpt]);//fasta file
    string bamfiletopen = string(argv[lastOpt+1]);//bam file
    //string fastaFile    = string(argv[argc-2]);//fasta file
    string region       = string(argv[lastOpt+2]);//fasta file

    cout<<lastOpt<<endl;
    cout<<fastaFile<<endl;
    cout<<bamfiletopen<<endl;
    cout<<region<<endl;
   
    for(int i=(lastOpt+3);i<(argc);i++){ //all but the last 3 args
	//cout<<"cont "<<string(argv[i])<<endl;
	filesAlFreq.push_back(string(argv[i]));	
    }




    ////////////////////////////////////
    // END Parsing arguments        //
    ////////////////////////////////////











    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    readIlluminaError(errFile,illuminaErrorsProb);

    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){
	    
    // 	    cout<<illuminaErrorsProb.s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }
    // return 1;
    
    ////////////////////////////////////////////////////////////////////////
    //
    // END  READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////











    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////
    readNucSubstitionFreq(deam5pfreqE,sub5p);
    readNucSubstitionFreq(deam3pfreqE,sub3p);
    readNucSubstitionFreq(deam5pfreqC,sub5pC);
    readNucSubstitionFreq(deam3pfreqC,sub3pC);



    // cout<<sub5p.size()<<endl;
    // cout<<sub3p.size()<<endl;

    // return 1;
    

    int defaultSubMatchIndex=0;

    for(int nuc1=0;nuc1<4;nuc1++){
    	for(int nuc2=0;nuc2<4;nuc2++){	    
    	    if(nuc1==nuc2)
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 1.0;
	    else
		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    	}
    	
    }


    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){	    
    // 	    cout<<sub5p[0].s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }

    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){	    
    // 	    cout<<sub3p[0].s[nuc2+nuc1*4]<<"\t";
    // 	}
    // 	cout<<endl;
    // }

    // return 1;
    
    ////////////////////////////////////////////////////////////////////////
    //
    // END  DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////
















    cout<<"test\t"<<bamfiletopen<<"\t"<<region<<endl;

     
    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
        cerr << "Could not open input BAM files " <<bamfiletopen<< endl;
        exit(1);
    }

    if ( !reader.OpenIndex(bamfiletopen+".bai") ) {
        cerr << "Could not open input index for BAM files " <<bamfiletopen+".bai"<< endl;
        exit(1);
    }

    const RefVector  references = reader.GetReferenceData();

    vector<string> tokens  = allTokens(region,':');
    if(tokens.size() != 2){
	cerr << "Invalid format for region:  " <<region<< endl;
	 exit(1);
    }

    vector<string> tokens2 = allTokens(tokens[1],'-');
    if(tokens2.size() != 2){
	cerr << "Invalid format for region:  " <<region<< endl;
	exit(1);
    }



    vector<MistarParser * > vectorOfMP;
    for(int i=(lastOpt+3);i<(argc);i++){ 
	// if(i==1 && string(argv[i]) == "-f"){
	//     force=true;
	//     continue;
	// }
	//cout<<"MP "<<string(argv[i])<<endl;
	if(!isFile(string(argv[i])+".tbi")){
	    cerr<<"Error: The allele count file must be tabix indexed: "<<string(argv[i])<<".tbi file not found"<<endl;
	    return 1;	
	}

	MistarParser * mp = new MistarParser(string(argv[i]),(string(argv[i])+".tbi"),tokens[0], destringify<int>(tokens2[0]) ,  destringify<int>(tokens2[1]));
	vectorOfMP.push_back(mp);
    }    

    // for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    // 	vectorOfMP[i]->repositionIterator(tokens[0], destringify<int>(tokens2[0]) ,  destringify<int>(tokens2[1]));
    // }


    int id = reader.GetReferenceID(tokens[0]);

    BamRegion regionbam (id, destringify<int>(tokens2[0]) , id, destringify<int>(tokens2[1]));
    
    reader.SetRegion(regionbam);


    Fasta fastaReference;
    if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	exit(1);
    }
	



    MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,&vectorOfMP,destringify<unsigned int>(tokens2[0]),  ignoreMQ);
    PileupEngine pileup;
    pileup.AddVisitor(cv);


    BamAlignment al;
    unsigned int numReads=0;
    cerr<<"Reading BAM file ..."<<endl;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<"name:\t"<<al.Name<<endl;
	numReads++;
	if(numReads !=0 && (numReads%100000)==0){
	    cerr<<"Read "<<thousandSeparator(numReads)<<" reads"<<endl;
	}

	if(al.IsMapped() && 
	   !al.IsFailedQC()){
	    // cerr<<al.Name<<endl;
	    pileup.AddAlignment(al);
	}
	    
    }
    cerr<<"...  done"<<endl;
    
    
    //clean up
    pileup.Flush();

    reader.Close();
    delete cv;


    return 0;
}

