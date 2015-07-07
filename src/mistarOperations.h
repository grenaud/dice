/*
 * mistarOperations
 * Date: Apr-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef mistarOperations_h
#define mistarOperations_h

#include <vector>
/* #include <mlpack/core.hpp> */

#include "MistarParser.h"
/* #include "GenomicRange.h" */

/* using namespace arma; */
using namespace std;

typedef struct{
    string name;
    uint64_t startIndexChr;
    uint64_t endIndexChr;
    uint64_t length;
} chrinfo;

//reading fasta index
/* void readFastaIndex(const string fastaIndex, */
/* 		    vector<chrinfo> & chrFound, */
/* 		    uint64_t & genomeLength); */

void initFiles(vector<MistarParser * > & vectorOfMP,
	       //bool & atLeastOneHasData,
	       vector<bool> & hasData,
	       vector<int> & popSizePerFile,
	       vector<AlleleRecords *> & vecAlleleRecords,
	       string & chr1,
	       unsigned int & coordCurrent,
	       bool printOnlyFirstPop=false); //if we print only the populations of the first file

bool sanityCheck(vector<MistarParser * > & vectorOfMP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 string & chr1,
		 unsigned int & coordCurrent,
		 string & chrcheck ,
		 char & refAllele,
		 bool force=false);

bool printAllele(vector<MistarParser * > & vectorOfMP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 string & chr1,
		 unsigned int & coordCurrent,
		 bool force=false);


bool hasCpG(vector<MistarParser * > & vectorOfMP,
	    vector<bool> & hasData,
	    vector<bool> & hasCoordinate,
	    vector<int> & popSizePerFile,
	    vector<AlleleRecords *> & vecAlleleRecords);

bool populateFreqVec(vector<MistarParser * > & vectorOfMP,
		     vector<bool> & hasData,
		     vector<bool> & hasCoordinate,
		     vector<int> & popSizePerFile,
		     vector<AlleleRecords *> & vecAlleleRecords,
		     string & chr1,
		     unsigned int & coordCurrent,
		     vector<double> * freqVec,
		     char & derAllele,
		     char & ancAllele,
		     bool force);

/* map< string, vector<GenomicRange> * > * readBEDSortedfile(string filetoread); */

#endif
