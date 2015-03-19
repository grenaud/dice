/*
 * MistarParser
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef MistarParser_h
#define MistarParser_h

#include <gzstream.h>
#include <sys/time.h> //for srand

#include "utils.h"
#include "SingleAllele.h"
#include "AlleleRecords.h"
#include "ReadTabix.h"

using namespace std;



/* typedef struct SingleAllele { */
/*     int refCount; */
/*     int altCount; */
/*     bool isCpg; */
/* } SingleAllele; */


/* typedef struct AlleleRecords { */
/*     string chr; */
/*     unsigned int coordinate; */
/*     char ref; */
/*     char alt; */
/*     vector<SingleAllele> * vectorAlleles; */
/* } AlleleRecords; */

/* inline string singleAlleleAsString(const SingleAllele & toprint){ */
/*     string toReturn =""+ */
/* 	stringify( toprint.refCount)+","+ */
/* 	stringify( toprint.altCount)+":"+ */
/* 	stringify( toprint.isCpg); */
/*     return toReturn; */
/* } */

/* inline string alleleRecordsAsString(const AlleleRecords & toprint){ */
/*     string toReturn =""+ */
/* 	toprint.chr+"\t"+ */
/* 	stringify(toprint.coordinate)+"\t"+ */
/* 	toprint.ref+","+ */
/* 	toprint.alt+"\t"+ */
/* 	vectorToString(*(toprint.vectorAlleles),"\t"); */
/* 	/\* for(unsigned int i=0;i<toprint.vectorAlleles->size();i++){        *\/ */
/* 	/\*     /\\* toReturn+= *\\/ *\/ */
/* 	/\*     /\\* 	toprint.vectorAlleles->at(i).toString(); *\\/ *\/ */
/* 	/\*     /\\* if( i!= (toprint.vectorAlleles->size()-1) ) *\\/ *\/ */
/* 	/\*     /\\* 	toReturn+="\t"; *\\/ *\/ */
/* 	/\* } *\/ */
/*     return toReturn; */
/* } */


class MistarParser{
 private:
    /* bool zipped; */
    /* istream * mystream; */
    vector<string> * populationNames;
    unsigned int numberPopulations;
    igzstream   * myFilezipped;
    /* ifstream    * myFile; */
    AlleleRecords * allRecToReturn;
    int numberOfTimesHasDataWasCalled;
    string header;
    string headerNoDefline;
    const vector<string> * dataToRead;
    unsigned int dataToReadInd;

    string defline;
    string currentline;

    int numbernew;
    int numberdel;

    ReadTabix * rt;
    void parseHeader(istream & in);
    bool getNextLine();

    bool stringMode;
    bool tabixMode;
    bool textMode;
 public:
    MistarParser(string filename);
    MistarParser(string file,string indexForFile,string chrName,int start,int end);
    MistarParser(const vector<string> * dataToRead,const vector<string> & populationNames_);
   
    MistarParser(const MistarParser&); // not implemented
    ~MistarParser();
    bool hasData();
    AlleleRecords  * getData();

    string getHeader(string prefix="");
    string getHeaderNoDefline(string prefix="");

    string getDefline();
    void repositionIterator(string chrName,int start,int end);

    const vector<string> *   getPopulationsNames() const ;
};




inline char sampleRandomRefAltAllele(char ref,char alt,int refCount,int altCount);

inline char sampleRandomRefAltAllele(char ref,char alt,int refCount,int altCount){

    if(refCount == 0 && altCount == 0 ){
	cerr<<"MistarParser::sampleRandomRefAltAllele() cannot sample when both allele counts are 0, find the null record in your input, exiting"<<endl;
	exit(1);
    }

    if(refCount == 0  ){//no reference, return alternative
	return alt;
    }

    if(altCount == 0  ){//no reference, return alternative
	return ref;
    }

    if(refCount ==  altCount ){//homozygous (kind of)
	if(randomBool())
	    return ref;
	else
	    return alt;
    }

    int totalCount = refCount+  altCount;
    int randnum=(callRand()%totalCount)+1;
    /* cout<<refCount<<" "<<altCount<<" randnum "<<randnum<<endl; */

    if(randnum<=refCount){
	/* cout<<"REF"<<endl; */
	return ref;
    }else{
	/* cout<<"ALT"<<endl; */
	return alt;
    }

    /* cerr<<"MistarParser::sampleRandomRefAltAllele() tell Gabriel to implement when the allele count is uneven, it should be sampled according to the count, exiting"<<endl; */
    /* exit(1); */
    
}

#endif
