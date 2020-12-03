#if defined(__CYGWIN__) 
#define atanl(X) atan(X)
#define logl(X) log(X)
#define sqrtl(X) sqrt(X)
#endif
  

 // #define DEBUGMST
// #define DEBUGERRORP
// #define DEBUGFREQ


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

#include "libgab.h"
#include "miscfunc.h"
#include "libnuc.h"

#include "ReconsReferenceBAM.h"
#include "MistarParser.h"
#include "mistarOperations.h"

using namespace BamTools;
using namespace std;
const long double PI  = atanl(1.0L)*4;   







char   offsetQual=33;

long double probMatch[MAXMAPPINGQUAL];
long double probMismatch[MAXMAPPINGQUAL];
// long double likeMatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
// long double likeMismatchMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

// long double likeMatchProb[MAXMAPPINGQUAL];
// long double likeMismatchProb[MAXMAPPINGQUAL];
// long double likeMatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];
// long double likeMismatchProbMQ[MAXMAPPINGQUAL][MAXMAPPINGQUAL];

long double probMapping[MAXMAPPINGQUAL];
long double probMismapping[MAXMAPPINGQUAL];



// probSubstition illuminaErrorsProb;
// vector<probSubstition> sub5p;
// vector<probSubstition> sub3p;
// vector<probSubstition> sub5pC;
// vector<probSubstition> sub3pC;

vector<substitutionRates> sub5p;
vector<substitutionRates> sub3p;
vector<substitutionRates> sub5pC;
vector<substitutionRates> sub3pC;

// probSubstition defaultSubMatch;
substitutionRates defaultSubMatch;

string dnaAlphabet="ACGT";
map<int, PHREDgeno> pos2phredgeno;

// map<string, long double > read2endoProb; //map seq id to probability that the read is endogenous using a deamination model
// long double read2endoProbInit=false;


// long double probCorrect(char b,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){
//     long double dnaProb[4];
//     long double probError  =0.0;
//     long double probCorrect=0.0;

//     for(int i=0;i<4;i++){
// 	int dinucIndex;
// 	char bModel=dnaAlphabet[i];
// 	if( isRev ){
// 	    dinucIndex =      baseResolved2int(complement(b)) *4+baseResolved2int(complement(bModel));
// 	}else{
// 	    dinucIndex =      baseResolved2int(           b)  *4+baseResolved2int(           bModel);
// 	}
// 	//cout<<"model ="<<bModel<<" "<<b<<" "<<dinucIndex<<" "<<probSubMatch->s[dinucIndex] <<endl;
// 	dnaProb[i]  = probMismapping[mpq]*0.25;
// 	if(bModel == b){
// 	    dnaProb[i] += probMapping[mpq]*( (probMatch[bq] * probSubMatch->s[dinucIndex]                          ) );
// 	    probCorrect+=dnaProb[i];
// 	}else{
// 	    dnaProb[i] += probMapping[mpq]*( (probMatch[bq] * probSubMatch->s[dinucIndex]  + probMismatch[bq]/3.0  ) );
// 	    probError+=dnaProb[i];
// 	}

// 	//cout<<"dnaProb ["<<i<<"]= "<<dnaProb[i]<<endl;
//     }
//     return (probCorrect/(probCorrect+probError));
//     //return probCorrect;
    
// }


// pair<long double,long double> probCorrectcubed(char b, char bModel,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){
//     // mpq = min(mpq,20);

//     int dinucIndex;
//     if( isRev ){
// 	dinucIndex =      baseResolved2int(complement(bModel)) *4+baseResolved2int(complement(b));
//     }else{
// 	dinucIndex =      baseResolved2int(           bModel)  *4+baseResolved2int(           b);
//     }


//     long double probSubDeam              = probSubMatch->s[dinucIndex];
//     long double probSameDeam             = 1.0-probSubDeam;
   
//     long double probCorrectAll3          = powl(probMatch[bq],3.0);
//     long double probIncorrectAll3        = 1.0-probCorrectAll3;

    
//     long double probCorrect2Params       = probCorrectAll3 * probSameDeam + probIncorrectAll3 * 0.5;
//     long double probIncorrect2Params     = probCorrectAll3 * probSubDeam  + probIncorrectAll3 * 0.5;

//     long double probCorrect2ParamsMQ     = probMapping[mpq]*  probCorrect2Params + probMismapping[mpq]*0.5;
//     long double probIncorrect2ParamsMQ   = probMapping[mpq]*probIncorrect2Params + probMismapping[mpq]*0.5;
    
//     pair<long double,long double>  toReturn;
//     toReturn.first  = probCorrect2ParamsMQ;
//     toReturn.second = probIncorrect2ParamsMQ;
//     return toReturn;
// }

//pair<long double,long double> probCorrect(char b, char bModel,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){

pair<long double,long double> probCorrect(char b, char bModel,bool isRev,int mpq,int bq, const substitutionRates * probSubMatch){
    //    mpq = min(mpq,20);

    int dinucIndex;
    if( isRev ){
	//dinucIndex =     baseResolved2int(complement(bModel)) *3+baseResolved2int(complement(b));
	dinucIndex =	 dimer2index(complement(bModel),complement(b));
    }else{
	//dinucIndex =      baseResolved2int(           bModel)  *3+baseResolved2int(           b);
	dinucIndex =	 dimer2index(           bModel ,           b);
    }

#ifdef DEBUGERRORP
    cout<<"probCorrect bModel\t"<<bModel<<"\tb\t"<<b<<"\t"<<dinucIndex<<"\t"<<isRev<<endl;
#endif

    // if(dinucIndex ==  7)
    //    mpq = min(mpq,10);
    long double probSubDeam              = probSubMatch->s[dinucIndex];
    long double probSameDeam             = 1.0-probSubDeam;
    //long double probSameDeam           = probSubMatch->s[dinucIndex];

    long double probCorrect              =     probMatch[bq];
    long double probIncorrect            = 1.0-probMatch[bq];

    
    long double probCorrect2Params       = probCorrect * probSameDeam + probIncorrect * 0.5;
    long double probIncorrect2Params     = probCorrect * probSubDeam  + probIncorrect * 0.5;

    long double probCorrect2ParamsMQ     = probMapping[mpq]*  probCorrect2Params + probMismapping[mpq]*0.5;
    long double probIncorrect2ParamsMQ   = probMapping[mpq]*probIncorrect2Params + probMismapping[mpq]*0.5;
    
    pair<long double,long double>  toReturn;
    toReturn.first  = probCorrect2ParamsMQ;
    toReturn.second = probIncorrect2ParamsMQ;


    //TO TEST, to remove
    // toReturn.first  = 1.0;
    // toReturn.second = 0.0;
    return toReturn;
}

// pair<long double,long double> probCorrect3(char b, char bModel,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){
//      // mpq = min(mpq,30);
//      // mpq = max(mpq,60);

//     int dinucIndex;
//     if( isRev ){
// 	dinucIndex =      baseResolved2int(complement(bModel)) *4+baseResolved2int(complement(b));
//     }else{
// 	dinucIndex =      baseResolved2int(           bModel)  *4+baseResolved2int(           b);
//     }


//     long double probSubDeam  = probSubMatch->s[dinucIndex];
//     long double probSameDeam = 1.0-probSubDeam;
   
//     long double probIncorrectAll3  = sqrtl(probMismatch[bq]);
//     long double probCorrectAll3    = 1.0-probIncorrectAll3;

    
//     long double probCorrect2Params   = probCorrectAll3 * probSameDeam + probIncorrectAll3 * 0.5;
//     long double probIncorrect2Params = probCorrectAll3 * probSubDeam  + probIncorrectAll3 * 0.5;

//     long double probCorrect2ParamsMQ     = probMapping[mpq]*probCorrectAll3   + probMismapping[mpq]*0.5;
//     long double probIncorrect2ParamsMQ   = probMapping[mpq]*probIncorrectAll3 + probMismapping[mpq]*0.5;
    
//     pair<long double,long double>  toReturn;
//     toReturn.first  = probCorrect2ParamsMQ;
//     toReturn.second = probIncorrect2ParamsMQ;
//     return toReturn;
// }


// pair<long double,long double> probCorrect4(char b, char bModel,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){
//      // mpq = min(mpq,30);
//      // mpq = max(mpq,60);

//     int dinucIndex;
//     if( isRev ){
// 	dinucIndex =      baseResolved2int(complement(bModel)) *4+baseResolved2int(complement(b));
//     }else{
// 	dinucIndex =      baseResolved2int(           bModel)  *4+baseResolved2int(           b);
//     }


//     long double probSubDeam  = probSubMatch->s[dinucIndex];
//     long double probSameDeam = 1.0-probSubDeam;
   
//     long double probIncorrectAll3  = sqrtl(probMismatch[bq]);
//     long double probCorrectAll3    = 1.0-probIncorrectAll3;

    
//     long double probCorrect2Params   = probCorrectAll3 * probSameDeam + probIncorrectAll3 * 0.5;
//     long double probIncorrect2Params = probCorrectAll3 * probSubDeam  + probIncorrectAll3 * 0.5;

//     long double probCorrect2ParamsMQ     = probMapping[mpq]*  probCorrect2Params + probMismapping[mpq]*0.5;
//     long double probIncorrect2ParamsMQ   = probMapping[mpq]*probIncorrect2Params + probMismapping[mpq]*0.5;
    
//     pair<long double,long double>  toReturn;
//     toReturn.first  = probCorrect2ParamsMQ;
//     toReturn.second = probIncorrect2ParamsMQ;
//     return toReturn;
// }



//pair<long double,long double> probCorrect(char b, char bModel,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){

    // long double dnaProb[4];
    // long double probError  =0.0;
    // long double probCorrect=0.0;



    // for(int i=0;i<4;i++){
    // 	int dinucIndex;
    // 	char bModel=dnaAlphabet[i];
    // 	if( isRev ){
    // 	    dinucIndex =      baseResolved2int(complement(b)) *4+baseResolved2int(complement(bModel));
    // 	}else{
    // 	    dinucIndex =      baseResolved2int(           b)  *4+baseResolved2int(           bModel);
    // 	}
    // 	//cout<<"model ="<<bModel<<" "<<b<<" "<<dinucIndex<<" "<<probSubMatch->s[dinucIndex] <<endl;
    // 	dnaProb[i]  = probMismapping[mpq]*0.25;
    // 	if(bModel == b){
    // 	    dnaProb[i] += probMapping[mpq]*( (probMatch[bq] * probSubMatch->s[dinucIndex]                          ) );
    // 	    probCorrect+=dnaProb[i];
    // 	}else{
    // 	    dnaProb[i] += probMapping[mpq]*( (probMatch[bq] * probSubMatch->s[dinucIndex]  + probMismatch[bq]/3.0  ) );
    // 	    probError+=dnaProb[i];
    // 	}

    // 	//cout<<"dnaProb ["<<i<<"]= "<<dnaProb[i]<<endl;
    // }
    // return (probCorrect/(probCorrect+probError));
    //return probCorrect;
    
//}


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
		  vector<MistarParser * > & vectorOfMP,
		  unsigned int coordFirst,
		  unsigned int coordLast,
		  const bool ignoreMQ,
		  vector<singleSite> * dataSitesVec,
		  bool useCpG,
		  // const long double contaminationPrior,
		  //const bool singleCont
		  bool * wasDataFound,
		  const int anchIDX,
		  const int admxIDX
		  )
    : PileupVisitor()
      // , m_references(references)
    , m_fastaReference(fastaReference)
    , m_vectorOfMP(vectorOfMP)
    , m_coordFirst(coordFirst)
    , m_coordLast(coordLast)
    , ignoreMQ(ignoreMQ)
    , m_dataSitesVec(dataSitesVec)
    , m_useCpG(useCpG)
    , m_anchIDX(anchIDX)
    , m_admxIDX(admxIDX)
    // , contaminationPrior(contaminationPrior)
    // , singleCont(singleCont)
  { 
     
      //      cout<<"constr size="<<vectorOfMP.size()<<endl;
      *wasDataFound=initFiles(vectorOfMP,
			      // atLeastOneHasData,
			      hasData,
			      popSizePerFile,
			      vecAlleleRecords,
			      chr1,
			      coordFirst,
			      false);
      
      //cout<<"constr"<<endl;
    hasCoordinate = vector<bool>(m_vectorOfMP.size(),true);//dummy value


  }
  ~MyPileupVisitor(void) { }
  
  // PileupVisitor interface implementation
public:
    
  void Visit(const PileupPosition& pileupData) {
      // for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){				
      // 	  cout<<"visit "<<i<<"\t"<<pileupData.PileupAlignments.size()<<endl;
      // }
      
    char referenceBase = 'N';
	    
    ///unsigned int posAlign = pileupData.Position+1;
    //    int posVector=int(pileupData.Position)%sizeGenome;
	    

    unsigned int posAlign = pileupData.Position+1;
    // int posVector=int(pileupData.Position)%sizeGenome;
	    
    // if( (posAlign%100) == 0){
    // 	cerr<<"pos  = "<<posAlign<<endl;
    // }
    // cout<<endl<<"pos = "<<posAlign<<endl;
    if(posAlign>m_coordLast)
	return ;

    if ( !m_fastaReference->GetBase(pileupData.RefId, posAlign-1, referenceBase ) ) {
	cerr << "bamtools convert ERROR: pileup conversion - could not read reference base from FASTA file at chr "<<pileupData.RefId<<" position "<<(posAlign-1) << endl;
	exit(1);
    }


    // BEGIN READING ALLELE FREQ

    bool stayLoop=true;
    bool samePos=false;
    //the following syncs the allele count files to make sure they are at the same positon

    while(stayLoop){

#ifdef DEBUGMST
	cout<<"posAlign "<<posAlign<<"\t"<<pileupData.PileupAlignments.size()<<endl;
#endif

	bool allHaveCoordinate=true;

	for(unsigned int i=0;i<m_vectorOfMP.size();i++){ 
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

	    samePos=true;
	    // printAllele(*m_vectorOfMP,
	    // 		hasData,
	    // 		hasCoordinate,
	    // 		popSizePerFile,
	    // 		vecAlleleRecords,
	    // 		chr1,
	    // 		posAlign,
	    // 		false);


	    // 	seekdata:
	    allHaveCoordinate=false;
	    stayLoop=false;
	    break;

#ifdef DEBUGFREQ
	    if(0)
	    for(unsigned int i=0;i<m_vectorOfMP.size();i++){ 
		 
		if(!hasData[i] ){
		    cerr<<"Invalid state"<<endl;
		    exit(1);
		}
		hasData[i]  =  m_vectorOfMP.at(i)->hasData();
		if(hasData[i]){
		    vecAlleleRecords[i] = m_vectorOfMP.at(i)->getData() ;
		}else{
		    stayLoop=false;
		    break;
		}
		
	    }
#endif

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
	    
	    
	    for(unsigned int i=0;i<m_vectorOfMP.size();i++){ 
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
		    hasData[i]  =   m_vectorOfMP[i]->hasData();
		    if(hasData[i]){
			vecAlleleRecords[i] = m_vectorOfMP.at(i)->getData() ;
		    }else{
			stayLoop=false;
			break;
		    }
		}
		
	    }

	    
	}//end different coord

	
    }//end stay loop
    // END READING ALLELE FREQ


    if(samePos){

	char derAllele='N';
	char ancAllele='N';
	vector<double> derFreq;

	if(    (pileupData.PileupAlignments.size()>0) ){
	    
	    bool populatedFreq= populateFreqVec(m_vectorOfMP,
						hasData,
						hasCoordinate,
						popSizePerFile,
						vecAlleleRecords,
						chr1,
						posAlign,
						&derFreq,
						derAllele,
						ancAllele,
						false);

	    
	    if(!populatedFreq){
		return ;
	    }

	    if(!m_useCpG){//if we care about CpGs
		if(hasCpG(m_vectorOfMP, //if one of the records has it
			  hasData,
			  hasCoordinate,
			  popSizePerFile,
			  vecAlleleRecords)){
		    return ;
		}
	    }


	    if(populatedFreq){

		
#ifdef DEBUGMST
		if(m_vectorOfMP.size()>1){		    
		    printAllele(m_vectorOfMP,
				hasData,
				hasCoordinate,
				popSizePerFile,
				vecAlleleRecords,
				chr1,
				posAlign,
				false);
		}else{
		    cout<<"alelle record "<<( *(vecAlleleRecords[0]) )<<endl;		    			
		}		
#endif

		if(vecAlleleRecords[0]->ref !=  referenceBase){
		    cerr<<"Error: the reference allele differs at  "<<chr1<<":"<<posAlign<<" allele count says : "<<vecAlleleRecords[0]->ref<<" bam file says "<<referenceBase<<endl;
		    exit(1);
		}

	
	    }
	}else{ //if has no alignments
	    return ;
	}

#ifdef DEBUGFREQ
	if(0){//testing

	    if(derFreq[0] == 0.0){//skip sites with fixed bases
		for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){				


		
		    if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
			pileupData.PileupAlignments[i].IsNextInsertion    &&
			(pileupData.PileupAlignments[i].InsertionLength    >0)){ //has insertion		    
			continue;//ignore
		    }

		    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
			pileupData.PileupAlignments[i].IsNextInsertion   &&
			(pileupData.PileupAlignments[i].InsertionLength == 0)){ //has deletion
			continue;
		    }

		    //base that was read
		    char b   = pileupData.PileupAlignments[i].Alignment.QueryBases[pileupData.PileupAlignments[i].PositionInAlignment];
		    //quality score
		    char q   = pileupData.PileupAlignments[i].Alignment.Qualities[pileupData.PileupAlignments[i].PositionInAlignment]-offsetQual;
		    int  m   = int(pileupData.PileupAlignments[i].Alignment.MapQuality);

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

		    //cout<<chr1<<":"<<posAlign<<"\t"<<derAllele<<"\t"<<ancAllele<<"\t"<<b<<"\t"<<int(q)<<"\t"<<m<<"\t"<<dist5p<<"\t"<<dist3p<<"\t"<<pileupData.PileupAlignments[i].Alignment.Name<<"\t"<<(b==ancAllele)<<endl;
	
		}

	    }
	}
#endif

	if( derAllele=='N' ||
	    ancAllele=='N' )
	    return ;


	if(m_admxIDX == -1){//2pop
	    if( derFreq[m_anchIDX] <= 0 || derFreq[m_anchIDX] >= 1.0)//skip sites with fixed bases for the anchor
	     	return ;
	}else{//3pop
	    double derFreqSum = ( (derFreq[m_anchIDX]+derFreq[m_admxIDX])/2.0 );
	    
	    if( derFreqSum <= 0 || derFreqSum >= 1.0)//skip sites with fixed bases for the combination of the anchor and admixed
	     	return ;

	}
	
	// for(unsigned int df=0;df<derFreq.size();df++){	    
	//     if(  m_anchIDX !=  int(df) )  //this index is the anchor, forego the polymorphic requirement
	// 	continue;

	//     if(derFreq[df] <= 0 || derFreq[df] >= 1.0)//skip sites with fixed bases
	//     	return ;
	// }

#ifdef DEBUGERRORP
	cout<<"samepos "<<chr1<<":"<<posAlign<<" derFreq "<<vectorToString(derFreq)<<"\td="<<derAllele<<"\ta="<<ancAllele<<endl;
	cout<<"derFreq: "<<vectorToString(derFreq,"\t")<<endl;
		// for(unsigned int k=0;k<singleBaseVec.size();k++){
		//     cout<<"\t"<<singleBaseVec[k].b<<"\t"<<singleBaseVec[k].bq;
		// }
	cout<<endl<<endl;
#endif
	//store 
	singleSite toAddSS;
	toAddSS.freqDerived = derFreq;
	toAddSS.ancAllele   = ancAllele;
	toAddSS.derAllele   = derAllele;
	toAddSS.chr         = chr1;
	toAddSS.coord       = posAlign;


	//	vector<singleBase> singleBaseVec;

	//re-iterate over each read to compute the likelihood for each insert (or pair of inserts)
	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){				


		
	    if( !pileupData.PileupAlignments[i].IsCurrentDeletion &&
		pileupData.PileupAlignments[i].IsNextInsertion    &&
		(pileupData.PileupAlignments[i].InsertionLength    >0)){ //has insertion		    
		continue;//ignore
	    }

	    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
		pileupData.PileupAlignments[i].IsNextInsertion   &&
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

	    if(b != derAllele &&
	       b != ancAllele) //skip tri-allelic sites	       
		return;
	
	    

	    // vector<singleBase> singleBaseVec;
	    singleBase sb;
	    sb.b    = b;
	    sb.bq   = int(q);
	    sb.mpq  = m;

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
		    		    
	    // probSubstition * probSubMatchToUseEndo = &defaultSubMatch ;
	    // probSubstition * probSubMatchToUseCont = &defaultSubMatch ;
	    substitutionRates * probSubMatchToUseEndo = &defaultSubMatch ;
	    substitutionRates * probSubMatchToUseCont = &defaultSubMatch ;


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

	    // for(int subi=0;subi<12;subi++){
	    // 	cout<<"subi"<<"\t"<<probSubMatchToUseEndo->s[subi]<<endl;
	    // }


	    //int dinucIndexModelObs;
	    // int dinucIndexDD;
	    // int dinucIndexDA;
	    // int dinucIndexAD;
	    // int dinucIndexAA;

	    //there are four possibilities:
	    // model       observation
	    // d           d
	    // d           a
	    // a           d
	    // a           a

	    // if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){

	    // 	dinucIndexDD =      baseResolved2int(complement(derAllele)) *4+baseResolved2int(complement(derAllele));
	    // 	dinucIndexDA =      baseResolved2int(complement(derAllele)) *4+baseResolved2int(complement(ancAllele));
	    // 	dinucIndexAD =      baseResolved2int(complement(ancAllele)) *4+baseResolved2int(complement(derAllele));
	    // 	dinucIndexAA =      baseResolved2int(complement(ancAllele)) *4+baseResolved2int(complement(ancAllele));

	    // }else{

	    // 	dinucIndexDD =      baseResolved2int(           derAllele)  *4+baseResolved2int(           derAllele);
	    // 	dinucIndexDA =      baseResolved2int(           derAllele)  *4+baseResolved2int(           ancAllele);
	    // 	dinucIndexAD =      baseResolved2int(           ancAllele)  *4+baseResolved2int(           derAllele);
	    // 	dinucIndexAA =      baseResolved2int(           ancAllele)  *4+baseResolved2int(           ancAllele);

	    // }

	    
	    // if(b== derAllele){//b is the observation and b is derived


	    // }else{

	    // 	if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){
	    // 	    dinucIndex =      baseResolved2int(complement(derAllele)) *4+baseResolved2int(complement(ancAllele));
	    // 	}else{
	    // 	    dinucIndex =      baseResolved2int(           derAllele)  *4+baseResolved2int(           ancAllele);
	    // 	}


	    // }
	    // END DEAMINATION COMPUTATION


	    // long double probDDeCorrect = 1.0-probSubMatchToUseEndo->s[dinucIndexDA];
	    // long double probDAeCorrect =     probSubMatchToUseEndo->s[dinucIndexDA];

	    // long double probADeCorrect =     probSubMatchToUseEndo->s[dinucIndexAD];
	    // long double probAAeCorrect = 1.0-probSubMatchToUseEndo->s[dinucIndexAD];

	    // long double probDDcCorrect = 1.0-probSubMatchToUseCont->s[dinucIndexDA];
	    // long double probDAcCorrect =     probSubMatchToUseCont->s[dinucIndexDA];

	    // long double probADcCorrect =     probSubMatchToUseCont->s[dinucIndexAD];
	    // long double probAAcCorrect = 1.0-probSubMatchToUseCont->s[dinucIndexAD];



	    // long double probDDe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDeCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;
	    // long double probDAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAeCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;

	    // long double probADe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADeCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;
	    // long double probAAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAeCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;

	    // long double probDDc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDcCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;
	    // long double probDAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAcCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;

	    // long double probADc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADcCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;
	    // long double probAAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAcCorrect) + (probMismatch[sb.bq])*0.5) ) + probMismapping[sb.mpq]*0.5;

	    // long double probDDe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDeCorrect) + 0.0                      ) ) + probMismapping[sb.mpq]*0.5;
	    // long double probDAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAeCorrect) + (probMismatch[sb.bq])    ) ) + probMismapping[sb.mpq]*0.5;

	    // long double probADe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADeCorrect) + (probMismatch[sb.bq])    ) ) + probMismapping[sb.mpq]*0.5;
	    // long double probAAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAeCorrect) + 0.0                      ) ) + probMismapping[sb.mpq]*0.5;

	    // long double probDDc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDcCorrect) + 0.0                      ) ) + probMismapping[sb.mpq]*0.5;
	    // long double probDAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAcCorrect) + (probMismatch[sb.bq])    ) ) + probMismapping[sb.mpq]*0.5;

	    // long double probADc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADcCorrect) + (probMismatch[sb.bq])    ) ) + probMismapping[sb.mpq]*0.5;
	    // long double probAAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAcCorrect) + 0.0                      ) ) + probMismapping[sb.mpq]*0.5;

	    //test 4
	    // long double probDDe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDeCorrect) + 0.0                      ) ) ;
	    // long double probDAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAeCorrect) + (probMismatch[sb.bq])    ) ) ;

	    // long double probADe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADeCorrect) + (probMismatch[sb.bq])    ) ) ;
	    // long double probAAe = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAeCorrect) + 0.0                      ) ) ;

	    // long double probDDc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDcCorrect) + 0.0                      ) ) ;
	    // long double probDAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAcCorrect) + (probMismatch[sb.bq])    ) ) ;

	    // long double probADc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probADcCorrect) + (probMismatch[sb.bq])    ) ) ;
	    // long double probAAc = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probAAcCorrect) + 0.0                      ) ) ;



	    // long double probDDe = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseEndo->s[dinucIndexDD]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexDD]);
	    // long double probDAe = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseEndo->s[dinucIndexDA]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexDA]);
	    // long double probADe = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseEndo->s[dinucIndexAD]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexAD]);
	    // long double probAAe = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseEndo->s[dinucIndexAA]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexAA]);

	    // long double probDDc = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseCont->s[dinucIndexDD]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexDD]);
	    // long double probDAc = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseCont->s[dinucIndexDA]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexDA]);
	    // long double probADc = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseCont->s[dinucIndexAD]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexAD]);
	    // long double probAAc = (likeMatchProb[sb.bq] * (0.0+probSubMatchToUseCont->s[dinucIndexAA]) + (1.0 - likeMatchProb[sb.bq])* illuminaErrorsProb.s[dinucIndexAA]);


	    
	    // long double probDDe =        probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);
	    // long double probDAe = 1.0 -  probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);

	    // long double probADe = 1.0 -  probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);
	    // long double probAAe =        probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);

	    // long double probDDc =        probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);
	    // long double probDAc = 1.0 -  probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);

	    // long double probADc = 1.0 -  probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);
	    // long double probAAc =        probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);


	    pair<long double,long double> probsub;
	    //                    obs       model
	    probsub = probCorrect(ancAllele,derAllele,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);


	    long double probDDe = probsub.first;
	    long double probDAe = probsub.second;
		
	    probsub = probCorrect(ancAllele,derAllele,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);

	    long double probDDc = probsub.first;
	    long double probDAc = probsub.second;

	    //                    obs       model
	    probsub = probCorrect(derAllele,ancAllele,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);

	    long double probAAe = probsub.first;
	    long double probADe = probsub.second;
	    
	    probsub = probCorrect(derAllele,ancAllele,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseCont);

	    long double probAAc = probsub.first;
	    long double probADc = probsub.second;



	    

#ifdef DEBUGERRORP	    

	    cout<<"qual "<<sb.bq<<" mq "<< sb.mpq<<endl;
	    int dinucIndex;
	    if( pileupData.PileupAlignments[i].Alignment.IsReverseStrand() ){
		//dinucIndex =      baseResolved2int(complement(derAllele)) *4+baseResolved2int(complement(ancAllele));
		dinucIndex =	  dimer2index(     complement(derAllele),  complement(ancAllele));
	    }else{
		//dinucIndex =      baseResolved2int(           derAllele)  *4+baseResolved2int(           ancAllele);
		dinucIndex =	  dimer2index(                derAllele,              ancAllele);
	    }

	    cout<<"D="<<derAllele<<",A="<<ancAllele<<"\t"<<probSubMatchToUseEndo->s[dinucIndex]<<endl;
	    for(int i=0;i<12;i++){
		cout<<i<<"="<<probSubMatchToUseEndo->s[i]<<endl;
	    }

	    // long double correct= probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDDeCorrect)                              ) ) + probMismapping[sb.mpq]*0.25;
	    // long double wrong  = probMapping[sb.mpq]*( (probMatch[sb.bq] * (probDAeCorrect) + (probMismatch[sb.bq])/3.0  ) ) + probMismapping[sb.mpq]*0.25;
	    //cout<<"test "<<correct/(correct+wrong*3.0)<<"\t"<<(1- (correct/(correct+wrong*3.0)))<<endl;
	    //long double probCorrect(char b,bool isRev,int mpq,int bq, const probSubstition * probSubMatch){
	    //probCorrect(b,pileupData.PileupAlignments[i].Alignment.IsReverseStrand(),sb.mpq,sb.bq,probSubMatchToUseEndo);


	    // cout<<"dde\t"<<probDDeCorrect<<endl;	    
	    // cout<<"dae\t"<<probDAeCorrect<<endl;	    
	    // cout<<"ade\t"<<probADeCorrect<<endl;	    
	    // cout<<"aae\t"<<probAAeCorrect<<endl;	    


	    cout<<"DDe"<<"\tdde\t"<<probDDe<<endl;	    
	    cout<<"DAe"<<"\tdae\t"<<probDAe<<endl;	    
	    cout<<"ADe"<<"\tade\t"<<probADe<<endl;	    
	    cout<<"AAe"<<"\taae\t"<<probAAe<<endl;	    

	    cout<<"DDc"<<"\tddc\t"<<probDDc<<endl;	    
	    cout<<"DAc"<<"\tdac\t"<<probDAc<<endl;	    
	    cout<<"ADc"<<"\tadc\t"<<probADc<<endl;	    
	    cout<<"AAc"<<"\taac\t"<<probAAc<<endl;	    	    

	    // cout<<dinucIndexDD<<"\tdde\t"<<probDDe<<endl;	    
	    // cout<<dinucIndexDA<<"\tdae\t"<<probDAe<<endl;	    
	    // cout<<dinucIndexAD<<"\tade\t"<<probADe<<endl;	    
	    // cout<<dinucIndexAA<<"\taae\t"<<probAAe<<endl;	    

	    // cout<<dinucIndexDD<<"\tddc\t"<<probDDc<<endl;	    
	    // cout<<dinucIndexDA<<"\tdac\t"<<probDAc<<endl;	    
	    // cout<<dinucIndexAD<<"\tadc\t"<<probADc<<endl;	    
	    // cout<<dinucIndexAA<<"\taac\t"<<probAAc<<endl;	    	    

#endif

	    // long double eest = 0.005;
	    // probDDe = 1-eest;
	    // probDAe = eest;
	    // probADe = eest;
	    // probAAe = 1-eest;

	    // probDDc = 1-eest;
	    // probDAc = eest;
	    // probADc = eest;
	    // probAAc = 1-eest;

	    //test 5
	    // probDAe=1.0-probDDe;
	    // probADe=1.0-probAAe;
	    // probDAc=1.0-probDDc;
	    // probADc=1.0-probAAc;


#ifdef DEBUGERRORP	    
	    cout<<"i="<<i<<"\t"<<posAlign<<"\tref="<<referenceBase<<"\tb_obs="<<b<<"\tq="<<sb.bq<<"\t"<<pileupData.PileupAlignments[i].Alignment.Name<<"\t"<<"\td="<<derAllele<<"\ta="<<ancAllele <<"\t"<<dist5p<<"\t"<<dist3p<<"\t"<<(pileupData.PileupAlignments[i].Alignment.IsReverseStrand()?"R":"F")<<endl;

	    // cout<<dinucIndexDD<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexDD])<<"\t"<<probSubMatchToUseCont->s[dinucIndexDD]<<endl;	    
	    // cout<<dinucIndexDA<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexDA])<<"\t"<<probSubMatchToUseCont->s[dinucIndexDA]<<endl;	    
	    // cout<<dinucIndexAD<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexAD])<<"\t"<<probSubMatchToUseCont->s[dinucIndexAD]<<endl;	    
	    // cout<<dinucIndexAA<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexAA])<<"\t"<<probSubMatchToUseCont->s[dinucIndexAA]<<endl;	    

	    // cout<<dinucIndexDD<<"\t"<<(1.0-probSubMatchToUseEndo->s[dinucIndexDA])<<"\t"<<probSubMatchToUseCont->s[dinucIndexDD]<<endl;	    
	    // cout<<dinucIndexDA<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexDA])<<"\t"<<probSubMatchToUseCont->s[dinucIndexDA]<<endl;	    
	    // cout<<dinucIndexAD<<"\t"<<(0.0+probSubMatchToUseEndo->s[dinucIndexAD])<<"\t"<<probSubMatchToUseCont->s[dinucIndexAD]<<endl;	    
	    // cout<<dinucIndexAA<<"\t"<<(1.0-probSubMatchToUseEndo->s[dinucIndexAD])<<"\t"<<probSubMatchToUseCont->s[dinucIndexAA]<<endl;	    


	 


	    // cout<<dinucIndexDD<<"\tdde\t"<<probDDe<<endl;	    
	    // cout<<dinucIndexDA<<"\tdae\t"<<probDAe<<endl;	    
	    // cout<<dinucIndexAD<<"\tade\t"<<probADe<<endl;	    
	    // cout<<dinucIndexAA<<"\taae\t"<<probAAe<<endl;	    

	    // cout<<dinucIndexDD<<"\tddc\t"<<probDDc<<endl;	    
	    // cout<<dinucIndexDA<<"\tdac\t"<<probDAc<<endl;	    
	    // cout<<dinucIndexAD<<"\tadc\t"<<probADc<<endl;	    
	    // cout<<dinucIndexAA<<"\taac\t"<<probAAc<<endl;	    

#endif


	    sb.probDDc = probDDc;
	    sb.probADc = probADc;
	    sb.probDAc = probDAc;
	    sb.probAAc = probAAc;

	    sb.probDDe = probDDe;
	    sb.probADe = probADe;
	    sb.probDAe = probDAe;
	    sb.probAAe = probAAe;


	    //singleBaseVec.push_back(sb);
	    if(b ==  derAllele){ //derived observation
		toAddSS.sitesBAMd.push_back(sb);
	    }else{
		if(b ==  ancAllele){ //ancestral observation
		    toAddSS.sitesBAMa.push_back(sb);
		}else{
		    cerr<<"Internal error, base "<<b<<" is neither ancestral or derived"<<endl;
		    exit(1);
		}
	    }

	
	    // vector<bool> hasData;
	    // vector<AlleleRecords *> vecAlleleRecords;
	    // vector<bool>  hasCoordinate (m_vectorOfMP->size(),true);//dummy value





	

	
	}//end for each read
	
	m_dataSitesVec->push_back(toAddSS);
	//put single base Vec somwhere
	
	// if(samePos && 
	// (pileupData.PileupAlignments.size()>0)){
	// 	cout<<samePos<<endl;

    }//samePos

    // cout<<"end visit"<<endl;



    
  }//end visit()
        
private:

    vector<bool> hasData;
    vector<int> popSizePerFile;
    vector<AlleleRecords *> vecAlleleRecords;
    string chr1;
    unsigned int coordCurrent;
    vector<bool>  hasCoordinate;


    unsigned int m_coordFirst;
    unsigned int m_coordLast;

    RefVector m_references;
    Fasta * m_fastaReference;
    vector<MistarParser * >  m_vectorOfMP;
    vector<singleSite> *     m_dataSitesVec;
    // vector<singlePosInfo> * m_infoPPos;
    // int sizeGenome;
    bool ignoreMQ;
    // long double contaminationPrior;
    // bool singleCont;
    bool m_useCpG;    
    //        ostream*  m_out;
    int m_anchIDX;
    int m_admxIDX;

};

























//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){


    for(int i=0;i<MAXMAPPINGQUAL;i++){
	long double correctBaseProb     = 1.0-powl(10.0,((long double)(i))/-10.0); //1-m
	long double incorrectBaseProb   =     powl(10.0,((long double)(i))/-10.0); //m

        probMatch[i]        = correctBaseProb;
        probMismatch[i]     = incorrectBaseProb;
    }


    //Adding mismapping probability
    for(int m=0;m<MAXMAPPINGQUAL;m++){
	long double correctMappingProb     = 1.0-powl(10.0,((long double)(m))/-10.0); //1-m
	long double incorrectMappingProb   =     powl(10.0,((long double)(m))/-10.0); //m
	
	probMapping[m]    = correctMappingProb;    //1-m
	probMismapping[m] = incorrectMappingProb;  //m
    }
// #ifdef DEBUG1
// 	cerr<<"m\t"<<m<<"\t"<<incorrectMappingProb<<"\t"<<correctMappingProb<<endl;
// #endif
	
//     	for(int i=0;i<2;i++){
//     	    likeMatchMQ[m][i]           = log(  correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0   )/log(10);         
//     	    likeMismatchMQ[m][i]        = log(  correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0   )/log(10);
//     	    likeMatchProbMQ[m][i]       = correctMappingProb*(1.0-pow(10.0,2.0/-10.0)    ) + incorrectMappingProb/4.0;
//     	    likeMismatchProbMQ[m][i]    = correctMappingProb*(    pow(10.0,2.0/-10.0)/3.0) + incorrectMappingProb/4.0;
//     	}


//     	//Computing for quality scores 2 and up
//     	for(int i=2;i<MAXMAPPINGQUAL;i++){
// 	    //  (1-m)(1-e) + m/4  = 1-m-e+me +m/4  = 1+3m/4-e+me
//     	    likeMatchMQ[m][i]         = log(  correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0    )/log(10);    
// 	    //  (1-m)(e/3) + m/4  = e/3 -me/3 + m/4
//     	    likeMismatchMQ[m][i]      = log(  correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0    )/log(10);    
	    
//     	    likeMatchProbMQ[m][i]           = correctMappingProb*(1.0-pow(10.0,i/-10.0)    ) + incorrectMappingProb/4.0;
//     	    likeMismatchProbMQ[m][i]        = correctMappingProb*(    pow(10.0,i/-10.0)/3.0) + incorrectMappingProb/4.0;
//     	}


// #ifdef DEBUG1
//     	for(int i=0;i<MAXMAPPINGQUAL;i++){
// 	    cerr<<"m\t"<<m<<"\t"<<i<<"\t"<<likeMatchMQ[m][i]<<"\t"<<likeMismatchMQ[m][i]<<"\t"<<likeMatchProbMQ[m][i]<<"\t"<<likeMismatchProbMQ[m][i]<<endl;
// 	}
// #endif

//     }

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
    // int sizeGenome=0;

    string outTabPrefix     = "";
    bool   outTabPrefixBool = false;

    // string outlog  = "/dev/stderr";
    // string nameMT  = "MT";

    // string outSeqC  = "";
    // string outLogC  = "";
    // bool   outSeqCflag  = false;
    // bool   outLogCflag  = false;


    // string nameMTC  = "MTc";
    // bool userWantsContProduced = false;
  
    // int minQual                = 0;
    bool ignoreMQ              = false;

    // string   outLog  ;
    // ofstream outLogFP;
    bool twoPopMode   = false;
    bool threePopMode = false;

    


    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
    initScores();
    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////


    // long double locatione=0.0;
    // long double locationc=0.0;
    // long double scalee   =0.0;
    // long double scalec   =0.0;



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
    // bool singleStrand           = false;
    // bool doubleStrand           = false;

    // vector<string> anchVec;
    string anchPop;
    //vector<string> admxVec;
    string admxPop;
    vector<string> contVec;
    

    vector<int> anchVecIDX;
    int         anchIDX=-1;
    vector<int> admxVecIDX;
    int         admxIDX=-1;
    vector<int> contVecIDX;

    bool flagTSTV=false;
    bool useCpG=false;

    const string usage=string("\nThis program takes an aligned BAM file for a mitonchondria and calls a\nconsensus for the endogenous material\n\n\t"+
			      string(argv[0])+			      
			      " [options] [fasta file] [bam file] [region or file with regions to use] [freq for pop1] [freq for pop2] ... "+"\n\n"+
			      
			      "\t\t"+"-2p" +"\t\t\t\t"+"Use 2pop mode (default: none)"+"\n"+
                              "\t\t"+"-3p" +"\t\t\t\t"+"Use 3pop mode (default: none)"+"\n"+


                              "\t\t"+"-o  [output prefix]" +"\t\t"+"Output prefix  (default: none )"+"\n"+
			      
			      "\n\tFiltering options:\n"+
			      "\t\t"+"-wcpg" +"\t\t\t\t"+"Allow CpG sites  (default:  "+boolStringify(useCpG)+" )"+"\n"+
			      "\n\tDeamination options:\n"+				      
			      "\t\t"+"-t" +"\t\t\t\t"+"Flag transitions/transversions  (default:  "+boolStringify(flagTSTV)+" )"+"\n"+

			      // "\tThese options mark potentially deaminated sites\n"+
			      // "\t\t"+"--single"+  "\t\t\t\t"+"Use single-stranded damage patterns, only C->T     (default: none)"+"\n"+
			      //  "\t\t"+"--double"+  "\t\t\t\t"+"Use double-stranded damage patterns, C->T and G->A (default: none)"+"\n"+

			      "\n\tPopulation options:\n"+ 
			      "\tThe name of the populations much be the same (case sensitive) in the frequency files\n"+
			      "\t\t"+"--anch"+  "\t\t\t\t"+"Population to use as anchor                        (default: none)"+"\n"+
			      "\t\t"+"--cont"+  "\t\t\t\t"+"Comma-separated list of contaminant populations    (default: all)"+"\n"+
			      "\t\t"+"--admx"+  "\t\t\t\t"+"Population to use as admixing population           (default: none)"+"\n"+
			      
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




        if(string(argv[i]) == "-wcpg" ){
	    useCpG   = true;
            continue;
        }


        if(string(argv[i]) == "-t" ){
	    flagTSTV   = true;
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
            outTabPrefix=string(argv[i+1]);
	    outTabPrefixBool=true;
            i++;
            continue;
        }


	if( string(argv[i]) == "--anch"  ){
	    anchPop= string(argv[i+1]) ;
            i++;
            continue;
        }

	if( string(argv[i]) == "--cont"  ){
	    contVec=allTokens( string(argv[i+1]) ,',' );
            i++;
            continue;
        }

	if( string(argv[i]) == "--admx"  ){
	    //admxVec=allTokens( string(argv[i+1]) ,',');
	    admxPop = string(argv[i+1]);
            i++;
            continue;
        }




	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    // if(output == outlog){
    // 	cerr<<"Error: The sequence output is the same as the log"<<endl;
    // 	return 1;	
    // }

    if(!outTabPrefixBool){

	cerr<<"Error: The output prefix is mandatory"<<endl;
	return 1;	

    }


    // outLogFP.open(outLog.c_str());

    // if (!outLogFP.is_open()){
    //     cerr << "Unable to write to output file "<<outLog<<endl;
    //     return 1;
    // }

    vector<string>  filesAlFreq;
    string fastaFile       = string(argv[lastOpt]);   //fasta file
    string bamfiletopen    = string(argv[lastOpt+1]); //bam file
    //string fastaFile    = string(argv[argc-2]);//fasta file
    string regionStr       = string(argv[lastOpt+2]); //fasta file


    //cout<<lastOpt<<endl;
    cout<<fastaFile<<endl;
    cout<<bamfiletopen<<endl;
    cout<<regionStr<<endl;
   
    for(int i=(lastOpt+3);i<(argc);i++){ //all but the last 3 args
	cout<<"cont "<<string(argv[i])<<endl;
	filesAlFreq.push_back(string(argv[i]));	
    }




    ////////////////////////////////////
    // END Parsing arguments        //
    ////////////////////////////////////













    for(int nuc=0;nuc<12;nuc++){
	defaultSubMatch.s[ nuc ] = 0.0;
	//     else
	// 	defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    	// }    	
    }


    if(twoPopMode == false && threePopMode==false){
	cerr<<"Either specify -2p or -3p"<<endl;
        return 1;
    }

    if(twoPopMode
       &&
       !admxPop.empty() ){
       //!admxVec.empty() ){
	cerr<<"Cannot specify --admx with -2p"<<endl;
        return 1;
    }

    if(threePopMode
       &&
       admxPop.empty() ){
       //!admxVec.empty() ){
	cerr<<"Must specify --admx with -3p"<<endl;
        return 1;
    }


    if( anchPop.empty() ){
	cerr<<"Must specify the anchor population with --anch"<<endl;
        return 1;
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
















    cout<<"test\t"<<bamfiletopen<<"\t"<<regionStr<<endl;

     
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
    vector<region> regionVec;

    if(isFile(regionStr)){
	
	
	igzstream regionFP;

	regionFP.open(regionStr.c_str(), ios::in);

	//    unsigned int counterCont=0;
	if (regionFP.good()){
	    vector<string> fields;
	    string line;
	
	    while(getline (regionFP,line)){
		region toadd; 
		vector<string> tokens  = allTokens(line,':');
		if(tokens.size() != 2){
		    cerr << "Invalid format for line:  " <<line<< ", the line has "<<tokens.size()<<" tokens if split wrt to ':' in region file:"<<regionStr<<endl;
		    exit(1);
		}

		vector<string> tokens2 = allTokens(tokens[1],'-');
		if(tokens2.size() != 2){
		    cerr << "Invalid format for line:  " <<line<< ", the line has "<<tokens2.size()<<" tokens if split wrt to '-' in region file:"<<regionStr<<endl;
		    exit(1);
		}

		toadd.id         =                   tokens[0];
		toadd.leftCoord  = destringify<int>(tokens2[0]);
		toadd.rightCoord = destringify<int>(tokens2[1]);
		regionVec.push_back(toadd);
	    }
		             	              
	    regionFP.close();
	}else{
	    cerr << "Unable to open file "<<regionStr<<endl;
	    exit(1);
	}


    }else{
	region toadd; 
	vector<string> tokens  = allTokens(regionStr,':');
	if(tokens.size() != 2){
	    cerr << "Invalid format for line:  " <<regionStr<< ", the line has "<<tokens.size()<<" tokens if split wrt to ':' in region string"<<endl;
	    exit(1);
	}

	vector<string> tokens2 = allTokens(tokens[1],'-');
	if(tokens2.size() != 2){
	    
	    cerr << "Invalid format for line:  " <<regionStr<< ", the line has "<<tokens2.size()<<" tokens if split wrt to '-' in region string"<<endl;
	    exit(1);
	}

	toadd.id         =                   tokens[0];
	toadd.leftCoord  = destringify<int>(tokens2[0]);
	toadd.rightCoord = destringify<int>(tokens2[1]);
	regionVec.push_back(toadd);
    }

    if(regionVec.size() == 0){
	cerr<<"Error: at least one region must be defined"<<endl;
	return 1;	
    }


    vector<MistarParser * > vectorOfMP;
    int numberOfPopulationForCont=0;
    vector<string> namesForPops;

    for(int i=(lastOpt+3);i<(argc);i++){ 
	// if(i==1 && string(argv[i]) == "-f"){
	//     force=true;
	//     continue;
	// }
	// cout<<"MP "<<string(argv[i])<<endl;
	if(!isFile(string(argv[i])+".tbi")){
	    cerr<<"Error: The allele count file must be tabix indexed: "<<string(argv[i])<<".tbi file not found"<<endl;
	    return 1;	
	}

	MistarParser * mp = new MistarParser(string(argv[i]),(string(argv[i])+".tbi"),regionVec[0].id,regionVec[0].leftCoord,regionVec[0].rightCoord);
	//cout<<"mp addr "<<mp->hasData()<<endl;
	vectorOfMP.push_back(mp);

	numberOfPopulationForCont+=(mp->getPopulationsNames()->size() -2);//minus the root and anc
	for(unsigned int n=2;n<mp->getPopulationsNames()->size();n++){
	    namesForPops.push_back( mp->getPopulationsNames()->at(n) );
	}

    }    
    
    if(vectorOfMP.size() == 0){
	cerr<<"Error: File with allele frequencies must be defined"<<endl;
	return 1;	
    }
    
    // for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    // 	vectorOfMP[i]->repositionIterator(tokens[0], destringify<int>(tokens2[0]) ,  destringify<int>(tokens2[1]));
    // }


    //should not be empty
    // if(anchVec.empty()){
    // 	for(unsigned int n=0;n<namesForPops.size();n++){
    // 	    anchVec.push_back(namesForPops[n]);
    // 	}
    // }

    if(contVec.empty()){
	for(unsigned int n=0;n<namesForPops.size();n++){
	    contVec.push_back(namesForPops[n]);
	}
    }


    //should not be empty if 3p
    // if(admxVec.empty()){
    // 	for(unsigned int n=0;n<namesForPops.size();n++){
    // 	    admxVec.push_back(namesForPops[n]);
    // 	}
    // }



   for(unsigned int i1=0;i1<namesForPops.size();i1++){
       //for(unsigned int i2=0;i2<anchVec.size();i2++){
       if(namesForPops[i1] == anchPop){
	   anchIDX            = i1;
	   anchVecIDX.push_back(i1);
       }       
       //}       
    }

   if(anchIDX == -1){
       cerr<<"Cannot find specified anchor population: "<<anchPop<<" among the populations"<<endl;
       return 1;
   }

    for(unsigned int i1=0;i1<namesForPops.size();i1++){
	for(unsigned int i2=0;i2<contVec.size();i2++){
	    if(namesForPops[i1] == contVec[i2]){
		contVecIDX.push_back(i1);
	    }
	}       
    }

    for(unsigned int i1=0;i1<namesForPops.size();i1++){
	//for(unsigned int i2=0;i2<admxVec.size();i2++){
	//if(namesForPops[i1] == admxVec[i2]){
	//	admxVecIDX.push_back(i1);
	//   }
	//}

	if(namesForPops[i1] == admxPop){
	    admxIDX            = i1;
	    admxVecIDX.push_back(i1);
	}
    }

    if(threePopMode)
	if(admxIDX == -1){
	    cerr<<"Cannot find specified admixed population: "<<admxIDX<<" among the populations"<<endl;
	    return 1;
	}


    // cout<<vectorToString(anchVecIDX)<<endl;
    // cout<<vectorToString(contVecIDX)<<endl;
    // cout<<vectorToString(admxVecIDX)<<endl;

    // return 1;

    // vector<int> anchVecIDX;
    // vector<int> admxVecIDX;
    // vector<int> contVecIDX;

    // for(unsigned int n=0;n<namesForPops.size();n++){
    // 	cout<<n<<"\t"<<namesForPops[n]<<endl;
    // }
    // for(unsigned int indexCont=0;indexCont<contVecIDX.size();indexCont++){
    // 	for(unsigned int indexAnchor=0;indexAnchor<anchVecIDX.size();indexAnchor++){
    // 	    for(unsigned int indexAdmix=0;indexAdmix<admxVecIDX.size();indexAdmix++){
    // 		cout<<contVecIDX[indexCont]<<"\t"<<anchVecIDX[indexAnchor]<<"\t"<<admxVecIDX[indexAdmix]<<endl;
    // 	    }
    // 	}
    // }
    // return 1;


    Fasta fastaReference;
    if ( !fastaReference.Open(fastaFile , fastaFile+".fai") ){ 
	cerr << "ERROR: failed to open fasta file " <<fastaFile<<" and index " << fastaFile<<".fai"<<endl;
	exit(1);
    }
    // cout<<"fine1"<<endl;






    vector<singleSite> *  dataSitesVec = new vector<singleSite>();

    for(unsigned int regionIdx=0;regionIdx<regionVec.size();regionIdx++){
	int id = reader.GetReferenceID( regionVec[regionIdx].id );

	cerr<<"Processing region :"<<"\t"<<regionVec[regionIdx].id<<"\t"<<regionVec[regionIdx].leftCoord<<"\t"<<regionVec[regionIdx].rightCoord<<endl;

	BamRegion regionbam ( id , regionVec[regionIdx].leftCoord , id, regionVec[regionIdx].rightCoord );    
	reader.SetRegion(regionbam);
	//cout<<"fine2 id "<<id<<endl;
	for(unsigned int mstIndex=0;mstIndex<vectorOfMP.size();mstIndex++){
	    vectorOfMP[mstIndex]->repositionIterator( regionVec[regionIdx].id , regionVec[regionIdx].leftCoord , regionVec[regionIdx].rightCoord );
	}
	//cerr<<"fine3 id "<<id<<endl;
	
	bool wasDataFound=true;

	MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,vectorOfMP, regionVec[regionIdx].leftCoord ,regionVec[regionIdx].rightCoord,  ignoreMQ,dataSitesVec , useCpG, &wasDataFound,anchIDX,admxIDX);
	//cerr<<"fine4 id "<<id<<endl;
	PileupEngine pileup;
	pileup.AddVisitor(cv);
	if(!wasDataFound){
	    cerr<<"No data was found in the frequency files"<<endl;
	    continue;
	}


	//cerr<<"fine5 id "<<id<<endl;
	BamAlignment al;
	unsigned int numReads=0;
	//cerr<<"Reading BAM file ..."<<endl;
	while ( reader.GetNextAlignment(al) ) {
	    //cout<<"name:\t"<<al.Name<<endl;
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
	// cerr<<"...  done"<<endl;
    
    
	//clean up
	pileup.Flush();
	delete cv;

    }

    reader.Close();
    cerr<<"done reading bam files"<<endl;
    //return 0;


    //if(outTabPrefixBool){
    if(twoPopMode){


	for(unsigned int indexCont=0;indexCont<contVecIDX.size();indexCont++){
	    for(unsigned int indexAnchor=0;indexAnchor<anchVecIDX.size();indexAnchor++){
	    //for(int indexAnchor=anchIDX;indexAnchor<=anchIDX;indexAnchor++){

		if(contVecIDX[indexCont] == anchVecIDX[indexAnchor]){	       

		    map<string,int> stringToPrint2Count;
		    for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){		
			//cout<<dataSitesVec->at(sitesIdx).ancAllele<<"\t"<<dataSitesVec->at(sitesIdx).derAllele<<"\t"<<endl;
			string tmpString = 
			    stringify(dataSitesVec->at(sitesIdx).sitesBAMa.size())+"\t"+
			    stringify(dataSitesVec->at(sitesIdx).sitesBAMd.size())+"\t"+
			    stringify(dataSitesVec->at(sitesIdx).freqDerived[ contVecIDX[indexCont] ]);
			
			if(flagTSTV){
			    tmpString+="\t"+stringify(isPotentialTransition(dataSitesVec->at(sitesIdx).ancAllele,dataSitesVec->at(sitesIdx).derAllele));
			}

			stringToPrint2Count[ tmpString  ]++;
			//"1"<<endl;
		    }


		    string filenameTab=outTabPrefix+"_Cont_Anch_"+namesForPops[ contVecIDX[indexCont] ]+".dice";
		    ofstream outTabFP;
		    
		    outTabFP.open(filenameTab.c_str());
		
		    if (!outTabFP.is_open()){
			cerr << "Unable to write to output file "<<filenameTab<<endl;
			return 1;
		    }

		    //cout<<"C "<<namesForPops[indexCont]<<"\tA "<<namesForPops[indexAnchor]<<endl;
		    outTabFP<<"Anc\tDer\tContAnchFreq\tNum";
		    if(flagTSTV)
			outTabFP<<"\tTs";

		    outTabFP<<endl;

		    for( map<string,int>::iterator itSTPC=stringToPrint2Count.begin(); 
			 itSTPC!=stringToPrint2Count.end(); ++itSTPC){
			//outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			if(flagTSTV){
			    //outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			    vector<string> s=allTokens(itSTPC->first,'\t');
			    vector<string> sToP;
			    for(unsigned int si=0;si<(s.size()-1);si++){
				sToP.push_back(s[si]);
			    }
			    sToP.push_back( stringify(itSTPC->second) );
			    sToP.push_back( s[s.size()-1] );
			    outTabFP << vectorToString(sToP,"\t") << endl;
			}else{
			    outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			}
		    }
		    outTabFP.close();

		}else{ //anchor is not cont

		    map<string,int> stringToPrint2Count;
		    for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){		

			string tmpString = 
			    stringify(dataSitesVec->at(sitesIdx).sitesBAMa.size())+"\t"+
			    stringify(dataSitesVec->at(sitesIdx).sitesBAMd.size())+"\t"+
			    stringify(dataSitesVec->at(sitesIdx).freqDerived[ anchVecIDX[indexAnchor] ])+"\t"+ 
			    stringify(dataSitesVec->at(sitesIdx).freqDerived[ contVecIDX[indexCont]   ]);

			if(flagTSTV){
			    tmpString+="\t"+stringify(isPotentialTransition(dataSitesVec->at(sitesIdx).ancAllele,dataSitesVec->at(sitesIdx).derAllele));
			}

			stringToPrint2Count[ tmpString ]++;
			//"1"<<endl;
		    }


		    string filenameTab=outTabPrefix+"_Cont_"+namesForPops[ contVecIDX[indexCont] ]+"_Anch_"+namesForPops[ anchVecIDX[indexAnchor] ]+".dice";
		    ofstream outTabFP;
		    
		    outTabFP.open(filenameTab.c_str());
		
		    if (!outTabFP.is_open()){
			cerr << "Unable to write to output file "<<filenameTab<<endl;
			return 1;
		    }

		    //cout<<"C "<<namesForPops[indexCont]<<"\tA "<<namesForPops[indexAnchor]<<endl;
		    outTabFP<<"Anc\tDer\tAnchFreq\tContFreq\tNum";
		    if(flagTSTV)
			outTabFP<<"\tTs";

		    outTabFP<<endl;

		    for( map<string,int>::iterator itSTPC=stringToPrint2Count.begin(); 
			 itSTPC!=stringToPrint2Count.end(); ++itSTPC){
			//outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			if(flagTSTV){
			    //outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			    vector<string> s=allTokens(itSTPC->first,'\t');
			    vector<string> sToP;
			    for(unsigned int si=0;si<(s.size()-1);si++){
				sToP.push_back(s[si]);
			    }
			    sToP.push_back( stringify(itSTPC->second) );
			    sToP.push_back( s[s.size()-1] );
			    outTabFP << vectorToString(sToP,"\t") << endl;
			}else{
			    outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
			}
		    }

		    outTabFP.close();

		}


	    }
	}
	//return 0;
    }else{
	if(threePopMode){

	    for(unsigned int indexCont=0;indexCont<contVecIDX.size();indexCont++){
		for(unsigned int indexAnchor=0;indexAnchor<anchVecIDX.size();indexAnchor++){
		    //for(int indexAnchor=anchIDX;indexAnchor<=anchIDX;indexAnchor++){
		   for(unsigned int indexAdmix=0;indexAdmix<admxVecIDX.size();indexAdmix++){
			
		       
		       //if(indexCont == indexAdmix){
		       if(contVecIDX[indexCont] == admxVecIDX[indexAdmix]){

			    map<string,int> stringToPrint2Count;
			    for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){		
				string tmpString = 
				    stringify(dataSitesVec->at(sitesIdx).sitesBAMa.size())+"\t"+
				    stringify(dataSitesVec->at(sitesIdx).sitesBAMd.size())+"\t"+
				    stringify(dataSitesVec->at(sitesIdx).freqDerived[  contVecIDX[indexCont]   ])+"\t"+ 
				    stringify(dataSitesVec->at(sitesIdx).freqDerived[  anchVecIDX[indexAnchor] ]);
				if(flagTSTV){
				    tmpString+="\t"+stringify(isPotentialTransition(dataSitesVec->at(sitesIdx).ancAllele,dataSitesVec->at(sitesIdx).derAllele));
				}

				stringToPrint2Count[ tmpString ]++;
				//"1"<<endl;
			    }
			    
			    
			    string filenameTab=outTabPrefix+"_ContAdmx_"+namesForPops[ contVecIDX[indexCont] ]+"_Anch_"+namesForPops[ anchVecIDX[indexAnchor] ]+".dice";
			    ofstream outTabFP;
		    
			    outTabFP.open(filenameTab.c_str());
		
			    if (!outTabFP.is_open()){
				cerr << "Unable to write to output file "<<filenameTab<<endl;
				return 1;
			    }

			    //cout<<"C "<<namesForPops[indexCont]<<"\tA "<<namesForPops[indexAnchor]<<endl;
			    outTabFP<<"Anc\tDer\tContAdmxFreq\tAnchFreq\tNum";
			    if(flagTSTV)
				outTabFP<<"\tTs";

			    outTabFP<<endl;

			    for( map<string,int>::iterator itSTPC=stringToPrint2Count.begin(); 
				 itSTPC!=stringToPrint2Count.end(); ++itSTPC){
				//outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
				if(flagTSTV){
				    //outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
				    vector<string> s=allTokens(itSTPC->first,'\t');
				    vector<string> sToP;
				    for(unsigned int si=0;si<(s.size()-1);si++){
					sToP.push_back(s[si]);
				    }
				    sToP.push_back( stringify(itSTPC->second) );
				    sToP.push_back( s[s.size()-1] );
				    outTabFP << vectorToString(sToP,"\t") << endl;
				}else{
				    outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
				}
			    }

			    outTabFP.close();
			    
			}else{//if contaminant not is the admixed pop
			 
			    map<string,int> stringToPrint2Count;
			    for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){		
				string tmpString = 
				    stringify(dataSitesVec->at(sitesIdx).sitesBAMa.size())+"\t"+
				    stringify(dataSitesVec->at(sitesIdx).sitesBAMd.size())+"\t"+
				    stringify(dataSitesVec->at(sitesIdx).freqDerived[ admxVecIDX[indexAdmix]  ])+"\t"+ 
				    stringify(dataSitesVec->at(sitesIdx).freqDerived[ anchVecIDX[indexAnchor] ])+"\t"+ 
				    stringify(dataSitesVec->at(sitesIdx).freqDerived[ contVecIDX[indexCont]   ]);

				if(flagTSTV){
				    tmpString+="\t"+stringify(isPotentialTransition(dataSitesVec->at(sitesIdx).ancAllele,dataSitesVec->at(sitesIdx).derAllele));
				}

				stringToPrint2Count[ tmpString ]++;
				//"1"<<endl;
			    }

			    string filenameTab=outTabPrefix+"_Admx_"+namesForPops[ admxVecIDX[indexAdmix] ]+"_Anch_"+namesForPops[ anchVecIDX[indexAnchor] ]+"_Cont_"+namesForPops[ contVecIDX[indexCont] ]+".dice";
			    ofstream outTabFP;
		    
			    outTabFP.open(filenameTab.c_str());
		
			    if (!outTabFP.is_open()){
				cerr << "Unable to write to output file "<<filenameTab<<endl;
				return 1;
			    }

			    //cout<<"C "<<namesForPops[indexCont]<<"\tA "<<namesForPops[indexAnchor]<<endl;
			    outTabFP<<"Anc\tDer\tAdmxFreq\tAnchFreq\tContFreq\tNum";
			    if(flagTSTV)
				outTabFP<<"\tTs";				
			    outTabFP<<endl;
			    
			    for( map<string,int>::iterator itSTPC=stringToPrint2Count.begin(); 
				 itSTPC!=stringToPrint2Count.end(); ++itSTPC){
				if(flagTSTV){
				    //outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
				    vector<string> s=allTokens(itSTPC->first,'\t');
				    vector<string> sToP;
				    for(unsigned int si=0;si<(s.size()-1);si++){
					sToP.push_back(s[si]);
				    }
				    sToP.push_back( stringify(itSTPC->second) );
				    sToP.push_back( s[s.size()-1] );
				    outTabFP << vectorToString(sToP,"\t") << endl;
				}else{
				    outTabFP << itSTPC->first << "\t" << (itSTPC->second) << endl;
				}
			    }

			    outTabFP.close();
			    
			}

		    
		    }
		}
	    }
	    
	}else{
	    cerr<<"Invalid state"<<endl;
	    return 1;
	}
    }


    return 0;
}

