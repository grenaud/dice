/*
 * diceBAM
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

#include "utils.h"
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
		  vector<singleSite> * dataSitesVec
		  // const long double contaminationPrior,
		  //const bool singleCont
		  )
    : PileupVisitor()
      // , m_references(references)
    , m_fastaReference(fastaReference)
    , m_vectorOfMP(vectorOfMP)
    , m_coordFirst(coordFirst)
    , m_coordLast(coordLast)
    , ignoreMQ(ignoreMQ)
    , m_dataSitesVec(dataSitesVec)
    // , contaminationPrior(contaminationPrior)
    // , singleCont(singleCont)
  { 
     
      //      cout<<"constr size="<<vectorOfMP.size()<<endl;
      initFiles(vectorOfMP,
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


	for(unsigned int df=0;df<derFreq.size();df++){
	    if(derFreq[df] <= 0 || derFreq[df] >= 1.0)//skip sites with fixed bases
	    	return ;
	}

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
	toAddSS.ancAllele = ancAllele;
	toAddSS.derAllele = derAllele;
	toAddSS.chr       = chr1;
	toAddSS.coord     = posAlign;


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


#ifdef DEBUGERRORP	    
	    cout<<"i="<<i<<"\t"<<posAlign<<"\tref="<<referenceBase<<"\tb_obs="<<b<<"\tq="<<sb.bq<<"\t"<<pileupData.PileupAlignments[i].Alignment.Name<<"\t"<<"\td="<<derAllele<<"\ta="<<ancAllele <<"\t"<<dist5p<<"\t"<<dist3p<<"\t"<<(pileupData.PileupAlignments[i].Alignment.IsReverseStrand()?"R":"F")<<endl;


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
    vector<singleSite> *  m_dataSitesVec;
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
    string output  = "/dev/stdout";
    // string outTabPrefix     = "";
    // bool   outTabPrefixBool = false;

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

    string   outLog  = "/dev/stdout";
    ofstream outLogFP;
    bool twoPopMode   = false;
    bool threePopMode = false;
    double step = 1000;
    int maxChains = 100000;


    //Constants
    long double innerdriftY =   0.16;
    long double innerdriftZ =   0.16;
    long double nC          =  100.0 ;
    long double nB          =  100.0 ;
    


    // Set lower boundaries for optimization algorithm
    // long double elower         = 0.00001;
    long double rlower         = 0.00001;
    long double tau_Clower     = 0.000001;
    long double tau_Alower     = 0.000001;
    long double admixratelower = 0.000001;
    long double admixtimelower = 0.05;

    // Set upper boundaries for optimization algorithm
    // long double eupper         = 0.1;
    long double rupper         = 0.5;
    long double tau_Cupper     = 1.0;
    long double tau_Aupper     = 1.0;
    long double admixrateupper = 0.5;
    long double admixtimeupper = 0.11;

    // long double e_i         ;
    long double r_i         ;
    long double tau_C_i     ;
    long double tau_A_i     ;
    long double admixrate_i ;
    long double admixtime_i ;

    // bool e_i_0         = false;
    bool r_i_0         = false;
    bool tau_C_i_0     = false;
    bool tau_A_i_0     = false;
    bool admixrate_i_0 = false;
    bool admixtime_i_0 = false;

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
    string anchPop;
    string admxPop;
    string contPop;

    int anchPopIDX;
    int admxPopIDX;
    int contPopIDX;


    const string usage=string("\nThis program takes an aligned BAM file for a mitonchondria and calls a\nconsensus for the endogenous material\n\n\t"+
			      string(argv[0])+			      
			      " [options] [fasta file] [bam file] [region or file with regions to use] [freq for pop1] [freq for pop2] ... "+"\n\n"+
			      
			      "\t\t"+"-2p" +"\t\t\t\t"+"Use 2pop mode (default: none)"+"\n"+
                              "\t\t"+"-3p" +"\t\t\t\t"+"Use 3pop mode (default: none)"+"\n"+

                              "\t\t"+"-o     [output log]" +"\t\t"+"Output log (default: stdout)"+"\n"+
                              // "\t\t"+"-otab  [output tab]" +"\t\t"+"Produce a tab delimited output, enter the prefix here but no MCMC computation (default: none)"+"\n"+
			      
                              "\n\tComputation options:\n"+
                              "\t\t"+"-s     [step]" +"\t\t\t"+"MCMC interval space step (default: "+stringify(step)+")"+"\n"+
                              "\t\t"+"-c     [#chains]" +"\t\t"+"Max. number of Markov chains (default: "+stringify(maxChains)+")"+"\n"+

                              "\n\tStarting values:\n"+
			      "\t\t"+"-e0     [error]"+"\t\t\t"+"Error rate         (default: random)"+"\n"+
			      "\t\t"+"-r0     [cont]" +"\t\t\t"+"Contamination rate (default: random)"+"\n"+
			      "\t\t"+"-tA0    [tauA]" +"\t\t\t"+"Tau ancient genome        (default: random)"+"\n"+
			      "\t\t"+"-tC0    [tauC]" +"\t\t\t"+"Tau anchor    (default: random)"+"\n"+
			      "\t\t"+"-aR0    [admR]" +"\t\t\t"+"Admixture time     (default: random)"+"\n"+
			      "\t\t"+"-aT0    [admT]" +"\t\t\t"+"Admixture rate     (default: random)"+"\n"+
                              
			      "\n\tRange for parameter values:\n"+
			      // "\t\t"+"-e     el,eh"+"\t\t\t"+"Error rate range          (default: "+stringify(elower)         +","+stringify(eupper)+" )"+"\n"+
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
			      
			      // "The BAM file has to be indexed and the region has to be in chr:start-end format\n\n"+
			      // "\n\tOutput options:\n"+	
			      // "\t\t"+"-out  [output]" +"\t\t"+"Output file (default: stdout)"+"\n"+
			      // "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
			      // "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
			      // "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQual)+") "+"\n"+

			      "\n\tDeamination options:\n"+				      
			      "\t\t"+"-deam5p  [.prof file]" +"\t\t"+"5p deamination frequency for the endogenous  (default: "+deam5pfreqE+")"+"\n"+
			      "\t\t"+"-deam3p  [.prof file]" +"\t\t"+"3p deamination frequency for the endogenous  (default: "+deam3pfreqE+")"+"\n"+
			      "\t\t"+"-deam5pc [.prof file]" +"\t\t"+"5p deamination frequency for the contaminant (default: "+deam5pfreqC+")"+"\n"+
			      "\t\t"+"-deam3pc [.prof file]" +"\t\t"+"3p deamination frequency for the contaminant (default: "+deam3pfreqC+")"+"\n"+


			      "\n\tPopulation options:\n"+ 
			      "\tThe name of the populations much be the same (case sensitive) in the frequency files\n"+
			      "\t\t"+"--anch"+  "\t\t\t\t"+"Anchor population         (default: all)"+"\n"+
			      "\t\t"+"--cont"+  "\t\t\t\t"+"Contaminant population    (default: all)"+"\n"+
			      "\t\t"+"--admx"+  "\t\t\t\t"+"Admixing populations      (default: all)"+"\n"+
			      
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

			      "\n\tComputation options:\n"+	
			      // "\t\t"+"-nomq" +"\t\t\t\t"+"Ignore mapping quality (default: "+booleanAsString(ignoreMQ)+")"+"\n"+
			      // "\t\t"+"-err" +"\t\t\t\t"+"Illumina error profile (default: "+errFile+")"+"\n"+
			      "\t\t"+"--phred64" +"\t\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+

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
    string cwdProg=getCWD(argv[0]);

    for(int i=1;i<(argc);i++){ 

	if(string(argv[i])[0] != '-'  ){
	    lastOpt=i;
	    break;
	}

	if(string(argv[i]) == "--phred64"  ){
	    offsetQual=64;
	    continue;
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

	// if( string(argv[i]) == "-otab"  ){
        //     outTabPrefix=string(argv[i+1]);
	//     outTabPrefixBool=true;
        //     i++;
        //     continue;
        // }

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

	if(string(argv[i]) == "-o" ){
	    output=string(argv[i+1]);
	    i++;
	    continue;
	}



	if( string(argv[i]) == "--anch"  ){
	    anchPop= string(argv[i+1]) ;
            i++;
            continue;
        }

	if( string(argv[i]) == "--cont"  ){
	    contPop= string(argv[i+1]) ;
            i++;
            continue;
        }

	if( string(argv[i]) == "--admx"  ){
	    admxPop= string(argv[i+1]) ;
            i++;
            continue;
        }




	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }


    if(twoPopMode
       &&
       !admxPop.empty() ){
	cerr<<"Cannot specify --admx with -2p"<<endl;
        return 1;
    }

    // if(output == outlog){
    // 	cerr<<"Error: The sequence output is the same as the log"<<endl;
    // 	return 1;	
    // }

    outLogFP.open(outLog.c_str());

    if (!outLogFP.is_open()){
        cerr << "Unable to write to output file "<<outLog<<endl;
        return 1;
    }

    vector<string>  filesAlFreq;
    string fastaFile       = string(argv[lastOpt]);   //fasta file
    string bamfiletopen    = string(argv[lastOpt+1]); //bam file
    //string fastaFile     = string(argv[argc-2]);    //fasta file
    string regionStr       = string(argv[lastOpt+2]); //fasta file


    cout<<lastOpt<<endl;
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











    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN READING ERROR PROFILE
    //
    ////////////////////////////////////////////////////////////////////////

    // readIlluminaError(errFile,illuminaErrorsProb);

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
    readNucSubstitionRatesFreq(deam5pfreqE,sub5p);
    readNucSubstitionRatesFreq(deam3pfreqE,sub3p);
    readNucSubstitionRatesFreq(deam5pfreqC,sub5pC);
    readNucSubstitionRatesFreq(deam3pfreqC,sub3pC);

    // for(int l=0;l<sub5p.size();l++){
    // 	cout<<l<<"";
    // 	for(int i=0;i<16;i++){
    // 	    cout<<"\t"<<sub5p[l].s[i]<<endl;
    // 	}
    // }

    // cout<<sub5p.size()<<endl;
    // cout<<sub3p.size()<<endl;

    // return 1;
    

    // int defaultSubMatchIndex=0;

    // for(int nuc1=0;nuc1<4;nuc1++){
    // 	for(int nuc2=0;nuc2<4;nuc2++){	    
    // 	    if(nuc1==nuc2)
    // 		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 1.0;
    // 	    else
    // 		defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    // 	}
    	
    // }

    for(int nuc=0;nuc<12;nuc++){
	defaultSubMatch.s[ nuc ] = 0.0;
	//     else
	// 	defaultSubMatch.s[ defaultSubMatchIndex++ ] = 0.0;
    	// }    	
    }


    if(twoPopMode == false && threePopMode==false){
	cerr<<"Either specify -2p or -3p "<<endl;
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
		    cerr << "Invalid format for line:  " <<line<< endl;
		    exit(1);
		}

		vector<string> tokens2 = allTokens(tokens[1],'-');
		if(tokens2.size() != 2){
		    cerr << "Invalid format for line:  " <<line<< endl;
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
	    cerr << "Invalid format for region:  " <<regionStr<< endl;
	    exit(1);
	}

	vector<string> tokens2 = allTokens(tokens[1],'-');
	if(tokens2.size() != 2){
	    cerr << "Invalid format for region:  " <<regionStr<< endl;
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
    vector<string> namesCont;

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
	    namesCont.push_back( mp->getPopulationsNames()->at(n) );
	}

    }    
    
    if(vectorOfMP.size() == 0){
	cerr<<"Error: File with allele frequencies must be defined"<<endl;
	return 1;	
    }



    if(!anchPop.empty()){
	for(unsigned int n=0;n<namesCont.size();n++){
	    if(anchPop==namesCont[n])
		anchPopIDX=int(n);
	}
    }

    if(!contPop.empty()){
	for(unsigned int n=0;n<namesCont.size();n++){
	    if(contPop==namesCont[n])
		contPopIDX=int(n);
	}
    }


    if(!admxPop.empty()){
	for(unsigned int n=0;n<namesCont.size();n++){
	    if(admxPop==namesCont[n])
		admxPopIDX=int(n);
	}
    }
    // for(unsigned int i=0;i<vectorOfMP.size();i++){ 
    // 	vectorOfMP[i]->repositionIterator(tokens[0], destringify<int>(tokens2[0]) ,  destringify<int>(tokens2[1]));
    // }






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
	// cout<<"fine3 id "<<id<<endl;
	
	MyPileupVisitor* cv = new MyPileupVisitor(references,&fastaReference,vectorOfMP, regionVec[regionIdx].leftCoord ,regionVec[regionIdx].rightCoord,  ignoreMQ,dataSitesVec);
	// cout<<"fine4 id "<<id<<endl;
	PileupEngine pileup;
	pileup.AddVisitor(cv);

	// cout<<"fine5 id "<<id<<endl;
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


    // if(outTabPrefixBool){
    // 	for(int indexCont=0;indexCont<numberOfPopulationForCont;indexCont++){
    // 	    for(int indexAnchor=0;indexAnchor<numberOfPopulationForCont;indexAnchor++){

    // 		string filenameTab=outTabPrefix+"_C_"+namesCont[indexCont]+"_A_"+namesCont[indexAnchor]+".dice";
    // 		ofstream outTabFP;

    // 		outTabFP.open(filenameTab.c_str());
		
    // 		if (!outTabFP.is_open()){
    // 		    cerr << "Unable to write to output file "<<filenameTab<<endl;
    // 		    return 1;
    // 		}

    // 		//cout<<"C "<<namesCont[indexCont]<<"\tA "<<namesCont[indexAnchor]<<endl;
    // 		if(indexCont == indexAnchor){	       
    // 		    outTabFP<<"Anc\tDer\tPanelFreq\tNum"<<endl;
    // 		    for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){		
    // 			outTabFP<<
    // 			    dataSitesVec->at(sitesIdx).sitesBAMa.size()<<"\t"<<
    // 			    dataSitesVec->at(sitesIdx).sitesBAMd.size()<<"\t"<<
    // 			    dataSitesVec->at(sitesIdx).freqDerived[indexCont]<<"\t"<<
    // 			    "1"<<endl;
    // 		    }
    // 		}else{
    // 		    //TODO , also do 3 pop
    // 		}
    // 		outTabFP.close();

    // 	    }
    // 	}
    // 	return 0;
    // }


    // cout<<"Anc\tDer\tPanelFreq\tNum"<<endl;
    // for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){
    // 	//cout<<dataSitesVec->at(sitesIdx).sitesBAMa.size()<<"\t"<<dataSitesVec->at(sitesIdx).sitesBAMd.size()<<"\t"<<vectorToString(dataSitesVec->at(sitesIdx).freqDerived,"\t")<<"\t"<<1<<endl;
    // 	cout<<dataSitesVec->at(sitesIdx).sitesBAMa.size()<<"\t"<<dataSitesVec->at(sitesIdx).sitesBAMd.size()<<"\t"<<vectorToString(dataSitesVec->at(sitesIdx).freqDerived,"\t")<<"\t"<<"\t"<<dataSitesVec->at(sitesIdx).chr<<"\t"<<dataSitesVec->at(sitesIdx).coord<<"\t"<<dataSitesVec->at(sitesIdx).ancAllele<<"\tanc"<<endl;
    // 	cout<<dataSitesVec->at(sitesIdx).sitesBAMa.size()<<"\t"<<dataSitesVec->at(sitesIdx).sitesBAMd.size()<<"\t"<<vectorToString(dataSitesVec->at(sitesIdx).freqDerived,"\t")<<"\t"<<"\t"<<dataSitesVec->at(sitesIdx).chr<<"\t"<<dataSitesVec->at(sitesIdx).coord<<"\t"<<dataSitesVec->at(sitesIdx).derAllele<<"\tder"<<endl;
    // }



    //vector Variables
   // if(!e_i_0)
   //     e_i       = randomLongDouble(elower,         eupper);

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


    // long double e_i_1;
    long double r_i_1;
    long double tau_C_i_1;
    long double tau_A_i_1;
    long double admixrate_i_1;
    long double admixtime_i_1;




   long double x_il;
   long double x_i_1l;
   //long double llik = LogFinalTwoPBAM( dataSitesVec, 0.05, 0.45, 0.5, indexCont, indexAnchor );

   outLogFP<<"chain"<<"\tlpos"<<"\tContRate"<<"\ttau_C"<<"\ttau_A"<<"\tadmixrate"<<"\tadmixtime\tacceptance"<<endl;

    // for(int indexCont=0;indexCont<numberOfPopulationForCont;indexCont++){
    // 	for(int indexAnchor=0;indexAnchor<numberOfPopulationForCont;indexAnchor++){
   //cout<<"C "<<namesCont[indexCont]<<"\tA "<<namesCont[indexAnchor]<<endl;
   
   if(twoPopMode){
       x_il     = LogFinalTwoPBAM( dataSitesVec, r_i, tau_C_i,tau_A_i, contPopIDX, anchPopIDX );
   }else{
       if(threePopMode){
	   //LogFinalThreeP(dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,admixrate_i_1,admixtime_i_1,innerdriftY,innerdriftZ,nC,nB,cwdProg,has5Cols );
	   x_il = LogFinalThreePBAM( dataSitesVec, r_i, tau_C_i,tau_A_i,admixrate_i_1,admixtime_i_1, innerdriftY,innerdriftZ,nC,nB,cwdProg, contPopIDX, anchPopIDX, admxPopIDX );
       }else{
	   cerr<<"Invalid state"<<endl;
	   return 1;
       }
   }

   //cout<<"init\t"<<"\t"<<std::setprecision(10)<<x_il<<"\t"<<e_i<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<"\t"<<admixrate_i<<"\t"<<admixtime_i<<endl;
   // cout<<x_il<<endl;
	    
   int accept=0;

   random_device rd;
   default_random_engine dre (rd());

   for(int chain=0;chain<=maxChains;chain++){
		
       long double partition= (long double)(step);


		
       // //e
       // normal_distribution<long double> distribution_e(e_i,     (eupper-elower)/partition  );
       // e_i_1      = distribution_e(dre);
       // // e_i_1      = e_i;

       // if(e_i_1 <= elower     ||  e_i_1 >= eupper     ){
       //     e_i_1      = e_i;
       //     //chain--;
       //     //continue;
       // }

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
      


		
       // if(chain!=0)
       outLogFP<<chain<<"\t"<<std::setprecision(10)<<x_il<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<"\t"<<admixrate_i<<"\t"<<admixtime_i<<"\t"<<double(accept)/double(chain)<<endl;
       //cout<<"currt\t"<<chain<<"\t"<<std::setprecision(10)<<x_il<<"\t"<<"\t"<<r_i<<"\t"<<tau_C_i<<"\t"<<tau_A_i<<"\t"<<admixrate_i<<"\t"<<admixtime_i<<"\t"<<double(accept)/double(chain)<<endl;
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
	   //x_i_1l    = LogFinalTwoP(  dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,                                                                  has4Cols );
	   x_i_1l        = LogFinalTwoPBAM(   dataSitesVec, r_i_1, tau_C_i_1,tau_A_i_1, contPopIDX, anchPopIDX );

       }else{
	   //x_i_1l    = LogFinalThreeP(dataToAdd,e_i_1,r_i_1,tau_C_i_1,tau_A_i_1,admixrate_i_1,admixtime_i_1,innerdriftY,innerdriftZ,nC,nB,cwdProg,has5Cols );
	   //x_il = LogFinalThreePBAM( dataSitesVec, r_i, tau_C_i,tau_A_i,admixrate_i_1,admixtime_i_1, innerdriftY,innerdriftZ,nC,nB,cwdProg, contPopIDX, anchPopIDX, admxPopIDX );
	   if(threePopMode){
	       x_i_1l    = LogFinalThreePBAM( dataSitesVec, r_i_1, tau_C_i_1,tau_A_i_1,admixrate_i_1,admixtime_i_1,innerdriftY,innerdriftZ,nC,nB,cwdProg, contPopIDX, anchPopIDX, admxPopIDX);
	   }else{
	       cerr<<"Invalid state"<<endl;
	       return 1;

	   }
		    
       }

       long double acceptance = min( (long double)(1.0)  , expl(x_i_1l-x_il) );

       //		cout<< "new   "<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;
       // outLogFP<< "ratio "<<std::setprecision(10)<<expl(x_i_1l-x_il)<<"\tnew "<<(x_i_1l)<<"\told "<<(x_il)<<"\t"<<(x_i_1l-x_il)<<"\t"<<acceptance<<endl;

       // outLogFP<<chain<<"p\t"<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;
       //cout<< "ratio "<<std::setprecision(10)<<expl(x_i_1l-x_il)<<"\tnew "<<(x_i_1l)<<"\told "<<(x_il)<<"\t"<<(x_i_1l-x_il)<<"\t"<<acceptance<<endl;

       //cout<<chain<<"p\t"<<std::setprecision(10)<<x_i_1l<<"\t"<<e_i_1<<"\t"<<r_i_1<<"\t"<<tau_C_i_1<<"\t"<<tau_A_i_1<<"\t"<<admixrate_i_1<<"\t"<<admixtime_i_1<<"\t"<<acceptance<<endl;

       if( (long double)(randomProb()) < acceptance){
	   // e_i           =  e_i_1;
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

   // cout<<"C"<<namesCont[indexCont]<<"\tA"<<namesCont[indexAnchor]<<endl;
   // long double llik = LogFinalTwoPBAM( dataSitesVec, 0.05, 0.45, 0.5, indexCont, indexAnchor );
   // cout<<llik<<endl;
   // 	}
   // }
   // for(unsigned int sitesIdx=0;sitesIdx<dataSitesVec->size();sitesIdx++){
	

   // 	cout<<dataSitesVec->at(sitesIdx).ancAllele<<"\t"<<dataSitesVec->at(sitesIdx).derAllele<<"\t"<<dataSitesVec->at(sitesIdx).sitesBAMd.size()<<"\t"<<dataSitesVec->at(sitesIdx).sitesBAMa.size()<<"\t"<<vectorToString(dataSitesVec->at(sitesIdx).freqDerived)<<endl;
   // }

   return 0;
}

