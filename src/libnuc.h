/*
 * libnuc
 * Date: Jan-27-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef libnuc_h
#define libnuc_h

/* #define DEBUGL1 */
/* #define DEBUGL2 */
/* #define DEBUGL3 */

using namespace std;


#include "utils.h"

typedef struct{
    int ancCount;
    int derCount;     
    long double panelFreqCont;
    long double panelFreqAdmx;
    long double panelFreqAnchor;
    int num;
} freqSite;

//! A method that computes the probability of seeing the derived allele for a given genotype model, to incorporate contamination and error rates
/*!
  \param i : dummy var 0 == homo anc, 1 = het, 2 == homo der
  \param r : cont. rate
  \param e : error. rate
  \param y : frequency of the derived allele for the contaminant
 
*/
inline long double qterm(int i,long double r,long double e,long double y){

#ifdef DEBUGL1
    cout<<"qterm begin\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;
#endif

    long double toreturn;
    if(i == 2){          //True genotype is homozygous derived
	//                   a:d->d          a:a->d           e:d->d
        toreturn=  (        r*y*(1.0-e) + r*(1.0-y)*e + (1.0-r)*(1.0-e)                     );
    }else{ 
        if(i == 1){      //True genotype is heterozygous
	    //                c:d->d         c:a->d          e:d->d                 e:a->d
	    toreturn= (     r*y*(1.0-e) + r*(1.0-y)*e + (1.0-r)*(1.0-e)/2.0 + (1.0-r)*e/2.0 );
        }else{

            if(i == 0){  // True genotype is homozygous ancestral
		//            a:d->d          a:a->d                           e:a->d
		toreturn= ( r*y*(1.0-e) + r*(1.0-y)*e                       + (1.0-r)*e     );
            }else{
		cerr<<"Internal error, wrong genotype term in qterm() "<<i<<endl;
		exit(1);
	    }
	}
    }

#ifdef DEBUGL1
    cout<<"qterm end\t"<<toreturn<<endl;
#endif

    return toreturn;
}

//! A method that computes Binomial probability of sampling ancestral and derived reads given q term
/*!
  \param a : ancestral count
  \param d : derived count
  \param i : dummy var 0 == homo anc, 1 = het, 2 == homo der
  \param r : cont. rate
  \param e : error. rate
  \param y : frequency of the derived allele for the contaminant
 */

//(a,d,i,r,e,freqcont)
//a = 
// Binomial probability of sampling ancestral and derived reads given q term
inline long double pad_given_irey(int a,int d,int i,long double r,long double e,long double y){
#ifdef DEBUGL2
    //cout<<"qterm begin\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;
    cout<<"pad_given_irey begin\t"<< a<<"\t"<<d<<"\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;
#endif


    long double qtermA   = qterm(i,r,e,y);
    long double toreturn = ( (long double)(nChoosek(a+d,d)) ) * (powl(qtermA,d))* (powl(1.0-qtermA,a) )  ; //defined in lib gab

#ifdef DEBUGL2
    cout<<"pad_given_irey res\ta="<< a<<"\t"<<d<<"\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;    
    cout<< (long double)(nChoosek(a+d,d)) <<"\t"<<(powl(qtermA,d)) <<"\t"<< (1.0-powl(qtermA,a)) <<endl;
    cout<<"Comb= "<<( (long double)(nChoosek(a+d,d)) )<<endl;
    cout<<qtermA<<endl;
    cout<<toreturn<<endl;
#endif

    return toreturn;
}



// Probability of ancient genotype given anchor population frequency and drift parameters   
inline long double Pgeno_given_ytau(int i,long double y,long double tau_C,long double tau_A){
    long double toreturn;
#ifdef DEBUGL2
    //cout<<"qterm begin\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;
    cout<<"Pgeno_given_ytau begin\t"<< i<<"\t"<<y<<"\t"<<tau_C<<"\t"<<tau_A<<endl;
#endif

    if(i == 0){           // Homozygous ancestral
        toreturn= (     1 - y*exp(      -1.0*tau_C)         - (0.5*y)*exp(-1.0*tau_A - 1.0*tau_C) + y*(y-0.5)*exp(-1.0*tau_A - 3.0*tau_C) );
    }else{
        if(i == 1){      // Heterozygous
            toreturn= (     y*exp(-1.0*tau_A-tau_C)   + y*(1.0-2.0*y)*exp(-1.0*tau_A - 3.0*tau_C)                                         );
        }else{
            if(i == 2){  // Homozygous derived
                toreturn= ( y*exp(-1.0*      tau_C)         - (0.5*y)*exp(-1.0*tau_A - 1.0*tau_C) + y*(y-0.5)*exp(-1.0*tau_A - 3.0*tau_C) );
            }else{
                cerr<<"Internal error, wrong genotype term in Pgeno_given_ytau() "<<i<<endl;
                exit(1);
            }
        }
    }

#ifdef DEBUGL2
    //cout<<"qterm begin\t"<<i<<"\t"<<r<<"\t"<<e<<"\t"<<y<<endl;
    cout<<"Pgeno_given_ytau begin\t"<< i<<"\t"<<y<<"\t"<<tau_C<<"\t"<<tau_A<<endl;
#endif
    return toreturn;

}

// Sum over each of the 3 types of genotypes for two-population method (no admixture)
inline long double pad_given_reytau(int a,int d,long double r,long double e,long double y,long double freqcont,long double tau_C,long double tau_A){

    double long sumResult=0.0;
    for(int i=0;i<=2;i++){
        sumResult += pad_given_irey(a,d,i,r,e,freqcont)*Pgeno_given_ytau(i,y,tau_C,tau_A);
    }

    return sumResult;
}

//Sum over each of the 3 types of genotypes for three-population method (incorporating admixture)
inline long double pad_given_reytau_dadi_threeP(int a,int d,long double r,long double e,long double y,long double z,long double freqcont, vector< vector<long double> * > * dadiTable,long double numhumy,long double numhumz){
    /* cout<<"pad_given_reytau_dadi_threeP1\t"<<a<<"\t"<<d<<endl; */

    long double  ycoord     = roundl( y * numhumy );
    long double  zcoord     = roundl( z * numhumz );
    long double  finalcoord = ycoord*(numhumy+1) + zcoord ; //TODO check for 
    double long sumResult=0.0;
    for(int i=0;i<=2;i++){
	/* cout<<"pad_given_reytau_dadi_threeP1\ts="<<dadiTable->size()<<"\ts="<<dadiTable->at(finalcoord)->size()<<"\t"<<ycoord<<"\t"<<zcoord<<"\t"<<finalcoord<<"\t"<<i<<endl; */
	/* cout<<"P "<<pad_given_irey(a,d,i,r,e,freqcont)<<endl; */
	/* cout<<"d "<<dadiTable->at(finalcoord)->at(i)<<endl; */
        sumResult += pad_given_irey(a,d,i,r,e,freqcont)*dadiTable->at(finalcoord)->at(i);
	//pad_given_irey(a,d,i,r,e,freqcont)*Pgeno_given_ytau(i,y,tau_C,tau_A);	             
    }

    /* cout<<"pad_given_reytau_dadi_threeP2\t"<< */
    /* 	a<<"\t"<< */
    /* 	d<<"\t"<< */
    /* 	r<<"\t"<< */
    /* 	e<<"\t"<< */
    /* 	y<<"\t"<< */
    /* 	z<<"\t"<< */
    /* 	freqcont<<"\t"<< */
    /* 	numhumy<<"\t"<< */
    /* 	numhumz<<endl; */
    /* cout<<sumResult<<endl; */
    //    result <- sum(sapply(c(0,1,2),function(i){Pad_given_irey(a,d,i,r,e,freqcont)*DadiTable[finalcoord,(i+1)]}))
    return sumResult ;
}

inline string ld2string(long double i){
    stringstream s;
    s << std::setprecision(10)<<i;
    return s.str();
}


inline vector< vector<long double>  * > * getDadiTableThreeP(long double tau_C,long double tau_A,long double admixrate,long double admixtime,long double innerdriftY,long double innerdriftZ,long double nC,long double nB,long double nA,const string & cwdProg){
    string cmd = "python  "+cwdProg+"/Dadi_three_pop_admix.py -c "+ld2string(tau_C)+" -a "+ld2string(tau_A)+" -x "+ld2string(admixrate)+" -t "+ld2string(admixtime)+" -y "+ld2string(innerdriftY)+" -z "+ld2string(innerdriftZ)+" -m "+ld2string(nC)+" -n "+ld2string(nA)+" -b "+ld2string(nB)+" ";
    //cerr<<cmd<<endl;
    string values = runCmdAndCaptureSTDOUTandSTDERR(cmd);
    
    vector<string>       numberS   = allTokens(values,'\n');
    vector<long double>  numbersLD;

    if(numberS[0] == "Traceback (most recent call last):"){
	cerr<<"An error has occurred, do you have all the python dependencies ? ex: dadi and numpy"<<endl<<"error:"<<endl<<values<<endl;
	exit(1);
    }

    for(unsigned int i=0;i<numberS.size();i++){

	if(!numberS[i].empty()){//avoid empty lines
	    long double numFound = destringify<long double>(numberS[i]) ;
	    //cout<<"#"<<numberS[i]<<"#"<<numFound<<"#"<<endl;
	    numbersLD.push_back( numFound );
	}
    }

    /* for(unsigned int i=0;i<numbersLD.size();i++){ */
    /* 	cout<<"#"<<(i+1)<<"#"<<numbersLD[i]<<"#"<<endl; */
    /* } */

    int    numberOfCols  = nA+1;
    double numberOfRows  = double(numbersLD.size())/(double(numberOfCols));

    if(floor(numberOfCols) != numberOfCols){
	cerr<<"Internal error, the number of numberOfCols in the dadi output is not a multiple of nA + 1 = "<<(nA+1)<<" in command "<<cmd<<endl;
    }
    
    //cout<<numberOfCols<<"\t"<<nA<<"\t"<<nB<<"\t"<<numberS.size()<<endl;
    vector< vector<long double>  * >  * toReturn = new vector< vector<long double> * >( int(numberOfRows) );
    
     for(int i=0;i<int(numberOfRows);i++){ 
	 toReturn->at(i) =  new vector<long double> ( numberOfCols );
     } 
     //cout<<"nr "<<numberOfRows<<" nc "<<numberOfCols<<endl;

     int numbersLDi=0;

     //copy to matrix
     for(int c=0;c<int(numberOfCols);c++){ 	     
	 for(int r=0;r<int(numberOfRows);r++){ 
	     toReturn->at(r)->at(c) = numbersLD[numbersLDi++];

	 }
     }

     //sum per col
     for(int r=0;r<int(numberOfRows);r++){ 
	 long double sumAllNum  = 0.0;

	 for(int c=0;c<int(numberOfCols);c++){ 	     
	     sumAllNum += toReturn->at(r)->at(c);
	     //cout<<"toReturn["<<r<<","<<c<<"] "<<toReturn->at(r)->at(c)<<" = numbersLD["<<numbersLDi<<"]"<<" = "<<numbersLD[numbersLDi]<<endl;
	 }

	 for(int c=0;c<int(numberOfCols);c++){ 	     
	     toReturn->at(r)->at(c) = toReturn->at(r)->at(c) / sumAllNum;
	     //sumAllNum += toReturn->at(r)->at(c);
	     //cout<<"toReturn["<<r<<","<<c<<"] "<<toReturn->at(r)->at(c)<<endl;
	 }


	 //cout<<endl;
     }

     

    return toReturn;
}

//  Log final posterior for two populations
inline long double LogFinalTwoP(vector<freqSite> * tableData,long double e,long double r,long double tau_C,long double tau_A,bool contequalanchor){
//     print(c(e,r,tau_C,tau_A))

	// Case where the anchor population is the same as the putative contaminant population (3rd column of data file)
    if(contequalanchor){
        long double sumterm=0.0;
        for(unsigned int indexSite=0;indexSite<tableData->size();indexSite++){

	    //result <- sum(apply(table,1,function(x){
	    //sumterm += log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[3],tau_C,tau_A))*x[4]
	    long double toaddToSum = log(pad_given_reytau(tableData->at(indexSite).ancCount,
							  tableData->at(indexSite).derCount,
							  r,e,
							  tableData->at(indexSite).panelFreqCont,
							  tableData->at(indexSite).panelFreqCont,
							  tau_C,tau_A))*tableData->at(indexSite).num;
	    //cout<<tableData->at(indexSite).ancCount<<"\t"<<tableData->at(indexSite).derCount<<"\t"<<tableData->at(indexSite).panelFreqCont<<"\t"<<tableData->at(indexSite).num<<"\t"<<toaddToSum<<endl;
	    sumterm+=toaddToSum;
	}
	return sumterm;
	// }))
    }else{
        long double sumterm=0.0;

	for(unsigned int indexSite=0;indexSite<tableData->size();indexSite++){

	    //result <- sum(apply(table,1,function(x){
	    //sumterm += log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[3],tau_C,tau_A))*x[4]
	    long double toaddToSum = log(pad_given_reytau(tableData->at(indexSite).ancCount,
							  tableData->at(indexSite).derCount,
							  r,e,
							  tableData->at(indexSite).panelFreqAnchor,
							  tableData->at(indexSite).panelFreqCont,
							  tau_C,
							  tau_A))*tableData->at(indexSite).num;
	    //cout<<tableData->at(indexSite).ancCount<<"\t"<<tableData->at(indexSite).derCount<<"\t"<<tableData->at(indexSite).panelFreqCont<<"\t"<<tableData->at(indexSite).num<<"\t"<<toaddToSum<<endl;
	    sumterm+=toaddToSum;
	}
	return sumterm;

	
    }

}



//  Log final posterior for three populations
inline long double LogFinalThreeP(vector<freqSite> * tableData,long double e,long double r,long double tau_C,long double tau_A,long double admixrate,long double admixtime,long double innerdriftY,long double innerdriftZ,long double nC,long double nB,const string & cwdProg,bool contequalanchor){
    //     print(c(e,r,tau_C,tau_A))
    //cout<<cwdProg<<"\t"<<contequalanchor<<endl;
    //exit(1);
    if( (e < 0)     || 
	(r < 0)     || 
	(tau_C < 0) || 
	(tau_A < 0) || 
	(admixrate < 0) || admixtime < 0){
        return(-1000000000000000);
    }else{
        long double nA = 2.0;
	vector< vector<long double> * > * dadiTable= getDadiTableThreeP(tau_C,
									tau_A,
									admixrate,
									admixtime,
									innerdriftY,
									innerdriftZ,
									nC,
									nB,
									nA,
									cwdProg);
	//cout<<"after getDadiTableThreeP"<<endl;
	long double sumterm=0.0;
	if(contequalanchor){


	    for(unsigned int indexSite=0;indexSite<tableData->size();indexSite++){
		//cout<<"LogFinalThreeP1 "<<indexSite<<endl;
		long double toaddToSum=0;

		toaddToSum = log(pad_given_reytau_dadi_threeP(tableData->at(indexSite).ancCount,
							      tableData->at(indexSite).derCount,
							      r,
							      e,
							      tableData->at(indexSite).panelFreqCont,
							      tableData->at(indexSite).panelFreqAnchor ,
							      tableData->at(indexSite).panelFreqCont,
							      dadiTable,
							      nC,
							      nB)
					     ) * tableData->at(indexSite).num;
		/* cout<<tableData->at(indexSite).ancCount<<"\t"<<tableData->at(indexSite).derCount<<"\t" */
		/*     <<"\t"<<r */
		/*     <<"\t"<<e */
		/*     <<"\t"<<tableData->at(indexSite).panelFreqCont */
		/*     <<"\t"<<tableData->at(indexSite).panelFreqAnchor  */
		/*     <<"\t"<<tableData->at(indexSite).panelFreqCont  */
		/*     <<"\t"<<nC */
		/*     <<"\t"<<nB<<"\t"<<tableData->at(indexSite).num<<"\t"<<toaddToSum<<endl; */
		//cout<<"LogFinalThreeP2 "<<indexSite<<endl;
		sumterm+=toaddToSum;
	    }
	    //cout<<"LogFinalThreeP done "<<endl;
	
	    
            /* result <- sum(apply(table,1,function(x){ */
            /*     sumterm <- log(Pad_given_reytau_dadi_threeP(x[1],x[2],r,e,x[3],x[4],x[3],DadiTable,nC,nB))*x[5] */
            /*     return(sumterm) */
            /* })) */
        }


	// Case where anchor populations are different from the putative contaminant population
	// Frequencies of anchor population 1 would be in the third column of data file
        // Frequencies of anchor population 2 would be in the fourth column of data file
        // Contaminant population frequencies would be in the fifth column of data file        
        else{
	    //TODO
	    //exit(1);

	    for(unsigned int indexSite=0;indexSite<tableData->size();indexSite++){
		//cout<<"LogFinalThreeP1 "<<indexSite<<endl;
		long double toaddToSum=0;

		//int a,int d,long double r,long double e,long double y,long double z,long double freqcont, vector< vector<long double> * > * dadiTable,long double numhumy,long double numhumz
		//int a,int d,long double r,long double e,long double y,long double z,long double freqcont, vector< vector<long double> * > * dadiTable,long double numhumy,long double numhumz

		toaddToSum = log(pad_given_reytau_dadi_threeP(tableData->at(indexSite).ancCount,
							      tableData->at(indexSite).derCount,
							      r,
							      e,
							      tableData->at(indexSite).panelFreqAdmx,
							      tableData->at(indexSite).panelFreqAnchor ,
							      tableData->at(indexSite).panelFreqCont,
							      dadiTable,
							      nC,
							      nB)
					     ) * tableData->at(indexSite).num;
		//cout<<"LogFinalThreeP2 "<<indexSite<<endl;
		sumterm+=toaddToSum;
	    }
	}


	//destroy table	
	for(unsigned int i=0;i<dadiTable->size();i++){
	
	    /* for(unsigned int j=0;j<dadiTable->at(i)->size();j++){ */
	    /* 	delete dadiTable->at(i)->at(j); */
	    /* 	//     cout<<dadiTable->at(i)->at(j).x <<","<<dadiTable->at(i)->at(j).y <<"\t"; */
	    /* } */
	    // cout<<endl;
	    delete dadiTable->at(i);
	}
	
	delete dadiTable;
	return sumterm;


    }

}


/* class libnuc{ */
/*  private: */

/*  public: */
/*     libnuc(); */
/*     libnuc(const libnuc & other); */
/*     ~libnuc(); */
/*     libnuc & operator= (const libnuc & other); */

/* }; */
#endif
