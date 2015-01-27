/*
 * libnuc
 * Date: Jan-27-2015 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef libnuc_h
#define libnuc_h

using namespace std;


#include "utils.h"

typedef struct{
    int ancCount;
    int derCount;     
    long double panelFreq;
    int num;
} freqSite;


// q term for binomial sampling from true genotype, to incorporate contamination and error rates
// i = genotype, 2 = homo anc, 1 = het anc/der, 0 = homo derived
// e = error 
// 
inline long double qterm(int i,long double r,long double e,long double y){

    if(i == 2){          //True genotype is homozygous ancestral
        return  (        r*y*(1.0-e) + r*(1.0-y)*e + (1.0-r)*(1.0-e)                     );
    }else{ 
        if(i == 1){      //True genotype is heterozygous
            return (     r*y*(1.0-e) + r*(1.0-y)*e + (1.0-r)*(1.0-e)/2.0 + (1.0-r)*e/2.0 );
        }else{

            if(i == 0){  // True genotype is homozygous derived
                return ( r*y*(1.0-e) + r*(1.0-y)*e                       + (1.0-r)*e     );
            }else{
		cerr<<"Internal error, wrong genotype term in qterm() "<<i<<endl;
		exit(1);
	    }
	}
    }

}

// Binomial probability of sampling ancestral and derived reads given q term
inline long double Pad_given_irey(int a,int d,int i,long double r,long double e,long double y){
    long double qtermA = qterm(i,r,e,y);
    return  ( (long double)(nChoosek(a+d,d)) ) * (powl(qtermA,d))* (1-powl(qtermA,a) )  ; //defined in lib gab
}



// Probability of ancient genotype given anchor population frequency and drift parameters   
inline long double Pgeno_given_ytau(int i,long double y,long double tau_C,long double tau_A){

    if(i == 0){           // Homozygous ancestral
        return (     1 - y*exp(      -1.0*tau_C)         - (0.5*y)*exp(-1.0*tau_A - 1.0*tau_C) + y*(y-0.5)*exp(-1.0*tau_A - 3.0*tau_C) );
    }else{
        if(i == 1){      // Heterozygous
            return (     y*exp(-1.0*tau_A-tau_C)   + y*(1.0-2.0*y)*exp(-1.0*tau_A - 3.0*tau_C)                                         );
        }else{
            if(i == 2){  // Homozygous derived
                return ( y*exp(-1.0*      tau_C)         - (0.5*y)*exp(-1.0*tau_A - 1.0*tau_C) + y*(y-0.5)*exp(-1.0*tau_A - 3.0*tau_C) );
            }else{
                cerr<<"Internal error, wrong genotype term in Pgeno_given_ytau() "<<i<<endl;
                exit(1);
            }
        }
    }
}

// Sum over each of the 3 types of genotypes for two-population method (no admixture)
inline long double Pad_given_reytau(int a,int d,long double r,long double e,long double y,long double freqcont,long double tau_C,long double tau_A){
    double long sumResult=0.0;
    for(int i=0;i<=2;i++){
        sumResult += Pad_given_irey(a,d,i,r,e,freqcont)*Pgeno_given_ytau(i,y,tau_C,tau_A);
    }
    // result <- sum(sapply(c(0,1,2),
    //                   function(i){})
    //            );
    return sumResult;

}

//  Log final posterior for two populations
long double LogFinalTwoP(vector<freqSite> * tableData,long double e,long double r,long double tau_C,long double tau_A,bool contequalanchor){
//     print(c(e,r,tau_C,tau_A))

	// Case where the anchor population is the same as the putative contaminant population (3rd column of data file)
    if(contequalanchor){
        long double sumterm=0.0;
        for(unsigned int indexSite=0;indexSite<tableData->size();indexSite++){
	    //result <- sum(apply(table,1,function(x){
	    //sumterm += log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[3],tau_C,tau_A))*x[4]
	    sumterm += log(Pad_given_reytau(tableData->at(indexSite).ancCount,
					    tableData->at(indexSite).derCount,
					    r,e,tableData->at(indexSite).panelFreq,tableData->at(indexSite).panelFreq,tau_C,tau_A))*tableData->at(indexSite).num;
	}
	return sumterm;
	// }))
    }else{
	// Case where the anchor population is different from the putative contaminant population
	// Anchor population frequencies are in third column of data file
	//Putative contaminant population frequencies are in fourth columnn of data file
	//TODO
        //          result <- sum(apply(table,1,function(x){
        //          sumterm <- log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[4],tau_C,tau_A))*x[5]
        //                   return(sumterm)
        // }))
	return 0;
	
    }

//     print(result)
//     if(is.na(result)){return(-1000000000000000)}
//     else{return(result)}
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
