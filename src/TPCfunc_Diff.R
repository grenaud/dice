library(ape)
library(caper)
library(geiger)
library(parallel)

# Probability of ancient genotype given anchor population frequency and drift parameters   
Pgeno_given_ytau <- function(i,y,tau_C,tau_A){
    # Homozygous ancestral
    if(i == 0){
        return( 1 - y*exp(-tau_C) - (1/2)*y*exp(-tau_A - tau_C) + y*(y-1/2)*exp(-tau_A - 3*tau_C) )
    }
    # Heterozygous
    else if(i == 1){
        return( y*exp(-tau_A-tau_C) + y*(1-2*y)*exp(-tau_A-3*tau_C) )
    }
    # Homozygous derived
    else if(i == 2){
        return( y*exp(-tau_C) - (1/2)*y*exp(-tau_A-tau_C) + y*(y-1/2)*exp(-tau_A-3*tau_C) )
    }
}

# q term for binomial sampling from true genotype, to incorporate contamination and error rates
qterm <- function(i,r,e,y){
    # True genotype is homozygous ancestral
    if(i == 2){
        result <- r*y*(1-e) + r*(1-y)*e + (1-r)*(1-e)
    }
    # True genotype is heterozygous
    else if(i == 1){
        result <- r*y*(1-e) + r*(1-y)*e + (1-r)*(1-e)/2 + (1-r)*e/2
    }
    # True genotype is homozygous derived
    else if(i == 0){
        result <- r*y*(1-e) + r*(1-y)*e + (1-r)*e
    }
    return(result)
}

# Binomial probability of sampling ancestral and derived reads given q term
Pad_given_irey <- function(a,d,i,r,e,y){
    qtermA <- qterm(i,r,e,y)
    result <- choose(a+d,d)*(qtermA^d)*((1-qtermA)^a)
    return(result)
}

# Sum over each of the 3 types of genotypes for two-population method (no admixture)
Pad_given_reytau <- function(a,d,r,e,y,freqcont,tau_C,tau_A){
    result <- sum(sapply(c(0,1,2),function(i){Pad_given_irey(a,d,i,r,e,freqcont)*Pgeno_given_ytau(i,y,tau_C,tau_A)}))
    return(result)
}

#Pad_given_reytau_dadi_twoP <- function(a,d,r,e,y,freqcont,DadiTable,numhumy){
#    ycoord <- round(y * numhumy, 0)
#    finalcoord <- ycoord + 1
#    result <- sum(sapply(c(0,1,2),function(i){Pad_given_irey(a,d,i,r,e,freqcont)*DadiTable[finalcoord,(i+1)]}))
#    return(result)
#}

# Sum over each of the 3 types of genotypes for three-population method (incorporating admixture)
Pad_given_reytau_dadi_threeP <- function(a,d,r,e,y,z,freqcont,DadiTable,numhumy,numhumz){
    ycoord <- round(y * numhumy, 0)
    zcoord <- round(z * numhumz, 0)
    finalcoord <- ycoord*(numhumy+1) + zcoord + 1
    result <- sum(sapply(c(0,1,2),function(i){Pad_given_irey(a,d,i,r,e,freqcont)*DadiTable[finalcoord,(i+1)]}))
    return(result)
}


# Uniform log prior for contamination rate
#LogP_r <- function(r,rmin,rmax){
#    if(r > rmin & r < rmax){ result <- 0 }
#    else{ result <- -1000000 }
#    return(result)
#}

# Uniform log prior for error rate
#LogP_e <- function(e,emin,emax){
#    if(e > emin & e < emax){ result <- 0 }
#    else{ result <- -1000000 }
#    return(result)
#}

# Uniform log prior for tau_A
#LogP_tau_A <- function(tau_A,tau_Amin,tau_Amax){
#    if(tau_A > tau_Amin & tau_A < tau_Amax){ result <- 0 }
#    else{ result <- -1000000 }
#    return(result)
#}

# Uniform log prior for tau_C
#LogP_tau_C <- function(tau_C,tau_Cmin,tau_Cmax){
#    if(tau_C > tau_Cmin & tau_C < tau_Cmax){ result <- 0 }
#    else{ result <- -1000000 }
#    return(result)
#}


#GetDadiTableTwoP <- function(tau_C,tau_A,admixrate,admixtime,nC,nA){
#    daditable <- system(paste("python ","Dadi_two_pop_admix.py -c ",tau_C," -a ",tau_A," -x ",admixrate," -t ",admixtime," -m ",nC," -n ",nA,sep=""),intern=T)
#    daditable <- as.numeric(daditable)
#    daditable <- matrix(daditable,ncol=(nA+1),byrow=TRUE)
#    daditable <- t(apply(daditable,1,function(x) { return( x / sum(x)) }))
#    return(daditable)
#}

# Call python script to obtain dadi table for three-population method
GetDadiTableThreeP <- function(tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,nA){
    daditable <- system(paste("python ","Dadi_three_pop_admix.py -c ",tau_C," -a ",tau_A," -x ",admixrate," -t ",admixtime," -y ",innerdriftY," -z ",innerdriftZ," -m ",nC," -n ",nA," -b ",nB,sep=""),intern=T)
    daditable <- as.numeric(daditable)
    daditable <- matrix(daditable,ncol=(nA+1),byrow=FALSE)
    daditable <- t(apply(daditable,1,function(x) { return( x / sum(x)) }))
    return(daditable)
}


# Log final posterior for two populations
LogFinalTwoP <- function(table,e,r,tau_C,tau_A,contequalanchor){
    print(c(e,r,tau_C,tau_A))

    # Case where the anchor population is the same as the putative contaminant population (3rd column of data file)
    if(contequalanchor == "TRUE"){
        result <- sum(apply(table,1,function(x){
            sumterm <- log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[3],tau_C,tau_A))*x[4]
            #sumterm <- sumterm + log(P_e(e)) + log(P_r(r)) + log(P_tau_C(tau_C)) + log(P_tau_A(tau_A))
            return(sumterm)
        }))
    }

    # Case where the anchor population is different from the putative contaminant population
    # Anchor population frequencies are in third column of data file
    # Putative contaminant population frequencies are in fourth columnn of data file
    else{
        result <- sum(apply(table,1,function(x){
            sumterm <- log(Pad_given_reytau(x[1],x[2],r,e,x[3],x[4],tau_C,tau_A))*x[5]
            return(sumterm)
        }))
    }
    
    print(result)
    if(is.na(result)){return(-1000000000000000)}
    else{return(result)}
}

# Log final posterior - dadi two pop (with admixture)
#LogFinalTwoPDadi <- function(table,e,r,tau_C,tau_A,admixrate,admixtime,nC,contequalanchor){
#
#    print(c(e,r,tau_C,tau_A,admixrate,admixtime))
#
#    if(e < 0 | r < 0 | tau_C < 0 | tau_A < 0 | admixrate < 0 | admixtime < 0){
#        return(-1000000000000000)
#    }
#    else{
#        # Create dadi table
#        nA = 2    
#        DadiTable <- GetDadiTableTwoP(tau_C,tau_A,admixrate,admixtime,nC,nA)
#    
#        if(contequalanchor == "TRUE"){
#            result <- sum(apply(table,1,function(x){
#                sumterm <- log(Pad_given_reytau_dadi_twoP(x[1],x[2],r,e,x[3],x[3],DadiTable,nC))*x[4]
#                return(sumterm)
#            }))
#        }
#        else{
#            result <- sum(apply(table,1,function(x){
#                sumterm <- log(Pad_given_reytau_dadi_twoP(x[1],x[2],r,e,x[3],x[4],DadiTable,nC))*x[5]
#                return(sumterm)
#            }))
#        }
#        
#        print(result)
#        if(is.na(result)){return(-1000000000000000)}
#        else{return(result)}
#    }
#}


# Log final posterior - dadi three pop (with admixture)
LogFinalThreePDadi <- function(table,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,contequalanchor){

    print(c(e,r,tau_C,tau_A,admixrate,admixtime))

    # If we are exploring a bad region of parameter space, return a very low log-likelihood
    if(e < 0 | r < 0 | tau_C < 0 | tau_A < 0 | admixrate < 0 | admixtime < 0){
        return(-1000000000000000)
    }
    else{
        # There are two individuals for the ancient sample (i.e. the sample is diploid)
        nA = 2
        # Create dadi table
        DadiTable <- GetDadiTableThreeP(tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,nA)
        
        # Case where one of the anchor populations is the same as the putative contaminant population (3rd column of data file)
        # The frequency of the other anchor population would be in the fourth column of data file
        if(contequalanchor == "TRUE"){
            result <- sum(apply(table,1,function(x){
                sumterm <- log(Pad_given_reytau_dadi_threeP(x[1],x[2],r,e,x[3],x[4],x[3],DadiTable,nC,nB))*x[5]
                return(sumterm)
            }))
        }
        # Case where anchor popuations are different from the putative contaminant population
        # Frequencies of anchor population 1 would be in the third column of data file
        # Frequencies of anchor population 2 would be in the fourth column of data file
        # Contaminant population frequencies would be in the fifth column of data file        
        else{
            result <- sum(apply(table,1,function(x){
                sumterm <- log(Pad_given_reytau_dadi_threeP(x[1],x[2],r,e,x[3],x[4],x[5],DadiTable,nC,nB))*x[6]
                return(sumterm)
            }))
        }
        
        print(result)
        if(is.na(result)){return(-1000000000000000)}
        else{return(result)}
    }
}
