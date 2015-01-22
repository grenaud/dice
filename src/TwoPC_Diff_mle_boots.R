library("bbmle")
source("TPCfunc_Diff.R")

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Name of data file
infilename <- args[1]
# Index of run (leave as NA)
i <- args[2]
# Number of bootstraps to run to obtain confidence intervals
bootsnum <- as.numeric(args[3])
# Set this to 1 if you don't want to run a point estimate but just want to append more bootstraps to existing file
# Otherwise, leave as NA
appendbool <- args[4]

if(i == "NA"){ i = NA }
if(appendbool == "NA"){ appendbool = NA }

# Load data table
table <- read.table(infilename,header=TRUE)

# Set lower boundaries for optimization algorithm
elower <- 0.00001
rlower <- 0.00001
tau_Clower <- 0.000001
tau_Alower <- 0.000001

# Set upper boundaries for optimization algorithm
eupper <- 0.1
rupper <- 0.5
tau_Cupper <- 1
tau_Aupper <- 1

# Only relevant if i != NA, otherwise ignore
if(is.na(i)){
    outfilename = gsub("simul_t","estimate_t",infilename)
}
if(!is.na(i)){
    outfilename = gsub("simul_t",paste("estimate_t",i,"_t",sep=""),infilename)
}

# If we are writing to the outfile for the first time...
if( is.na(appendbool) | appendbool != 1){
    # This will be the function that will spit out the minus log-likelihood for the first optimzation
    # DO NOT MOVE IT FROM HERE
    OptimFunc <- function(e,r,tau_C,tau_A){
        result <- -LogFinalTwoP(table,e,r,tau_C,tau_A,TRUE)
        #print(result)
        return(result)
    }

    # We'll run the algorithm 5 times from different starting points to ensure convergence is achieved
    # We'll keep the run with the best log-likelihood
    # Note that these are not bootstrap replicates! We'll do those later in the script
    bestlik <- -Inf
    for(rep in seq(1,5)){ 

        # Set starting values for optimization algorithm - just pick a random number from within the bounds
        estart <- runif(1,elower,eupper)
        rstart <- runif(1,rlower,rupper)
        tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
        tau_Astart <- runif(1,tau_Alower,tau_Aupper)
        
        # Run optimization algorithm
        temppost <- mle2(OptimFunc, method="L-BFGS-B",start=list(e=estart,r=rstart,tau_C=tau_Cstart,tau_A=tau_Astart),lower=list(e=elower,r=rlower,tau_C=tau_Clower,tau_A=tau_Alower),upper=list(e=eupper,r=rupper,tau_C=tau_Cupper,tau_A=tau_Aupper),control=list(maxit=10^6))
        # Obtain -log-likelihood
        temploglik <- -attributes(summary(temppost))$m2logL

        # If current estimate is better than current best estimate, replace for current best estimate
        if( temploglik > bestlik ){
            postestimate <- temppost
            loglik <- temploglik
            bestlik <- temploglik
        }
    }
    print(attributes(summary(postestimate))$coef[,1])
    # Write best estimate to file
    write(c(attributes(summary(postestimate))$coef[,1],loglik),file=outfilename,sep="\t",append=FALSE)
}

# Function to obtain bootstrap samples from data
boots <- function(table){
    total <- sum(table[,4])
    probvec <- table[,4] / total
    sampidx <- sample(seq(1,length(table[,4])),total,replace=TRUE,prob=probvec)
    sampcounts <- table(sampidx)
    sampidx <- as.numeric(names(sampcounts))
    bootstable <- cbind(table[sampidx,seq(1,3)],as.vector(sampcounts))
    colnames(bootstable) <- c("Anc","Der","PanelFreq","Num")
    return(bootstable)
}

# Run bootstrap estimates
sapply(seq(1,bootsnum),function(x){

    # Obtain bootstrap sample
    bootstable <- boots(table)

    # Run optimization on bootstrapped sample 
    OptimBoots <- function(e,r,tau_C,tau_A){
        result <- -LogFinalTwoP(bootstable,e,r,tau_C,tau_A,TRUE)
        return(result)
    }
    estart <- runif(1,elower,eupper)
    rstart <- runif(1,rlower,rupper)
    tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
    tau_Astart <- runif(1,tau_Alower,tau_Aupper)
    postestimate <- mle2(OptimBoots, method="L-BFGS-B",start=list(e=estart,r=rstart,tau_C=tau_Cstart,tau_A=tau_Astart),lower=list(e=elower,r=rlower,tau_C=tau_Clower,tau_A=tau_Alower),upper=list(e=eupper,r=rupper,tau_C=tau_Cupper,tau_A=tau_Aupper),control=list(maxit=10^6))

    # Obtain log-likelihood from boostrap
    loglik <- -attributes(summary(postestimate))$m2logL
    # Write point estimates from bootstrap to file
    write(c(attributes(summary(postestimate))$coef[,1],loglik),file=outfilename,sep="\t",append=TRUE)
})

