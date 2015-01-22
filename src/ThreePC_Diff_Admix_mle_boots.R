library("bbmle")
source("TPCfunc_Diff.R")

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Name of file containing data
infilename <- args[1]
# ID number (leave as NA)
i <- args[2]
# Number of bootstraps to perform for calculating confidence intervals
bootsnum <- as.numeric(args[3])
# Set this to 1 if you don't want to run a point estimate but just want to append more bootstraps to existing file
# Otherwise, leave as NA
appendbool <- args[4]

if(i == "NA"){ i = NA }
if(appendbool == "NA"){ appendbool = NA }

DadiTable <- NA

# These are the drift times specific to the two most closely related (i.e. modern human) samples
# The three-pop admixture inference method assumes we know this from some other study
# Calculation of these drift times would not require an archaic human so it's ok to have as given
innerdriftY <- 0.16
innerdriftZ <- 0.16
# Number of individuals from the two samples
nC <- 20
nB <- 20

# Load data table
table <- read.table(infilename,header=TRUE)

# Set lower boundaries for optimization algorithm
elower <- 0.00001
rlower <- 0.00001
tau_Clower <- 0.000001
tau_Alower <- 0.000001
admixratelower <- 0.000001
admixtimelower <- 0.05

# Set upper boundaries for optimization algorithm
eupper <- 0.1
rupper <- 0.5
tau_Cupper <- 1
tau_Aupper <- 1
admixrateupper <- 0.5
admixtimeupper <- 0.11

# Only relevant if i != NA, otherwise ignore
if(is.na(i)){
    outfilename = gsub("simul_t","estimate_dadi_t",infilename)
}
if(!is.na(i)){
    outfilename = gsub("simul_t",paste("estimate_dadi_",i,"_t",sep=""),infilename)
}

# OPTIMIZATION
# Note that, unlike the two-populaion method, the optimzation here is done in two steps
# because doing it in one step led to inaccurate estimates
# This is because here there are 6 parameters insted of 4, which makes the likelihood surface harder to explore

# If we are writing to the outfile for the first time...
if( is.na(appendbool) | appendbool != 1){

    # This will be the function that will spit out the minus log-likelihood for the first optimzation
    # DO NOT MOVE IT FROM HERE
    OptimFunc <- function(e,r,tau_C,tau_A,admixrate,admixtime){
        result <- -LogFinalThreePDadi(table,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,TRUE)
        #print(result)
        return(result)
    }
    
    # Set starting values for first optimization algorithm - just pick a random number from within the bounds
    estart <- runif(1,elower,eupper)
    rstart <- runif(1,rlower,rupper)
    tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
    tau_Astart <- runif(1,tau_Alower,tau_Aupper)
    admixratestart <- runif(1,admixratelower,admixrateupper)
    admixtimestart <- runif(1,admixtimelower,admixtimeupper)
    
    # Run first opitimzation (L-BFGS-B algorithm)
    postestimate <- mle2(OptimFunc, method="L-BFGS-B",start=list(e=estart,r=rstart,tau_C=tau_Cstart,tau_A=tau_Astart,admixrate=admixratestart,admixtime=admixtimestart),lower=list(e=elower,r=rlower,tau_C=tau_Clower,tau_A=tau_Alower,admixrate=admixratelower,admixtime=admixtimelower),upper=list(e=eupper,r=rupper,tau_C=tau_Cupper,tau_A=tau_Aupper,admixrate=admixrateupper,admixtime=admixtimeupper),control=list(maxit=10^6))

    # Obtain the point estiamtes
    pointest <- attributes(summary(postestimate))$coef[,1]

    # Fix the error and contamination parameters to be those obtained from the first optimization
    e <- pointest[1]
    r <- pointest[2]

    # The rest of the parameters will be used as starting parameters for the second optimization
    tau_Cstart <- pointest[3]
    tau_Astart <- pointest[4]
    admixratestart <- pointest[5]
    admixtimestart <- pointest[6]

    # This will be the function that will spit out the -log-likelihood for the second optimization
    # DO NOT MOVE IT FROM HERE
    OptimFuncB <- function(tau_C,tau_A,admixrate,admixtime){
        result <- -LogFinalThreePDadi(table,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,TRUE)
        return(result)
    }

    # Run second optimization
    postestimateB <- mle2(OptimFuncB, method="L-BFGS-B",start=list(tau_C=tau_Cstart,tau_A=tau_Astart,admixrate=admixratestart,admixtime=admixtimestart),lower=list(tau_C=tau_Clower,tau_A=tau_Alower,admixrate=admixratelower,admixtime=admixtimelower),upper=list(tau_C=tau_Cupper,tau_A=tau_Aupper,admixrate=admixrateupper,admixtime=admixtimeupper),control=list(maxit=10^6))

    # Obtain point estimates
    pointest <- attributes(summary(postestimateB))$coef[,1]
    # Merge with previously obtained error and contamination estimates
    pointest <- c(e,r,pointest[1],pointest[2],pointest[3],pointest[4])
    # Obtain log-likelihood
    loglik <- -attributes(summary(postestimateB))$m2logL
    # Write to file
    write(paste(c(pointest,loglik),collapse="\t"),file=outfilename,append=FALSE)
    print(pointest)
}

# Function to obtain bootstrap samples from data
lastrow <- dim(table)[2]
boots <- function(table){
    total <- sum(table[,lastrow])
    probvec <- table[,lastrow] / total
    sampidx <- sample(seq(1,length(table[,lastrow])),total,replace=TRUE,prob=probvec)
    sampcounts <- table(sampidx)
    sampidx <- as.numeric(names(sampcounts))
    bootstable <- cbind(table[sampidx,seq(1,(lastrow-1))],as.vector(sampcounts))
    if(lastrow == 5){ colnames(bootstable) <- c("Anc","Der","PanelFreqA","PanelFreqB","Num") }
    if(lastrow == 6){ colnames(bootstable) <- c("Anc","Der","PanelFreqA","PanelFreqB","ContPanelFreq","Num") }
    return(bootstable)
}

# Run bootstrap estimates
sapply(seq(1,bootsnum),function(x){
    # Obtain bootstrap sample
    bootstable <- boots(table)

    # Run first optimization on bootstrapped sample
    OptimBoots <- function(e,r,tau_C,tau_A,admixrate,admixtime){
        result <- -LogFinalThreePDadi(bootstable,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,TRUE)
        return(result)
    }
    estart <- runif(1,elower,eupper)
    rstart <- runif(1,rlower,rupper)
    tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
    tau_Astart <- runif(1,tau_Alower,tau_Aupper)
    admixratestart <- runif(1,admixratelower,admixrateupper)
    admixtimestart <- runif(1,admixtimelower,admixtimeupper)
    postestimate <- mle2(OptimBoots, method="L-BFGS-B",start=list(e=estart,r=rstart,tau_C=tau_Cstart,tau_A=tau_Astart,admixrate=admixratestart,admixtime=admixtimestart),lower=list(e=elower,r=rlower,tau_C=tau_Clower,tau_A=tau_Alower,admixrate=admixratelower,admixtime=admixtimelower),upper=list(e=eupper,r=rupper,tau_C=tau_Cupper,tau_A=tau_Aupper,admixrate=admixrateupper,admixtime=admixtimeupper),control=list(maxit=10^6))    
    pointest <- attributes(summary(postestimate))$coef[,1]

    # Run second optimization on bootstrapped sample
    e <- pointest[1]
    r <- pointest[2]
    tau_Cstart <- pointest[3]
    tau_Astart <- pointest[4]
    admixratestart <- pointest[5]
    admixtimestart <- pointest[6]
    OptimBootsB <- function(tau_C,tau_A,admixrate,admixtime){
        result <- -LogFinalThreePDadi(bootstable,e,r,tau_C,tau_A,admixrate,admixtime,innerdriftY,innerdriftZ,nC,nB,TRUE)
        return(result)
    }
    postestimateB <- mle2(OptimBootsB, method="L-BFGS-B",start=list(tau_C=tau_Cstart,tau_A=tau_Astart,admixrate=admixratestart,admixtime=admixtimestart),lower=list(tau_C=tau_Clower,tau_A=tau_Alower,admixrate=admixratelower,admixtime=admixtimelower),upper=list(tau_C=tau_Cupper,tau_A=tau_Aupper,admixrate=admixrateupper,admixtime=admixtimeupper),control=list(maxit=10^6))

    # Obtain point estimates from bootstrap
    pointest <- attributes(summary(postestimateB))$coef[,1]
    pointest <- c(e,r,pointest[1],pointest[2],pointest[3],pointest[4])
    # Obtain log-likelihood from bootstrap
    loglik <- -attributes(summary(postestimateB))$m2logL
    # Write to file
    write(paste(c(pointest,loglik),collapse="\t"),file=outfilename,append=TRUE)
    print(pointest)
})
     
