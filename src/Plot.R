contprior <- c(0.00001,0.5)
errorprior <- c(0.00001,0.1)
driftAprior <- c(0.00001,1)
driftBprior <- c(0.00001,1)

add.error.bars <- function(X,Y,SE,w,col=1){
    print(X)
    print(Y)
    print(SE)
    X0 = X; Y0 = (Y-SE); X1 =X; Y1 = (Y+SE);
    print(Y0)
    arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col);
}


add.conf.int <- function(X,Y,low,high,w,col=1){
    for(i in seq(1,length(low))){
        if( (high[i] - low[i]) < 0.001){
            high[i] <- high[i] + 0.001
            low[i] <- low[i] - 0.001
        }
    }
    print(high-low)
    X0 = X; Y0 = low; X1 =X; Y1 = high;
    arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col);
}





par(mfrow=c(5,4))
sapply(c(0.05,0.1,0.25), function(cont){
    sapply(c(0.5,0.25,0.1), function(split){
        sapply(c(0.001),function(error){
            estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff/estimate_t10_s",split,"_a0_n10000_m30_c",cont,"_e",error,".txt",sep=""))
            contest <- estimatetab[2]
            errorest <- estimatetab[1]
            driftAest <- estimatetab[3]
            driftBest <- estimatetab[4]
            plot(cont,error,xlim=contprior,ylim=errorprior,log="y",main=paste("cont=",cont,", error=",error,", driftC=",split,", driftA=",split,sep=""),cex.main=0.9)
            points(contest,errorest,pch=4,col="red")
            plot(split,split,xlim=driftAprior,ylim=driftBprior)
            points(driftAest,driftBest,pch=4,col="red")
        })
    })
})




# Diffusion - Bootstrap confidence intervals plot
errorpar <- c(0.001); admix <- 0; meancov <- 30; humsamp <- ""
#errorpar <- c(0.1); admix <- 0; meancov <- 30; humsamp <- ""
#errorpar <- c(0.001); admix <- 0.05; meancov <- 30; humsamp <- ""
#errorpar <- c(0.001); admix <- 0; meancov <- 3; humsamp <- ""
#errorpar <- c(0.001); admix <- 0; meancov <- 30; humsamp <- "_lowhum"
bootsnum <- 20
par(mfrow=c(3,3))
sapply(c(0.05,0.1,0.25), function(cont){
    sapply(c(0.5,0.25,0.1), function(split){
        sapply(errorpar,function(error){
            if(admix == 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff/estimate_t10_s",split,"_a",admix,"_n10000_m",meancov,"_c",cont,"_e",error,humsamp,".txt",sep=""))
            }
            if(admix > 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff/estimate_t10_s",split,"_a",admix,"_b0.0375_n10000_m",meancov,"_c",cont,"_e",error,humsamp,".txt",sep=""))
            }
            allest <- unlist(estimatetab[1,])
            allboots <- estimatetab[seq(2,(bootsnum+1)),]
            lowboots <- apply(allboots,2, function(x) { quantile(x,0.025) } )
            highboots <- apply(allboots,2, function(x) { quantile(x,0.975) } )
            numpar <- seq(1,4)
            plot(numpar,allest[numpar],main=paste("cont = ",cont,", error = ",error,", drift C = ",split,", drift A = ",split," admix rate = ",admix,sep=""),xaxt="n",ylim=c(0,1.1),ylab="Parameter value",xlab="",cex=2,cex.lab=1.5)
            add.conf.int(numpar,allest[numpar],lowboots[numpar],highboots[numpar],0.1,col="black")
            legend("topleft",c("True","Estimated"),col=c("red","black"),pch=c(4,1),cex=1.5)
            axis(1,at=numpar,label=c("error","contamination","drift C", "drift A"),cex.axis=1.5)
            points(1,error,pch=4,col="red",cex=3)
            points(2,cont,pch=4,col="red",cex=3)
            points(3,split,pch=4,col="red",cex=3)
            points(4,split,pch=4,col="red",cex=3)
        })
    })
})



# Diffusion - Error bar plot
errorpar <- c(0.001); admix <- 0; meancov <- 30; humsamp <- ""
#errorpar <- c(0.1); admix <- 0; meancov <- 30; humsamp <- ""
#errorpar <- c(0.001); admix <- 0.05; meancov <- 30; humsamp <- ""
#errorpar <- c(0.001); admix <- 0; meancov <- 3; humsamp <- ""
#errorpar <- c(0.001); admix <- 0; meancov <- 30; humsamp <- "_lowhum"

par(mfrow=c(3,3))
sapply(c(0.05,0.1,0.25), function(cont){
    sapply(c(0.5,0.25,0.1), function(split){
        sapply(errorpar,function(error){
            if(admix == 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff/estimate_t10_s",split,"_a",admix,"_n10000_m",meancov,"_c",cont,"_e",error,humsamp,".txt",sep=""),nrows=2)
            }
            if(admix > 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff/estimate_t10_s",split,"_a",admix,"_b0.0375_n10000_m",meancov,"_c",cont,"_e",error,humsamp,".txt",sep=""),nrows=2)
            }
            allest <- unlist(estimatetab[1,])
            allstderr <- 15*unlist(estimatetab[2,])
            numpar <- seq(1,4)
            plot(numpar,allest,main=paste("cont = ",cont,", error = ",error,", drift C = ",split,", drift A = ",split," admix rate = ",admix,sep=""),xaxt="n",ylim=c(0,1.1),ylab="Parameter value",xlab="",cex=2,cex.lab=1.5)
            #print(allstderr)
            #add.error.bars(numpar,allest,allstderr,0.1,col="black")
            legend("topleft",c("True","Estimated"),col=c("red","black"),pch=c(4,1),cex=1.5)
            axis(1,at=numpar,label=c("error","contamination","drift C", "drift A"),cex.axis=1.5)
            points(1,error,pch=4,col="red",cex=3)
            points(2,cont,pch=4,col="red",cex=3)
            points(3,split,pch=4,col="red",cex=3)
            points(4,split,pch=4,col="red",cex=3)
        })
    })
})


# Diffusion Three Pops - Bootstrap confidence intervals plot

errorpar <- c(0.001)
admix <- 0.2

par(mfrow=c(2,3))
sapply(c(0.05,0.1,0.25), function(cont){
    sapply(c(0.5,0.25), function(split){
        sapply(errorpar,function(error){
            
            numpar <- seq(1,5)
            plot(numpar,c(error,cont,split - 0.16,split,admix),col="red",cex=2,pch=4,xaxt="n",ylim=c(0,1),ylab="Parameter value",xlab="",main="",cex.lab=1.5,cex.axis=1.5)
            axis(1,at=numpar,label=c("error","contamination","drift C","drift A", "admixture rate"),cex.axis=1.2)
            legend("topright",c("True","Estimated"),col=c("red","black"),pch=c(4,1),cex=2)
            bootsnum <- 4
            if(admix == 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff_3/estimate_dadi_t10_s",split,"_a",admix,"_n10000_m30_c",cont,"_e",error,".txt",sep=""))
                print(dim(estimatetab))
            }
            if(admix > 0){
                estimatetab <- read.table(paste("~/TwoPopCont/simul_cont_diff_3/estimate_dadi_t10_s",split,"_a",admix,"_b0.08_n10000_m30_c",cont,"_e",error,".txt",sep=""))
            }
            allest <-estimatetab[1,][numpar]
            allboots <- estimatetab[seq(2,(bootsnum+1)),numpar]
            lowboots <- apply(allboots,2, function(x) { quantile(x,0.025) } )
            highboots <- apply(allboots,2, function(x) { quantile(x,0.975) } )
            add.conf.int(numpar,allest,lowboots,highboots, 0.1, col="black")
            print(lowboots)
            print(highboots)
            print(numpar)
            points(numpar,allest,pch=1,cex=2)
        })
    })
})
