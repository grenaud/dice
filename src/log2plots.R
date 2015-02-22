#!/usr/bin/env Rscript

#log2plots.R [nucular log] [output suffix] <burnins>

args=(commandArgs(TRUE))

burnins<-10000;

if(length(args)==3){
  burnins<-args[3];
}

data <- read.table(args[1]);


pdf(paste(args[2],".it.pdf",sep=""));
plot(data[,1],data[,2],xlab="MCMC iteration",ylab="posterior prob.",main="Posterior probability wrt iteration",col="darkred");
dev.off();

data<-data[-c(1:burnins),];

pdf(paste(args[2],".e.pdf",sep=""));
plot(density(data[,3]),xlab="Error rate",ylab="Density",main="Post. prob. for error",col="darkblue");
dev.off();


pdf(paste(args[2],".r.pdf",sep=""));
plot(density(data[,4]),xlab="r",ylab="Density",main="Post. prob. for cont. rate (r)",col="darkgreen");
dev.off();

pdf(paste(args[2],".tauC.pdf",sep=""));
plot(density(data[,5]),xlab=paste(expression("tau"),"_C",sep=""),ylab="Density",main="Post. prob. for contaminant lineage drift",col="darkorange");
dev.off();


pdf(paste(args[2],".tauA.pdf",sep=""));
plot(density(data[,6]),xlab=paste(expression("tau"),"_A",sep=""),ylab="Density",main="Post. prob. for ancient lineage drift",col="purple");
dev.off();


    

#     density(dataa$V1),xlab="xlab",ylab="ylab",main="main",col="red");
#legend("topright", title="sample", c("data"), lty=c(1), lwd=c(2.5),col=c("red")) 
