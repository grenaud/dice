#!/usr/bin/env Rscript
library(coda);
#log2plots.R [dice logs] 

args=(commandArgs(TRUE))

burnins<-10000;

maxFound<- -1.0*Inf;
maxFoundF<-"none";

for(i in 1:length(args)){
  write(paste("file:",args[i],sep=" "),file="/dev/stdout");
  data <- read.table(args[i],header=T);
  m<-max(data$lpos);
  if(m>maxFound){
    maxFound<-m;
    maxFoundF<-args[i];    
  }

  #print(args[i]);
  #print(maxFound);
  #print(maxFoundF);
  data<-data[-c(1:burnins),];

  e<-HPDinterval(as.mcmc(data$error), prob=0.95);
  c<-HPDinterval(as.mcmc(data$ContRate), prob=0.95);
  tC<-HPDinterval(as.mcmc(data$tau_C), prob=0.95);
  tA<-HPDinterval(as.mcmc(data$tau_A), prob=0.95);

write(paste("error:",e[1],e[2],sep=" "),file="/dev/stdout");
write(paste("contamination:",c[1],c[2],sep=" "),file="/dev/stdout");
write(paste("tauC:",tC[1],tC[2],sep=" "),file="/dev/stdout");
write(paste("tauA:",tA[1],tA[2],sep=" "),file="/dev/stdout");

}

write(paste("max found with ",maxFoundF,sep=""),file="/dev/stdout");
data <- read.table(maxFoundF,header=T);

data<-data[-c(1:burnins),];

e<-HPDinterval(as.mcmc(data$error), prob=0.95);
c<-HPDinterval(as.mcmc(data$ContRate), prob=0.95);
tC<-HPDinterval(as.mcmc(data$tau_C), prob=0.95);
tA<-HPDinterval(as.mcmc(data$tau_A), prob=0.95);

write(paste("error:",e[1],e[2],sep=" "),file="/dev/stdout");
write(paste("contamination:",c[1],c[2],sep=" "),file="/dev/stdout");
write(paste("tauC:",tC[1],tC[2],sep=" "),file="/dev/stdout");
write(paste("tauA:",tA[1],tA[2],sep=" "),file="/dev/stdout");

                                        #
#
#
#pdf(paste(args[2],".it.pdf",sep=""));
#plot(data[,1],data[,2],xlab="MCMC iteration",ylab="posterior prob.",main="Posterior probability wrt iteration",col="darkred");
#dev.off();
#
#data<-data[-c(1:burnins),];
#
#pdf(paste(args[2],".e.pdf",sep=""));
#plot(density(data[,3]),xlab="Error rate",ylab="Density",main="Post. prob. for error",col="darkblue");
#dev.off();
#
#
#pdf(paste(args[2],".r.pdf",sep=""));
#plot(density(data[,4]),xlab="r",ylab="Density",main="Post. prob. for cont. rate (r)",col="darkgreen");
#dev.off();
#
#pdf(paste(args[2],".tauC.pdf",sep=""));
#plot(density(data[,5]),xlab=paste(expression("tau"),"_C",sep=""),ylab="Density",main="Post. prob. for contaminant lineage drift",col="darkorange");
#dev.off();
#
#
#pdf(paste(args[2],".tauA.pdf",sep=""));
#plot(density(data[,6]),xlab=paste(expression("tau"),"_A",sep=""),ylab="Density",main="Post. prob. for ancient lineage drift",col="purple");
#dev.off();
#
#
#    
#
##     density(dataa$V1),xlab="xlab",ylab="ylab",main="main",col="red");
##legend("topright", title="sample", c("data"), lty=c(1), lwd=c(2.5),col=c("red")) 
#
