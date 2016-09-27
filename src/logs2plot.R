#!/usr/bin/env Rscript-3.2.0
library(coda);
#require(reshape2)
require(ggplot2)

args=(commandArgs(TRUE))


burnins<-10000;

maxFound<- -1.0*Inf;
maxFoundF<-"none";
maxFoundIND<- -1;

erangemin<-1;
erangemax<-0;

crangemin<-1;
crangemax<-0;

tauCrangemin<-10;
tauCrangemax<-0;

tauArangemin<-10;
tauArangemax<-0;

outputprefix<-"outlog";
anch<-"ANCH";

nameCont<-list();
eMe<-list();
eMin<-list();
eMax<-list();

cMe<-list();
cMin<-list();
cMax<-list();

tauC<-list();
tauCMin<-list();
tauCMax<-list();

tauA<-list();
tauAMin<-list();
tauAMax<-list();

llik<-list();
llikMin<- 0;
llikMax<- -1.0*Inf



outputprefix<-"";
anch<-"";

for(i in 1:length(args)){
  write(paste("file:",args[i],sep=" "),file="/dev/stderr");
  data <- read.table(args[i],header=T);
  fields<-strsplit(args[i],"_")[[1]];
  #write(  str(length(strsplit(args[i],"_")[[1]])) ,file="/dev/stdout");
  
  if(length(fields) == 4){

    outputprefix<-fields[1];

    fields2<-strsplit(fields[4],".",fixed=TRUE)[[1]];    
    anch<-fields2[1];
    write(paste("using output prefix:",outputprefix,sep=" "),file="/dev/stderr");
    write(paste("using anchor:",       anch,sep=" "),        file="/dev/stderr");
    write(paste("using contam:",       anch,sep=" "),        file="/dev/stderr");
    
    nameCont<-append(nameCont,fields2[1]);

  }else{
    if(length(fields) == 5){
      
      if(length(outputprefix)==1){
        outputprefix <- fields[1];
        fields2<-strsplit(fields[5],".",fixed=TRUE)[[1]];    
        anch         <- fields2[1];
        write(paste("using output prefix:",outputprefix,sep=" "),        file="/dev/stderr");
        write(paste("using anchor:",       anch,        sep=" "),        file="/dev/stderr");
        write(paste("using contam:",       fields[3],   sep=" "),        file="/dev/stderr");
    
      }else{
        fields2<-strsplit(fields[5],".",fixed=TRUE)[[1]];    
        if(outputprefix != fields[1]){
          write(paste("The file ",args[i],"does not match the expected prefix=",outputprefix,sep=" "),file="/dev/stderr");
          quit();
        }
        write(paste("using contam:",  fields[3],sep=" "),        file="/dev/stderr");
        
        if( anch != fields2[1] ){
          write(paste("The file ",args[i],"does not match the expected anchor=",anch,sep=" "),file="/dev/stderr");
          quit();
        }
      }
      nameCont<-append(nameCont,fields[3]);
      
    }
    
  }

  m<-max(data$lpos);
  if(m>maxFound){
    maxFound<-m;
    maxFoundF<-args[i];
    maxFoundIND <- i;
  }

  #print(args[i]);
  #print(maxFound);
  #print(maxFoundF);
  data<-data[-c(1:burnins),];

  e<-HPDinterval(as.mcmc(data$error), prob=0.95);

  eMe<-append(eMe,mean(data$error));
  eMin<-append(eMin,e[1]);
  eMax<-append(eMax,e[2]);

  erangemin<-min(  erangemin , e[1]);
  erangemax<-max(  erangemin , e[2]);

  
  c<-HPDinterval(as.mcmc(data$ContRate), prob=0.95);


  cMe<-append(cMe,mean(data$ContRate));
  cMin<-append(cMin,c[1]);
  cMax<-append(cMax,c[2]);

  #write(paste("error:",e[1],e[2],sep=" "),file="/dev/stdout");
  crangemin<-min(  crangemin , c[1]);
  crangemax<-max(  crangemax , c[2]);

  tC<-HPDinterval(as.mcmc(data$tau_C), prob=0.95);

  tauC<-append(tauC,mean(data$tau_C));
  tauCMin<-append(tauCMin,tC[1]);
  tauCMax<-append(tauCMax,tC[2]);


  tauCrangemin<-min(  tauCrangemin , tC[1]);
  tauCrangemax<-max(  tauCrangemax , tC[2]);

  tA<-HPDinterval(as.mcmc(data$tau_A), prob=0.95);

  tauA<-append(tauA,mean(data$tau_A));
  tauAMin<-append(tauAMin,tA[1]);
  tauAMax<-append(tauAMax,tA[2]);

  tauArangemin<-min(  tauArangemin , tA[1]);
  tauArangemax<-max(  tauArangemax , tA[2]);

  llik<-append(llik,m);
  llikMin<-min(llikMin,m);
  llikMax<-max(llikMax,m);


}




#l<-list( 'nameCont'=nameCont,'e'=eMe,'eMin'=eMin,'eMax'=eMax, 'c'=cMe,'cMin'=cMin,'cMax'=cMax, 'tauC'=tauC,'tauCMin'=tauCMin,'tauCMax'=tauCMax, 'tauA' = tauA,'tauAMin'=tauAMin,'tauAMax'=tauAMax);

nameCont <- unlist( nameCont );
eMe <- unlist( eMe );
eMin <- unlist( eMin );
eMax <- unlist( eMax );
cMe <- unlist( cMe );
cMin <- unlist( cMin );
cMax <- unlist( cMax );
tauC <- unlist( tauC );
tauCMin <- unlist( tauCMin );
tauCMax <- unlist( tauCMax );
tauA <- unlist( tauA );
tauAMin <- unlist( tauAMin );
tauAMax <- unlist( tauAMax );
llik <- unlist( llik );
#llikMin <- unlist( llikMin );
#llikMax <- unlist( llikMax );

d<-data.frame(nameCont,eMe,eMin,eMax,cMe,cMin,cMax,tauC,tauCMin,tauCMax,tauA,tauAMin,tauAMax , llik   );
colnames(d)<-c( 'nameCont','e','eMin','eMax', 'c','cMin','cMax', 'tauC','tauCMin','tauCMax', 'tauA' ,'tauAMin','tauAMax', 'llik' );

ds <- d[order(-d$llik) ,];
ds$nameCont <- factor(ds$nameCont, levels = ds$nameCont);



pdf(paste(outputprefix,"_e.pdf",sep=""));
ewiggle<-erangemax/20
ggplot(ds, aes( x=nameCont,y=e )) + coord_cartesian(ylim = c(max(erangemin-ewiggle,0),erangemax+ewiggle))+geom_bar(position= position_dodge(width = 0.9), stat="identity",fill="cyan") +geom_errorbar( aes(ymax = eMax, ymin= eMin), position= position_dodge(width = 0.9), lwd=0.7,width=0.4) + theme_bw() + xlab("\nPopulation code for the contaminant source") + ylab("estimated error\n") +  ggtitle("Estimated error for each contaminant population\n") + theme(axis.text.x = element_text(angle = 90, hjust = 1));
dev.off();


pdf(paste(outputprefix,"_c.pdf",sep=""));

cwiggle<-crangemax/10

ggplot(ds, aes( x=nameCont,y=c )) + coord_cartesian(ylim = c(max(crangemin-cwiggle,0),crangemax+cwiggle))+geom_bar(position=position_dodge(width = 0.9), stat="identity",fill="brown1") +geom_errorbar( aes(ymax = cMax, ymin= cMin), position= position_dodge(width = 0.9), lwd=0.7,width=0.4) + theme_bw() + xlab("\nPopulation code for the contaminant source") + ylab("Contamination rate\n") +  ggtitle("Estimated present-day contamination\nfor each contaminant population\n") + theme(axis.text.x = element_text(angle = 90, hjust = 1));

dev.off();


pdf(paste(outputprefix,"_tauC.pdf",sep=""));

tauCwiggle<-tauCrangemax/10

ggplot(ds, aes( x=nameCont,y=tauC )) + coord_cartesian(ylim = c(max(tauCrangemin-tauCwiggle,0),tauCrangemax+tauCwiggle))+geom_bar(position=position_dodge(width = 0.9), stat="identity",fill="chartreuse") +geom_errorbar( aes(ymax = tauCMax, ymin= tauCMin), position= position_dodge(width = 0.9), lwd=0.7,width=0.4) + theme_bw() + xlab("\nPopulation code for the contaminant source") + ylab("Drift time for the contaminant\n") +  ggtitle("Estimated drift time for the contamination\npopulation for each contaminant population\n") + theme(axis.text.x = element_text(angle = 90, hjust = 1));
dev.off();

pdf(paste(outputprefix,"_tauA.pdf",sep=""));

tauAwiggle<-tauArangemax/20

ggplot(ds, aes( x=nameCont,y=tauA )) + coord_cartesian(ylim = c(max(tauArangemin-tauAwiggle,0),tauArangemax+tauAwiggle))+geom_bar(position=position_dodge(width = 0.9), stat="identity",fill="yellow2") +geom_errorbar( aes(ymax = tauAMax, ymin= tauAMin), position= position_dodge(width = 0.9), lwd=0.7,width=0.4) + theme_bw() + xlab("\nPopulation code for the contaminant source") + ylab("Drift time for the archaic sample\n") +  ggtitle("Estimated drift time for the archaic\npopulation for each contaminant population\n") + theme(axis.text.x = element_text(angle = 90, hjust = 1));

dev.off();

