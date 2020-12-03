/*
 * mistarOperations
 * Date: Apr-01-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "mistarOperations.h"




//reading fasta index
// void readFastaIndex(const string fastaIndex,
// 		    vector<chrinfo> & chrFound,
// 		    uint64_t & genomeLength){
//     ifstream myFile;
//     string line;

//     myFile.open(fastaIndex.c_str(), ios::in);
//     if (myFile.is_open()){
// 	while ( getline (myFile,line)){
// 	    vector<string> fields = allTokens(line,'\t');
// 	    chrinfo toadd;

// 	    toadd.name         =fields[0];
// 	    toadd.startIndexChr=genomeLength+1;
// 	    toadd.length       =string2uint(fields[1]);
// 	    toadd.endIndexChr  =genomeLength+toadd.length;
// 	    chrFound.push_back(toadd);
// 	    genomeLength+=toadd.length;
// 	    // cout<<"gl "<<genomeLength<<endl;
// 	}
//     }else{
// 	cerr<<"Cannot open fasta index  "<<fastaIndex<<endl;
// 	exit(1);
//     }
//     myFile.close();

// }

bool initFiles(vector<MistarParser * > & vectorOfMP,
	       //bool & atLeastOneHasData,
	       vector<bool> & hasData,
	       vector<int> & popSizePerFile,
	       vector<AlleleRecords *> & vecAlleleRecords,
	       string & chr1,
	       unsigned int & coordCurrent,
	       bool printOnlyFirstPop){
    // cout<<"#chr\tcoord\tREF,ALT\troot\tanc\t";


    //   atLeastOneHasData=false;
   hasData = vector<bool>(vectorOfMP.size(),false);

   //   cerr<<"initFiles " <<endl;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	hasData[i] = vectorOfMP[i]->hasData()  ;
	//	cerr<<i<<"\t"<<hasData[i]<<endl;
	if(!hasData[i]){
	    cerr<<"File #"<<(i+1)<<" does not have any data, exiting"<<endl;
	    // exit(1);    
	    return false;
	}// else{
	//     //atLeastOneHasData=true;
	// }
    }

    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	unsigned int nonPop = vectorOfMP[i]->getPopulationsNames()->size();
	vector<string> pops;
	popSizePerFile.push_back(nonPop);

	for(unsigned int j=2;j<nonPop;j++){
	    pops.push_back( vectorOfMP[i]->getPopulationsNames()->at(j));
	}
	//cout<<vectorToString(pops,"\t");
	if(printOnlyFirstPop)
	    break;

	if( i!=(vectorOfMP.size() -1)){
	    cout<<"\t";
	}       	
    }

    // cout<<"\n";



    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	vecAlleleRecords.push_back( vectorOfMP[i]->getData() );
	if(i==0){
	    chr1          = vecAlleleRecords[i]->chr;
	    coordCurrent  = vecAlleleRecords[i]->coordinate;
	}else{
	    if(chr1 != vecAlleleRecords[i]->chr ){
		cerr<<"Chromosomes differ between "<<chr1<<" and "<< vecAlleleRecords[i]->chr <<endl;
		exit(1);    	
	    }
	    coordCurrent  = min(coordCurrent,vecAlleleRecords[i]->coordinate);
	}	
    }
    return true;

}

bool sanityCheck(vector<MistarParser * > & vectorOfMP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 string & chr1,
		 unsigned int & coordCurrent,
		 string & chrcheck ,
		 char & refAllele,
		 bool force  ){
  

    //cerr<<"size =  "<<vectorOfMP.size()<<"\tref="<<refAllele<<endl;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	//cerr<<i<<" "<<hasData[i]<<" "<<hasCoordinate[i]<<"\ "<<vecAlleleRecords[i]->ref<<endl;
	if(hasData[i] && hasCoordinate[i]){

	    if(refAllele == '\0'){
		chrcheck   = vecAlleleRecords[i]->chr;
		refAllele  = vecAlleleRecords[i]->ref;
		if( chrcheck   != chr1){
		    cerr<<"Chromosomes differ between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    exit(1);   
		}

	    }else{
		if( chrcheck   != vecAlleleRecords[i]->chr){
		    cerr<<"Chromosomes differ between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    exit(1);   
		}

		if( refAllele  != vecAlleleRecords[i]->ref){
		    cerr<<"The reference allele differs between "<<(* vecAlleleRecords[0])<<" and "<<(*vecAlleleRecords[i])<<endl;
		    if(force)
			return false;
		    else
			exit(1);  
		}
	    }
	}
    }

    //cerr<<"size =  "<<vectorOfMP.size()<<"\tref="<<refAllele<<endl;
    if(refAllele == '\0'){
	cerr<<"The reference allele could not be determined at coordinate "<<coordCurrent<<endl;	
	exit(1);  
    }

    return true;

}














bool hasCpG(vector<MistarParser * > & vectorOfMP,
	    vector<bool> & hasData,
	    vector<bool> & hasCoordinate,
	    vector<int> & popSizePerFile,
	    vector<AlleleRecords *> & vecAlleleRecords){

    bool oneIsCpG=false;

    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){


	    for(int k=0;k<popSizePerFile[i];k++){
		
		oneIsCpG = oneIsCpG || (vecAlleleRecords[i]->vectorAlleles->at(k).getIsCpg() );

	    }
	}
    }

    return oneIsCpG;
}


bool populateFreqVec(vector<MistarParser * > & vectorOfMP,
		     vector<bool> & hasData,
		     vector<bool> & hasCoordinate,
		     vector<int> & popSizePerFile,
		     vector<AlleleRecords *> & vecAlleleRecords,
		     string & chr1,
		     unsigned int & coordCurrent,
		     vector<double> * freqVec,
		     char & derAllele,
		     char & ancAllele,
		     bool force){

    //sanity checks
    string chrcheck = "";  //= vecAlleleRecords[0]->chr;
    char refAllele  = '\0'; // = vecAlleleRecords[0]->ref;
    //cerr<<"populateFreqVec"<<endl;
    bool isSane=sanityCheck(vectorOfMP,hasData,hasCoordinate,vecAlleleRecords,chr1,coordCurrent,chrcheck,refAllele,force);
    if(!isSane)
	return false;
    
    //ancestral info
    SingleAllele chimp;
    SingleAllele anc;
    bool chimpSet=false;
    bool ancSet  =false;
	    
    //determining new alternative allele
    char newAlt = 'N';
    vector<SingleAllele> toPrint;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){

	    if( !isResolvedDNA(newAlt)  && //is 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) ){ //not 'N'
		newAlt = vecAlleleRecords[i]->alt;
	    }
		
	    if( isResolvedDNA(newAlt)                   && //not 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) && //not 'N'
		vecAlleleRecords[i]->alt != newAlt){       //differ
		//tri-allelic
		//goto seekdata;
		return false;
	    }
	}
    }
	    
    bool ancIsDefined  = false;
    bool rootIsDefined = false;
    bool ancIsRef      = false;
    bool rootIsRef     = false;

    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){
	    //chimp
	    if(!vecAlleleRecords[i]->vectorAlleles->at(0).alleleCountNull()){
		rootIsDefined=true;

		if(chimpSet){
		    if(chimp != vecAlleleRecords[i]->vectorAlleles->at(0)){
			cerr<<"Discrepancy in chimp info between "<<(* vecAlleleRecords[i])<<" and "<<(chimp)<<endl;
			if(force)
			    return false;
			else
			    exit(1);  
		    }
		}else{
		    chimpSet=true;
		    chimp=vecAlleleRecords[i]->vectorAlleles->at(0);
		    
		    if( (chimp.getRefCount() !=0 ) &&  (chimp.getAltCount() !=0 ) ){
			cerr<<"More than one root allele for "<<(* vecAlleleRecords[i])<<" and "<<(chimp)<<endl;
			exit(1);
		    }

		    rootIsRef= (chimp.getRefCount() !=0);
		    

		    
		}
	    }

	    //anc
	    if(!vecAlleleRecords[i]->vectorAlleles->at(1).alleleCountNull()){
		ancIsDefined=true;
		if(ancSet){
		    if(anc != vecAlleleRecords[i]->vectorAlleles->at(1)){
			cerr<<"Disprency in ancestral info between "<<(* vecAlleleRecords[i])<<" and "<<(anc)<<endl;
			if(force)
			    return false;
			else
			    exit(1);  
		    }
		}else{
		    ancSet=true;
		    anc=vecAlleleRecords[i]->vectorAlleles->at(1);

		    if( (anc.getRefCount() !=0 ) &&  (anc.getAltCount() !=0 ) ){
			cerr<<"More than one anc allele for "<<(* vecAlleleRecords[i])<<" and "<<(anc)<<endl;
			exit(1);
		    }

		    ancIsRef= (anc.getRefCount() !=0);
		    
		}
	    }

	}
    }
    //cout<<endl;
	     

    // 	printnewrecord:
    //cout<<chr1<<"\t"<<coordCurrent<<"\t"<<refAllele<<","<<newAlt<<"\t"<<chimp<<"\t"<<anc<<"\t";
    
    if(!rootIsDefined  || !ancIsDefined)
	return false;


    if(ancIsRef){
	derAllele = newAlt;
	ancAllele = refAllele;
    }else{ //anc is the alt	
	derAllele = refAllele;
	ancAllele = newAlt;
    }

    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){

	    double freq;
	    for(int k=2;k<popSizePerFile[i];k++){

		if( (vecAlleleRecords[i]->vectorAlleles->at(k).getRefCount() == 0 ) &&
		    (vecAlleleRecords[i]->vectorAlleles->at(k).getAltCount() == 0 ) ){
		    return false;
		}
		
		if(ancIsRef){

		    freq= double(vecAlleleRecords[i]->vectorAlleles->at(k).getAltCount() )
			      /
			double(vecAlleleRecords[i]->vectorAlleles->at(k).getRefCount() + vecAlleleRecords[i]->vectorAlleles->at(k).getAltCount()  );

		}else{ //anc is the alt

		    freq= double(vecAlleleRecords[i]->vectorAlleles->at(k).getRefCount() )
			      /
			double(vecAlleleRecords[i]->vectorAlleles->at(k).getRefCount() + vecAlleleRecords[i]->vectorAlleles->at(k).getAltCount()  );

		}
		freqVec->push_back(freq);
	    }
		    
	}else{//has not data

	    return false;

	}
    }
    //cout<<vectorToString(toPrint,"\t")<<endl;

    return true;
}








bool printAllele(vector<MistarParser * > & vectorOfMP,
		 vector<bool> & hasData,
		 vector<bool> & hasCoordinate,
		 vector<int> & popSizePerFile,
		 vector<AlleleRecords *> & vecAlleleRecords,
		 string & chr1,
		 unsigned int & coordCurrent,
		 bool force){

    //sanity checks
    string chrcheck = "";  //= vecAlleleRecords[0]->chr;
    char refAllele  = '\0'; // = vecAlleleRecords[0]->ref;
    bool isSane=sanityCheck(vectorOfMP,hasData,hasCoordinate,vecAlleleRecords,chr1,coordCurrent,chrcheck,refAllele,force);
    if(!isSane)
	return false;
    
    //ancestral info
    SingleAllele chimp;
    SingleAllele anc;
    bool chimpSet=false;
    bool ancSet  =false;
	    
    //determining new alternative allele
    char newAlt = 'N';
    vector<SingleAllele> toPrint;
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){

	    if( !isResolvedDNA(newAlt)  && //is 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) ){ //not 'N'
		newAlt = vecAlleleRecords[i]->alt;
	    }
		
	    if( isResolvedDNA(newAlt)                   && //not 'N'
		isResolvedDNA(vecAlleleRecords[i]->alt) && //not 'N'
		vecAlleleRecords[i]->alt != newAlt){       //differ
		//tri-allelic
		//goto seekdata;
		return false;
	    }
	}
    }
	    
	  
    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){
	    //chimp
	    if(!vecAlleleRecords[i]->vectorAlleles->at(0).alleleCountNull()){
		if(chimpSet){
		    if(chimp != vecAlleleRecords[i]->vectorAlleles->at(0)){
			cerr<<"Disprency in chimp info between "<<(* vecAlleleRecords[i])<<" and "<<(chimp)<<endl;
			if(force)
			    return false;
			else
			    exit(1);  
		    }
		}else{
		    chimpSet=true;
		    chimp=vecAlleleRecords[i]->vectorAlleles->at(0);
		}
	    }

	    //anc
	    if(!vecAlleleRecords[i]->vectorAlleles->at(1).alleleCountNull()){
		if(ancSet){
		    if(anc != vecAlleleRecords[i]->vectorAlleles->at(1)){
			cerr<<"Disprency in ancestral info between "<<(* vecAlleleRecords[i])<<" and "<<(anc)<<endl;
			if(force)
			    return false;
			else
			    exit(1);  
		    }
		}else{
		    ancSet=true;
		    anc=vecAlleleRecords[i]->vectorAlleles->at(1);
		}
	    }

	}
    }
    //cout<<endl;
	     

    // 	printnewrecord:
    cout<<chr1<<"\t"<<coordCurrent<<"\t"<<refAllele<<","<<newAlt<<"\t"<<chimp<<"\t"<<anc<<"\t";
	     

    for(unsigned int i=0;i<vectorOfMP.size();i++){ 
	if(hasData[i] && hasCoordinate[i]){
	    for(int k=2;k<popSizePerFile[i];k++){
		toPrint.push_back(vecAlleleRecords[i]->vectorAlleles->at(k));
	    }
		    
	}else{

	    for(int k=2;k<popSizePerFile[i];k++){
		SingleAllele t;
		toPrint.push_back(t);
	    }

	}
    }
    cout<<vectorToString(toPrint,"\t")<<endl;

    return true;
}




// //! Method to read a sorted bed file (copy pasted from vcfcompute.cpp (did not integrate well with libgab.h))
// /*!
//  *
//  * This method checks for the records being ordered

//   \param filetoread : String with the full path to the file to read
//   \return           : Return(head) a pointer to a map where the key is the chromosome name and the value a vector of genomic ranges
//   \sa  readBEDSortedfile()
// */

// map< string, vector<GenomicRange> * > * readBEDSortedfile(string filetoread){
//     //vector<GenomicRange> toReturn;
//     map< string, vector<GenomicRange> * > * toReturn= new map< string, vector<GenomicRange> *>();
//     igzstream myFile;
//     myFile.open(filetoread.c_str(), ios::in);
//     string line;
//     unsigned int     lastEndCoord = 0;
//     string           lastChrname  = "###";

//     if (myFile.good()){
// 	while ( getline (myFile,line)){	    
// 	    //cout<<line<<endl;
// 	    if(line.empty()){
// 		continue;
// 	    }
// 	    string       chrName;
// 	    unsigned int startCoord;
// 	    unsigned int endCoord;
// 	    vector<string> temp=allTokens(line,'\t');
// 	    if(temp.size() != 3){
// 		cerr << "Error in mistarfilter in readBEDSortedfile(): following line does not have 3 fields"<< line<<endl;
// 		exit(1);		
// 	    }
// 	    chrName     = destringify<string>(temp[0]);
// 	    startCoord  = destringify<unsigned int>(temp[1])+1; //the left coordinate is zero based
// 	    endCoord    = destringify<unsigned int>(temp[2]);

// 	    if(lastChrname != chrName){//new chr found
// 		//check for previously found		
// 		lastChrname  = chrName;
// 		lastEndCoord = endCoord;
// 		if(toReturn->find(chrName) == toReturn->end() ){//not previously found
// 		    //cout<<"new chr "<<chrName<<endl;
// 		    (*toReturn)[chrName]=new vector<GenomicRange>();
// 		}else{
// 		    cerr << "Cannot have multiple records on the same chromosome at different parts of the file "<< line<<endl;
// 		    exit(1);
// 		}
// 	    }else{//stay on same chr
// 		if(startCoord < lastEndCoord ){
// 		    cerr << "Problem with line =  "<<line<<" the start of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
// 		    exit(1);
// 		}

// 		if(endCoord   < lastEndCoord ){
// 		    cerr << "Problem with line =  "<<line<<" the end of the coordinate lesser than the end of the previous record "<<lastEndCoord<<endl;
// 		    exit(1);		    
// 		}
// 	    }

// 	    GenomicRange toadd (chrName,startCoord,endCoord);	    
// 	    (*toReturn)[chrName]->push_back(toadd);
// 	    //cout<<(*toReturn)[chrName]->size()<<endl;
// 	}
// 	myFile.close();
//     }else{
// 	cerr << "Unable to open file "<<filetoread<<endl;
// 	exit(1);
//     }

//     return toReturn;
// }


