/*
 * MistarParser
 * Date: Jan-25-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "MistarParser.h"



MistarParser::MistarParser(string file,string indexForFile,string chrName,int start,int end){
    //for rand() ins SingleAllele
    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );


    rt =new ReadTabix (file,indexForFile,chrName,start,end);
    string headertemp=rt->getHeader();

    istringstream is(headertemp);
  
    //reading the header
    numberPopulations=0;
    populationNames=new vector<string>();

    parseHeader(is);

    if(defline.empty()){
	cerr << "Error: MistarParser cannot get definition line"  <<endl;
	exit(1);
    }

    numberOfTimesHasDataWasCalled=-1;
    tabixMode  = true;
    textMode   = false;
    stringMode = false;
}


MistarParser::MistarParser(const vector<string> * dataToRead_,const vector<string> & populationNames_){
    //populationNames = populationNames_;
    populationNames = new vector<string>( populationNames_);
    numberPopulations = populationNames->size();
    dataToRead      = dataToRead_;
    header="";
    headerNoDefline="";

    numbernew=0;
    numberdel=0;
    dataToReadInd=0;
    numberOfTimesHasDataWasCalled=-1;

    stringMode = true;
    tabixMode  = false;
    textMode   = false;
}

MistarParser::MistarParser(string filename){
    header="";
    headerNoDefline="";

    defline="";


    numbernew=0;
    numberdel=0;

    myFilezipped=new igzstream();
    myFilezipped->open( filename.c_str() );

    //for rand() ins SingleAllele
    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );


    if ( ! myFilezipped->good()) {
	cerr << "Error: MistarParser failed to open file=" << filename <<endl;
	exit(1);;
    }

    //reading the header
    numberPopulations=0;
    populationNames=new vector<string>();

    parseHeader(*myFilezipped);

    if(defline.empty()){
	cerr << "Error: MistarParser cannot get definition line" << filename <<endl;
	exit(1);
    }

    numberOfTimesHasDataWasCalled=-1;
    
    tabixMode  = false;
    textMode   = true;
    stringMode = false;
}

MistarParser::~MistarParser(){
    //cerr<<"destructor MistarParser"<<endl;
    if(numberOfTimesHasDataWasCalled == 0){//last called was getData
	//delete(allRecToReturn->vectorAlleles);
	delete(allRecToReturn);
	numberdel++;
    }

    if(tabixMode){
	delete rt; //calling the destructor
    }


    delete(populationNames);
    if(textMode){
	delete(myFilezipped);
    }
    // else
    // 	delete(myFile);
}

void MistarParser::parseHeader(istream & in){
    bool firstLine=true;
    string line;

    while(getline ( in,line)){
	//cout<<"line "<<line<<endl;
	if(line[0] == '#'){
	    // cout<<line;
	    if(firstLine){
		if(line != "#MISTAR"){
		    cerr << "Error: MistarParser first line must be #MISTAR found: " << line <<endl;
		    exit(1);	    
		}		
		firstLine=false;
		continue;
	    }

	    
	    if(strBeginsWith(line, "#chr")){
		defline=line;
		vector<string> fields=allTokens(line,'\t');
		if(fields[0] != "#chr")   { cerr<<"Field #1 of header must be #chr ";    exit(1); }
		if(fields[1] != "coord")  { cerr<<"Field #2 of header must be coord ";   exit(1); }
		if(fields[2] != "REF,ALT"){ cerr<<"Field #3 of header must be REF,ALT "; exit(1); }
		if(fields[3] != "root")   { cerr<<"Field #4 of header must be root ";    exit(1); }
		if(fields[4] != "anc")    { cerr<<"Field #5 of header must be anc ";     exit(1); }

		for(unsigned int i=3;i<fields.size();i++){
		    populationNames->push_back(fields[i]);
		    numberPopulations++;
		}
		header+=line+"\n";

		break;
	    }else{
		header+=line+"\n";
		headerNoDefline+=line+"\n";
	    }
	    
	}else{
	    cerr << "Error: MistarParser cannot get header"  <<endl;
	    exit(1);
	}
    }
}


string MistarParser::getHeader(string prefix){
    vector<string> fields=allTokens(header,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}



string MistarParser::getHeaderNoDefline(string prefix){
    vector<string> fields=allTokens(headerNoDefline,'\n');
    vector<string> toreturn;
    for(unsigned int i=0;i<fields.size();i++){
	if(!fields[i].empty())
	    toreturn.push_back(prefix+fields[i]);
    }

    return vectorToString(toreturn,"\n");
}




string MistarParser::getDefline(){
    return defline;
}

void MistarParser::repositionIterator(string chrName,int start,int end){
    
    if(!tabixMode){
	cerr<<"The subroutine repositionIterator can only be called on objects constructed using tabix " <<endl;
	exit(1);	
    }
    numberOfTimesHasDataWasCalled=-1;

    rt->repositionIterator(chrName,start,end);
}

bool MistarParser::hasData(){

    if(numberOfTimesHasDataWasCalled!=-1){
	//	cout<<"delete"<<endl;
	//cerr<<"del "<<allRecToReturn<<endl;
	//delete(allRecToReturn->vectorAlleles);
	delete(allRecToReturn);
	numberdel++;
    }else{
	numberOfTimesHasDataWasCalled=0;
    }
    
    numberOfTimesHasDataWasCalled++;
    //    string line;
    //if(getline ( *myFilezipped,line)){
    if(getNextLine()){
	numbernew++;
	allRecToReturn                = new AlleleRecords();
	//	cerr<<"new "<<allRecToReturn<<endl;
	//allRecToReturn->vectorAlleles = new vector<SingleAllele>();
	// cout<<"currentline "<<currentline<<endl;
	
	vector<string> fields=allTokens(currentline,'\t');

	if(fields.size() != (numberPopulations+3)){
	    cerr << "Error: MistarParser the following line should have "<<(numberPopulations+3)<<" fields " << currentline <<endl;
	    exit(1);	   
	}
	if(fields[2].length() != 3){
	    cerr << "Error: MistarParser the following line " << currentline <<" does not have 2 comma separated alleles"<<endl;
	    exit(1);	   
	}

	allRecToReturn->chr        =                           fields[0];
	allRecToReturn->coordinate = destringify<unsigned int>(fields[1]);
	allRecToReturn->ref        =                           fields[2][0];
	allRecToReturn->alt        =                           fields[2][2];
	if(allRecToReturn->ref == allRecToReturn->alt){
	    cerr << "Error: MistarParser the following line " << currentline <<" the reference is equal to the alt allele, exiting"<<endl;
	    exit(1);	   
	}

	allRecToReturn->vectorAlleles = new vector<SingleAllele>();
	for(unsigned int i=3;i<fields.size();i++){
	    unsigned int indexComma=0;
	    unsigned int indexColon=0;
	    for(unsigned int k=0;k<fields[i].size();k++){
		if(fields[i][k]==',')
		    indexComma=k;
		if(fields[i][k]==':')
		    indexColon=k;
	    }

	    if(indexComma == 0 || indexColon == 0 ){
		cerr << "Error: MistarParser problem with the following line " << currentline <<" cannot get allele count"<<endl;
		exit(1);	   
	    }
	    
	    SingleAllele sa (destringify<int>( fields[i].substr(0,indexComma)),
			     destringify<int>( fields[i].substr(indexComma+1,indexColon)),
			     destringify<bool>(fields[i].substr(indexColon+1))   );

	    allRecToReturn->vectorAlleles->push_back(sa);
 	}

	if( allRecToReturn->vectorAlleles->size() != numberPopulations){
	    cerr << "Error: MistarParser problem with the following line " << currentline <<" number of allele count read is not "<<numberPopulations<<endl;
		exit(1);	   	    
	}

	return true;

    }else{//if has no data

	if(textMode){
	    myFilezipped->close();
	}

	// else
	//     myFile->close();
	return false;
    }
    return false;
}


bool MistarParser::getNextLine(){
    if(tabixMode){
	return rt->readLine(currentline);
    }

    if(textMode){
	return getline ( *myFilezipped,currentline);
    }

    if(stringMode){
	//cout<<dataToReadInd<<"\t"<<dataToRead->size()<<endl;
	if(dataToReadInd<dataToRead->size()){
	    currentline = dataToRead->at(dataToReadInd++);
	    return true;
	}else{
	    return false;
	}
	
    }
    
    cerr<<"Invalid state in MistarParser::getNextLine()"<<endl;
    exit(1);
    return false;
}



AlleleRecords * MistarParser::getData(){
    
    if(numberOfTimesHasDataWasCalled != 1){
	cerr<<"The subroutine hasData must have been called once prior to calling getData it was called:  "<<numberOfTimesHasDataWasCalled<<" times " <<endl;
	exit(1);
    }
    numberOfTimesHasDataWasCalled=0;


    return allRecToReturn;
}





const vector<string> * MistarParser::getPopulationsNames() const{
    return populationNames;
}



