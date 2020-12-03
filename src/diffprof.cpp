/*
 * diffprof
 * Date: Jul-16-2019 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "libgab.h"

using namespace std;

int main (int argc, char *argv[]) {
    bool absVal=false;
    
    string usage=string(""+string(argv[0])+" <options>  [prof file 1]  [prof file 2]  "+
                        "\nThis program reads 2 .prof files produced by bam2prof\n"+
                        "and outputs the difference on STDOUT\n"+
			"\n\n\tOther options:\n"+
                        "\t\t"+"-abs\t\t\tUse absolute value (Default: "+stringify( absVal )+")\n"+
                        "\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
        cerr << "Usage "<<usage<<endl;
        return 1;       
    }

    
    for(int i=1;i<(argc-2);i++){

        if(string(argv[i]) == "-abs"  ){
            absVal=true;
            continue;
        }


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;       
    }

    
    string line1;
    string line2;
    
    igzstream myFile1;
    igzstream myFile2;
    
    string filename1 = string(argv[ argc-2 ]);
    string filename2 = string(argv[ argc-1 ]);

    myFile1.open(filename1.c_str(), ios::in);
    myFile2.open(filename2.c_str(), ios::in);
    
    if (!myFile1.good()){
	cerr << "Unable to open file "<<filename1<<endl;
	return 1;
    }
    if (!myFile2.good()){
	cerr << "Unable to open file "<<filename2<<endl;
	return 1;
    }
    
    //headers
    getline (myFile1,line1);
    getline (myFile2,line2);
    
    while ( getline (myFile1,line1)){
	if( getline (myFile2,line2) ){
	    // cout<<line1<<"#"<<line2<<endl;
	    vector<string> fields1 = allTokens(line1, '\t');
	    vector<string> fields2 = allTokens(line2, '\t');
	    if( fields1.size() != fields2.size() ){
		cerr<<"Line:"<<line1<<" from  "<<filename1<<" and "<<line2<<" from "<<filename2<<" have a different # of fields "<<fields1.size()<<" vs. "<<fields2.size() << endl;
		return 1;
	    }

	    vector<double> toprint;
	    for(unsigned int i=0;i<fields1.size();i++){
		double d1  =  destringify<double>( fields1[i] );
		double d2  =  destringify<double>( fields2[i] );
		double d   = d1-d2;
		if(absVal)
		    d=abs(d);
		
		toprint.push_back( d );
	    }
	    cout<<vectorToString(toprint,"\t")<<endl;
	    
	}else{
	    cerr << "Uneven length between  "<<filename1<<" "<<filename2<<endl;
	    return 1;	    
	}
	
    }
    
    if( getline (myFile1,line1) ){
	cerr << "Uneven length between  "<<filename1<<" "<<filename2<<endl;
	return 1;	    	
    }
    
    myFile1.close();
    myFile2.close();


    return 0;
}

