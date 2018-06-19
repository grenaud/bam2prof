/*
 * ReconsReferenceHTSLIB
 * Date: Oct-03-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef ReconsReferenceHTSLIB_h
#define ReconsReferenceHTSLIB_h

#include <vector>
#include <iostream>
#include <string>

#include <stdlib.h> //for exit()
extern "C" {
    //#include "tabix.h"
    //#include "bam.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "bam.h"

#include "samtools.h"
#include "sam_opts.h"

}
#define bam_get_n_cigar_op(b) ((b)->core.n_cigar)
	
/* #include "api/BamMultiReader.h" */
/* #include "api/BamReader.h" */
/* #include "api/BamWriter.h" */
/* #include "api/BamAux.h" */

using namespace std;

typedef struct{
    char bp;
    int offset;
} mdField;

static int asciiOffsetZero=48;
static char DUMMYCHAR='#';




//! To convert the MD first into a vector of mdField structs
/*!
 *
 * This function converts the MD field into a vector of mdField structs
 * we skip deletions in the read (ins in reference)
 *
  \return vector<mdField> a vector of mdField structures
*/

inline vector<mdField> mdString2Vector(const string & mdFieldToParse){
    vector<mdField> toReturn;
    int i=0;
    // int addToOffset=0;
    mdField toadd;
    

    toadd.offset=0;
    toadd.bp='N';

    while(int(mdFieldToParse.length()) != i){
	if(isdigit(mdFieldToParse[i])){
	    toadd.offset=toadd.offset*10+(int(mdFieldToParse[i])-asciiOffsetZero);
	}else{
	    //deletions in read (insertion in reference)
	    if(mdFieldToParse[i] == '^'){
		if(toadd.offset != 0){
		    toadd.bp=DUMMYCHAR;
		    toReturn.push_back(toadd);
		    toadd.offset=0;
		    toadd.bp='N';
		}

		i++;
		mdField toadd2;
		toadd2.offset=0;
		toadd2.bp='^';
		while(isalpha(mdFieldToParse[i])){
		    i++;
		    toadd2.offset++;
		}
		toReturn.push_back(toadd2);
		i--;
	    }else{
		toadd.bp=mdFieldToParse[i];
		toReturn.push_back(toadd);

		toadd.offset=0;
		toadd.bp='N';
	    }

	}
	i++;
    }
    return toReturn;
}

/* string reconstructRef(const BamAlignment  * al); */
pair< string, vector<int> > reconstructRefWithPosHTS(const bam1_t   * b);

/* int numberOfDeletions(const BamAlignment  * al); */

#endif
