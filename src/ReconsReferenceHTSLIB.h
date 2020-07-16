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
#include "htslib/kstring.h"
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




/* string reconstructRef(const BamAlignment  * al); */
void reconstructRefWithPosHTS(const bam1_t   * b,pair< kstring_t *, vector<int> > &);

/* int numberOfDeletions(const BamAlignment  * al); */

#endif
