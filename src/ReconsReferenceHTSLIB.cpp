/*
 * ReconsReferenceHTSLIB
 * Date: Oct-03-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 * modified by tsk 7july 2020
 */

#include "ReconsReferenceHTSLIB.h"



//! To convert the MD first into a vector of mdField structs
/*!
 *
 * This function converts the MD field into a vector of mdField structs
 * we skip deletions in the read (ins in reference)
 *
  \return vector<mdField> a vector of mdField structures
*/

void  mdString2Vector(const char * mdFieldToParse,vector<mdField> &toReturn){
  toReturn.clear();
    int i=0;
    // int addToOffset=0;
    mdField toadd;
    

    toadd.offset=0;
    toadd.bp='N';

    while(int(strlen(mdFieldToParse)) != i){//tsk compiler warning
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
}


//! This function returns a pair with the string representation of the reference and the vector of the positions on the ref
/*!
 *
 * This function parses the MD field and CIGAR and reconstructs the reference
 * useful for detecting mutations
 *
  \return A pair with the string representation of the reference and the vector of the positions on the ref
*/
static int maxlength=1024;
static char *reconstructedTemp=(char*)calloc(maxlength,1);
void  reconstructRefWithPosHTS(const bam1_t   * b,pair< kstring_t *, vector<int> > &smart){
  smart.first->l = 0; 
  smart.second.clear();
  memset(reconstructedTemp,0,maxlength);
  
  static vector<mdField> parsedMD;
    //initialize
    // int editDist=-1;


    //skip unmapped
    if( ((b)->core.flag&BAM_FUNMAP) != 0 ){
	cerr<<"The function reconstructRefWithPosOnReadHTS()  cannot be called for unmapped reads"<<endl;
	exit(1);
    }
    //cerr<<"reconstructRefWithPosOnReadHTS"<<endl;
    //get relevant data
    // if(!al->GetTag("NM",editDist)){
    // 	cerr<<"Cannot get NM tag from "<<al->Name<<endl;
    // 	exit(1);
    // }
    // if(!al->GetTag("MD",mdFieldString)){
    // 	cerr<<"ReconsReferenceHTSLIB: Cannot get MD tag from "<<al->Name<<endl;
    // 	exit(1);
    // }
    uint8_t *mdptr = bam_aux_get(b, "MD");
    //cerr<<"reconstructRefWithPosOnReadHTS "<<(mdptr+1)<<endl;
    // cerr<<"rg1 "<<rgptr<<endl;
    // cout<<"isize "<<isize<<endl;
            
    if(mdptr==NULL){
	cerr<<"ReconsReferenceHTSLIB: Cannot get MD tag from "<<bam_get_qname(b)<<endl;
	exit(1);
    }

    int32_t   n_cigar_op = bam_get_n_cigar_op(b);
    uint32_t *cigar      = bam_get_cigar(b);
    int at =0;
    for(int32_t i = 0; i < n_cigar_op; i++){
	char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);	
	memset(reconstructedTemp+at,opchr,oplen);
	at += oplen;
    }
    //exit(1);
    /*    
    vector<CigarOp> cigarData=al->CigarData;
    for(unsigned int i=0;i<cigarData.size();i++){
	reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
    }
    */

    //get a vector representation of the MD field	

    mdString2Vector((char *)mdptr+1,parsedMD);

    
    
    //int initialPositionControl=al->Position;
    int initialPositionControl=b->core.pos;

    //combine the CIGAR and MD into one single string
    int mdVectorIndex=0;

    for(unsigned int i=0;i<strlen(reconstructedTemp);i++){
	if(reconstructedTemp[i] == 'M' ){ //only look at matches and indels	    
		
	    if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

		    if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
			i--;
		    }else{
			kputc(parsedMD[mdVectorIndex].bp,smart.first);
			smart.second.push_back(initialPositionControl++);
		    }
		    mdVectorIndex++;
		}else{ //wait until we reach a mismatch
		  kputc(reconstructedTemp[i],smart.first);
		    parsedMD[mdVectorIndex].offset--;
		    smart.second.push_back(initialPositionControl++);
		}

		//skipping all the positions with deletions on the read
		//if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		while( (mdVectorIndex<int(parsedMD.size())) &&
		       (parsedMD[mdVectorIndex].bp == '^' ) ){ 
		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
		    mdVectorIndex++;
		}
		    
	    }else{
	      kputc(reconstructedTemp[i],smart.first);
	      smart.second.push_back(initialPositionControl++);
	    }
	}else{
	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
	      kputc(reconstructedTemp[i],smart.first);
	      smart.second.push_back(initialPositionControl);
	    }
	}
    }

    if(int(smart.first->l) != b->core.l_qseq){
	cerr << "Could not recreate the sequence for read "<<bam_get_qname(b)  << endl;
	exit(1);
    }

    if(smart.second.size() != smart.first->l){
	cerr << "Could not determine the positions for the read "<<bam_get_qname(b) << endl;
	exit(1);
    }
}
