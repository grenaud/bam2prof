/*
 * ReconsReferenceHTSLIB
 * Date: Oct-03-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "ReconsReferenceHTSLIB.h"


// int numberOfDeletions(const BamAlignment  * al){
//     int toReturn=0;
//     vector<CigarOp> cigarData=al->CigarData;
//     for(unsigned int i=0;i<cigarData.size();i++){
// 	if(cigarData[i].Type == 'D')
// 	    toReturn+=cigarData[i].Length;
// 	//reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
//     }
//     return toReturn;
// }

//! This function returns a string representation of the reference
/*!
 *
 * This function parses the MD field and CIGAR and reconstructs the reference
 * useful for detecting mutations
 *
  \return A string representation of the reference
*/
/*
string reconstructRef(const BamAlignment  * al){
    pair< string, vector<int> > toReturn = reconstructRefWithPos(al);
    return toReturn.first;

    // //initialize
    // // int editDist=-1;
    // string mdFieldString="";
    // string reconstructed="";
    // string reconstructedTemp="";

    // //skip unmapped
    // if(!al->IsMapped()){
    // 	cerr<<"The function reconstructRef()  cannot be called for unmapped reads"<<endl;
    // 	exit(1);
    // }

    // //get relevant data
    // // if(!al->GetTag("NM",editDist)){
    // // 	cerr<<"Cannot get NM tag from "<<al->Name<<endl;
    // // 	exit(1);
    // // }
    // if(!al->GetTag("MD",mdFieldString)){
    // 	cerr<<"ReconsReferenceHTSLIB: Cannot get MD tag from "<<al->Name<<endl;
    // 	exit(1);
    // }
	
    // vector<CigarOp> cigarData=al->CigarData;
    // for(unsigned int i=0;i<cigarData.size();i++){
    // 	reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
    // }


    // //get a vector representation of the MD field	

    // vector<mdField> parsedMD=mdString2Vector(mdFieldString);

    // vector<int> positionsOnControl;
    // int initialPositionControl=al->Position;

    // //combine the CIGAR and MD into one single string
    // int mdVectorIndex=0;

    // for(unsigned int i=0;i<reconstructedTemp.size();i++){
    // 	if(reconstructedTemp[i] == 'M' ){ //only look at matches and indels	    
		
    // 	    if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

    // 		if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

    // 		    if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
    // 			i--;
    // 		    }else{
    // 			reconstructed+=parsedMD[mdVectorIndex].bp;
    // 			positionsOnControl.push_back(initialPositionControl++);
    // 		    }
    // 		    mdVectorIndex++;
    // 		}else{ //wait until we reach a mismatch
    // 		    reconstructed+=reconstructedTemp[i];
    // 		    parsedMD[mdVectorIndex].offset--;
    // 		    positionsOnControl.push_back(initialPositionControl++);
    // 		}


    // 		//skipping all the positions with deletions on the read
    // 		//if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches
    // 		while( (mdVectorIndex<int(parsedMD.size())) &&
    // 		       parsedMD[mdVectorIndex].bp == '^'){ 
    // 		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
    // 		    mdVectorIndex++;
    // 		}
    // 		//}
		    
    // 	    }else{
    // 		reconstructed+=reconstructedTemp[i];
    // 		positionsOnControl.push_back(initialPositionControl++);
    // 	    }
    // 	}else{
    // 	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
    // 		reconstructed+=reconstructedTemp[i];
    // 		positionsOnControl.push_back(initialPositionControl);
    // 	    }
    // 	}
    // }

    // if(reconstructed.size() != al->QueryBases.size()){
    // 	cerr << "Could not recreate the sequence for read "<<al->Name << endl;
    // 	exit(1);
    // }

    // if(positionsOnControl.size() != reconstructed.size()){
    // 	cerr << "Could not determine the positions for the read "<<al->Name << endl;
    // 	exit(1);
    // }


    // return reconstructed;

}
*/



//! This function returns a pair with the string representation of the reference and the vector of the positions on the ref
/*!
 *
 * This function parses the MD field and CIGAR and reconstructs the reference
 * useful for detecting mutations
 *
  \return A pair with the string representation of the reference and the vector of the positions on the ref
*/
pair< string, vector<int> >  reconstructRefWithPosHTS(const bam1_t   * b){
    //initialize
    // int editDist=-1;
    string mdFieldString="";
    string reconstructed="";
    string reconstructedTemp="";

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
            
    if(mdptr){
	mdFieldString = string( (const char*)(mdptr+1));
    }else{
	cerr<<"ReconsReferenceHTSLIB: Cannot get MD tag from "<<bam_get_qname(b)<<endl;
	exit(1);
    }

    int32_t   n_cigar_op = bam_get_n_cigar_op(b);
    uint32_t *cigar      = bam_get_cigar(b);

    for(int32_t i = 0; i < n_cigar_op; i++){
	char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);
	reconstructedTemp+=string(oplen,opchr);
    }
    //exit(1);
    /*    
    vector<CigarOp> cigarData=al->CigarData;
    for(unsigned int i=0;i<cigarData.size();i++){
	reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
    }
    */

    //get a vector representation of the MD field	

    vector<mdField> parsedMD=mdString2Vector(mdFieldString);

    vector<int> positionsOnControl;
    
    //int initialPositionControl=al->Position;
    int initialPositionControl=b->core.pos;

    //combine the CIGAR and MD into one single string
    int mdVectorIndex=0;

    for(unsigned int i=0;i<reconstructedTemp.size();i++){
	if(reconstructedTemp[i] == 'M' ){ //only look at matches and indels	    
		
	    if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

		    if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
			i--;
		    }else{
			reconstructed+=parsedMD[mdVectorIndex].bp;
			positionsOnControl.push_back(initialPositionControl++);
		    }
		    mdVectorIndex++;
		}else{ //wait until we reach a mismatch
		    reconstructed+=reconstructedTemp[i];
		    parsedMD[mdVectorIndex].offset--;
		    positionsOnControl.push_back(initialPositionControl++);
		}

		//skipping all the positions with deletions on the read
		//if(mdVectorIndex<int(parsedMD.size())){ //still have mismatches

		while( (mdVectorIndex<int(parsedMD.size())) &&
		       (parsedMD[mdVectorIndex].bp == '^' ) ){ 
		    initialPositionControl+=parsedMD[mdVectorIndex].offset;
		    mdVectorIndex++;
		}
		    
	    }else{
		reconstructed+=reconstructedTemp[i];
		positionsOnControl.push_back(initialPositionControl++);
	    }
	}else{
	    if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
		reconstructed+=reconstructedTemp[i];
		positionsOnControl.push_back(initialPositionControl);
	    }
	}
    }

    if(reconstructed.size() != b->core.l_qseq){
	cerr << "Could not recreate the sequence for read "<<bam_get_qname(b)  << endl;
	exit(1);
    }

    if(positionsOnControl.size() != reconstructed.size()){
	cerr << "Could not determine the positions for the read "<<bam_get_qname(b) << endl;
	exit(1);
    }


    return pair< string, vector<int> >(reconstructed,positionsOnControl);
}
