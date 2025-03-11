// utility.cpp
#include "utilityFunctions.h"

const int offset=0;
int numberOfCycles;
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

//A=0,C=1,G=2,T=3
char refToChar[256] = {
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char toIndex[4][4]={
  {0,1,2,3},
  {4,5,6,7},
  {8,9,10,11},
  {12,13,14,15}
};
//a->t,c->g,g->c,t->a
char com[4] = {3,2,1,0};


double returnRatioFS(int num,int denom,double errorToRemove,bool failsafe){
  if(failsafe){
    if(denom==0){
      return (0.0-errorToRemove);
    }else{
      return ((double(num)/double(denom))-errorToRemove);
    }
  }else{
    return ((double(num)/double(denom)-errorToRemove));
  }
};
			    
//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const bam1_t  * b, char *reconstructedReference, const vector<int> &  reconstructedReferencePos, const int & minQualBase, string & refFromFasta, const bam_hdr_t *h, void *bed,bool mask, bool ispaired, bool isfirstpair, std::vector<std::vector<unsigned int>>& typesOfDimer5p, std::vector<std::vector<unsigned int>>& typesOfDimer3p, std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle){ // ,int firstCycleRead,int increment

    char refeBase;
    char readBase;
    int  qualBase;
    //Checking if the 5' is deaminated
    bool isDeam5pS=false; //C->T 5'
    bool isDeam3pS=false; //C->T 3'
    bool isDeam5pD=false; //C->T 5'
    bool isDeam3pD=false; //G->A 3'

    int i;
    //cerr<<"read  "<<bam_get_qname(b)<<endl;
    if(ispaired){ //since we cannot evaluate the 5' ends or 3' ends
	goto iterateLoop;
    }


    i=0; //5p for forward str, 3p for reverse
    //cerr<<"bed1 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<endl;
    if(bed && bed_overlap(bed, h->target_name[b->core.tid], reconstructedReferencePos[i], reconstructedReferencePos[i] + 1)==int(mask)) 	goto eval3pdeam;
    //cerr<<"bed2 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<endl;
    
	refeBase=toupper(reconstructedReference[i]);
    
    //readBase=toupper(         al.QueryBases[i]);
    readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
    //qualBase=int(             al.Qualities[i])-offset;
    qualBase=int(             bam_get_qual(b)[i])-offset;  

    if(qualBase < minQualBase)
	goto eval3pdeam;

    if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	goto eval3pdeam;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }

    refeBase = refToChar[refeBase];
    readBase = refToChar[readBase];
    if( refeBase!=4  && readBase!=4 ){
	// if(al.IsReverseStrand()){ //need to take the complement
	if( bam_is_reverse(b) ){
	    refeBase=com[refeBase];
	    readBase=com[readBase];
	}
	
	if(refeBase == 1 && readBase == 3 ){ //C->T

	    //if(al.IsReverseStrand()){ //3'
	    if( bam_is_reverse(b) ){
		isDeam3pS=true;		
	    }else{                    //5'
		isDeam5pS=true;
		isDeam5pD=true;
	    }
	}


	if(refeBase == 2 && readBase == 0 ){ //G->A

	    //if(al.IsReverseStrand()){ //3'
	    if( bam_is_reverse(b) ){
		isDeam3pD=true;		
	    }else{                    //5'
	    }
	}
	   
    }


 eval3pdeam:
    //i=int(al.QueryBases.size())-1; //3p for forward str, 5p for reverse
    i=b->core.l_qseq-1;
    //cerr<<"bed3 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<endl;
    if(bed && bed_overlap(bed, h->target_name[b->core.tid], reconstructedReferencePos[i], reconstructedReferencePos[i] + 1) == int(mask)) 	goto iterateLoop;
    //cerr<<"bed4 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<endl;
    
    refeBase=toupper(reconstructedReference[i]);
    // readBase=toupper(         al.QueryBases[i]);
    // qualBase=int(              al.Qualities[i])-offset;
    readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
    qualBase=int(             bam_get_qual(b)[i])-offset;  
    
    if(qualBase < minQualBase)
	goto iterateLoop;
    
    if(refeBase == 'S' || refeBase == 'I'){ //don't care about soft clipped or indels
	goto iterateLoop;
    }
    
    if(refeBase == 'M'){//match
	refeBase =  readBase;
    }
    refeBase = refToChar[refeBase];
    readBase = refToChar[readBase];
    if( refeBase!=4  && readBase!=4 ){
	//if(al.IsReverseStrand()){ //need to take the complement
	if( bam_is_reverse(b) ){
	    refeBase=com[refeBase];
	    readBase=com[readBase];
	}
	
	if(refeBase == 1 &&
	   readBase == 3 ){ //C->T

	    // if(al.IsReverseStrand()){ //5'
	    if( bam_is_reverse(b) ){

		isDeam5pS=true;
		isDeam5pD=true;		
	    }else{                    //3'
		isDeam3pS=true;
	    }
	}

	if(refeBase == 2 &&
	   readBase == 0 ){ //G->A

	    // if(al.IsReverseStrand()){ //5'
	    if( bam_is_reverse(b) ){
    
	    }else{                    //3'
		isDeam3pD=true;
	    }
	}
	   


    }

 iterateLoop:
    
    char refBaseFromFasta      = 'N';
    char refBaseFromFastaPrev  = 'N';
    char refBaseFromFastaNext  = 'N';
    int j=0;
    //for(i=0;i<int(al.QueryBases.size());i++,j++){
    for(i=0;i<int(b->core.l_qseq);i++,j++){
        // cout<<i<<endl;
        //cerr<<"bed5 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<" "<<i<<endl;
        if(bed && bed_overlap(bed, h->target_name[b->core.tid], reconstructedReferencePos[j], reconstructedReferencePos[j] + 1) == int(mask)) 	continue;
        //cerr<<"bed6 "<<bed<<" "<<bam_get_qname(b)<<" "<<h->target_name[b->core.tid]<<" "<<reconstructedReferencePos[i]<<" "<<(reconstructedReferencePos[i] + 1)<<" "<<i<<endl;

        refeBase=toupper(reconstructedReference[j]);

        // readBase=toupper(          al.QueryBases[i]);
        // qualBase=int(              al.Qualities[i])-offset;

        readBase=toupper( alphabetHTSLIB[ bam_seqi(bam_get_seq(b),i) ] ); //b->core.l_qseq[i]);
        qualBase=int(             bam_get_qual(b)[i])-offset;  
    
        //cout<<i<<"\t"<<qualBase<<"\t"<<minQualBase<<endl;
        //cout<<"-"<<i<<"\t"<<qualBase<<"\t"<<minQualBase<<endl;
        //cout<<"i="<<i<<" j="<<j<<" "<< refeBase<<" "<<readBase<<" "<<refFromFasta[j+1]<<endl;
        //cerr<<"readBase "<<readBase<<" refeBase "<<refeBase<<endl;
        if( refeBase == 'S'){ //don't care about soft clipped or indels	    
            j--;
            continue;
        }


        
        if( refeBase == 'I'){ //don't care about soft clipped or indels
        //i--;
        continue;
        }


        if(refeBase == 'D'){//deletion
            //j++;
            i--;
            continue;
        }

        if(qualBase < minQualBase)
            continue;
        
        if(refeBase == 'M'){//match
            refeBase =  readBase;

            if(!refFromFasta.empty()){
            refBaseFromFasta         = refFromFasta[j+1];
            refBaseFromFastaPrev     = refFromFasta[j  ];
            refBaseFromFastaNext     = refFromFasta[j+2];		
            if(refeBase != refBaseFromFasta){
                cerr<<"Discrepency#1 for "<<bam_get_qname(b)<<" where the reference base at position "<<i<<" "<<refeBase<<" "<<refBaseFromFasta<<endl;
                exit(1);
            }

            }
            
        
            
        }else{
            if(!refFromFasta.empty()){
            refBaseFromFasta         = refFromFasta[j+1];
            refBaseFromFastaPrev     = refFromFasta[j  ];
            refBaseFromFastaNext     = refFromFasta[j+2];		
            if(refeBase != refBaseFromFasta){
                cerr<<"Discrepency#2 for "<<bam_get_qname(b)<<" where the reference base at position "<<i<<" "<<refeBase<<" "<<refBaseFromFasta<<endl;
                exit(1);
            }

            }

        }

        // cout<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<endl;
        refeBase = refToChar[refeBase];
        readBase = refToChar[readBase];
        if( refeBase!=4  && readBase!=4 ){
            int dist5p=i;
            //int dist3p=int(al.QueryBases.size())-1-i;
            int dist3p=b->core.l_qseq-1-i;
            
            //  if(al.IsReverseStrand()){ //need to take the complement
            if( bam_is_reverse(b) ){
            refeBase=com[refeBase];
            readBase=com[readBase];
            //dist5p=int(al.QueryBases.size())-1-i;
            dist5p=int(b->core.l_qseq)-1-i;
            dist3p=i;
            }

            if(dist5p > MAXLENGTH ||
            dist3p > MAXLENGTH ){
            cerr<<"Molecule found "<<bam_get_qname(b)<<" with length greater than limit"<<endl;
            exit(1);
            }
            

            //mismatches[cycleToUse]++;
            if( !ispaired ||  isfirstpair){
            //cerr<<"increase 5p"<<endl;
            typesOfDimer5p[dist5p][toIndex[refeBase][readBase]]++;
            }

            if( !ispaired || !isfirstpair){
            //cerr<<"increase 3p"<<endl;
            typesOfDimer3p[dist3p][toIndex[refeBase][readBase]]++;
            }
            
            if(!refFromFasta.empty()){
            if(
                ( (refBaseFromFasta     == 'C' && refBaseFromFastaNext == 'G') && !bam_is_reverse(b) ) //!al.IsReverseStrand() )
                ||
                ( (refBaseFromFastaPrev == 'C' && refBaseFromFasta     == 'G') &&  bam_is_reverse(b) ) //al.IsReverseStrand() )
            ){
                //cout<<"   CPG: "<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<" ref:"<<refeBase<<" read:"<<readBase<<" "<<al.IsReverseStrand()<<" same="<<(refeBase==readBase)<<endl;
                if( !ispaired ||  isfirstpair)
                typesOfDimer5p_cpg[dist5p][toIndex[refeBase][readBase]]++;
                if( !ispaired || !isfirstpair)
                typesOfDimer3p_cpg[dist3p][toIndex[refeBase][readBase]]++;

            }else{
                if( isResolvedDNA(refBaseFromFasta)                               &&
                isResolvedDNA(refBaseFromFastaPrev)                           &&
                isResolvedDNA(refBaseFromFastaNext)                           &&
                !(refBaseFromFasta     == 'C' && refBaseFromFastaNext == 'G') &&
                !(refBaseFromFastaPrev == 'C' && refBaseFromFasta     == 'G')
                ){

                //cout<<"nonCPG: "<<refBaseFromFastaPrev<<" "<<refBaseFromFasta<<" "<<refBaseFromFastaNext<<" ref:"<<refeBase<<" read:"<<readBase<<" "<<al.IsReverseStrand()<<" same="<<(refeBase==readBase)<<endl;
                if( !ispaired ||  isfirstpair)
                    typesOfDimer5p_noncpg[dist5p][toIndex[refeBase][readBase]]++;
                if( !ispaired || !isfirstpair)
                    typesOfDimer3p_noncpg[dist3p][toIndex[refeBase][readBase]]++;
                
                }
            }
            }

            if(isDeam5pS){
                if( !ispaired || !isfirstpair)
                    typesOfDimer3pSingle[dist3p][toIndex[refeBase][readBase]]++;
            }

            if(isDeam3pS){
                if( !ispaired ||  isfirstpair)
                    typesOfDimer5pSingle[dist5p][toIndex[refeBase][readBase]]++;
            }

            if(isDeam5pD){
                if( !ispaired || !isfirstpair)
                    typesOfDimer3pDouble[dist3p][toIndex[refeBase][readBase]]++;
            }

            if(isDeam3pD){
                if( !ispaired ||  isfirstpair)
                    typesOfDimer5pDouble[dist5p][toIndex[refeBase][readBase]]++;
            }

        }
    }
}


double dbl2log(const double d,bool phred){
    double t= -10.0*(log(d)/log(10.0));
    // if(d == 0){
    // 	t = 
    // }
    if(phred)
	return t;
    else 
	return d;
}


void countSubsPerRef(bool genomeFileB, IndexedGenome* genome, const bam1_t  * b, std::pair<kstring_t*, std::vector<int>>& reconstructedReference, const int & minQualBase, string & refFromFasta, string & refFromFasta_, const bam_hdr_t *h, void *bed,bool mask, bool ispaired, bool isfirstpair, std::vector<std::vector<unsigned int>>& typesOfDimer5p, std::vector<std::vector<unsigned int>>& typesOfDimer3p, std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle){

	reconstructRefWithPosHTS(b,reconstructedReference);
	if(genomeFileB){
		//string ch = refData[al.RefID].RefName;
		string ch = h->target_name[b->core.tid];
		

		if(genome->name2index.find(ch) == genome->name2index.end()){
		cerr<<"Cannot find chr "<<ch<<endl;
		}else{
		//cout<<"found"<<endl;
		}
		faidx1_t & findx=genome->name2index[ch];
	
	
		unsigned int lengthToExtract = reconstructedReference.first->l;
		for(unsigned int i=0;i<reconstructedReference.first->l;i++){
		if(reconstructedReference.first->s[i] == 'I')
			lengthToExtract--;		
		}
		//int startPos = al.Position;
		int startPos = b->core.pos;
		if(startPos!=0)
		startPos--;
		else
		return;
	
		refFromFasta_ = genome->returnStringCoord(&findx,startPos,(lengthToExtract+2));
		refFromFasta = "";
		refFromFasta=refFromFasta_[0];
		int j=1;
		for(unsigned int i=0;i<reconstructedReference.first->l;i++){		
			if(reconstructedReference.first->s[i] == 'I'){
				refFromFasta+="I";
			}else{
				refFromFasta+=refFromFasta_[j++];
			}
		}
		refFromFasta+=refFromFasta_[ refFromFasta_.size() -1 ];
    }
	
	increaseCounters(b,reconstructedReference.first->s, reconstructedReference.second, minQualBase, refFromFasta, h, bed, mask, ispaired, isfirstpair, typesOfDimer5p, typesOfDimer3p, typesOfDimer5p_cpg, typesOfDimer3p_cpg, typesOfDimer5p_noncpg, typesOfDimer3p_noncpg, typesOfDimer5pDouble, typesOfDimer3pDouble, typesOfDimer5pSingle,typesOfDimer3pSingle); //start cycle numberOfCycles-1

}

vector<vector<unsigned int>> initializeDimerVectors(int maxLength, int innerSize) {
    vector<vector<unsigned int>> dimerVector;
    for (int i = 0; i < maxLength; ++i) {
        dimerVector.push_back(std::vector<unsigned int>(innerSize, 0));
    }
    return dimerVector;
}



// void generateDamageProfile( const std::string& outDir, const std::string& bamfiletopen, const std::string& refId, int lengthMaxToPrint, bool dpFormat, bool hFormat, bool allStr, bool singAnddoubleStr, bool doubleStr, bool singleStr, bool endo, bool genomeFileB, bool cpg, double errorToRemove, bool failsafe, bool phred, const std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer5p, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer3p, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg) {
//     std::string file5p, file3p;
//     std::string file5pDefault, file3pDefault;  // These should be defined or passed as parameters

//     // Creating the outdir
//     if (outDir != "/dev/stdout") {
//         std::string command = "mkdir -p " + outDir;
//         int result = system(command.c_str());
//         if (result != 0) {
//             std::cerr << "Failed to create output directory: " << outDir << std::endl;
//         }
//     } 

//     // Retrieving the basename of the bam file that is read
//     std::string bamfiletopenBase = bamfiletopen.substr(bamfiletopen.find_last_of("/\\") + 1);
//     std::string::size_type const p(bamfiletopenBase.find_last_of('.'));
//     std::string file_base = bamfiletopenBase.substr(0, p);

//     // Set up file names
//     if (outDir == "/dev/stdout") {
//         file5p = file_base + "_" + refId + "_5p.prof";
//         file3p = file_base + "_" + refId + "_3p.prof";
//     } else {
//         file5p = outDir + "/" + file_base + "_" + refId + "_5p.prof";
//         file3p = outDir + "/" + file_base + "_" + refId + "_3p.prof";
//     }

//     // Open 5' file
//     std::ofstream file5pFP(file5p.c_str());
//     if (!file5pFP.is_open()) {
//         std::cerr << "Unable to write to 5p file " << file5p << std::endl;
//         return;
//     }

//     // Write header for 5' file
//     if (dpFormat) file5pFP << "\t";
//     if (hFormat) file5pFP << "pos\t";
//     file5pFP << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << std::endl;

//     // Choose the appropriate 5' dimer type
//     const std::vector<std::vector<unsigned int>>* typesOfDimer5pToUse;
//     if (endo) {
//         typesOfDimer5pToUse = doubleStr ? &typesOfDimer5pDouble : &typesOfDimer5pSingle;
//     } else {
//         typesOfDimer5pToUse = &typesOfDimer5p;
//     }
//     if (genomeFileB) {
//         typesOfDimer5pToUse = cpg ? &typesOfDimer5p_cpg : &typesOfDimer5p_noncpg;
//     }

//     // Write 5' profile data
//     for (int l = 0; l < lengthMaxToPrint; l++) {
//         if (dpFormat) file5pFP << l << "\t";
//         if (hFormat) file5pFP << printIntAsWhitePaddedString(l, int(log10(lengthMaxToPrint)) + 1) << "\t";

//         for (int n1 = 0; n1 < 4; n1++) {   
//             int totalObs = 0;
//             for (int n2 = 0; n2 < 4; n2++) {   
//                 totalObs += (*typesOfDimer5pToUse)[l][4 * n1 + n2];
//             }

//             for (int n2 = 0; n2 < 4; n2++) {   
//                 if (n1 == n2) continue;
                
//                 double ratio = MAX(0.0, returnRatioFS((*typesOfDimer5pToUse)[l][4 * n1 + n2], totalObs, errorToRemove, failsafe));
                
//                 if (allStr || 
//                     (singAnddoubleStr && ((n1 == 1 && n2 == 3) || (n1 == 2 && n2 == 0))) ||
//                     (doubleStr && n1 == 1 && n2 == 3) ||
//                     (singleStr && n1 == 1 && n2 == 3)) {
                    
//                     if (dpFormat) file5pFP << dbl2log(ratio, phred) << " [0..0]";
//                     else if (hFormat) file5pFP << printDoubleAsWhitePaddedString(ratio, 1, 5);
//                     else file5pFP << dbl2log(ratio, phred);
//                 } else {
//                     if (dpFormat) file5pFP << (phred ? "-Inf" : "0.0") << " [0..0]";
//                     else if (hFormat) file5pFP << printDoubleAsWhitePaddedString(0.0, 1, 5);
//                     else file5pFP << (phred ? "-Inf" : "0.0");
//                 }

//                 if (!(n1 == 3 && n2 == 2)) file5pFP << "\t";
//             }
//         }
//         file5pFP << std::endl;
//     }

//     file5pFP.close();

//     // Open 3' file
//     std::ofstream file3pFP;
//     if (outDir == "/dev/stdout") {
//         file3pFP.open(file3p.c_str(), std::ofstream::out | std::ofstream::app);
//     } else {
//         file3pFP.open(file3p.c_str());
//     }

//     if (!file3pFP.is_open()) {
//         std::cerr << "Unable to write to 3p file " << file3p << std::endl;
//         return;
//     }

//     // Write header for 3' file
//     if (dpFormat) file3pFP << "\t";
//     if (hFormat) file3pFP << "pos\t";
//     file3pFP << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << std::endl;

//     // Choose the appropriate 3' dimer type
//     const std::vector<std::vector<unsigned int>>* typesOfDimer3pToUse;
//     if (endo) {
//         typesOfDimer3pToUse = doubleStr ? &typesOfDimer3pDouble : &typesOfDimer3pSingle;
//     } else {
//         typesOfDimer3pToUse = &typesOfDimer3p;
//     }
//     if (genomeFileB) {
//         typesOfDimer3pToUse = cpg ? &typesOfDimer3p_cpg : &typesOfDimer3p_noncpg;
//     }

//     // Write 3' profile data
//     for (int le = 0; le < lengthMaxToPrint; le++) {
//         int l = le;
//         if (dpFormat) {
//             l = lengthMaxToPrint - 1 - le;
//             file3pFP << (l == 0 ? "" : "-") << l << "\t";
//         }
//         if (hFormat) {
//             file3pFP << printIntAsWhitePaddedString(l, int(log10(lengthMaxToPrint)) + 1) << "\t";
//         }

//         for (int n1 = 0; n1 < 4; n1++) {   
//             int totalObs = 0;
//             for (int n2 = 0; n2 < 4; n2++) {   
//                 totalObs += (*typesOfDimer3pToUse)[l][4 * n1 + n2];
//             }

//             for (int n2 = 0; n2 < 4; n2++) {   
//                 if (n1 == n2) continue;
                
//                 double ratio = MAX(0.0, returnRatioFS((*typesOfDimer3pToUse)[l][4 * n1 + n2], totalObs, errorToRemove, failsafe));
                
//                 if (allStr || 
//                     (singAnddoubleStr && ((n1 == 1 && n2 == 3) || (n1 == 2 && n2 == 0))) ||
//                     (doubleStr && n1 == 2 && n2 == 0) ||
//                     (singleStr && n1 == 1 && n2 == 3)) {
                    
//                     if (dpFormat) file3pFP << dbl2log(ratio, phred) << " [0..0]";
//                     else if (hFormat) file3pFP << printDoubleAsWhitePaddedString(ratio, 1, 5);
//                     else file3pFP << dbl2log(ratio, phred);
//                 } else {
//                     if (dpFormat) file3pFP << (phred ? "-Inf" : "0.0") << " [0..0]";
//                     else if (hFormat) file3pFP << printDoubleAsWhitePaddedString(0.0, 1, 5);
//                     else file3pFP << (phred ? "-Inf" : "0.0");
//                 }

//                 if (!(n1 == 3 && n2 == 2)) file3pFP << "\t";
//             }
//         }
//         file3pFP << std::endl;
//     }

//     file3pFP.close();
// }

void generateDamageProfile( const std::string& outDir, const std::string& bamfiletopen, const std::string& refId, int lengthMaxToPrint, bool dpFormat, bool hFormat, bool allStr, bool singAnddoubleStr, bool doubleStr, bool singleStr, bool endo, bool genomeFileB, bool cpg, double errorToRemove, bool failsafe, bool phred, const std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer5p, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer3p, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, uint64_t mapped) {
    std::string file5p, file3p;
    std::string file5pDefault, file3pDefault;
    std::string mappedStr = std::to_string(mapped);


    // Creating the output directory
    if (outDir != "/dev/stdout") {
        std::string command = "mkdir -p " + outDir;
        int result = system(command.c_str());
        if (result != 0) {
            std::cerr << "Failed to create output directory: " << outDir << std::endl;
        }
    }

    // Retrieving the basename of the bam file
    std::string bamfiletopenBase = bamfiletopen.substr(bamfiletopen.find_last_of("/\\") + 1);
    std::string::size_type const p(bamfiletopenBase.find_last_of('.'));
    std::string file_base = bamfiletopenBase.substr(0, p);

    // Set up file names
    file5p = outDir + "/" + file_base + "_" + refId + "_n" + mappedStr + "_5p.prof";
    file3p = outDir + "/" + file_base + "_" + refId + "_n" + mappedStr + "_3p.prof";

    // Open 5' file
    std::ofstream file5pFP(file5p.c_str());
    if (!file5pFP.is_open()) {
        std::cerr << "Unable to write to 5p file " << file5p << std::endl;
        return;
    }

    // Write header for 5' file
    if (dpFormat) file5pFP << "\t";
    if (hFormat) file5pFP << "pos\t";
    file5pFP << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << std::endl;

    // Choose the appropriate 5' dimer type
    const std::vector<std::vector<unsigned int>>* typesOfDimer5pToUse = &typesOfDimer5p;
    if (endo) {
        typesOfDimer5pToUse = doubleStr ? &typesOfDimer5pDouble : &typesOfDimer5pSingle;
    }
    if (genomeFileB) {
        typesOfDimer5pToUse = cpg ? &typesOfDimer5p_cpg : &typesOfDimer5p_noncpg;
    }

    // Write 5' profile data
    for (int l = 0; l < lengthMaxToPrint; l++) {
        if (dpFormat) file5pFP << l << "\t";
        if (hFormat) file5pFP << printIntAsWhitePaddedString(l, int(log10(lengthMaxToPrint)) + 1) << "\t";

        for (int n1 = 0; n1 < 4; n1++) {   
            int totalObs = 0;
            for (int n2 = 0; n2 < 4; n2++) {   
                totalObs += (*typesOfDimer5pToUse)[l][4 * n1 + n2];
            }

            for (int n2 = 0; n2 < 4; n2++) {   
                if (n1 == n2) continue;

                double ratio = MAX(0.0, returnRatioFS((*typesOfDimer5pToUse)[l][4 * n1 + n2], totalObs, errorToRemove, failsafe));

                // Abort if any NaN is found
                if (std::isnan(ratio)) {
                    file5pFP.close();
                    std::remove(file5p.c_str());
                    std::remove(file3p.c_str());  // Ensure both files are deleted
                    return;  // Stop execution immediately
                }

                if (allStr || (singAnddoubleStr && ((n1 == 1 && n2 == 3) || (n1 == 2 && n2 == 0))) || 
                    (doubleStr && n1 == 1 && n2 == 3) || (singleStr && n1 == 1 && n2 == 3)) {
                    
                    file5pFP << dbl2log(ratio, phred);
                } else {
                    file5pFP << (phred ? "-Inf" : "0.0");
                }

                if (!(n1 == 3 && n2 == 2)) file5pFP << "\t";
            }
        }
        file5pFP << std::endl;
    }

    file5pFP.close();

    // Open 3' file
    std::ofstream file3pFP;
    if (outDir == "/dev/stdout") {
        file3pFP.open(file3p.c_str(), std::ofstream::out | std::ofstream::app);
    } else {
        file3pFP.open(file3p.c_str());
    }

    if (!file3pFP.is_open()) {
        std::cerr << "Unable to write to 3p file " << file3p << std::endl;
        return;
    }

    // Write header for 3' file
    if (dpFormat) file3pFP << "\t";
    if (hFormat) file3pFP << "pos\t";
    file3pFP << "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G" << std::endl;

    // Choose the appropriate 3' dimer type
    const std::vector<std::vector<unsigned int>>* typesOfDimer3pToUse;
    if (endo) {
        typesOfDimer3pToUse = doubleStr ? &typesOfDimer3pDouble : &typesOfDimer3pSingle;
    } else {
        typesOfDimer3pToUse = &typesOfDimer3p;
    }
    if (genomeFileB) {
        typesOfDimer3pToUse = cpg ? &typesOfDimer3p_cpg : &typesOfDimer3p_noncpg;
    }

    // Write 3' profile data
    for (int le = 0; le < lengthMaxToPrint; le++) {
        int l = le;
        if (dpFormat) {
            l = lengthMaxToPrint - 1 - le;
            file3pFP << (l == 0 ? "" : "-") << l << "\t";
        }
        if (hFormat) {
            file3pFP << printIntAsWhitePaddedString(l, int(log10(lengthMaxToPrint)) + 1) << "\t";
        }

        for (int n1 = 0; n1 < 4; n1++) {   
            int totalObs = 0;
            for (int n2 = 0; n2 < 4; n2++) {   
                totalObs += (*typesOfDimer3pToUse)[l][4 * n1 + n2];
            }

            for (int n2 = 0; n2 < 4; n2++) {   
                if (n1 == n2) continue;
                
                double ratio = MAX(0.0, returnRatioFS((*typesOfDimer3pToUse)[l][4 * n1 + n2], totalObs, errorToRemove, failsafe));

                // Abort if any NaN is found
                if (std::isnan(ratio)) {
                    file3pFP.close();
                    std::remove(file5p.c_str());
                    std::remove(file3p.c_str());  // Ensure both files are deleted
                    return;  // Stop execution immediately
                }
                
                if (allStr || 
                    (singAnddoubleStr && ((n1 == 1 && n2 == 3) || (n1 == 2 && n2 == 0))) ||
                    (doubleStr && n1 == 2 && n2 == 0) ||
                    (singleStr && n1 == 1 && n2 == 3)) {
                    
                    if (dpFormat) file3pFP << dbl2log(ratio, phred) << " [0..0]";
                    else if (hFormat) file3pFP << printDoubleAsWhitePaddedString(ratio, 1, 5);
                    else file3pFP << dbl2log(ratio, phred);
                } else {
                    if (dpFormat) file3pFP << (phred ? "-Inf" : "0.0") << " [0..0]";
                    else if (hFormat) file3pFP << printDoubleAsWhitePaddedString(0.0, 1, 5);
                    else file3pFP << (phred ? "-Inf" : "0.0");
                }

                if (!(n1 == 3 && n2 == 2)) file3pFP << "\t";
            }
        }
        file3pFP << std::endl;
    }

    file3pFP.close();
}
