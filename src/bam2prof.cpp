#include <iostream>
#include <vector>
#include <cstring>
#include <cstring>
#include <set>
#include <ctype.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <sys/types.h>

#include <sys/stat.h>
#include <fcntl.h>
#include "libgab.h"

extern "C" {
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "bam.h"
#include "samtools.h"
#include "sam_opts.h"
#include "bedidx.h"
}
#include "ReconsReferenceHTSLIB.h"

#define bam_is_reverse(b)     (((b)->core.flag&BAM_FREVERSE)    != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)

#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)

#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_sec(b) || bam_is_supp(b) )
#define bam_mqual(b)          ((b)->core.qual)

using namespace std;


const int offset=0;
int numberOfCycles;
string alphabetHTSLIB = "NACNGNNNTNNNNNNN";

#define MAXLENGTH 1000


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



/* it is the structure created by samtools faidx */
typedef struct faidx1_t {
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
}faidx1_t,*FaidxPtr;
/**
 * wrapper for a mmap, a fileid and some faidx indexes
 */
class IndexedGenome
{
private:
    /* used to get the size of the file */
    struct stat buf;
    /* genome fasta file file descriptor */
    int fd;
    /* the mmap (memory mapped) pointer */
    char *mapptr;
    /** reads an fill a string */
    bool readline(gzFile in,string& line)
    {
	if(gzeof(in)) return false;
	line.clear();
	int c=-1;
	while((c=gzgetc(in))!=EOF && c!='\n') line+=(char)c;
	return true;
    }
	    
public:
    /* maps a chromosome to the samtools faidx index */
    map<string,faidx1_t> name2index;

    /** constructor 
     * @param fasta: the path to the genomic fasta file indexed with samtools faidx
     */
    IndexedGenome(const char* fasta):fd(-1),mapptr(NULL)
    {
	string faidx(fasta);
	//cout<<fasta<<endl;
	string line;
	faidx+=".fai";
	/* open *.fai file */
	//cout<<faidx<<endl;
	ifstream in(faidx.c_str(),ios::in);
	if(!in.is_open()){
	    cerr << "cannot open " << faidx << endl;
	    exit(EXIT_FAILURE);
	}
	/* read indexes in fai file */
	while(getline(in,line,'\n'))
	    {
		if(line.empty()) continue;
		const char* p=line.c_str();
		char* tab=(char*)strchr(p,'\t');
		if(tab==NULL) continue;
		string chrom(p,tab-p);
		++tab;
		faidx1_t index;
		if(sscanf(tab,"%ld\t%ld\t%d\t%d",
			  &index.len, &index.offset, &index.line_blen,&index.line_len
		)!=4)
		    {
			cerr << "Cannot read index in "<< line << endl;
			exit(EXIT_FAILURE);
		    }
		/* insert in the map(chrom,faidx) */
		name2index.insert(make_pair(chrom,index));
	    }
	/* close index file */
	in.close();

	/* get the whole size of the fasta file */
	if(stat(fasta, &buf)!=0)
	    {
		perror("Cannot stat");
		exit(EXIT_FAILURE);
	    }
			
	/* open the fasta file */
	fd = open(fasta, O_RDONLY);
	if (fd == -1)
	    {
		perror("Error opening file for reading");
		exit(EXIT_FAILURE);
	    }
	/* open a memory mapped file associated to this fasta file descriptor */
	mapptr = (char*)mmap(0, buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (mapptr == MAP_FAILED)
	    {
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	    }
    }
    /* destructor */
    ~IndexedGenome()
    {
	/* close memory mapped map */
	if(mapptr!=NULL && munmap(mapptr,buf.st_size) == -1)
	    {
		perror("Error un-mmapping the file");
	    }
	/* dispose fasta file descriptor */
	if(fd!=-1) close(fd);
    }
			

    /* return the base at position 'index' for the chromosome indexed by faidx */
    string returnStringCoord(const FaidxPtr faidx,int64_t index, unsigned int length){

	int64_t index2=index;
	// int64_t st=index2;
	// int64_t en=index2+length;
	string strToReturn="";
	
	for(unsigned int j=0;j<length;j++){ //for each char
	    long pos= faidx->offset +
		index2 / faidx->line_blen * faidx->line_len +
		index2 % faidx->line_blen
		;
	    //cout<<char(toupper(mapptr[pos]));
	    strToReturn+=char(toupper(mapptr[pos]));
	    index2++;
	}
	
	return strToReturn;
    }//end returnStringCoord

};
// vector<unsigned int>    matches;
// vector<unsigned int> mismatches;

vector< vector<unsigned int> > typesOfDimer5p; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p; //3' deam rates

vector< vector<unsigned int> > typesOfDimer5p_cpg; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p_cpg; //3' deam rates
vector< vector<unsigned int> > typesOfDimer5p_noncpg; //5' deam rates
vector< vector<unsigned int> > typesOfDimer3p_noncpg; //3' deam rates


vector< vector<unsigned int> > typesOfDimer5pDouble; //5' deam rates when the 3' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer3pDouble; //3' deam rates when the 5' is deaminated according to a double str.
vector< vector<unsigned int> > typesOfDimer5pSingle; //5' deam rates when the 3' is deaminated according to a single str.
vector< vector<unsigned int> > typesOfDimer3pSingle; //3' deam rates when the 5' is deaminated according to a single str.



inline double returnRatioFS(int num,int denom,double errorToRemove,bool failsafe=false){
  if(failsafe){
    if(denom==0){
      return (0.0-errorToRemove);
    }else{
      return ((double(num)/double(denom))-errorToRemove);
    }
  }else{
    return ((double(num)/double(denom)-errorToRemove));
  }
}
			    
//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const   bam1_t  * b,char *reconstructedReference,const vector<int> &  reconstructedReferencePos,const int & minQualBase,const string & refFromFasta, const bam_hdr_t *h,void *bed,bool mask,bool ispaired,bool isfirstpair){ // ,int firstCycleRead,int increment

    char refeBase;
    char readBase;
    int  qualBase;
    //    int cycleToUse=firstCycleRead;
    // cout<<"name "<<al.Name<<endl;
    // cout<<"firstCycleRead "<<firstCycleRead<<endl;
    // cout<<"increment      "<<increment<<endl;
    // cout<<refFromFasta<<endl;
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

int main (int argc, char *argv[]) {

    string file5p="/dev/stdout";
    string file3p="/dev/stdout";
    bool endo=false;

    bool allStr   =true;
    bool singleStr=false;
    bool doubleStr=false;
    bool singAnddoubleStr=false;

    int lengthMaxToPrint = 5;
    int minQualBase      = 0;
    int minLength        = 35;
    bool dpFormat=false;
    bool hFormat=false;
    double errorToRemove=0.0;
    bool phred=false;
    string genomeFile;
    bool   genomeFileB=false;
    IndexedGenome* genome=NULL;
    bool cpg=false;
    string bedfilename;
    void *bed = 0; // BED data structure

    bool bedF=false;
    bool mask=false;
    bool failsafe=false;
    bool paired=false;
    bool quiet=false;

    string usage=string(""+string(argv[0])+" <options>  [in BAM file]"+
			"\nThis program reads a BAM file and produces a deamination profile for the\n"+
			"5' and 3' ends\n"+

			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+

			"\n\n\tOther options:\n"+
			"\t\t"+"-minq\t\t\tRequire the base to have at least this quality to be considered (Default: "+stringify( minQualBase )+")\n"+
			"\t\t"+"-minl\t\t\tRequire the fragment/read to have at least this length to be considered (Default: "+stringify( minLength )+")\n"+
			"\t\t"+"-endo\t\t\tRequire the 5' end to be deaminated to compute the 3' end and vice-versa (Default: "+stringify( endo )+")\n"+
			"\t\t"+"-length\t[length]\tDo not consider bases beyond this length  (Default: "+stringify(lengthMaxToPrint)+" ) \n"+
			"\t\t"+"-err\t[error rate]\tSubstract [error rate] from the rates to account for sequencing errors  (Default: "+stringify(errorToRemove)+" ) \n"+
			"\t\t"+"-log\t\t\tPrint substitutions on a PHRED logarithmic scale  (Default: "+stringify(phred)+" ) \n"+
			"\t\t"+"-bed\t[bed file]\tOnly consider positions in the bed file  (Default: "+booleanAsString( bedF )+" ) \n"+
			"\t\t"+"-mask\t[bed file]\tMask these positions in the bed file    (Default: "+booleanAsString( mask )+" ) \n"+
			"\t\t"+"-paired\t\t\tAllow paired reads    (Default: "+booleanAsString( paired )+" ) \n"+


			"\n\n\tYou can specify either one of the two:\n"+
			"\t\t"+"-single\t\t\tUse the deamination profile of a single strand library  (Default: "+booleanAsString( singleStr )+")\n"+
			"\t\t"+"-double\t\t\tUse the deamination profile of a double strand library  (Default: "+booleanAsString( doubleStr )+")\n"+
			"\n\tor specify this option:\n"+
			"\t\t"+"-both\t\t\tReport both C->T and G->A regardless of stand  (Default: "+booleanAsString( singAnddoubleStr )+")\n"+

			"\n\n\tOutput options:\n"+
			"\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5p)+")\n"+
			"\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3p)+")\n"+
			"\t\t"+"-dp\t\t\tOutput in damage-patterns format (Default: "+booleanAsString(dpFormat)+")\n"+
			"\t\t"+"-h\t\t\tMore human readible output (Default: "+booleanAsString(hFormat)+")\n"+
			"\t\t"+"-q\t\t\tDo not print why reads are skipped (Default: "+booleanAsString(quiet)+")\n"+
		       
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "--help") )
    ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }

    
    for(int i=1;i<(argc-1);i++){ //all but the last 3 args


        if(string(argv[i]) == "-dp"  ){
            dpFormat=true;
            continue;
        }

        if(string(argv[i]) == "-log"  ){
            phred=true;
            continue;
        }


        if(string(argv[i]) == "-paired"  ){
            paired=true;
            continue;
        }


        if(string(argv[i]) == "-bed"  ){
            bedfilename=string(argv[i+1]);
	    bedF=true;
	    i++;
            continue;
        }

        if(string(argv[i]) == "-mask"  ){
            bedfilename=string(argv[i+1]);
	    mask=true;
	    i++;
            continue;
        }

	if(string(argv[i]) == "-failsafe"  ){
	  failsafe=true;
	  i++;
	  continue;
        }

        if(string(argv[i]) == "-h"  ){
            hFormat=true;
            continue;
        }

        if(string(argv[i]) == "-q"  ){
            quiet=true;
            continue;
        }

        if(string(argv[i]) == "-minq"  ){
            minQualBase=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-minl"  ){
            minLength=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-fa"  ){
	    genomeFile=string(argv[i+1]);
	    genomeFileB=true;
            i++;
            continue;
        }

        if(string(argv[i]) == "-cpg"  ){
	    cpg=true;
            continue;
        }

        if(string(argv[i]) == "-length"  ){
            lengthMaxToPrint=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-err"  ){
            errorToRemove=destringify<double>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-5p" ){
	    file5p = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
	    file3p = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-endo" ){
	    endo   = true;
            continue;
        }

        if(string(argv[i]) == "-both" ){
	    //doubleStr=true;

	    allStr           = false;
	    singleStr        = false;
	    doubleStr        = false;
	    singAnddoubleStr = true;
            continue;
        }


        if(string(argv[i]) == "-single" ){

	    allStr    = false;
	    singleStr = true;
	    doubleStr = false;

            continue;
        }

        if(string(argv[i]) == "-double" ){
	    //doubleStr=true;

	    allStr    = false;
	    singleStr = false;
	    doubleStr = true;

            continue;
        }


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }

    if(  endo &&  paired ){
	cerr<<"Error: cannot specify both -endo and -paired"<<endl;
	return 1;
    }

    
    if(  bedF &&  mask ){
	cerr<<"Error: cannot specify both -bed and -mask"<<endl;
	return 1;
    }

    
    if(phred && hFormat){
	cerr<<"Error: cannot specify both -log and -h"<<endl;
	return 1;
    }

    if(dpFormat && hFormat){
	cerr<<"Error: cannot specify both -dp and -h"<<endl;
	return 1;
    }

    if(endo){
	if(singAnddoubleStr){
	    cerr<<"Error: cannot use -singAnddoubleStr with -endo"<<endl;
	    return 1;
	}
	
	if( !singleStr &&
	    !doubleStr ){
	    cerr<<"Error: you have to provide the type of protocol used (single or double) when using endogenous"<<endl;
	    return 1;
	}

    }
    if(!bedfilename.empty()){
	bed = bed_read(bedfilename.c_str()); 
    }

    typesOfDimer5p       = vector< vector<unsigned int> >();
    typesOfDimer3p       = vector< vector<unsigned int> >();
    typesOfDimer5p_cpg   = vector< vector<unsigned int> >();
    typesOfDimer3p_cpg   = vector< vector<unsigned int> >();
    typesOfDimer5p_noncpg= vector< vector<unsigned int> >();
    typesOfDimer3p_noncpg= vector< vector<unsigned int> >();
    
    typesOfDimer5pDouble = vector< vector<unsigned int> >();
    typesOfDimer3pDouble = vector< vector<unsigned int> >();
    typesOfDimer5pSingle = vector< vector<unsigned int> >();
    typesOfDimer3pSingle = vector< vector<unsigned int> >();

    for(int l=0;l<MAXLENGTH;l++){
	//for(int i=0;i<16;i++){
	typesOfDimer5p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5p_cpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p_cpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5p_noncpg.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3p_noncpg.push_back( vector<unsigned int> ( 16,0 ) );


	typesOfDimer5pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pDouble.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer5pSingle.push_back( vector<unsigned int> ( 16,0 ) );
	typesOfDimer3pSingle.push_back( vector<unsigned int> ( 16,0 ) );

	//}
    }
    
    if(genomeFileB){
	genome=new IndexedGenome(genomeFile.c_str());
	cerr<<genomeFile<<" mapped into memory"<<endl;
    }
    
    string bamfiletopen = string( argv[ argc-1 ] );
    // string deambam      = string( argv[ argc-2 ] );
    // string nondeambam   = string( argv[ argc-1 ] );

    samFile  *fp;
    bam1_t    *b;
    bam_hdr_t *h;

    string refFromFasta_;
    string refFromFasta;
    
    fp = sam_open_format(bamfiletopen.c_str(), "r", NULL); 
    if(fp == NULL){
	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
	return 1;
    }

    h = sam_hdr_read(fp);
    if(h == NULL){
        cerr<<"Could not read header for "<<bamfiletopen<<endl;
        return 1;
    }
    b = bam_init1();
    pair< kstring_t *, vector<int> >  reconstructedReference;
    reconstructedReference.first =(kstring_t *) calloc(sizeof(kstring_t),1);
    reconstructedReference.first->s =0 ;
    reconstructedReference.first->l =    reconstructedReference.first->m =0;
    while(sam_read1(fp, h, b) >= 0){
	if(bam_is_unmapped(b) ){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" unmapped"<<endl;
	    continue;
	}
	if(bam_is_failed(b) ){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" failed"<<endl;
	    continue;
	}
	if(b->core.l_qseq < minLength){
	    if(!quiet)
		cerr<<"skipping "<<bam_get_qname(b)<<" too short"<<endl;
	    continue;
	}
	bool ispaired    = bam_is_paired(b);
	bool isfirstpair = bam_is_read1(b);
	
	if(!paired){    
	    if( ispaired   ){
		if(!quiet)
		    cerr<<"skipping "<<bam_get_qname(b)<<" is paired (can be considered using the -paired flag"<<endl;
		continue;
	    }
	}
	//cerr<<"keeping "<<bam_get_qname(b)<<endl;
	//cout<<bam_get_qname(b)<<endl;
	// 	#define bam_is_paired(b)    (((b)->core.flag&BAM_FPAIRED) != 0)
	// #define bam_is_qcfailed(b)  (((b)->core.flag&BAM_FQCFAIL)     != 0)
	// #define bam_is_rmdup(b)     (((b)->core.flag&BAM_FDUP)        != 0)
	// #define bam_is_failed(b)    ( bam_is_qcfailed(b) || bam_is_rmdup(b) )
	// #define bam_mqual(b)        ((b)->core.qual)

	
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
		continue;
	
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
	
	increaseCounters(b,reconstructedReference.first->s, reconstructedReference.second,minQualBase,refFromFasta,h,bed,mask,ispaired,isfirstpair); //start cycle numberOfCycles-1
    }
    
    bam_destroy1(b);
    sam_close(fp);
    
    ofstream file5pFP;
    file5pFP.open(file5p.c_str());

    if (!file5pFP.is_open()){
	cerr << "Unable to write to 5p file "<<file5p<<endl;
	exit(1);
    }




    //cout<<"cycle\tmatches\tmismatches\tmismatches%\tA>C\tA>C%\tA>G\tA>G%\tA>T\tA>T%\tC>A\tC>A%\tC>G\tC>G%\tC>T\tC>T%\tG>A\tG>A%\tG>C\tG>C%\tG>T\tG>T%\tT>A\tT>A%\tT>C\tT>C%\tT>G\tT>G%"<<endl;
    if(dpFormat)
	file5pFP<<"\t";
    if(hFormat)
	file5pFP<<"pos\t";
    file5pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
  

    vector< vector<unsigned int> > * typesOfDimer5pToUse;

    if(endo){
	if(doubleStr)
	    typesOfDimer5pToUse = &typesOfDimer5pDouble;
	else
	    typesOfDimer5pToUse = &typesOfDimer5pSingle;
    }else{
	typesOfDimer5pToUse     = &typesOfDimer5p;
    }
    
    if(genomeFileB){
	if(cpg)
	    typesOfDimer5pToUse     = &typesOfDimer5p_cpg;
	else
	    typesOfDimer5pToUse     = &typesOfDimer5p_noncpg;
    }
    
    for(int l=0;l<lengthMaxToPrint;l++){
	if(dpFormat)
	    file5pFP<<l<<"\t";

	if(hFormat)
	    file5pFP<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";
	
	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer5pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
		      file5pFP<<dbl2log(MAX(0.0,returnRatioFS( (*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]";
		    else
			if(hFormat)
			  file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS( (*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5);
			else
			  file5pFP<<dbl2log( MAX(0.0,returnRatioFS( (*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred);  
		}else{ 
		    if(singAnddoubleStr){
			if(         (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
				file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]";
			    else
				if(hFormat)
				  file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),1,5);
				else
				  file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
			} else { 
			    if(dpFormat)
				file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
			    else
				if(hFormat)
				    file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				else
				    file5pFP<<(phred?"-Inf":"0.0")<<"";
			}

		    }else{
			if(doubleStr){
			    //          C        T
			    if(         n1==1 && n2==3  ) { 
				if(dpFormat)
				    file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]";
				else
				    if(hFormat)
					file5pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5); 
				    else
					file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
			    } else { 
				if(dpFormat)
				    file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				else
				    if(hFormat)
					file5pFP<<printDoubleAsWhitePaddedString( 0.0,1,5);
				    else					
					file5pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					  file5pFP<<printDoubleAsWhitePaddedString(MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),1,5);
					else
					    file5pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer5pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
					    
				} else { 
				    if(dpFormat)
					file5pFP<<(phred?"-Inf":"0.0")<<" [0..0]";
				    else
					if(hFormat)
					    file5pFP<<printDoubleAsWhitePaddedString(0.0,1,5);
					else										    
					    file5pFP<<(phred?"-Inf":"0.0"); 
				}
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file5pFP<<"\t";
	    }


	}
	file5pFP<<endl;
    }


    file5pFP.close();

    ofstream file3pFP;
    if(file3p == "/dev/stdout"){
	file3pFP.open(file3p.c_str(), ofstream::out | ofstream::app);
    }else{
	file3pFP.open(file3p.c_str());
    }

    if (!file3pFP.is_open()){
	cerr << "Unable to write to 3p file "<<file3p<<endl;
	exit(1);
    }

    if(dpFormat)
	file3pFP<<"\t";
    if(hFormat)
	file3pFP<<"pos\t";

    file3pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;


    vector< vector<unsigned int> > * typesOfDimer3pToUse;

    if(endo){
	if(doubleStr)
	    typesOfDimer3pToUse = &typesOfDimer3pDouble;
	else
	    typesOfDimer3pToUse = &typesOfDimer3pSingle;
    }else{
	typesOfDimer3pToUse     = &typesOfDimer3p;
    }

    if(genomeFileB){
	if(cpg)
	    typesOfDimer3pToUse     = &typesOfDimer3p_cpg;
	else
	    typesOfDimer3pToUse     = &typesOfDimer3p_noncpg;
    }
    

    for(int le=0;le<lengthMaxToPrint;le++){

	int l=le;
	if(dpFormat){
	    l=lengthMaxToPrint-1-le;
	    if(l==0)
		file3pFP<<""<<l<<"\t";	    
	    else
		file3pFP<<"-"<<l<<"\t";	    
	}

	if(hFormat){
	    //l=lengthMaxToPrint-1-le;
	    // if(l==0)
	    // 	file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    // else
	    // 	file3pFP<<"-"<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint)))<<"\t";	    
	    file3pFP<<""<<printIntAsWhitePaddedString(l,int(log10(lengthMaxToPrint))+1)<<"\t";	    
	}

	for(int n1=0;n1<4;n1++){   
	    int totalObs=0;
	    for(int n2=0;n2<4;n2++){   
		totalObs+=(*typesOfDimer3pToUse)[l][4*n1+n2];
	    }

	    for(int n2=0;n2<4;n2++){   
		if(n1==n2)
		    continue;
		if(allStr){
		    if(dpFormat)
		      file3pFP<<dbl2log( MAX(0.0,returnRatioFS( (*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]";
		    else
			if(hFormat)
			  file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS( (*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5);
			else
			  file3pFP<<dbl2log( MAX(0.0,returnRatioFS( (*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred);
		}else{ 
		    if(singAnddoubleStr){			
			if(   (n1==1 && n2==3) || (n1==2 && n2==0 )  ) { 
			    if(dpFormat)
			      file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]"; 
			    else
				if(hFormat)
				  file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5);
				else
				  file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
			} else { 
			    if(dpFormat)
				file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
			    else
				if(hFormat)
				    file3pFP<<printDoubleAsWhitePaddedString( 0.0 ,1,5);
				else
				    file3pFP<<(phred?"-Inf":"0.0"); 
			}

		    }else{
			if(doubleStr){
			    //          G        A
			    if(         n1==2 && n2==0  ) { 
				if(dpFormat)
				    file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5);
				    else					
					file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
			    } else { 				
				if(dpFormat)
				    file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				else
				    if(hFormat)
					file3pFP<<printDoubleAsWhitePaddedString(  0.0 ,1,5);
				    else	
					file3pFP<<(phred?"-Inf":"0.0"); 
			    }
			}else{ 
			    if(singleStr){
				//      C        T
				if(     n1==1 && n2==3  ) { 
				    if(dpFormat)
					file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred)<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)) ,1,5);
					else					
					    file3pFP<<dbl2log( MAX(0.0,returnRatioFS((*typesOfDimer3pToUse)[l][4*n1+n2],totalObs,errorToRemove,failsafe)),phred); 
				} else { 				    
				    if(dpFormat)
					file3pFP<<(phred?"-Inf":"0.0")<<" [0..0]"; 
				    else
					if(hFormat)
					    file3pFP<<printDoubleAsWhitePaddedString(  0.0 , 1,5);
					else	
					    file3pFP<<(phred?"-Inf":"0.0"); 

				}
				
			    }
			}
		    }
		}

		
		if(!(n1 ==3 && n2 == 2 ))
		    file3pFP<<"\t";
	    }


	}
	file3pFP<<endl;
    }


    //cerr<<"TEST1"<<endl;
    file3pFP.close();
    ///cerr<<"TEST2"<<endl;
    if (bed) bed_destroy(bed);
    //cerr<<"TEST3"<<endl;
    return 0;
}

