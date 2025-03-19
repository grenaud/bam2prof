// utility.h
#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <cstring>
#include <set>
#include <ctype.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <fstream>
#include <cstdio> 

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

// louis was here, get the ref id
#define bam_ref_id(b) 		((b)->core.tid)

#define MAXLENGTH 1000

using namespace std;

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



// Function declarations
double returnRatioFS(int num,int denom,double errorToRemove,bool failsafe=false);
			    
//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(const bam1_t  * b, char *reconstructedReference, const vector<int> &  reconstructedReferencePos, const int & minQualBase, string & refFromFasta, const bam_hdr_t *h, void *bed,bool mask, bool ispaired, bool isfirstpair, std::vector<std::vector<unsigned int>>& typesOfDimer5p, std::vector<std::vector<unsigned int>>& typesOfDimer3p, std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle);

double dbl2log(const double d,bool phred);

void countSubsPerRef(bool genomeFileB, IndexedGenome* genome, const bam1_t  * b, std::pair<kstring_t*, std::vector<int>>& reconstructedReference, const int & minQualBase, string & refFromFasta, string & refFromFasta_, const bam_hdr_t *h, void *bed,bool mask, bool ispaired, bool isfirstpair, std::vector<std::vector<unsigned int>>& typesOfDimer5p, std::vector<std::vector<unsigned int>>& typesOfDimer3p, std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle);

vector<vector<unsigned int>> initializeDimerVectors(int maxLength, int innerSize);

void generateDamageProfile( const std::string& outDir, const std::string& bamfiletopen, const std::string& refId, int lengthMaxToPrint, bool dpFormat, bool hFormat, bool allStr, bool singAnddoubleStr, bool doubleStr, bool singleStr, bool endo, bool genomeFileB, bool cpg, double errorToRemove, bool failsafe, bool phred, const std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer5p, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer3p, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, uint64_t mapped);

#endif // UTILITY_H
