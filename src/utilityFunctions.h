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
    /* size in bytes of the fasta file */
    int64_t fileSize;
    /* if true map the fasta a window of windowBases bases at a time, else map all of it on first use */
    bool windowed;
    int64_t windowBases;
    /* bases kept mapped before the requested coordinate when a new window is mapped, so reads that
       step slightly backwards do not trigger a remap */
    static const int64_t windowMargin=1000000;

    /* a memory mapped stretch of the fasta file: file bytes [lo,hi) are at ptr (lo is page-aligned),
       which hold bases [startBase,endBase) of the chromosome whose faidx offset is chrOffset */
    struct Window {
	char *ptr;
	int64_t lo, hi;
	uint64_t chrOffset;
	int64_t startBase, endBase;
	Window():ptr(NULL),lo(0),hi(0),chrOffset(0),startBase(0),endBase(0){}
	bool covers(int64_t first, int64_t last) const { return ptr!=NULL && first>=lo && last<hi; }
    };
    /* the window being read from, and the one prefetched after it (windowed mode only) */
    Window cur, nxt;

    /* file byte offset of base 'index' (0-based) of the chromosome described by faidx */
    static int64_t byteOffset(const faidx1_t * faidx, int64_t index){
	return (int64_t)faidx->offset + index / faidx->line_blen * faidx->line_len + index % faidx->line_blen;
    }

    static void unmapWindow(Window & w){
	if(w.ptr!=NULL && munmap(w.ptr,w.hi-w.lo) == -1){
	    perror("Error un-mmapping the file");
	}
	w=Window();
    }

    void unmapAll(){
	unmapWindow(cur);
	unmapWindow(nxt);
    }

    /* maps file bytes [lo,hi) (lo is rounded down to a page boundary) into w; if prefetch is set the
       kernel is asked to start reading the pages in right away */
    void mapBytes(Window & w, int64_t lo, int64_t hi, bool prefetch){
	int64_t page=sysconf(_SC_PAGESIZE);
	lo=(lo/page)*page;
	if(hi>fileSize) hi=fileSize;
	char * p = (char*)mmap(0, hi-lo, PROT_READ, MAP_SHARED, fd, lo);
	if (p == MAP_FAILED)
	    {
		close(fd);
		perror("Error mmapping the file");
		exit(EXIT_FAILURE);
	    }
	if(prefetch) madvise(p, hi-lo, MADV_WILLNEED);
	w.ptr=p;
	w.lo=lo;
	w.hi=hi;
    }

    /* maps bases [startBase,endBase) of the chromosome into w */
    void mapBases(Window & w, const faidx1_t * faidx, int64_t startBase, int64_t endBase, bool prefetch){
	mapBytes(w, byteOffset(faidx,startBase), byteOffset(faidx,endBase-1)+1, prefetch);
	w.chrOffset=faidx->offset;
	w.startBase=startBase;
	w.endBase=endBase;
    }

    /* makes sure the bases [index,index+length) of the chromosome are mapped and returns the window
       holding them. In windowed mode, once reads reach the middle of the current window the following
       window is mapped (and prefetched) so it is ready by the time reads cross into it; when they do, the
       old window is released, so at most two windows are mapped at once. */
    const Window * ensureMapped(const faidx1_t * faidx, int64_t index, unsigned int length){
	int64_t first=byteOffset(faidx,index);
	int64_t last =byteOffset(faidx,index+length-1);
	if(!windowed){
	    if(!cur.ptr) mapBytes(cur,0,fileSize,false);
	    return &cur;
	}
	if(!cur.covers(first,last)){
	    if(nxt.covers(first,last)){ //crossed into the prefetched window
		unmapWindow(cur);
		cur=nxt;
		nxt=Window();
	    }else{ //first lookup, or a jump (new chromosome, unsorted reads): start over around the request
		unmapAll();
		int64_t startBase=index-windowMargin;
		if(startBase<0) startBase=0;
		int64_t endBase=startBase+windowBases;
		if(endBase<index+(int64_t)length) endBase=index+length;
		if(endBase>faidx->len) endBase=faidx->len;
		mapBases(cur,faidx,startBase,endBase,false);
	    }
	}
	//halfway through the current window: bring in the next one
	if(nxt.ptr==NULL && cur.chrOffset==faidx->offset && cur.endBase<faidx->len &&
	   index >= cur.startBase + (cur.endBase-cur.startBase)/2){
	    int64_t endBase=cur.endBase+windowBases;
	    if(endBase>faidx->len) endBase=faidx->len;
	    mapBases(nxt,faidx,cur.endBase,endBase,true);
	}
	return &cur;
    }
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
    IndexedGenome(const char* fasta):fd(-1),fileSize(0),windowed(false),windowBases(10000000)
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
	fileSize=buf.st_size;
	/* the fasta is memory mapped lazily, on the first lookup, see setWindowed() */
    }

    /** chooses how the fasta is memory mapped. Call before the first lookup.
     * @param sorted: true if the reads come in coordinate order, the reference is then mapped
     *                windowBases bases at a time and the previous window is released; false maps the whole fasta
     * @param bases: window size in bases
     */
    void setWindowed(bool sorted, int64_t bases=10000000)
    {
	unmapAll();
	windowed=sorted;
	windowBases=bases;
    }
    /* destructor */
    ~IndexedGenome()
    {
	/* close memory mapped maps */
	unmapAll();
	/* dispose fasta file descriptor */
	if(fd!=-1) close(fd);
    }
			

    /* return the base at position 'index' for the chromosome indexed by faidx */
    string returnStringCoord(const FaidxPtr faidx,int64_t index, unsigned int length){

	if(length==0) return "";
	if(byteOffset(faidx,index+length-1) >= fileSize){
	    return string(length,'N'); //past the end of the file, nothing to map
	}
	const Window * w=ensureMapped(faidx,index,length);

	string strToReturn="";
	strToReturn.reserve(length);
	for(unsigned int j=0;j<length;j++){ //for each char
	    int64_t pos=byteOffset(faidx,index+j);
	    strToReturn+=char(toupper(w->ptr[pos-w->lo]));
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

//counts the reference nucleotide composition flanking a read's aligned span (mapDamage-style), needs -fa
void countBaseCompositionFlank(IndexedGenome* genome, const bam1_t * b, const bam_hdr_t *h, int aroundFlank, bool ispaired, bool isfirstpair, std::vector<std::vector<unsigned int>>& baseComp5pFlank, std::vector<std::vector<unsigned int>>& baseComp3pFlank);

void generateBaseCompositionProfile( const std::string& outDir,
				      const std::string& bamfiletopen,
				      const std::string& refId,
				      int lengthMaxToPrint,
				      int aroundFlank,
				      bool genomeFileB,
				      const std::vector<std::vector<unsigned int>>& typesOfDimer5p,
				      const std::vector<std::vector<unsigned int>>& typesOfDimer3p,
				      const std::vector<std::vector<unsigned int>>& baseComp5pFlank,
				      const std::vector<std::vector<unsigned int>>& baseComp3pFlank,
				      uint64_t mapped);

void generateDamageProfile( const std::string& outDir,
			    const std::string& file5pparam,
			    const std::string& file3pparam,
			    const std::string& bamfiletopen, const std::string& refId, int lengthMaxToPrint, bool dpFormat, bool hFormat, bool allStr, bool singAnddoubleStr, bool doubleStr, bool singleStr, bool endo, bool genomeFileB, bool cpg, double errorToRemove, bool failsafe, bool phred, const std::vector<std::vector<unsigned int>>& typesOfDimer5pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer5pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer5p, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer5p_noncpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3pSingle, const std::vector<std::vector<unsigned int>>& typesOfDimer3pDouble, const std::vector<std::vector<unsigned int>>& typesOfDimer3p, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_cpg, const std::vector<std::vector<unsigned int>>& typesOfDimer3p_noncpg, uint64_t mapped);

#endif // UTILITY_H
