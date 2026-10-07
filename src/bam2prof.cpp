#include "utilityFunctions.h" // Include the utility header

#define OPENMP
#ifdef OPENMP
#include <omp.h>
#endif

// Originally written by Gabriel Renaud & Thorfinn Korneliussen
// Modified and extended by Louis Kraft


int main (int argc, char *argv[]) {

    string file5pDefault="/dev/stdout";
    string file3pDefault="/dev/stdout";
    string outDir="/dev/stdout";
    string outDirUnsafe="/dev/stdout";
    vector<string> refIdsList={};
    vector<string> bamFiles={};
    
    bool endo=false;

    bool allStr   =true;
    bool singleStr=false;
    bool doubleStr=false;
    bool singAnddoubleStr=false;

    int lengthMaxToPrint = 5;
    int minQualBase      = 0;
    int minLength        = 35;
    int numAlns        = 10000000;
    bool dpFormat=false;
    bool hFormat=false;
    double errorToRemove=0.0;
    bool phred=false;
    string genomeFile;
    bool genomeFileB=false;
    IndexedGenome* genome=NULL;
    bool cpg=false;
    string bedfilename;
    void *bed = 0; // BED data structure

    bool bedF=false;
    bool mask=false;
    bool failsafe=false;
    bool metaMode=false;
    bool paired=false;
    bool quiet=true;
    bool classicMode=true;
    double precisionThresh=0.0;
    double precisionConverge=0.01;
    int stepsizeConverge=500;
    unsigned int convergeUntil=100000;
    //string refId;

    bool isizeB=false;
    bool isizeAllPairs=false; //-is-allpaired: count every read1 of a pair, not only properly paired ones
    string isizeFile;
    map<int32_t,uint64_t> insertSizeCounts;       //all fragments: properly paired + merged/single-end
    map<int32_t,uint64_t> insertSizeCountsPaired; //properly paired fragments only (every pair with -is-allpaired)
    map<int32_t,uint64_t> insertSizeCountsMerged; //merged / single-end molecules only

    bool compFlag=false;
    int aroundFlank=10;
    vector< vector<unsigned int> > baseComp5pFlank; //nt composition upstream  of the read's 5' end (from the reference, needs -fa)
    vector< vector<unsigned int> > baseComp3pFlank; //nt composition downstream of the read's 3' end (from the reference, needs -fa)

	//#define DEBUG

    string usage=string(""+string(argv[0])+" <options> <mode> [in BAM file]"+
			"\nThis program reads a BAM file and produces a deamination profile for the\n"+
			"5' and 3' ends\n"+

			"\n\nPlease provide sorted bam files (by reference name).\n"+

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
			"\t\t"+"-meta\t\t\tOne Profile for each unique reference    (Default: "+booleanAsString( metaMode )+" ) \n"+
			"\t\t"+"-classic\t\tOne Profile per bam file (default)   (Default: "+booleanAsString( classicMode )+" ) \n"+
			"\t\t"+"-precision\t\tSet minimum precision for substitution frequency computation (Default: All alignments [= 0.0]; Speed up by setting precision to either 0.01, 0.001, ... ) \n"+
			"\t\t"+"-minAligned\t\tNumber of aligned sequences after which substitution patterns are checked for converging (Default: "+stringify( numAlns )+")\n"+
			"\t\t"+"-ref-id\t\t\tSpecify reference ID; if multiple references: Provide comma seperated list (no spaces!) ( Default: Not Set ) \n"+
			"\t\t"+"-is\t[output file]\tAlso compute the fragment size distribution over the same reads used for the profile and write it there as \"count length\" per line, sorted by length. Properly paired fragments are reported as abs(TLEN) of read1 (requires -paired), merged / single-end molecules as their SEQ length. The same counts are also written separately to \"<file>.properly_paired\" and \"<file>.merged\" (Default: not computed) \n"+
			"\t\t"+"-is-allpaired\t\tWith -is, count every paired read1 as a fragment, not only properly paired ones; the paired counts are then written to \"<file>.paired\" instead of \"<file>.properly_paired\" (Default: "+booleanAsString( isizeAllPairs )+")\n"+
			"\t\t"+"-comp\t\t\tAlso compute a base composition profile (A/C/G/T frequency per position) written as _5p_comp.prof/_3p_comp.prof next to the substitution profiles. Without -fa, only positions inside the fragment are reported; with -fa, it also reports "+stringify(aroundFlank)+" bp of reference sequence flanking the fragment on either side  (Default: "+booleanAsString( compFlag )+") \n"+
			"\t\t"+"-around\t[N]\t\tNumber of reference bp to report outside the fragment for -comp; only used together with -fa (Default: "+stringify(aroundFlank)+")\n"+

			"\n\n\tYou can specify either one of the two:\n"+
			"\t\t"+"-single\t\t\tUse the deamination profile of a single strand library  (Default: "+booleanAsString( singleStr )+")\n"+
			"\t\t"+"-double\t\t\tUse the deamination profile of a double strand library  (Default: "+booleanAsString( doubleStr )+")\n"+
			"\n\tor specify this option:\n"+
			"\t\t"+"-both\t\t\tReport both C->T and G->A regardless of stand  (Default: "+booleanAsString( singAnddoubleStr )+")\n"+

			"\n\n\tOutput options:\n"+
			"\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5pDefault)+")\n"+
			"\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3pDefault)+")\n"+
			"\t\t"+"-o\t[output dir]\tOutput Directory for all matrices (Default: "+stringify(outDir)+")\n"+			
			"\t\t"+"-dp\t\t\tOutput in damage-patterns format (Default: "+booleanAsString(dpFormat)+")\n"+
			"\t\t"+"-h\t\t\tMore human readible output (Default: "+booleanAsString(hFormat)+")\n"+
			"\t\t"+"-q\t\t\tDo not print why reads are skipped. Turn On [1] or Off [0] (Default: "+booleanAsString(quiet)+")\n"+

			"\n\n\tExpert options (Identifiying zig-zag profiles):\n"+
			"\t\t"+"-minConverge\t\tSet minimum threshold for substitution frequency convergence (Default: "+stringify( precisionConverge)+")\n"+
			"\t\t"+"-stepsConverge\t\tCheck every step-size number of aligned fragments for convergence (Default: "+stringify( stepsizeConverge )+")\n"+
		       
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

	if(string(argv[i]) == "-meta"  ){
	    metaMode=true;
	    continue;
        }

	if(string(argv[i]) == "-classic"  ){
	    classicMode=true;
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
            quiet=destringify<int>(argv[i+1]);
			i++;
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

		if(string(argv[i]) == "-precision"  ){
            precisionThresh=destringify<double>(argv[i+1]);
			i++;
            continue;
        }

        if(string(argv[i]) == "-minAligned"  ){
            numAlns=destringify<int>(argv[i+1]);
            i++;
            continue;
        }


		if(string(argv[i]) == "-minConverge"  ){
            precisionConverge=destringify<double>(argv[i+1]);
			i++;
            continue;
        }

        if(string(argv[i]) == "-stepsConverge"  ){
            stepsizeConverge=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-convergeUntil"  ){
            convergeUntil=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-is-allpaired" ){
	    isizeAllPairs = true;
            continue;
        }

        if(string(argv[i]) == "-is" ){
	    isizeFile = string(argv[i+1]);
	    isizeB    = true;
	    i++;
            continue;
        }

        if(string(argv[i]) == "-comp" ){
	    compFlag=true;
            continue;
        }

        if(string(argv[i]) == "-around" ){
	    aroundFlank=destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-ref-id" ){
            string refIdList = string(argv[i+1]);  // Capture the comma-separated list
            stringstream ss(refIdList);
            string item;

            while (std::getline(ss, item, ',')) {
                refIdsList.push_back(item);  // Add each ref ID to the vector
            }

	    	i++;
            continue;
        }

        if(string(argv[i]) == "-5p" ){
	    file5pDefault = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
	    file3pDefault = string(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "-o" ){
	    outDir = string(argv[i+1]);
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

    if ( classicMode && metaMode ){
	std::cerr << "Error: cannot specify both -classic and -meta )" << std::endl;
	return 1;
    }

    if ( !classicMode && !metaMode ){
	std::cerr << "Error: need to specify mode: -classic or -meta )" << std::endl;
	return 1;
    }

    if ( classicMode && (refIdsList.size() != 0 || metaMode) ){
	std::cerr << "Error: cannot specify none: At least one must be specified -classic OR -meta )" << std::endl;
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

    if(compFlag && aroundFlank<0){
	cerr<<"Error: -around cannot be negative"<<endl;
	return 1;
    }

    if(compFlag && !genomeFileB){
	cerr<<"Warning: -comp without -fa will only report base composition inside the fragment, not the flanking reference bases"<<endl;
    }

    if(genomeFileB){
	genome=new IndexedGenome(genomeFile.c_str());
    }
    
    // Define the output path for the profiles that have not converged
    if (!outDir.empty() && outDir.back() == '/'){
	outDir.pop_back();
    } 
    
    outDirUnsafe = outDir + "_notConverged";
    string outDirSwap = outDir;
    
    // Creating the output directory
    if ( outDir != "/dev/stdout" ) {
        std::string command = "mkdir -p " + outDir;
        int result = system(command.c_str());
        if (result != 0) {
            std::cerr << "Failed to create output directories: " << outDir << std::endl;
        }
    }
    if ( metaMode && outDir != "/dev/stdout") {
        std::string command2 = "mkdir -p " + outDirUnsafe;
        int result2 = system(command2.c_str());
        if (result2 != 0) {
            std::cerr << "Failed to create output directories: " << outDirUnsafe << std::endl;
        }
    }

    // string bamfilelist = string( argv[ argc-1 ] );
    // stringstream ss(bamfilelist);
    // string bamitem;
    // while (std::getline(ss, bamitem, ',')) {
    // 	bamFiles.push_back(bamitem);  // Add each ref ID to the vector
    // }
    
	
    string bamPath = string( argv[ argc-1 ] );
    bool fromStdin = (bamPath == "-" || bamPath == "/dev/stdin"); //piped input: "-" is htslib's name for stdin
    string bamfiletopen = fromStdin ? string("stdin") : bamPath;   //used to name the output files
    
    bam_hdr_t *h;
    samFile  *fp;
    hts_idx_t *idx;
    
    fp = sam_open_format(bamPath.c_str(), "r", NULL); 
    if(fp == NULL){
	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
	return 1;
    }

    h = sam_hdr_read(fp);
    if(h == NULL){
	cerr<<"Could not read header for "<<bamfiletopen<<endl;
	return 1;
    }

    if(genomeFileB){
	//coordinate sorted BAM: the reads only move forward along each chromosome, so map the reference
	//10Mb at a time instead of all of it; otherwise map the whole fasta
	kstring_t so = KS_INITIALIZE;
	bool sorted = (sam_hdr_find_tag_hd(h, "SO", &so) == 0 && so.s != NULL && string(so.s) == "coordinate");
	ks_free(&so);
	genome->setWindowed(sorted);
	if(sorted){
	    cerr<<genomeFile<<": BAM is coordinate sorted, reference will be memory mapped 10Mb at a time"<<endl;
	}else{
	    cerr<<genomeFile<<": BAM is not flagged as coordinate sorted, whole reference will be memory mapped"<<endl;
	}
    }

    // Load the index for the BAM file
    // Without an index (piped input, or a BAM that is not coordinate sorted) the reads are processed
    // front to back in one pass instead of one chromosome at a time
    idx = fromStdin ? NULL : sam_index_load(fp, bamPath.c_str());
    bool streaming = (idx == NULL);
    if(streaming){
	if(!classicMode){
	    std::cerr << "Error: -meta and -ref need an indexed, coordinate sorted BAM file; without an index only -classic can be used" << std::endl;
	    return 1;
	}
	std::cerr << (fromStdin ? "Reading from stdin" : "No index found for "+bamPath) << ", processing the reads sequentially" << std::endl;
    }
    
    
    //std::set<int32_t> refIdSet;
    std::set<std::string> refNameSet;
    
    for (const auto& refName : refIdsList) {
	refNameSet.insert(refName);
	//refIdSet.insert(sam_hdr_name2tid(h, refName.c_str()));
    }


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
    
    // Then we initialize a new vectors to count:
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

    baseComp5pFlank = vector< vector<unsigned int> >();
    baseComp3pFlank = vector< vector<unsigned int> >();
    for(int l=0;l<aroundFlank;l++){
	baseComp5pFlank.push_back( vector<unsigned int> ( 4,0 ) );
	baseComp3pFlank.push_back( vector<unsigned int> ( 4,0 ) );
    }

    // Initiating the early stop rule for the "classic" mode:

    vector< vector<unsigned int> > typesOfDimer5pTmp; //5' deam rates
    vector< vector<unsigned int> > typesOfDimer3pTmp; //3' deam rates

    vector< vector<unsigned int> > typesOfDimer5p_cpgTmp; //5' deam rates
    vector< vector<unsigned int> > typesOfDimer3p_cpgTmp; //3' deam rates
    vector< vector<unsigned int> > typesOfDimer5p_noncpgTmp; //5' deam rates
    vector< vector<unsigned int> > typesOfDimer3p_noncpgTmp; //3' deam rates

    vector< vector<unsigned int> > typesOfDimer5pDoubleTmp; //5' deam rates when the 3' is deaminated according to a double str.
    vector< vector<unsigned int> > typesOfDimer3pDoubleTmp; //3' deam rates when the 5' is deaminated according to a double str.
    vector< vector<unsigned int> > typesOfDimer5pSingleTmp; //5' deam rates when the 3' is deaminated according to a single str.
    vector< vector<unsigned int> > typesOfDimer3pSingleTmp; //3' deam rates when the 5' is deaminated according to a single str.
    
    // Then we initialize a new vectors to count:
    typesOfDimer5pTmp       = typesOfDimer5p;
    typesOfDimer3pTmp       = typesOfDimer3p;
    typesOfDimer5p_cpgTmp   = typesOfDimer5p_cpg;
    typesOfDimer3p_cpgTmp   = typesOfDimer3p_cpg;
    typesOfDimer5p_noncpgTmp= typesOfDimer5p_noncpg;
    typesOfDimer3p_noncpgTmp= typesOfDimer3p_noncpg;
    
    typesOfDimer5pDoubleTmp = typesOfDimer5pDoubleTmp;
    typesOfDimer3pDoubleTmp = typesOfDimer3pDoubleTmp;
    typesOfDimer5pSingleTmp = typesOfDimer5pSingleTmp;
    typesOfDimer3pSingleTmp = typesOfDimer3pSingleTmp;
    
    uint64_u totalMapped = 0;
    
    bool stopEarly = false;
    bool isConvergedOnce = false;

    unsigned int numAlnsSafeRange = convergeUntil;
    unsigned int numAlnsSafe = stepsizeConverge;
    float critThresh = precisionConverge;

    // Iterate over each reference 
    // streaming: a single pass over every read, in file order, instead of one iterator per chromosome
    for (int i = 0; i < (streaming ? 1 : h->n_targets); i++) {

	unsigned int processedAlns = 1;

	//std::cerr << "iterating targets" << std::endl;
	const char* refName = streaming ? "classic" : h->target_name[i];
	std::string refNameStr(refName);

	if ( !classicMode && refIdsList.size() > 0 ){
	    if (refNameSet.find(refNameStr) == refNameSet.end()) {
		//std::cerr << "inside continue" << std::endl;
		continue;
	    }
	}
	
	if ( !classicMode ){
	    // Then we initialize a new vector to count:
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

	    baseComp5pFlank = vector< vector<unsigned int> >();
	    baseComp3pFlank = vector< vector<unsigned int> >();
	    for(int l=0;l<aroundFlank;l++){
		baseComp5pFlank.push_back( vector<unsigned int> ( 4,0 ) );
		baseComp3pFlank.push_back( vector<unsigned int> ( 4,0 ) );
	    }
	}

	hts_itr_t *iter = streaming ? NULL : sam_itr_queryi(idx, i, 0, h->target_len[i]);
	bam1_t *b = bam_init1();
	
	// Check if iter is null
	if (!streaming && iter == NULL) {
	    std::cerr << "Could not create iterator for target " << h->target_name[i] << std::endl;
	    continue;
	}

	string refFromFasta_;
	string refFromFasta;
	
	pair< kstring_t *, vector<int> >  reconstructedReference;
	reconstructedReference.first =(kstring_t *) calloc(sizeof(kstring_t),1);
	reconstructedReference.first->s =0 ;
	reconstructedReference.first->l =    reconstructedReference.first->m =0;
	

	// Get the number of mapped reads for the current reference 'i'
	uint64_t mapped = 0, unmapped = 0;
	if(!streaming){
	    hts_idx_get_stat(idx, i, &mapped, &unmapped);
	    totalMapped += mapped;
	}
	
	bool isConvergedLow = false;  // Flag to indicate if all changes are small; Also tells us if converged at all or not
	
	stopEarly = false; // for any mode but references where at least 10 Mio have aligned PER REFERENCE
		
	// Iterate over the BAM records
	while ( (streaming ? sam_read1(fp, h, b) : sam_itr_next(fp, iter, b)) >= 0) {
	    
	    if ( processedAlns >= 10000 ){
		unsigned int numAlnsSafe = 1000;
		float critThresh = precisionConverge;
	    }
	    
	    bool critRange = (processedAlns <= numAlnsSafeRange);
	    if(bam_is_unmapped(b)){
		if(!quiet)
		    cerr<<"skipping "<<bam_get_qname(b)<<" unmapped"<<endl;
		continue;
	    }
	    if(streaming) totalMapped++; //same count the index reports: every read not flagged unmapped
	    if(bam_is_failed(b)){
		if(!quiet)
		    cerr<<"skipping "<<bam_get_qname(b)<<" failed"<<endl;
		continue;
	    }
	    if(b->core.l_qseq < minLength){
		if(!quiet)
		    cerr<<"skipping "<<bam_get_qname(b)<<" too short"<<endl;
		continue;
	    }
	    bool ispaired = bam_is_paired(b);
	    bool isfirstpair = bam_is_read1(b);
	    
	    if(!paired){    
		if(ispaired){
		    if(!quiet)
			cerr<<"skipping "<<bam_get_qname(b)<<" is paired (can be considered using the -paired flag)"<<endl;
		    continue;
		}
	    }
	    
	    // fragment size distribution (mirrors the standalone insertsize/insize tool)
	    if(isizeB){
		if(ispaired){
		    if(isfirstpair && (isizeAllPairs || (b->core.flag & BAM_FPROPER_PAIR))){ //properly paired fragments only, unless -is-allpaired
			int32_t isize = b->core.isize;
			if(isize != 0){ //skip mates on different contigs / unmapped mates
			    insertSizeCounts[ abs(isize) ]++;
			    insertSizeCountsPaired[ abs(isize) ]++;
			}
		    }//read2 is skipped: each fragment is counted once, via read1
		}else{ //merged (collapsed) or single-end molecule: the read is the whole fragment
		    insertSizeCounts[ b->core.l_qseq ]++;
		    insertSizeCountsMerged[ b->core.l_qseq ]++;
		}
	    }

	    // base composition of the reference flanking the fragment ( needs -fa)
	    if(compFlag && genomeFileB){
		countBaseCompositionFlank(genome, b, h, aroundFlank, ispaired, isfirstpair, baseComp5pFlank, baseComp3pFlank);
	    }

	    // update matrix
	    countSubsPerRef(genomeFileB, genome, b, reconstructedReference, minQualBase, refFromFasta, refFromFasta_, h, bed, mask, ispaired, isfirstpair, typesOfDimer5p, typesOfDimer3p, typesOfDimer5p_cpg, typesOfDimer3p_cpg, typesOfDimer5p_noncpg, typesOfDimer3p_noncpg, typesOfDimer5pDouble, typesOfDimer3pDouble, typesOfDimer5pSingle, typesOfDimer3pSingle);

	    
	    // Early stop mechanism: Check every numAlns alignments
	    if (processedAlns % numAlns == 0 || ( critRange && processedAlns % numAlnsSafe == 0 && metaMode ) ) {
		
		
		// std::cerr << "critRange\t" << critRange << std::endl;
		// std::cerr << "isConvergedLow\t" << isConvergedLow << std::endl;
		// std::cerr << "processedAlns\t" << processedAlns << std::endl;
		// std::cerr << "refNameStr\t" << refNameStr << std::endl;
		// std::cerr << "\n" << std::endl;
		
		// Initialize the vector to store changes
		std::vector<std::vector<double>> difference5p(MAXLENGTH, std::vector<double>(16, 0.0));
		std::vector<std::vector<double>> difference3p(MAXLENGTH, std::vector<double>(16, 0.0));
		
		// Choose the appropriate 5' dimer type for current and temporary matrices
		const std::vector<std::vector<unsigned int>>* dimer5pToUse;
		const std::vector<std::vector<unsigned int>>* dimer5pTmpToUse;
		if (endo) {
		    dimer5pToUse = doubleStr ? &typesOfDimer5pDouble : &typesOfDimer5pSingle;
		    dimer5pTmpToUse = doubleStr ? &typesOfDimer5pDoubleTmp : &typesOfDimer5pSingleTmp;
		} else {
		    dimer5pToUse = &typesOfDimer5p;
		    dimer5pTmpToUse = &typesOfDimer5pTmp;
		}
		if (genomeFileB) {
		    dimer5pToUse = cpg ? &typesOfDimer5p_cpg : &typesOfDimer5p_noncpg;
		    dimer5pTmpToUse = cpg ? &typesOfDimer5p_cpgTmp : &typesOfDimer5p_noncpgTmp;
		}
		
		// Same for 3'
		const std::vector<std::vector<unsigned int>>* dimer3pToUse;
		const std::vector<std::vector<unsigned int>>* dimer3pTmpToUse;
		if (endo) {
		    dimer3pToUse = doubleStr ? &typesOfDimer3pDouble : &typesOfDimer3pSingle;
		    dimer3pTmpToUse = doubleStr ? &typesOfDimer3pDoubleTmp : &typesOfDimer3pSingleTmp;
		} else {
		    dimer3pToUse = &typesOfDimer3p;
		    dimer3pTmpToUse = &typesOfDimer3pTmp;
		}
		if (genomeFileB) {
		    dimer3pToUse = cpg ? &typesOfDimer3p_cpg : &typesOfDimer3p_noncpg;
		    dimer3pTmpToUse = cpg ? &typesOfDimer3p_cpgTmp : &typesOfDimer3p_noncpgTmp;
		}
		
		
		// Assign based on input arguments (endo, doubleStr, etc.)
		// For brevity, the same logic as previously described can be used here
		
		// Loop through each position
		for (int l = 0; l < MAXLENGTH; ++l) {
		    // Calculate differences for 5' end
		    for (int n1 = 0; n1 < 4; ++n1) {
			int totalObsCurrent5p = 0, totalObsTmp5p = 0;
			
			// Get total observations for each combination
			for (int n2 = 0; n2 < 4; ++n2) {
			    totalObsCurrent5p += (*dimer5pToUse)[l][4 * n1 + n2];
			    totalObsTmp5p += (*dimer5pTmpToUse)[l][4 * n1 + n2];
			}
			
			// Calculate ratio differences for each match/mismatch type
			for (int n2 = 0; n2 < 4; ++n2) {
			    double currentRatio5p = (totalObsCurrent5p == 0) ? static_cast<double>(0) : returnRatioFS((*dimer5pToUse)[l][4 * n1 + n2], totalObsCurrent5p, errorToRemove, failsafe);
			    double tmpRatio5p = (totalObsTmp5p == 0) ? static_cast<double>(0) : returnRatioFS((*dimer5pTmpToUse)[l][4 * n1 + n2], totalObsTmp5p, errorToRemove, failsafe);
			    
			    // Store the absolute difference in ratios
			    difference5p[l][4 * n1 + n2] = fabs(currentRatio5p - tmpRatio5p);
			}
		    }
		    
		    // Repeat the same process for the 3' end
		    for (int n1 = 0; n1 < 4; ++n1) {
			int totalObsCurrent3p = 0, totalObsTmp3p = 0;
			
			// Get total observations for each combination
			for (int n2 = 0; n2 < 4; ++n2) {
			    totalObsCurrent3p += (*dimer3pToUse)[l][4 * n1 + n2];
			    totalObsTmp3p += (*dimer3pTmpToUse)[l][4 * n1 + n2];
			}
			
			// Calculate ratio differences for each match/mismatch type
			for (int n2 = 0; n2 < 4; ++n2) {
			    
			    double currentRatio3p = (totalObsCurrent3p == 0) ? static_cast<double>(0) : returnRatioFS((*dimer3pToUse)[l][4 * n1 + n2], totalObsCurrent3p, errorToRemove, failsafe);
			    double tmpRatio3p = (totalObsTmp3p == 0) ? static_cast<double>(0) : returnRatioFS((*dimer3pTmpToUse)[l][4 * n1 + n2], totalObsTmp3p, errorToRemove, failsafe);
			    
			    // Store the absolute difference in ratios
			    difference3p[l][4 * n1 + n2] = fabs(currentRatio3p - tmpRatio3p);
			}
		    }
		}
		
		int seq_len = b->core.l_qseq;
		int max_l   = std::min(seq_len, MAXLENGTH);
		
		// 1) collect the l‐indices: first 5 and last 5 (but not exceeding seq_len)
		std::vector<int> l_positions;
		for (int j = 0; j < 5 && j < seq_len; ++j)
		    l_positions.push_back(j);
		for (int j = std::max(5, seq_len - 5); j < seq_len; ++j)
		    l_positions.push_back(j);
		
		// Assume convergence across ALL l_positions and i==7,8
		bool allConverged = true;
		bool allConvergedStopEarly = true;
		
		for (int l : l_positions) {
		    for (int i : {7, 8}) {
			// double thresh = critRange ? critThresh : precisionThresh;
			
			// if *any* single check fails, we’re not fully converged
			if ( critRange ) { // is entered only in meta mode
			    if (!(difference5p[l][i] < critThresh && difference3p[l][i] < critThresh))
				{
				    // std::cerr << "critRange\t" << critRange << std::endl;
				    // std::cerr << "isConvergedLow\t" << isConvergedLow << std::endl;
				    // std::cerr << "processedAlns\t" << processedAlns << std::endl;
				    // std::cerr << "difference5p[l][i]\t" << difference5p[l][i] << std::endl;
				    // std::cerr << "difference3p[l][i]\t" << difference3p[l][i] << std::endl;
				    // std::cerr << "l\t" << l << std::endl;
				    // std::cerr << "i\t" << i << std::endl;
				    // std::cerr << "refNameStr\t" << refNameStr << std::endl;
				    // std::cerr << "\n" << std::endl;
				    allConverged = false;
				    break;
				}
			}
			else {
			    if (!(difference5p[l][i] < precisionThresh && difference3p[l][i] < precisionThresh))
				{
				    allConvergedStopEarly = false;
				    break;
				}
			}
		    }
		    if (!allConverged) break;  // no need to keep scanning
		    if (!allConvergedStopEarly) break;
		}
		
		// now set your flag if—and only if—all passed
		if (critRange && allConverged) { // can only be true in meta mode
		    // std::cerr << "allConverged\t" << allConverged << std::endl;
		    isConvergedLow = true;
		    isConvergedOnce = true;
		    // break;
		}
		
		if (!critRange && allConvergedStopEarly) {
		    stopEarly = true;
		    // std::cerr << "stopEarly\t" << stopEarly << std::endl;
		    break;
		}
		
		typesOfDimer5pTmp = typesOfDimer5p;
		typesOfDimer3pTmp = typesOfDimer3p;
		typesOfDimer5p_cpgTmp = typesOfDimer5p_cpg;
		typesOfDimer3p_cpgTmp = typesOfDimer3p_cpg;
		typesOfDimer5p_noncpgTmp = typesOfDimer5p_noncpg;
		typesOfDimer3p_noncpgTmp = typesOfDimer3p_noncpg;
		typesOfDimer5pDoubleTmp = typesOfDimer5pDouble;
		typesOfDimer3pDoubleTmp = typesOfDimer3pDouble;
		typesOfDimer5pSingleTmp = typesOfDimer5pSingle;
		typesOfDimer3pSingleTmp = typesOfDimer3pSingle;
	    }
	    
	    processedAlns++;
	}
	
	if ( metaMode ){
	    if ( isConvergedLow ){
		generateDamageProfile(outDir,
				      file5pDefault,
				      file3pDefault,
				      bamfiletopen, refNameStr, lengthMaxToPrint, dpFormat, hFormat, 
				      allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
				      genomeFileB, cpg, errorToRemove, failsafe, phred, 
				      typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
				      typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
				      typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p, 
				      typesOfDimer3p_cpg, typesOfDimer3p_noncpg, mapped);
		if(compFlag){
		    generateBaseCompositionProfile(outDir, bamfiletopen, refNameStr, lengthMaxToPrint, aroundFlank, genomeFileB, typesOfDimer5p, typesOfDimer3p, baseComp5pFlank, baseComp3pFlank, mapped);
		}
	    }
	    else{
		generateDamageProfile(outDirUnsafe,
				      file5pDefault,
				      file3pDefault,
				      bamfiletopen, refNameStr, lengthMaxToPrint, dpFormat, hFormat, 
				      allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
				      genomeFileB, cpg, errorToRemove, failsafe, phred, 
				      typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
				      typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
				      typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p,
				      typesOfDimer3p_cpg, typesOfDimer3p_noncpg, mapped);
		if(compFlag){
		    generateBaseCompositionProfile(outDirUnsafe, bamfiletopen, refNameStr, lengthMaxToPrint, aroundFlank, genomeFileB, typesOfDimer5p, typesOfDimer3p, baseComp5pFlank, baseComp3pFlank, mapped);
		}
	    }
	}

        // Clean up
        hts_itr_destroy(iter);
        bam_destroy1(b);
	
	//std:cerr << "processedAlns\t" << processedAlns << std::endl;
    }
    
    if ( classicMode ){
	// if ( isConvergedOnce || totalMapped >= numAlnsSafeRange ){
	generateDamageProfile(outDir,
			      file5pDefault,
			      file3pDefault,
			      bamfiletopen, "classic", lengthMaxToPrint, dpFormat, hFormat, 
			      allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
			      genomeFileB, cpg, errorToRemove, failsafe, phred, 
			      typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
			      typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
			      typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p,
			      typesOfDimer3p_cpg, typesOfDimer3p_noncpg, totalMapped);
	if(compFlag){
	    generateBaseCompositionProfile(outDir, bamfiletopen, "classic", lengthMaxToPrint, aroundFlank, genomeFileB, typesOfDimer5p, typesOfDimer3p, baseComp5pFlank, baseComp3pFlank, totalMapped);
	}
	// }
	// else{
	// 	generateDamageProfile(outDirUnsafe, bamfiletopen, "classic", lengthMaxToPrint, dpFormat, hFormat, 
	// 						allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
	// 						genomeFileB, cpg, errorToRemove, failsafe, phred, 
	// 						typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
	// 						typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
	// 						typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p, 
	// 						typesOfDimer3p_cpg, typesOfDimer3p_noncpg, totalMapped);
	// 	}
    }

    if(isizeB){
	std::ofstream isizeFP(isizeFile.c_str());
	if(!isizeFP){
	    std::cerr << "Could not open insert size output file " << isizeFile << std::endl;
	    return 1;
	}
	for(const auto & lengthAndCount : insertSizeCounts){
	    isizeFP << lengthAndCount.second << "\t" << lengthAndCount.first << "\n";
	}
	isizeFP.close();
	const pair<string, const map<int32_t,uint64_t>*> perClass[2] = { make_pair(isizeFile+(isizeAllPairs ? ".paired" : ".properly_paired"), &insertSizeCountsPaired),
									 make_pair(isizeFile+".merged", &insertSizeCountsMerged) };
	for(const auto & fileAndCounts : perClass){
	    std::ofstream classFP(fileAndCounts.first.c_str());
	    if(!classFP){
		std::cerr << "Could not open insert size output file " << fileAndCounts.first << std::endl;
		return 1;
	    }
	    for(const auto & lengthAndCount : *fileAndCounts.second){
		classFP << lengthAndCount.second << "\t" << lengthAndCount.first << "\n";
	    }
	}
    }

    // Clean up
    hts_idx_destroy(idx);
    bam_hdr_destroy(h);
    sam_close(fp);



    ///cerr<<"TEST2"<<endl;
    if (bed) bed_destroy(bed);
    //cerr<<"TEST3"<<endl;
    return 0;
}

