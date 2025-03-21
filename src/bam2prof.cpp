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
    bool quiet=false;
	bool classicMode=false;
	double precisionThresh=0.0;
	//string refId;

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
			"\t\t"+"-classic\t\tOne Profile per bam file    (Default: "+booleanAsString( classicMode )+" ) \n"+
			"\t\t"+"-precision\t\tSet minimum precision for substitution frequency computation (Default: All alignments [= 0.0]; Speed up by setting precision to either 0.001, 0.0001, 0.00001, ... ) \n"+
			"\t\t"+"-minAligned\t\tNumber of aligned sequences after which substitution patterns are checked for converging (Default: "+stringify( numAlns )+")\n"+

			"\t\t"+"-ref-id\t\t\tSpecify reference ID; if multiple references: Provide comma seperated list (no spaces!)    ( Default: Not Set ) \n"+

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


    if(genomeFileB){
	genome=new IndexedGenome(genomeFile.c_str());
	cerr<<genomeFile<<" mapped into memory"<<endl;
    }


	// string bamfilelist = string( argv[ argc-1 ] );
    // stringstream ss(bamfilelist);
    // string bamitem;
	// while (std::getline(ss, bamitem, ',')) {
	// 	bamFiles.push_back(bamitem);  // Add each ref ID to the vector
	// }

	
	string bamfiletopen = string( argv[ argc-1 ] );

	bam_hdr_t *h;
	samFile  *fp;
	hts_idx_t *idx;

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

	// Load the index for the BAM file
	idx = sam_index_load(fp, bamfiletopen.c_str());
	if(idx == NULL){
		std::cerr << "Could not load index for " << bamfiletopen << std::endl;
		return 1;
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

	// Iterate over each reference 
	for (int i = 0; i < h->n_targets; i++) {
		unsigned int processedAlns = 1;

		//std::cerr << "iterating targets" << std::endl;
		const char* refName = h->target_name[i];
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
		}

		hts_itr_t *iter = sam_itr_queryi(idx, i, 0, h->target_len[i]);
		bam1_t *b = bam_init1();

		// Check if iter is null
		if (iter == NULL) {
			std::cerr << "Could not create iterator for target " << h->target_name[i] << std::endl;
			continue;
		}

		string refFromFasta_;
		string refFromFasta;
		
		pair< kstring_t *, vector<int> >  reconstructedReference;
		reconstructedReference.first =(kstring_t *) calloc(sizeof(kstring_t),1);
		reconstructedReference.first->s =0 ;
		reconstructedReference.first->l =    reconstructedReference.first->m =0;


        // long unsigned int mapped = 0;
		// for (int map = 0; map < sam_hdr_nref(h); ++map) {
        //     uint64_t u, v;
        //     hts_idx_get_stat(idx, map, &u, &v);
		// 	mapped += u;
        // }
		// std::cerr << "bamfiletopen\t" << bamfiletopen << std::endl;
		// std::cerr << "refNameStr\t" << refNameStr << std::endl;
		// std::cerr << "mapped\t" << mapped << std::endl;
		// std::cerr << "unmapped\t" << unmapped << std::endl;
		// Get the number of mapped reads for the current reference 'i'
		uint64_t mapped = 0, unmapped = 0;
		hts_idx_get_stat(idx, i, &mapped, &unmapped);
		totalMapped += mapped;
		// std::cerr << "bamfiletopen\t" << bamfiletopen << std::endl;
		// std::cerr << "refNameStr\t" << refNameStr << std::endl;
		// std::cerr << "mapped\t" << mapped << std::endl;
		// std::cerr << "unmapped\t" << unmapped << std::endl;


		// Iterate over the BAM records
		while (sam_itr_next(fp, iter, b) >= 0) {
			
			if(bam_is_unmapped(b)){
				if(!quiet)
					cerr<<"skipping "<<bam_get_qname(b)<<" unmapped"<<endl;
				continue;
			}
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

			// update matrix
			countSubsPerRef(genomeFileB, genome, b, reconstructedReference, minQualBase, refFromFasta, refFromFasta_, h, bed, mask, ispaired, isfirstpair, typesOfDimer5p, typesOfDimer3p, typesOfDimer5p_cpg, typesOfDimer3p_cpg, typesOfDimer5p_noncpg, typesOfDimer3p_noncpg, typesOfDimer5pDouble, typesOfDimer3pDouble, typesOfDimer5pSingle, typesOfDimer3pSingle);


			// bool getDamProf = true;
			// if ( mapped <= 10000 ){
			// 	numAlns = 100;
			// 	precisionThresh = 0.01;
			// 	getDamProf = false;
			// }






			// Early stop mechanism: Check every numAlns alignments
			if (processedAlns % numAlns == 0) {

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


				bool stopEarly = false;  // Flag to indicate if all changes are small

				for (int l = 0; l < MAXLENGTH; ++l) {
					for (int i = 0; i < 16; ++i) {
						if ( precisionThresh != 0 ){
							if (( difference5p[l][i] != 0 && difference5p[l][i] < precisionThresh ) && ( difference3p[l][i] != 0 && difference3p[l][i] < precisionThresh)) {
								//std::cerr << "difference5p[l][i] " << difference5p[l][i] << std::endl;
								//std::cerr << "difference3p[l][i] " << difference3p[l][i] << std::endl;
								stopEarly = true;  // If any change exceeds 0.001, set flag to false
								break;  // Break out early if we find a significant change
							}
						}
					}
					if (stopEarly) break;  // Exit the outer loop if we already know there's a significant change
				}

				// If no significant changes were detected across all positions and types, we can stop further iterations
				if (stopEarly) {
					//std::cout << "No significant change detected. Stopping iterations." << std::endl;
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

		if ( !classicMode ){
			generateDamageProfile(outDir, bamfiletopen, refNameStr, lengthMaxToPrint, dpFormat, hFormat, 
								allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
								genomeFileB, cpg, errorToRemove, failsafe, phred, 
								typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
								typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
								typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p, 
								typesOfDimer3p_cpg, typesOfDimer3p_noncpg, mapped);
		}

        // Clean up
        hts_itr_destroy(iter);
        bam_destroy1(b);
    }


	if ( classicMode ){
		generateDamageProfile(outDir, bamfiletopen, "classic", lengthMaxToPrint, dpFormat, hFormat, 
							allStr, singAnddoubleStr, doubleStr, singleStr, endo, 
							genomeFileB, cpg, errorToRemove, failsafe, phred, 
							typesOfDimer5pSingle, typesOfDimer5pDouble, typesOfDimer5p, 
							typesOfDimer5p_cpg, typesOfDimer5p_noncpg, 
							typesOfDimer3pSingle, typesOfDimer3pDouble, typesOfDimer3p, 
							typesOfDimer3p_cpg, typesOfDimer3p_noncpg, totalMapped);
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

