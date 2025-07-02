version 1.0

workflow hilow {
    meta {
        title: "StJude HiC Looping Organization Workflow"
        summary: "HiC parallel data processing, Peak and loop Calling Pipeline"
        description: "An all-in-one pipeline which takes FASTQs and perform mapping by HiC-Pro, filtering, merging, normalization, converting to .hic files, peaks and loops analysis."
        version: "1.1"
        details: {
            citation: "",
            contactEmail: [
                "jie.lu@stjude.org",
                "brian.abraham@stjude.org",
                "modupeore.adetunji@stjude.org",
            ],
            contactOrg: "St Jude Children's Research Hospital",
            contactUrl: "",
            upstreamLicenses: "MIT",
            upstreamUrl: "",
            whatsNew: "",
        }
        outputs: {
            hichipjsonfile: "JSON file for HiChIP results",
            hicjsonfile: "JSON file for HiC results",
            hicpro_tar: "Tar file containing HiC-Pro results",
            hicfile: "HiC file in .hic format",
            qc_outputzip: "QC output files in zip format",
            active_qcfiles: "Active QC files from HiChIP analysis",
            peaksbw: "BigWig file of HiChIP peaks",
            peaksbdg: "BedGraph file of HiChIP peaks",
            peaksbed: "BED file of HiChIP peaks",
            peakslog: "Log file of HiChIP peak calling",
            peaksreport: "Report file of HiChIP peak calling",
            fithichip_zip: "Zip file containing FitHiChIP results",
            active_fithichip: "Active FitHiChIP results zip file",
            anchor_bw: "BigWig file of loop anchors",
            pp_bw: "ProteinPaint BigWig file",
            pp_hic: "ProteinPaint HiC file",
            pp_peaks: "ProteinPaint peaks BED file",
            pp_tbi: "ProteinPaint TBI index for peaks BED file",
            bedloop_1: "BED file for loop type 1",
            bedloop_2: "BED file for loop type 2",
            bedloop_3: "BED file for loop type 3",
            tbiloop_1: "TBI index for loop type 1 BED file",
            tbiloop_2: "TBI index for loop type 2 BED file",
            tbiloop_3: "TBI index for loop type 3 BED file",
            fastqchtml_1: "Array of FastQC HTML files for R1 reads",
            fastqchtml_2: "Array of FastQC HTML files for R2 reads",
            fastqczip_1: "Array of FastQC zip files for R1 reads",
            fastqczip_2: "Array of FastQC zip files for R2 reads",
            qcreports: "QC report files from HiC analysis",
            active_qcreports: "Active QC report files from HiChIP analysis",
            active_promoters: "Active promoters BED file from GTF annotation processing",
            active_generegion: "Active gene regions BED file from GTF annotation processing",
        }
    }

    parameter_meta {
        fastqfiles_R1: {
            description: "One or more sample R1 FASTQs",
            help: "Upload zipped R1 FASTQ files.",
            patterns: [
                "*_R1.fq.gz",
                "*_R1.fastq.gz",
            ],
        }
        fastqfiles_R2: {
            description: "One or more sample R2 FASTQs",
            help: "Upload zipped R2 FASTQ files.",
            patterns: [
                "*_R2.fq.gz",
                "*_R2.fastq.gz",
            ],
        }
        genomefragment: {
            description: "restriction fragment file in BED format",
            help: "Input genome fragments digested from restriction enzyme(s).",
            patterns: [
                "*.bed",
                "*.bed.gz",
            ],
        }
        genomename: {
            description: "Reference Genome Name",
            help: "Reference genome name. Default: hg19. Examples: hg19, mm10, mm9",
            example: "hg19",
        }
        pp_directory: {
            description: "ProteinPaint directory",
            help: "preferred ProteinPaint directory name, the directory contains the visualization files and required json. Directory should be copied to the ProteinPaint '*/tp' directory.",
        }
        sampleid: {
            description: "Sample ID",
            help: "Sample ID for the output files. If not defined, the sample ID will be derived from the FASTQ file names.",
            example: "Sample1",
        }
        reference: {
            description: "Reference FASTA file",
            patterns: [
                "*.fa",
                "*.fasta",
                "*.fa.gz",
                "*.fasta.gz",
            ],
        }
        blacklist: {
            description: "Blacklist file in BED format",
            help: "If defined, blacklist regions listed are excluded after reference alignment.",
            patterns: [
                "*.bed",
                "*.bed.gz",
            ],
        }
        chromsizes: {
            description: "Chromosome sizes file in tab delimited format",
            help: "If defined, genome chromosome size file.",
        }
        loopsAnchor: {
            description: "Optional file containing pre-defined set of loop anchors, only the first 3 cols are relevant (Chr,Start, End). 1D peak calling step will skipped if provided",
            patterns: [
                "*.bed"
            ],
        }
        promoters: {
            description: "Optional file containing pre-defined set of promoters(4kb regions around TSS), in the format of chr,start,end,TranscriptID|GeneID. LoopAnchor will be generated by merging 1D HiChIP peak with top 2/3 promoters by HiChIP coverage."
        }
        geneAnnotation: {
            description: "Input gene annotation file from RefSeq or GENCODE (.gtf).",
            patterns: [
                "*.gtf",
                "*.gtf.gz",
            ],
        }
        bowtie2index: {
            description: "bowtie v2 index files (*.bt2)",
            help: "If not defined, bowtie v2 index files are generated, will take a longer compute time.",
            patterns: [
                "*.bt2"
            ],
        }
        ligationsite: {
            description: "Restriction enzyme ligation site used by HiC-Pro"
        }
        ashichip: {
            description: "Boolean to indicate if HiChIP peak calling should be performed",
            help: "If set to True, HiChIP peak calling will be performed, otherwise only HiC analysis will be performed.",
            default: true,
        }
    }

    input {
        # Required Inputs
        Array[File] fastqfiles_R1
        Array[File] fastqfiles_R2
        File? genomefragment
        String genomename = "hg19"
        # Optional Inputs
        String pp_directory = "browser_files"
        String? sampleid
        File? reference
        File? blacklist
        File? chromsizes
        File? loopsAnchor
        File? promoters
        File? geneAnnotation
        Array[File]? bowtie2index
        # HiC and HiChIP mapping using HiC-Pro
        String? ligationsite
        # Limit HiC analysis
        Boolean ashichip = true
    }

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 1 ----------- ###
    ### ------ Pre-process Analysis Files ------ ###
    ### ---------------------------------------- ###

    # Required Parameters but not for user inputs
    String hicpro_output = "HiCPro_out"
    Float FDR = 0.4
    Int LowerThreshold = 5000
    Int UpperThreshold = 4000000
    Int IntType = 3
    Float loopThreshold_1 = 0.01
    Float loopThreshold_2 = 0.05
    Float loopThreshold_3 = 0.1

    # IntType is 1 or 3. Default 3
    if (IntType != 1 && IntType != 3) {
        call raise_exception as error_wrong_IntType { input: msg = "Type of interaction reported by FitHiChIP. Choices: 1 (peak to peak) or 3 (peak to all). Default = 3"}
    }

    if (ashichip && (!defined(genomefragment) || !defined(ligationsite))){
        call raise_exception as error_missing_data { input: msg = "Genome Fragment file and/or Ligation sites are not specified for HiChIP peak calling, set 'ashichip' to False" }
    }

    # Generating INDEX files
    # 1. Bowtie INDEX files if not provided
    if ( !defined(bowtie2index) ||
    length(select_first([bowtie2index, []])) != 6 ) {
        # create bowtie index when not provided
        call bowtie2_index { input:
            genomename = genomename,
            reference = select_first([reference, ""]),
        }
    }

    Array[File] actual_bowtie2_index = select_first([bowtie2_index.bowtie_indexes, bowtie2index])

    # 2. Create FASTA index and chrom sizes files
    if (! defined(chromsizes) ){
        call chrfaidx { input:
            reference = select_first([reference, ""])
        }
    }
    File actual_chromsizes = select_first([chrfaidx.chromsizes, chromsizes])

    # 3. Create promoters file from gtf as chr,start,end,TranscriptID|GeneID
    if ( defined(geneAnnotation) ){
        call promgeneregion { input:
            annotation = select_first([geneAnnotation, ""])
        }
    }

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 2 ----------- ###
    ### ------ Perform HiC-Pro Alignments ------ ###
    ### ---------------------------------------- ###

    Int number_fastqs = length(fastqfiles_R1)
    scatter (fastqfile_R1 in fastqfiles_R1) {
        call fastqsnumber as extrapolate_R1 { input:
            fastqfile = fastqfile_R1,
            lengthfastqs = number_fastqs,
        }
        call splitfastqs as validsplit_R1 { input:
            nreads = extrapolate_R1.batch_size,
            fastqfile = fastqfile_R1,
        }
        call fastqc as fastqc_R1 { input: fastq = fastqfile_R1 }
    }

    scatter (fastqfile_R2 in fastqfiles_R2) {
        call fastqsnumber as extrapolate_R2 { input:
            fastqfile = fastqfile_R2,
            lengthfastqs = number_fastqs,
        }
        call splitfastqs as validsplit_R2 { input:
            nreads = extrapolate_R2.batch_size,
            fastqfile = fastqfile_R2,
        }
        call fastqc as fastqc_R2 { input: fastq = fastqfile_R2 }
    }

    Array[Pair[File, File]] all_fastqfiles = zip(flatten(validsplit_R1.outputfiles), flatten(validsplit_R2.outputfiles))

    Array[Int] files = range(length(all_fastqfiles))
    scatter (i in files) {
        call hicpro_align as step1hicpro { input:
            hicpro_out = hicpro_output,
            fastqfile_R1 = all_fastqfiles[i].left,
            fastqfile_R2 = all_fastqfiles[i].right,
            ligationsite = ligationsite,
            genomefragment = genomefragment,
            genomename = genomename,
            chromsizes = actual_chromsizes,
            bowtie2index = actual_bowtie2_index,
            iteration = i + 1,
        }
    }
    call hicpro_merge as step2hicpro { input:
        hicpro_out = hicpro_output,
        fastqfiles_R1 = flatten(validsplit_R1.outputfiles),
        fastqfiles_R2 = flatten(validsplit_R2.outputfiles),
        hicpro_align = step1hicpro.hicpro_align,
        ligationsite = ligationsite,
        genomefragment = genomefragment,
        genomename = genomename,
        chromsizes = actual_chromsizes,
        bowtie2index = actual_bowtie2_index,
    }

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 3 ----------- ###
    ### ------- Filter Blacklist Regions ------- ###
    ### ---------------------------------------- ###

    if ( defined(blacklist) ) {
        File blacklist_file = select_first([blacklist, ""])

        call filterblklist { input:
            count_fastqs = step2hicpro.count_fastqs,
            hicpro_out = hicpro_output,
            hicpro_result = step2hicpro.hicpro_merge,
            blacklist = blacklist_file,
            chromsizes = actual_chromsizes,
        }
    }

    File final_hicpro = select_first([filterblklist.hicpro_filtered, step2hicpro.hicpro_merge])

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 4 ----------- ###
    ### ----------- Create .hic file ----------- ###
    ### ---------------------------------------- ###

    call compressfiles { input:
        inputfile = final_hicpro
    }

    call converthic { input:
        count_fastqs = step2hicpro.count_fastqs,
        hicpro_out = hicpro_output,
        pp_directory = pp_directory,
        hicpro_result = final_hicpro,
        chromsizes = actual_chromsizes,
        genomename = genomename,
        sampleid = sampleid,
    }

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 5 ----------- ###
    ### ------------- Peak Calling ------------- ###
    ### ---------------------------------------- ###

    if (ashichip) {
        if ( !defined(loopsAnchor) ){
            call oneDpeaks { input:
                count_fastqs = step2hicpro.count_fastqs,
                FDR = FDR,
                hicpro_out = hicpro_output,
                pp_directory = pp_directory,
                hicpro_result = final_hicpro,
                genomefragment = genomefragment,
                chromsizes = actual_chromsizes,
                genomename = genomename,
                sampleid = sampleid,
            }
        }

        if ( defined(promoters) || defined(geneAnnotation) ) {
            call activeRegionsMerge { input:
                    hicpro_out = hicpro_output,
                    loopBed = select_first([loopsAnchor, oneDpeaks.oneDpeaksbed, ""]),
                    promoters = select_first([promoters, promgeneregion.promoters, ""]),
                    genomename = genomename,
                    hicpro_result = final_hicpro,
            }
        }

        ### ---------------------------------------- ###
        ### ------------ S E C T I O N 6 ----------- ###
        ### ------------- Loop Calling ------------- ###
        ### ---------------------------------------- ###

        call twoDloops { input:
            count_fastqs = step2hicpro.count_fastqs,
            FDR = FDR,
            hicpro_out = hicpro_output,
            pp_directory = pp_directory,
            loopBed = select_first([loopsAnchor, oneDpeaks.oneDpeaksbed, ""]),
            loopsAnchor = loopsAnchor,
            hicpro_result = final_hicpro,
            chromsizes = actual_chromsizes,
            LowDistThr = LowerThreshold,
            UppDistThr = UpperThreshold,
            IntType = IntType,
            loopThreshold_1 = loopThreshold_1,
            loopThreshold_2 = loopThreshold_2,
            loopThreshold_3 = loopThreshold_3,
            genomename = genomename,
            sampleid = sampleid,
        }

        if ( defined(promoters) || defined(geneAnnotation) ){
            call twoDloops as active_twoDloops { input:
                count_fastqs = step2hicpro.count_fastqs,
                FDR = FDR,
                hicpro_out = hicpro_output,
                pp_directory = pp_directory,
                loopBed = select_first([activeRegionsMerge.combined_promoters, ""]),
                promoters = select_first([promoters, promgeneregion.promoters]),
                hicpro_result = final_hicpro,
                chromsizes = actual_chromsizes,
                LowDistThr = LowerThreshold,
                UppDistThr = UpperThreshold,
                IntType = IntType,
                loopThreshold_1 = loopThreshold_1,
                loopThreshold_2 = loopThreshold_2,
                loopThreshold_3 = loopThreshold_3,
                genomename = genomename,
                sampleid = sampleid,
            }

            call qcreport as active_qcreport { input:
                sampleid = sampleid,
                genomename = genomename,
                activeregion = true,
                qczip = select_first([active_twoDloops.qc_output, ""]),
            }
        }

        File final_pp_peaks = select_first([twoDloops.pp_peaks, twoDloops.pp_peaks2, ""])
        File final_pp_tbi = select_first([twoDloops.pp_tbi,twoDloops.pp_tbi2, ""])

        call createjson as hichipjson { input:
            pp_bw = select_first([oneDpeaks.pp_bw,twoDloops.anchor_bw]),
            pp_hic = converthic.pp_hic,
            pp_peaks = final_pp_peaks,
            bedloop_1 = twoDloops.bedloop_1,
            bedloop_2 = twoDloops.bedloop_2,
            bedloop_3 = twoDloops.bedloop_3,
            pp_directory = pp_directory,
            IntType = IntType,
            genomename = genomename,
            sampleid = sampleid,
        }

        call qcreport { input:
            sampleid = sampleid,
            genomename = genomename,
            qczip = select_first([twoDloops.qc_output, ""]),
        }
    }  # end if hichip

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 7 ----------- ###
    ### ----------- Create JSON File ----------- ###
    ### ---------------------------------------- ###

    if (! ashichip) {
        call createjson as hicjson { input:
            pp_hic = converthic.pp_hic,
            pp_directory = pp_directory,
            IntType = IntType,
            genomename = genomename,
            sampleid = sampleid,
        }

        call qcreport as hic_report { input:
            sampleid = sampleid,
            genomename = genomename,
            qczip = step2hicpro.qc_output,
        }
    }

    ### ---------------------------------------- ###
    ### ------------ S E C T I O N 8 ----------- ###
    ### --------------- The End ---------------- ###
    ### ---------------------------------------- ###

    output {
        File? hichipjsonfile = hichipjson.json
        File? hicjsonfile = hicjson.json
        File? hicpro_tar = compressfiles.zipfile
        File? hicfile = converthic.hic_file
        File? qc_outputzip = select_first([twoDloops.qc_output, step2hicpro.qc_output])
        File? active_qcfiles = active_twoDloops.qc_output
        File? peaksbw = oneDpeaks.oneDpeaksbw
        File? peaksbdg = oneDpeaks.oneDpeaksbdg
        File? peaksbed = oneDpeaks.oneDpeaksbed
        File? peakslog = oneDpeaks.oneDpeakslog
        File? peaksreport = oneDpeaks.oneDpeaksreport
        File? fithichip_zip = twoDloops.twoDloops_zip
        File? active_fithichip = active_twoDloops.twoDloops_zip
        File? anchor_bw = twoDloops.anchor_bw
        File? pp_bw = oneDpeaks.pp_bw
        File? pp_hic = converthic.pp_hic
        File? pp_peaks = final_pp_peaks
        File? pp_tbi = final_pp_tbi
        File? bedloop_1 = twoDloops.bedloop_1
        File? bedloop_2 = twoDloops.bedloop_2
        File? bedloop_3 = twoDloops.bedloop_3
        File? tbiloop_1 = twoDloops.tbiloop_1
        File? tbiloop_2 = twoDloops.tbiloop_2
        File? tbiloop_3 = twoDloops.tbiloop_3
        Array[File?] fastqchtml_1 = fastqc_R1.htmlfile
        Array[File?] fastqchtml_2 = fastqc_R2.htmlfile
        Array[File?] fastqczip_1 = fastqc_R1.zipfile
        Array[File?] fastqczip_2 = fastqc_R2.zipfile
        File? qcreports = select_first([hic_report.qcreport, qcreport.qcreport])
        File? active_qcreports = active_qcreport.qcreport
        File? active_promoters = promgeneregion.promoters
        File? active_generegion = promgeneregion.generegion
    }
}

### ---------------------------------------- ###
### ------------ TASKS UTILIZED ------------ ###
### ---------------------------------------- ###

task promgeneregion {
    meta {
        description: "Generate promoter and gene regions from GTF annotation file"
        outputs: {
            promoters: "File containing promoter regions in BED format",
            generegion: "File containing gene regions in BED format",
        }
    }

    parameter_meta {
        annotation: {
            description: "Input gene annotation file from RefSeq or GENCODE (.gtf).",
            patterns: [
                "*.gtf",
                "*.gtf.gz",
            ],
        }
        output_promoters: "File name for output file containing promoter regions"
        output_generegion: "File name for output file containing  gene regions"
        outdir: "Directory to store output files during task execution"
    }

    input {
        File annotation
        String output_promoters = "PromoterRegions.bed"
        String output_generegion = "GeneRegions.bed"
        String outdir = "RegionFiles"
    }

    String base = basename(annotation, ".gz")

    command <<<
        set -euo pipefail
        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}

        sed -i 's/ /\t/g' ~{base}
        
        #get promoterregion
        awk -F\\t '{ if($3 == "transcript" && $1 !~ "chrM")
            if($7 =="+")
                print $1 "++" $4-2000 "++" $4+2000 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3)
            else if($7=="-")
                print $1 "++" $5-2000 "++" $5+2000 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3) }' \
            ~{base} \
            | grep -v "unassigned" \
            | sort \
            | uniq \
            | sed s/\+\+/\\t/g >| aaa
        
        
        mkdir -p ~{outdir}
        sort -k1,1 -k2,2n aaa > ~{outdir}/~{output_promoters}
            && rm -rf aaa

        #get whole gene
        awk -F\\t '{ if($3 == "transcript" && $1 !~ "chrM") 
            print $1 "++" $4 "++" $5 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3)}' \
            ~{base} \
            | grep -v "unassigned" \
            | sort \
            | uniq \
            | sed s/\+\+/\\t/g >| aaa

        sort -k1,1 -k2,2n aaa > ~{outdir}/~{output_generegion}
            && rm -rf aaa
    >>>

    output {
        File promoters = "~{outdir}/~{output_promoters}"
        File generegion = "~{outdir}/~{output_generegion}"
    }

    runtime {
        memory: "5 GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task fastqsnumber {
    meta {
        description: "Calculate number of reads in FASTQ files in order to yield a split into 30 jobs or less"
        outputs: {
            batch_size: "Number of reads per split FASTQ file",
            total_reads: "Total number of reads in the input FASTQ file",
        }
    }

    parameter_meta {
        fastqfile: {
            description: "Input FASTQ file",
            patterns: [
                "*.fastq.gz",
                "*.fq.gz",
            ],
        }
        lengthfastqs: {
            description: "Number of FASTQ files provided by the user for read 1 or read 2",
        }
        jobs: {
            description: "The maximum number of jobs to run",
        }
    }

    # Get number of fastq reads to split into 30 jobs or less
    input {
        File fastqfile
        Int lengthfastqs
        Int jobs = 30
    }

    Int minimum_batch_size = jobs * 1000000

    command <<<
        # Compute number of threads to dedicate to this FASTQ file
        maxthreads=$(printf %.0f $(echo "~{jobs} / ~{lengthfastqs}" | bc)); echo $maxthreads
        # Get the count of lines in the FASTQ file
        count=$(zcat ~{fastqfile} | wc -l); echo $count
        # Compute the batches in millions of reads
        reads=$(printf %.0f $(echo "$count/(4*$maxthreads)" | bc)); echo $reads
        if [ $reads -le ~{minimum_batch_size} ]; then
            reads=~{minimum_batch_size};
        fi

        echo "$count" > total_reads.txt
        echo $reads > batch_size.txt
    >>>

    output {
        Int batch_size = read_int("batch_size.txt")
        File total_reads = "total_reads.txt"
    }

    runtime {
        memory: "5 GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task splitfastqs {
    meta {
        description: "Split FASTQ files into smaller chunks for HiC-Pro processing"
        outputs: {
            outputfiles: "Array of split FASTQ files"
        }
    }

    parameter_meta {
        fastqfile: {
            description: "Input FASTQ file to be split",
            patterns: [
                "*.fastq.gz",
                "*.fq.gz",
            ],
        }
        nreads: {
            description: "Number of reads to per split",
        }
        output_location: {
            description: "Subdirectory to store split FASTQ files",
        }
        prefix: {
            description: "Prefix for the split FASTQ file(s). The extension `.fastq` will be added.",
        }
    }

    # Split fastq files into chunks for running HiC-Pro
    input {
        File fastqfile
        Int nreads = 30000000
        String output_location = "_split/fastq"
        String prefix = basename(fastqfile)
    }

    command <<<
        mkdir -p ~{output_location}
        cd ~{output_location}
        prefix=~{prefix}
        prefix=${prefix%.f*q*}

        # Making sure paired reads are R1 and R2.
        if [[ $prefix != *"_R1_0"* ]]; then
            if [[ $prefix != *"_R2_0"* ]]; then
                if [[ $prefix == *"_1"* ]]; then
                    prefix=${prefix/_1/_R1}
                elif [[ $prefix == *"_2"* ]]; then
                    prefix=${prefix/_2/_R2}
                fi
            fi
        fi

        nlines=$(echo "~{nreads} * 4" | bc)

        if [[ "~{fastqfile}" == *"gz" ]]; then
            zcat ~{fastqfile} | split -l $nlines -d - temp
        else
            split -l $nlines -d ~{fastqfile} temp
        fi

        files=$(ls temp* | sort)
        for each in $(ls temp* | sort); do mv $each ${each:4}_$prefix.fastq; done
    >>>

    output {
        Array[File] outputfiles = glob("~{output_location}/*.fastq")
    }

    runtime {
        memory: "5 GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task hicpro_align {
    meta {
        description: "Perform HiC-Pro step 1: mapping and filtering of FASTQ files"
        outputs: {
            hicpro_align: "Tar file containing HiC-Pro alignment results"
        }
    }

    parameter_meta {
        fastqfile_R1: "File containing read ones"
        fastqfile_R2: "File containing read twos"
        chromsizes: "Chromosome sizes file in tab delimited format"
        genomefragment: "Restriction fragment file in BED format"
        genomename: "Reference genome name"
        ligationsite: "Restriction enzyme ligation site used by HiC-Pro"
        bowtie2index: "Array of bowtie v2 index files (*.bt2)"
        hicpro_out: "Output directory for HiC-Pro results"
        iteration: "Iteration number for the current FASTQ file"
        ncpu: "Number of CPU cores to use for processing"
        memory_gb: "Memory in GB to allocate for the task"
    }

    input {
        File fastqfile_R1
        File fastqfile_R2
        File chromsizes
        File? genomefragment
        String genomename
        String? ligationsite
        Array[File] bowtie2index
        String hicpro_out = "HiCProOut"
        Int iteration = 0
        Int ncpu = 4
        Int memory_gb = 5
    }

    File actual_genomefragment = select_first([genomefragment, ""])

    command <<<
        set -euo pipefail
        # Update Config File
        mkdir temp
        pwd=$(pwd)

        config=~{basename(fastqfile_R1)}.tar
        cp /data/config-hicpro.txt config-hicpro.txt

        newfastqname=~{sub(basename(fastqfile_R1),".fastq|.fq|.fastq.gz|.fq.gz","")}
        newfastqname_edit="${newfastqname/_R1/}"

        # Allow for Genome Fragment to not be specified
        if [ -f ~{actual_genomefragment} ]; then
            fragment_filename=~{basename(actual_genomefragment, ".gz")}
            GFragment=$pwd/$fragment_filename
            gunzip -c ~{actual_genomefragment} > $GFragment \
                || ln -sf ~{actual_genomefragment} $GFragment
            genomefragment=$GFragment
        else
            genomefragment=""
        fi

        # Allow for Ligation Site to not be specified
        ligationsite=~{ligationsite}

        # Make sure Min_Cis_Dist is set when "Fragments" and "Ligation Sites" are not specified
        if [[ ${#ligationsite} -le 1 ]] && [[ ${#assignfragment} -le 1 ]] ; then
            sed -i "s/MIN_CIS_DIST\ \=/MIN_CIS_DIST\ \=\ 1000/" config-hicpro.txt
        fi

        sed -i "s#Xgenomefragment#${genomefragment}#;s#GATCGATC#${ligationsite}#" config-hicpro.txt

        index_path=~{bowtie2index[0]}; index_path=${index_path%/*}

        # Changes to the config file
        sed -i "s#Xbowtieindex#${index_path}#;s#Xchromsizes#~{chromsizes}#;s#Xgenome#~{genomename}#" config-hicpro.txt

        # Copy Files to fastq location
        mkdir -p _split/fastq
        cp ~{fastqfile_R1} ~{fastqfile_R2} _split/fastq
        echo `date`

        newcount=~{iteration}
        echo "this is sleep time $newcount"

        # sleep $newcount #$(( ( RANDOM % 10 )  + 1 ))
        # Generating SUB steps
        /HiC-Pro_3.1.0/bin/HiC-Pro -i _split -o ~{hicpro_out} -c config-hicpro.txt -p

        # Execute Alignment Step
        if [ -d ~{hicpro_out} ]; then
            cd ~{hicpro_out}
            bash HiCPro_step1_.sh

            if [ -f $pwd/~{hicpro_out}/hic_results/data/fastq/$newfastqname_edit"_"~{genomename}.bwt2pairs.validPairs ]; then
                # removing not neccessary files
                rm rawdata HiCPro_step1_.sh HiCPro_step2_.sh inputfiles_.txt config-hicpro.txt

                # Compress Alignment Results
                cd $pwd
                mkdir ~{basename(fastqfile_R1)}
                mv ~{hicpro_out} ~{basename(fastqfile_R1)}
                tar -cpf ~{basename(fastqfile_R1)}.tar ~{basename(fastqfile_R1)}

            else
                echo -e "hilow error: Can not find valid pairs. Failed to align FASTQs using HiCPro."
                exit 1
            fi

        else
            echo -e "hilow error: ~{hicpro_out} folder was not created"
            exit 1
        fi
    >>>

    output {
        File hicpro_align = "~{basename(fastqfile_R1)}.tar"
    }

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: ncpu
    }
}

task hicpro_merge {
    meta {
        description: "Perform HiC-Pro step 2: merging and normalization of HiC-Pro results"
        outputs: {
            hicpro_merge: "Merged and normalized HiC-Pro results in tar format",
            qc_output: "QC files in zip format",
            count_fastqs: "Number of split FASTQs processed",
        }
    }

    parameter_meta {
        fastqfiles_R1: " Array of read one FASTQ files"
        fastqfiles_R2: " Array of read two FASTQ files"
        hicpro_align: "Array of HiC-Pro alignment tar files"
        chromsizes: "Chromosome sizes file in tab delimited format"
        genomefragment: "Restriction fragment file in BED format"
        genomename: "Reference genome name"
        ligationsite: "Restriction enzyme ligation site used by HiC-Pro"
        bowtie2index: "Array of bowtie v2 index files (*.bt2)"
        hicpro_out: "Output directory name for HiC-Pro results"
    }

    # Perform HiC-Pro step 2: merging and normalization
    input {
        Array[File] fastqfiles_R1
        Array[File] fastqfiles_R2
        Array[File] hicpro_align
        File chromsizes
        File? genomefragment
        String genomename
        String? ligationsite
        Array[File] bowtie2index
        String hicpro_out = "HiCProOut"
    }

    File actual_genomefragment = select_first([genomefragment, ""])
    Int memory_gb = ceil(length(fastqfiles_R1) * 8 )  # memory is times 3 of number of split fastqs (which = <= 90)
    Int int_fastqs = length(fastqfiles_R1)  # number of split fastqs

    command <<<
        # Update Config File
        pwd=$(pwd)
        cp /data/config-hicpro.txt config-hicpro.txt

        # Allow for Genome Fragment to not be specified
        if [ -f ~{actual_genomefragment} ]; then
            fragment_filename=~{basename(actual_genomefragment, ".gz")}
            GFragment=$pwd/$fragment_filename
            gunzip -c ~{actual_genomefragment} > $GFragment \
                || ln -sf ~{actual_genomefragment} $GFragment
            genomefragment=$GFragment
        else
            genomefragment=""
        fi

        # Allow for Ligation Site to not be specified
        ligationsite=~{ligationsite}

        # Make sure Min_Cis_Dist is set when "Fragments" and "Ligation Sites" are not specified
        if [[ ${#ligationsite} -le 1 ]] && [[ ${#assignfragment} -le 1 ]] ; then
            sed -i "s/MIN_CIS_DIST\ \=/MIN_CIS_DIST\ \=\ 1000/" config-hicpro.txt
        fi
        sed -i "s#Xgenomefragment#${genomefragment}#;s#GATCGATC#${ligationsite}#" config-hicpro.txt
        
        index_path=~{bowtie2index[0]}; index_path=${index_path%/*}

        # Changes to the config file
        sed -i "s#Xbowtieindex#${index_path}#;s#Xchromsizes#~{chromsizes}#;s#Xgenome#~{genomename}#" config-hicpro.txt

        # Copy Files to fastq location
        mkdir -p _split/fastq
        cp -rf ~{sep=" " fastqfiles_R1} ~{sep=" " fastqfiles_R2} _split/fastq
        echo `date`

        # Generating SUB steps
        /HiC-Pro_3.1.0/bin/HiC-Pro -i _split -o ~{hicpro_out} -c config-hicpro.txt -p
        file_array=()
        mkdir -p AAA ; mv ~{hicpro_out} AAA
        echo `date`
        for TAR in ~{sep=" " hicpro_align}; do
            tar -xpf $TAR
            TAR=${TAR##*/}
            file_array+=("${TAR%.tar}/~{hicpro_out}")
        done

        # Organize HiCPro File
        cp -R ${file_array[@]} AAA/~{hicpro_out} ./
        cd ~{hicpro_out}
        bash HiCPro_step2_.sh

        # Compress Merge Results
        echo `date`
        cd $pwd
        # rm -rf ${file_array} AAA/~{hicpro_out}
        # rm -rf ~{basename(hicpro_out)}/bowtie_results
        # rm -rf ~{basename(hicpro_out)}/rawdata
        tar -cpf ~{basename(hicpro_out)}.tar ~{basename(hicpro_out)}

        # Export QCfiles
        mkdir -p QCfiles
        cp HiCPro_out/hic_results/pic/fastq/plot* QCfiles
        cp HiCPro_out/hic_results/stats/fastq/* QCfiles
        zip -9qr QCfiles.zip QCfiles
        # rm -rf $pwd/~{basename(hicpro_out)}

    >>>

    output {
        File hicpro_merge = "~{basename(hicpro_out)}.tar"
        File qc_output = "QCfiles.zip"
        Int count_fastqs = int_fastqs
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task oneDpeaks {
    meta {
        description: "Call 1D peaks from HiC-Pro results using HiChIP-Peaks"
        outputs: {
            oneDpeaksbw: "BigWig file of 1D peaks",
            oneDpeaksbdg: "BedGraph file of 1D peaks",
            oneDpeaksbed: "BED file of 1D peaks",
            oneDpeakslog: "Log file of 1D peak calling",
            oneDpeaksreport: "Report PDF of 1D peak calling",
            pp_bw: "ProteinPaint BigWig file for visualization",
        }
    }

    parameter_meta {
        count_fastqs: "Number of FASTQ files used in the HiC-Pro run"
        hicpro_result: "HiC-Pro result archive containing the processed data"
        genomefragment: "Genome fragment file for peak calling (optional)"
        FDR: "False Discovery Rate threshold for peak calling"
        chromsizes: "Chromosome sizes file for BigWig conversion"
        genomename: "Name of the genome used in the analysis"
        pp_directory: "Directory to store ProteinPaint output files"
        sampleid: "Sample identifier for naming output files"
        sampleid_m: "Modified sample identifier for output file naming"
        hicpro_out: "Directory where HiC-Pro results are stored"
        outdir: "Output directory for peak calling results"
        ncpu: "Number of CPU cores to use for peak calling"
    }

    # Call 1D peaks from  HiC-Pro result using HiChIP-Peaks
    input {
        Int count_fastqs
        File hicpro_result
        File? genomefragment
        Float FDR=0.4
        File chromsizes
        String genomename
        String pp_directory
        String? sampleid
        String sampleid_m = "~{sampleid + '.' + genomename + '_'}"
        String hicpro_out = "HiCProOut"
        String outdir = "PeakOut"+ FDR
        Int ncpu = 4
    }

    File actual_genomefragment = select_first([genomefragment, ""])
    Int memory_gb = if (count_fastqs < 5) then ceil(5 * 3) else ceil(count_fastqs * 3)

    command <<<
        pwd=$(pwd)
        mkdir -p ~{outdir} ~{pp_directory}
        tar -xpf ~{hicpro_result}

        # Allow for Genome Fragment to not be specified
        if [ -f ~{actual_genomefragment} ]; then
            fragment_filename=~{basename(actual_genomefragment, ".gz")}
            GFragment=$pwd/$fragment_filename
            gunzip -c ~{actual_genomefragment} > $GFragment \
                || ln -sf ~{actual_genomefragment} $GFragment
            genomefragment=$GFragment
        else
            genomefragment=""
        fi

        echo this is genomefragment $genomefragment
        peak_call -i ~{hicpro_out}/hic_results/data/fastq -o ~{outdir} -r $genomefragment -f ~{FDR} -a ~{chromsizes}

        cd ~{outdir}
        cp fastqbedgraph.bdg ~{sampleid_m}fastqbedgraph.bdg
        cp fastqpeaks.bed ~{sampleid_m}fastqpeaks.bed
        cp fastqlog.log ~{sampleid_m}fastqlog.log
        cp fastqreport.pdf ~{sampleid_m}fastqreport.pdf
        LC_COLLATE=C sort -k1,1 -k2,2n ~{sampleid_m}fastqbedgraph.bdg > ~{sampleid_m}fastqbedgraph.sorted.bdg
        bedGraphToBigWig ~{sampleid_m}fastqbedgraph.sorted.bdg ~{chromsizes} ~{sampleid_m}fastqbedgraph.bw

        #proteinpaint directory
        cp ~{sampleid_m}fastqbedgraph.bw $pwd/~{pp_directory}

    >>>

    output {
        File oneDpeaksbw = "~{outdir}/~{sampleid_m}fastqbedgraph.bw"
        File oneDpeaksbdg = "~{outdir}/~{sampleid_m}fastqbedgraph.bdg"
        File oneDpeaksbed = "~{outdir}/~{sampleid_m}fastqpeaks.bed"
        File oneDpeakslog = "~{outdir}/~{sampleid_m}fastqlog.log"
        File oneDpeaksreport = "~{outdir}/~{sampleid_m}fastqreport.pdf"
        File pp_bw = "~{pp_directory}/~{sampleid_m}fastqbedgraph.bw"
    }

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: ncpu
    }
}

task filterblklist {
    meta {
        description: "Filter HiC-Pro results based on BlackList regions"
        outputs: {
            hicpro_filtered: "Filtered HiC-Pro results archive",
        }
    }

    parameter_meta {
        count_fastqs: "Number of FASTQ files"
        hicpro_result: "HiC-Pro result archive"
        blacklist: "List of regions to filter out"
        chromsizes: "Chromosome sizes file"
        slop: "Slop size for filtering regions"
        hicpro_out: "Output directory for HiC-Pro results"
    }

    # Filter based on BlackList regions
    input {
        Int count_fastqs
        File hicpro_result
        File blacklist
        File chromsizes
        Int slop = 50
        String hicpro_out = "HiCProOut"
    }

    Int memory_gb = ceil((count_fastqs * 3) + 10)
    String base = basename(blacklist, ".gz")

    command <<<
        set -euo pipefail
        pwd=$(pwd)
        tar -xpf ~{hicpro_result}

        all_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.allValidPairs
        orig_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.AVP.Orig
        mv $all_valid_pairs $orig_valid_pairs

        cd $pwd/~{hicpro_out}/hic_results/data/fastq

        gunzip -c ~{blacklist} > ~{base} \
           || ln -sf ~{blacklist} ~{base}

        #left bed
        awk -F "\t" -v OFS="\t" '{print $2,$3,$3,$1}' $orig_valid_pairs|slopBed -i stdin -b ~{slop} -g ~{chromsizes}|intersectBed -a stdin -b ~{base} -u > left.bed

        #right bed
        awk -F "\t" -v OFS="\t" '{print $5,$6,$6,$1}' $orig_valid_pairs|slopBed -i stdin -b ~{slop} -g ~{chromsizes}|intersectBed -a stdin -b ~{base} -u >right.bed


        cat <(cut -f 4 left.bed) <(cut -f 4 right.bed)|sort -u > filter.pair

        python <<CODE
        from collections import defaultdict
        f=open("filter.pair")
        blackID=defaultdict(int)
        while True:
            line=f.readline()
            if not line:
                break
            cols=line.strip().split("\t")
            ID=cols[0]
            blackID[ID]=1

        f2=open("fastq.AVP.Orig")
        outfile1="fastq.allValidPairs"
        outfile2="fastq.AVP.removed"
        of1=open(outfile1,'w')
        of2=open(outfile2,'w')
        while True:
            line=f2.readline()
            if not line:
                break
            cols=line.strip().split("\t")
            id=cols[0]
            if id in blackID:
                of2.write(line)
            else:
                of1.write(line)
        CODE

        all=`wc -l $orig_valid_pairs |cut -d " " -f 1`
        filtered=`wc -l $pwd/~{hicpro_out}/hic_results/data/fastq/fastq.AVP.removed|cut -d " " -f 1`
        percentage=`echo "$filtered/$all*100"|bc -l`
        echo -e "$filtered read pairs are filtered out from a total of $all, accounting for ${percentage} %\n"

        cd $pwd
        tar -cpf ~{basename(hicpro_out)}.tar ~{hicpro_out}
    >>>

    output {
        File hicpro_filtered = "~{basename(hicpro_out)}.tar"
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task twoDloops {
    meta {
        description: "Call 2D loops from HiC-Pro result using FitHiChIP"
        outputs: {
            twoDloops_zip: "Output zip file containing FitHiChIP results",
            pp_peaks: "ProteinPaint peaks file in BigWig format",
            pp_tbi: "ProteinPaint peaks index file in TBI format",
            pp_peaks2: "ProteinPaint peaks file in BigWig format for second set of loops",
            pp_tbi2: "ProteinPaint peaks index file in TBI format for second set of loops",
            anchor_bw: "Anchor BigWig file for FitHiChIP results",
            bedloop_1: "Bed file for the first set of loops",
            bedloop_2: "Bed file for the second set of loops",
            bedloop_3: "Bed file for the third set of loops",
            tbiloop_1: "TBI index file for the first set of loops",
            tbiloop_2: "TBI index file for the second set of loops",
            tbiloop_3: "TBI index file for the third set of loops",
            qc_output: "QC output zip file containing FitHiChIP results",
        }
    }

    parameter_meta {
        count_fastqs: "Number of FASTQ files"
        hicpro_result: "HiC-Pro result archive"
        genomename: "Name of the genome"
        FDR: "False discovery rate"
        IntType: "Interaction type"
        LowDistThr: "Lower distance threshold"
        UppDistThr: "Upper distance threshold"
        chromsizes: "Chromosome sizes file"
        loopBed: "Loop regions BED file"
        input_loopBed: "Input loop BED file"
        bw_loopBed: "BigWig file for loop regions"
        hicpro_out: "HiC-Pro output directory"
        loopsAnchor: "File containing loop anchors"
        promoters: "File containing promoter regions"
        loopThreshold_1: "Loop threshold for the first set of loops"
        loopThreshold_2: "Loop threshold for the second set of loops"
        loopThreshold_3: "Loop threshold for the third set of loops"
        loopThresholds: "Array of loop thresholds"
        pp_directory: "Directory for ProteinPaint output"
        sampleid: "Sample ID for the output files"
        sampleid_m: "Modified sample ID for output files"
        loopName: "Name for the loop output files"
        loopOut: "Output directory for FitHiChIP results"
        qcfiles: "QC files directory"
    }

    # Call 2D loops from  HiC-Pro result using FitHiChIP
    input {
        Int count_fastqs
        File hicpro_result
        String genomename
        Float FDR = 0.4
        Int IntType = 3
        Int LowDistThr = 5000
        Int UppDistThr = 4000000
        File chromsizes
        File loopBed
        String input_loopBed = basename(loopBed)
        String bw_loopBed = basename(loopBed,".bed")
        String hicpro_out = "HiCProOut"
        File? loopsAnchor
        File? promoters
        Float loopThreshold_1 = 0.01
        Float loopThreshold_2 = 0.05
        Float loopThreshold_3 = 0.1
        Array[Float] loopThresholds = [loopThreshold_1,loopThreshold_2,loopThreshold_3]
        String pp_directory
        String? sampleid
        String sampleid_m = "~{sampleid + '.' + genomename + '_'}"
        String loopName = if defined(loopsAnchor) then IntType + "." + FDR + "." + basename(loopBed,".bed") + "." + LowDistThr + "." + UppDistThr else if defined(promoters) then IntType + "." + FDR + "." + basename(loopBed,".bed") + "." + LowDistThr + "." + UppDistThr else IntType + "." + FDR + "." + LowDistThr + "." + UppDistThr
        String loopOut = if defined(loopsAnchor) then "LoopOut_FitHiChIP." + IntType + "." + FDR + "." + basename(loopBed,".bed") + "." + LowDistThr + "." + UppDistThr else if defined(promoters) then "LoopOut_FitHiChIP." + IntType + "." + FDR + "." + basename(loopBed,".bed") + "." + LowDistThr + "." + UppDistThr else "LoopOut_FitHiChIP." + IntType + "." + FDR + "." + LowDistThr + "." + UppDistThr
        String qcfiles = if defined(promoters) then "active_QCfiles" else "QCfiles"
    }

    Int memory_gb = ceil(count_fastqs * 3)

    command <<<
        pwd=$(pwd)
        mkdir -p ~{loopOut}/Anchor ~{pp_directory}
        tar -xpf ~{hicpro_result}

        InputFile=~{input_loopBed}
        Flag=${InputFile##*.}

        if [ -f "~{loopsAnchor}" ]; then
            cd $pwd/~{pp_directory}
            awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' ~{loopsAnchor} > ~{bw_loopBed}'.bdg'

            LC_COLLATE=C sort -k1,1 -k2,2n ~{bw_loopBed}'.bdg' > ~{bw_loopBed}'.sorted.bdg'
            bedGraphToBigWig ~{bw_loopBed}'.sorted.bdg' ~{chromsizes} ~{bw_loopBed}'.bw'
        fi

        cd $pwd/~{loopOut}
        provided_sampleid = ~{sampleid_m}
        if [[ -n $provided_sampleid ]] && [[ $InputFile != "~{sampleid_m}"* ]]; then
            Output=Anchor/~{sampleid_m}$InputFile
        else
            Output=Anchor/$InputFile
        fi

        if [[ $Flag == "bedpe" ]];then
            awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' ~{loopBed} | sort -k1,1 -k2,2n > $Output
        elif [[ $Flag == "bed" ]];then
            sort -k1,1 -k2,2n ~{loopBed} |cut -f 1-3  > $Output
        elif [[ $Flag == "bw" ]];then
            cp ~{loopBed} $Output
        fi

        bgzip -f $Output
        tabix -f -p bed $Output.gz

        # Copy to proteinpaint directory
        cp $Output.gz{,.tbi} $pwd/~{pp_directory}

        config=~{sampleid + "_"}config.txt

        cat > $config <<EOF
        #Interval==$pwd/~{hicpro_out}/hic_results/matrix/fastq/raw/5000/fastq_5000_abs.bed
        #Matrix=$pwd/~{hicpro_out}/hic_results/matrix/fastq/raw/5000/fastq_5000.matrix
        ValidPairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.allValidPairs
        Bed=
        PeakFile=~{loopBed}
        OutDir=$pwd/~{loopOut}
        IntType=~{IntType}
        BINSIZE=5000
        LowDistThr=~{LowDistThr}
        UppDistThr=~{UppDistThr}
        UseP2PBackgrnd=0
        BiasType=1
        MergeInt=1
        QVALUE=0.01
        ChrSizeFile=~{chromsizes}
        PREFIX=FitHiChIP
        OverWrite=1
        EOF

        # Run FitHiChIP
        /FitHiChIP/FitHiChIP_HiCPro.sh -C $config

        if [ ~{IntType} -eq 3 ];then
            LoopDir=$pwd/~{loopOut}/FitHiChIP_Peak2ALL_b5000_L~{LowDistThr}_U~{UppDistThr}/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr
        elif [ ~{IntType} -eq 1 ];then
            LoopDir=$pwd/~{loopOut}/FitHiChIP_Peak2Peak_b5000_L~{LowDistThr}_U~{UppDistThr}/Coverage_Bias/FitHiC_BiasCorr
        else
            echo "Error!! IntType must be 1 (Peak2Peak) or 3 (Peak2All)"
            exit 1
        fi

        cd $LoopDir
        LoopThr_list=(~{sep=" " loopThresholds}) && count=~{length(loopThresholds)}

        LoopAll=$LoopDir/FitHiChIP.interactions_FitHiC.bed

        echo -e "Get significant loops on $count FDR threholds"

        for Qvalue in ${LoopThr_list[@]}; do
            Output=$LoopDir/~{sampleid_m}"Qvalue."$Qvalue.~{loopName}.bed
            Output2=$LoopDir/~{sampleid_m}"Qvalue."$Qvalue.~{loopName}.bedpe
            awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7}' $LoopAll | sort -k1,1 -k2,2n > $Output2
            awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' $LoopAll | sort -k1,1 -k2,2n > $Output
            bgzip -f $Output
            tabix -f -p bed $Output.gz

            # Copy to specified Proteinpaint directory name
            cp $Output.gz{,.tbi} $pwd/~{pp_directory}

            echo -e "processed Qvalue $Qvalue and stored results at $Output"
        done

        cd $pwd
        zip -1qr ~{loopOut}.zip ~{loopOut}

        mkdir -p ~{qcfiles}
        cp $pwd/~{hicpro_out}/hic_results/stats/fastq/* ~{qcfiles}
        cp ~{loopBed} $LoopAll $LoopDir/FitHiChIP.interactions_FitHiC_Q0.01.bed ~{qcfiles}
        zip -9qr ~{qcfiles}.zip ~{qcfiles}

    >>>

    output {
        File twoDloops_zip = "~{loopOut}.zip"
        File? pp_peaks = "~{pp_directory}/~{sampleid_m}~{input_loopBed}.gz"
        File? pp_tbi = "~{pp_directory}/~{sampleid_m}~{input_loopBed}.gz.tbi"
        File? pp_peaks2 = "~{pp_directory}/~{input_loopBed}.gz"
        File? pp_tbi2 = "~{pp_directory}/~{input_loopBed}.gz.tbi"
        File? anchor_bw = "~{pp_directory}/~{bw_loopBed}.bw"
        File? bedloop_1 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_1}.~{loopName}.bed.gz"
        File? bedloop_2 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_2}.~{loopName}.bed.gz"
        File? bedloop_3 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_3}.~{loopName}.bed.gz"
        File? tbiloop_1 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_1}.~{loopName}.bed.gz.tbi"
        File? tbiloop_2 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_2}.~{loopName}.bed.gz.tbi"
        File? tbiloop_3 = "~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_3}.~{loopName}.bed.gz.tbi"
        File? qc_output = "~{qcfiles}.zip"
    }

    runtime {
        memory: ceil(memory_gb * 3) + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 4
    }
}

task activeRegionsMerge {
    meta {
        description: "Merge promoter and enhancer regions to create active coordinate regions"
        outputs: {
            combined_promoters: "Combined promoter and enhancer regions in BED format",
            merge_bam: "Merged BAM file after removing duplicates",
        }
    }

    parameter_meta {
        hicpro_result: "HiC-Pro result file in tar format"
        hicpro_out: "Output directory for HiC-Pro results"
        promoters: "Promoter regions in BED format"
        loopBed: "Loop regions in BED format"
        genomename: "Name of the genome used in the analysis"
        sampleid: "Sample ID for the data"
        sampleid_m: "Modified Sample ID for output files"
        comb_promoters: "Combined promoter file name"
        outdir: "Output directory for active coordinate regions"
        memory_gb: "Memory allocated for the task in GB"
        ncpu: "Number of CPUs allocated for the task"
    }

    input {
        File hicpro_result
        String hicpro_out = "HiCProOut"
        File promoters
        File loopBed
        String genomename
        String? sampleid
        String sampleid_m = "~{sampleid + '.' + genomename + '_'}"
        String comb_promoters = sub(basename(promoters),".bed",".LoopAnchors.Enhs.combined.sort.bed")
        String outdir = "ActiveCoordRegionFiles"
        Int memory_gb = 20
        Int ncpu = 8
    }

    command <<<
        pwd=$(pwd)
        tar -xpf ~{hicpro_result}
        BamDir=$pwd/~{hicpro_out}/hic_results/data/fastq

        mkdir -p ~{outdir}
        cd ~{outdir}
        #merge BAM
        samtools merge -@ ~{ncpu} ~{sampleid_m}merged.bam $BamDir/*.bam

        #rmdup BAM
        samtools rmdup -S ~{sampleid_m}merged.bam ~{sampleid_m}merged.rmdup.bam

        prom=~{basename(promoters,".bed")}
        sort -k1,1 -k2,2n -k4,4 ~{promoters} |awk -F'\t' -v OFS="\t" '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$4 }END{ for(i in a) print i,a[i] }'| sort -k1,1 -k2,2n > $prom".unique.bed"

        #get read counts
        intersectBed -c -a $prom".unique.bed" -b ~{sampleid_m}"merged.rmdup.bam" > $prom"_Signal"

        #compute top 2/3
        length=$(wc -l $prom"_Signal" | awk -F\  '{print $1}'); 
        keepN=$(echo "$(($length * 2 / 3))"); 

        sort -k5,5nr $prom"_Signal" | head -n $keepN | sort -k1,1 -k2,2n | cut -f 1-3 > $prom"_Signal.keep."$keepN".3cols.bed"

        #get non-promoter region
        cut -f 1-3 ~{loopBed} | sort -k1,1 -k2,2n > ~{basename(loopBed)}
        bedtools subtract -a ~{basename(loopBed)} -b $prom".unique.bed" | sort > Enhancers.Subtract.4kbPromoters.sort.bed

        cat $prom"_Signal.keep."$keepN".3cols.bed" Enhancers.Subtract.4kbPromoters.sort.bed | sort -k1,1 -k2,2n > ~{comb_promoters}

    >>>

    output {
        File combined_promoters = "~{outdir}/~{comb_promoters}"
        File merge_bam = "~{outdir}/~{sampleid_m}merged.rmdup.bam"
    }

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: ncpu
    }
}

task converthic {
    meta {
        description: "Convert HiC-Pro results to .hic format for visualization"
        outputs: {
            hic_file: "Converted HiC file in .hic format",
            pp_hic: "ProteinPaint compatible HiC file",
        }
    }

    parameter_meta {
        count_fastqs: "Number of FASTQ files processed"
        hicpro_result: "HiC-Pro result file in tar format"
        chromsizes: "Chromosome sizes file in tab delimited format"
        hicpro_out: "Output directory for HiC-Pro results"
        hic_out: "Output directory for .hic files"
        genomename: "Name of the genome used in the analysis"
        pp_directory: "Directory for ProteinPaint output"
        sampleid: "Sample ID for the data (optional)"
        sampleid_m: "Modified Sample ID for output files (default: '~{sampleid + '.' + genomename + '_'}')"
    }

    # Convert HiC-Pro result to .hic file
    input {
        Int count_fastqs
        File hicpro_result
        File chromsizes
        String hicpro_out = "HiCProOut"
        String hic_out = "HicFile"
        String genomename
        String pp_directory
        String? sampleid
        String sampleid_m = "~{sampleid + '.' + genomename + '_'}"
    }

    Int memory_gb = ceil(count_fastqs * 15)

    command <<<
        pwd=$(pwd)
        mkdir -p ~{hic_out} ~{pp_directory}
        tar -xpf ~{hicpro_result}
        all_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.allValidPairs

        cd ~{hic_out}
        /HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i $all_valid_pairs -g ~{chromsizes} -j /usr/local/bin/juicer_tools_1.22.01.jar

        # Copy to proteinpaint directory
        cp fastq.allValidPairs.hic $pwd/~{pp_directory}/~{sampleid_m}fastq.allValidPairs.hic
    >>>

    output {
        File hic_file = "~{hic_out}/fastq.allValidPairs.hic"
        File pp_hic = "~{pp_directory}/~{sampleid_m}fastq.allValidPairs.hic"
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task bowtie2_index {
    meta {
        description: "Generate Bowtie2 indexes for a reference genome"
        outputs: {
            bowtie_indexes: "Array of Bowtie2 index files",
        }
    }

    parameter_meta {
        reference: "Reference genome file in FASTA format (can be gzipped)"
        genomename: "Name of the genome for the index"
        memory_gb: "Memory allocated for the task in GB"
    }

    # Generate bowtie2 indexes
    input {
        File reference
        String genomename
        Int memory_gb = 20
    }

    String base = basename(reference, ".gz")

    command <<<
        gunzip -c ~{reference} > ~{base} \
           || ln -sf ~{reference} ~{base}

        bowtie2-build --threads 20 ~{base} ~{genomename}
    >>>

    output {
        Array[File] bowtie_indexes = glob("*.bt2")
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 20
    }
}

task chrfaidx {
    meta {
        description: "Generate chromosome sizes and faidx file from a reference genome"
        outputs: {
            faidx_file: "FAIDX file for the reference genome",
            chromsizes: "Chromosome sizes in tab delimited format",
        }
    }

    parameter_meta {
        reference: "Reference genome file in FASTA format (can be gzipped)"
        memory_gb: "Memory allocated for the task in GB"
    }

    input {
	    File reference
        Int memory_gb = 5
    }

    String base = basename(reference, ".gz")

    command <<<
        gunzip -c ~{reference} > ~{base} \
           || ln -sf ~{reference} ~{base}

        samtools faidx ~{base} -o ~{base}.fai
        cut -f1,2 ~{base}.fai > ~{base}.tab
    >>>

    output {
        File faidx_file = "~{base}.fai"
        File chromsizes = "~{base}.tab"
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task raise_exception {
    meta {
        description: "Utility task to raise an exception with a custom error message"
        outputs: {
            error_msg: "Error message indicating the reason for failure",
        }
    }

    parameter_meta {
        msg: "Custom error message to be displayed"
    }

    input {
        String msg
    }

    command <<<
        echo -e "\n* Error: ~{msg}\n" >&2
        exit 2
    >>>

    output {
        String error_msg = "~{msg}"
    }

    runtime {
        memory: "5 GB"
        maxRetries: 0
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task createjson {
    meta {
        description: "Create a JSON file for ProteinPaint visualization"
        outputs: {
            json: "JSON file for ProteinPaint",
        }
    }

    parameter_meta {
        genomename: "Name of the genome"
        sampleid: "Sample ID for the data"
        sampleid_m: "Modified Sample ID for output files"
        pp_directory: "Directory for ProteinPaint output"
        IntType: "Type of interaction (1 for Peak2Peak, 3 for Peak2All)"
        pp_bw: "BigWig file for coverage"
        pp_hic: "HiC file for visualization"
        pp_peaks: "Peaks file in BED format"
        bedloop_1: "First loop file in BED format"
        bedloop_2: "Second loop file in BED format"
        bedloop_3: "Third loop file in BED format"
    }

    input {
        String genomename
        String? sampleid
        String sampleid_m = "~{sampleid + '.' + genomename + '_'}"
        String pp_directory
        Int IntType = 3
        File? pp_bw
        File pp_hic
        File? pp_peaks
        File? bedloop_1
        File? bedloop_2
        File? bedloop_3
    }

    Array[File] bedloops = select_all([bedloop_1,bedloop_2,bedloop_3])

    command <<<
        mkdir -p ~{pp_directory}
        cd ~{pp_directory}

        cat > ~{sampleid_m}proteinpaint.json <<EOF
        [
        EOF
        cat >> ~{sampleid_m}proteinpaint.json <<EOF
            {
                "type":"hicstraw",
                "name":"~{sampleid + '-'}hic",
                "mode_arc":false,
                "mode_hm":true,
                "file":"~{pp_directory}/~{basename(pp_hic)}"
        EOF

        # if there is peak data
        pp_bw=~{pp_bw}
        pp_peaks=~{pp_peaks}
        if [ -f "~{pp_peaks}" ]; then
            cat >> ~{sampleid_m}proteinpaint.json <<EOF
            },
            {
                "type":"bigwig",
                "name":"~{sampleid + '-'}Coverage",
                "file":"~{pp_directory}/${pp_bw##*/}",
                "scale":{
                    "auto": 1
                },
                "height":50
            },
            {
                "type":"bedj",
                "name":"~{sampleid + '-'}Anchor",
                "file":"~{pp_directory}/${pp_peaks##*/}",
                "stackheight":14,
                "stackspace":1
        EOF
        fi

        # Looping files
        for loop in ~{sep=" " bedloops};
        do
            newloop=${loop##*Qvalue.}
            Qvalue=(${newloop//.~{IntType}/ })
            cat >> ~{sampleid_m}proteinpaint.json <<EOF
            },
            {
                "type":"hicstraw",
                "name":"~{sampleid + '-'}$Qvalue",
                "mode_arc":true,
                "mode_hm":false,
                "bedfile":"~{pp_directory}/${loop##*/}"
        EOF
        done

        # hic file
        cat >> ~{sampleid_m}proteinpaint.json <<EOF
            }
        ]
        EOF


    >>>

    output {
        File json = "~{pp_directory}/~{sampleid_m}proteinpaint.json"
    }

    runtime {
        memory: "5 GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task compressfiles {
    meta {
        description: "Compress input files into a gzipped format"
        outputs: {
            zipfile: "Gzipped file of the input file",
        }
    }

    parameter_meta {
        inputfile: "Input file to be compressed"
    }

    input {
        File inputfile
    }

    command <<<
        cp ~{inputfile} ~{basename(inputfile)}
        gzip ~{basename(inputfile)}
    >>>

    output {
        File zipfile = "~{basename(inputfile)}.gz"
    }

    runtime {
        memory: "5 GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}

task fastqc {
    meta {
        description: "Run FastQC on a given FASTQ file"
        outputs: {
            htmlfile: "HTML file containing FastQC report",
            zipfile: "ZIP file containing FastQC results",
        }
    }

    parameter_meta {
        fastq: "Input FASTQ file"
        prefix: "Prefix for output files"
        default_location: "Default location for FastQC output files"
        memory_gb: "Memory in GB allocated for the task"
    }

    input {
        File fastq
        String prefix = sub(basename(fastq),".fastq.gz|.fq.gz","")
        String default_location="FastQC"
        Int memory_gb = 5
    }

    command <<<
        ln -s ~{fastq} ~{sub(basename(fastq),".bam$",".bam.bam")}

        echo ~{prefix}
        mkdir -p ~{default_location}

        fastqc \
            -o ~{default_location} \
            ~{sub(basename(fastq),".bam$",".bam.bam")}
    >>>

    output {
        File htmlfile = "~{default_location}/~{prefix}_fastqc.html"
        File zipfile = "~{default_location}/~{prefix}_fastqc.zip"
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/fastqc:v0.11.9"
        cpu: 1
    }
}

task qcreport {
    meta {
        description: "Generate QC report from HiC-Pro output"
        outputs: {
            qcreport: "Text file containing QC statistics from HiC-Pro output"
        }
    }

    parameter_meta {
        qczip: "Input zip file containing HiC-Pro QC files"
        activeregion: "Boolean to indicate if active region analysis is performed"
        genomename: "Name of the genome used"
        sampleid: "Sample ID for naming output files"
        sampleid_m: "Modified sample ID for output file naming"
        memory_gb: "Memory in GB allocated for the task"
    }

    input {
        File qczip
        Boolean activeregion = false
        String genomename
        String? sampleid
        String sampleid_m = if activeregion
            then "active_"
            else "" + "~{sampleid + "." + genomename + "_"}"
        Int memory_gb = 10
    }

    command <<<
        pwd=$(pwd)
        unzip "~{qczip}"
        cd "~{basename(qczip,".zip")}"

        python <<CODE
        import os
        FILES = ["fastq_allValidPairs.mergestat", "fastq.mpairstat", "fastq_R1.mmapstat", "fastq_R2.mmapstat"]
        VARIABLE = ["valid_interaction", "cis_shortRange", "cis_longRange", "Total_pairs_processed", "Reported_pairs", "total_R2", "mapped_R2", "total_R1", "mapped_R1"]
        RESULTS = {}
        for eachfile in FILES:
            with open(eachfile) as f:
                for line in f:
                    array = line.split()
                    if array[0] in VARIABLE:
                        RESULTS[array[0]] = int(array[1])

        PERCENTAGES = {}
        COUNTS = {}
        VERDICT = {}
        REP_VAR = ["R1_aligned", "R2_aligned", "valid_interactionPairs", "cis_shortRange", "cis_longRange"]

        if os.path.isfile('FitHiChIP.interactions_FitHiC.bed'):
            with open('FitHiChIP.interactions_FitHiC_Q0.01.bed') as fithichip:
                LOOPS_SIGNIFICANT = len(fithichip.readlines())-1
            with open('FitHiChIP.interactions_FitHiC.bed') as fithichip:
                LOOPS = len(fithichip.readlines())-1
        if os.path.isfile("~{sampleid_m}fastqpeaks.bed"):
            with open("~{sampleid_m}fastqpeaks.bed") as peaks:
                PEAKS = len(peaks.readlines())-1

        PERCENTAGES["R1_aligned"] = round((RESULTS["mapped_R1"]*100)/RESULTS["total_R1"])
        PERCENTAGES["R2_aligned"] = round((RESULTS["mapped_R2"]*100)/RESULTS["total_R2"])
        PERCENTAGES["valid_interactionPairs"] = round((RESULTS["valid_interaction"]*100)/RESULTS["Reported_pairs"])
        PERCENTAGES["cis_shortRange"] = round((RESULTS["cis_shortRange"]*100)/RESULTS["valid_interaction"])
        PERCENTAGES["cis_longRange"] = round((RESULTS["cis_longRange"]*100)/RESULTS["valid_interaction"])
        COUNTS["R1_aligned"] = RESULTS["mapped_R1"]
        COUNTS["R2_aligned"] = RESULTS["mapped_R2"]
        COUNTS["valid_interactionPairs"] = RESULTS["valid_interaction"]
        COUNTS["cis_shortRange"] = RESULTS["cis_shortRange"]
        COUNTS["cis_longRange"] = RESULTS["cis_longRange"]

        for aligned in ["R1_aligned", "R2_aligned"]:
            if PERCENTAGES[aligned] > 80:
                VERDICT[aligned] = "GOOD"
            else:
                VERDICT[aligned] = "BAD"

        if PERCENTAGES["valid_interactionPairs"] > 50:
            VERDICT["valid_interactionPairs"] = "GOOD"
        else:
            VERDICT["valid_interactionPairs"] = "BAD"

        if PERCENTAGES["cis_shortRange"] > 50:
            VERDICT["cis_shortRange"] = "BAD"
        elif PERCENTAGES["cis_shortRange"] > 30:
            VERDICT["cis_shortRange"] = "MARGINAL"
        else:
            VERDICT["cis_shortRange"] = "GOOD"

        if PERCENTAGES["cis_longRange"] > 40:
            VERDICT["cis_longRange"] = "GOOD"
        elif PERCENTAGES["cis_longRange"] > 20:
            VERDICT["cis_longRange"] = "MARGINAL"
        else:
            VERDICT["cis_longRange"] = "BAD"

        REPORT = open("~{sampleid_m}QCreport.txt", 'w')
        REPORT.write("STAT\tCOUNTS\tPERCENTAGE\tVERDICT\n")
        REPORT.write("Total_pairs_processed\t" + str(RESULTS["Total_pairs_processed"]) + "\n")
        for var in REP_VAR:
            REPORT.write(var + "\t" + str(COUNTS[var]) + "\t" + str(PERCENTAGES[var]) + "\t" + VERDICT[var] + '\n')

        if os.path.isfile("~{sampleid_m}fastqpeaks.bed"):
            REPORT.write("peaks\t" + str(PEAKS) + "\n")
        if os.path.isfile('FitHiChIP.interactions_FitHiC.bed'):
            REPORT.write("loops\t" + str(LOOPS) + "\n")
            REPORT.write("loops_significant\t" + str(LOOPS_SIGNIFICANT) + "\n")
        CODE
        cp "~{sampleid_m}QCreport.txt" "$pwd"
    >>>

    output {
        File qcreport = "~{sampleid_m}QCreport.txt"
    }

    runtime {
        memory: memory_gb + " GB"
        maxRetries: 1
        docker: "ghcr.io/stjude/abralab/hilow:v1.0"
        cpu: 1
    }
}
