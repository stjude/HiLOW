version 1.0

workflow hilow {
    String pipeline_ver = 'v1.0.0'

    meta {
        title: 'StJude HiC Looping Organization Workflow'
        summary: 'HiC parallel data processing, Peak and loop Calling Pipeline'
        description: 'An all-in-one pipeline which takes FASTQs and perform mapping by HiC-Pro, filtering, merging, normalization, converting to .hic files, peaks and loops analysis.'
        version: '1.0'
        details: {
            citation: '',
            contactEmail: ['jie.lu@stjude.org','brian.abraham@stjude.org','modupeore.adetunji@stjude.org'],
            contactOrg: "St Jude Children's Research Hospital",
            contactUrl: "",
            upstreamLicenses: "MIT",
            upstreamUrl: '',
            whatsNew: ''
        }
    }

    input {
        # Required Inputs
        Array[File] fastqfiles_R1
        Array[File] fastqfiles_R2
        File? genomefragment
        String genomename = 'hg19'

        # Optional Inputs
        String pp_directory = "browser_files"
        String? sampleid
        File? reference
        File? blacklist
        File? chromsizes
        Array[File]? bowtie2index

        # HiC and HiChIP mapping using HiC-Pro
        String? ligationsite

        # Limit HiC analysis
        Boolean ashichip = true

        # Peak calling parameters using HiChIP_peaks
        Float FDR=0.4

        # Loop calling parameters using FitHiChip
        File? loopsAnchor
        Int LowerThreshold=5000
        Int UpperThreshold=4000000 
        Int IntType=3
        String loopThreshold='0.01,0.05,0.1'
        Float loopThreshold_1=0.01
        Float loopThreshold_2=0.05
        Float loopThreshold_3=0.1
    }

    parameter_meta {
        fastqfiles_R1: {
            description: 'One or more sample R1 FASTQs',
            help: 'Upload zipped R1 FASTQ files.',
            patterns: ["*_R1.fq.gz", "*_R1.fastq.gz"]
        }
        fastqfiles_R2: {
            description: 'One or more sample R2 FASTQs',
            help: 'Upload zipped R2 FASTQ files.',
            patterns: ["*_R2.fq.gz", "*_R2.fastq.gz"]
        }
        genomefragment: {
            description: 'restriction fragment file in BED format',
            help: 'Input genome fragments digested from restriction enzyme(s).',
            patterns: ["*.bed", "*.bed.gz"]
        }
        genomename: {
            description: 'Reference Genome Name',
            help: 'Reference genome name. Default: hg19. Examples: hg19, mm10, mm9',
            example: 'hg19'
        }
        hicpro_output: {
            description: 'HiC-Pro output name',
            help: 'Input preferred HiC-Pro output name.',
            example: 'HiCPro_out'
        }
        reference: {
            description: 'Reference FASTA file',
            patterns: ["*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz"]
        }
        blacklist: {
            description: 'Blacklist file in BED format',
            help: 'If defined, blacklist regions listed are excluded after reference alignment.',
            patterns: ["*.bed", "*.bed.gz"]
        }
        chromsizes: {
            description: 'Chromosome sizes file in tab delimited format',
            help: 'If defined, genome chromosome size file.',
        }
        pp_directory: {
            description: 'ProteinPaint directory',
            help: 'preferred ProteinPaint directory name, the directory contains the visualization files and required json. Directory should be copied to the ProteinPaint "*/tp" directory.'
        }
        bowtie2index: {
            description: 'bowtie v2 index files (*.bt2)',
            help: 'If not defined, bowtie v2 index files are generated, will take a longer compute time.',
            patterns: ["*.bt2"]
        }
        loopsAnchor: {
            description: 'Optional file containing pre-defined set of loop anchors, only the first 3 cols are relevent (Chr,Start, End). 1D peak calling step will skipped if provided',
            patterns: ["*.bed"]
        }
        FDR: {
            description: 'False Discovery Rate used to call Peaks using HiChIP-peaks.',
            help: 'False discovery rate'
        }
        ligationsite: {
            description: 'Restriction enzyme ligation site used by HiC-Pro'
        }
        LowerThreshold: {
            description: 'Lower distance threshold of interaction between two intervals. Interactions below this distance threshold will not be considered for statistical significance.'
        }
        UpperThreshold: {
            description: 'Upper distance threshold of interaction between two intervals. Interactions above this distance threshold will not be considered for statistical significance.'
        }
        IntType: {
            description: 'Type of interaction (foreground) reported by FitHiChIP.'
        }
        loopThreshold_1: {
            description: 'FDR threshold 1 used to call loops using FitHiChIP.',
            example: '0.01'
        }
        loopThreshold_2: {
            description: 'FDR threshold 2 used to call loops using FitHiChIP.',
            example: '0.1'
        }
        loopThreshold_3: {
            description: 'FDR threshold 3 used to call loops using FitHiChIP.',
            example: '0.01'
        }
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 1 ----------- ###
### ------ Pre-process Analysis Files ------ ###
### ---------------------------------------- ###

    # HiCPro directory
    String hicpro_output = 'HiCPro_out'
    
    # IntType is 1 or 3. Default 3
    if (IntType != 1 && IntType != 3) {
        call raise_exception as error_wrong_IntType { input: msg = 'Type of interaction reported by FitHiChIP. Choices: 1 (peak to peak) or 3 (peak to all). Default = 3'}
    }

    if (ashichip && (!defined(genomefragment) || !defined(ligationsite))){
        call raise_exception as error_missing_data { input: msg = 'Genome Fragment file and/or Ligation sites are not specified for HiChIP peak calling, set "ashichip" to False' }
    }
    
    # if (ashichip && !defined(genomefragment)) {
    #     call raise_exception as error_missing_fragments  { 
    #         input: msg = 'Genome Fragments file is not specified for HiChiP peak calling, set ashichip to False'
    #     }
    # }
    # if (ashichip && !defined(ligationsite)) {
    #     call raise_exception as error_missing_ligationsite  { 
    #         input: msg = 'Ligation sites are not specified for HiChiP peak calling, set ashichip to False'
    #     }
    # }

    # Generating INDEX files
    String string_reference = ""
    #File reference_m = select_first[(reference,string_reference)]
    #1. Bowtie INDEX files if not provided
    if ( !defined(bowtie2index) ) {
        # create bowtie index when not provided
        call bowtie2_index {
            input :
                genomename=genomename,
                reference=select_first([reference,string_reference])
        }
    }
    #2. Make sure indexes are six else build indexes
    if ( defined(bowtie2index) ) {
        # check total number of bowtie indexes provided
        Array[String] string_bowtie2index = [1] #buffer to allow for bowtie_index optionality
        Array[File] int_bowtie2index = select_first([bowtie2index, string_bowtie2index])
        if ( length(int_bowtie2index) != 6 ) {
            # create bowtie index if 6 index files aren't provided
            call bowtie2_index as bowtie2_index_2 {
                input :
                    genomename=genomename,
                    reference=select_first([reference,string_reference])
            }
        }
    }
    Array[File] actual_bowtie2_index = select_first([bowtie2_index.bowtie_indexes, bowtie2_index_2.bowtie_indexes, bowtie2index])
    
    if (! defined(chromsizes) ){
        call chrfaidx {
            # create FASTA index and chrom sizes files
            input :
                reference=select_first([reference,string_reference])
        }
    }
    File actual_chromsizes = select_first([chrfaidx.chromsizes,chromsizes])

### ---------------------------------------- ###
### ------------ S E C T I O N 2 ----------- ###
### ------ Perform HiC-Pro Alignments ------ ###
### ---------------------------------------- ###

    Int number_fastqs = length(fastqfiles_R1)
    scatter (fastqfile_R1 in fastqfiles_R1) { 
        call fastqsnumber as extrapolate_R1 { input : fastqfile=fastqfile_R1, lengthfastqs=number_fastqs }
        call splitfastqs as validsplit_R1 { input : nreads=extrapolate_R1.cal_reads, fastqfile=fastqfile_R1 } 
    }

    scatter (fastqfile_R2 in fastqfiles_R2) { 
        call fastqsnumber as extrapolate_R2 { input : fastqfile=fastqfile_R2, lengthfastqs=number_fastqs }
        call splitfastqs as validsplit_R2 { input : nreads=extrapolate_R2.cal_reads, fastqfile=fastqfile_R2 } 
    }

    Array[Pair[File, File]] all_fastqfiles = zip(flatten(validsplit_R1.outputfiles), flatten(validsplit_R2.outputfiles))

    scatter (eachfastq in all_fastqfiles) {
        call hicpro_align as step1hicpro {
            input :
                hicpro_out=hicpro_output,
                fastqfile_R1=eachfastq.left,
                fastqfile_R2=eachfastq.right,
                ligationsite=ligationsite,
                genomefragment=genomefragment,
                genomename=genomename,
                chromsizes=actual_chromsizes,
                bowtie2index=actual_bowtie2_index
        }
    }
    call hicpro_merge as step2hicpro {
        input :
            hicpro_out=hicpro_output,
            fastqfiles_R1=flatten(validsplit_R1.outputfiles),
            fastqfiles_R2=flatten(validsplit_R2.outputfiles),
            hicpro_align=step1hicpro.hicpro_align,
            ligationsite=ligationsite,
            genomefragment=genomefragment,
            genomename=genomename,
            chromsizes=actual_chromsizes,
            bowtie2index=actual_bowtie2_index
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 3 ----------- ###
### ------- Filter Blacklist Regions ------- ###
### ---------------------------------------- ###

    if ( defined(blacklist) ) {
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_file = select_first([blacklist, string_blacklist])

        call filterblklist {
            input :
                count_fastqs=step2hicpro.count_fastqs,
                hicpro_out=hicpro_output,
                hicpro_result=step2hicpro.hicpro_merge,
                blacklist=blacklist_file,
                chromsizes=actual_chromsizes
        }
    }

    File final_hicpro = select_first([filterblklist.hicpro_filtered, step2hicpro.hicpro_merge])

### ---------------------------------------- ###
### ------------ S E C T I O N 4 ----------- ###
### ----------- Create .hic file ----------- ###
### ---------------------------------------- ###

    call converthic {
        input :
            count_fastqs=step2hicpro.count_fastqs,
            hicpro_out=hicpro_output,
            pp_directory=pp_directory,
            hicpro_result=final_hicpro,
            chromsizes=actual_chromsizes,
            genomename=genomename,
            sampleid=sampleid
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 5 ----------- ###
### ------------- Peak Calling ------------- ###
### ---------------------------------------- ###

    String string_peaks = "" #buffer to allow for not running peaks and looping analysis
    if (ashichip) {
        if ( !defined(loopsAnchor) ){
            call oneDpeaks {
                input :
                    count_fastqs=step2hicpro.count_fastqs,
                    FDR=FDR,
                    hicpro_out=hicpro_output,
                    pp_directory=pp_directory,
                    hicpro_result=final_hicpro,
                    genomefragment=genomefragment,
                    chromsizes=actual_chromsizes,
                    genomename=genomename,
                    sampleid=sampleid
            }
        }

### ---------------------------------------- ###
### ------------ S E C T I O N 6 ----------- ###
### ------------- Loop Calling ------------- ###
### ---------------------------------------- ###

        call twoDloops {
            input :
                count_fastqs=step2hicpro.count_fastqs,
                FDR=FDR,
                hicpro_out=hicpro_output,
                pp_directory=pp_directory,
                loopBed=select_first([loopsAnchor, oneDpeaks.oneDpeaksbed,string_peaks]),
                loopsAnchor=loopsAnchor,
                hicpro_result=final_hicpro,
                chromsizes=actual_chromsizes,
                LowDistThr=LowerThreshold,
                UppDistThr=UpperThreshold,
                IntType=IntType,
                loopThreshold_1=loopThreshold_1,
                loopThreshold_2=loopThreshold_2,
                loopThreshold_3=loopThreshold_3,
                genomename=genomename,
                sampleid=sampleid
        }

        call createjson as hichipjson {
            input:
                pp_bw=select_first([oneDpeaks.pp_bw,twoDloops.anchor_bw]),
                pp_hic=converthic.pp_hic,
                pp_peaks=twoDloops.pp_peaks,
                bedloop_1=twoDloops.bedloop_1,
                bedloop_2=twoDloops.bedloop_2,
                bedloop_3=twoDloops.bedloop_3,
                pp_directory=pp_directory,
                IntType=IntType,
                genomename=genomename,
                sampleid=sampleid
        }
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 7 ----------- ###
### ----------- Create JSON File ----------- ###
### ---------------------------------------- ###

    if (! ashichip) {
        call createjson as hicjson {
            input:
                pp_hic=converthic.pp_hic,
                pp_directory=pp_directory,
                IntType=IntType,
                genomename=genomename,
                sampleid=sampleid
        }
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 8 ----------- ###
### --------------- The End ---------------- ###
### ---------------------------------------- ###

    output {
        File? hichipjsonfile = hichipjson.out_json
        File? hicjsonfile = hicjson.out_json
        File? hicpro_zip = select_first([filterblklist.hicpro_filtered_zip, step2hicpro.hicpro_zip])
        File? hicfile = converthic.hic_file
        File? peaksbw = oneDpeaks.oneDpeaksbw
        File? peaksbdg = oneDpeaks.oneDpeaksbdg
        File? peaksbed = oneDpeaks.oneDpeaksbed
        File? peakslog = oneDpeaks.oneDpeakslog
        File? peaksreport = oneDpeaks.oneDpeaksreport
        File? fithichip_tar = twoDloops.twoDloops_tar
        File? fithichip_zip = twoDloops.twoDloops_zip
        File? anchor_bw = twoDloops.anchor_bw
        File? pp_bw = oneDpeaks.pp_bw
        File? pp_hic = converthic.pp_hic
        File? pp_peaks = twoDloops.pp_peaks
        File? pp_tbi = twoDloops.pp_tbi
        File? bedloop_1 = twoDloops.bedloop_1
        File? bedloop_2 = twoDloops.bedloop_2
        File? bedloop_3 = twoDloops.bedloop_3
        File? tbiloop_1 = twoDloops.tbiloop_1
        File? tbiloop_2 = twoDloops.tbiloop_2
        File? tbiloop_3 = twoDloops.tbiloop_3
    }
}

### ---------------------------------------- ###
### ------------ TASKS UTILIZED ------------ ###
### ---------------------------------------- ###

task fastqsnumber {
    # Get number of fastq reads to split into 30 jobs or less
    input {
        File fastqfile
        Int lengthfastqs

        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 5
    }
    
    command <<<
        maxthreads=$(printf %.0f $(echo "30 / ~{lengthfastqs}" | bc)); echo $maxthreads
        count=$(zcat ~{fastqfile} | wc -l); echo $count
        # echo "thesting"
        # ccc=$(echo "$count/(4*$maxthreads)" | bc)
        # echo $ccc
        # echo "end test"
        reads=$(printf %.0f $(echo "$count/(4*$maxthreads*1000000)" | bc)); echo $reads
        if [ $reads -le 10 ]; then 
            nreads=10000000;
        else 
            nreads=$(echo "$reads * 1000000" | bc) 
        fi

        echo "$count" | bc > totalreads.txt
        echo $nreads > nreads.txt
    >>>
    
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        Int cal_reads = read_int("nreads.txt")
        File total_reads = "totalreads.txt"
    }
}

task splitfastqs {
    # Split fastq files into chunks for running HiC-Pro
    input {
        File fastqfile
        Int nreads = 30000000
        String output_location = '_split/fastq'
        String prefix = basename(fastqfile)

        Int max_retries = 1
        Int ncpu = 1
        Int memory_gb = 5
    }

    command <<<
        mkdir -p ~{output_location}
        cd ~{output_location}
        prefix=~{prefix}
        prefix=${prefix%.f*q*}

        #making sure paired reads are R1 and R2.
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

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        Array[File] outputfiles = glob("~{output_location}/*.fastq")
    }
}

task hicpro_align {
    # Perform HiC-Pro step 1: mapping and filtering
    input {
        File fastqfile_R1
        File fastqfile_R2
        File chromsizes
        File? genomefragment
        String genomename
        String? ligationsite
        Array[File] bowtie2index
        String hicpro_out = 'HiCProOut'

        Int max_retries = 1
        Int ncpu = 4
        Int memory_gb = 10
    }

    command <<<
        # Update Config File
        mkdir temp
        pwd=$(pwd)

        config=~{basename(fastqfile_R1)}.tar
        cp /data/config-hicpro.txt config-hicpro.txt

        newfastqname=~{sub(basename(fastqfile_R1),'.fastq|.fq|.fastq.gz|.fq.gz','')}
        newfastqname_edit="${newfastqname/_R1/}"

        # Allow for Genome Fragment to not be specified
        assignfragment=~{genomefragment}
        if [[ ${#assignfragment} -ge 1 ]] && [ -f $assignfragment ]; then
            fragment_filename=${assignfragment##*/}
            GFragment=$pwd/${fragment_filename%.gz}
            if [[ "~{genomefragment}" == *"gz" ]]; then
                gunzip -c ~{genomefragment} > ${fragment_filename%.gz}
            else
                ln -s ~{genomefragment} ${fragment_filename%.gz}
            fi
            genomefragment=$GFragment; genomefragment="${genomefragment//\//\\/}"
        else
            genomefragment="";
        fi

        # Allow for Ligation Site to not be specified
        ligationsite=~{ligationsite}
        
        # Make sure Min_Cis_Dist is set when "Fragments" and "Ligation Sites" are not specified
        if [[ ${#ligationsite} -le 1 ]] && [[ ${#assignfragment} -le 1 ]] ; then 
            sed -i "s/MIN_CIS_DIST\ \=/MIN_CIS_DIST\ \=\ 1000/" config-hicpro.txt
            sed -i "s/Xgenomefragment/${genomefragment}/" config-hicpro.txt
            sed -i "s/GATCGATC/${ligationsite}/" config-hicpro.txt
        else
            sed -i "s/Xgenomefragment/${genomefragment}/" config-hicpro.txt
            sed -i "s/GATCGATC/${ligationsite}/" config-hicpro.txt
        fi
        
        index_path=~{bowtie2index[0]}; index_path=${index_path%/*}; index_path="${index_path//\//\\/}"
        chromsizes=~{chromsizes}; chromsizes="${chromsizes//\//\\/}"

        #Changes to the config file
        sed -i "s/Xbowtieindex/${index_path}/" config-hicpro.txt
        sed -i "s/Xchromsizes/${chromsizes}/" config-hicpro.txt
        sed -i "s/Xgenome/~{genomename}/" config-hicpro.txt

        # Copy Files to fastq location
        mkdir -p _split/fastq
        mv ~{fastqfile_R1} ~{fastqfile_R2} _split/fastq
        echo `date`

        getcount=${pwd%/*}
        newcount=${getcount##*/shard-}
        echo "this is sleep time $newcount"

        sleep $newcount #$(( ( RANDOM % 10 )  + 1 ))
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
                echo -e 'hilow error: Can not find valid pairs. Failed to align FASTQs using HiCPro.'
            fi
        
        else
            echo -e 'hilow error: ~{hicpro_out} folder was not created'
        fi
    >>>

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File hicpro_align = "~{basename(fastqfile_R1)}.tar"
    }
}


task hicpro_merge {
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
        String hicpro_out = 'HiCProOut'

        Int max_retries = 1
        Int ncpu = 4
    }

    Int memory_gb = ceil(length(fastqfiles_R1) * 8 ) # memory is times 3 of number of split fastqs (which = <= 90)
    Int int_fastqs = length(fastqfiles_R1) #number of split fastqs

    command <<<
        # Update Config File
        pwd=$(pwd)
        cp /data/config-hicpro.txt config-hicpro.txt

        # Allow for Genome Fragment to not be specified
        assignfragment=~{genomefragment}
        if [[ ${#assignfragment} -ge 1 ]] && [ -f $assignfragment ]; then
            fragment_filename=${assignfragment##*/}
            GFragment=$pwd/${fragment_filename%.gz}
            if [[ "~{genomefragment}" == *"gz" ]]; then
                gunzip -c ~{genomefragment} > ${fragment_filename%.gz}
            else
                ln -s ~{genomefragment} ${fragment_filename%.gz}
            fi
            genomefragment=$GFragment; genomefragment="${genomefragment//\//\\/}"
        else
            genomefragment="";
        fi

        # Allow for Ligation Site to not be specified
        ligationsite=~{ligationsite}
        
        # Make sure Min_Cis_Dist is set when "Fragments" and "Ligation Sites" are not specified
        if [[ ${#ligationsite} -le 1 ]] && [[ ${#assignfragment} -le 1 ]] ; then 
            sed -i "s/MIN_CIS_DIST\ \=/MIN_CIS_DIST\ \=\ 1000/" config-hicpro.txt
            sed -i "s/Xgenomefragment/${genomefragment}/" config-hicpro.txt
            sed -i "s/GATCGATC/${ligationsite}/" config-hicpro.txt
        else
            sed -i "s/Xgenomefragment/${genomefragment}/" config-hicpro.txt
            sed -i "s/GATCGATC/${ligationsite}/" config-hicpro.txt
        fi
        
        index_path=~{bowtie2index[0]}; index_path=${index_path%/*}; index_path="${index_path//\//\\/}"
        chromsizes=~{chromsizes}; chromsizes="${chromsizes//\//\\/}"

        #Changes to the config file
        sed -i "s/Xbowtieindex/${index_path}/" config-hicpro.txt
        sed -i "s/Xchromsizes/${chromsizes}/" config-hicpro.txt
        sed -i "s/Xgenome/~{genomename}/" config-hicpro.txt

        # Copy Files to fastq location
        mkdir -p _split/fastq
        cp -rf ~{sep=' ' fastqfiles_R1} ~{sep=' ' fastqfiles_R2} _split/fastq
        echo `date`

        # Generating SUB steps
        /HiC-Pro_3.1.0/bin/HiC-Pro -i _split -o ~{hicpro_out} -c config-hicpro.txt -p
        file_array=()
        mkdir -p AAA ; mv ~{hicpro_out} AAA
        echo `date`
        for TAR in ~{sep=' ' hicpro_align}; do
            tar -xpf $TAR
            TAR=${TAR##*/}
            file_array+=("${TAR%.tar}/~{hicpro_out}")
        done

        # Organize HiCPro File
        cp -R ${file_array[@]} AAA/~{hicpro_out} ./
        #rm -rf ${file_array} AAA/~{hicpro_out}
        cd ~{hicpro_out}
        bash HiCPro_step2_.sh

        # Compress Merge Results
        echo `date`
        cd $pwd
        tar -cpf ~{basename(hicpro_out)}.tar ~{basename(hicpro_out)}
        zip -9qr ~{basename(hicpro_out)}.zip ~{basename(hicpro_out)}
        echo `date`
        #rm -rf $pwd/~{basename(hicpro_out)}


    >>>

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File hicpro_merge = "~{basename(hicpro_out)}.tar"
        File hicpro_zip = "~{basename(hicpro_out)}.zip"
        Int count_fastqs = int_fastqs
    }
}

task oneDpeaks {
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
        String sampleid_m = if defined(sampleid) then sampleid + '.' + genomename + '_' else ""
        String hicpro_out = 'HiCProOut'
        String outdir = 'PeakOut'+ FDR

        Int max_retries = 1
        Int ncpu = 4
    }

    Int memory_gb = count_fastqs

    command <<<
        pwd=$(pwd)
        mkdir -p ~{outdir} ~{pp_directory}
        tar -xpf ~{hicpro_result}

        # Allow for Genome Fragment to not be specified
        assignfragment=~{genomefragment}
        if [ -f $assignfragment ]; then
            fragment_filename=${assignfragment##*/}
            GFragment=$pwd/${fragment_filename%.gz}
            if [[ "~{genomefragment}" == *"gz" ]]; then
                gunzip -c ~{genomefragment} > ${fragment_filename%.gz}
            else
                ln -s ~{genomefragment} ${fragment_filename%.gz}
            fi
            genomefragment=$GFragment;
        else
            genomefragment="";
        fi

        echo this is genomefragment $genomefragment
        peak_call -i ~{hicpro_out}/hic_results/data/fastq -o ~{outdir} -r $genomefragment -f ~{FDR} -a ~{chromsizes}

        cd ~{outdir}
        LC_COLLATE=C sort -k1,1 -k2,2n fastqbedgraph.bdg > fastqbedgraph.sorted.bdg
        bedGraphToBigWig fastqbedgraph.sorted.bdg ~{chromsizes} fastqbedgraph.bw

        #proteinpaint directory
        cp fastqbedgraph.bw $pwd/~{pp_directory}/~{sampleid_m}fastqbedgraph.bw
        
    >>>

    runtime {
        memory: ceil(memory_gb * 2.5) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File oneDpeaksbw = "~{outdir}/fastqbedgraph.bw"
        File oneDpeaksbdg = "~{outdir}/fastqbedgraph.bdg"
        File oneDpeaksbed = "~{outdir}/fastqpeaks.bed"
        File oneDpeakslog = "~{outdir}/fastqlog.log"
        File oneDpeaksreport = "~{outdir}/fastqreport.pdf"
        File pp_bw = "~{pp_directory}/~{sampleid_m}fastqbedgraph.bw"
    }
}

task filterblklist {
    # Filter based on BlackList regions
    input {
        Int count_fastqs
        File hicpro_result
        File blacklist
        File chromsizes
        Int slop=50
        String hicpro_out = 'HiCProOut'

        Int max_retries = 1
        Int ncpu = 1
    }

    Int memory_gb = ceil((count_fastqs * 3) + 10)

    command <<<
        pwd=$(pwd)
        tar -xpf ~{hicpro_result}

        all_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.allValidPairs
        orig_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.AVP.Orig
        mv $all_valid_pairs $orig_valid_pairs

        cd $pwd/~{hicpro_out}/hic_results/data/fastq

        if [[ "~{blacklist}" == *"gz" ]]; then
            gunzip -c ~{blacklist} > ~{sub(basename(blacklist),'.gz','')}
        else
           ln -s ~{blacklist} ~{sub(basename(blacklist),'.gz','')}
        fi

        #left bed
        awk -F "\t" -v OFS="\t" '{print $2,$3,$3,$1}' $orig_valid_pairs|slopBed -i stdin -b ~{slop} -g ~{chromsizes}|intersectBed -a stdin -b ~{sub(basename(blacklist),'.gz','')} -u > left.bed

        #right bed
        awk -F "\t" -v OFS="\t" '{print $5,$6,$6,$1}' $orig_valid_pairs|slopBed -i stdin -b ~{slop} -g ~{chromsizes}|intersectBed -a stdin -b ~{sub(basename(blacklist),'.gz','')} -u >right.bed
        

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
        tar -cpf ~{hicpro_out}.tar ~{hicpro_out}
        zip -9qr ~{hicpro_out}.zip ~{hicpro_out}
        #rm -rf ~{hicpro_out}
    >>>

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File hicpro_filtered = "~{hicpro_out}.tar"
        File hicpro_filtered_zip = "~{hicpro_out}.zip"
    }
}

task twoDloops {
    # Call 2D loops from  HiC-Pro result using FitHiChIP
    input {
        Int count_fastqs
        File hicpro_result
        String genomename
        Float FDR=0.4
        Int IntType=3
        Int LowDistThr=5000
        Int UppDistThr=4000000
        File chromsizes
        File loopBed
        String input_loopBed = basename(loopBed)
        String bw_loopBed = basename(loopBed,'.bed')
        String hicpro_out = 'HiCProOut'
        File? loopsAnchor
        Float? loopThreshold_1=0.01
        Float? loopThreshold_2=0.05
        Float? loopThreshold_3=0.1
        Array[Float] loopThresholds = select_all([loopThreshold_1,loopThreshold_2,loopThreshold_3])

        String pp_directory
        String? sampleid
        String sampleid_m = if defined(sampleid) then select_first([sampleid, genomename])+ '.' else ""

        #String outdir = 'PeakOut'+ FDR
        String loopName = if defined(loopsAnchor) then IntType + '.' + FDR + '.' + basename(loopBed,'.bed') + '.' + LowDistThr + '.' + UppDistThr else IntType + '.' + FDR + '.' + LowDistThr + '.' + UppDistThr
        String loopOut = if defined(loopsAnchor) then 'LoopOut_FitHiChIP.' + IntType + '.' + FDR + '.' + basename(loopBed,'.bed') + '.' + LowDistThr + '.' + UppDistThr else 'LoopOut_FitHiChIP.' + IntType + '.' + FDR + '.' + LowDistThr + '.' + UppDistThr

        Int max_retries = 1
        Int ncpu = 4
    }

    Int memory_gb = count_fastqs

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
        if [[ $Flag == "bedpe" ]];then
            Output=Anchor/~{sampleid_m}$InputFile
            awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' ~{loopBed} | sort -k1,1 -k2,2n > $Output
        elif [[ $Flag == "bed" ]];then
            Output=Anchor/~{sampleid_m}$InputFile
            sort -k1,1 -k2,2n ~{loopBed} |cut -f 1-3  > $Output
        elif [[ $Flag == "bw" ]];then
            Output=Anchor/~{sampleid_m}$InputFile
            cp ~{loopBed} $Output
        fi

        bgzip -f $Output
        tabix -f -p bed $Output.gz

        #copy to proteinpaint directory
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

        #run FitHiChIP
        /FitHiChIP/FitHiChIP_HiCPro.sh -C $config

        if [ ~{IntType} -eq 3 ];then
            LoopDir=$pwd/~{loopOut}/FitHiChIP_Peak2ALL_b5000_L~{LowDistThr}_U~{UppDistThr}/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr
        elif [ ~{IntType} -eq 1 ];then
            LoopDir=$pwd/~{loopOut}/FitHiChIP_Peak2Peak_b5000_L~{LowDistThr}_U~{UppDistThr}/Coverage_Bias/FitHiC_BiasCorr
        else
            echo "Error!! IntType must be 1 (Peak2Peak) or 3 (Peak2All)"
            exit;
        fi

        cd $LoopDir
        LoopThr_list=(~{sep=' ' loopThresholds}) && count=~{length(loopThresholds)}

        LoopAll=$LoopDir/FitHiChIP.interactions_FitHiC.bed

        echo -e "Get significant loops on $count FDR threholds"
        
        for Qvalue in ${LoopThr_list[@]}; do
            Output=$LoopDir/$sampleid"Qvalue."$Qvalue.~{loopName}.bed
            Output2=$LoopDir/$sampleid"Qvalue."$Qvalue.~{loopName}.bedpe
            awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7}' $LoopAll | sort -k1,1 -k2,2n > $Output2
            awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' $LoopAll | sort -k1,1 -k2,2n > $Output
            bgzip -f $Output
            tabix -f -p bed $Output.gz

            #proteinpaint directory
            cp $Output.gz{,.tbi} $pwd/~{pp_directory}
            
            echo -e "processed Qvalue $Qvalue and stored results at $Output"
        done

        cd $pwd
        tar -cpf ~{loopOut}.tar ~{loopOut}
        zip -9qr ~{loopOut}.zip ~{loopOut}
        #rm -rf $pwd/~{loopOut}

    >>>

    runtime {
        memory: ceil(memory_gb * 3) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File twoDloops_tar = "~{loopOut}.tar"
        File twoDloops_zip = "~{loopOut}.zip"
        File pp_peaks="~{pp_directory}/~{sampleid_m}~{input_loopBed}.gz"
        File pp_tbi="~{pp_directory}/~{sampleid_m}~{input_loopBed}.gz.tbi"
        File? anchor_bw = "~{pp_directory}/{bw_loopBed}.bw"
        File? bedloop_1="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_1}.~{loopName}.bed.gz"
        File? bedloop_2="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_2}.~{loopName}.bed.gz"
        File? bedloop_3="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_3}.~{loopName}.bed.gz"
        File? tbiloop_1="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_1}.~{loopName}.bed.gz.tbi"
        File? tbiloop_2="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_2}.~{loopName}.bed.gz.tbi"
        File? tbiloop_3="~{pp_directory}/~{sampleid_m}Qvalue.~{loopThreshold_3}.~{loopName}.bed.gz.tbi"
    }
}

task converthic {
    # Convert HiC-Pro result to .hic file
    input {
        Int count_fastqs
        File hicpro_result
        File chromsizes
        String hicpro_out = 'HiCProOut'
        String hic_out = 'HicFile'
        String genomename
        String pp_directory
        String? sampleid
        String sampleid_m = if defined(sampleid) then sampleid + '.' + genomename + '_' else ""

        Int max_retries = 1
        Int ncpu = 1
    }

    Int memory_gb = ceil(count_fastqs * 5) ###wait for Cas9 to be completed. It was 15

    command <<<
        pwd=$(pwd)
        mkdir -p ~{hic_out} ~{pp_directory}
        tar -xpf ~{hicpro_result}
        all_valid_pairs=$pwd/~{hicpro_out}/hic_results/data/fastq/fastq.allValidPairs

        cd ~{hic_out}
        /HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i $all_valid_pairs -g ~{chromsizes} -j /usr/local/bin/juicer_tools_1.22.01.jar

        #proteinpaint directory
        cp fastq.allValidPairs.hic $pwd/~{pp_directory}/~{sampleid_m}fastq.allValidPairs.hic
    >>>

    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File hic_file = "~{hic_out}/fastq.allValidPairs.hic"
        File pp_hic = "~{pp_directory}/~{sampleid_m}fastq.allValidPairs.hic"
    }
}

task bowtie2_index {
    # Generate bowtie2 indexes
    input {
        File reference
        String genomename

        Int memory_gb = 20
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        if [[ "~{reference}" == *"gz" ]]; then
            gunzip -c ~{reference} > ~{sub(basename(reference),'.gz','')}
        else
           ln -s ~{reference} ~{sub(basename(reference),'.gz','')}
        fi

        bowtie2-build --threads ~{ncpu} ~{sub(basename(reference),'.gz','')} ~{genomename}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        Array[File] bowtie_indexes = glob("*.bt2")
    }
}

task chrfaidx {
    # Generate genome chromosome sizes in tab delimited format
    input {
	    File reference

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        if [[ "~{reference}" == *"gz" ]]; then
            gunzip -c ~{reference} > ~{sub(basename(reference),'.gz','')}
        else
           ln -s ~{reference} ~{sub(basename(reference),'.gz','')}
        fi

        samtools faidx ~{sub(basename(reference),'.gz','')} -o ~{sub(basename(reference),'.gz','')}.fai
        cut -f1,2 ~{sub(basename(reference),'.gz','')}.fai > ~{sub(basename(reference),'.gz','')}.tab
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: ncpu
    }
    output {
        File faidx_file = "~{sub(basename(reference),'.gz','')}.fai"
        File chromsizes = "~{sub(basename(reference),'.gz','')}.tab"
    }
}

task raise_exception {
    # Raise Exception
    input {
        String msg
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        exit 2
    }
    runtime {
        memory: "5 GB"
        maxRetries: 0
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: 1
    }
    output {
        String error_msg = '${msg}'
    }
}

task createjson {
    input {
        String genomename
        String? sampleid
        String sampleid_m = if defined(sampleid) then select_first([sampleid,genomename]) + '.' + genomename + '_' else ""
        String pp_directory
        Int IntType=3
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
        
        # if there is peak data
        pp_bw=~{pp_bw}
        pp_peaks=~{pp_peaks}
        if [ -f ~{pp_peaks} ]; then
            cat >> ~{sampleid_m}proteinpaint.json <<EOF
            {
                "type":"bigwig",
                "name":"~{sampleid + '.'}Coverage",
                "file":"~{pp_directory}/${pp_bw##*/}",
                "scale":{
                    "auto": 1
                },
                "height":50
            },
            {
                "type":"bedj",
                "name":"~{sampleid + '.'}Anchor",
                "file":"~{pp_directory}/${pp_peaks##*/}",
                "stackheight":14,
                "stackspace":1
            },
        EOF
        
        fi
        
        #looping files
        for loop in ~{sep=' ' bedloops};
        do
            newloop=${loop##*Qvalue.}
            Qvalue=(${newloop//.~{IntType}/ })
            cat >> ~{sampleid_m}proteinpaint.json <<EOF
            {
                "type":"hicstraw",
                "name":"~{sampleid + '.'}$Qvalue",
                "mode_arc":true,
                "mode_hm":false,
                "bedfile":"~{pp_directory}/${loop##*/}"
            },
        EOF
        done

        # hic file
        cat >> ~{sampleid_m}proteinpaint.json <<EOF
            {
                "type":"hicstraw",
                "name":"~{sampleid + '.'}hic",
                "mode_arc":false,
                "mode_hm":true,
                "file":"~{pp_directory}/~{basename(pp_hic)}"
            }
        ]
        EOF


    >>>
    runtime {
        memory: "5 GB"
        maxRetries: 0
        docker: 'ghcr.io/stjude/abralab/hilow:v1.0'
        cpu: 1
    }
    output {
        File out_json = "~{pp_directory}/~{sampleid_m}proteinpaint.json"
    }
}
