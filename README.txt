HiChIP/HiC Looping Organization Workflow (Hi-LOW): one-click identification of high-confidence HiChIP loops

Introduction
========================

This is all-in-one pipleline performing the following analysis:

1. Splitting large fastq files into chunks of 10M reads

2. Mapping and filtering reads in parallel (HiC-Pro step1)

3. Merging mapped reads and normalizing contact matrix (HiC-Pro step2)

4. Call 1D ChIP peaks from HiC-Pro result  (macs2)

5. Call 2D loops from HiC-Pro result (FitHiChIP)

6. Converting allValid Pairs to .hic files (jucer_tools)

7. Copy .hic files to appropriate location (/rgs01/resgen/legacy/gb_customTracks/tp) and generate JSON file for visualization by ProteinPaint

8. Generate bigwig file for 1D peaks and related JSON file for visualization by ProteinPaint

9. Generate BEDPE file for 2D loops  and related JSON file for visualization by ProteinPaint

Usage
========================
usage : bash HiC-Pro_parallel_auto.sh -i InputDir -o OutDir [-c ConfigFile] [-f flag] [-n CPU_step1] [-g Gsize] [-j JuicerJar] [-p ProteinPaintPath] [-R refGenome] [-r Restrition] [-l LowDistThr] [-u UppDistThr] [-I IntType] [-t LoopThr] [-h]
Use option -h for more information
    Options:
        -i  InputDir            InputDir as requrest for HiC-Pro input, the input directory is in the format of WD/SampleID, the fastq file are orgnised as WD/SampleID/fastq/*.{fq,fastq,fastq.gz}
        -d  ScriptDir           directory for all HiChIP pipeline scripts, default to be "/rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/HiChIP_pipeline/script". External user should change this to their own installation directory
        -o  OutDir              user specified directory for all output of pipeline
        -c  ConfigFile          Configuration file for HiChIP/Hi-C processing
        -f  Flag                HiC data flag, the default data type is HiChIP
        -n  CPU_step1           Number of CPUs used per mapping job in HiC-Pro step 1. (total CPUs needed = Number of CPU per job * Number of parallel jobs)
        -g  Gsize               genome size file used during HiC-Pro processing
        -j  JuicerJar           path to juicer_tools.jar file, used to convert HiC-Pro ressult to .hic file
        -p  ProteinPaintPath    path to store files for visulization using ProteinPaint
        -H  HiCProBaseDir       main directory of HiCPro output
        -R  refGenome           Reference genome string used for displaying track in UCSC. Default is 'hs19' for human chromosome. For mouse, specify 'mm10' or 'mm9'
        -r  RestrictionF        Absolute path of restriction fragment file used by HiC-Pro and HiChIP-Peaks.
        -l  LowDistThr          Lower distance threshold of interaction between two intervals (CIS). Default: 5000 (indicates 5Kb). Interactions below this distance threshold will not be considered for statistical significance.
        -u  UppDistThr          Upper distance threshold of interaction between two intervals (CIS). Default: 4000000 (indicates 4Mb). Interactions above  this distance threshold will not be considered for statistical significance.
        -I  IntType             Type of interaction (foreground) reported by FitHiChIP.1) peak to peak; 2) peak 2 non peak; 3) peak to all (default); 4) all to all; 5) everything from 1 to 4
        -t  LoopThr             A list of FDR threshold used to calling loops using FitHiChIP, provide in the format of "0.01 0.05 0.1".
        -M  FDR                 FDR used to call 1D peaks by HiChIP peaks
        -L  ReadLength          Length of reads for the HiC-pro generated reads. Default 75 *** (for potential furthur usage)
        -a  loopAnchor          Optional file containing pre-defined set of loop anchors, only the first 3 cols are relevent (Chr,Start,End). 1D peak calling step will skipped if provided.
        -s  Slop                Increase the Valid Pair entry by the same number base pairs in each direction,Default=50, meaning 100bp for both paired read.
        -b  BlackList           Optional file containing a list of region need to be filtered out on the level of allValidPairs. With each of the paired reads being set to be 2*Slop (-s), pairs would be filtered out if either read overlap 1bp with any black listed regions.

