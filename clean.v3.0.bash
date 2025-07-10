
###Clean HiChIP pipeline result

# Could runt script in two ways
# 1. Provide the absolute path of the Hi-Low output file directory as argument, e.g.
#   bash /rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/HiChIP_pipeline/script/clean.sh /rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/Jie/archive/HiChIP/Ying/Nalm6-AM-CTCF-HiChIP_2Covaris-G_112020.pipeline_output
# 2. run the script within the HiChIP output directory and without argument, e.g.,
#    cd /rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/Jie/archive/HiChIP/Ying/Nalm6-AM-CTCF-HiChIP_2Covaris-G_112020.pipeline_output
#    bash /rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/HiChIP_pipeline/script/clean.sh

InputDir=${1:-$PWD}

echo target dirctory  to clean: $InputDir

## Initial size
echo "#####"initial outputsize: $(du -hs $InputDir|cut -f 1)


# remove splited fastq input files

target=$(find $InputDir/*_split -maxdepth 0 -type d -exec echo {} \;)

if [[ -d $target ]];then
    s=$(du -hs $target|cut -f 1)
    echo remove file $target $s
    rm -fr $target
else
    echo "directory for split fastq not found!!"
fi

# remove bowtie results of HiC-Pro

target2=$(find $InputDir/HiCPro_out/bowtie_results -maxdepth 0 -type d -exec echo {} \;)

if [[ -d $target2 ]]; then
    s=$(du -hs $target2|cut -f 1)
    echo remove file $target2 $s
    rm -fr $target2
else
    echo "directory for bowtie results not found!!"
fi

# remove bam and validPaires output for splited fastq from HiC-Pro
if target4="$(find $InputDir/HiCPro_out/hic_results/data/fastq/*.{bam,validPairs,allValidPairs,Orig} -type f -exec echo {} \;)";then
    s=$(du -hs $target4|cut -f 1)
    echo remove file $target4 $s
    rm $target4
else
    echo "bam and/or validPaires output for splited fastq not found!!"
fi

: <<'END'
# remove matrix files from HiC-Pro

target6=$(find $InputDir/HiCPro_out/hic_results/matrix/fastq -maxdepth 0 -type d -exec echo {} \;)

if [[ -d $target6 ]]; then
    s=$(du -hs $target6|cut -f 1)
    echo remove file $target6 $s
    rm -fr $target6
else
    echo "matrix files from HiC-Pro not found!!"
fi
END

# remove bedgraph files of MACS2

if target3="$(find $InputDir/PeakOut.0.4/*.{bdg,bam} -type f -exec echo {} \;)"; then
    s=$(du -hs $target3|cut -f 1)
    echo remove file $target3 $s
    rm $target3
else
    echo "bedgraph files and/or merged bam files not found!!"
fi

# remove inputs of MACS2

#if target5="$(find $InputDir/PeakOut/MACS2_input.bed -type f -exec echo {} \;)"; then
#    s=$(du -hs $target5|cut -f 1)
#    echo remove file $target5 $s
#    rm $target5
#else
#    echo "MACS input  not found!!"
#fi

## size after remove files
echo "#####" size after cleaning: $(du -hs $InputDir|cut -f 1)
