#!/bin/bash
#Promoter list, format: Chr/tStart/tEnd/tRefSeqID|GeneSymbol
PromoterFile=${1:-'/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/new_mm10.clean.4kbproms.bed'}
#bam file of promoter signal measurment, e.g., H3K27ac ChIP-Seq
#BamFile=${2:-'/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie_mm10/SRR1613251.bam'}
BamDir=$2
WD=$3
HiChIPPeaks=$WD/fastqpeaks.bed
#bam file of Contorl for promoter signal measurment, e.g., H3K27ac ChIP-Seq Input
#ControlBamFile=${3:-'/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie_mm10/SRR1613249.bam'}
#/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_MED1.sorted.bam
#/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_Input.sorted.bam

module load bedtools/2.30.0 samtools/1.17

BamFile=$WD/merged.bam

cd $WD

echo -e "$(date)"
echo "################Merge HiCPro generated bams######################"
if [ ! -f $BamFile ]
then
    start_time="$(date -u +%s.%N)"
    samtools merge -@ 8 $BamFile  $BamDir/*.bam
    end_time="$(date -u +%s.%N)"
    elapsed="$(bc <<<"$end_time-$start_time")"
    echo -e "Total of $elapsed seconds elapsed for process of merging HicPro bam files\n"
else
    echo -e "Detecting merged bam files at location: $BamFile. Skipping...\n"
fi

echo -e "$(date)"
echo "################Remove PCR duplicates from merged bam file######################"
if [ ! -f ${BamFile/.bam/.rmdup.bam} ]
then
    start_time="$(date -u +%s.%N)"
    samtools rmdup -S $BamFile ${BamFile/.bam/.rmdup.bam.temp}
    samtools sort -@ 8 ${BamFile/.bam/.rmdup.bam.temp} -o ${BamFile/.bam/.rmdup.bam}
    end_time="$(date -u +%s.%N)"
    elapsed="$(bc <<<"$end_time-$start_time")"
    echo -e "Total of $elapsed seconds elapsed for process of removing PCR duplicates and sorting bam file: $BamFile\n"
else
    echo -e "Detecting duplicates removed and sorted  bam files at location: ${BamFile/.bam/.rmdup.bam}. Skipping...\n"
fi


echo -e "$(date)"
echo "################Get active promoters######################"
base1=$(basename $PromoterFile)
id1=${base1/.bed/}
if [ ! -f $id1.unique.activePromoter.bed ]
then
    start_time="$(date -u +%s.%N)"
    sort -k1,1 -k2,2n -k4,4 $PromoterFile |awk -F'\t' -v OFS="\t" '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$4 }END{ for(i in a) print i,a[i] }'|sort -k1,1 -k2,2n > $id1.unique.bed
    totalOrig=$(cat $PromoterFile|wc -l)
    total=$(cat $id1.unique.bed|wc -l)
    keepN=$(bc <<< "$total*2/3")

    #echo $total, $keepN
    intersectBed -c -a $id1.unique.bed -b ${BamFile/.bam/.rmdup.bam} > $id1.Signal
    #intersectBed -c -a $PromoterFile -b $ControlBamFile > $id3.Signal
    #paste $id2.Signal  $id3.Signal|awk -F '\t' -v OFS='\t' '{norm= ($5>=$10 ? $5-$10 : 0)} {print $1,$2,$3,$4,norm}' > $id2.Norm.Signal
    sort -k5,5nr  $id1.Signal|head -n $keepN > $id1.keep.$keepN.Signal
    sort -k1,1 -k2,2n $id1.keep.$keepN.Signal|cut -f 1-4 > $id1.unique.activePromoter.bed
    end_time="$(date -u +%s.%N)"
    elapsed="$(bc <<<"$end_time-$start_time")"
    echo -e "Total of $elapsed seconds elapsed for process of quantifying promoter HiChIP coverage: $id1.unique.bed\n"
    echo -e "$totalOrig promoters are collapsed into $total unique ones, among which $keepN are selected as active promoters based on bam file: ${BamFile/.bam/.rmdup.bam}\n"
else
    echo -e "Detecting active promoter files at location: $WD/$id1.unique.bed. Skipping...\n"
fi

echo -e "$(date)"
echo "################Generate Loop Anchor######################"
if [ ! -f $id1.HiChIPPeak.combined.loopAnchor.bed ]
then
    cat <(cut -f 1-3 $HiChIPPeaks|awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$1"_"$2"_"$3}')  <(cut -f 1-4 $id1.unique.activePromoter.bed)|sort -k1,1 -k2,2n|mergeBed -i stdin -c 4 -o collapse -delim ";" > $id1.HiChIPPeak.combined.bed
    cut -f 1-3 $id1.HiChIPPeak.combined.bed > $id1.HiChIPPeak.combined.loopAnchor.bed
    echo -e "Loop anchors for FitHiChIP generated: $WD/$id1.HiChIPPeak.combined.loopAnchor.bed"
else
    echo -e "Detecting Loop anchors for FitHiChIP: $WD/$id1.HiChIPPeak.combined.loopAnchor.bed. Skipping...\n"
fi


