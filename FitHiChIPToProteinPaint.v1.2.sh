SampleID=$1
Input=$2
LoopThr_list=(${3//,/ })

InputDir=$(dirname $Input)
temp=${InputDir/*\/LoopOut_FitHiChIP./}
name=${temp/\/*/}
InputFile=$(basename $Input)

#echo $Input, $InputFile, $SampleID, ${LoopThr_list[@]}
module load tabix

for Qvalue in ${LoopThr_list[@]}; do
    Output=$InputDir/$SampleID.Qvalue.$Qvalue.$name.bed
    Output2=$InputDir/$SampleID.Qvalue.$Qvalue.$name.bedpe
    awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7}' $Input|sort -k1,1 -k2,2n > $Output2
    awk -F '\t' -v OFS='\t' -v Q=$Qvalue '((NR>1) && ($NF<Q)) {print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' $Input|sort -k1,1 -k2,2n > $Output
    bgzip -f $Output
    tabix -f -p bed $Output.gz
    echo -e "processed Qvalue $Qvalue and stored results at $Output"

done
