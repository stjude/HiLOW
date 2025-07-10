InputBED=$1
Flag=$2
Dir=$3
Name=$4

name_temp=$(basename $InputBED)

InputDir=$(dirname $InputBED)

: ${Flag:=bed}
: ${Dir:=$InputDir}

mkdir -p $Dir
if [[ $Flag == "bedpe" ]];then
    name=${name_temp%.bedpe}
    : ${Name:=$name}
    Output=$Dir/$Name.ppt.bedpe
elif [[ $Flag == "bed" ]]; then
    name=${name_temp%.bed}
    : ${Name:=$name}
    Output=$Dir/$Name.ppt.bed
elif [[ $Flag == "bw" ]]; then
    name=${name_temp%.bw}
    : ${Name:=$name}
    Output=$InputBED
fi

echo $InputBED, $Dir, $Name, $Flag, $Output

module load tabix

ProteinPaintPath="/rgs01/resgen/legacy/gb_customTracks/tp/abrahgrp/jlu"
JSONPath=${ProteinPaintPath#/rgs01/resgen/legacy/gb_customTracks/tp/}

if [[ $Flag != "bw" ]]; then
    if [[ $Flag == "bedpe" ]];then
        awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7"\n"$4,$5,$6,$1,$2,$3,$7}' $InputBED|sort -k1,1 -k2,2n > $Output
    elif [[ $Flag == "bed" ]];then
        sort -k1,1 -k2,2n $InputBED|cut -f 1-3  > $Output
    fi

    bgzip -f $Output
    tabix -f -p bed $Output.gz

    echo -e "processed bedfile  $InputBED  and stored results at $Output.gz"
    cp $Output.gz $Output.gz.tbi $ProteinPaintPath

fi

JSONFile=$Dir/$Name.json

if [[ $Flag == "bedpe" ]];then
    cat > $JSONFile <<EOF
    {
        "type": "hicstraw",
        "name": "$Name",
        "mode_arc": true,
        "mode_hm": false,
        "bedfile": "$JSONPath/$Name.ppt.bedpe.gz",
        "enzyme": "MboI"
    }
EOF

elif [[ $Flag == "bed" ]];then
    cat > $JSONFile <<EOF
    {
        "type":"bedj",
        "name":"$Name",
        "file":"$JSONPath/$Name.ppt.bed.gz",
        "stackheight":14,
        "stackspace":1
    }
EOF

elif [[ $Flag == "bw" ]];then
    cp $InputBED $ProteinPaintPath
    cat > $JSONFile <<EOF
    {
        "type":"bigwig",
        "name":"$Name",
        "file":"$JSONPath/$Name.bw",
        "scale": {
             "auto": 1
         },
         "height": 50
     }
EOF
fi

