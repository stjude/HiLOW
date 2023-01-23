#!/bin/bash
## Test WDL script using samplefiles provided.

# CROMWELL on St. Jude hpc user directory
sjhpc_cromwell="/home/madetunj/Modupe/publiclink/cromwell-52.jar"
lsf_config="/home/madetunj/Modupe/publiclink/lsf.conf"
script="../hilow.wdl"

input="./inputs.json"
option="./options.json"

# STD OUT and ERR files
logout="wdlhilow_out"
logerr="wdlhilow_err"

#script syntax
wdlscript="java -Dconfig.file=$lsf_config \
    -jar $sjhpc_cromwell \
    run $script \
    --inputs $input \
    --options $option"
    bsub -P watcher -q compbio \
      -R "rusage[mem=10000]" \
      -J wdlhilow \
      -o $logout \
      -e $logerr \
      -N $wdlscript
