#!/bin/bash

declare -a ARRAY

script_path=$1
input_path=$2
output_path=$input_path
rowbarcode=$3

#Read file names
samples=$4

exec 10<&0
exec < $4
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of cells to be processed:" ${#ARRAY[@]}
nsamples=${#ARRAY[@]}

let count=0

outpath_std=$output_path'stdout_stderr/'
mkdir -p $outpath_std

while [ $count -lt $nsamples ]; do

    FILE=${ARRAY[$count]}
    fastqfile=$input_path$FILE'_R2_001.fastq.gz'

    $script_path'demultiplex/demultiplex_run.sh' $script_path $input_path $output_path $rowbarcode $FILE $fastqfile

    sleep 5s
    
    ((count++))

done

