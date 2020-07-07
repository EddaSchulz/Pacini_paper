#!/bin/bash

declare -a ARRAY

script_path=$1
input_path=$2
output_path=$3

#STAR
mismatches=$4
indicespath=$5
multimap=$6

#Read file names
samples=$7
snp_file=$8

exec 10<&0
exec < $7
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of time points to be processed:" ${#ARRAY[@]}
ntime=${#ARRAY[@]}

let count=0

STAR_output=$output_path'STAR/'
mkdir -p $STAR_output

while [ $count -lt $ntime ]; do

   TimePoint=${ARRAY[$count]}
   fastq_path=$input_path$TimePoint'/'
   out_time=$STAR_output$TimePoint'/'
   
   mkdir -p $out_time
   $script_path'CompleteLaunch_STAR_AS.sh' $fastq_path $out_time $script_path $mismatches $indicespath $multimap $snp_file $TimePoint
    
   sleep 60s
    
   ((count++))

done

