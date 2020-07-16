#!/bin/bash

declare -a ARRAY

script_path=$1
input_path=$2
sample_names=$3
unique_only=$4
reference=$5
annotation=$6
rowbarcode=$7
output_path=$8

#Read file names
unaligned_path="${10}"

exec 10<&0
exec < $9
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of time points to be processed:" ${#ARRAY[@]}
ntime=${#ARRAY[@]}

let count=0

mkdir -p $input_path'MergedBAM'
cd $input_path'MergedBAM/'

while [ $count -lt $ntime ]; do

    echo Processing time point ${ARRAY[$count]} ..
    echo $count over $ntime

    time=${ARRAY[$count]}
    unaligned=$unaligned_path$time'/'
    aligned=$input_path'STAR/'$time'/'
    includeMM='false'
    paired='false'

    $script_path'counting/counting_run.sh' $script_path $input_path $time $reference $unaligned $aligned $includeMM $paired $unique_only $annotation $rowbarcode $output_path $sample_names

    sleep 60s
    
    ((count++))

done

