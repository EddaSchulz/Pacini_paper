#!/bin/bash

script_path=$1
input_path=$2
time=$3
reference=$4
unaligned=$5
aligned=$6
includeMM=$7
paired=$8
unique_only=$9
annotation="${10}"
rowbarcode="${11}"
output_path="${12}"

#Read file names
samplenames="${13}"
wholegene_path="${14}"
intron_path="${15}"

exec 10<&0
exec < $samplenames
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    #echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of cells to be processed:" ${#ARRAY[@]}
ncells=${#ARRAY[@]}

let count=0

while [ $count -lt $ncells ]; do
    
    samplename=${ARRAY[$count]}
    echo Processing sample $samplename from time point $time
    unaligned_file=$unaligned$samplename'_sorted.bam'
aligned_file=$aligned$samplename'_DemFilt/Aligned.sortedByCoord.out.bam'

    output=$unaligned$samplename

    name='RemoveMisleading_'$time

    # Remove misleading reads and repeat counting
    echo Remove misleading

$script_path'Merge_Tag_RemoveMisleading_AS_intronexon_removemisleading.sh' $samplename $unique_only $rowbarcode $time $script_path $wholegene_path $intron_path

    ((count++))
done
