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
    output=$unaligned$samplename    # to extend as _merged and _tagged

    name='Counting_nomisleading_'$time

    aligned_file1=$aligned$samplename'_DemFilt/Aligned.sortedByCoord.out.genome1.bam'
    aligned_file2=$aligned$samplename'_DemFilt/Aligned.sortedByCoord.out.genome2.bam'
    g1='genome1'
    g2='genome2'

    # Remove misleading reads and repeat counting

    out=$output_path'DGE/'$time'/'$samplename'_unique_restrUMI_dge.txt'
    out1=$output_path'DGE/'$time'/'$samplename'_genome1_unique_restrUMI_dge.txt'
    out2=$output_path'DGE/'$time'/'$samplename'_genome2_unique_restrUMI_dge.txt'

$script_path'counting/counting_noASmisleading_run.sh' $reference $unaligned_file $aligned_file2 $output $includeMM $paired $samplename $unique_only $annotation $input_path $rowbarcode $output_path $time $g2 $script_path 

    ((count++))

done

