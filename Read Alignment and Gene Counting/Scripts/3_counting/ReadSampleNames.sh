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

    output=$unaligned$samplename
    name='Counting_SeqHT_'$time

        ### Merge and Tag alignments with DGE --> notAS
    final_file=$output_path'DGE/'$time'/'$samplename'_unique_UMI_dge.txt'

    if [ ! -f $final_file ]

    then
	$script_path'Merge_Tag_Alignments.sh' $reference $unaligned_file $aligned_file $output $includeMM $paired $samplename $unique_only $annotation $input_path $rowbarcode $output_path $time $script_path

	sleep 1s
    else
	   echo 'notAS file exists already...'
    fi

    
    ### Merge and Tag alignments with DGE --> AS

    aligned_file1=$aligned$samplename'_DemFilt/Aligned.sortedByCoord.out.genome1.bam'
    aligned_file2=$aligned$samplename'_DemFilt/Aligned.sortedByCoord.out.genome2.bam'

    name1='Counting_SeqHT_'$time'_BL6'
    name2='Counting_SeqHT_'$time'_CastEiJ'
    name3='Counting_nomisleading'

    g1='genome1'
    g2='genome2'


    # B6    
    out=$output_path'DGE/'$time'/'$samplename'_'$g1'_unique_UMI_dge.txt'
    if [ ! -e $out ]
    then
	echo Processing to produce B6 reads
	$script_path'Merge_Tag_Alignments_AS.sh' $reference $unaligned_file $aligned_file1 $output $includeMM $paired $samplename $unique_only $annotation $input_path $rowbarcode $output_path $time $g1 $script_path	

	sleep 1s
    else
	echo 'B6 file is already present...'
    fi

    # Cast
    out=$output_path'DGE/'$time'/'$samplename'_'$g2'_unique_UMI_dge.txt'
    if [ ! -e $out ]
    then
	echo Processing to produce Cast reads
	$script_path'Merge_Tag_Alignments_AS.sh' $reference $unaligned_file $aligned_file2 $output $includeMM $paired $samplename $unique_only $annotation $input_path $rowbarcode $output_path $time $g2 $script_path

	sleep 1s
    else
	echo 'Cast file is already present...'
    fi
        

    ((count++))

done

