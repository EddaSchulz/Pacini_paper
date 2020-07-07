#!/bin/bash

declare -a ARRAY

inpath=$1
outpath=$2
softwpath=$3
mismatches=$4
indicespath=$5
multimap=$6
snp_file=$8
time=$9

exec 10<&0
exec < $7
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of cells to be aligned:" ${#ARRAY[@]}
nsamples=${#ARRAY[@]}

outpath_std=$outpath'stdout_stderr/'
mkdir -p $outpath_std

let count=0

while [ $count -lt $nsamples ]; do

    FILE=$(echo ${ARRAY[$count]} | sed 's#.*/##g' | sed 's#.fastq##g')
    FILE1=$inpath$FILE'.fastq'

    echo This is FILE: $FILE
    echo This is FILE1: $FILE1
    echo inp: $inpath
    echo out: $outpath
    echo soft: $softwpath
    echo mismatches: $mismatches
    echo indices: $indicespath
    echo multimap: $multimap
    echo This is multimap: $multimap
    

    name='scRNASeq_SeqHT_AS_'$time

    if [ ! -e $output$FILE ]
    then
	$softwpath'runSTAR_AS.sh' $FILE1 $inpath $outpath $mismatches $indicespath $multimap $snp_file $softwpath
    else
	echo file exists already!
    fi
    
    ((count++))

done

