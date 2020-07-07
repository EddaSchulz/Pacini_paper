#!/bin/bash

inpath=$1
outpath=$2
softwpath=$3
mismatches=$4
indicespath=$5
multimap=$6
snp_file=$7
time=$8

echo The fastq files are stored at $inpath ...
echo The output text files are then saved at $outpath ...
echo The .sh files can be found at the directory $softwpath ...
echo The maximum number of mismatches is $mismatches ...
echo The indices are stored at $indicespath ...
echo The maximum number of sequences to which a read can align is specified to be $multimap ...
echo The SNP file is stored at $snp_file ...

cd $inpath

ncells=$( ls *.fastq | wc -l )

echo There are $ncells different cells to analyze

ls $inpath*.fastq > $outpath'dir_fastq.txt'

$softwpath'STAR_cluster_AS.sh' $inpath $outpath $softwpath $mismatches $indicespath $multimap $outpath'dir_fastq.txt' $snp_file $time


