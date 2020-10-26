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

#Xist 5' quantification
output_path_5prime=$9
snp_file_5prime="${10}"

# summarize input files
echo The .sh files can be found at the directory $script_path ...
echo The fastq files are stored at $input_path ...
echo The alignment output is saved at $output_path ...
echo The maximum number of mismatches is $mismatches ...
echo The indices are stored at $indicespath ...
echo The maximum number of alignments are $multimap ...
echo The fastq file names are stored at $samples ...
echo The SNP file is stored at $snp_file ...
echo The alignment output restricted to Xist5prime SNPs is saved at $output_path_5prime ...
echo The Xist5prime SNP file is stored at $snp_file_5prime ...

# Read and store file names
exec 10<&0
exec < $7
let count=0

while read LINE; do
    ARRAY[$count]=$LINE
    echo ${ARRAY[$count]}
    ((count++))
done

echo "Number of samples to be processed:" ${#ARRAY[@]}
nsamples=${#ARRAY[@]}

let count=0
mkdir -p $output_path

while [ $count -lt $nsamples ]; do

   FILE=${ARRAY[$count]}
   echo $FILE
   
   $script_path'alignment/alignment_star_bulk.sh' \
   $input_path \
   $script_path \
   $indicespath \
   $mismatch \
   $multimap \
   $FILE \
   $output_path \
   $snp_file \
   $output_path_5prime \
   $snp_file_5prime
   
   sleep 5s
    
   ((count++))

done

