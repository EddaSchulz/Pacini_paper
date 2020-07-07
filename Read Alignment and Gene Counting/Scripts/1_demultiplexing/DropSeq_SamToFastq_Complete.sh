#!/bin/bash

inpath=$1
samplenames=$2
timepoint=$3

exec 10<&0
exec < $2

declare -a SAMPLE
let count=0

while read LINE; do
    SAMPLE[$count]=$LINE
    echo ${SAMPLE[$count]}
    ((count++))
done

echo The number of samples is ${#SAMPLE[@]}
nsamples=${#SAMPLE[@]}

### Sam to Fastq with Picard

out=$inpath'Demultiplexed_Filtered_Fastq/'$timepoint'/'
mkdir -p $out

let count=0

while [ $count -lt $nsamples ]; do

    echo Transforming back to fastq sample ${SAMPLE[$count]} ..
    echo $count over $nsamples	
input=$inpath'MergedBAM/DemplxBAM/'$timepoint'/'${SAMPLE[$count]}'_sorted.bam'
    output=$out${SAMPLE[$count]}'_DemFilt.fastq'
    java -jar $script_path’/picard-2.7.1/picard.jar’ SamToFastq \
	     INPUT=$input \
	     FASTQ=$output 
  	     
    ((count++))
done

