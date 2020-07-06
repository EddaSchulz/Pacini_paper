#!/bin/bash

input=$1
output=$2
rowbarcode=$3
column=$4

exec 10<&0
exec < $3

declare -a SAMPLE
let count=0

while read LINE; do
    SAMPLE[$count]=$LINE
    echo ${SAMPLE[$count]}
    ((count++))
done

echo The number of Cells is ${#SAMPLE[@]}
nsamples=${#SAMPLE[@]}

let count=0
let row=1

tempsam=/scratch/local2/temp_$$.sam

if [[ $column =~ ^-?[0-9]+$ ]]
 while [ $count -lt $nsamples ]; do
	echo Demultiplexing Cell with barcode ${SAMPLE[$count]} ..
	echo $row over $nsamples

	pattern='XC:Z:'${SAMPLE[$count]}
	echo Demultiplexing...
	samtools view -H $input > $tempsam
	samtools view $input | grep -E $pattern >> $tempsam
	samtools view -bS $tempsam > $output'/Col_'$column'_Row_'$row'.bam'
	samtools sort $output'/Col_'$column'_Row_'$row'.bam' -o $output'/Col_'$column'_Row_'$row'_sorted.bam'
	samtools index $output'/Col_'$column'_Row_'$row'_sorted.bam'
	rm $tempsam
	((count++))
	((row++))
    done
else
    echo Demultiplexing Tube Control...
    samtools view -H $input > $tempsam
    samtools view $input >> $tempsam
    samtools view -bS $tempsam > $output'/Col_'$column'.bam'
    samtools sort $output'/Col_'$column'.bam' -o $output'/Col_'$column'_sorted.bam'
    samtools index $output'/Col_'$column'_sorted.bam'
    rm $tempsam
fi

