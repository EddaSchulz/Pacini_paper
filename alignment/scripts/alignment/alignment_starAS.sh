#!/bin/bash

FASTQ=$1
inpath=$2
outpath=$3
mismatches=$4
indicespath=$5
multimap=$6
snp_file=$7
softwpath=$8

echo Start processing file $FASTQ

### Define path to STAR aligner
STAR=$softwpath'STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR'

### Define the number of threads
threads=5

### Alignment
FILE=$(echo $FASTQ | sed 's#.*/##g' | sed 's#.fastq##g')

echo Make directory $FILE
mkdir -p $outpath$FILE;

$STAR --runThreadN $threads \
	--genomeDir $indicespath \
	--genomeLoad NoSharedMemory \
      --readFilesIn $FASTQ \
	--outFileNamePrefix $outpath$FILE'/' \
      --quantMode TranscriptomeSAM GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
      --outWigType bedGraph \
	--outWigStrand Unstranded \
      --outFilterMismatchNmax $mismatches \
	--outFilterMultimapNmax $multimap  \
      --limitBAMsortRAM 20000000000 \
      --alignEndsType EndToEnd \
      --outSAMattributes NH HI NM MD

### Index BAM files for visualization on IGV
samtools index $outpath$FILE'/Aligned.sortedByCoord.out.bam'

### SNPsplit

echo Splitting alignment into parental genomes...

alignment=$outpath$FILE'/Aligned.sortedByCoord.out.bam'

$softwpath'/SNPsplit_v0.3.2/SNPsplit --snp_file $snp_file \
	     --conflicting \
	     $alignment

echo Analysis on cell $FILE completed.

