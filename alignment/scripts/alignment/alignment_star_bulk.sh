#!/bin/bash

input_path=$1
script_path=$2
indicespath=$3
mismatch=$4
multimap=$5
sample=$6
output_path=$7
snp_file=$8
output_path_5prime=$9
snp_file_5prime="${10}"

# Define read 1 and 2
R1=$input_path$sample'_R1_001.fastq.gz'
R2=$input_path$sample'_R2_001.fastq.gz'

# Define output folder for each sample
id=$(echo $sample | sed 's#^.*\-A[1-6]-##')
out=$output_path$id"/"
mkdir -p $out

# Define path to STAR aligner and number of threads
STAR=$script_path'STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR'
threads=5

# Alignment

$STAR --runThreadN $threads \
	--genomeDir $indicespath \
	--genomeLoad NoSharedMemory \
      --readFilesIn $R1 $R2 \
	--outFileNamePrefix $out \
      --quantMode TranscriptomeSAM GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
      --outWigType bedGraph \
	--outWigStrand Unstranded \
      --outFilterMismatchNmax $mismatch \
	--outFilterMultimapNmax $multimap  \
      --limitBAMsortRAM 20000000000 \
      --alignEndsType EndToEnd \
      --outSAMattributes NH HI NM MD \
      --readFilesCommand zcat

# Index BAM files for visualization on IGV
samtools index $out'Aligned.sortedByCoord.out.bam'

# SNPsplit: all SNPs
echo Splitting alignment into parental genomes...
alignment=$out'Aligned.sortedByCoord.out.bam'
$softwpath'/SNPsplit_v0.3.2/SNPsplit' --snp_file $snp_file \
	     --conflicting \
	     $alignment
	     
# SNPsplit: Xist 5' SNPs
out=$output_path_5prime$id"/"
mkdir -p $out
cp $alignment $out
$softwpath'/SNPsplit_v0.3.2/SNPsplit' --snp_file $snp_file_5prime \
	     --conflicting \
	     $alignment
	     
echo Analysis on sample $id is completed.