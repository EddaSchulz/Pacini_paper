#!/bin/bash

script_path=$1
input_path=$2
output_path=$3
rowbarcode=$3
samplename=$4
fastqfile=$5

columnid=$( echo $samplename | sed 's#.*COL##g' | sed 's#_S.*##g' )
timepoint=$( echo $samplename | sed 's#.*_SC##g' | sed 's#-COL.*##g' )



#########################################################################
### 1) Fastq to SAM

echo '1) fastq to SAM file...'

outpath=$output_path'UnalignedBAM/'
mkdir -p $outpath

READ1=$input_path$samplename'_R1_001.fastq.gz'
READ2=$input_path$samplename'_R2_001.fastq.gz'
OUTFILE=$outpath$samplename'_unBAM.bam'
SAMPLENAME=$samplename
LIBRARYNAME='Column_'$columnid

java -jar $script_path’/picard-2.7.1/picard.jar’ FastqToSam \
	 FASTQ=$READ1 \
	 FASTQ2=$READ2 \
	 OUTPUT=$OUTFILE \
	 SAMPLE_NAME=$SAMPLENAME \
	 LIBRARY_NAME=$LIBRARYNAME

sleep 5s


############################################################################ 2) Demultiplex Unaligned BAM files

echo '2.1) Cell - Demultiplex Unaligned BAM files...'

outpath=$output_path'CellandMolecule_Demultiplex/'
mkdir -p $outpath

input=$OUTFILE
OUTFILE=$outpath$samplename'_Cell_DMPX.bam'
SUMMARY=$outpath$samplename'_Cell_DMPX_summary.txt'

java -jar $script_path’/Drop-seq_tools-1.12/jar/dropseq.jar’ TagBamWithReadSequenceExtended \
	 INPUT=$input \
	 OUTPUT=$OUTFILE \
	 SUMMARY=$SUMMARY \
	 BASE_RANGE=1-6 \
	 BASE_QUALITY=10 \
	 BARCODED_READ=1 \
	 DISCARD_READ=False \
	 TAG_NAME=XC \
	 NUM_BASES_BELOW_QUALITY=1 

sleep 5s

##### Molecular demultiplexed

echo '2.2) Molecular - Demultiplex Unaligned BAM files...'

input=$OUTFILE
OUTFILE=$outpath$samplename'_CellMolecular_DMPX.bam'
SUMMARY=$outpath$samplename'_CellMolecular_DMPX_summary.txt'

java -jar $script_path’/Drop-seq_tools-1.12/jar/dropseq.jar’ TagBamWithReadSequenceExtended \
	 INPUT=$input \
	 OUTPUT=$OUTFILE \
	 SUMMARY=$SUMMARY \
	 BASE_RANGE=7-11 \
	 BASE_QUALITY=10 \
	 BARCODED_READ=1 \
	 DISCARD_READ=True \
	 TAG_NAME=XM \
	 NUM_BASES_BELOW_QUALITY=1 

sleep 5s


############################################################################ 3) Filter out tagged reads

echo '3.1) Filter out tagged reads...'

outpath=$output_path'CellandMolecule_Demultiplex/'

input=$OUTFILE
OUTFILE=$outpath$samplename'_Filtered_DMPX.bam'

java -jar $script_path’/Drop-seq_tools-1.12/jar/dropseq.jar’ FilterBAM \
	 TAG_REJECT=XQ \
	 INPUT=$input \
	 OUTPUT=$OUTFILE

sleep 5s


#################


echo '3.2) Filter out tagged BAM...'

outpath=$output_path'MergedBAM/DemplxBAM/'$timepoint
mkdir -p $outpath

input=$OUTFILE

$script_path'demultiplex/demultiplex_TaggedBAM.sh' $input $outpath $rowbarcode $columnid

sleep 5s


#################

echo '3.2) Filter out tagged BAM...'

combinations=$input_path'trimmed_sample_'$columnid'.txt'

outpath=$output_path'MergedBAM/DemplxBAM/'$timepoint
mkdir -p $outpath

input=$OUTFILE

$script_path' demultiplex/demultiplex_samTofastq.sh' $output_path $combinations $timepoint
