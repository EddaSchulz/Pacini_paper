#!/bin/bash

echo Reading files...

reference=$1
unaligned_file=$2
aligned_file=$3
output=$4
includeMM=$5
paired=$6
samplename=$7
unique_only=$8
annotation=$9
input_path="${10}"
rowbarcode="${11}"
output_path="${12}"
time="${13}"
script_path="${14}"

############### Merge file

echo '1) Merge files...'

#output_merge=$output'_merged.bam'
mkdir -p $output_path'MergedBAM/DemplxBAM/'$time'/'
output_merge=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_merged.bam'

java -jar $script_path'picard-2.7.1/picard.jar' MergeBamAlignment \
	 REFERENCE_SEQUENCE=$reference \
	 UNMAPPED_BAM=$unaligned_file \
	 ALIGNED_BAM=$aligned_file \
	 OUTPUT=$output_merge \
	 INCLUDE_SECONDARY_ALIGNMENTS=$includeMM \
	 PAIRED_RUN=$paired

################ Tag file

echo '2) Tag files...'

tag='GE'
#output_tag=$output'_exontagged.bam'
output_tag=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged.bam'

java -jar $script_path'Drop-seq_tools-1.12/jar/dropseq.jar' TagReadWithGeneExon \
	 I=$output_merge \
	 O=$output_tag \
	 ANNOTATIONS_FILE=$annotation \
	 TAG=$tag \
	 ALLOW_MULTI_GENE_READS='true' \
	 USE_STRAND_INFO='false'


################ Count features

echo '3) Count features and define DGE list...'

input=$output_tag
#output=$input_path'DGE/'$time'/'
output=$output_path'DGE/'$time'/'
mkdir -p $output

echo The output directory is $output

output_file=$output$samplename

if [ $unique_only == 'TRUE' ]
then

    echo Unique Only option...
    
    tempfile=/scratch/local2/temp_$$.sam
    tempfile_second=/scratch/local2/temp_$$_2.sam
    

    unique_bam=$( echo $input | sed 's#_exontagged.bam##g' )
    sample=$( echo $input | sed 's#_exontagged.bam##g' | sed 's#.*/##g' )

    samtools view -H $input > $tempfile
    samtools view $input >> $tempfile
    samtools view $tempfile -S -h -q 255 > $tempfile_second
    samtools view -bS $tempfile_second > $unique_bam'_exontagged_unique.bam'
    samtools sort $unique_bam'_exontagged_unique.bam' -o $unique_bam'_exontagged_unique_sorted.bam'
    samtools index $unique_bam'_exontagged_unique_sorted.bam'

    rm $tempfile
    rm $tempfile_second

    input=$unique_bam'_exontagged_unique_sorted.bam'
    output_file=$output_file'_unique'
fi

### not UMI counts

echo '3.1) notUMI - Count features...'

output_temp=$output_file'_notUMI_dge.txt'
summary_temp=$output_file'_notUMI_summary.txt'

java -jar $script_path'Drop-seq_tools-1.12/jar/dropseq.jar' DigitalExpression \
	 I=$input \
	 O=$output_temp \
	 SUMMARY=$summary_temp \
	 OUTPUT_READS_INSTEAD=true \
	 CELL_BC_FILE=$rowbarcode \
	 EDIT_DISTANCE=0 \
	 READ_MQ=0 \
	 MIN_BC_READ_THRESHOLD=0 \
	 USE_STRAND_INFO=true \
	 RARE_UMI_FILTER_THRESHOLD=0.0

### UMI counts

echo '3.2) UMI - Count features...'

output_temp=$output_file'_UMI_dge.txt'
summary_temp=$output_file'_UMI_summary.txt'

java -jar $script_path'Drop-seq_tools-1.12/jar/dropseq.jar' DigitalExpression \
	 I=$input \
	 O=$output_temp \
	 SUMMARY=$summary_temp \
	 OUTPUT_READS_INSTEAD=false \
	 CELL_BC_FILE=$rowbarcode \
	 EDIT_DISTANCE=0 \
	 READ_MQ=0 \
	 MIN_BC_READ_THRESHOLD=0 \
	 USE_STRAND_INFO=true \
	 RARE_UMI_FILTER_THRESHOLD=0.0



