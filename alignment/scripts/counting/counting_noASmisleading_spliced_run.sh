#!/bin/bash

samplename=$1
unique_only=$2
rowbarcode=$3
time=$4
script_path=$5
wholegene_path=$6
intron_path=$7

############### Define file with reads to keep

for output_path in $wholegene_path $intron_path	   
do

    echo '1) Identifying reads to keep...'

    file_notAS=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged.bam'
    file1=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged_genome1.bam'
    file2=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged_genome2.bam'

    # define whole gene reads as the reference
    w_notAS=$wholegene_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged.bam'
    w1=$wholegene_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged_genome1.bam'
    w2=$wholegene_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_exontagged_genome2.bam'

    temp=/scratch/local2/temp_$$.txt
    temp_b6=/scratch/local2/temp_$$_b6.txt
    temp_cast=/scratch/local2/temp_$$_cast.txt

    # B6
    samtools view $w1 | grep "GE:Z:" > $temp
    awk '{printf "%s", $0 "\tB6\n"}' $temp > $temp_b6
    rm $temp

    # Cast
    samtools view $w2 | grep "GE:Z:" > $temp
    awk '{printf "%s", $0 "\tCast\n"}' $temp > $temp_cast

    # Combine files
    cat $temp_b6 $temp_cast > $output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_reads.txt'
    rm $temp_b6
    rm $temp_cast


    ################# Subset BAM files to selected reads

    echo '2) Subset BAM files to selected reads, and keep UA reads...'

    outpath=$output_path'MergedBAM/DemplxBAM/'$time'/'
    reads=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_reads.txt'
    outfile_keep=$output_path'MergedBAM/DemplxBAM/'$time'/'$samplename'_keepreads.txt'

    $script_path'counting/counting_noASmisleading_run.R' $outpath $samplename $reads $outfile_keep

    toremove=$( wc -l < $outfile_keep )

    sleep .5s

    for file in $file_notAS $file1 $file2
    do
	echo Processing file $file...

	outfile=$( echo $file | sed 's#exontagged#keep#g' )
	id=$( echo $file | sed 's#.*/##g' | sed 's#.bam##g' | sed 's#_exontagged##g' )

	if [ $toremove == 0 ]
	then
	    echo 'Nothing to remove...'
	    cp $file $outfile
	else
	    echo 'Remove misleading reads...'
	   java -jar $script_path'picard-2.7.1/picard.jar' FilterSamReads \
		 INPUT=$file \
		 FILTER=excludeReadList \
		 READ_LIST_FILE=$outfile_keep \
		 OUTPUT=$outfile \
		 WRITE_READS_FILES=true \
		 VERBOSITY=INFO \
		 QUIET=false \
		 VALIDATION_STRINGENCY=STRICT \
		 COMPRESSION_LEVEL=5 \
		 MAX_RECORDS_IN_RAM=500000 \
		 CREATE_INDEX=false \
		 CREATE_MD5_FILE=false \
		 GA4GH_CLIENT_SECRETS=client_secrets.json
	fi
	

	# subset to UA reads only

	outfile_unique=$( echo $outfile | sed 's#keep#keep_unique#g' )

	if [ $unique_only == 'TRUE' ]
	then
	    tempfile=/scratch/local2/temp_$$.sam
	    tempfile_second=/scratch/local2/temp_$$_2.sam

	    unique=$( echo $outfile | sed 's#.bam##g' )

	    samtools view -H $outfile > $tempfile
	    samtools view $outfile >> $tempfile
	    samtools view $tempfile -S -h -q 255 > $tempfile_second
	    samtools view -bS $tempfile_second > $unique'_unique.bam'
	    samtools sort $unique'_unique.bam' -o $unique'_unique_sorted.bam'
	    samtools index $unique'_unique_sorted.bam'

	    rm $tempfile
	    rm $tempfile_second

	    input=$unique'_unique_sorted.bam'
	fi


	# not UMI counts

	out_dge=$output_path'DGE/'$time'/'$id'_unique_restrnotUMI_dge.txt'
	summary_dge=$output_path'DGE/'$time'/'$id'_unique_restrnotUMI_summary.txt'

	java -jar $script_path'Drop-seq_tools-1.12/jar/dropseq.jar' DigitalExpression \
	     I=$input \
	     O=$out_dge \
	     SUMMARY=$summary_dge \
	     OUTPUT_READS_INSTEAD=true \
	     CELL_BC_FILE=$rowbarcode \
	     EDIT_DISTANCE=0 \
	     READ_MQ=0 \
	     MIN_BC_READ_THRESHOLD=0 \
	     USE_STRAND_INFO=true \
	     RARE_UMI_FILTER_THRESHOLD=0.0


	# Count unique UMIs

	out_dge=$output_path'DGE/'$time'/'$id'_unique_restrUMI_dge.txt'
	summary_dge=$output_path'DGE/'$time'/'$id'_unique_restrUMI_summary.txt'

	java -jar $script_path'Drop-seq_tools-1.12/jar/dropseq.jar' DigitalExpression \
	     I=$input \
	     O=$out_dge \
	     SUMMARY=$summary_dge \
	     OUTPUT_READS_INSTEAD=false \
	     CELL_BC_FILE=$rowbarcode \
	     EDIT_DISTANCE=0 \
	     READ_MQ=0 \
	     MIN_BC_READ_THRESHOLD=0 \
	     USE_STRAND_INFO=true \
	     RARE_UMI_FILTER_THRESHOLD=0.0

    done
    
done

