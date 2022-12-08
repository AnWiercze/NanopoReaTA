#!/usr/bin/env bash


###### Function to extract read length from fastq files 
echo "--------------------------------------------------"
echo "#################### READ LENGTHS ################"
echo "--------------------------------------------------"

fastq_file="$1" # $input_dir/*/*/*/*.fastq.gz
sample="$2" # ERRxyz
output_dir="$3"

mkdir -p $output_dir
if [ ! -f $output_dir/"$sample"_read_lengths_pass.txt ]
then
    echo "Length" > $output_dir/"$sample"_read_lengths_pass.txt
fi
#if [ ! -f $output_dir/"$sample"_read_lengths_fail.txt ]; then
#    echo "Length" > $output_dir/"$sample"_read_lengths_fail.txt
#fi

### Only passed reads will be processed
#if [[ "$fastq_file" == *"pass"* ]]; then
if [[ "$fastq_file" == *".gz" ]]
then
    echo "HERE" 
    zcat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_pass.txt
else
    echo "HERE2"
    cat $fastq_file | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_pass.txt
fi
#fi
#if [[ "$fastq_file" == *"fail"* ]]; then
#    if [[ "$fastq_file" == *".gz" ]]; then
#        zcat $j | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_fail.txt
#    else
#        cat $j | awk 'NR%4==2' | awk '{ print length }' >> $output_dir/"$sample"_read_lengths_fail.txt
#    fi
#fi

