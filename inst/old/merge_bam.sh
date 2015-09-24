#!/usr/bin/env bash
echo "#usage: merge_bam.sh bam_dir output_dir output.bam"
echo "#usage: merge *.bam under bam_dir into output_dir/output.bam"
if [ -z $3 ]; then exit; fi

echo "#input bam_files: find $1/ -name \"*.bam\""
bam_files=($(find $1/ -name "*.bam"))
echo "#first bam_file is ${bam_files[0]}"
header_sam=$2/$3.sam
echo "#header bam_file is $2/$3.sam" 
echo "samtools view -H ${bam_files[0]} >$header_sam"
echo "#all bamfiles are ${bam_files[*]}"
bam_input=${bam_files[*]}
echo "samtools merge -h $header_sam $2/$3 $bam_input"
