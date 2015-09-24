#!/usr/bin/env bash
if [ -z $1 ]; then 
	echo "usage: vcf2bed.sh your.vcf"
	echo "bed format"
	echo "chr \t start \t end \t type \t size"
	exit
fi
file=$1
cat $file | grep -v "#" | grep -o ";END=[0-9]\+" | grep -o "[0-9]\+" >$file.col3.txt #end point, must to ;END= otherwise conflict with CIEND
cat $file | grep -v "#" | cut -f 1 >$file.col1.txt   #chr
cat $file | grep -v "#" | cut -f 2 >$file.col2.txt   #start point
cat $file | grep -v "#" | cut -f 5 | grep -o "[^><]*" >$file.col4.txt
paste $file.col3.txt $file.col2.txt | awk '{a=$2; print $1-a}' >$file.col5.txt 
cat $file | grep -v "#" | cut -f 4- >$file.col6.txt #
paste $file.col1.txt $file.col2.txt $file.col3.txt $file.col4.txt $file.col5.txt $file.col6.txt >${file%.vcf}.bed
rm -f $file.col*.txt
