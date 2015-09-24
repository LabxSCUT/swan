#!/usr/bin/env bash

usage="usage: split_lib.sh input.bam"
info="info: split_lib.sh splits input.bam by libraries defined in read group"
dep="dep: split_lib.sh depends on samtools"

if [ -z $1 ]; then
  echo "usage error!"
  echo $usage
  echo $info
  echo $dep
  exit
fi
bamfile=$1
#bamopre=`basename $bamfile .bam`
bamopre=${bamfile%.bam}
liblist=$bamopre.lib
if [ -e $liblist ]; then
  echo "overwrite old lib list!"
  rm -f $liblist
fi

lb=`samtools view -H $bamfile | grep @RG | grep -o 'LB:\S*' | awk -F ':' '{print $2}' | uniq | sort | tr '\n' ' '`
#note some libriries may have multiple read groups, we combine them
echo "libs:" $lb
let i=1
for lib in `echo $lb`; do
  echo samtools view -bh -l $lib -o $bamopre.lib$i.bam $bamfile
  samtools view -bh -l $lib -o $bamopre.lib$i.bam $bamfile
  echo "$bamopre.lib$i.bam $lib" >>$liblist 
  let i=i+1
done
