#!/usr/bin/env bash

usage="usage: split_vcf.sh input.vcf chr_format"
info="info: split_vcf.sh splits input.vcf by chromosomes"
dep="dep: split_vcf.sh depends on samtools"

declare -a chr_1_22=($(seq 1 22))
#echo ${chr_1_22[@]}
declare -a chr_X_Y=("X" "Y")
#echo ${chr_X_Y[@]}
declare -a chr_noM=( "${chr_1_22[@]}" "${chr_X_Y[@]}" )
echo "#chr_noM=${chr_noM[@]}"

if [ -z $1 ]; then
  echo "#usage error!"
  echo "#$usage"
  echo "#$info"
  echo "#$dep"
  exit
fi
vcffile=$1
chrfmt=$2

for chr in ${chr_noM[@]}; do
  vcfopre=${vcffile%.vcf}.$chr.vcf
  echo "cat $vcffile | grep \"^$chrfmt$chr\s\" | grep \"<DUP>\|<DEL>\|<INS>\|<INV>\" >$vcfopre"
  echo "vcf2bed.sh $vcfopre"
done
