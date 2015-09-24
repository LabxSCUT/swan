#!/bin/bash
if [ -z $1 ]; then
  echo "usage: plot_bed.sh output.dir{$1} bed.file{$2} data.dir/scan.prefix{$3} data.dir/trunk.prefix {$4} chr_format{$5}"
	echo "example: plot_bed.sh example example/example.bed example/example.lib1 example/example.lib1"
	#echo "example: ./plot_bed.sh set1_lCd set1_lCd/somatic.synthetic.challenge.set1.lCd.bed set1/synthetic.challenge.set1.tumor.v2.rg1.{swan.txt.gz,trunk?.RData} set1/synthetic.challenge.set1.tumor.v2.RData"
  #echo 'for file in `ls set1/*tumor*.lCd.del.bed`; do cat $file>>set1_lCd/synthetic.challenge.set1.tumor.v2.lCd.del.bed; done;'
  #echo 'for file in `ls set1/*normal*.lCd.del.bed`; do cat $file>>set1_lCd/synthetic.challenge.set1.tumor.v2.lCd.del.bed; done;'
  #echo "bedtools intersect -v -a set1_lCd/synthetic.challenge.set1.tumor.v2.lCd.del.bed -b set1_lCd/synthetic.challenge.set1.normal.v2.lCd.del.bed >set1_lCd/somatic.synthetic.challenge.set1.lCd.bed"
  exit
fi

flank=10000
chr_format=$5
if [ -z $chr_format ]; then chr_format=""; else chr_format="chr"; fi
if [ ! -e $1 ]; then mkdir $1; fi;

opre=`basename $2 .bed`
bpre=`basename $3`
#echo "grep -v '^#' $1/$2"
#for line in `grep -v '^#' $1/$2 | cat`; do
#$lines=`sed -e '/^#/d' $1/$2`

cat $2 | while read -r line; do
  #echo $line
  [[ $line = \#* ]] && continue #skipping comment lines
  chr=`echo $line | awk {'print $1'}`
  #echo $chr
  start=`echo $line | awk {'print $2'}`
  #echo $start
  end=`echo $line | awk {'print $3'}`
  #echo $end
  type=`echo $line | grep -oP '[[:space:]]DUP[[:space:]]|[[:space:]]DEL[[:space:]]|[[:space:]]INS[[:space:]]|[[:space:]]INV[[:space:]]|[[:space:]]CNV[[:space:]]|[[:space:]]BND[[:space:]]' | grep -oP 'DUP|DEL|INS|INV|CNV|BND'`
  #echo $type
  if [ -z $type ]; then
    type="SVTYPE"
  fi
  flank_start=`expr $start - $flank`
  flank_end=`expr $end + $flank`
  echo "swan_dbg.R -c $chr -u $flank_start -v $flank_end -t $start:$end -o $1/$opre.$bpre.$type.$chr_format$chr.$flank_start-$flank_end -r $4.$chr_format$chr -e $3.$chr_format$chr $3.$chr_format$chr.swan.txt.gz"
done
  
