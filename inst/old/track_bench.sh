#!/bin/bash
#echo "# run_bench.sh truth.bed{$1} attemp.bed{$2}"
hit=`bedtools intersect -a $1 -b $2 -u | wc -l`
call=`cat $2 | wc -l`
target=`cat $1 | wc -l`
#echo "# bedtools intersect -a $1 -b $2 -u | wc -l"
#echo "# cat $2 | wc -l"
#echo "# cat $1 | wc -l"
declare -i hit
declare -i call
declare -i target
prec=`echo "scale=3;$hit/$call" | bc`
recall=`echo "scale=3;$hit/$target" | bc`
accuracy=`echo "scale=3;2*$prec*$recall/($prec+$recall)" | bc`
echo $hit $call $target $accuracy $prec $recall
