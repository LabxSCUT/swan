#!/usr/bin/env bash
echo "usage: run_sub.sh sample_dir config_dir pbs_letter"
echo "now processing $1/$2/*.$3.pbs"

for sub in $(ls -d $1/$2); do
  echo $sub
  for pbsfile in $(ls -f $sub/*.$3.pbs); do
    bamfile=${pbsfile%.$3.pbs}.bam
    echo $bamfile
      if  [ ! -e $bamfile ]; then
        echo "-------------------"
        echo "processing $pbsfile"
        echo "-------------------"
        ssa.py $pbsfile
        #. $pbsfile
        sleep 1
      fi
  done
done
echo "all submitted"
