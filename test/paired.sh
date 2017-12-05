  #!/bin/bash

  #0 prompt the correct usage 
  if [ -z $2 ]; then
    echo "usage:"
    echo "  paired.sh $SWAN_BIN one|two|all"
    exit
  fi

	SWAN_BIN=$1
	NLIB=$2

  #1.  checking options
  echo "testing SWAN binaries in $SWAN_BIN in the case of $NLIB libraries"

  #2   see details on https://bitbucket.org/charade/swan/wiki/Example#rst-header-a-simple-example-run-for-one-sample
  #2.1 single sample single read group:
  #    run swan_stat
  if [ $NLIB != "two" ]; then
  echo "testing paired sample single read group pipeline..."
  $SWAN_BIN/swan_stat -c a example/normal.lib1.bam
  $SWAN_BIN/swan_stat -c a example/tumor.lib1.bam
  #    run swan_scan
  $SWAN_BIN/swan_scan -c a -n example/example.gap.bed example/example.fna example/normal.lib1.bam 
  $SWAN_BIN/swan_scan -c a -n example/example.gap.bed example/example.fna example/tumor.lib1.bam 
  #    run sclip_scan
  $SWAN_BIN/sclip_scan -c a -o example/normal_tumor.lib1 -n example/example.gap.bed example/example.fna example/normal.lib1.bam:example/tumor.lib1.bam 
  #    run swan_join
  $SWAN_BIN/swan_join -c a -t example/normal.lib1.stat:example/tumor.lib1.stat -i example/normal.lib1:example/tumor.lib1 -j example/normal.lib1:example/tumor.lib1 -m example/normal.lib1:example/tumor.lib1 -l example/normal_tumor.lib1.sclip.RData example/example.fna example/normal.lib1.bam:example/tumor.lib1.bam
  fi

  if [ $NLIB != "one" ]; then
  echo "testing paired sample multiple read group pipeline..."
  #2.2 single sample multiple read group pipeline..."
  #    mlib bam file is merged lib1 and lib2 bam files
  #    run swan_stat on unmerged bamfiles
  $SWAN_BIN/swan_stat -c a -m example/normal.mlib example/normal.lib1.bam,example/normal.lib2.bam
  $SWAN_BIN/swan_stat -c a -m example/tumor.mlib example/tumor.lib1.bam,example/tumor.lib2.bam
  #    run swan_scan on unmerged bamfiles get SWAN scan tracks in .swan.txt.gz file
  $SWAN_BIN/swan_scan -o example/normal.mlib -c a -n example/example.gap.bed example/example.fna example/normal.lib1.bam,example/normal.lib2.bam 
  $SWAN_BIN/swan_scan -o example/tumor.mlib -c a -n example/example.gap.bed example/example.fna example/tumor.lib1.bam,example/tumor.lib2.bam 
  #    run sclip_scan on merged bamfiles to get SWAN soft-sclipping cluster remapping in .RData file
  $SWAN_BIN/sclip_scan -c a -o example/normal_tumor.mlib -n example/example.gap.bed example/example.fna example/normal.mlib.bam:example/tumor.mlib.bam
  #    run swan_join on unmerged bamfiles to join the calls
  $SWAN_BIN/swan_join -c a -t example/normal.mlib.stat:example/tumor.mlib.stat -i example/normal.mlib:example/tumor.mlib -j example/normal.mlib:example/tumor.mlib -m example/normal.mlib:example/tumor.mlib -l example/normal_tumor.mlib.sclip.RData example/example.fna example/normal.lib1.bam,example/normal.lib2.bam:example/tumor.lib1.bam,example/tumor.lib2.bam
  fi

  #NOTE:
  #add -q to above commands to increase verbosity 
  #add -a to above commands to get additional debug output
  #the script is simply for testing functionality of SWAN
  #there is no scientific value within the analysis
  #Any questions or comments, please send to lixia@stanford.edu
  #Also check SWAN's bitbucket site for latest updates
  #Thanks you for using SWAN

  echo "done testing $SWAN_BIN"
