  #!/bin/bash
  
  #0 prompt the correct usage
  if [ -z $2 ]; then
    echo "usage:"
    echo "  single.sh $SWAN_BIN one|two|all"
    exit
  fi

  SWAN_BIN=$1
  NLIB=$2

  #1.  checking options
  echo "testing SWAN binaries in $SWAN_BIN in the case of $NLIB libraries"

  #2.   see details on https://bitbucket.org/charade/swan/wiki/Example#rst-header-a-simple-example-run-for-paired-sample
  if [ $NLIB != "two" ]; then
  echo "testing single sample single read group pipeline..."
  #2.1 single sample single read group:
  #    run swan_stat
  $SWAN_BIN/swan_stat -c a example/example.lib1.bam
  #    run swan_scan
  $SWAN_BIN/swan_scan -c a -n example/example.gap.bed example/example.fna example/example.lib1.bam 
  #    run sclip_scan
  $SWAN_BIN/sclip_scan -c a -o example/example.lib1 -n example/example.gap.bed example/example.fna example/example.lib1.bam 
  #    run swan_join
  $SWAN_BIN/swan_join -c a -i example/example.lib1 -j example/example.lib1 -l example/example.lib1.sclip.RData -m example/example.lib1 example/example.fna example/example.lib1.bam
  fi

  if [ $NLIB != "one" ]; then
  echo "testing single sample multiple read group pipeline..."
  #2.2 single sample multiple read group pipeline..."
  #    mlib bam file is merged lib1 and lib2 bam files
  #    run swan_stat on unmerged bamfiles
  $SWAN_BIN/swan_stat -c a -m example/example.mlib example/example.lib1.bam,example/example.lib2.bam
  #    run swan_scan on unmerged bamfiles get SWAN scan tracks in .swan.txt.gz file
  $SWAN_BIN/swan_scan -o example/example.mlib -c a -n example/example.gap.bed example/example.fna example/example.lib1.bam,example/example.lib2.bam 
  #    run sclip_scan on merged bamfiles to get SWAN soft-sclipping cluster remapping in .RData file
  $SWAN_BIN/sclip_scan --vcf -c a -o example/example.mlib -n example/example.gap.bed example/example.fna example/example.mlib.bam 
  #    run swan_join on unmerged bamfiles to join the calls
  $SWAN_BIN/swan_join -c a -t example/example.mlib.stat -i example/example.mlib -j example/example.mlib -l example/example.mlib.sclip.RData -m example/example.mlib example/example.fna example/example.lib1.bam,example/example.lib2.bam
  fi

  #NOTE:
  #add -q to above commands to increase verbosity 
  #add -a to above commands to get additional debug output
  #add -q to above commands to increase verbosity
  #add -a to above commands to get additional debug output
  #the script is simply for testing functionality of SWAN
  #there is no scientific value within the analysis
  #Any questions or comments, please send to lixia@stanford.edu
  #Also check SWAN's bitbucket site for latest updates
  #Thanks you for using SWAN

  echo "done testing $SWAN_BIN"
