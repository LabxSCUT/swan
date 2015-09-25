#!/bin/bash
  
  echo "testing SWAN binaries in $SWAN_BIN"
  if [ -z $1 ]; then
		echo "usage:"
		echo "  paired.sh one|two|all"
	fi
  #add -a -q to SWAN commands to increase verbosity and debug output

  #$SWAN_BIN/swan_stat
  ### this is an example to demonstrate SWAN use for single sample ###  

  #0. OBSOLETE to run this example please download and install svengine
  #echo rm -rf master.zip charade-svengine-*
  #rm -rf master.zip charade-svengine-*
  #echo wget https://bitbucket.org/charade/svengine/get/master.zip
  ##wget https://bitbucket.org/charade/svengine/get/master.zip
  #echo wget https://bitbucket.org/charade/svengine/get/devel.zip -O master.zip #for developer use only
  #wget https://bitbucket.org/charade/svengine/get/devel.zip -O master.zip #for developer use only
  #echo unzip master.zip 
  #unzip master.zip 
  #echo 'export PDIR=`ls -d charade-svengine-*/`' 
  #export PDIR=`ls -d charade-svengine-*/` 
  #echo cd $PDIR
  #cd $PDIR
  #echo easy_install . 
  #easy_install . 
  #cd ..
  
  #1. OBSOLETE use svengine to generate a small size WGS data with embedded structure variant
  #these input files such as .meta, .fna, .gap.bed, .par are provided
  #echo 'fasforge -o example/example -m example/example.meta example/example.bed example/example.par example/example.fna'
  #fasforge -f 10000 -o example/example -m example/example.meta example/example.gap.bed example/example.par example/example.fna
  #echo 'for file in `ls example/example.lib*.bam`; do mv $file ${file%.bam}.us.bam; samtools sort ${file%.bam}.us.bam ${file%.bam}; rm -f ${file%.bam}.us.bam; samtools index $file; done;'
  #for file in `ls example/example.lib*.bam`; do mv $file ${file%.bam}.us.bam; samtools sort ${file%.bam}.us.bam ${file%.bam}; rm -f ${file%.bam}.us.bam; samtools index $file; done;
  #merge mlib files
  #samtools merge example/example.mlib.bam example/example.lib1.bam example/example.lib2.bam
  #samtools index example/example.mlib.bam
  #for file in `ls example/example.lib*.bam`; do $SWAN_BIN/swan_stat -c a $file; done
  
  #2. run swan_stat to get hold of the sequencing parameters and statistics
  #    see details on https://bitbucket.org/charade/swan/wiki/Example#rst-header-a-simple-example-run-for-one-sample
  if [ $1 != "two" ]; then
  echo "testing paired sample single read group pipeline..."
  #2.1 single sample single read group:
  #    run swan_stat
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

  if [ $1 != "one" ]; then
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

  echo "done testing $SWAN_BIN"
  
  #6. OBSOLETE check the recovery of embedded events
  #cat example/example.bed | wc -l
  #bedtools intersect -a -u  example/example.bed -b example/example.lib1.conf.bed | wc -l
  
  #7. OBSOLETE plot the spiked-in and found region to visualize event characteristics
  #../inst/plot_bed example example/example.bed example/example.lib1 example/example.lib1
  #../inst/plot_bed example example/example.lib1.conf.bed example/example.lib1 example/example.lib1
  #../inst/swan_plot -f bed -o example/example example/example.bed example/example.lib1.bam #all true region
  #../inst/swan_plot -f bed -o example/example.lib1.conf example/example.lib1.conf.bed example/example.lib1.bam #all true region
  
  #8. Any questions or comments, please send to lixia@stanford.edu
  #   Also check SWAN's bitbucket site for latest updates
  #   Thanks you for using SWAN
