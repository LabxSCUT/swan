#!/usr/bin/env bash
if [ -z $3 ]; then
  echo "#usage: swan_run.sh data_dir spX:spY stage dry chr_format"
  echo "#aim: a script to ease use of swan pipeline on human genomics"
  echo "#input: data_dir/spX.bam | data_dir/spX.bam:data_dir/spY.bam | data_dir/{spX.1.bam,spX.2.bam,...}:data_dir{spY.1.bam,spY.2.bam,...}"
  echo "#working: PWD is working directory where intermediatary results are"
  echo "#data_dir: data_dir is the directory where bam files are"
  echo "#chromosome: reference is formated as chr_format{$5}+seq_index"
  exit
fi

env_value="#$ -v LD_LIBRARY_PATH=$HOME/usr/lib:$HOME/usr/lib64\n#$ -v PATH=/usr/local/bin:/bin:/usr/bin:$HOME/scripts" #\n#$ -v BOOST_HOME=$HOME/usr"
data_dir=$1
swan_stage=$3
chr_format=""
dry_mode=""
if [ -z $4 ]; then dry_mode=$4; fi
if [ -z $5 ]; then chr_format=$5; fi
declare -a chr_1_22=($(seq 1 22))
#echo ${chr_1_22[@]}
declare -a chr_X_Y=("X" "Y")
#echo ${chr_X_Y[@]}
declare -a chr_noM=( "${chr_1_22[@]}" "${chr_X_Y[@]}" )
echo "#chr_noM=${chr_noM[@]}"
if [ -z $chr_format ]; then # formatted as X, use 1kg
  ref_file="$HOME/ws/hg/hg19/human_g1k_v37.fasta"
  gap_file="$HOME/ws/hg/hg19/human_g1k_v37.fasta.gaphc.bed"
else
  ref_file="$HOME/ws/hg/hg19/human_g1k_v37.ucsc.fasta"
  gap_file="$HOME/ws/hg/hg19/human_g1k_v37.ucsc.fasta.gaphc.bed"
fi
sp_files=($(echo $2 | tr ":" "\t")) #0-indexed array of samples
echo "#data_dir=$data_dir"
echo "#sp_files=${sp_files[*]}"
echo "#swan_stage=$swan_stage"
echo "#dry_mode=$dry_mode"
echo "#chr_format=$chr_format"
let start_time=$SECONDS
declare -a jobs
declare -a names
let cnt=0
run_cmd=""

xbam=($(echo ${sp_files[0]} | tr "," "\t"))          #components of spX
xbam1=${xbam[0]}             #first of spX
ybam=($(echo ${sp_files[1]} | tr "," "\t"))          #components of spY
ybam1=${ybam[0]}             #first of spY
echo "#xbam=${xbam[*]}"
echo "#xbam=${#xbam[@]}"
echo "#xbam1=$xbam1"
echo "#ybam=${ybam[*]}"
echo "#xbam=${#ybam[@]}"
echo "#ybam1=$ybam1"
echo "#PWD=$PWD"

if [ ${#xbam[@]} == 1 ]; then
  xbam_files=$data_dir/${xbam1}
else
  xbam_files="echo ${sp_files[0]} | awk -F , '{ for (i=1;i<=NF;i++){ printf \"$data_dir\"\"/\"\$i\",\" }}'"
  xbam_files=`eval $xbam_files`
  xbam_files=${xbam_files%","}
fi
if [ ${#ybam[@]} == 1 ]; then
  ybam_files=$data_dir/${ybam1}
else
  ybam_files="echo ${sp_files[1]} | awk -F , '{ for (i=1;i<=NF;i++){ printf \"$data_dir\"\"/\"\$i\",\" }}'"
  ybam_files=`eval $ybam_files`
  ybam_files=${ybam_files%","}
fi
echo "#xbam_files=$xbam_files"
echo "#ybam_files=$ybam_files"

if [ ${#xbam[@]} == 1 ]; then
  xprefix=${xbam1%".bam"}
else
  xprefix="echo $xbam1 | awk -F . '{ for (i=1;i<=NF-2;i++){ printf \$i\".\" } }'" #add \ to $ sign is important to preserve it from expanding in .sh
  #echo $xprefix
  xprefix=`eval $xprefix` # xprefix would be exactly runnable on cmd line and eval takes its return, easy for debug
  xprefix=${xprefix%"."}
fi
if [ ${#ybam[@]} == 1 ]; then
  yprefix=${ybam1%".bam"}
else
  yprefix="echo $ybam1 | awk -F . '{ for (i=1;i<=NF-2;i++){ printf \$i\".\" } }'"
  yprefix=$($yprefix)
  yprefix=${yprefix%"."}
fi
echo "#xprefix=$xprefix"
echo "#yprefix=$yprefix"

### all chr ###
if [ $swan_stage == "swstat" ]; then
  run_name=$xprefix
  run_cmd="swan_stat.R -o $xprefix $data_dir/$xbam1"
  pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=32gb,vmem=32gb,pmem=32gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  jobs[$cnt]=$pbs_cmd
  names[$cnt]=$run_name
  let cnt=$cnt+1
fi

if [ $swan_stage == "scscanA" ]; then
  if [ ! -e "$PWD/${xprefix}_scscanA" ]; then
    mkdir "$PWD/${xprefix}_scscanA" #create subdir named $bam_file to save RData files
  fi
  run_cmd="sclip_scan.R -a -c A -n $gap_file -o $PWD/${xprefix}_scscanA/${xprefix} $ref_file $xbam_files";
  run_name=${xprefix}.chrA
  pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  jobs[$cnt]=$pbs_cmd
  names[$cnt]=$run_name
  let cnt=$cnt+1
fi

### chr-wise ###
for chr in ${chr_noM[@]}; do

  #replaced sclip_scan2, sclip_call2, sclip_event with sclip_scan
  #if [ $swan_stage == "1sccall2" ]; then
  #  if [ ! -e "$data_dir/$xprefix-1sccall2" ]; then
  #    echo mkdir "$data_dir/$xprefix-1sccall2"
  #    mkdir "$data_dir/$xprefix-1sccall2" #create subdir for RData files
  #  fi
  #  trunk_num=30
  #  for ti in $(seq 1 $trunk_num); do #default 3 30
  #    run_cmd="sclip_call2.R -n $gap_file -o $data_dir/$xprefix-1sccall2/$xprefix -c $trunk_num:$ti $ref_file $data_dir/$xprefix-scscan2"
  #    run_name=$xprefix.$ti
  #    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -l h_vmem=30G\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=8gb,vmem=8gb,pmem=8gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #    jobs[$cnt]=$pbs_cmd
  #    names[$cnt]=$run_name
  #    let cnt=$cnt+1
  #  done
  #  break #we don't need to loop over chr here, break out
  #fi
  
  #if [ $swan_stage == "sccall" ]; then #not available for multibam yet
  #  if [ ! -e "$data_dir/$xprefix-sccall" ]; then
  #    echo mkdir "$data_dir/$xprefix-sccall"
  #    mkdir "$data_dir/$xprefix-sccall" #create subdir for RData files
  #  fi
  #  trunk_num=30
  #  for ti in $(seq 1 $trunk_num); do
  #    run_cmd="sclip_call.R -d $data_dir/$yprefix-scscan -n $gap_file -o $data_dir/$xprefix-sccall/$xprefix -c $trunk_num:$ti $ref_file $data_dir/$xprefix-scscan"
  #    run_name=$xprefix.$ti
  #    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=8gb,vmem=8gb,pmem=8gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #    jobs[$cnt]=$pbs_cmd
  #    names[$cnt]=$run_name
  #    let cnt=$cnt+1
  #  done
  #  break #we don't need loop over chr here, break out
  #fi

  #if [ $swan_stage == "sccall" ]; then #not available for multibam yet
  #  if [ ! -e "$data_dir/$xprefix-sccall" ]; then
  #    echo mkdir "$data_dir/$xprefix-sccall"
  #    mkdir "$data_dir/$xprefix-sccall" #create subdir for RData files
  #  fi
  #  trunk_num=10
  #  for ti in $(seq 1 $trunk_num); do  #default 3 and 30
  #    run_cmd="sclip_call2.R  -d $data_dir/$yprefix-scscan -n $gap_file -o $data_dir/$xprefix-sccall/$xprefix -c $trunk_num:$ti $ref_file $data_dir/$xprefix-scscan"
  #    run_name=$xprefix.$ti
  #    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=8gb,vmem=8gb,pmem=8gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #    jobs[$cnt]=$pbs_cmd
  #    names[$cnt]=$run_name
  #    let cnt=$cnt+1
  #  done
  #  break #we don't need chr loop here, break out
  #fi

  #if [ $swan_stage == "sceventTN" ]; then
  #  run_cmd="sclip_events.R -a $data_dir/$xbam1.bam -b $data_dir/$ybam1.bam -s $data_dir/$xprefix.stat -t $data_dir/$yprefix.stat -o $data_dir/$xprefix.sceventTN -p $data_dir/$xprefix.sclip.event $ref_file $gap_file $data_dir/$xprefix-$yprefix-sccallTN"
  #  run_name=$xprefix
  #  pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #  jobs[$cnt]=$pbs_cmd
  #  names[$cnt]=$run_name
  #  let cnt=$cnt+1
  #  break #we don't need chr loop here, break out
  #fi

  #if [ $swan_stage == "1scevent" ]; then
  #  run_name=$xprefix
  #  run_cmd="sclip_events.R -a $data_dir/$xbam1.bam -s $data_dir/$xprefix.stat -o $data_dir/$xprefix.sclip -p $data_dir/$xprefix.sclip $ref_file $gap_file $data_dir/$xprefix-1sccall2 "
  #  pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=4\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #  jobs[$cnt]=$pbs_cmd
  #  names[$cnt]=$run_name
  #  let cnt=$cnt+1
  #  break
  #fi

  if [ $swan_stage == "seqcbs" ]; then
    run_cmd="seqcbs_scan.R -a -q -c $chr_format$chr $ref_file $ybam_files:$xbam_files -o $data_dir/$yprefix.$chr:$data_dir/$xprefix.$chr"
    run_name=$xprefix.$chr
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=20gb,vmem=20gb,pmem=20gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  #if [ $swan_stage == "sccallTN" ]; then
  #  if [ ! -e "$data_dir/$xprefix-sccallTN" ]; then
  #    echo mkdir "$data_dir/$xprefix-sccallTN"
  #    mkdir "$data_dir/$xprefix-sccallTN" #create subdir for RData files
  #  fi
  #  trunk_num=30
  #  for ti in $(seq 1 $trunk_num); do #default 3 30
  #    run_cmd="sclip_call2.R -d $data_dir/$yprefix-scscan2N -n $gap_file -o $data_dir/$xprefix-$yprefix-sccallTN/$xprefix -c $trunk_num:$ti $ref_file $data_dir/$xprefix-scscan2T"
  #    run_name=$xprefix.$ti
  #    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -l h_vmem=30G\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l mem=8gb,vmem=8gb,pmem=8gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
  #    jobs[$cnt]=$pbs_cmd
  #    names[$cnt]=$run_name
  #    let cnt=$cnt+1
  #  done
  #  break #we don't need to loop over chr here, break out
  #fi

  if [ $swan_stage == "swjoinTN_ql" ]; then #join seqcbs (q) + sclip (l)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -l $data_dir/$xprefix.sclip.RData -k $data_dir/$xprefix.$chr.seqcbs.txt $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.ql.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.ql.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_q" ]; then #join only seqcbs (q)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -k $data_dir/$xprefix.$chr.seqcbs.txt $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.cb.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.cb.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_l" ]; then #join only sclip (l)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -c $chr_format$chr -l $data_dir/$xprefix.sclip.RData -o $data_dir/$xprefix.$chr $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.sc.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.sc.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_g" ]; then #join only big deletion (g)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -c $chr_format$chr -j $data_dir/$yprefix.$chr.bigd.txt:$data_dir/$xprefix.$chr.bigd.txt -o $data_dir/$xprefix.$chr $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.bd.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.bd.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_s" ]; then #join only disc duplication/translocation (s)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -c $chr_format$chr -m $data_dir/$yprefix.$chr.disc.txt:$data_dir/$xprefix.$chr.disc.txt -o $data_dir/$xprefix.$chr $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.di.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.di.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_C" ]; then #join only lCd (C)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -u track=lCd,method=theo,thresh=level3,sup=100,gap=100 -o $data_dir/$xprefix.$chr $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.cd.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.cd.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_D" ]; then #join only lDl+lDr (D)
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -u track=lDl+lDr,method=theo,thresh=level3,sup=100,gap=100 -o $data_dir/$xprefix.$chr $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.dx.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.dx.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_glqCD" ]; then #join all cbs,sc,swan results in matched case, letter from a-zA-Z
    run_name=$xprefix.$chr
    run_cmd="swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -j $data_dir/$yprefix.$chr.swan.bigd.txt:$data_dir/$xprefix.$chr.swan.bigd.txt -k $data_dir/$xprefix.$chr.seqcbs.txt -l $data_dir/$xprefix.sclip.RData $ref_file $data_dir/$ybam:$data_dir/$xbam "
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_dream" ]; then #special join and confirm for dream challenge
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -j $data_dir/$yprefix.$chr.bigd.txt:$data_dir/$xprefix.$chr.bigd.txt -k $data_dir/$xprefix.$chr.seqcbs.txt -l $data_dir/$xprefix.sclip.RData -m $data_dir/$yprefix.$chr.disc.txt:$data_dir/$xprefix.$chr.disc.txt -u \"\" $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.bed $data_dir/$xprefix.$chr.dr.bed; mv $data_dir/$xprefix.$chr.vcf $data_dir/$xprefix.$chr.dr.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_C" ]; then #join only lCd
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -u track=lCd,method=theo,thresh=level1,sup=100,gap=100 $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.swd.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.swd.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_D" ]; then #join only lDx
    run_name=$xprefix.$chr
    run_cmd="(swan_join.R -a -q -o $data_dir/$xprefix.$chr -c $chr_format$chr -i $data_dir/$yprefix.$chr.swan.txt.gz:$data_dir/$xprefix.$chr.swan.txt.gz -u track=lDl+lDr,method=theo,thresh=level3,sup=100,gap=100 $ref_file $data_dir/$ybam:$data_dir/$xbam; mv $data_dir/$xprefix.$chr.raw.bed $data_dir/$xprefix.$chr.swx.bed; mv $data_dir/$xprefix.$chr.raw.vcf $data_dir/$xprefix.$chr.swx.vcf)"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_gsCD" ]; then
    run_name=$xprefix.$chr
    run_cmd="swan_join.R -a -q -o $PWD/$xprefix.$chr -c $chr_format$chr -t $PWD/$yprefix.stat:$PWD/$xprefix.stat -i $PWD/$yprefix.$chr.swan.txt.gz:$PWD/$xprefix.$chr.swan.txt.gz -j $PWD/$yprefix.$chr.bigd.txt:$PWD/$xprefix.$chr.bigd.txt -m $PWD/$xprefix.$chr.disc.txt:$PWD/$yprefix.$chr.disc.txt $ref_file $data_dir/$ybam:$data_dir/$xbam "
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoinTN_sCD" ]; then
    run_name=$xprefix.$chr
    run_cmd="swan_join.R -a -q -o $PWD/$xprefix.$chr -c $chr_format$chr -t $PWD/$yprefix.stat:$PWD/$xprefix.stat -i $PWD/$yprefix.$chr.swan.txt.gz:$PWD/$xprefix.$chr.swan.txt.gz -m $PWD/$xprefix.$chr.disc.txt:$PWD/$yprefix.$chr.disc.txt $ref_file $data_dir/$ybam:$data_dir/$xbam "
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi


  if [ $swan_stage == "swjoin" ]; then # CD for swan, g for bigd, s for disc, l for clip
    run_name=$xprefix.$chr
    run_cmd="swan_join.R -a -q -o $PWD/$xprefix.$chr -c $chr_format$chr -i $PWD/$xprefix.$chr.swan.txt.gz -j $PWD/$xprefix.$chr.swan.bigd.txt -l $PWD/$xprefix.sclip.RData -m $PWD/$xprefix.$chr.disc.txt $ref_file $xbam_files"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swjoin_glsCD" ]; then # CD for swan, g for bigd, s for disc, l for clip
    run_name=$xprefix.$chr
    run_cmd="swan_join.R -a -q -o $PWD/$xprefix.$chr -c $chr_format$chr -i $PWD/$xprefix.$chr.swan.txt.gz -j $PWD/$xprefix.$chr.swan.bigd.txt -l $PWD/$xprefix.sclip.RData -m $PWD/$xprefix.$chr.disc.txt $ref_file $xbam_files"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swscan" ]; then #default
    run_name=$xprefix.$chr
    run_cmd="swan_scan.R -c $chr_format$chr -n $gap_file -o $PWD/$xprefix.$chr $ref_file $xbam_files"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -l h_vmem=30G\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "swscan_EF" ]; then #fast (F) and low bigd bound (E)
    run_name=$xprefix.$chr
    run_cmd="swan_scan.R -f fast -e 1200 -c $chr_format$chr -n $gap_file -o $PWD/$xprefix.$chr $ref_file $xbam_files"
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -l h_vmem=40G\n$env_value\n#$ -cwd\n#$ -N j$run_name\n#$ -l h_vmem=30G\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=40gb,vmem=40gb,pmem=40gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "scscanTN" ]; then
    if [ ! -e "$PWD/${xprefix}_${yprefix}_scscanTN" ]; then
      mkdir "$PWD/${xprefix}_${yprefix}_scscanTN" #create subdir named $bam_file to save RData files
    fi
    #intrestingly -i 3,7 -j 30,70 produce more clusters than -i 5,7 -j 50,70
    run_cmd="sclip_scan.R -a -i 5,7 -j 50,70 -c $chr_format$chr -n $gap_file -o $PWD/${xprefix}_${yprefix}_scscanTN/${xprefix}_${yprefix}.$chr $ref_file $ybam_files:$xbam_files";
    run_name=${xprefix}_${yprefix}.chr$chr
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "scscanT" ]; then
    if [ ! -e "$PWD/${xprefix}_${yprefix}_scscanT" ]; then
      mkdir "$PWD/${xprefix}_${yprefix}_scscanT" #create subdir named $bam_file to save RData files
    fi
    run_cmd="sclip_scan.R -a -i 7 -j 70 -c $chr_format$chr -n $gap_file -o $PWD/${xprefix}_${yprefix}_scscanT/${xprefix}_${yprefix}.$chr $ref_file $ybam_files:$xbam_files";
    run_name=${xprefix}_${yprefix}.chr$chr
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "scscanN" ]; then
    if [ ! -e "$PWD/${xprefix}_${yprefix}_scscanN" ]; then
      mkdir "$PWD/${xprefix}_${yprefix}_scscanN" #create subdir named $bam_file to save RData files
    fi
    run_cmd="sclip_scan.R -a -i 3 -j 30 -c $chr_format$chr -n $gap_file -o $PWD/${xprefix}_${yprefix}_scscanN/${xprefix}_${yprefix}.$chr $ref_file $ybam_files:$xbam_files";
    run_name=${xprefix}_${yprefix}.chr$chr
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  if [ $swan_stage == "scscan" ]; then
    if [ ! -e "$PWD/${xprefix}_scscan" ]; then
      mkdir "$PWD/${xprefix}_scscan" #create subdir named $bam_file to save RData files
    fi
    run_cmd="sclip_scan.R -a -c $chr_format$chr -n $gap_file -o $PWD/${xprefix}_scscan/${xprefix}.$chr $ref_file $xbam_files"
    run_name=${xprefix}.chr$chr
    pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=16gb,vmem=16gb,pmem=16gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
    jobs[$cnt]=$pbs_cmd
    names[$cnt]=$run_name
    let cnt=$cnt+1
  fi

  ### why this section is needed? ###
  pbs_cmd=""
  for sp in $sp_files; do
    bam_files=$(echo $sp | tr "," "\t")
    for bam_file in $bam_files; do
      out_prefix=`basename $bam_file .bam`
      if [ $swan_stage == "svplot" ]; then #plot whole genome
        run_name=$out_prefix.$chr
        raw_bed="$data_dir/$out_prefix.$chr.sv.bed"
        x_ticks=`cat $raw_bed | cut -f 2 | tr "\n" ","`
        y_ticks=`cat $raw_bed | cut -f 3 | tr "\n" ","`
        x_ticks=${x_ticks%,}
        y_ticks=${y_ticks%,}
        run_cmd="swan_plot.R -a -q -l lCd,lDl,lDr,HAF,HAR -c $chr_format$chr -t $x_ticks:$y_ticks -r $data_dir/$out_prefix.$chr.rg1 $data_dir/$out_prefix.$chr.swan.txt.gz"
        #swan_plot.R -u 16146206 -v 16166107 -t 16148206:16164107 -l lDl,lDl,HAF,HAR -r dream_set3/set3.tumor.22.rg1 dream_set3/set3.tumor.22.swan.txt.gz
        pbs_cmd="#!/bin/bash\n#$ -q bigram\n#$ -cwd\n#$ -N j$run_name\n#$ -S /bin/bash\n#$ -j y\n#$ -m e\n#$ -M xiac\n#PBS -q main\n#PBS -d $PWD\n#PBS -N j$run_name\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m e\n#PBS -M lixia\n#PBS -l walltime=300:00:00\n#PBS -l nodes=1:ppn=8\n#PBS -l mem=8gb,vmem=8gb,pmem=8gb\n$run_cmd >$run_name.$swan_stage.log 2>&1\n"
      fi
      if [ -n "$pbs_cmd" ]; then
        #echo $run_cmd
        #echo "this shouldn't show up"
        jobs[$cnt]=$pbs_cmd
        names[$cnt]=$run_name
        let cnt=$cnt+1
      fi
    done
  done
done

if [ $cnt == 0 ]; then
  echo "$swan_stage not found!"
  exit
fi

echo "#total jobs: $cnt"
echo "#total jobs: ${#jobs[@]}"
let cnt=0
for job in $(seq 1 ${#jobs[@]}); do
  run_name=${names[$cnt]}
  pbs_cmd=${jobs[$cnt]}
  #echo "#run_name $run_name"
  #echo "#stage $swan_stage"
  #echo "#pbs_cmd $pbs_cmd"
  echo -e $pbs_cmd >$run_name.$swan_stage.pbs
  #if [ $dry_mode -eq "run" ]; then
  #  echo "#ssa.py $run_name.$swan_stage.pbs"
  #   ssa.py $run_name.$swan_stage.pbs
  #fi
  let cnt=$cnt+1
done
let elapse_time=$SECONDS-$start_time
echo "all $cnt submitted, use $elapse_time s"
