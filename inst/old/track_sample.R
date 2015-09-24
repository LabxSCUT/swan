#!/usr/bin/env Rscript
suppressMessages(library(Biostrings))
suppressMessages(library(Rsamtools))
library(swan)
setwd("~/work/test/")

## change following
fa_dir = file.path(getwd(),"2013Sep18_fa")
sv_tag = "d5000000"
seq_name = "11"
sample_dir = file.path(getwd(),paste("2014Jun17_sample",sep="_"))
## mixing rates
r_range=c(1,0.5,0.25,0.05,0.00)
print(c("r length:",length(r_range),"values:",r_range))
## coverage
c_range=c(5,20,50,80)
print(c("c length:",length(c_range),"values:",c_range))
## replicates
p_size=3
p_seed=sample(1:10000, size=p_size, replace=FALSE)
p_seed=c(0001,0002,0003)
p0_size=0
p0_seed=sample(1:10000, size=p0_size, replace=FALSE)
cat(c("p0 length:",length(p0_seed),"\n"))
cat(c("p length:",length(p_seed),"\n"))
## isize
is_range=c(300,350,400)
print(c("is length:",length(is_range),"values:",is_range))
## isize_sd
min_sd=0; max_sd=2; sd_step=0.05
sd_range=t(floor(as.matrix(2**seq(min_sd,max_sd,1)*0.05)%*%(is_range)))
print(c("sd length:",length(sd_range),"values:",sd_range))
## read length
rl_range=c(75,100,150)
print(c("rl length:",length(rl_range),"values:",rl_range))

#tag = rxx_cxx_lxx_ixx_sxx_pxx
#glob all .fa files in specified dir
fa_files=Sys.glob(file.path(fa_dir,paste("*",sv_tag,"*.fas",sep="")))
if(length(fa_files)==0) cat("WARNING: no fa files selected \n")
dir.create(sample_dir, showWarnings = FALSE)
pbs_prefix="#!/bin/bash
#PBS -N %s
#PBS -o %s
#PBS -S /bin/bash
#PBS -j oe
#PBS -m e
#PBS -M lixia
#PBS -l walltime=299:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l mem=%d000mb,vmem=%d000mb,pmem=%d000mb
cd %s"
ref_file="~/hg/hg19/human_g1k_v37.fasta"
sampe_max_is=100000
#wgsim -d isize -s issd -N readnum -1 lenR1 -2 lenR2 -r mutrate -R indelrate \
#   -S rseed <in.ref.fa> <out.read1.fq> <out.read2.fq>
#wgsim_cmd="wgsim -r 0.0010 -R 0.15 -d %d -s %d -N %d -1 %d -2 %d -S %d %s %s %s"
wgsim_cmd="dwgsim -r 0.0010 -R 0.15 -y 0 -d %d -s %d -N %d -1 %d -2 %d -z %d %s %s"
merge_cmd="cat %s %s >%s"
#bwa index -a bwtsw $ref.fa
#bwa_index="#bwa index -a bwtsw %s"
bwa_index="sleep 1"
#bwa aln -q 15 -f $sam_1 $reference_fasta $fastq_file
#Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
bwa_aln="bwa aln -f %s %s %s &"
bwa_mem="bwa mem -M -t %s %s %s %s | samtools view -Sbh - | samtools sort - %s"
#bwa sampe -s -a $max_insert_size -f $sam $ref $sai1 $sai2 $fq1 $fq2; remove -s to enable sw-align of unmapped reads
bwa_sampe="bwa sampe -a %d -f %s %s %s %s %s %s" #use -s to disable SW for unmapped mate
#samtools view -bSu $sam_2 | samtools sort -n -o - samtools_nsort_tmp | 
#  samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - samtools_csort_tmp |
#  samtools fillmd -u - $reference_fasta > $bam_file
#samtl_cmd="samtools view -bSu %s | samtools sort -n -o - %s | samtools fixmate /dev/stdin /dev/stdout | samtools sort -o - %s | samtools fillmd -b - %s > %s"
samtl_cmd="samtools view -bSu %s | samtools sort -m %s -o - samtools_nsort_tmp > %s"
#this can be shortened
samtl_idx="samtools index %s"
samtl_stat="samtools flagstat %s"
chk_exit="exit_code=$?; if [[ $exit_code != 0 ]]; then echo $exit_code; exit $exit_code; fi"
sep_cmd="echo --"
set_time="SECONDS=0"
chk_time='echo "that took approximately $SECONDS seconds"'
chk_mem_cmd='peak_mem.sh $!'

total_sample_per_rset=length(fa_files)*length(p_seed)
total_sample_per_r0set=length(p0_seed)
print(paste("total_sample_per_rset=",total_sample_per_rset))
print(paste("total_sample_per_r0set=",total_sample_per_r0set))
total_set=length(c_range)*length(rl_range)*length(is_range)*length(sd_range)*
  ((length(r_range)-1)*total_sample_per_rset+total_sample_per_r0set)
print(paste("total_set=",total_set))
cat("p_seed",p_seed,"\n")
cat("p0_seed",p0_seed,"\n")

mem1=32
ppn1=4

#will produce
#suppose for fa files indexed
#sv_tag/sample_tag/sv_tag.fa.sample_tag.pseed.pbs
#sv_tag/sample_tag/sv_tag.fa.sample_tag.pseed.abd
#sv_tag/sample_tag/sv_tag.fa.sample_tag.pseed.bam

for(fa_file in fa_files){
  fa_prefix=basename(fa_file)
  #ref_file = paste(c(head(unlist(strsplit(fa_file,"\\.")),-1),"fna"),collapse=".")
  fa = FaFile(fa_file)  ## index is some.fa.gz.
  sved = scanFa(fa)
  ref = FaFile(ref_file)
  rg = scanFa(ref,param=scanFaIndex(ref))
  ref_len = length(rg[[seq_name]])
  for(i in seq(1,length(r_range))){
    for(j in seq(1,length(c_range)))
      for(k in seq(1,length(rl_range)))
        for(l in seq(1,length(is_range)))
          for(m in seq(1,length(sd_range[l,]))){
            key=c("r","c","rl","is","sd")
            value=c(r_range[i]*10000,c_range[j],rl_range[k],
                    is_range[l],sd_range[l,m])
            sample_tag=paste(paste(key,value,sep=""),collapse="_")
            sub_dir=file.path(sample_dir,sample_tag)
            if(file.exists(sub_dir)) unlink(sub_dir, recursive=T, force=T)
            dir.create(sub_dir, showWarnings = FALSE)
            if(i==0)
              p_range=p0_seed
            else
              p_range=p_seed
            for(n in seq(1,length(p_range))){
              pseed = paste("ps",p_range[n],sep="")
              fq_prefix = paste(fa_prefix,sample_tag,pseed,sep=".")
              log_file = file.path(sub_dir, paste(fq_prefix,"w","log",sep="."))
              pbs_file = file.path(sub_dir, paste(fq_prefix,"w","pbs",sep="."))
              abd_file = file.path(sub_dir, paste(fq_prefix,"abd",sep="."))
              abd_tmp=paste( paste(c(names(sved),names(rg)),c(1-r_range[i],r_range[i]),sep="\t")
                            ,collapse="\n")
              cat(abd_tmp,file=abd_file)
              #cat("writing ",fq_prefix,".* ...\n")
              ref_fq_file = paste(fq_prefix,"ref",sep=".")
              ref_bfast_file = paste(fq_prefix,"ref","bfast","fastq",sep=".")
              ref_muta_file = paste(fq_prefix,"ref","mutations","txt",sep=".")
              ref_fq1_file = paste(fq_prefix,"ref","bwa","read1","fastq",sep=".")
              ref_fq2_file = paste(fq_prefix,"ref","bwa","read2","fastq",sep=".")
              sved_fq_file = paste(fq_prefix,"sved",sep=".")
              sved_muta_file = paste(fq_prefix,"sved","mutations","txt",sep=".")
              sved_bfast_file = paste(fq_prefix,"sved","bfast","fastq",sep=".")
              sved_fq1_file = paste(fq_prefix,"sved", "bwa","read1","fastq",sep=".")
              sved_fq2_file = paste(fq_prefix,"sved", "bwa","read2","fastq",sep=".")
              fq1_file = paste(fq_prefix,"fq1",sep=".")
              fq2_file = paste(fq_prefix,"fq2",sep=".")
              sai1_file = paste(fq_prefix,"sai1",sep=".")
              sai2_file = paste(fq_prefix,"sai2",sep=".")
              sort_file = paste(fq_prefix,"sort",sep=".")
              sam_file = paste(fq_prefix,"sam",sep=".")
              bam_file = paste(fq_prefix,"bam",sep=".")
              pbs_cmd_tmp=gettextf(pbs_prefix, fq_prefix, log_file,
                                    ppn1, mem1, mem1, mem1, sub_dir)
              ref_num_reads = round((1-r_range[i])*(c_range[j]*ref_len)/rl_range[k]/2)
              #cat("ref_num_reads=", ref_num_reads,"\n")
              sved_num_reads = round(r_range[i]*(c_range[j]*ref_len)/rl_range[k]/2)
              #cat("sved_num_reads=", sved_num_reads,"\n")
              #wgsim_cmd_ref=gettextf(wgsim_cmd, is_range[l], sd_range[l,m], ref_num_reads,
              #                       rl_range[k], rl_range[k], p_range[n], ref_file, 
              #                       ref_fq1_file, ref_fq2_file)
              #wgsim_cmd_sved=gettextf(wgsim_cmd, is_range[l], sd_range[l,m], sved_num_reads,
              #                        rl_range[k], rl_range[k], p_range[n], fa_file, 
              #                        sved_fq1_file, sved_fq2_file)
              x=(1-r_range[i])*(c_range[j])*ref_len/(rl_range[k])/2
              y=(r_range[i])*(c_range[j])*ref_len/(rl_range[k])/2
              cat("r=",r_range[i],"c=",c_range[j],"rl=", rl_range[k], "ref_len=", ref_len, "\n")
              cat("ref=", x, "sved=", y,"\n")
              cat("ref=", round(x), "sved=", round(y),"\n")
              cat("sved=", ref_fq_file, "\n")
              cat("ref_file=", ref_file, "\n")
              if(ref_num_reads==0) {
                wgsim_cmd_ref=gettextf("touch %s %s",ref_fq1_file,ref_fq2_file)
              } else { #wgsim: -d isize=insert_size; dwgsim: -d isize=insert_size-2*rl
                wgsim_cmd_ref=gettextf(wgsim_cmd, is_range[l]-2*rl_range[k], sd_range[l,m], ref_num_reads,
                                     rl_range[k], rl_range[k], p_range[n], ref_file, 
                                     ref_fq_file)
              }
              if(sved_num_reads==0) {
                wgsim_cmd_sved=gettextf("touch %s %s",sved_fq1_file,sved_fq2_file)
              } else {
                wgsim_cmd_sved=gettextf(wgsim_cmd, is_range[l]-2*rl_range[k], sd_range[l,m], sved_num_reads,
                                       rl_range[k], rl_range[k], p_range[n], fa_file, 
                                       sved_fq_file)
              }
              merge_cmd_fq1=gettextf(merge_cmd, ref_fq1_file, sved_fq1_file, fq1_file)
              merge_cmd_fq2=gettextf(merge_cmd, ref_fq2_file, sved_fq2_file, fq2_file)
              bwa_index_tmp=gettextf(bwa_index, ref_file)
              bwa_mem_tmp=gettextf(bwa_mem,ppn1, ref_file,fq1_file,fq2_file,substr(bam_file,1,nchar(bam_file)-4))
              bwa_aln1_tmp="sleep 1"
              bwa_aln2_tmp="sleep 1"
              bwa_sampe_tmp="sleep 1"
              samtl_cmd_tmp="sleep 1"
              #bwa_aln1_tmp=gettextf(bwa_aln, sai1_file, ref_file, fq1_file)
              #bwa_aln2_tmp=gettextf(bwa_aln, sai2_file, ref_file, fq2_file)
              #bwa_sampe_tmp=gettextf(bwa_sampe, sampe_max_is, sam_file, ref_file,
              #                       sai1_file, sai2_file, fq1_file, fq2_file)
              #samtl_cmd_tmp=gettextf(samtl_cmd, sam_file, format(mem1*1000^3,scientific=F), bam_file)
              samtl_stat_tmp=gettextf(samtl_stat, bam_file)
              samtl_idx_tmp=gettextf(samtl_idx, bam_file)
              rm_cmd_tmp=gettextf("rm -f %s %s %s %s %s %s %s %s %s %s %s %s %s", sai1_file, sai2_file, 
                                  ref_fq1_file, ref_fq2_file, sved_fq1_file, sved_fq2_file,
                                  sam_file, fq1_file, fq2_file, 
                                  ref_bfast_file, sved_bfast_file, ref_muta_file, sved_muta_file)
              sh_tmp=paste(pbs_cmd_tmp, 
                           sep_cmd, set_time, wgsim_cmd_ref,  chk_exit, chk_time,
                           sep_cmd, set_time, wgsim_cmd_sved, chk_exit, chk_time,  
                           sep_cmd, set_time, merge_cmd_fq1, chk_exit, chk_time, 
                           sep_cmd, set_time, merge_cmd_fq2, chk_exit, chk_time,  
                           sep_cmd, set_time, bwa_index_tmp, chk_exit, chk_time,  
                           sep_cmd, set_time, bwa_mem_tmp, chk_exit, chk_time,  
                           sep_cmd, set_time, bwa_aln1_tmp, chk_exit, chk_time, 
                           sep_cmd, set_time, bwa_aln2_tmp, chk_exit, chk_time,  
                           sep_cmd, set_time, bwa_sampe_tmp, chk_exit, chk_time, 
                           sep_cmd, set_time, samtl_cmd_tmp, chk_exit, chk_time, 
                           sep_cmd, set_time, samtl_idx_tmp, chk_exit, chk_time, 
                           sep_cmd, set_time, samtl_stat_tmp, chk_exit, chk_time, 
                           sep_cmd, rm_cmd_tmp, sep="\n")
              cat(sh_tmp, file=pbs_file)
              #quit() #for debug
            }  #p_range
          }    #sd_range[m,]
    }          #r_range
}              #fa_files

#Read names explained
#Read names are of the form:
#@<#1>_<#2>_<#3>_<#4>_<#5>_<#6>_<#7>_<#8>:<#9>:<#10>_<#11>:<#12>:<#13>_<#14>
#1. contig name (chromsome name)
#2. start end 1 (zero-based)
#3. start end 2 (zero-based)
#4. strand end 1 (0 - forward, 1 - reverse)
#5. strand end 2 (0 - forward, 1 - reverse)
#6. random read end 1 (0 - from the mutated reference, 1 - random)
#7. random read end 2 (0 - from the mutated reference, 1 - random)
#8. number of sequencing errors end 1 (color errors for colorspace)
#9. number of SNPs end 1
#10. number of indels end 1
#11. number of sequencing errors end 2 (color errors for colorspace)
#number of SNPs end 2
#number of indels end 2
#read number (unique within a given contig/chromsome)
cat("warnings if any:\n")
warnings()

