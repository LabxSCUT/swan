#!/usr/bin/env Rscript
gtk=proc.time()[3] #global time keeper
#### Externalities ####
suppressMessages(library(Biobase))
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(swan))
suppressMessages(library(parallel))
options(warn=2)

## UNDO:
#source("~/software/swan/R/libswan.R")
#source("~/software/swan/R/svlib.r")

#MIN.READS.PER.CLUSTER=5
#MIN.BASES.PER.CLUSTER=50
SC.PROPORTION.MISMATCH.THRESH=0.05
MIN.GAP.DISTANCE=10000
#SC.MIN.CONSENSUS.LENGTH=30

#### CommandLines ####
option_list <- list(
  #all options cann't default to NA, no multiline help, -g not usable  !!!!!
  #options the same cross samples and read groups
  make_option(c("-c", "--chromosomeName"), default="11",
              help="chromosome to scan [default %default]"),
  make_option(c("-n", "--gap"), default="None",
              help="gap/N locations of hg19 in ucsc format [default %default]"),
  make_option(c("-i", "--minReadPerCluster"), type="integer", default=3,
              help="minimal number of reads per cluster, [default %default]"),
  make_option(c("-j", "--minBasePerCluster"), type="integer", default=30,
              help="minimal number of total bases per cluster, [default %default]"),
  make_option(c("-u", "--scanStart"), type="integer", default=1,
              help="1-indexed scan start, [default %default]"),
  make_option(c("-v", "--scanEnd"), type="integer", default=max_chr_len,
              help="1-indexed scan end, [default %default]"),
  make_option(c("-z", "--trunkSize"), type="integer", default=1000000,
              help="trunk size for scanning bamfile [default %default]"),
  make_option(c("-o", "--spout"), default="input",
              help="sample output prefix, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]")
)
parser=OptionParser(usage="%prog [options] sp.rg1.bam,sp.rg2.bam,sp.rg3.bam...", option_list=option_list)
cat("-Info: invoking command:",commandArgs(),"\n")
args=commandArgs(trailingOnly=TRUE)
cmd=parse_args(parser,args,print_help_and_exit=TRUE,positional_arguments=TRUE)
if(length(cmd$args)!=1){ print_help(parser); quit(); }
verbose=cmd$options$noQuiet #true to turn on
inputs=cmd$args[1]
sample_tags=strsplit(inputs,split=':')[[1]]
n_sample=length(sample_tags); rg_files=strsplit(sample_tags,split=",")
trunk_size=cmd$options$trunkSize  
seqname=cmd$options$chromosomeName
MIN.READS.PER.CLUSTER=cmd$options$minReadPerCluster
MIN.BASES.PER.CLUSTER=cmd$options$minBasePerCluster
scan_start=max(1,cmd$options$scanStart)
scan_end=cmd$options$scanEnd
gap=read_gap(cmd$options$gap)
if(cmd$options$spout=='input') {
  sp_prefix=gsub("\\.\\w+$","",sapply(rg_files,lcPrefix),perl=T)
} else {
  sp_prefix=strsplit(cmd$options$spout,split=':')[[1]]
}
files=rg_files[[1]]  ### Check: is this right?
debug=cmd$options$debug
cat("files=",files,"\n")
## UNDO: set this back
#debug=TRUE
#files= "../../synthetic/set1.v2/normal/ae16ceb3-ce31-4648-840c-66f3c5d180a6/synthetic.challenge.set1.normal.v2.bam"
#gap=NULL
#scan_start=1
#trunk_size=500000

header=scanBamHeader(files[1])
chrnames=names(header[[1]]$targets)  # find out what the chromosomes are called.
idx = which(chrnames==seqname)
if(length(idx)==0){
    cat("\n\nCan not find ",seqname," among target names in bam file.\n",sep="")
    cat("Target names: ",chrnames,"\n")
    cat("Your input: ",seqname,"\n")
    quit()
}
n_trunks=ceiling((scan_end-scan_start+1)/trunk_size);
cat("scan_start=",nsf(scan_start),"scan_end=",nsf(scan_end),"n_trunks=",n_trunks,"\n")

rcpp_core_sclip = function(ti) {
  .Call("core_sclip", ti, n_trunks,scan_start,scan_end,trunk_size,files,seqname,gap, MIN.READS.PER.CLUSTER, MIN.BASES.PER.CLUSTER, SC.PROPORTION.MISMATCH.THRESH, MIN.GAP.DISTANCE, PACKAGE="swan");
}

if(debug) {
  sclip_res=lapply(seq(n_trunks), rcpp_core_sclip)
} else { #accelerate using mclapply
  sclip_res=mclapply(seq(n_trunks), rcpp_core_sclip)
}

var_names=ls(sclip_res[[1]])
var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },sclip_res[[1]])
all_RData=lapply(seq_along(var_names),new_merge_RData,var_names,var_classes,sclip_res)


#for(di in seq_along(all_RData)) { all_RData[[di]]=all_RData[[di]][[1]] }
names(all_RData)=var_names
#names(all_RData)
#all_Rdata
#class(sclip_res[[1]][["scL"]])
#class(sclip_res[[1]][["scR"]])
#lapply(all_RData,length)
#lapply(all_RData,dim)
#sum(sapply(seq_along(sclip_res),function(x,tag){length(sclip_res[[x]][[tag]])},"scL"))
#sum(sapply(seq_along(sclip_res),function(x,tag){length(sclip_res[[x]][[tag]])},"scR"))

var_names=ls(sclip_res[[1]])
var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },sclip_res[[1]])
if(debug) {
  cat("scL","scR","\n")
  for(res in sclip_res) { cat(length(res$scL),length(res$scR),"\n") }
  all_RData=lapply(seq_along(var_names),new_merge_RData,var_names,var_classes,sclip_res)
}
for(di in seq_along(all_RData)) { all_RData[[di]]=all_RData[[di]][[1]] }
names(all_RData)=var_names
cat("scL","scR","\n")
cat(length(all_RData$scL),length(all_RData$scR),"\n")

#add in chrom information; so sc1L and sc1R always on the same chromosome
for(i in seq_along(all_RData$scL)){ all_RData$scL[[i]]$chrom=seqname }   
for(i in seq_along(all_RData$scR)){ all_RData$scR[[i]]$chrom=seqname }

#save the results here

sclip_file=paste(sp_prefix[1],seqname,"RData",sep=".")
dir.create(dirname(sclip_file),recursive=TRUE,showWarnings=FALSE)
save(MIN.READS.PER.CLUSTER,MIN.BASES.PER.CLUSTER,
  SC.PROPORTION.MISMATCH.THRESH,MIN.GAP.DISTANCE,
  all_RData,file=sclip_file)
cat("output to",sclip_file,"\n")
cat("done sclip_scan",taggie(gtk),"\n")
#sp_par=as.data.frame(list(chr=seqname,start=scan_start,end=scan_end,trunk_size=trunk_size,min_reads=MIN.READS.PER.CLUSTER,min_bases=MIN.BASES.PER.CLUSTER,max_mismatch=SC.PROPORTION.MISMATCH.THRESH,min_gapdist=MIN.GAP.DISTANCE))
#parf=paste(sp_prefix[2],"sclip","par","txt",sep=".")
#write.table(format(sp_par,digits=2,scientific=F),parf,quote=F,row.names=F,sep="\t")
cat("warnings if any:\n")
warnings()
