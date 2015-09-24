#!/usr/bin/env Rscript
#Rprof(filename="Rprof.out", line.profiling=TRUE) ##Eable to Profile
gtk=proc.time()[3] #global time keeper
#### Externalities ####
suppressMessages(library(Biobase))
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(IRanges))
suppressMessages(library(swan))
suppressMessages(library(seqCBS))
gmk=get_gmk(Sys.getpid())

#### CommandLines ####
option_list <- list(
  #all options cann't default to NA, no multiline help, -g not usable  !!!!!
  #options the same cross samples and read groups
  make_option(c("-c", "--chromosomeName"), default="11",
              help="chromosome to scan [default %default]"), 
  make_option(c("-u", "--scanStart"), type="integer", default=1,
              help="1-indexed scan start, [default %default]"),
  make_option(c("-v", "--scanEnd"), type="integer", default=max_chr_len,
              help="1-indexed scan end, [default %default]"),
  make_option(c("-w", "--windowWidth"), type="integer", default=width_fail,
              help="scan widow witdh, must be an integer >0 [default %default]"),
  make_option(c("-s", "--stepSize"), type="integer", default=10,
              help="scan window step size for the scan [default %default]"),
  make_option(c("-n", "--gap"), default="None",
              help="gap/N locations of hg19 in ucsc format [default %default]"),
  make_option(c("-o", "--spout"), default="input",
              help="sample output prefix, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)
parser=OptionParser(usage="%prog [options] ref_file spY.rg1.bam,spY.rg2.bam,...:spX.rg1.bam,spX.rg2.bam,...", option_list=option_list)
cat("-Info: invoking command:",commandArgs(),"\n")
args=commandArgs(trailingOnly=TRUE)
cmd=parse_args(parser,args,print_help_and_exit=TRUE,positional_arguments=TRUE)
#seqcbs_what=c("pos","mapq"); MAPQ.THRESH=0; seqcbs_pos0=c(); seqcbs_read0=IRanges()

if(length(cmd$args)!=2){ print_help(parser); quit(); }
ref_file=cmd$args[1]
inputs=cmd$args[2]
sample_tags=strsplit(inputs,split=':')[[1]]
n_sample=length(sample_tags); rg_files=strsplit(sample_tags,split=",")
if(cmd$options$spout=='input') { 
  sp_prefix=gsub("\\.\\w+$","",sapply(rg_files,lcPrefix),perl=T)
  #rg_prefix=rg_files
} else {
  sp_prefix=strsplit(cmd$options$spout,split=':')[[1]]
  #rg_prefix=lapply(seq_along(sp_prefix),function(i,a,b){paste(a[i],".rg",seq_along(b[[i]]),sep="")},sp_prefix,rg_files)
}
#if(!length(sp_prefix)>=2) stop("need at lease input two bam files and outputs\n")
#print(rg_files)
#print(rg_prefix)
#print(sp_prefix)
verbose=cmd$options$noQuiet #true to turn on
debug=cmd$options$debug
if(debug) options(warn=2)
if(cmd$options$windowWidth<=0) stop("error window width, retry")
width=cmd$options$windowWidth

seq_name=cmd$options$chromosomeName
if(cmd$options$stepSize<=0) stop("error sliding block step size, retry")
stepsize=cmd$options$stepSize
gap=read_gap(cmd$options$gap)
if(verbose) cat("-Info: gap_file=",cmd$options$gap,"\n")
## -u and -v
ref=FaFile(ref_file)
rg=scanFa(ref, param=scanFaIndex(ref)) #this returns a DNAStringSet
ref_len=length(rg[[seq_name]]) #[]->DNAStringSet, support width; [[]]->DNAString, support length
list[ref_head,ref_tail]=ref_ends(rg[[seq_name]])
rm(rg); tgc=gc();
start=max(1,cmd$options$scanStart)
end=min(ref_len,cmd$options$scanEnd)
win_pos=seq(start,end,stepsize)
n_wins=length(win_pos)

cat("-Info: doing scan,",taggie(gtk),"\n")
apply_par=list(rg_files=rg_files, sp_prefix=sp_prefix, gap=gap, seq_name=seq_name, 
               start=start, end=end, width=width, stepsize=stepsize,
               ref_head=ref_head, ref_tail=ref_tail, debug=debug, verbose=verbose, gtk=gtk)
do.call(apply_seqcbs, apply_par)
cat("=Info: done,",taggie(gtk),"\n")
cat("=warnings if any:\n")
warnings()
