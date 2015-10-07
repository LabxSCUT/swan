#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
#Rprof(filename="Rprof.out", line.profiling=TRUE) ##Eable to Profile
#### Externalities ####
library(swan)
for(p in c("optparse")) myrequire(p)
for(p in c("Biobase","GenomicRanges","Rsamtools","IRanges")) myrequire(p,repo="Bioc")
#suppressMessages(library(Biobase))
#suppressMessages(library(optparse))
#suppressMessages(library(GenomicRanges))
#suppressMessages(library(Rsamtools))
#suppressMessages(library(IRanges))

#	tryCatch({
#    source("~/setup/swan/R/libswan.R")
#    source("~/setup/swan/R/svlib.r")
#    source("~/setup/swan/R/ThresholdFunctions.R")
#    dyn.load("~/setup/swan/src/swan.so")
#    #suppressMessages(library(swan))
#    cat("-Info: loaded from ~/setup/swan package!\n")
#  }, error=function(e) {
#    cat("-Warn: cant load from ~/setup/swan package!\n")
#    suppressMessages(library(swan))
#    cat("-Info: loaded from installed swan package!\n")
#  })

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
  make_option(c("-r", "--mixingRate"), default=mixing_rate_fail, #can't default to NA
              help="mixing rates, [default %default]"),
  make_option(c("-w", "--windowWidth"), type="integer", default=width_fail,
              help="scan widow witdh, must be an integer >0 [default %default]"),
  make_option(c("-g", "--lwWindowWidth"), type="integer", default=lw_width_fail,
              help="Lw scan widow witdh, must be an integer >0 [default %default]"),
  make_option(c("-s", "--stepSize"), type="integer", default=10,
              help="scan window step size for the scan [default %default]"),
  make_option(c("-n", "--gap"), default="",
              help="gap/N locations of hg19 in ucsc format [default %default]"),
  make_option(c("-k", "--stat"), action="store_true",default=FALSE,
              help="provide precomputed stat file to disable tracks, [default %default]"),
  #options may be different cross samples and read groups
  make_option(c("-x", "--propClip"), default="learn",
           help="required aligned length soft isize; e.g. RL+1=>none, 0=>all [default %default, .5xRL]"),
  make_option(c("-y", "--hangClip"), default="0",
          help="required aligned length soft hang; e.g. RL+1=>all, 0=>none [default %default, RL-20]"),
  make_option(c("-b", "--coverageMean"), default="learn",
              help="coverage mean, [default %default]"),
  make_option(c("-l", "--readLength"), default="learn",
              help="read length, [default %default]"),
  make_option(c("-i", "--insertSize"), default="learn",
              help="biological insert size mean,sdR,sdL [default %default]"),
  make_option(c("-m", "--marginDelta"), default="learn",
              help="margin/delta size [default %default, IS+6*ISSD]"),
  make_option(c("-e", "--bigDel"), default="learn",
              help="big deletion size [default %default, IS+3*ISSD]"),
  make_option(c("-p", "--probHang"), default="learn",
              help="global probablity seeing hang read [default %default]"),
  make_option(c("-d", "--probSoft"), default="learn",
              help="global probablity seeing soft read [default %default]"),
  make_option(c("-t", "--otherOpt"), default="smallDel=20,smallIns=20,maxInsert=learn,multiCore=1",
              help="other options [default %default]"),
  #options for computing tuning
  make_option(c("-z", "--trunkSize"), type="integer", default=1000000, #can't default to NA
              help="trunk size for processing scanning bamfile, for 50x within 8G mem use, must be multiples of stepsize -s and blocksize -k [default %default]"),
  make_option(c("-o", "--spout"), default="input",
              help="sample output prefix, [default %default]"),
  make_option(c("-f", "--fastSave"), default="normal",
              help="compute fast, can use normal, fast or super [default %default]"),
  make_option(c("-j", "--memSave"), action="store_true", default=FALSE,
              help="save memory, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)
parser=OptionParser(usage="%prog [options] ref_file sp.rg1.bam,sp.rg2.bam,sp.rg3.bam,... ", option_list=option_list)
cat("-Info: swan_scan vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args=commandArgs(trailingOnly=TRUE)
cmd=parse_args(parser,args,print_help_and_exit=TRUE,positional_arguments=TRUE)
#cat("-Info: input options:\n")
#print(cmd$options)
if(length(cmd$args)!=2){ print_help(parser); quit(); }
ref_file=cmd$args[1]
inputs=cmd$args[2]
debug=cmd$options$debug
verbose=cmd$options$noQuiet

if(debug){ options(warn=2) } else { options(warn=-1) }
gtk=proc.time()[3]
gmk=get_gmk(Sys.getpid())

usestat=cmd$options$stat
spout_text=cmd$options$spout
mem_save=cmd$options$memSave
fastSave=cmd$options$fastSave
trunk_size=cmd$options$trunkSize
seq_name=cmd$options$chromosomeName
stepsize=cmd$options$stepSize
mixing_rate=cmd$options$mixingRate
width=cmd$options$windowWidth
lw_width=cmd$options$lwWindowWidth
gap_text=cmd$options$gap
scan_start=cmd$options$scanStart
scan_end=cmd$options$scanEnd
Delta_global=cmd$options$marginDelta
RL_global=cmd$options$readLength
isize_global=cmd$options$insertSize
isize_sd_global=-1
coverage_global=cmd$options$coverageMean
p_global=cmd$options$probHang
q_global=cmd$options$probSoft
hang_clip_global=cmd$options$hangClip
prop_clip_global=cmd$options$propClip
bigDel_global=cmd$options$bigDel
other_opt_text=cmd$options$otherOpt

### debug options
#ref_file="../real/human_g1k_v37.ucsc.fasta";inputs="NA12878.chr22.16M-21M.bam";verbose=TRUE;spout_text="input";debug=TRUE;mem_save=FALSE;fastSave="normal";trunk_size=5000000;seq_name="chr22";stepsize=10;mixing_rate=0.1;width=10;scan_start=16000001;scan_end=21000000;Delta_global="learn";RL_global="learn";isize_global="learn";coverage_global="learn";p_global="learn";q_global="learn";hang_clip_global="learn";prop_clip_global="learn";bigDel_global="learn";gap_text="../data/human_g1k_v37.ucsc.fasta.gap.txt.gz";usestat=TRUE;other_opt_text="smallDel=20,smallIns=20"

### global constants ###
what=c("pos","mpos","isize","strand","qwidth","cigar","qname","flag")
pattern="(?<key>[a-z]+)(?<value>[0-9]+)"

### parsing positional arguments ###
sample_tags=strsplit(inputs,split=':')[[1]]
n_sample=length(sample_tags); rg_files=strsplit(sample_tags,split=",") #rg[[sp]][rg] -> bam #now we only allow n_sample==1
if(n_sample!=1) stop("only one sample allowed!")
if(spout_text=='input') {  #if -o is not set
  sp_prefix=gsub("\\.\\w+$","",sapply(rg_files,lcPrefix),perl=T) #get common prefix of rg_files; 
} else {  #if -o is set
  sp_prefix=strsplit(spout_text,split=':') #sample prefixes
}
if(length(sp_prefix)!=length(rg_files)) stop("=Error: # of output samples not equal to input samples\n") 
#for(i in seq_along(rg_files)) stopifnot(length(rg_files[[i]])==length(rg_prefix[[i]])) #make sure in/output prefixes consistent
block_size=trunk_size

### parsing global scan options ###
other_opt=parse_opt2(other_opt_text)[[1]]
smallDel_global=other_opt[["smallDel"]]
smallIns_global=other_opt[["smallIns"]]
maxInsert_global=other_opt[["maxInsert"]]
multiCore=as.integer(other_opt[["multiCore"]])
if(stepsize<=0) stop("error sliding block step size, retry")
if(mixing_rate>1 | mixing_rate<0) stop("error mixing rate, retry")
if(width<=0) stop("error window width, retry")
if(lw_width<=0) stop("error window width, retry")
gap=read_gap(gap_text)
if(verbose) cat("-Info: gap_file=",gap_text,"\n")
ref=FaFile(ref_file)
rg=scanFa(ref, param=scanFaIndex(ref)) #this returns a DNAStringSet
if(seq_name %in% c("A","a","ALL","all","All")){
	seq_name = names(rg)[seq_len(min(length(rg),24))]
	cat("all chrs:",seq_name,"\n")
	ref_len = sapply(seq_name, function(x) { length(rg[[x]]) })
	start = rep(1,length(ref_len))
	end = ref_len
} else if(grepl(",",seq_name)){
	seq_name = strsplit(seq_name,",")[[1]]
	cat("all chrs:",seq_name,"\n")
	ref_len = lapply(seq_name, function(x) { length(rg[[x]]) })
	start = rep(1,length(ref_len))
	end = ref_len
} else {
	cat("all chrs:",seq_name,"\n")
	ref_len = length(rg[[seq_name]]) #[]->DNAStringSet, support width; [[]]->DNAString, support length
	start = max(1,scan_start); end=min(ref_len,scan_end)
}
### apply swan scan ###
# apply_scan will take proper arguments and do the scan job
# scan from one sample will be merged

apply_par = list(rg_files=rg_files,
                 what=what,isize_global=isize_global,
                 coverage_global=coverage_global,
                 RL_global=RL_global,gap=gap,trunk_size=trunk_size,
                 width=width,lw_width=lw_width,mixing_rate=mixing_rate,stepsize=stepsize,block_size=block_size,
                 Delta_global=Delta_global,bigDel_global=bigDel_global,
                 smallDel_global=smallDel_global,smallIns_global=smallIns_global,maxInsert_global=maxInsert_global,
                 p_global=p_global,q_global=q_global,
                 hang_clip_global=hang_clip_global,prop_clip_global=prop_clip_global,
                 fast=fastSave,usestat=usestat,debug=debug,verbose=verbose,gtk=gtk) 

cat("-Info: doing scan,",taggie(gtk),"\n")
apply_chr = function(i,arg,chr,start,end,rg_files,sp_prefix){ #apply one chromosome
	#arg[["sp_prefix"]] = paste(sp_prefix, paste("chr",gsub(x=chr[i],pattern="chr",replacement=""),sep=""), sep=".")
	arg[["sp_prefix"]] = paste(sp_prefix, chr[i], sep=".")
  arg[["rg_prefix"]] = lapply(seq_along(arg[["sp_prefix"]]),function(i,a,b){paste(a[i],".rg",seq_along(b[[i]]),sep="")},arg[["sp_prefix"]],rg_files)
 	#rg_prefix=lapply(rg_files,function(x) { x=gsub(".bam$","",x); if(length(x)==1) paste(x,"rg1",sep="."); x }) 
	cat("-Info: apply_chr:", chr[i],", rg_files:\n"); print(rg_files)
	cat("-Info: apply_chr:", chr[i],", sp_prefix:\n"); print(arg[["sp_prefix"]])
	cat("-Info: apply_chr:", chr[i],", rg_prefix:\n"); print(arg[["rg_prefix"]])
	arg[["seq_name"]]=chr[i]; arg[["start"]]=start[i]; arg[["end"]]=end[i]
	tmp=do.call(apply_scan,arg)
}
if(!debug){
  tmp=lapply(seq_along(seq_name),apply_chr,apply_par,seq_name,start,end,rg_files,sp_prefix)
} else {
  tmp=mclapply(seq_along(seq_name),apply_chr,apply_par,seq_name,start,end,rg_files,sp_prefix,mc.cores=multiCore)
}
cat("=Info: done,",taggie(gtk),"\n")
if(debug) { cat("=Info: warnings if any\n"); warnings() }
cat("please find outputs in: {",paste(sp_prefix,collapse=',',sep=':'),"} . {",paste(seq_name,sep="",collapse=","),"} . {swan.txt.gz, bigd.txt, disc.txt, anch.txt and par.txt} files","\n",sep=" ")
cat("SWANDONE\n")
