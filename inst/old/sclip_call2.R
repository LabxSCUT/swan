#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
gtk=proc.time()[3]
suppressMessages(library(Biostrings))
suppressMessages(library(swan))
suppressMessages(library(parallel))
gmk=get_gmk(Sys.getpid())
options(warn=2)

####CommandLines####
option_list2 <- list(
    make_option(c("-d", "--contdir"), default="none",
              help="contrast directory, [default %default]"),
    make_option(c("-n", "--gapfile"), default="none",
              help="gap annotation file, [default %default]"),          
    make_option(c("-o", "--outd"), default="sclip_call",
              help="output files prefix, [default %default]"),
    make_option(c("-c", "--trunk"), default="1:1",
              help="number of trunks:trunk number, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] reffile clusterdir", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cat("-Info: sclip_call2.R vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
#quit if not valid cmd args
if(length(cmd$args)!=2){ 
  print_help(parser)
  quit()
}
temp=strsplit(cmd$options$trunk,split=":")[[1]]
ntrunks=as.numeric(temp[1]);trunknum=as.numeric(temp[2])
cat("ntrunks=",ntrunks,", trunknum=",trunknum,"\n")
if(trunknum>ntrunks){
    cat("-c option must be ntrunks:trunknum, where trunknum<=ntrunks.")
    quit()
}
ref_file=cmd$args[1]
clusterdir1=cmd$args[2]
gapfile = cmd$options$gapfile
if(gapfile=="none") gapfile=NULL
clusterdir0 = cmd$options$contdir
if(clusterdir0=="none") clusterdir0=NULL
outprefix=cmd$options$outd
rdatfile = paste(outprefix,".",trunknum,".RData",sep="")
dir.create(dirname(rdatfile),recursive=TRUE,showWarnings=FALSE)
cat("clusterdir1=",clusterdir1,", clusterdir0=",clusterdir0,"\n")
cat("Writing to ",rdatfile,"\n")

#### comment this:
#clusterdir1 = "sclip_scan"
#clusterdir0 = NULL
#rdatfile="sclip_calls.RData"
#gapfile="~/software/swan/data/human_g1k_v37.fasta.gap.txt.gz"
#ref_file="~/structure/xiac/hg/hg19/human_g1k_v37.fasta"
#source("~/software/swan/R/libswan.R")
#source("~/software/swan/R/svlib.r")
debug=TRUE

###############################################
# Tuning parameters
###############################################

SC.MIN.CONSENSUS.LENGTH=30
MAX.MATCHES=10
MAX.DIST.CONSENSUS.CLUSTER=100 # Maximum distance allowed between consensus and SC cluster in destination.
REMAP.WIN=1000 # window for remapping of destination consensus back to home.

###############################################
#  Get the reference template.
###############################################
rg = .Call("scanFa", ref_file, PACKAGE="swan")
setClass("BWA")
tmp = setMethod("vmatchPattern", "BWA",
  function(pattern, subject,
           max.mismatch=1, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
           algorithm="auto") {
    if (is(pattern, "XString"))
      pattern = toString(pattern)
    # NOTE: this function only supports max.mismatch = 1 for now.
    stopifnot(max.mismatch == 1)
    a = .Call("matchPattern", pattern, subject, PACKAGE="swan")
    lapply(a, function (x) { IRanges(start=x[[1]], width=x[[2]]) })
  }
)

###############################################
# Combine sclip clusters over all chromosomes.
# Load all bam files in directory.
###############################################
### Control directory
if(!is.null(clusterdir0)){
    files=dir(clusterdir0)
    sel=grep(".RData",files,fixed=TRUE)
    load(file.path(clusterdir0,files[sel[1]]))
    scL=all_RData$scL
    scR=all_RData$scR
    if(length(sel)>=2){ #only needed if >2 sscan files
      for(i in 2:length(sel)){
        load(file.path(clusterdir0,files[sel[i]]))
        scL = c(scL, all_RData$scL)
        scR = c(scR, all_RData$scR)
      }
    }
    scL0 = scL
    scR0 = scR
    sctabL0 = summarizeClusters(scL0,includeNBases=TRUE)
    sctabR0 = summarizeClusters(scR0,includeNBases=TRUE)
    cat(length(scL0),"left clip clusters loaded for control.\n")
    cat(length(scR0),"right clip clusters loaded for control.\n")
    cat("done loading control\n")
}
### Case directory
files=dir(clusterdir1)
sel=grep(".RData",files,fixed=TRUE)
envs = lapply(sel, function(i) {
  env = new.env()
  load(file.path(clusterdir1, files[i]), env)
  env
})
scL = do.call(c, lapply(envs, function(env) { env$all_RData$scL }))
scR = do.call(c, lapply(envs, function(env) { env$all_RData$scR }))
rm(envs)
cat("done loading control\n")
scL1 = scL
scR1 = scR
sctabL1 = summarizeClusters(scL1,includeNBases=TRUE)
sctabR1 = summarizeClusters(scR1,includeNBases=TRUE)

cat(length(scL1),"left clip clusters loaded for case.\n")
cat(length(scR1),"right clip clusters loaded for case.\n")

###############################################
# Filtering
###############################################

### Remove centromere and gap regions.
gap=read_gap(gapfile)
seqnames=unique(gap$chrom)
throwL1 = rep(FALSE,length(scL1)); throwR1 = rep(FALSE,length(scR1))
throwL0 = rep(FALSE,length(scL1)); throwR0 = rep(FALSE,length(scR1))
for(s in seqnames){
    sel2 = which(gap$chrom==s)
    thisgap = IRanges(start=gap$chromStart[sel2],end=gap$chromEnd[sel2])
    sel = which(sctabL1$Chrom == s); pos = IRanges(start=sctabL1$Position[sel],width=1)
    co = countOverlaps(pos,thisgap); throwL1[which(co>0)]=TRUE
    sel = which(sctabR1$Chrom == s); pos = IRanges(start=sctabR1$Position[sel],width=1)
    co = countOverlaps(pos,thisgap); throwR1[which(co>0)]=TRUE
    if(!is.null(clusterdir0)){
        sel = which(sctabL0$Chrom == s); pos = IRanges(start=sctabL0$Position[sel],width=1)
        co = countOverlaps(pos,thisgap); throwL0[which(co>0)]=TRUE
        sel = which(sctabR0$Chrom == s); pos = IRanges(start=sctabR0$Position[sel],width=1)
        co = countOverlaps(pos,thisgap); throwR0[which(co>0)]=TRUE
    }
}
scL1=scL1[!throwL1]; scR1=scR1[!throwR1]; sctabL1=sctabL1[!throwL1,]; sctabR1=sctabR1[!throwR1,]
cat(length(scL1),"left clip clusters left for case after gap removal.\n")
cat(length(scR1),"right clip clusters left for case after gap removal.\n")
if(!is.null(clusterdir0)){
    scL0=scL0[!throwL0]; scR0=scR0[!throwR0]; sctabL0=sctabL0[!throwL0,]; sctabR0=sctabR0[!throwR0,]
    cat(length(scL0),"left clip clusters left for control after gap removal.\n")
    cat(length(scR0),"right clip clusters left for control after gap removal.\n")
}


### Remove clusters found in contrast.
if(!is.null(clusterdir0)){
    seqnames=unique(sctabL1$Chrom)
    throwL1=rep(FALSE,length(scL1));throwR1=rep(FALSE,length(scR1))
    for(s in seqnames){
        sel1 = which(sctabL1$Chrom==s)
        sel0 = which(sctabL0$Chrom==s)
        mm = match(sctabL1$Position[sel1],sctabL0$Position[sel0])
        throwL1[sel1[!is.na(mm)]]=TRUE
        sel1 = which(sctabR1$Chrom==s)
        sel0 = which(sctabR0$Chrom==s)
        mm = match(sctabR1$Position[sel1],sctabR0$Position[sel0])
        throwR1[sel1[!is.na(mm)]]=TRUE
    }
    scL1=scL1[!throwL1]; scR1=scR1[!throwR1]; sctabL1=sctabL1[!throwL1,]; sctabR1=sctabR1[!throwR1,]
    cat(length(scL1),"left clip clusters left for case after removing contrast.\n")
    cat(length(scR1),"right clip clusters left for case after removing contrast.\n")
}

###############################################
# Remap right-clipped clusters.
###############################################
sc1 = scR1
sctab1 = sctabR1
sctab1m = sctabL1
sc1m = scL1
idprefix="SCR"
mateidprefix="SCL"

trunksize=ceiling(length(sc1)/ntrunks) #here ntrunks=1
if(trunknum<ntrunks){
  trunkR = trunksize*(trunknum-1)+c(1:trunksize)
} else {
  if(length(sc1)==0) trunkR=c() else trunkR=c((trunksize*(trunknum-1)+1):length(sc1))
  #trunkR = c((trunksize*(trunknum-1)+1):length(sc1))
}

remap_resR=list()
if(length(trunkR)>0){
	cat("Processing trunkR ",trunkR[1]," to ",trunkR[length(trunkR)],".\n",sep="")
	if(debug) {
  	remap_resR=lapply(trunkR,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=rg,doLeftClip=FALSE,idprefix=idprefix,mateidprefix=mateidprefix)
	} else { #accelerate using mclapply
  remap_resR=mclapply(trunkR,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=rg,doLeftClip=FALSE,idprefix=idprefix,mateidprefix=mateidprefix)
	}
}
save(trunkR, remap_resR,file=rdatfile)

###############################################
# Remap left-clipped clusters.
###############################################
sc1=scL1
sctab1=sctabL1
sctab1m = sctabR1
sc1m = scR1
idprefix = "SCL"
mateidprefix="SCR"
trunksize=ceiling(length(sc1)/ntrunks) #here ntrunks=1
if(trunknum<ntrunks){
  trunkL = trunksize*(trunknum-1)+c(1:trunksize)
} else {
  if(length(sc1)==0) trunkL=c() else trunkL=c((trunksize*(trunknum-1)+1):length(sc1))
}
#trunksize=floor(length(sc1)/(ntrunks-1))
#if(trunknum<ntrunks){
#    trunkL = trunksize*(trunknum-1)+c(1:trunksize)
#} else {
#    trunkL = c((trunksize*(trunknum-1)+1):length(sc1))
#}
remap_resL=list()
if(length(trunkL)>0){
	cat("Processing trunkL ",trunkL[1]," to ",trunkL[length(trunkL)],".\n",sep="")
	if(debug) {
  	remap_resL=lapply(trunkL,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=rg,doLeftClip=TRUE,idprefix=idprefix,mateidprefix=mateidprefix)
	} else { #accelerate using mclapply
  	remap_resL=mclapply(trunkL,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=rg,doLeftClip=TRUE,idprefix=idprefix,mateidprefix=mateidprefix)
	}
}
save(remap_resR,remap_resL,scL1,scR1,trunkL,trunkR,file=rdatfile)
cat("=Info: done,",taggie(gtk),"\n")
cat("warnings if any:\n")
warnings()

