#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
#### Externalities ####
suppressMessages(library(Biobase))
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(parallel))
suppressMessages(library(Biostrings))
suppressMessages(library(swan))

#  tryCatch({
#    source("~/setup/swan/R/libswan.R")
#    source("~/setup/swan/R/svlib.r")
#    source("~/setup/swan/R/ThresholdFunctions.R")
#    dyn.load("~/setup/swan/src/swan.so")
#    cat("-Info: loaded from ~/setup/swan package!\n")
#  }, error=function(e) {
#    cat("-Warn: cant load from ~/setup/swan package!\n")
#    suppressMessages(library(swan))
#    cat("-Info: loaded from installed swan package!\n")
#  })

#MIN.READS.PER.CLUSTER=5
#MIN.BASES.PER.CLUSTER=50
#SC.MIN.CONSENSUS.LENGTH=30
SC.PROPORTION.MISMATCH.THRESH=0.05
MIN.GAP.DISTANCE=10000
SC.MIN.CONSENSUS.LENGTH=30
MAX.MATCHES=10
MAX.DIST.CONSENSUS.CLUSTER=100 # Maximum distance allowed between consensus and SC cluster in destination.
REMAP.WIN=1000 # window for remapping of destination consensus back to home.

#### CommandLines ####
option_list <- list(
  #all options cann't default to NA, no multiline help, -g not usable  !!!!!
  #options the same cross samples and read groups
  make_option(c("-c", "--chromosomeName"), default="a",
              help="chromosome to scan [default %default]"),
  make_option(c("-n", "--gapfile"), default="none",
              help="gap/N locations of hg19 in ucsc format [default %default]"),
  make_option(c("-i", "--minReadPerCluster"), default="3:5",
              help="minimal number of reads per cluster, [spY]:spX [default %default]"),
  make_option(c("-j", "--minBasePerCluster"), default="30:50",
              help="minimal number of total bases per cluster, [spY]:spX [default %default]"),
  make_option(c("-u", "--scanStart"), type="integer", default=1,
              help="1-indexed scan start, [default %default]"),
  make_option(c("-v", "--scanEnd"), type="integer", default=max_chr_len,
              help="1-indexed scan end, [default %default]"),
  make_option(c("-z", "--trunkSize"), type="integer", default=1000000,
              help="trunk size for scanning bamfile [default %default]"),
  make_option(c("-d", "--contdir"), default="none",
              help="contrast directory, [default %default]"),
  make_option(c("-r", "--sample"), default="spX,INFO,MIX,DESCRIPTION",
              help="mannual override of spX information [default %default]"),
  make_option(c("-t", "--stat"), default="",
              help="stat inputs: [spY.stat:]spX.stat;
                    [spY.stat:]spX.stat implicitly assumed, must explicit supplied for multi-lib bams \n"),
  make_option(c("-e", "--delthresh"), default=0.8,
              help="foldchange threshold for deletion events, [default %default]"),
  make_option(c("-k", "--dupthresh"), default=1.2,
             help="foldchange threshold for duplication events, [default %default]"),
  make_option(c("-m", "--maxfc"), type="integer", default=20000000,
             help="maximum region size for fold change check (due to memory considerations), [default %default]"),
  make_option(c("-l", "--minfc"), type="integer", default=10000,
             help="minimum region size for fold change check for del/dup calls, [default %default]"),
  make_option(c("-f", "--fineconf"), action="store_true", default=FALSE,
    help="fine conf mode and .bam is assumed for all inputs, see manual [default %default]"),
  make_option(c("-b", "--minGapPair"), type="integer", default=25,
             help="A breakpoint and its mate must be separated by at least this value, [default %default]"),
  make_option(c("-y", "--insparam"), default="0:0",
             help="parameters for calling insertions, 0 means to estimate from data [default %default]"),
  make_option(c("-x", "--hotspot"), default="10000:3",
             help="setting for hotspot filtering, [default %default]"),
  make_option(c("-g", "--gapdist"), type="integer", default=100000,
             help="setting for gaps (centromere or telomere) filtering, [default %default]"),
  make_option(c("-s", "--savedRData"), default="none",
             help="whehter there is saved RData to load [default %default]"),
  make_option(c("-o", "--spout"), default="input",
              help="sample output prefix, [default %default]"),
  make_option(c("-p", "--plot"), default="null",
             help="file for diagnostic plots, [default %default]"),
  make_option(c("--vcf"), action="store_true", default=FALSE,
             help="output VCF file, [default %default]"),
  make_option(c("--nobam"), action="store_true", default=FALSE,
             help="Use bam file for calling, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]")
)

parser=OptionParser(usage="%prog [options] ref_file [spY.rg1.bam,spY.rg2.bam]:spX.rg1.bam,spX.rg2.bam", option_list=option_list)
cat("-Info: sclip_scan vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args=commandArgs(trailingOnly=TRUE)
cmd=parse_args(parser,args,print_help_and_exit=TRUE,positional_arguments=TRUE)
if(length(cmd$args)!=2){ print_help(parser); quit(); }
verbose=cmd$options$noQuiet #true to turn on
debug=cmd$options$debug

if(debug){ options(warn=2) } else { options(warn=-1) }
gtk=proc.time()[3]
gmk=get_gmk(Sys.getpid())

fineconf=cmd$options$fineconf
inputs=cmd$args[2]; ref_file=cmd$args[1]; sample_tags=strsplit(inputs,split=':')[[1]]
n_sample=length(sample_tags); rg_files=strsplit(sample_tags,split=","); n_sp=n_sample
outprefix=cmd$options$spout
if(outprefix=='input') {
  rg_prefix=rg_files
  rg_prefix=lapply(rg_prefix,function(x) { gsub(".bam$","",x) })
  sp_prefix=gsub("\\.\\w+$","",sapply(rg_files,lcPrefix),perl=T)
  outprefix=sp_prefix[[n_sp]]
} 
text_stat=cmd$options$stat
stat_files=lapply(rg_files,function(x){ gsub(".bam$",".stat",x) }) #here, we need .stat file for coverage, rl
if(text_stat!=""){
  stat_files=as.list(strsplit(text_stat,split=':')[[1]])
}
nstat=sapply(stat_files,length)
if(max(nstat)>1){
  stop("Error Input: must provide stat file by -t if multiple lib bams \n")
}
statfile=stat_files[[n_sp]]; constatfile=NULL;
if(n_sp>1) { constat=stat_files[[1]] }
trunk_size=cmd$options$trunkSize  
seq_name=cmd$options$chromosomeName
MIN.READS.PER.CLUSTER=as.integer(strsplit(cmd$options$minReadPerCluster,":")[[1]])
MIN.BASES.PER.CLUSTER=as.integer(strsplit(cmd$options$minBasePerCluster,":")[[1]])
#when default is given for one sample analysis, use the second as real options
if(n_sp==1 && length(MIN.READS.PER.CLUSTER)>1) { MIN.READS.PER.CLUSTER=MIN.READS.PER.CLUSTER[2] }
if(n_sp==1 && length(MIN.BASES.PER.CLUSTER)>1) { MIN.BASES.PER.CLUSTER=MIN.BASES.PER.CLUSTER[2] }
scan_start=max(1,cmd$options$scanStart)
scan_end=cmd$options$scanEnd
gap_file=cmd$options$gapfile
gap=read_gap(gap_file)
ref=FaFile(ref_file)
ref_seq=scanFa(ref, param=scanFaIndex(ref)) #this returns a DNAStringSet #rg = .Call("scanFa", ref_file, PACKAGE="swan")

###################

if(seq_name %in% c("A","a","ALL","all","All")){
  seq_name = names(ref_seq)[seq_len(min(length(ref_seq),24))]
  cat("all chrs:",seq_name,"\n")
  ref_len = sapply(seq_name, function(x) { length(ref_seq[[x]]) })
  start = rep(1,length(ref_len))
  end = as.vector(ref_len)
  chrname = NULL
} else {
  ref_len = length(ref_seq[[seq_name]]) #[]->DNAStringSet, support width; [[]]->DNAString, support length
  start = max(1,scan_start); end=min(ref_len,scan_end)
  chrname = seq_name
}

sample_id=strsplit(cmd$options$sample,split=",")[[1]][1]
sample_info=matrix(strsplit(cmd$options$sample,split=",")[[1]],ncol=4)
contig_info=matrix(c("seqname","ref_len","md5_of_seq","human"),ncol=4)
outputVCF=if(cmd$options$vcf) paste(outprefix,"sclip","vcf",sep=".") else NULL
noBam=cmd$options$nobam
gapdist=cmd$options$gapdist
plotpf=ifelse(cmd$options$plot=="null",paste(outprefix,"sclip",sep="."),cmd$options$plot)

bamfile1=rg_files[[n_sp]]; bamfile0=NULL
sclip_opt=list(); mean_cvg1=list(); rl1=rep(NA,length(bamfile1));
cvg1=rep(NA,length(bamfile1));

for(ix in seq_along(ref_seq)){
  sn=names(ref_seq)[ix]; mean_cvg1[[sn]]=0; idx=1;
  for(bam_file in bamfile1){
    cat("-Info: calculating coverage...",bam_file,",",sn,",")
    #this is assumed that genome are evenly covered
    #it is user's responsibility to provide a .stat file by running swan_stat on the bam input file
    rl1[idx]=read.table(stat_files[[n_sp]], header=T)$rl[1]
    mean_cvg_tmp=sum(read.table(stat_files[[n_sp]], header=T)$cvg)
    if(fineconf) { # actually load bam file by fine mode
      nreads=count_bam(bam_file,sn);
      mean_cvg_tmp=nreads*rl1[idx]/ref_len[ix]
    }
    cat(mean_cvg_tmp,"\n")
    idx=idx+1
    mean_cvg1[[sn]]=mean_cvg1[[sn]]+mean_cvg_tmp
  }
}
if(n_sp>1){#
  bamfile0=rg_files[[1]]; mean_cvg0=list();  rl0=rep(NA,length(bamfile0)); cvg0=rep(NA,length(bamfile0));
  for(ix in seq_along(ref_seq)){
    sn=names(ref_seq)[ix]; mean_cvg0[[sn]]=0; idx=1;
    for(bam_file in bamfile0){
      cat("-Info: calculating coverage...",bam_file,",",sn,",")
      rl0[idx]=read.table(stat_files[[1]], header=T)$rl[1]
      mean_cvg_tmp=sum(read.table(stat_files[[1]], header=T)$cvg)
      if(fineconf) { # actually load bam file by fine mode
         nreads=count_bam(bam_file,sn);
         mean_cvg_tmp=nreads*rl0[idx]/ref_len[ix]
      }
      cat(mean_cvg_tmp,"\n")
      idx=idx+1
      mean_cvg0[[sn]]=mean_cvg0[[sn]]+mean_cvg_tmp
    }
  }
}
coverage1=mean(unlist(mean_cvg1),na.rm=T)
rl=median(rl1,na.rm=T)
scalefac=1
if(n_sp>1){
 coverage0=mean(unlist(mean_cvg0),na.rm=T)
 scalefac=coverage1/coverage0
}
### if we got saved RData, use that ###
savedRData=cmd$options$savedRData
fromRData=FALSE
reg_sclip=list()
if(savedRData!="none") {
  sclip_load=new.env()
  tryCatch( {load(savedRData,envir=sclip_load)}, error=function(e){cat("\n==Warn:",savedRData,"empty, nonexist or errornous formatted\n")})
  bndtab = sclip_load[["bndtab"]]
  scL1 = sclip_load[["scL1"]]
  scR1 = sclip_load[["scR1"]]
  resL = sclip_load[["resL"]]
  resR = sclip_load[["resR"]]
  fromRData = TRUE
} else {
### show what we set for scan ###

cat("-Info: rg_files\n"); print(rg_files)
cat("MIN.READS.PER.CLUSTER=",MIN.READS.PER.CLUSTER,"\n")
cat("MIN.BASES.PER.CLUSTER=",MIN.BASES.PER.CLUSTER,"\n")

#######################
### sclip_scan     ####
#######################

### Suko's Work which I still need to understand ###
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

scan_chr=function(sn,start,end,files,mrpc,mbpc){
 
  cat("sn=",sn,"files=",files,"\n")
  n_trunks=ceiling((end-start+1)/trunk_size);
  cat("scan_start=",nsf(start),"scan_end=",nsf(end),"n_trunks=",n_trunks,"\n")
  
  rcpp_core_sclip = function(ti) {
    .Call("core_sclip", ti, n_trunks, start, end, trunk_size, files, sn, gap, mrpc, mbpc, SC.PROPORTION.MISMATCH.THRESH, MIN.GAP.DISTANCE, PACKAGE="swan");
  }
  
  if(!debug) { #mcapply of small memories are assumed safe and faster then lapply
    sclip_res=mclapply(seq(n_trunks), rcpp_core_sclip)
  } else { #accelerate using mclapply
    sclip_res=mclapply(seq(n_trunks), rcpp_core_sclip)
  }
  
  #paradimg to combine data from list(list(name1, name2),list(name1,name2)) to list(name1,name2)
  var_names=ls(sclip_res[[1]]) #names
  var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },sclip_res[[1]])
  all_RData=lapply(seq_along(var_names),new_merge_RData,var_names,var_classes,sclip_res)
  for(di in seq_along(all_RData)) { all_RData[[di]]=all_RData[[di]][[1]] }
  names(all_RData)=var_names
  cat("scL","scR","\n")
  cat(length(all_RData$scL),length(all_RData$scR),"\n")
 #done the paradigm and check if number is correct
  
  #add in chrom information; sc1L and sc1R
  for(i in seq_along(all_RData$scL)){ all_RData$scL[[i]]$chrom=sn }   
  for(i in seq_along(all_RData$scR)){ all_RData$scR[[i]]$chrom=sn }
  
  all_RData
}

scan_RData=list() #list is not in-place mutable
for(i in seq_len(n_sp)){
 scan_tmp=lapply(seq_along(seq_name),function(x){ scan_chr(seq_name[x],start[x],end[x],rg_files[[i]],MIN.READS.PER.CLUSTER[i],MIN.BASES.PER.CLUSTER[i]) })
 #paradimg to combine data from list(list(name1, name2),list(name1,name2)) to list(name1,name2)
  var_names=ls(scan_tmp[[1]]) #names of c(scL,scR) 
  var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },scan_tmp[[1]])
  tmp_data=lapply(seq_along(var_names),new_merge_RData,var_names,var_classes,scan_tmp)
  for(di in seq_along(tmp_data)) { tmp_data[[di]]=tmp_data[[di]][[1]] }
 names(tmp_data)=var_names
 #done the paradigm and check if number is correct
 scan_RData[[i]]=tmp_data
}

#######################
### sclip_call     ####
#######################

###############################################
# Combine sclip clusters over all chromosomes.
# Load all bam files in directory.
###############################################
### Case directory
scL1 = scan_RData[[n_sp]]$scL
scR1 = scan_RData[[n_sp]]$scR 
cat(length(scL1),"left clip clusters loaded for case.\n")
cat(length(scR1),"right clip clusters loaded for case.\n")
sctabL1 = summarizeClusters(scL1,includeNBases=TRUE)
sctabR1 = summarizeClusters(scR1,includeNBases=TRUE)
### Control directory
if(n_sp>1){
 scL0 = scan_RData[[1]]$scL
 scR0 = scan_RData[[1]]$scR 
  cat(length(scL0),"left clip clusters loaded for control.\n")
  cat(length(scR0),"right clip clusters loaded for control.\n")
  sctabL0 = summarizeClusters(scL0,includeNBases=TRUE)
  sctabR0 = summarizeClusters(scR0,includeNBases=TRUE)
}

###############################################
# Filtering
###############################################
### Remove clusters found in contrast.
if(n_sp>1){
  gapchrnames=unique(sctabL1$Chrom)
  throwL1=rep(FALSE,length(scL1));throwR1=rep(FALSE,length(scR1))
  for(s in gapchrnames){
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
### Remove centromere and gap regions.
gap=read_gap(gap_file)
gapchrnames = unique(gap$chrom)
throwL1 = rep(FALSE,length(scL1)); throwR1 = rep(FALSE,length(scR1))
throwL0 = rep(FALSE,length(scL1)); throwR0 = rep(FALSE,length(scR1))
for(s in gapchrnames){
  sel2 = which(gap$chrom==s)
  thisgap = IRanges(start=gap$chromStart[sel2],end=gap$chromEnd[sel2])
  sel = which(sctabL1$Chrom == s); pos = IRanges(start=sctabL1$Position[sel],width=1)
  co = countOverlaps(pos,thisgap); throwL1[which(co>0)]=TRUE
  sel = which(sctabR1$Chrom == s); pos = IRanges(start=sctabR1$Position[sel],width=1)
  co = countOverlaps(pos,thisgap); throwR1[which(co>0)]=TRUE
  if(n_sp>1){
    sel = which(sctabL0$Chrom == s); pos = IRanges(start=sctabL0$Position[sel],width=1)
    co = countOverlaps(pos,thisgap); throwL0[which(co>0)]=TRUE
    sel = which(sctabR0$Chrom == s); pos = IRanges(start=sctabR0$Position[sel],width=1)
    co = countOverlaps(pos,thisgap); throwR0[which(co>0)]=TRUE
  }
}
scL1=scL1[!throwL1]; scR1=scR1[!throwR1]; sctabL1=sctabL1[!throwL1,]; sctabR1=sctabR1[!throwR1,]
cat(length(scL1),"left clip clusters left for case after gap removal.\n")
cat(length(scR1),"right clip clusters left for case after gap removal.\n")
if(n_sp>1){
  scL0=scL0[!throwL0]; scR0=scR0[!throwR0]; sctabL0=sctabL0[!throwL0,]; sctabR0=sctabR0[!throwR0,]
  cat(length(scL0),"left clip clusters left for control after gap removal.\n")
  cat(length(scR0),"right clip clusters left for control after gap removal.\n")
}

# Remap right-clipped clusters.
###############################################
sc1 = scR1
sctab1 = sctabR1
sctab1m = sctabL1
sc1m = scL1
idprefix="SCR"
mateidprefix="SCL"
ntrunks=1
trunknum=1

trunksize=ceiling(length(sc1)/ntrunks) #here ntrunks=1
if(trunknum<ntrunks){
  trunkR = trunksize*(trunknum-1)+c(seq_len(trunksize))
} else {
  if(length(sc1)==0) trunkR=c() else trunkR=c((trunksize*(trunknum-1)+1):length(sc1))
  #trunkR = c((trunksize*(trunknum-1)+1):length(sc1))
}

remap_resR=list()
if(length(trunkR)>0){
  if(debug) cat("Processing trunkR ",trunkR[1]," to ",trunkR[length(trunkR)],".\n",sep="")
  if(!debug) {
    remap_resR=mclapply(trunkR,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=ref_seq,doLeftClip=FALSE,idprefix=idprefix,mateidprefix=mateidprefix)
  } else { #accelerate using mclapply
    remap_resR=lapply(trunkR,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=ref_seq,doLeftClip=FALSE,idprefix=idprefix,mateidprefix=mateidprefix)
  }
}

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
  trunkL = trunksize*(trunknum-1)+c(seq_len(trunksize))
} else {
  if(length(sc1)==0) trunkL=c() else trunkL=c((trunksize*(trunknum-1)+1):length(sc1))
}
remap_resL=list()
if(length(trunkL)>0){
  if(debug) cat("Processing trunkL ",trunkL[1]," to ",trunkL[length(trunkL)],".\n",sep="")
  if(!debug) {
    remap_resL=mclapply(trunkL,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=ref_seq,doLeftClip=TRUE,idprefix=idprefix,mateidprefix=mateidprefix)
  } else { #accelerate using mclapply
    remap_resL=lapply(trunkL,remapSoftclipCluster,sc=sc1,sctab=sctab1,scm=sc1m,sctabm=sctab1m,templateseq=ref_seq,doLeftClip=TRUE,idprefix=idprefix,mateidprefix=mateidprefix)
  }
}
#remap_resR, remap_resL, scL1, scR1, trunkL, trunkR

######################
### sclip_events.R ###
######################

insparam= as.numeric(strsplit(cmd$options$insparam,":")[[1]])
MAX.GAP.INSERTION=insparam[1]
MIN.TOTAL.READS.INS=insparam[2]
if(MAX.GAP.INSERTION==0) MAX.GAP.INSERTION=NA
if(MIN.TOTAL.READS.INS==0) MIN.TOTAL.READS.INS=NA
MIN.TOTAL.READS.INS=qbinom(0.1, qpois(0.2,coverage1), 0.5)
MIN.TOTAL.READS.INS=max(MIN.TOTAL.READS.INS, 6)

#if(np){
#  stattab=read.table(statfile,header=TRUE,sep="\t")
#  rl=stattab$rl
#  if(is.na(MIN.TOTAL.READS.INS)){
#    coverage=NULL
#    if(!is.null(stattab$cvg)) coverage=stattab$cvg
#    if(!is.null(stattab$coverage)) coverage=stattab$coverage
#    if(!is.null(coverage)){
#    }
#  }
#} else {
#  rl=100
#}

MIN.GAP.BETWEEN.BNDPAIR = cmd$options$minGapPair
MAX.FOLDCHANGE.DEL = cmd$options$delthresh
MIN.READS.DELDUP.CALL = 30
MIN.FOLDCHANGE.DUP = cmd$options$dupthresh
hotparam= as.numeric(strsplit(cmd$options$hotspot,":")[[1]])
DIST.TO.NEXT.IN.HOTSPOT=hotparam[1]
HOTSPOT.CLUSTER.SIZE=hotparam[2]
MAX.FC.CHECK=cmd$options$maxfc
MIN.FC.CHECK=cmd$options$minfc

#if(!noBam && is.null(bamfile1)){
#  cat("Bam file(s) need to be specified unless running in noBam option.\n")
#  quit()
#}
#if(!is.null(bamfile0) && statfile=="none"){
#  cat("Stat file needs to be specified if contrast bam given (use swan_stat.R).\n")
#  quit()
#}
#if(!is.null(bamfile0) && constatfile=="none"){
#  cat("Contrast stat file needs to be specified if contrast bam given (use swan_stat.R).\n")
#  quit()
#}

if(!is.null(bamfile1)) bamfile1=strsplit(bamfile1,split=",")[[1]]
if(!is.null(bamfile0)) bamfile0=strsplit(bamfile0,split=",")[[1]]

#if(statfile!="none"){
#    stattab=read.table(statfile,header=TRUE,sep="\t")
#    rl=stattab$rl
#} else {
#    rl=100
#}

MIN.GAP.BETWEEN.CALLS = rl
MAX.DISTANCE.SOFTCLIP.PAIR = rl

if(!is.null(statfile) && !is.null(constatfile)){
  stattab=read.table(statfile,header=TRUE,sep="\t")
  constattab=read.table(constatfile,header=TRUE,sep="\t")
  scalefac=stattab$cvg/constattab$cvg
  if(is.null(scalefac)) scalefac=1
} else {
  scalefac=1 
}
cat("scalefac=",scalefac,"\n")

resR = remap_resR
resL = remap_resL
#files=dir(calldir)
#sel=grep(".RData",files,fixed=TRUE)
#load(file.path(calldir,files[sel[1]]))
#cat("Loading data from ",file.path(calldir,files[sel[1]]),"\n")
#resL=vector("list",length(scL1))
#resR=vector("list",length(scR1))
#resL[trunkL] = remap_resL
#resR[trunkR] = remap_resR
#if(length(sel)>=2){ #only needed if >2 sscan files
#  for(i in 2:length(sel)){
#    cat("Loading data from ",file.path(calldir,files[sel[i]]),"\n")
#    load(file.path(calldir,files[sel[i]]))
#    resR[trunkR] = remap_resR
#    resL[trunkL] = remap_resL
#  }
#}

############################################## 
# Organize bnds and call events.
############################################## 
# ------- Step 1: Make a data frame of all bnd entries: chr1,pos1,goRight1,chr2,pos2,goRight2
chr1=rep("",0); pos1=rep(0,0);goRight1=rep(0,0);chr2=rep("",0); pos2=rep(0,0);goRight2=rep(0,0);
cid=rep(0,0);mcid=rep(0,0)
for(i in seq_along(resR)){
    if(length(resR[[i]]$bndlist)>0){
        for(j in seq_along(resR[[i]]$bndlist)){
            bndl=resR[[i]]$bndlist[[j]];
            chr1=c(chr1,bndl$chr1);pos1=c(pos1,bndl$pos1);goRight1=c(goRight1,bndl$goRight1);
            chr2=c(chr2,bndl$chr2);pos2=c(pos2,bndl$pos2);goRight2=c(goRight2,bndl$goRight2);
            cid=c(cid,bndl$cid);mcid=c(mcid,bndl$mcid)            
        }
    }
}
cat(length(chr1)," entries found using right soft-clip clusters.\n",sep="")
for(i in seq_along(resL)){
    if(length(resL[[i]]$bndlist)>0){
        for(j in seq_along(resL[[i]]$bndlist)){
            bndl=resL[[i]]$bndlist[[j]];
            chr1=c(chr1,bndl$chr1);pos1=c(pos1,bndl$pos1);goRight1=c(goRight1,bndl$goRight1);
            chr2=c(chr2,bndl$chr2);pos2=c(pos2,bndl$pos2);goRight2=c(goRight2,bndl$goRight2);
            cid=c(cid,bndl$cid);mcid=c(mcid,bndl$mcid)            
        }
    }
}
bndtab=data.frame(chr1,pos1,goRight1,chr2,pos2,goRight2,cid,mcid)
bndtab$chr1 = paste(bndtab$chr1)
bndtab$chr2 = paste(bndtab$chr2)
cat(nrow(bndtab)," total entries found.\n",sep="")
#if seq_name==all, chrname=NULL; else chrname=seq_name 
if(!is.null(chrname)){
    sel=which(bndtab$chr1==chrname | bndtab$chr2==chrname)
    bndtab=bndtab[sel,]
    cat(nrow(bndtab)," total entries found mapping to target ",chrname,"\n",sep="")
} else {
    cat(nrow(bndtab)," total entries found mapping to genome.\n",sep="")
}

# ------- Step 2a: Sort, each breakend has smaller (chr,pos), followed by larger (chr,pos)
cat("Sorting each breakend pair so that (chr1,pos1) << (chr2,pos2) ...\n")
print(nrow(bndtab))
for(i in seq_len(nrow(bndtab))){
    if(i %% 100 == 0) cat("\t\tDone with ",i," of ",nrow(bndtab),"\n")
    if(!posIncreasing(bndtab$chr1[i],bndtab$pos1[i],bndtab$chr2[i],bndtab$pos2[i])){
        chr2 = bndtab$chr1[i]; pos2 = bndtab$pos1[i]; goRight2 = bndtab$goRight1[i]; mcid=bndtab$cid[i]
        bndtab$chr1[i]=bndtab$chr2[i]; bndtab$pos1[i] = bndtab$pos2[i]; 
        bndtab$goRight1[i]=bndtab$goRight2[i]; bndtab$cid[i]=bndtab$mcid[i]
        bndtab$chr2[i]=chr2; bndtab$pos2[i]=pos2; 
        bndtab$goRight2[i]=goRight2; bndtab$mcid[i]=mcid
    }
}

if(nrow(bndtab)>0){
  # ------- Step 2b: Sort, the rows of the table should increase in (chr1,pos1)
  cat("Sorting rows of table so that breakend pairs are increasing...\n")
  o=order(bndtab$chr1, bndtab$pos1)
  bndtab=bndtab[o,]

  if(debug){
    plotFile=paste(plotpf,"fig1","pdf",sep=".")
    cat("saving diagnostic plots to",plotFile,"\n")
    pdf(plotFile)
    plot(bndtab$pos1[!is.na(bndtab$pos1)],pch=seq_along(sort(unique(bndtab$chr1))))
    dev.off()
  }

  # ------- Step 3: Various filtering

  bndtab.old=bndtab
  # ------- Step 3a: Hot spot filtering.
  # filter hotspots based on extreme clustering of soft-clip clusters.
  cat("Filtering out hotspots where soft-clip breakends cluster.\nClusters defined as ",HOTSPOT.CLUSTER.SIZE," breakends with less than ",DIST.TO.NEXT.IN.HOTSPOT," between adjacent entries.\n",sep="")
  # only apply this when nrow(bndtab)>1
  if(nrow(bndtab)>1){
    toNext=rep(NA,nrow(bndtab)-1)
    for(i in seq_len(nrow(bndtab)-1)){
      toNext[i] = bndtab$pos1[i+1]-bndtab$pos1[i]
    }
    runst=0; runend=-1
    throw=rep(FALSE,nrow(bndtab))
    for(i in seq_len(nrow(bndtab)-1)){
      if(toNext[i]<DIST.TO.NEXT.IN.HOTSPOT && toNext[i]>=0){
        #cat(i,"toNext=",toNext[i],": part of run, st=",runst,", end=",runend,"\n")
        if(runst==0) runst=i;
        runend=i
      } else {
        #cat(i,"toNext=",toNext[i],": not in run, st=",runst,", end=",runend,"\n")
        if(runend-runst>HOTSPOT.CLUSTER.SIZE) throw[runst:runend]=TRUE
        runst=0; runend=-1
      }
    }
    if(debug){
      plotFile=paste(plotpf,"fig2","pdf",sep=".")
      pdf(plotFile)
      par(mfrow=c(2,1))
      plot(bndtab$pos1[!is.na(bndtab$pos1)],pch=seq_along(sort(unique(bndtab$chr1))),main="Distribution of breakpoints",ylab="pos1")
      plot(pmax(pmin(c(toNext,NA),DIST.TO.NEXT.IN.HOTSPOT),0), col=throw+1, main="Distance to next bnd",ylab=paste("pmax(pmin(Dist,",DIST.TO.NEXT.IN.HOTSPOT,"),0)",sep=""),pch=seq_along(sort(unique(bndtab$chr1))))
      dev.off()
    }
    cat(sum(throw)," bnd entries thrown due to being within hotspots.\n",sep="")
    bndtab = bndtab[!throw,]
  }

  # ------- Step 3b: Get rid of duplicates in table.  
  if(nrow(bndtab)>1){
    cat("Throwing away duplicates...\n")
    isDuplicate<-function(entree1,entree2,DELTA){
      if(entree1$chr1==entree2$chr1 && entree1$chr2==entree2$chr2){
        if(entree1$goRight1==entree2$goRight1 && entree1$goRight2==entree2$goRight2){
            if(abs(entree1$pos1-entree2$pos1)<DELTA && abs(entree1$pos2-entree2$pos2)<DELTA){
                return(TRUE)
            }
        }
      } 
      return(FALSE)
    }

    throw=rep(FALSE,nrow(bndtab))
    for(j in 2:nrow(bndtab)){
      if(j %% 10000==0) cat("\t\tDone with ",j," rows out of ",nrow(bndtab),".\n") 
      for(i in (j-1):1){
        if (bndtab[i,]$chr1 != bndtab[j,]$chr1 || abs(bndtab[i,]$pos1 - bndtab[j,]$pos1) >= MIN.GAP.BETWEEN.CALLS) {
          break  # We can break early because we sort the table by (chr1,pos1).
        }
        if (isDuplicate(bndtab[i,],bndtab[j,],MIN.GAP.BETWEEN.CALLS)) {
          #cat("Found Duplicate!\n")
          throw[j]=TRUE
          break
        }
      }
    }
    cat(sum(throw)," bnd entries thrown due to being within ",MIN.GAP.BETWEEN.CALLS," of other calls.\n",sep="")
    bndtab = bndtab[!throw,]
    if(debug) {
      plotFile=paste(plotpf,"fig4","pdf",sep=".")
      pdf(plotFile)
      plot(bndtab$pos1[!is.na(bndtab$pos1)],pch=seq_along(sort(unique(bndtab$chr1))))
      dev.off()
    }
  }

  # -------- Step 3c: Get rid of entries in table where pos1-pos2<MIN.GAP.BETWEEN.CALLS.
  cat("Throwing away breakend pairs separated by less than ",MIN.GAP.BETWEEN.BNDPAIR," ...\n")
  if(nrow(bndtab)>1){
    throw=rep(FALSE,nrow(bndtab))
    for(i in seq_len(nrow(bndtab))){
      if(bndtab$chr1[i]==bndtab$chr2[i] && bndtab$pos2[i]-bndtab$pos1[i]<MIN.GAP.BETWEEN.BNDPAIR){
        throw[i]=TRUE
      }
    }
    cat(sum(throw)," bnd entries thrown due to bnd pair being within ",MIN.GAP.BETWEEN.BNDPAIR,".\n",sep="")
    bndtab =bndtab[!throw,]
    if(debug) {
      plotFile=paste(plotpf,"fig5","pdf",sep=".")
      pdf(plotFile)
      plot(bndtab$pos1[!is.na(bndtab$pos1)],pch=seq_along(sort(unique(bndtab$chr1))))
      dev.off()
    }
  }

  # ------- Step 4: Create extra columns EventType, PartnerBnds=other bnds that participate.
  EventType = rep(NA,nrow(bndtab))
  EventStart = rep(NA,nrow(bndtab)) # In case of TRA,DUP.INS, EventStart, EventEnd are site of insertion.
  EventEnd = rep(NA,nrow(bndtab))
  Distance=rep(NA,nrow(bndtab))
  DonorStart=rep(NA,nrow(bndtab))
  DonorEnd=rep(NA,nrow(bndtab))
  DonorChr=rep(NA,nrow(bndtab))
  EventChr=rep(NA,nrow(bndtab))
  PartnerBnds=rep("",nrow(bndtab))
  FoldChange=rep(NA,nrow(bndtab))

  # ------- Step 5: Treat double-double clusters: 
  #         pairs of breakends (A,Am), (B,Bm) 
  #         where goRightA=FALSE, goRightB=TRUE 
  #         where posB-posA < DELTA.
  isPair<-function(chrA,posA,goRightA,chrB,posB,goRightB,DELTA){
    if(chrA==chrB){
        if(goRightA && !goRightB){
            d=posA-posB
            if(d>=0 && d<DELTA) return(TRUE)
        } 
        if(!goRightA && goRightB){
            d=posB-posA
            if(d>=0 && d<DELTA) return(TRUE)
        }
    } 
    return(FALSE)
  }
  nbnd = nrow(bndtab)
  bndPair=matrix(nrow=0,ncol=4)
  if(nbnd>1){
    for(i in seq_len(nbnd-1)){
    if(i %% 10000==0) cat("\t\tDone with ",i," out of ",nbnd,".\n") 
   
        for(j in (i+1):nbnd){
            pair11 = isPair(bndtab$chr1[i],bndtab$pos1[i],bndtab$goRight1[i],bndtab$chr1[j],bndtab$pos1[j],bndtab$goRight1[j],DELTA=MIN.GAP.BETWEEN.CALLS)
            pair12 = isPair(bndtab$chr1[i],bndtab$pos1[i],bndtab$goRight1[i],bndtab$chr2[j],bndtab$pos2[j],bndtab$goRight2[j],DELTA=MIN.GAP.BETWEEN.CALLS)
            pair21 = isPair(bndtab$chr2[i],bndtab$pos2[i],bndtab$goRight2[i],bndtab$chr1[j],bndtab$pos1[j],bndtab$goRight1[j],DELTA=MIN.GAP.BETWEEN.CALLS)
            pair22 = isPair(bndtab$chr2[i],bndtab$pos2[i],bndtab$goRight2[i],bndtab$chr2[j],bndtab$pos2[j],bndtab$goRight2[j],DELTA=MIN.GAP.BETWEEN.CALLS)
            if(pair11) bndPair=rbind(bndPair,c(i,1,j,1))
            if(pair12) bndPair=rbind(bndPair,c(i,1,j,2))
            if(pair21) bndPair=rbind(bndPair,c(i,2,j,1))
            if(pair22) bndPair=rbind(bndPair,c(i,2,j,2))
        }
    }
  }
  # Go through breakend pairs, and classify.
  cat("There are ",nrow(bndPair)," (A,Am);(B,Bm) pairs where goRightA=FALSE, goRightB=TRUE, and posB-posA<",MIN.GAP.BETWEEN.CALLS,".\n")
  if(nrow(bndPair)>=1){
    for(i in seq_len(nrow(bndPair))){
        # Organize breakend pairs into the following form:
        # BreakendA=(A,Am), BreakendB=(B,Bm)
        # posA and posB are adjacent.
        # By definition, goRightA=FALSE, goRightB=TRUE
        if(bndPair[i,2]==1){
            chrA=bndtab$chr1[bndPair[i,1]];chrAm=bndtab$chr2[bndPair[i,1]]
            posA=bndtab$pos1[bndPair[i,1]];posAm=bndtab$pos2[bndPair[i,1]]
            goRightA=bndtab$goRight1[bndPair[i,1]];goRightAm=bndtab$goRight2[bndPair[i,1]]
        } else {
            chrA=bndtab$chr2[bndPair[i,1]];chrAm=bndtab$chr1[bndPair[i,1]]
            posA=bndtab$pos2[bndPair[i,1]];posAm=bndtab$pos1[bndPair[i,1]]
            goRightA=bndtab$goRight2[bndPair[i,1]];goRightAm=bndtab$goRight1[bndPair[i,1]]    
        }
        if(bndPair[i,4]==1){
            chrB=bndtab$chr1[bndPair[i,3]];chrBm=bndtab$chr2[bndPair[i,3]]
            posB=bndtab$pos1[bndPair[i,3]];posBm=bndtab$pos2[bndPair[i,3]]
            goRightB=bndtab$goRight1[bndPair[i,3]];goRightBm=bndtab$goRight2[bndPair[i,3]]
        } else {
            chrB=bndtab$chr2[bndPair[i,3]];chrBm=bndtab$chr1[bndPair[i,3]]
            posB=bndtab$pos2[bndPair[i,3]];posBm=bndtab$pos1[bndPair[i,3]]
            goRightB=bndtab$goRight2[bndPair[i,3]];goRightBm=bndtab$goRight1[bndPair[i,3]]    
        }
        if(goRightA && !goRightB){
            # Switch A and B.
            ochrA=chrA;ochrAm=chrAm;oposA=posA;oposAm=posAm;ogoRightA=goRightA;ogoRightAm=goRightAm
            chrA=chrB;chrAm=chrBm;posA=posB;posAm=posBm;goRightA=goRightB;goRightAm=goRightBm
            chrB=ochrA;chrBm=ochrAm;posB=oposA;posBm=oposAm;goRightB=ogoRightA;goRightBm=ogoRightAm
        } 
        # Is this an inversion?
        # Mates of A,B are also adjacent, goRightAm=FALSE,goRightBm=TRUE.
        if(!goRightAm && goRightBm && chrA==chrAm && chrB==chrBm && abs(posBm-posAm)<MAX.DISTANCE.SOFTCLIP.PAIR){
            EventType[bndPair[i,1]]="INV"; EventType[bndPair[i,3]]="INV"
            EventChr[bndPair[i,1]]=chrA; EventStart[bndPair[i,1]]=bndtab$pos1[bndPair[i,1]]; EventEnd[bndPair[i,1]]=bndtab$pos2[bndPair[i,1]]
            EventChr[bndPair[i,3]]=chrA; EventStart[bndPair[i,3]]=bndtab$pos1[bndPair[i,3]]; EventEnd[bndPair[i,3]]=bndtab$pos2[bndPair[i,3]]
            PartnerBnds[bndPair[i,1]] = bndPair[i,3]; PartnerBnds[bndPair[i,3]]=bndPair[i,1]
            cat("Found Inversion on target ",chrA,".\n",sep="")
        } else {
            # Otherwise, do mates of A,B bracket region on chromosome?
            if(chrAm==chrBm){
                res=NULL
                endtype=""
                if(goRightAm && !goRightBm && posAm<posBm){
                    # donor sequence not inverted.
                    if(!noBam && posBm-posAm < MAX.FC.CHECK) res=getCoverageFoldchange(chrAm,posAm,posBm,scalefac=scalefac,bamfile1,bamfile0) 
                    endtype=""
                    EventChr[bndPair[i,1]]=chrA; EventChr[bndPair[i,3]]=chrA;
                    EventStart[bndPair[i,1]]=posA; EventEnd[bndPair[i,1]]=posB;
                    EventStart[bndPair[i,3]]=posA; EventEnd[bndPair[i,3]]=posB;
                    DonorStart[bndPair[i,1]]=posAm; DonorEnd[bndPair[i,1]]=posBm;
                    DonorStart[bndPair[i,3]]=posAm; DonorEnd[bndPair[i,3]]=posBm;
                }
                if(!goRightAm && goRightBm && posAm>posBm){
                    # donor sequence inverted.
                    if(!noBam && posBm-posAm < MAX.FC.CHECK) res=getCoverageFoldchange(chrAm,posBm,posAm,scalefac=scalefac,bamfile1,bamfile0) 
                    endtype=".INV"
                    EventChr[bndPair[i,1]]=chrA; EventChr[bndPair[i,3]]=chrA;
                    EventStart[bndPair[i,1]]=posA; EventEnd[bndPair[i,1]]=posB;
                    EventStart[bndPair[i,3]]=posA; EventEnd[bndPair[i,3]]=posB;
                    DonorStart[bndPair[i,1]]=posBm; DonorEnd[bndPair[i,1]]=posAm;
                    DonorStart[bndPair[i,3]]=posBm; DonorEnd[bndPair[i,3]]=posAm;
                }
                PartnerBnds[bndPair[i,1]] = bndPair[i,3]
                PartnerBnds[bndPair[i,3]] = bndPair[i,1]

                if(!is.null(res)){
                    #print(res$fc)
                    if(length(res$fc)==0) res$fc=1 #To Nancy: why res$fc sometimes return 0 length
                    if(res$fc>MIN.FOLDCHANGE.DUP && res$contrastcov>MIN.READS.DELDUP.CALL){
                        EventType[bndPair[i,1]] = paste("INS.DUP",endtype,sep="")
                        EventType[bndPair[i,3]] = paste("INS.DUP",endtype,sep="")
                        FoldChange[bndPair[i,1]] = res$fc
                        FoldChange[bndPair[i,3]] = res$fc
                        cat("Found ",EventType[bndPair[i,1]]," on target ",chrA,"\n")
                    } else {
                #if(is.null(res$fc)) res$fc=-1  #TODO: res$fc
                        EventType[bndPair[i,1]] = paste("TRA",endtype,sep="")
                        EventType[bndPair[i,3]] = paste("TRA",endtype,sep="")
                        FoldChange[bndPair[i,1]] = res$fc
                        FoldChange[bndPair[i,3]] = res$fc
                        cat("Found ",EventType[bndPair[i,1]]," on target ",chrA,".\n",sep="")
                    }
                } else {
                    EventType[bndPair[i,1]] = paste("DOM.INS",endtype,sep="")
                    EventType[bndPair[i,3]] = paste("DOM.INS",endtype,sep="")
                    cat("Found ",EventType[bndPair[i,1]]," on target ",chrA,".\n",sep="")
                }
            }
        }
    }
  }

  # ------- Step 6: Consider all bnds that are not yet classified as an event:
  ### Step 6a: Are chr1=chr2?  If so: is this a deletion, tandem duplication, or a missed inversion event?  (see slide)
  ### Look for deletions.
  rm(chr1,pos1,goRight1,chr2,pos2,goRight2,cid,mcid) # this allows attaching.
  attach(bndtab)
  
  sel = which(is.na(EventType) & (chr1==chr2) & !goRight1 & goRight2)
  fc=rep(1,length(sel)); nullcov=rep(0,length(sel)); len=pos2[sel]-pos1[sel]
  if(length(sel)>0){
        for(i in seq_along(sel)){    
            if(!noBam && len[i]<MAX.FC.CHECK){
                if(debug) cat("Getting fold change in coverage for region ",i,": ",paste(bndtab[sel[i],],collapse=" "),"\n",sep="")
                res=getCoverageFoldchange(chr1[sel[i]],pos1[sel[i]],pos2[sel[i]],scalefac=scalefac,bamfile1,bamfile0); 
                fc[i]=res$fc
                nullcov[i]=res$contrastcov
            } else {
                fc[i]=NA
                nullcov[i]=NA
            }       
        }
        
        sel2 = which(len<MIN.FC.CHECK | (fc<MAX.FOLDCHANGE.DEL & nullcov>MIN.READS.DELDUP.CALL))
        cat(length(sel2),"events called as DEL from single bnd entries (MAX.FOLDCHANGE.DEL=",MAX.FOLDCHANGE.DEL,", MIN.READS.DELDUP.CALL=",MIN.READS.DELDUP.CALL,", MIN.FC.CHECK=",MIN.FC.CHECK,")\n")
        cat(length(sel)-length(sel2),"breakpoints in this orientation did not get called.\n")
        EventType[sel[sel2]] = "DEL"
        EventStart[sel[sel2]] = pos1[sel[sel2]]; EventEnd[sel[sel2]]=pos2[sel[sel2]]; EventChr[sel[sel2]]=chr1[sel[sel2]]
        FoldChange[sel]=fc
  } 

  ### Look for tandem duplications.
  sel = which(is.na(EventType) & (chr1==chr2) & goRight1 & !goRight2)
  fc=rep(1,length(sel)); nullcov=rep(0,length(sel)); len=pos2[sel]-pos1[sel]
  if(length(sel)>0){
    for(i in seq_along(sel)) {
         if(!noBam && len[i] < MAX.FC.CHECK){
            cat("Getting fold change in coverage for region ",i,": ",paste(bndtab[sel[i],],collapse=" "),"\n",sep="")
            res=getCoverageFoldchange(chr1[sel[i]],pos1[sel[i]],pos2[sel[i]],scalefac=scalefac,bamfile1,bamfile0); 
            fc[i]=res$fc; nullcov[i]=res$contrastcov
         } else {
            fc[i]=NA
            nullcov[i]=NA
         }
    }
    sel2 = which(len<MIN.FC.CHECK | (fc>MIN.FOLDCHANGE.DUP & nullcov>MIN.READS.DELDUP.CALL))
    cat(length(sel2),"events called as TANDEM.DUP from single bnd entries (MIN.FOLDCHANGE.DUP=",MIN.FOLDCHANGE.DUP,", MIN.READS.DELDUP.CALL=",MIN.READS.DELDUP.CALL,", MIN.FC.CHECK=",MIN.FC.CHECK,")\n")
    cat(length(sel)-length(sel2),"breakpoints in this orientation did not get called.\n")
    EventType[sel[sel2]] = "TANDEM.DUP"
    EventStart[sel[sel2]] = pos1[sel[sel2]]; EventEnd[sel[sel2]]=pos2[sel[sel2]]; EventChr[sel[sel2]]=chr1[sel[sel2]]
    FoldChange[sel]=fc
  }
  ### Step 6b: chr1 not the same as chr2:  Event type is inter-chromosomal, report as is.
  sel = which(EventType=="NA" & chr1!=chr2)
  EventType[sel] = "INTERCHROMOSOMAL"
  Distance = bndtab$pos2-bndtab$pos1; Distance[bndtab$chr2!=bndtab$chr1]=NA
  bndtab=data.frame(bndtab,Distance,EventType,PartnerBnds,FoldChange,EventStart,EventEnd,DonorStart,DonorEnd)
} else {
  cat("no events...\n")
  EventType = rep("NA",nrow(bndtab))
  EventStart = rep(NA,nrow(bndtab)) # In case of TRA,DUP.INS, EventStart, EventEnd are site of insertion.
  EventEnd = rep(NA,nrow(bndtab))
  Distance=rep(NA,nrow(bndtab))
  DonorStart=rep(NA,nrow(bndtab))
  DonorEnd=rep(NA,nrow(bndtab))
  DonorChr=rep(NA,nrow(bndtab))
  EventChr=rep(NA,nrow(bndtab))
  PartnerBnds=rep("",nrow(bndtab))
  FoldChange=rep(NA,nrow(bndtab))
  bndtab=data.frame(bndtab,Distance,EventType,PartnerBnds,FoldChange,EventStart,EventEnd,DonorStart,DonorEnd)
}

# ------- Look for putative insertions indicated by LClip-RClip pair.
scL1tab = summarizeClusters(scL1,includeNBases=FALSE)
scR1tab = summarizeClusters(scR1,includeNBases=FALSE)
if(!is.null(chrname)){
    scL1tab=scL1tab[which(scL1tab$Chrom==chrname),]
    scR1tab=scR1tab[which(scR1tab$Chrom==chrname),]
}
hotspotsL=findHotspots(scL1tab$Position,DIST.TO.NEXT.IN.HOTSPOT=DIST.TO.NEXT.IN.HOTSPOT,HOTSPOT.CLUSTER.SIZE=HOTSPOT.CLUSTER.SIZE)
#cex=c(0.2,2)
#plot(scL1tab$Position, col=hotspotsL+1, cex=cex[hotspotsL+1])
hotspotsR=findHotspots(scR1tab$Position,DIST.TO.NEXT.IN.HOTSPOT=10000,HOTSPOT.CLUSTER.SIZE=3)
#cex=c(0.2,2)
#plot(scR1tab$Position, col=hotspotsR+1, cex=cex[hotspotsR+1])
scL1tab = scL1tab[!hotspotsL,]
scR1tab = scR1tab[!hotspotsR,]
Rmate=rep(NA,nrow(scL1tab))
cat("Finding L-R clusters that may indicate insertion ...\n")
minD=rep(NA,nrow(scL1tab))
setNegativeToNA<-function(x,thresh=-1){
    x[which(x<thresh)]=NA
    x
}
for(i in seq_len(nrow(scL1tab))){
    if(debug) if(i %%100==0) cat("\t\t",i," out of ",nrow(scL1tab),"\n",sep="")
    x = scL1tab$Position[i]
    inds = which(scR1tab$Chrom==scL1tab$Chrom[i] & (scR1tab$Position-x)>= -1 & (scR1tab$Position-x)<=rl)
    minD[i]=tryCatch( max(-Inf,min(setNegativeToNA(scR1tab$Position[scR1tab$Chrom==scL1tab$Chrom[i]]-x),na.rm=TRUE)),warning=function(w) { NA }, error=function(e) { NA } ) #DONE: scR1tab$Chrom==scL1tab$Chrom[i] could be nothing, warning captured
    if(length(inds)>0){
        Rmate[i] = inds[which.min(abs(scR1tab$Position[inds]-x))]
    }
}
if(is.na(MAX.GAP.INSERTION)){ 
    MAX.GAP.INSERTION = min(25, max(quantile(minD, 0.1,na.rm=T),10))  # between 10 and 25 bp
    cat("Estimated MAX.GAP.INSERTION from data: ",MAX.GAP.INSERTION,"\n")
}
inds=which(!is.na(Rmate))
if(length(inds)>0){
    LRpairs=cbind(inds,Rmate[inds])
    cat("There are a total of ",nrow(LRpairs),"Lclip-Rclip pairs.\n")
    #plot(scL1tab$Position[LRpairs[,1]])
    scpairs=data.frame(Chr=scL1tab$Chrom[LRpairs[,1]],Lclip=scL1tab$Position[LRpairs[,1]],Rclip=scR1tab$Position[LRpairs[,2]],
                        Lreads=scL1tab$nReads[LRpairs[,1]],Rreads=scR1tab$nReads[LRpairs[,2]],
                        cid=row.names(scL1tab)[LRpairs[,1]],mcid=row.names(scR1tab)[LRpairs[,2]],
                        pmmL=scL1tab$ProportionMismatch[LRpairs[,1]],pmmR=scR1tab$ProportionMismatch[LRpairs[,2]])
    d=scpairs$Rclip-scpairs$Lclip
    totreads=scpairs$Rreads+scpairs$Lreads
    maxreads=pmax(scpairs$Rreads,scpairs$Lreads)
    maxppr=pmax(scpairs$pmmL,scpairs$pmmR)
    if(is.na(MIN.TOTAL.READS.INS)){ 
        MIN.TOTAL.READS.INS = quantile(totreads, 0.75)  # 75-th percentile.
        cat("Estimated MIN.TOTAL.READS.INS from data: ",MIN.TOTAL.READS.INS,"\n")
    }
    
  if(debug){
     plotFile=paste(plotpf,"fig3","pdf",sep=".")
     pdf(plotFile)
     par(mfrow=c(2,1))
     if(!all(d==0)) {
       hist(d,100,col="cornflowerblue", xlab="Distance Between Lclip and Rclip Cluster", ylab="Number of Pairs",main=paste(round(100*length(d)/length(Rmate), digits=1),"% Pairs Within Read Length"))
     }
     pLessThresh=rep(0,rl); for(i in seq_len(rl)){ pLessThresh[i]=sum(totreads[d==i]>MIN.TOTAL.READS.INS)/sum(d==i)}
     if(!all(is.na(pLessThresh))) {
       plot(pLessThresh, xlab="Distance Between Lclip and Rclip Cluster", ylab=paste("P(Total Supporting Reads>",MIN.TOTAL.READS.INS,")",sep=""),type="b")
     }
     dev.off()
  }
    
    cat("Keeping calls with:\n\tat least ",MIN.TOTAL.READS.INS," total reads,\n\tat most ",
        MAX.GAP.INSERTION," gap between clusters.\n",sep="")
    throw=totreads<MIN.TOTAL.READS.INS  | d>MAX.GAP.INSERTION
    scpairs=scpairs[!throw,]
    if(nrow(scpairs)>0){
        instab=data.frame(chr1=scpairs$Chr,pos1=scpairs$Lclip,goRight1=rep(0,nrow(scpairs)),chr2=scpairs$Chr,pos2=scpairs$Rclip,goRight2=rep(1,nrow(scpairs)),
                    cid=scpairs$cid,mcid=scpairs$mcid,
                    Distance=scpairs$Rclip-scpairs$Lclip,EventType=rep("INS",nrow(scpairs)),PartnerBnds=rep(NA,nrow(scpairs)),FoldChange=rep(NA,nrow(scpairs)),
                    EventStart=scpairs$Lclip,EventEnd=scpairs$Rclip,DonorStart=rep(NA,nrow(scpairs)),DonorEnd=rep(NA,nrow(scpairs)))
        
        bndtab=rbind(bndtab,instab)
    }
}

#filter out if any breakpoint within a limit of distance to centromere or telomere 
#TODO does this breaks PartnerBnds? the other end may still remain
cat("bndtab=",nrow(bndtab),"\n")
throw=rep(FALSE,nrow(bndtab))
for(i in seq_len(nrow(bndtab))){
  #print(distanceFromGap(bndtab$pos1[i],bndtab$chr1[i],gap))
  if((distanceFromGap(bndtab$pos1[i],bndtab$chr1[i],gap)<gapdist)||(distanceFromGap(bndtab$pos2[i],bndtab$chr2[i],gap)<gapdist))
  throw[i]=TRUE
}
cat(sum(throw)," bnd entries thrown due to within ",nsf(gapdist)," to centromere or telomere.\n",sep="")
bndtab=bndtab[!throw,]
cat("final N of bndtab=",nrow(bndtab),"\n")

} #end of non saved RData loop

event_file=paste(outprefix,"sclip","RData",sep=".")
cat("outputRDATA=",event_file,"\n")
dir.create(dirname(event_file),recursive=TRUE,showWarnings=FALSE)
sclip_opt$rl=rl1
save(bndtab,resL,resR,scL1,scR1,reg_sclip,fromRData,file=event_file) #does saved this
reg_sclip=bnd2vcf(bndtab,scL1,scR1,rg_files[[n_sp]],ref_seq,sclip_opt)
reg_sclip=lapply(seq_len(length(reg_sclip)),prep_reg,reg_sclip)
reg_sclip[sapply(reg_sclip, is.null)] <- NULL
save(bndtab,resL,resR,scL1,scR1,reg_sclip,fromRData,file=event_file)
#file.copy(statfile,paste(outprefix,"par","txt",sep="."))

if(!is.null(outputVCF)) { #this is always disabled now dueto the complicated changes to bnd2vcf
  cat("outputVCF=",outputVCF,"\n")
  vcf_file=paste(outprefix,"sclip","vcf",sep=".")
  list[vcf_raw,vcf_meta,bed_raw]=vcf_merge(reg_sclip,sample_id,sample_info,contig_info,ref_seq,ref_file)
  if(nrow(vcf_raw)!=0) { 
    sortnames=c("CHROM","POS")
    vcf_raw[do.call("order", vcf_raw[sortnames]), ]
    #sv_order=sort(vcf_raw$POS,vcf_raw$CHROM,index.return=T)$ix
    #vcf_raw=vcf_raw[sv_order,]
    for(vi in seq_len(nrow(vcf_raw))) vcf_raw$ID[vi]=paste(sample_id,vi,vcf_raw$ID[vi],sep=".")
    cat("==Info: raw sorted!\n")
    write_com(vcf_raw,comment=vcf_meta,filename=vcf_file,gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
    cat("==Info:total ",nrow(vcf_raw),"raw calls, ",vcf_file,"written!\n")
  }
}
cat("=Info: done,",taggie(gtk),"\n")
if(debug) { cat("=Info: warnings if any\n"); warnings() }
cat("please find outputs in:",event_file,"and",outputVCF,"\n",sep=" ")
cat("SCLIPDONE\n")
