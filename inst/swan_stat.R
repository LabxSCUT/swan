#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
####Externalities####
library(swan)
for(p in c("optparse")) myrequire(p)
for(p in c("Biobase","GenomicRanges","GenomeInfoDb","Rsamtools")) myrequire(p,repo="Bioc")
#suppressMessages(library(optparse))
#suppressMessages(library(GenomicRanges))
#suppressMessages(library(Rsamtools))
#suppressMessages(library(Biobase))
#suppressMessages(library(GenomeInfoDb))

#this has to come early than options because of option defaults
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

####CommandLines####
option_list2 <- list(
  make_option(c("-x", "--xmax"), type="integer", default=2000,
              help="x limit on histogram, [default %default]"),
  make_option(c("-y", "--yield"), type="integer", default=1000000,
              help="number of reads to be sampled for stats [default %default]"),
  make_option(c("-c", "--chrname"), default="A",
              help="use all if bam have all chromosomes or use the chr name if only have one, or chr name separated by ',' if multiple [default]"),
  make_option(c("-o", "--oprefix"), default="none",
              help="bam-wise stat output prefixs [default %default] "),
  make_option(c("-m", "--mprefix"), default="none",
              help="merged stat output prefixs [default %default] "),
  make_option(c("-s", "--step"), type="integer", default=10,
              help="bin step size on histogram [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)
parser <- OptionParser(usage = "%prog [options] rg1.bam,rg2.bam,...", option_list=option_list2)
cat("-Info: swan_stat vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=1){ print_help(parser); quit() }
verbose=cmd$options$noQuiet
debug=cmd$options$debug
if(debug){ options(warn=2) } else { options(warn=-1) }
ptm=proc.time() #time keeper
xmax=cmd$options$xmax
binsize=cmd$options$step
yield=cmd$options$yield
chrname=cmd$options$chrname
what=c("pos","mpos","isize","strand","qwidth","cigar","qname","flag")
max_chr_len=3e+8
bam_files=strsplit(cmd$args,",")[[1]]
oprefixs=strsplit(cmd$options$oprefix,",")[[1]]; #if libwise oprefix is specified, will use specified
if(cmd$options$oprefix=="none" || length(oprefixs)!=length(bam_files)) {
  oprefixs=gsub(".bam$","",bam_files) #oprefixs=libwise statfiles, mprefix=merged statfiles
}
cat("==oprefixs:",oprefixs,"\n")
if(cmd$options$mprefix=="none") {
  mprefix=lcPrefix(oprefixs)
} else {
  mprefix=cmd$options$mprefix
} 
cat("==mprefix:",mprefix,"\n")
bam_tmps=paste(oprefixs,"stat.bam",sep=".")
cat("==bam_tmps:",bam_tmps,"\n")

myseqinfo=function(bam_file){
  raw=system(sprintf("samtools view -H %s | grep @SQ",bam_file),intern=T)
  mat=as.matrix(as.data.frame(strsplit(raw,"\t")))
  chr=gsub("SN:","",mat[2,])
  chr_len=as.numeric(gsub("LN:","",mat[3,]))
  seq=Seqinfo(seqnames=chr, seqlengths=chr_len)
  return(seq)
}

mycountbam=function(bam_file,chr=""){
  raw=system(sprintf("samtools idxstats %s",bam_file),intern=T)
  mat=as.data.frame(t(as.matrix(as.data.frame(strsplit(raw,"\t")))),stringsAsFactors=F) #all characters
  #set_colClass(mat,c("character","numeric","numeric","numeric"))
  colnames(mat)=c("seqname","seqlength","mapped","unmapped")
  if(chr!=""){
    nreads=as.numeric(mat$mapped[mat$seqname==chr])
  } else {
    nreads=sum(as.numeric(mat$mapped))
  }
  return(nreads)
}
  
mtable=data.frame()
for(ni in seq_along(bam_files)){
  bam_file=bam_files[ni]
  bam_tmp=bam_tmps[ni]
  cat("==bam_file:",bam_file,"\n")
  oprefix=oprefixs[ni]
  #bam_full=BamFile(bam_file) #TODO: significant time spend here, reduce
  bam_header=myseqinfo(bam_file)
  if(chrname %in% c("a","A","all","ALL","All")){
    chr_length=sum(as.numeric(seqlengths(bam_header)))
    nreads=mycountbam(bam_file)
  } else {
    chr_length=as.numeric(seqlengths(bam_header)[seqnames(bam_header)==chrname])
    nreads=mycountbam(bam_file,chr=chrname)
  }
  cat("reads:", nreads, "\n")
  cat("sampled :", yield, "reads\n")
  frac=min(yield/nreads,1)
  cat("chr_name:", chrname,"\n")
  cat("chr_len:", chr_length,"\n")
  cat(sprintf("samtools view -bh -s %f %s >%s",frac,bam_file,bam_tmp),"\n")
  system(sprintf("samtools view -bh -s %f %s >%s",frac,bam_file,bam_tmp)) 
  cat(sprintf("samtools index %s",bam_tmp),"\n") 
  system(sprintf("samtools index %s",bam_tmp)) 
  bam_sample=BamFile(bam_tmp) #TODO: significant time spend here, reduce
  param <- ScanBamParam(what=what,
                      flag=scanBamFlag(isDuplicate=F,isNotPassingQualityControls=F))
  stats <- scanBam(bam_sample, param = param)[[1]]
  rl_mean=mean(stats[["qwidth"]],na.rm=TRUE)
  if(is.na(rl_mean)) rl_mean=0
  isize_mean=mean(abs(stats[["isize"]]),na.rm=TRUE)
  isize_sd=sd(abs(stats[["isize"]]),na.rm=TRUE)
  isize_med=median(abs(stats[["isize"]]),na.rm=TRUE)
  isize_sd_mad=mad(abs(stats[["isize"]]),na.rm=TRUE)
  n_plus=sum(stats[["strand"]]=="+",na.rm=T)
  n_minus=sum(stats[["strand"]]=="-",na.rm=T)
  p_left=sum(stats[["flag"]] %in% hang_flags & stats[["strand"]]=="+",na.rm=T)/(n_plus+1)
  p_right=sum(stats[["flag"]] %in% hang_flags & stats[["strand"]]=="-",na.rm=T)/(n_minus+1)
  q_left=sum(grepl(tail_pattern_new,stats[["cigar"]],perl=TRUE),na.rm=T)/(n_plus+1)
  q_right=sum(grepl(head_pattern_new,stats[["cigar"]],perl=TRUE),na.rm=T)/(n_minus+1)
  cat("rate of plus read with mate hang",p_left,"\n")
  cat("rate of minus read with mate hang",p_right,"\n")
  cat("rate of plus read with tail clip",q_left,"\n")
  cat("rate of minus read with head clip",q_right,"\n")
  cat("read length:", rl_mean,"\n")
  cat("mean isize (red):", isize_mean,"\n")
  cat("sd isize (red):", isize_sd,"\n")
  cat("median isize (green):", isize_med,"\n")
  cat("sd by mad isize (green):", isize_sd_mad,"\n")
  isize_V=stats[["isize"]][stats[["isize"]]>0&!is.na(stats[["isize"]])]
  cat("find isize by mode, also seediagnostics in ",bam_file,".xxx.hist.pdf \n",sep="")
  if(is.na(isize_mean)) isize_nancy=list(target.isize=0,sdR=0,sdL=0,good.for.deletion=0,good.for.insertion=0) else isize_nancy=isizePar(isize_V,MAX.ISIZE=NA,doplot=TRUE,method=2,opre=oprefix)
  cat("use nancy's robust method\n")
  mean_cvg=nreads*rl_mean/chr_length
  if(is.na(mean_cvg)) mean_cvg=0
  cat("mean_cvg:", mean_cvg,"\n")
  file.remove(bam_tmp)
  file.remove(paste(bam_tmp,"bai",sep="."))
  print(isize_nancy)
  #cat(paste(rl_mean, mean_cvg, p_left, p_right, q_left, q_right, nreads, sep=","),"\n")
  par=as.data.frame(list(rl=rl_mean,cvg=mean_cvg,p_left=p_left,p_right=p_right,q_right=q_right,q_left=q_left,
  nreads=nreads,
  is=isize_nancy$target.isize,sdR=isize_nancy$sdR,sdL=isize_nancy$sdL,
  lCd=isize_nancy$good.for.deletion,
  lCi=isize_nancy$good.for.insertion,
  lDl=p_left<=p_fail,lDr=p_right<=p_fail,lSl=q_left<=q_fail,lSr=q_right<=q_fail))
  print(par)
  mtable=rbind(mtable,par)
  parf=paste(oprefix,"stat",sep=".")
  write.table(format(par,digits=2,scientific=F),parf,quote=F,row.names=F,sep="\t")
  cat("written",parf,"\n")
  if(debug) save.image(paste(parf,"RData",sep="."))
}
if(mprefix!="none") {
  mfile=paste(mprefix,"stat",sep=".") 
  write.table(format(mtable,digits=2,scientific=F),mfile,quote=F,row.names=F,sep="\t")
  cat("written",mfile,"\n")
}
elapsed=proc.time()-ptm
cat("bam_stat took",elapsed[3],"seconds.\n")
if(debug) { cat("warnings if any\n"); warnings() }
cat("STATDONE\n")
