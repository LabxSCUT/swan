#!/usr/bin/env Rscript
gtk=proc.time()[3]
suppressMessages(library(optparse))
suppressMessages(library(swan))
suppressMessages(library(Rsamtools))
#suppressMessages(library(Cairo))
options(warn=2)
gmk=get_gmk(Sys.getpid())

parse.one <- function(res, result) {
  m <- do.call(rbind, lapply(seq_along(res), function(i) {
    if(result[i] == -1) return("")
    st <- attr(result, "capture.start")[i, ]
    substring(res[i], st, st + attr(result, "capture.length")[i, ] - 1)
  }))
  colnames(m) <- attr(result, "capture.names")
  m
}

option_list <- list(
  make_option(c("-f", "--format"), default="bed",
              help="bed/vcf format, [default: %default] \n"),
  make_option(c("-a", "--bam0"), default="none",
              help="control bam file, [default: %default] \n"),
  make_option(c("-b", "--alt"), default="DEL",
              help="alternative, [default: %default] \n"),
  make_option(c("-o", "--outdir"), default="plot",
              help="output prefix, [default: %default] \n")
	)	
parser <- OptionParser( 
  usage="%prog [options] vcfFile spX.bam", 
  option_list=option_list)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
print(cmd$args)
if(length(cmd$args)!=2){ print_help(parser); quit(); }

varfile=cmd$args[1]
bamfile1=strsplit(cmd$args[2],",")[[1]]
format=cmd$options$format
chrs=cmd$options$chrname
alt=cmd$options$alt
print(cmd$options)
outdir=cmd$options$outdir
if(cmd$options$bam0=="none") bamfile0=NULL else bamfile0=strsplit(cmd$options$bam0,",")[[1]]

if(!file.exists(outdir)) dir.create(file.path(getwd(),outdir))
regions=read.table(varfile,sep="\t",header=FALSE,stringsAsFactors=F)
if(format=="vcf"){
	names(regions)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","V10","V11")[1:ncol(regions)]
	regions$ALT=paste(regions$ALT)
	isSNP = nchar(paste(regions$REF))==1 & nchar(paste(regions$ALT))==1
	if(length(isSNP)>0) regions=regions[!isSNP,]
	isMSK = (regions$ALT=="<MSK>")
	if(length(isMSK)>0) regions=regions[!isMSK,]
	isIGN = (regions$ALT=="<IGN>")
	if(length(isIGN)>0) regions=regions[!isIGN,]
	svlen_pat="SVLEN=(?<size>[0-9]+)"
	len=get_al_omp(regions$INFO,svlen_pat)
	regions$LEN=len
} else {
	if(ncol(regions)==3) regions=cbind(regions, list(ALT=alt))
	names(regions)=c("CHROM","POS","END","ALT",rep("V",ncol(regions)-4))
	regions$LEN=regions$END-regions$POS
}
ord = order(regions$CHROM, regions$POS)
regions.plot=regions[ord,]
header=scanBamHeader(bamfile1[1])
chrnames=names(header[[1]]$targets)  # find out what the chromosomes are called.
if(nrow(regions.plot)>1){
  for(ix in 1:nrow(regions.plot)){
		#print(regions.plot[ix,])
    cat("Processing region ",ix," out of ", nrow(regions.plot),"...\n",sep="")
		#print(len); print(regions.plot$POS[ix]);
    st=max(1,as.integer(regions.plot$POS[ix])-200); ed=st+abs(regions.plot$LEN[ix])+400
		if(is.na(st)) { print(regions.plot$POS[ix]); next }
		if(is.na(ed)) { print(regions.plot$LEN[ix]); next }
    cat(regions.plot[ix,ncol(regions.plot)-1],"\n")
    pur_found=gregexpr("PURITY=[-+]?([0-9]*\\.[0-9]+|[0-9]+)",regions.plot[ix,ncol(regions.plot)-1],perl=T)[[1]]   
    pur=as.vector(parse.one(regions.plot[ix,ncol(regions.plot)-1],pur_found))
    chr=which(chrnames==regions.plot$CHROM[ix])
    buffer=min(max(floor((ed-st)*0.5),500),5000)
    st0=max(1,st-buffer); ed0=ed+buffer
    #filename=paste(outdir,"/Region_Chr",chr,"_ST",st0,"_ED",ed0,"_PU",pur,"_",regions.plot$ALT[ix],"_LEN",ed-st,".png",sep="")
    filename=paste(outdir,"/Region_Chr",chr,"_ST",st0,"_ED",ed0,"_PU",pur,"_",regions.plot$ALT[ix],"_LEN",ed-st,".pdf",sep="")
    label=filename
    #label=paste("SV Length: ",ed-st,", SV Type: ",regions.plot$ALT[ix])
    cat("Region size: ",ed0-st0," bases, output to ", filename,"\n",sep="")
    #CairoPNG(filename, height=6,width=4, units="in", res=1200, pointsize=4)
		#par(mar = c(6,6,6,6))
    pdf(filename, height=24,width=12)
		par(mar = c(1,1,1,1))
    cat("plotting: ",bamfile1,bamfile0,chr,st0,ed0,st,ed,label,"\n")
    res=plotRegions(bamfile1=bamfile1,bamfile0=bamfile0,chr=chr,st0=st0,ed0=ed0,st=st,ed=ed,label=label,covymax=5)
    dev.off()
  }
}
