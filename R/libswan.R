#TODO: cleanup code for paper submission
myrequire = function(pkg, repo="CRAN", ...){
  cat("requiring package", pkg, "\n")
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    if(repo!="CRAN"){
      source("http://bioconductor.org/biocLite.R")
      biocLite(pkg,...)
    } else {
      install.packages(pkg,repo="http://cran.us.r-project.org",...)
    }
  })
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    stop(pkg," was not installed and cannot install on the fly!\n")
  })
}

#TODO: do documentation and fix cross-reference
#TODO: to pass R CMD CHECK
#suppressMessages(library(ggplot2))
#suppressMessages(library(gridExtra))
#suppressMessages(library(seqCBS))
for(p in c("Rcpp","RcppArmadillo")) myrequire(p,type="source")
for(p in c("hash","sets","plyr","stringr","data.table","robustbase")) myrequire(p)
for(p in c("Biostrings","Rsamtools","IRanges","GenomicRanges")) myrequire(p,repo="Bioc")
#suppressMessages(library(hash))
#suppressMessages(library(sets))
#suppressMessages(library(plyr))
#suppressMessages(library(stringr))
#suppressMessages(library(data.table))
#suppressMessages(library(robustbase))
#suppressMessages(library(Rcpp))
#suppressMessages(library(RcppArmadillo))
#suppressMessages(library(Biostrings))
#suppressMessages(library(Rsamtools))
#suppressMessages(library(IRanges))
#suppressMessages(library(GenomicRanges))
Rcpp::loadModule("pairmap_module", TRUE)
#CODE: how to load an external module
#Rcpp::loadModule("hashmap_module", TRUE) 
#Rcpp::loadModule("hashset_module", TRUE)
#options(warn=2) #warning will error
#dynamic unload of swan library
#library(swan)
#detach("package:swan", unload = TRUE)
#library.dynam()
#library.dynam.unload("swan", system.file(package = "swan"))


#CODE: remove null items from list
#lst[sapply(lst, is.null)] <- NULL

#CODE: to generate documentation for swan
#R -e 'library(devtools);document()' 
#edit man/swan-package.Rd for change of package descriptions

#NOTE: fail safe fool proof values
min_trunk_MPRs=100
mixing_rate_fail=0.75 # now in default, use 0.75 instead of 0.5 to boost signal
width_fail=100 #  now in default, use 200 instead of 100 to boost signal
lw_width_fail=1000 #
delta_fail=100000 # maximum possible delta preset, this decides overlap region
flank_fail=10000 # minimum flank|ligation to avoid edge effect
#NOTE: maxInsert > bigDel > Delta > minInsert, as designed
smallDel_fail=20
smallIns_fail=20
Delta_sd=6
bigDel_sd=10     #over stringent?
#bigDel_sd=3     #good with Broad and BCM, broke down with WashU
maxInsert_fail=50000 #define the max insert allowed to have lCd contribution
maxbigDel_fail=1000000 #define the possible upper limit of bigDel, 1M
minbigDel_fail=300 #define the possible lower limit of bigDel, 0.3k
maxDelta_fail=2000 #define the possible upper of delta
minDelta_fail=1000 #define the possible lower of delta
minInsert_fail=30 #define the min insert allowed to have lCd contribution
p_fail=0.01 # we typically see this less than 1/100 in a good library 
q_fail=0.01 # we typically see this less than 1/100 in a good library
max_chr_len=3e+8 #  aovid integer burst, should be larger than any single chr_len
prop_clip_fail=0.5 # if less than prop_clip_fail*RL aligned => not suitable for isize inference
soft_cut_fail=20 # if more than soft_cut_fail bases soft-clipped => consider non-chance clipping
fy_cap_fail=20  # single read score cap of fy ratio for lCd 
mem_save=FALSE
mvsum_fail=21
alpha_level=c(0.1,0.05,0.01,0.005,0.001) # emperical levels 
complex_pattern="[SH]"                        #pattern 
softclip_pattern="[S]"                        #pattern 
simple_pattern="[0-9]+M"                        #pattern 
alr_pattern="(?<size>[0-9]+)[MX=I]"      # actual aligned length
alih_pattern="^(?<size>[0-9+])I"           # exclude ^I and I$ things
alit_pattern="(?<size>[0-9+])I$"
alg_pattern="(?<size>[0-9]+)[MX=DN]"     # actual aligned length to genome
ald_pattern="[^D](?<size>[0-9]+)[D][0-9]"          # read embed deletion size
ali_pattern="[^I](?<size>[0-9]+)[I][0-9]"          # read embed insertion size
#NOTE: Some bams do have xIyM... and ...xMyI patterns, what does that mean? Is it different from S
head_pattern="^(?<size>[0-9]+)S" 
tail_pattern="(?<size>[0-9]+)S$"
head_pattern_new="^(?<size1>[0-9]+)[SHI]|(?<size2>[0-9]+)S[0-9]+[MIDNPX=]"
tail_pattern_new="(?<size1>[0-9]+)[SHI]$|[0-9]+[MIDNPX=](?<size2>[0-9]+)S"
#NOTE: for flag summary see manual/sam_output.pdf
#Let's for now exclude reads without proper pair marking, i.e. 0x2 is not set
#current understanding of the direction and position relationship
#1. forward == read is from the reference strand (5->3),
#2.  reverse == read is from the reverse complementary of reference strand (5->3),
#3. whether or not the read is forward/reverse, SEQ always give the 5->3, rc if necessary
#4. in shotgun sequencing, first/second is not correlated with forward/reverse
#5. in target sequencing like OS-seq, first read is always forward while second read is always reverse?
#conc_flags=c(99,147,97,145,83,163,81,161) #this is all concordant in directions and relative position, regardless of isize
impp_flags=c(97,145,81,161) #inprop pair only cares FR pair with TLEN<=0 or TLEN>big or RNEXT is another chr
conc_flags=c(99,147,83,163) #prop pair cares FR pair with TLEN proper
hang_flags=c(73,137,121,185,105,169,89,153) 
#NOTE: 73=105 1m+, 89=121 1m-, 137=169 2m+, 153 == 185 2m-
#1m+ is 1st read mapped to forward, try decode 1m-, 2m+, 2m- similarily
#89, 153, 105, 169 never exists though, it is con
disc_flags=c(65,129,113,177,67,115,131,179) #same direction with correct/incorrect isize 
umap_flags=c(77,141,133,165,181,101,117,69) #both/itself unmapped
rmfw_flags=c(73,105,169,137,99,163,67,131,161,97,65,129) #all read mapped forward    
rmrv_flags=c(89,121,153,185,83,147,115,179,81,145,113,177) #all read mapped reverse
#NOTE: disc pair cares FF pair or RR pair whether RNEXT is either *,=,chr
#67,115,131,179 is now included in disc_flags according to the breakdown found on inet
#insert size less than RL is mostlikely due to some erroroneous library; now we take this part out in isize calculation
scan_cols_preset=c("start","lW","lCd","lCi","lDr","lDl","lSr","lSl","cvg","cCd","cCi","cDr","cDl")
scan_cols_preset=c(scan_cols_preset,"ins","del","HAF","HAR") #additional scan columns, must be the last 4
#NOTE: window start, sore for Lw, Lcd, Lci, Ldr, lDl, lSl (when sc are included), coverage,
#read count for Lcd, Lci, Ldr, Ldl 
#cigar of ins, del, read count for forward read hanging and reverse read hanging
out_cols_preset=scan_cols_preset

#NOTE: global constants by Nancy, need to clean them up by METHOD.WHAT.FOR
MIN.READS.PER.CLUSTER=5 # Minimum reads in soft clipping cluster.
SC.FROM.BK.THRESH=1000   # Maximum distance of soft clipping to breakpoint 
SC.PROPORTION.MISMATCH.THRESH = 0.01
SC.LEN.THRESH = 30  # Minium length of soft clipped sequence for remapping.
SC.MIN.BUFFER=10000
SC.MAX.BUFFER=50000
gainRelCNTHRESH = 1.2
lossRelCNTHRESH = 0.8
MIN.COV.SOFTCLIP.OVERRIDE = 2 #minimum reads coverage
MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE = 2
NUM.READS.FOR.CI = 100
CI.LEN = 500
REMAP.RELCN.RATIO.THRESH = 2
PRECISE.BKPOINT.OFFSET.ALLOWANCE = 20 # Should be some value that fits well within a read.
BIC.OVERSHOOT=10
LCD.MIN.DEL=1000 #deletion >1000bp will be found by lcd

#NOTE: some utility functions
asc <- function(x) { strtoi(charToRaw(x),16L) } # check if x contains a wholenumber
qscore <- function(x) { asc(x)-asc("!") }       # convert letter quality to qscore
pscore <- function(x) { 10**(-0.1*qscore(x)) }  # convert qscore to error probability
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol #
get_gtk=function(gtk){ proc.time()[3]-gtk }     # get global time
get_gmk=function(pid){ tryCatch( {paste("(rss,vsz)",system(paste("ps -p",pid,"-o rss,vsz | tail -n 1"),intern=T,ignore.stderr=T)) }, error = function(e) { "gmk err" } ) } # get global memory
taggie=function(gtk){ paste("gtk=",get_gtk(gtk),"s,gmk=",suppressMessages(get_gmk(Sys.getpid())),"kb",sep="") } # format gtk and gmk to attachable tag

#FUNC: options parser 1
parse_opt=function(opt_text,fs=c(",","%"),as="="){
  if(opt_text=="") return(list())
  opt_tmp=strsplit(strsplit(opt_text,split=fs[1])[[1]],split=as)
  opt_names=sapply(opt_tmp,"[",1)
  opt_values=strsplit(sapply(opt_tmp,"[",2),split=fs[2])
  names(opt_values)=opt_names
  opt_values
} #parse_opt("a=1%2,b=2%3")

#FUNC: options parser 2
parse_opt2=function(opt_text,fs=c("_",","),as="="){
  if(opt_text=="") return(list())
  opt_tmp=strsplit(opt_text,split=fs[1])[[1]]
  opt_values=lapply(opt_tmp,parse_opt,fs=c(fs[2],fs[1]),as=as)
  opt_values
} #parse_opt2("track=lCd,thresh=level3,tele=100_track=lDl+lDr,tele=20")

#FUNC: argument splitter
split_args=function(arg_text,guide_lst,arg_num=1){
  args=lapply(seq_len(arg_num),function(i,x){x},guide_lst)
  if(grepl("^[0-9]",arg_text)){ #1_2,3_4:5_6; 1,2:3
    tmp=strsplit(strsplit(arg_text,split=':')[[1]],split=",")
    for(a in seq_len(arg_num))
      for(s in seq_along(tmp)){
        args[[a]][[s]]=rep(NA,length(guide_lst[[s]]))
        if(length(guide_lst[[s]])==length(tmp[[s]])){
          for(r in seq_along(tmp[[s]]))
            args[[a]][[s]][r]=tmp[[s]][[r]][a]
        }else{
            args[[a]][[s]]=rep(tmp[[s]][[1]][a],length(guide_lst[[s]]))
        } 
      }
  } else {
    for(a in seq_along(args))
      for(s in seq_along(args[[a]]))
        args[[a]][[s]]=rep(arg_text,length(guide_lst[[s]]))
  }
  args
} #list[arg1,arg2]=split_args(arg_text,test_lst,2)

#FUNC: set argument
set_args=function(arg_text,guide_lst,arg_num=1){
  args=lapply(seq_len(arg_num),function(i,x){x},guide_lst)
  if(grepl("^[0-9]",arg_text)){ #1_2,3_4:5_6; 1,2:3
    for(a in seq_len(arg_num))
      for(s in seq_along(guide_lst)){
        args[[a]][[s]]=rep(NA,length(guide_lst[[s]]))
        for(r in seq_along(guide_lst[[s]]))
          args[[a]][[s]][r]=arg_text
      }
  } else {
    for(a in seq_along(args))
      for(s in seq_along(args[[a]]))
        args[[a]][[s]]=rep(arg_text,length(guide_lst[[s]]))
  }
  args
}

#OBSO: functions for swan_join ###
#all these call functions will take its input files and output a list of SV events
#the break point has to be as accurate as possible and 
#in load_bam buffer region is automatically loaded base on sv_type, sv_length and ci's
#where each list item has following required field:
#sv_start, sv_end, sv_start_ci_low, sv_start_ci_high, sv_end_ci_low, sv_end_ci_high, sv_type, sv_mtd, sv_val
#each above field is a integer reflect best guess from overlaping with softclipping
#sv_mtd_aux, sv_val_aux, sv_start_aux, sv_end_aux
#each above field is a vector and will extend if ovlp is found with sclip cluster
#sv_remap, sv_bam_cvg, sv_bam_hang, sv_bam_soft, sv_bam_strad
#each above field is a boolean and marks the confirmation received for this sv
#if overlap with sclip clusters, we keep the sclip cluster one as the main entry
#remove those ovelap entries from their own list, append method and score to sclip entry

#FUNC: test if ends (including ci) of SV overlap 
ovlap_sv_end=function(chr1,pos1,pos1_ci_low,pos1_ci_high,chr2,pos2,pos2_ci_low,pos2_ci_high){
  st1=pos1+pos1_ci_low;ed1=pos1+pos1_ci_high;st2=pos2+pos2_ci_low;ed2=pos2+pos2_ci_high;
  return((chr1==chr2)&&((st1>=st2&st1<=ed2)||(ed1>=st2&ed1<=ed2)||(st2>=st1&st2<=ed1)||(ed2>=st1&ed2<=ed1)))
}

#FUNC: try read a table file and return NULL if not successful
safe_read=function(filename,...){
  res=tryCatch(read.table(filename,...),error=function(e){
    cat("Warn:",filename,"empty or unreadable!\n"); NULL
  })
  res
}

#FUNC: read in a gap file or return NULL
read_gap=function(gapfile){
  if(gapfile=="" | gapfile=="none" | gapfile=="None" | is.null(gapfile)){
    gap=NULL
  } else {
    gap=read.table(file=gapfile,comment.char='#',stringsAsFactors=F,sep='\t',header=F,colClass=NA)
    gap=if(grepl(gapfile,".txt.gz")) gap[,2:4] else gap[,1:3]
    gap=set_colClass(gap,c("character",rep("numeric",2)))
    colnames(gap)=c("chrom","chromStart","chromEnd")
  }
  gap
}

#' @title call_seqcbs
#'
#' @description \code{call_seqcbs} \cr
#' If called, \code{call_seqcbs} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_seqcbs(seqname="1",
#'  seqcbs_file="spX.seqcbs.txt",
#'  seqcbs_parf="spX.seqcbs.par.txt",
#'  seqcbs_opt="learn",
#'  rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_seqcbs=function(seqname,seqcbs_file,seqcbs_parf,seqcbs_opt,rg_files){ 
	#seqcbs_parf should contain scalefac
  reg_seqcbs=list(); reg_seqcbs_good=list(); idx=1; jdx=1; seqcbs_win_size=1000; par_seqcbs=NULL
  tryCatch({tab_seqcbs=read.table(seqcbs_file,header=T)},error=function(e){cat("Warn:",seqcbs_file,"empty, nonexist or erroneous formatted\n"); return(reg_seqcbs)})
  tryCatch({par_seqcbs=read.table(seqcbs_parf,header=T)},error=function(e){cat("Warn:",seqcbs_parf,"empty, nonexist or erroneous formatted\n"); })
  scale_frac=if(!is.null(par_seqcbs)) par_seqcbs$cvgrxy else 1 #cvgrxy is in seqcbs.par.txt
  minstat = if(!is.null(seqcbs_opt$minstat)) as.numeric(seqcbs_opt$minstat) else 0
  sup = if(!is.null(seqcbs_opt$sup)) as.numeric(seqcbs_opt$sup) else 0
  gap = if(!is.null(seqcbs_opt$gap)) as.numeric(seqcbs_opt$gap) else 0
  expand = if(!is.null(seqcbs_opt$expand)) as.numeric(seqcbs_opt$expand) else 0

  rsize=rep(0,nrow(tab_seqcbs))
  for(ri in seq_along(rsize)){
    rsize[ri]=tab_seqcbs$cptR[ri]-tab_seqcbs$cptL[ri]
  }
  rsize_sort=sort(rsize,index.return=T)$ix
  
  #sort seq cbs regions in order of size; 
  #if seqbsequent regions overlap with previous regions; throw those regions  
  #print(tab_seqcbs[rsize_sort,])
  runique=IRanges(); throw=rep(FALSE,nrow(tab_seqcbs))
  for(ri in rsize_sort){
    rtmp=IRanges(start=tab_seqcbs$cptL[ri],end=tab_seqcbs$cptR[ri])
    if(sum(countOverlaps(rtmp,runique))>=1)
      throw[ri]=TRUE   
    else
      runique=c(runique,rtmp)
  }
  tab_seqcbs=tab_seqcbs[!throw,]  
  #print(throw)
  sel = which(tab_seqcbs$stat>minstat)
  #print(tab_seqcbs)
  pile_track=c(IRanges::as.vector(IRanges::coverage(IRanges(start=tab_seqcbs$cptL[sel],end=tab_seqcbs$cptR[sel]))),rep(0,1000000)) #add some buffer on edge
  tmp_track=pile_track>0
  #cat("tmp_track="); print(tmp_track)
  cbs_run=filter_calls(tmp_track,supnum=sup,gaplim=gap+2*expand) 
  cat("sup=",sup,"gap=",gap,"expand=",expand,"gap+2*expand=",gap+2*expand,"\n")
  for(ri in rsize_sort){
    if(!is.na(runValue(cbs_run)[ri])&runValue(cbs_run)[ri]>0) {
      res=getCoverageFoldchange(seqname,start(cbs_run)[ri]-expand,end(cbs_run)[ri]+expand,scale_frac,rg_files[[2]],rg_files[[1]]) #first normal, then tumor
      fc=res$fc
      #if(length(fc)==0) fc=1 #FIXME: what does numeric(0) means here? foced to 1
      if(fc>seqcbs_opt$good) { #TODO: scalfac?
        reg_seqcbs_good[[jdx]]=list(sv_chrname=seqname, sv_start=start(cbs_run)[ri]-expand, sv_end=end(cbs_run)[ri]+expand, sv_start_ci_low=-seqcbs_win_size, sv_start_ci_high=seqcbs_win_size, sv_end_ci_low=-seqcbs_win_size, sv_end_ci_high=seqcbs_win_size, sv_type="DUP", sv_alt="<DUP>",sv_mtd="cbs", sv_val=fc, sv_remap=FALSE, sv_filter=TRUE, sv_qual=6, sv_id="cbs", sv_cvg=-1, sv_imprecise=TRUE)
        jdx=jdx+1
      } else {
        if(fc>1) sv_type="DUP" else sv_type="DEL"
        reg_seqcbs[[idx]]=list(sv_chrname=seqname, sv_start=start(cbs_run)[ri]-expand, sv_end=end(cbs_run)[ri]+expand, sv_start_ci_low=-seqcbs_win_size, sv_start_ci_high=seqcbs_win_size, sv_end_ci_low=-seqcbs_win_size, sv_end_ci_high=seqcbs_win_size, sv_type=sv_type, sv_alt=paste("<",sv_type,">",sep=""),sv_mtd="cbs", sv_val=fc, sv_remap=FALSE, sv_filter=FALSE, sv_qual=NA, sv_id="cbs", sv_cvg=-1, sv_imprecise=TRUE)
        idx=idx+1
      }
    }
  }
  return(list(reg_seqcbs,reg_seqcbs_good))
}

#' @title bnd2vcf
#'
#' @description \code{bnd2vcf} \cr
#' If called, \code{bnd2vcf} will return a list of sv_event converted from bndtab
#'
#' @param bndtab: break point table
#' @param scL: soft clipping cluster on Left
#' @param scR: soft clipping cluster on Left
#' @param bam_files: bam files possibly read group separated
#' @param ref_seq: set of reference sequences
#' @param sclip_opt: options for sclip call
#' @param rle_cvg: rle object of coverage (inactive)
#' @param drop_na: whether NA type should be dropped
#'
#' @note
#'
#' @examples
#' bnd2vcf(bndtab=bndtab,
#' scL=scL,scR=scR,
#' bam_files=c("spX.rg1.bam","spX.rg2.bam"),
#' ref_seq=ref_seq,
#' sclip_opt="learn",
#' rel_cvg=rle(), drop_na=F
#' )
#'
bnd2vcf=function(bndtab,scL,scR,bam_files,ref_seq,sclip_opt,rle_cvg=NULL,drop_na=FALSE){

  if(nrow(bndtab)==0) return(list())
  bndtab$EventType=as.character(bndtab$EventType) #convert to character
  bndtab$ChrChrMatch= bndtab$chr1 == bndtab$chr2   #T if and only if chr1==chr2
  ci_sclip=100

  for(i in seq_len(nrow(bndtab))){ #assume single entry for a reciprocal pair; 
    if( grepl("DEL",bndtab$EventType[i]) && na_equal(bndtab$ChrChrMatch[i],TRUE) ) { 
      bndtab$EventType[i]="DEL"
    } else if( grepl("DUP",bndtab$EventType[i]) && na_equal(bndtab$ChrChrMatch[i],TRUE) ) {
      bndtab$EventType[i]="DUP"
    } else if( grepl("INV",bndtab$EventType[i]) && na_equal(bndtab$ChrChrMatch[i],TRUE) ) {
      bndtab$EventType[i]="INV"
    } else if( grepl("INS",bndtab$EventType[i]) && na_equal(bndtab$ChrChrMatch[i],TRUE) ) {
      bndtab$EventType[i]="INS"
    } else if( grepl("INTERCHROMOSOMAL",bndtab$EventType[i]) ) {
      bndtab$EventType[i]="TRA"
    } else if( grepl("TRA",bndtab$EventType[i]) ) {
      bndtab$EventType[i]="TRA"
    } else if( grepl("TRP",bndtab$EventType[i]) ) {
      bndtab$EventType[i]="TRA"
    } else {
      bndtab$EventType[i]="TRA"
    }
  } #cleaning up type and chromosome inconsistencies, consider all as TRA by default

  idx=1;reg_sclip=list(); mappedbuffer=15; clippedbuffer=5
  #cat("scL:",length(scL),"\n"); cat("scR:",length(scR),"\n")
  purity=rep(NA,nrow(bndtab)); nalt=rep(NA,nrow(bndtab));
  totreads=rep(NA,nrow(bndtab)); nref=rep(NA,nrow(bndtab));
  if(!is.null(bam_files)){
    for(sx in seq_len(nrow(bndtab))){
      bam1 = load_bam3(1,list(list(chrom=bndtab$chr1[sx],start=bndtab$pos1[sx],end=bndtab$pos1[sx])),bam_files,sclip_opt$rl,list(what=c("pos","cigar")))
      cig1=matrix(nrow=length(bam1$pos),ncol=5)
      for(i in seq_along(bam1$pos)){
         cig1[i,] = SoftClipFromCigar(bam1$cigar[i])
      }
      readst=bam1$pos
      readed=readst+cig1[,3]-1
      readRanges=IRanges(start=as.vector(na.omit(readst)),end=as.vector(na.omit(readed)))
      if(bndtab$goRight1[sx]==1){
        leftbuffer=clippedbuffer; rightbuffer=mappedbuffer
      } else {
        leftbuffer=mappedbuffer; rightbuffer=clippedbuffer
      }
      nref1 = countOverlaps(IRanges(start=bndtab$pos1[sx]-leftbuffer,end=bndtab$pos1[sx]+rightbuffer),readRanges,type="within")
      bam2 = load_bam3(1,list(list(chrom=bndtab$chr2[sx],start=bndtab$pos2[sx],end=bndtab$pos2[sx])),bam_files,sclip_opt$rl,list(what=c("pos","cigar")))
      cig2=matrix(nrow=length(bam2$pos),ncol=5)
      for(i in seq_along(bam2$pos)){
          cig2[i,] = SoftClipFromCigar(bam2$cigar[i])
      }
      readst=bam2$pos
      readed =readst+cig2[,3]-1
      readRanges2=IRanges(start=as.vector(na.omit(readst)),end=as.vector(na.omit(readed)))
      if(bndtab$goRight2[sx]==1){
          leftbuffer=clippedbuffer; rightbuffer=mappedbuffer
      } else {
          leftbuffer=mappedbuffer; rightbuffer=clippedbuffer
      }
      nref2 = countOverlaps(IRanges(start=bndtab$pos2[sx]-leftbuffer,end=bndtab$pos2[sx]+rightbuffer),readRanges2,type="within")
      ###TODO, why wee need this here?
      cid=as.integer(bndtab$cid[sx]); mcid=as.integer(bndtab$mcid[sx])
      if(bndtab$EventType[sx] %in% c("INS")){
        scn1=scL[[cid]]$nreads; pos1=scL[[cid]]$position 
        scn2=scR[[mcid]]$nreads; pos2=scR[[mcid]]$position 
      } else {
        scn1=if(bndtab$goRight1[sx]==1) scL[[cid]]$nreads else scR[[cid]]$nreads
        pos1=if(bndtab$goRight1[sx]==1) scL[[cid]]$position else scR[[cid]]$position
        scn2=if(bndtab$goRight2[sx]==1) scL[[mcid]]$nreads else scR[[mcid]]$nreads
        pos2=if(bndtab$goRight2[sx]==1) scL[[mcid]]$position else scR[[mcid]]$position
      }
      nref[sx]=nref1+nref2
      nalt[sx]=scn1+scn2
      totreads[sx]=nalt[sx]+nref[sx]
      purity[sx]=nalt[sx]/totreads[sx]
    }
  }
  
  for(i in seq_len(nrow(bndtab))){ #assume single entry for a reciprocal pair; cleaning up inconsistencies 
    cid=as.integer(bndtab$cid[i]); mcid=as.integer(bndtab$mcid[i])
    tryCatch({
      if(bndtab$EventType[i] %in% c("INS")){
        scn1=scL[[cid]]$nreads; pos1=scL[[cid]]$position 
        scn2=scR[[mcid]]$nreads; pos2=scR[[mcid]]$position 
      } else {
        scn1=if(bndtab$goRight1[i]==1) scL[[cid]]$nreads else scR[[cid]]$nreads
        pos1=if(bndtab$goRight1[i]==1) scL[[cid]]$position else scR[[cid]]$position
        scn2=if(bndtab$goRight2[i]==1) scL[[mcid]]$nreads else scR[[mcid]]$nreads
        pos2=if(bndtab$goRight2[i]==1) scL[[mcid]]$position else scR[[mcid]]$position
      }
    }, error=function(e) { 
      cat("scL1:",bndtab$goRight1[i],"\n"); cat("scL2:",bndtab$goRight2[i],"\n")
      cat("cid:",cid,"\n"); cat("mcid:",mcid,"\n")
      cat("chk",pos1,"=",bndtab$pos1[i],"\n")
      cat("chk",pos2,"=",bndtab$pos2[i],"\n")
      cat("scn1",scn1,"scn2",scn2,"\n")
      scn1=NA; pos1=NA; scn2=NA; pos2=NA 
    })
    if(bndtab$EventType[i] %in% c("DEL","DUP","INV")) { #report full event
      sv_alt=bndtab$EventType[i]; 
      reg_sclip[[idx]]=list(sv_chrname=bndtab$chr1[i],sv_start_chr=bndtab$chr1[i],sv_start=bndtab$pos1[i],sv_end_chr=bndtab$chr2[i],sv_end=bndtab$pos2[i],sv_start_ci_low=-ci_sclip,sv_start_ci_high=ci_sclip,sv_end_ci_low=-ci_sclip,sv_end_ci_high=ci_sclip,sv_type=sv_alt,sv_mtd="sclip",sv_val=bndtab$FoldChange[i],sv_remap=TRUE,sv_filter=TRUE,sv_qual=6,sv_id="sclip",sv_cvg=FALSE,sv_imprecise=FALSE,sv_alt=sv_alt,sv_freq=purity[i],sv_dep=totreads[i])  
      reg_sclip[[idx]]$conf_mtd=rep(FALSE,4); names(reg_sclip[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_sclip[[idx]]$conf_val=rep(NA,4); reg_sclip[[idx]]$conf_p=rep(NA,4);
      reg_sclip[[idx]]$conf_type=rep("NA",4); reg_sclip[[idx]]$conf_freq=rep(NA,4); 
      reg_sclip[[idx]]$conf_st=rep(NA,4); reg_sclip[[idx]]$conf_ed=rep(NA,4);
      reg_sclip[[idx]]$conf_mtd[3]=TRUE;
      reg_sclip[[idx]]$conf_val[3]=nalt[i]
      reg_sclip[[idx]]$conf_type[3]=sv_alt;
      reg_sclip[[idx]]$conf_freq[3]=purity[i];
      reg_sclip[[idx]]$conf_st[3]=bndtab$pos1[i];
      reg_sclip[[idx]]$conf_ed[3]=bndtab$pos2[i];
      idx=idx+1
    } else if (bndtab$EventType[i] %in% c("INS")) { #report events
      sv_alt=bndtab$EventType[i]; 
      reg_sclip[[idx]]=list(sv_chrname=bndtab$chr1[i],sv_start_chr=bndtab$chr1[i],sv_start=bndtab$pos1[i],sv_end_chr=bndtab$chr2[i],sv_end=bndtab$pos2[i],sv_start_ci_low=-ci_sclip,sv_start_ci_high=ci_sclip,sv_end_ci_low=-ci_sclip,sv_end_ci_high=ci_sclip,sv_type=sv_alt,sv_mtd="sclip",sv_val=bndtab$FoldChange[i],sv_remap=TRUE,sv_filter=TRUE,sv_qual=6,sv_id="sclip",sv_cvg=FALSE,sv_imprecise=FALSE,sv_alt=sv_alt,sv_freq=purity[i],sv_dep=totreads[i])  
      reg_sclip[[idx]]$conf_mtd=rep(FALSE,4); names(reg_sclip[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_sclip[[idx]]$conf_val=rep(NA,4); reg_sclip[[idx]]$conf_p=rep(NA,4);
      reg_sclip[[idx]]$conf_type=rep("NA",4); reg_sclip[[idx]]$conf_freq=rep(NA,4)
      reg_sclip[[idx]]$conf_st=rep(NA,4); reg_sclip[[idx]]$conf_ed=rep(NA,4);
      reg_sclip[[idx]]$conf_mtd[3]=TRUE;
      reg_sclip[[idx]]$conf_val[3]=nalt[i]
      reg_sclip[[idx]]$conf_type[3]=sv_alt;
      reg_sclip[[idx]]$conf_freq[3]=purity[i];
      reg_sclip[[idx]]$conf_st[3]=bndtab$pos1[i];
      reg_sclip[[idx]]$conf_ed[3]=bndtab$pos2[i];
      idx=idx+1
    } else { #report breakpoint
      goRight1=bndtab$goRight1[i]; goRight2=bndtab$goRight2[i];
      if(is.na(bndtab$EventType[i])&&drop_na) { cat("gR1",bndtab$goRight1[i],"gR2",bndtab$goRight2[i],"\n"); next }
      #cat(bndtab$chr1[i],bndtab$chr1[i],"\n")
      ref_base1=toString(subseq(ref_seq[[bndtab$chr1[i]]],bndtab$pos1[i],bndtab$pos1[i]))
      ref_base2=toString(subseq(ref_seq[[bndtab$chr2[i]]],bndtab$pos2[i],bndtab$pos2[i]))
      if(goRight1 && !goRight2) {  #gR is whether reference base on right
        sv1_alt=paste("]",bndtab$chr2[[i]],":",bndtab$pos2[[i]],"]",ref_base1,sep="") 
        sv2_alt=paste(ref_base2,"[",bndtab$chr1[[i]],":",bndtab$pos1[[i]],"[",sep="") 
      } else if(!goRight1 && goRight2) {
        sv1_alt=paste(ref_base1,"[",bndtab$chr2[i],":",bndtab$pos2[i],"[",sep="")
        sv2_alt=paste("]",bndtab$chr1[i],":",bndtab$pos1[i],"]",ref_base2,sep="")
      } else if(goRight1 && goRight2) {
         sv1_alt=paste("[",bndtab$chr2[i],":",bndtab$pos2[i],"[",ref_base1,sep="") 
        sv2_alt=paste("[",bndtab$chr1[i],":",bndtab$pos1[i],"[",ref_base2,sep="") 
      } else { 
        sv1_alt=paste(ref_base1,"]",bndtab$chr2[i],":",bndtab$pos2[i],"]",sep="") 
        sv2_alt=paste(ref_base2,"]",bndtab$chr1[i],":",bndtab$pos1[i],"]",sep="") 
      }
      reg_sclip[[idx]]=list(sv_chrname=bndtab$chr1[i],sv_start_chr=bndtab$chr1[i],sv_start=bndtab$pos1[i],sv_end_chr=bndtab$chr1[i],sv_end=bndtab$pos1[i],sv_start_ci_low=-ci_sclip,sv_start_ci_high=ci_sclip,sv_end_ci_low=-ci_sclip,sv_end_ci_high=ci_sclip,sv_type="TRA",sv_mtd="sclip",sv_val=bndtab$FoldChange[i],sv_remap=TRUE,sv_filter=TRUE,sv_qual=6,sv_id="sclip",sv_cvg=FALSE,sv_imprecise=FALSE,sv_alt=sv1_alt,sv_freq=purity[i],sv_dep=totreads[i])  
      reg_sclip[[idx]]$conf_mtd=rep(FALSE,4); names(reg_sclip[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_sclip[[idx]]$conf_val=rep(NA,4); reg_sclip[[idx]]$conf_p=rep(NA,4);
      reg_sclip[[idx]]$conf_type=rep("NA",4); reg_sclip[[idx]]$conf_freq=rep(NA,4)
      reg_sclip[[idx]]$conf_st=rep(NA,4); reg_sclip[[idx]]$conf_ed=rep(NA,4);
      reg_sclip[[idx]]$conf_mtd[3]=TRUE;
      reg_sclip[[idx]]$conf_val[3]=nalt[i]
      reg_sclip[[idx]]$conf_type[3]=sv1_alt;
      reg_sclip[[idx]]$conf_freq[3]=purity[i];
      reg_sclip[[idx]]$conf_st[3]=bndtab$pos1[i];
      idx=idx+1
      reg_sclip[[idx]]=list(sv_chrname=bndtab$chr2[i],sv_start_chr=bndtab$chr2[i],sv_start=bndtab$pos2[i],sv_end_chr=bndtab$chr2[i],sv_end=bndtab$pos2[i],sv_start_ci_low=-ci_sclip,sv_start_ci_high=ci_sclip,sv_end_ci_low=-ci_sclip,sv_end_ci_high=ci_sclip,sv_type="TRA",sv_mtd="sclip",sv_val=bndtab$FoldChange[i],sv_remap=TRUE,sv_filter=TRUE,sv_qual=6,sv_id="sclip",sv_cvg=FALSE,sv_imprecise=FALSE,sv_alt=sv2_alt,sv_freq=purity[i],sv_dep=totreads[i])  
      reg_sclip[[idx]]$conf_mtd=rep(FALSE,4); names(reg_sclip[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_sclip[[idx]]$conf_val=rep(NA,4); reg_sclip[[idx]]$conf_p=rep(NA,4);
      reg_sclip[[idx]]$conf_type=rep("NA",4); reg_sclip[[idx]]$conf_freq=rep(NA,4)
      reg_sclip[[idx]]$conf_st=rep(NA,4); reg_sclip[[idx]]$conf_ed=rep(NA,4);
      reg_sclip[[idx]]$conf_mtd[3]=TRUE;
      reg_sclip[[idx]]$conf_val[3]=nalt[i]
      reg_sclip[[idx]]$conf_type[3]=sv2_alt;
      reg_sclip[[idx]]$conf_freq[3]=purity[i];
      reg_sclip[[idx]]$conf_st[3]=bndtab$pos2[i];
      idx=idx+1
    }
    #TODO: what does a quality mean here?
    #TODO: how to write partner information?
  }
  reg_sclip
}

#cat("loaded line 500...\n")

#' @title call_sclip
#'
#' @description \code{call_sclip} \cr
#' If called, \code{call_sclip} will return one list of sv_event structure
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_sclip(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_sclip=function(sclip_file,sclip_opt,bam_files,ref_seq,rle_cvg=NULL) {
  sclip_load=new.env(); reg_sclip=list(); rl=sclip_opt$rl
	cat(sclip_file)
  tryCatch( {load(sclip_file,envir=sclip_load)}, 
						error=function(e){cat("Warn:",sclip_file,"empty, nonexist or errornous formatted\n"); return(reg_sclip)})
  reg_sclip=list(); idx=1 #we now assume reg_sclip is done in sclip_call now!
  all_sclip=sclip_load$reg_sclip
  for(reg in all_sclip){ #overwrite proxy by sclip confirmation
    if(reg$sv_start_chr %in% sclip_opt$seq_name) {
      reg_sclip[[idx]]=reg; idx=idx+1
    }
  }
  throw=rep(FALSE,length(reg_sclip))
  for(idx in seq_along(reg_sclip)){
    throw[idx]=FALSE
    if(reg_sclip[[idx]]$sv_type %in% c("DUP","TANDEM.DUP") && abs(reg_sclip[[idx]]$sv_start-reg_sclip[[idx]]$sv_end)<rl)
      throw[idx]=TRUE #size filter
  }
  reg_sclip=reg_sclip[!throw]
  return(reg_sclip)
}

#' @title call_cvg
#'
#' @description \code{call_cvg} \cr
#' If called, \code{call_cvg} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_cvg(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_cvg=function(seqname,swan_file,swan_par,ovrd_file,swan_opt,coverage) { #return dup_track and del_track
  #print(swan_par)
  stepsize=swan_par$stepsize[1];start=swan_par$start[1];end=swan_par$end[1];
  width=swan_par$lw_width[1]; if(length(width)==0) width=1000
  stepsize=swan_par$stepsize[1];avg_win=width/stepsize
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score="cvg"; reg_cvg=list();
  if(!(score %in% names(swan_opt))) { return(reg_cvg) } 
  list[del_thresh,del_base]=smart_level(swan_opt[[score]]$del)
  list[dup_thresh,dup_base]=smart_level(swan_opt[[score]]$dup)
  list[sup_thresh,sup_base]=smart_level(swan_opt[[score]]$sup)
  list[gap_thresh,gap_base]=smart_level(swan_opt[[score]]$gap)
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score_track=as.numeric(filter(score_tmp[["cvg"]],rep(1./avg_win,avg_win)))
  cvg_track=score_tmp[["cvg"]]
  cvg_expect=mean(coverage)
  del_thresh_track=rep(del_thresh*cvg_expect,length(win_pos))
  dup_thresh_track=rep(dup_thresh*cvg_expect,length(win_pos))
  cvg_del_track=score_track<=del_thresh_track
  cvg_dup_track=score_track>=dup_thresh_track
  cvg_del_run=filter_calls(cvg_del_track,supnum=sup_thresh/stepsize,gaplim=gap_thresh/stepsize) 
  cvg_dup_run=filter_calls(cvg_dup_track,supnum=sup_thresh/stepsize,gaplim=gap_thresh/stepsize) 
  #print(Rle(score_track))
  #print(Rle(cvg_track))
  #print(Rle(cvg_expect))
  #print(Rle(del_thresh_track))
  #print(Rle(dup_thresh_track))
  #print(cvg_del_run)
  #print(cvg_dup_run)
  idx=1; isize_w=round(max(swan_par$isize,na.rm=T)/stepsize)
  for(ri in seq_len(nrun(cvg_del_run))){
    if(!is.na(runValue(cvg_del_run)[ri])&runValue(cvg_del_run)[ri]>0) {
      w_start=start(cvg_del_run)[ri];w_end=end(cvg_del_run)[ri]
      sv_cvg=if(max(cvg_track[w_start:w_end],0,na.rm=T)>=3*sum(coverage,na.rm=T)) TRUE else FALSE
      #reg_cvg[[idx]]=list(w_start=w_start,w_end=w_end,sv_cvg=sv_cvg)
      reg_cvg[[idx]]=list(sv_chrname=seqname,sv_start=win_pos[w_start],sv_end=win_pos[w_end],sv_start_ci_low=-1*width,sv_start_ci_high=width,sv_end_ci_low=-1*width,sv_end_ci_high=width, sv_type="DEL",sv_mtd="cvg",sv_val=mean(score_track[w_start:w_end]),sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_id="cvg",sv_cvg=sv_cvg,sv_imprecise=TRUE)
      reg_cvg[[idx]]$conf_mtd=rep(FALSE,4); names(reg_cvg[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_cvg[[idx]]$conf_val=rep(NA,4); reg_cvg[[idx]]$conf_p=rep(NA,4);
      reg_cvg[[idx]]$conf_type=rep("",4); reg_cvg[[idx]]$conf_freq=rep(NA,4)
      reg_cvg[[idx]]$conf_st=rep(NA,4); reg_cvg[[idx]]$conf_ed=rep(NA,4);
      idx=idx+1
    }
  }
  for(ri in seq_len(nrun(cvg_dup_run))){
    if(!is.na(runValue(cvg_dup_run)[ri])&runValue(cvg_dup_run)[ri]>0) {
      w_start=start(cvg_dup_run)[ri];w_end=end(cvg_dup_run)[ri]
      sv_cvg=if(max(cvg_track[(w_start-isize_w):(w_end+isize_w)],0,na.rm=T)>=3*sum(coverage,na.rm=T)) TRUE else FALSE
      reg_cvg[[idx]]=list(w_start=w_start,w_end=w_end,sv_cvg=sv_cvg)
      reg_cvg[[idx]]=list(sv_chrname=seqname,sv_start=win_pos[w_start],sv_end=win_pos[w_end],sv_start_ci_low=-width,sv_start_ci_high=width,sv_end_ci_low=-width,sv_end_ci_high=width, sv_type="DUP",sv_mtd="cvg",sv_val=mean(score_track[w_start:w_end]),sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_id="cvg",sv_cvg=sv_cvg,sv_imprecise=TRUE)
      reg_cvg[[idx]]$conf_mtd=rep(FALSE,4); names(reg_cvg[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
      reg_cvg[[idx]]$conf_val=rep(NA,4); reg_cvg[[idx]]$conf_p=rep(NA,4);
      reg_cvg[[idx]]$conf_type=rep("",4); reg_cvg[[idx]]$conf_freq=rep(NA,4)
      reg_cvg[[idx]]$conf_st=rep(NA,4); reg_cvg[[idx]]$conf_ed=rep(NA,4);
      idx=idx+1
    }
  }
  return(reg_cvg)
}

#' @title geno_sclip
#'
#' @description \code{geno_sclip} \cr
#' If called, \code{geno_sclip} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' geno_sclip(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
geno_sclip=function(si,sclip_set,bam_files,seqname,scan_par,sc1L,sc1R,r_lst){
  #list[si,sclip_set,bam_files,seqname,scan_par,sc1L,sc1R,sctab1L,sctab1R,r_lst]=list(3,sclip_set,case_files,seqname,swan_par[[2]],sc1L,sc1R,sctab1L,sctab1R,r_lst)
  #sv_start and sv_end is the exact break point
  reg=NULL
  if(length(sclip_set[[si]]$sv_target)!=0) { 
    #can we assume registered reads are from mutant others are from reference?
    reg=sclip_set[[si]]$sv_target; #idx=as.integer(reg$sv_id) #
    reg_st_br=reg$sv_start; reg_ed_br=reg$sv_end;
    reg_st_idx=which.min(c(sc1L[[si]]$position,sc1R[[si]]$position))
    reg_ed_idx=which.max(c(sc1L[[si]]$position,sc1R[[si]]$position))
    reg_st_sc=c(sc1L[[si]]$position,sc1R[[si]]$position)[reg_st_idx]
    reg_ed_sc=c(sc1L[[si]]$position,sc1R[[si]]$position)[reg_ed_idx]
    bam_sclip_st=DataFrame(); bam_sclip_ed=DataFrame()
    #cat(reg_st_sc,reg_ed_sc)
    what=c("pos","qwidth")
    #print(bam_files)
    for(bi in seq_len(length(bam_files))){ #bi=1
      seq_info=GRanges(seqname,IRanges(start=reg_st_sc-500,end=reg_ed_sc+500))
      #print(DataFrame(allFunction(seq_info,bam_files[bi],what=what))$qwidth)
      rl=median(DataFrame(allFunction(seq_info,bam_files[bi],what=what))$qwidth,na.rm=T)
      #print(c(reg_st_sc,rl))
      #print(c(reg_ed_sc,rl))
      #print(IRanges(start=reg_st_sc-rl,end=reg_st_sc))
      #print(IRanges(start=reg_ed_sc-rl,end=reg_ed_sc))
      seq_info_st=GRanges(seqname,IRanges(start=reg_st_sc-rl,end=reg_st_sc))
      seq_info_ed=GRanges(seqname,IRanges(start=reg_ed_sc-rl,end=reg_ed_sc))
      #print(seq_info_st); print(seq_info_ed)
      bam_tmp1=DataFrame(allFunction(seq_info_st,bam_files[bi],what=what))
      bam_sclip_st=IRanges::rbind(bam_sclip_st,bam_tmp1)
      bam_tmp2=DataFrame(allFunction(seq_info_ed,bam_files[bi],what=what))
      bam_sclip_ed=IRanges::rbind(bam_sclip_ed,bam_tmp2) 
    }
    reg_st=if(reg_st_sc==1) sc1L[[si]] else sc1R[[si]]
    reg_ed=if(reg_ed_sc==1) sc1R[[si]] else sc1L[[si]]
    reg_st_sc_n=reg_st$consmat$nseq[1]
    reg_ed_sc_n=reg_ed$consmat$nseq[1]
    reg_st_n=nrow(bam_sclip_st); reg_ed_n=nrow(bam_sclip_ed)
    lk_st=sapply(r_lst,function(x){dbinom(reg_st_sc_n,reg_st_n,x)})
    lk_ed=sapply(r_lst,function(x){dbinom(reg_ed_sc_n,reg_ed_n,x)})
    rj=which.max(c(lk_st,lk_ed)); rmax=r_lst[rj%%3]; reg$freq=rmax
  }
  reg
} 

#FUNC: determine threshold level from option input smartly
smart_level=function(opt){
  tryCatch({if(grepl("level",opt)) {
      opt_thresh=as.integer(substr(opt,nchar(opt),nchar(opt))); opt_num=FALSE
    } else {
      opt_thresh=as.numeric(opt); opt_num=TRUE 
    }}, error=function(e){ stop(paste("==Error: error parsing",opt,", possibly missing required opt?\n")) })
  list(opt_thresh,opt_num)
}

#FUNC: return track peaks with at least supnum window width after combining peaks within gaplim
filter_calls=function(merge_calls,supnum=supnum,gaplim=gaplim) {
  merge_run=Rle(merge_calls)
  while(length(which(runLength(merge_run)<supnum&!is.na(runValue(merge_run))&runValue(merge_run)))>0)
    runValue(merge_run)[which(runLength(merge_run)<supnum&!is.na(runValue(merge_run))&runValue(merge_run))]=FALSE
  while(length(which(runLength(merge_run)<gaplim&!is.na(runValue(merge_run))&!runValue(merge_run)))>0)
    runValue(merge_run)[which(runLength(merge_run)<gaplim&!is.na(runValue(merge_run))&!runValue(merge_run))]=TRUE
  return(merge_run)
}

#' @title call_lcd
#'
#' @description \code{call_lcd} \cr
#' If called, \code{call_seqcbs} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_seqcbs(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_lcd=function(seqname,swan_file,swan_par,ovrd_file,swan_opt,coverage) { 
  rl=swan_par$rl[1]; w=swan_par$w[1]; sdR=swan_par$isize_sdR[1]
  stepsize=swan_par$stepsize[1];start=swan_par$start[1];end=swan_par$end[1]; #bonferroni=0.05/(end-start)
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score="lCd"; reg_lcd=list()
  if(!(score %in% names(swan_opt))) { return(reg_lcd) }
  list[gap_thresh,gap_base]=smart_level(swan_opt[[score]]$gap)
  list[sup_thresh,sup_base]=smart_level(swan_opt[[score]]$sup)
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score_track=score_tmp[[score]]; cvg_track=score_tmp[["cvg"]]
  K=round(100000/(2*stepsize))*2+1
  Z=score_track
  sel=which(is.na(Z))
  Z[sel]=median(Z,na.rm=TRUE)
  Zmed = runmed(Z,k=K)
  ZAD=abs(Z-Zmed)
  Zmad = runmed(ZAD,k=K)
  deviation_track = (Z-Zmed)/pmax(Zmad,median(Zmad))
  thresh_track=thresh_score(score,idx,score_track,swan_par,ovrd_file,swan_opt[[score]]$method,swan_opt[[score]]$thresh,seqname)
  cat("---Info: score=",score,"threshold mean=",mean(thresh_track,na.rm=T),"\n")
  if(length(score_track)!=length(thresh_track)) stop("==Error: thresh and score track length must be equal!")
  tmp_track=score_track>thresh_track
  #supnum determine the smallest deletion can be trusted from lcd signal, min(isize)+3*min(isize_sdR)
  #gaplim determine the max random disruption in runs of lcd singal, bases 99% chance at least MIN.READS.DEL reads will drop into the gap
  #MIN.READS.DEL=5; PROB.AT.LEAST.MIN.READS.DEL=(1-bonferroni); 
  #for(size in seq(10,GAP.LIM.MAX,10)) if(qbinom(PROB.AT.LEAST.MIN.READS.DEL,size,sum(swan_par$lambda),lower.tail=F)>MIN.READS.DEL) { gaplim=size; break; }
  #supnum=qnorm(1-bonferroni,min(swan_par$isize),min(swan_par$isize_sdR)) #optimistically using min
  lcd_run=filter_calls(tmp_track,supnum=sup_thresh/stepsize,gaplim=gap_thresh/stepsize) 
  #connect calls gap<gaplim and report calls sup>supnum (too many flase positives, remove marginal first)
  #zero_len=runLength(lcd_run)[!runValue(lcd_run)&!is.na(runValue(lcd_run))]
  #one_len=runLength(lcd_run)[runValue(lcd_run)&!is.na(runValue(lcd_run))]
  idx=1; ci_lcd=rl; ci_w=ci_lcd/stepsize #isize_w=round(max(swan_par$isize,na.rm=T)/stepsize)
  for(ri in seq_len(nrun(lcd_run))){
    if(!is.na(runValue(lcd_run)[ri])&runValue(lcd_run)[ri]>0) {
      w_start=start(lcd_run)[ri];w_end=end(lcd_run)[ri]
      sv_cvg=mean(cvg_track[(w_start-ci_w):(w_end+ci_w)],na.rm=T) #store coverage
      reg_lcd[[idx]]=list(sv_chrname=seqname,sv_start_chr=seqname,sv_start=win_pos[w_start]+1.9*rl,sv_end_chr=seqname,sv_end=win_pos[w_end]-.9*rl,sv_start_ci_low=-ci_lcd,sv_start_ci_high=ci_lcd,sv_end_ci_low=-ci_lcd,sv_end_ci_high=ci_lcd,sv_type="DEL",sv_mtd="lcd",sv_val=max(score_track[w_start:w_end]),sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_stat1=max(deviation_track[w_start:w_end]),sv_stat2=-1,sv_id="lcd",sv_cvg=sv_cvg,sv_imprecise=TRUE)
      idx=idx+1
    }
  }
  throw=rep(FALSE,length(reg_lcd))
  for(idx in seq_along(reg_lcd)){
    if(abs(reg_lcd[[idx]]$sv_end-reg_lcd[[idx]]$sv_start)<sdR) throw[idx]=TRUE #size filter, throw if too small
    if(nna_and_geq(reg_lcd[[idx]]$sv_cvg,2*coverage)) throw[idx]=TRUE  #coverage filter, throw if in high coverage region
  }
  reg_lcd=reg_lcd[!throw]
  return(reg_lcd)
}

#' @title call_haf_har
#'
#' @description \code{call_haf_har} \cr
#' If called, \code{call_haf_har} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_haf_har(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_haf_har=function(seqname,swan_file,scan_par,ovrd_file,scan_opt,coverage){ #call based on forward hanning and reverse hanging pairs
  #list[swan_file,scan_par,ovrd_file,scan_opt,seqname]=list(swan_files[2],swan_par[[2]],ovrd_files[2],swan_opt,seqname)
  #list[swan_file,scan_par,ovrd_file,scan_opt,seqname]=list(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt,seqname)
  score1="HAF"; score2="HAR"; 
  stepsize=scan_par$stepsize[1];start=scan_par$start[1];end=scan_par$end[1];bonferroni=0.05/(end-start)
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score=paste(score1,score2,sep="+"); reg_haf_har=list()
  if(!(score %in% names(scan_opt))) { return(reg_haf_har) }
  list[gap_thresh,gap_base]=smart_level(scan_opt[[score]]$gap)
  list[sup_thresh,sup_base]=smart_level(scan_opt[[score]]$sup)
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score1_track=score_tmp[[score1]] #lDl
  score2_track=score_tmp[[score2]] #lDr
  cvg_track=score_tmp[["cvg"]]
  thresh1_track=thresh_score(score1,idx,score1_track,scan_par,ovrd_file,scan_opt[[score]]$method,scan_opt[[score]]$thresh,seqname)
  cat("--Info: score=",score1,"threshold mean=",mean(thresh1_track,na.rm=T),"\n")
  thresh2_track=thresh_score(score2,idx,score2_track,scan_par,ovrd_file,scan_opt[[score]]$method,scan_opt[[score]]$thresh,seqname)
  cat("--Info: score=",score2,"threshold mean=",mean(thresh2_track,na.rm=T),"\n")
  haf_track=score1_track>thresh1_track;
  har_track=score2_track>thresh2_track;
  haf_run=filter_calls(haf_track, supnum=sup_thresh/stepsize, gaplim=gap_thresh/stepsize)
  har_run=filter_calls(har_track, supnum=sup_thresh/stepsize, gaplim=gap_thresh/stepsize)
  haf_vec=IRanges::as.vector(haf_run); haf_vec[is.na(haf_vec)]=F; haf_run=Rle(haf_vec) #convert NA to FALSE, safe if only call for INS 
  har_vec=IRanges::as.vector(har_run); har_vec[is.na(har_vec)]=F; har_run=Rle(har_vec)
  haf_start=start(haf_run)[runValue(haf_run)]
  har_end=end(har_run)[runValue(har_run)]
  tele_dist=sapply(haf_start,function(x){  b=har_end[har_end>=x]; if(length(b)>0) b[1]-x else NA })
  tele_peak=sapply(haf_start,function(x){  b=har_end[har_end>=x]; if(length(b)>0) b[1] else NA })
  dist_distr=tele_dist
  if(length(dist_distr)>0){
    dist_peak=median(dist_distr,na.rm=T); #find the median in steps
    dist_sd=mad(dist_distr,na.rm=T) #find sd in steps by taking < 10*median
    dist_delta=(dist_peak-2*dist_sd)*stepsize #expected lower bound of peak distance in bp
  } else {
    dist_delta=0  #also in bp
  }
  cat("dist_delta=",dist_delta,"isize",scan_par$isize,"isize_sdR",scan_par$isize_sdR)
  if(dist_delta>(3*scan_par$isize+10*scan_par$isize_sdR) | dist_delta<3*scan_par$rl) #we have wired libraries, use read length based
    tele_thresh=3*scan_par$rl
  else
    tele_thresh=dist_delta
  cat("=Info: tele_thresh:",tele_thresh,"bp \n")
  reg_tmp=list(); idx=1; 
  for(i in seq_along(dist_distr))
    if(dist_distr[i]<=(tele_thresh/stepsize) & dist_distr[i]>0 & !is.na(dist_distr[i])) {
      reg_tmp[[idx]] = list(start=haf_start[i],end=tele_peak[i],tele=FALSE); idx=idx+1
    }  #in steps
  rl_w=round(max(scan_par$rl,na.rm=T)/stepsize)
  reg_tmp=lapply(reg_tmp,function(r) { nr=r;if(max(cvg_track[(r$start-rl_w):(r$end+rl_w)],0,na.rm=T)>=3*sum(coverage,na.rm=T)) nr$sv_cvg=TRUE else nr$sv_cvg=FALSE; return(nr) })
  ci_hax=max(scan_par$rl); idx=1; reg_haf_har=list()
  for(idx in seq_len(length(reg_tmp))){
    reg_haf_har[[idx]]=list(sv_chrname=seqname, sv_start_chr=seqname, sv_start=win_pos[reg_tmp[[idx]]$start]+100, sv_end_chr=seqname, sv_end=win_pos[reg_tmp[[idx]]$end], sv_start_ci_low=-ci_hax, sv_start_ci_high=ci_hax, sv_end_ci_low=-ci_hax, sv_end_ci_high=ci_hax, sv_type=if(reg_tmp[[idx]]$tele) "DEL" else "INS", sv_mtd="hax", sv_val=max(pmax(score1_track[reg_tmp[[idx]]$start:reg_tmp[[idx]]$end], score2_track[reg_tmp[[idx]]$start:reg_tmp[[idx]]$end])), sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_id="hax",sv_cvg=reg_tmp[[idx]]$sv_cvg, sv_imprecise=TRUE) #INS here is not reliable, maybe small del
  }
  reg_haf_har
}

#' @title call_ldr_ldl
#'
#' @description \code{call_ldr_ldl} \cr
#' If called, \code{call_ldr_ldl} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_ldr_ldl(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_ldr_ldl=function(seqname,swan_file,scan_par,ovrd_file,scan_opt,coverage){ 
  #seqname="2"; swan_file="example/example.lib1.2.swan.txt.gz";
  #scan_par=swan_par[[1]][[2]]; ovrd_file=NULL; scan_opt=swan_opt[[1]] 
  #coverage=mean_cvg1[[seqname]]
  #Thoughts of improving:
  #INV: if ldr starts early than ldl and ldl overlap with ldr for a INV length 
  rl=scan_par$rl[1]; is=scan_par$isize[1]; w=scan_par$w[1]
  score1="lDr"; score2="lDl";     #ldr signal leads ldl signal 
  stepsize=scan_par$stepsize[1];start=scan_par$start[1];end=scan_par$end[1];bonferroni=0.05/(end-start)
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score=paste(score1,score2,sep="+"); reg_ldr_ldl=list()
  if(!(score %in% names(scan_opt))) { return(reg_ldr_ldl) }
  #print(scan_opt[[score]])
  list[gap_thresh,gap_is_base]=smart_level(scan_opt[[score]]$gap)
  list[sup_thresh,sup_is_base]=smart_level(scan_opt[[score]]$sup)
  if(!is.null(scan_opt[[score]]$back)){
    list[back_thresh,back_is_base]=smart_level(scan_opt[[score]]$back)
  } else {
    back_thresh=2*w; back_is_base=TRUE
  }
  if(!is.null(scan_opt[[score]]$tele)){
    list[tele_thresh,tele_is_base]=smart_level(scan_opt[[score]]$tele)
  } else {
    tele_thresh=1.5*w; tele_is_base=TRUE
  }
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score1_track=score_tmp[[score1]] #lDr, step based
  score2_track=score_tmp[[score2]] #lDl, step based
  cvg_track=score_tmp[["cvg"]]
  thresh1_track=thresh_score(score1,idx,score1_track,scan_par,ovrd_file,scan_opt[[score]]$method,scan_opt[[score]]$thresh,seqname)
  cat("---Info: score=",score1,"threshold mean=",mean(thresh1_track,na.rm=T),"\n")
  thresh2_track=thresh_score(score2,idx,score2_track,scan_par,ovrd_file,scan_opt[[score]]$method,scan_opt[[score]]$thresh,seqname)
  cat("---Info: score=",score2,"threshold mean=",mean(thresh2_track,na.rm=T),"\n")
  ldr_track=score1_track>thresh1_track;
  ldl_track=score2_track>thresh2_track;
  ldr_run=filter_calls(ldr_track, supnum=sup_thresh/stepsize, gaplim=gap_thresh/stepsize)
  ldl_run=filter_calls(ldl_track, supnum=sup_thresh/stepsize, gaplim=gap_thresh/stepsize)
  ldr_vec=IRanges::as.vector(ldr_run);ldr_vec[is.na(ldr_vec)]=FALSE;ldr_run=Rle(ldr_vec) 
  ldl_vec=IRanges::as.vector(ldl_run);ldl_vec[is.na(ldl_vec)]=FALSE;ldl_run=Rle(ldl_vec)
  ldr_peak=unlist(sapply(seq(nrun(ldr_run)),function(x){if(runValue(ldr_run)[x]) which.max(score1_track[start(ldr_run)[x]:end(ldr_run)[x]])+start(ldr_run)[x]-1}))
  ldl_peak=unlist(sapply(seq(nrun(ldl_run)),function(x){if(runValue(ldl_run)[x]) which.max(score2_track[start(ldl_run)[x]:end(ldl_run)[x]])+start(ldl_run)[x]-1}))
  tele_dist=sapply(ldr_peak,function(x){b=ldl_peak[ldl_peak>=x-back_thresh]; if(length(b)>0) b[1]-x else NA })
  tele_peak=sapply(ldr_peak,function(x){b=ldl_peak[ldl_peak>=x-back_thresh]; if(length(b)>0) b[1] else NA })
  #save.image("single.RData")
  summary(tele_dist)
  cat("=Info: tele_thresh:",tele_thresh,"bp \n")
  reg_tmp=list(); idx=1; 
  for(i in seq_along(tele_dist))
    if(tele_dist[i]<=(tele_thresh/stepsize) & !is.na(tele_dist[i])) {
      reg_tmp[[idx]] = list(start=ldr_peak[i],end=tele_peak[i],tele=TRUE); idx=idx+1 #ldr is start, ldl is end
    }  #in steps
  rl_w=round(max(scan_par$rl,na.rm=T)/stepsize)
  ci_ldx=max(scan_par$rl); idx=1; reg_ldr_ldl=list()
  for(idx in seq_len(length(reg_tmp))){
    sv_type=ifelse(reg_tmp[[idx]]$start-reg_tmp[[idx]]$end>-w,"INS","DEL")
    sv_wstart=min(reg_tmp[[idx]]$start,reg_tmp[[idx]]$end)
    sv_wend=max(reg_tmp[[idx]]$start,reg_tmp[[idx]]$end)
    reg_ldr_ldl[[idx]]=list(sv_chrname=seqname,sv_start_chr=seqname,sv_start=win_pos[sv_wstart]+rl,sv_end_chr=seqname,sv_end=win_pos[sv_wend]+rl,sv_start_ci_low=-ci_ldx,sv_start_ci_high=ci_ldx,sv_end_ci_low=-ci_ldx, sv_end_ci_high=ci_ldx, sv_type=sv_type, sv_mtd="ldx", sv_val=max(pmax(score1_track[sv_wstart:sv_wend], score2_track[sv_wstart:sv_wend])), sv_remap=FALSE, sv_filter=NA, sv_qual=NA, sv_id="ldx",sv_cvg=mean(cvg_track[(sv_wstart-rl_w):(sv_wend+rl_w)],na.rm=T), sv_imprecise=TRUE) #INS here is not reliable, maybe small del
  }
  reg_ldr_ldl
}

#' @title call_del
#'
#' @description \code{call_del} \cr
#' If called, \code{call_del} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_del(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_del=function(seqname,swan_file,swan_par,swan_opt,rle_cvg=NULL) { 
  #list[swan_file,swan_parf,ovrd_file,swan_opt,seqname]=list(swan_files[2],swan_parfs[2],ovrd_files[2],swan_opt,seqname)
  #list[swan_file,swan_parf,ovrd_file,swan_opt,seqname]=list(swan_files[1],swan_parfs[1],ovrd_files[1],swan_opt,seqname)
  #swan_par=read.table(swan_parf,header=TRUE)
  stepsize=swan_par$stepsize[1];start=swan_par$start[1];end=swan_par$end[1]; #bonferroni=0.05/(end-start)
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score="del"; reg_del=NULL
  if(!(score %in% names(swan_opt))) { return(reg_del) }
  list[cvg_thresh,cvg_base]=smart_level(swan_opt[[score]]$cvg) #cvg=5
  list[sup_thresh,sup_base]=smart_level(swan_opt[[score]]$sup) #sup=50
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=tryCatch({read_com(swan_file,comment_char="#",colClasses="numeric")},error=function(e){cat("==Warn:",swan_file,"empty, broken or nonexist\n"); stop(); list(comment="",score_tmp=data.frame())})
  #list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score_track=score_tmp[[score]]; thresh=cvg_thresh; cvg_track=score_tmp[["cvg"]]; 
  tmp_track=score_track>=thresh
  del_run=filter_calls(tmp_track,supnum=sup_thresh/stepsize,gaplim=0) 
  reg_tmp=list(); idx=1; ci_del=stepsize
  for(ri in seq_len(nrun(del_run))){
    if(!is.na(runValue(del_run)[ri])&runValue(del_run)[ri]>0) {
      w_start=start(del_run)[ri];w_end=end(del_run)[ri]; sv_cvg=mean(cvg_track[w_start:w_end],na.rm=T)
      reg_tmp[[idx]]=list(w_start=w_start-1,w_end=w_end+1,sv_cvg=sv_cvg)
      idx=idx+1
    }
  }
  #we take median sdR of insert size as ci for lcd call
  for(idx in seq_len(length(reg_tmp))){
    #sv_dep=(runValue(rle_cvg[win_pos[reg_tmp[[idx]]$w_start]])+runValue(rle_cvg[win_pos[reg_tmp[[idx]]$w_end]]))/2
    sv_dep=mean(cvg_track[reg_tmp[[idx]]$w_start:reg_tmp[[idx]]$w_end])
    sv_val=max(score_track[reg_tmp[[idx]]$w_start:reg_tmp[[idx]]$w_end])
    reg_del[[idx]]=list(sv_chrname=seqname,sv_start_chr=seqname,sv_start=win_pos[reg_tmp[[idx]]$w_start],sv_end_chr=seqname,sv_end=win_pos[reg_tmp[[idx]]$w_end],sv_start_ci_low=-ci_del,sv_start_ci_high=ci_del,sv_end_ci_low=-ci_del,sv_end_ci_high=ci_del, sv_type="DEL",sv_mtd="del",sv_val=sv_val,sv_remap=TRUE,sv_filter=NA,sv_qual=NA,sv_id="del",sv_cvg=reg_tmp[[idx]]$sv_cvg,sv_imprecise=TRUE,sv_freq=sv_val/sv_dep,sv_dep=sv_dep)
    reg_del[[idx]]$conf_mtd=rep(FALSE,4); names(reg_del[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
    reg_del[[idx]]$conf_val=rep(NA,4); reg_del[[idx]]$conf_p=rep(NA,4); 
    reg_del[[idx]]$conf_type=rep("NA",4); reg_del[[idx]]$conf_freq=rep(NA,4)
    reg_del[[idx]]$conf_st=rep(NA,4); reg_del[[idx]]$conf_ed=rep(NA,4);
    reg_del[[idx]]$conf_type[3]="DEL"
    reg_del[[idx]]$conf_mtd[3]=TRUE
    reg_del[[idx]]$conf_freq[3]=sv_val/sv_dep
  }
  reg_del
}

#cat("loaded line 1000...\n")

#' @title call_ins
#'
#' @description \code{call_ins} \cr
#' If called, \code{call_ins} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_ins(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_ins=function(seqname,swan_file,swan_par,swan_opt,rle_cvg=NULL) { 
  #list[swan_file,swan_parf,ovrd_file,swan_opt,seqname]=list(swan_files[2],swan_parfs[2],ovrd_files[2],swan_opt,seqname)
  #list[swan_file,swan_parf,ovrd_file,swan_opt,seqname]=list(swan_files[1],swan_parfs[1],ovrd_files[1],swan_opt,seqname)
  #swan_par=read.table(swan_parf,header=TRUE)
  stepsize=swan_par$stepsize[1];start=swan_par$start[1];end=swan_par$end[1]; #bonferroni=0.05/(end-start)
  win_pos=seq(start,end,stepsize) #reproduce the scanning windows
  score="ins";reg_ins=NULL; 
  if(!(score %in% names(swan_opt))) { return(reg_ins) }
  list[cvg_thresh,cvg_base]=smart_level(swan_opt[[score]]$cvg)
  list[sup_thresh,sup_base]=smart_level(swan_opt[[score]]$sup)
  cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
  list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
  score_track=score_tmp[[score]]; thresh=cvg_thresh; cvg_track=score_tmp[["cvg"]];
  tmp_track=score_track>=thresh
  ins_run=filter_calls(tmp_track,supnum=sup_thresh/stepsize,gaplim=0) 
  reg_tmp=list(); idx=1; ci_ins=stepsize
  for(ri in seq_len(nrun(ins_run))){
    if(!is.na(runValue(ins_run)[ri])&runValue(ins_run)[ri]>0) {
      w_start=start(ins_run)[ri];w_end=end(ins_run)[ri]; sv_cvg=mean(cvg_track[w_start:w_end],na.rm=T)
      reg_tmp[[idx]]=list(w_start=w_start-1,w_end=w_end+1,sv_cvg=sv_cvg)
      idx=idx+1
    }
  }
  #cat("br0\n")
  #we take median sdR of insert size as ci for lcd call
  #print(reg_tmp)
  #print(length(cvg_track))
  #print(length(score_track))
  for(idx in seq_along(reg_tmp)){
    #cat(idx,"\n")
    #if(idx==2) print(reg_tmp[[idx]])
    #sv_dep=runValue(rle_cvg[win_pos[reg_tmp[[idx]]$w_start]])
    sv_dep=mean(cvg_track[reg_tmp[[idx]]$w_start:reg_tmp[[idx]]$w_end],na.rm=T)
    sv_val=max(score_track[reg_tmp[[idx]]$w_start:reg_tmp[[idx]]$w_end])
    sv_freq=sv_val/sv_dep
    reg_ins[[idx]]=list(sv_chrname=seqname,sv_start_chr=seqname,sv_start=win_pos[reg_tmp[[idx]]$w_start],sv_end_chr=seqname,sv_end=win_pos[reg_tmp[[idx]]$w_end],sv_start_ci_low=-ci_ins,sv_start_ci_high=ci_ins,sv_end_ci_low=-ci_ins,sv_end_ci_high=ci_ins, sv_type="INS",sv_mtd="ins",sv_val=sv_val,sv_remap=TRUE,sv_filter=NA,sv_qual=NA,sv_id="ins",sv_cvg=reg_tmp[[idx]]$sv_cvg,sv_imprecise=TRUE,sv_freq=sv_freq,sv_dep=sv_dep)
    reg_ins[[idx]]$conf_mtd=rep(FALSE,4); names(reg_ins[[idx]]$conf_mtd)=c("CV","HA","SC","ST")
    reg_ins[[idx]]$conf_val=rep(NA,4); reg_ins[[idx]]$conf_p=rep(NA,4); 
    reg_ins[[idx]]$conf_type=rep("",4); reg_ins[[idx]]$conf_freq=rep(NA,4)
    reg_ins[[idx]]$conf_st=rep(NA,4); reg_ins[[idx]]$conf_ed=rep(NA,4);
    reg_ins[[idx]]$conf_type[3]="INS"
    reg_ins[[idx]]$conf_mtd[3]=TRUE
    reg_ins[[idx]]$conf_freq[3]=sv_freq
  }
  #cat("br1\n")
  reg_ins
}

#' @title call_bigd
#'
#' @description \code{call_bigd} \cr
#' If called, \code{call_bigd} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_bigd(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_bigd=function(seqname,bigd_file,bigd_par,bigd_opt,swan_file=NULL,coverage=NULL) { 
  rl=bigd_par$rl[1]; stepsize=bigd_par$stepsize[1]; #swan_par is bigd_par
  if(!is.null(swan_file)) {
    cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
    list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
    cvg_track=score_tmp[["cvg"]]
  }
  bigd_call=list(); idx=1; bigd_lst=data.frame()
  minmpr=max(1,as.integer(bigd_opt$minmpr[1]),na.rm=T) 
  maxins=min(10000000,as.integer(bigd_opt$maxins[1]),na.rm=T)
  tmp_lst=tryCatch({read.table(bigd_file,header=T)},error=function(e){cat("==Warn:",bigd_file,"empty or nonexist\n"); data.frame()}) 
  ovlap_lst=c()
  if(nrow(tmp_lst)>0) {
    colnames(tmp_lst)=c("lstart","lend","rstart","rend","support")
    bigd_lst=tmp_lst[tmp_lst$support>minmpr & abs(tmp_lst$lstart-tmp_lst$rstart)<maxins,]
    #sv_start=round((bigd_lst$lstart+bigd_lst$lend)/2); 
    #sv_end=round((bigd_lst$rstart+bigd_lst$rend)/2); 
    sv_start=bigd_lst$lend
    sv_end=bigd_lst$rstart
    bigd_ranges=IRanges(start=pmin(sv_start,sv_end),end=pmax(sv_start,sv_end))
    bigd_ovlap=IRanges::as.matrix(IRanges::findOverlaps(bigd_ranges,bigd_ranges)) #query subject
    for(i in seq_len(nrow(bigd_ovlap))){ #assuming no tiling, if ovlap keep the small one
      if(bigd_ovlap[i,2]>bigd_ovlap[i,1]) ovlap_lst=c(ovlap_lst,bigd_ovlap[i,2])
    }
  }
  #we take lend as sv_start, rstart as sv_end
  for(bi in seq_len(nrow(bigd_lst))){
    if(!(bi %in% ovlap_lst)){
      #sv_start=round((bigd_lst[bi,]$lstart+bigd_lst[bi,]$lend)/2)
      #sv_end=round((bigd_lst[bi,]$rstart+bigd_lst[bi,]$rend)/2)
      sv_start=bigd_lst[bi,]$lend; sv_wstart=round(sv_start/stepsize)
      sv_end=bigd_lst[bi,]$rstart; sv_wend=round(sv_end/stepsize)
      if(abs(sv_start-sv_end)>maxins) next
      if(!is.null(swan_file)) bigd_cvg=mean(cvg_track[sv_wstart:sv_wend],na.rm=T) else bigd_cvg=-1;
      bigd_call[[idx]]=list(sv_start_chr=seqname, sv_start=sv_start-.4*rl, sv_end_chr=seqname, sv_end=sv_end+.5*rl, sv_chrname=seqname,sv_start_ci_low=bigd_lst[bi,]$lstart-sv_start, sv_start_ci_high=bigd_lst[bi,]$lend-sv_start, sv_end_ci_low=bigd_lst[bi,]$rstart-sv_end, sv_end_ci_high=bigd_lst[bi,]$rend-sv_end, sv_type="DEL", sv_mtd="bigd", sv_val=bigd_lst[bi,]$support, sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_id="bigd",sv_cvg=bigd_cvg,sv_imprecise=TRUE)
      idx=idx+1
    }#end if big_lst
  }#end bi
  bigd_call
}

#' @title call_disc
#'
#' @description \code{call_disc} \cr
#' If called, \code{call_disc} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' call_disc(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
call_disc=function(seqname,disc_file,disc_par,disc_opt,swan_file=NULL,coverage=NULL) { #NOTE this doesn't call translocation
  rl=disc_par$rl[1]; is=disc_par$isize[1]; stepsize=disc_par$stepsize[1]
  if(!is.null(swan_file)) {
    cat("reading in", swan_file, "if the program halt here, likely the file was truncated, need regenerate\n")
    list[comment,score_tmp]=read_com(swan_file,comment_char="#",colClasses="numeric")
    cvg_track=score_tmp[["cvg"]]
  }
  disc_call=list(); idx=1; disc_lst=data.frame()
  minmpr=max(1,as.integer(disc_opt$minmpr[1]),na.rm=T) 
  maxins=min(10000000,as.integer(disc_opt$maxins[1]),na.rm=T)
  tmp_lst=tryCatch({read.table(disc_file,header=T)},error=function(e){cat("==Warn:",disc_file,"empty or nonexist\n"); data.frame()}) 
  ovlap_lst=c()
  if(nrow(tmp_lst)>0) {
    colnames(tmp_lst)=c("lstart","lend","rstart","rend","support")
    disc_lst=tmp_lst[tmp_lst$support>minmpr & abs(tmp_lst$lstart-tmp_lst$rstart)<maxins,]
    sv_start=disc_lst$lstart
    sv_end=disc_lst$rend
    disc_ranges=IRanges(start=pmin(sv_start,sv_end),end=pmax(sv_start,sv_end))
    disc_ovlap=IRanges::as.matrix(IRanges::findOverlaps(disc_ranges,disc_ranges)) #query subject
    for(i in seq_len(nrow(disc_ovlap))){ #assuming no tiling, if ovlap keep the small one
      if(disc_ovlap[i,2]>disc_ovlap[i,1]) ovlap_lst=c(ovlap_lst,disc_ovlap[i,2])
    }
  }
  #we take mid of lstart and lend as sv_start, mid of rstart and rend as sv_end
  #should use lstart as start, rend as end here (haven't implemented)
  for(bi in seq_len(nrow(disc_lst))){
    if(!(bi %in% ovlap_lst)){
      sv_left=round((disc_lst[bi,]$lstart+disc_lst[bi,]$lend)/2)
      sv_right=round((disc_lst[bi,]$rstart+disc_lst[bi,]$rend)/2)
      if(sv_left>sv_right) { #left and right reversed
        sv_start=sv_right; sv_start_ci_low=disc_lst[bi,]$rstart-sv_right; sv_start_ci_high=disc_lst[bi,]$rend-sv_right
        sv_end=sv_left; sv_end_ci_low=disc_lst[bi,]$lstart-sv_left; sv_end_ci_high=disc_lst[bi,]$lend-sv_left
      } else {
        sv_start=sv_left; sv_start_ci_low=disc_lst[bi,]$lstart-sv_left; sv_start_ci_high=disc_lst[bi,]$lend-sv_left
        sv_end=sv_right; sv_end_ci_low=disc_lst[bi,]$rstart-sv_right; sv_end_ci_high=disc_lst[bi,]$rend-sv_right
      }
      if(abs(sv_start-sv_end)>maxins) next
      disc_call[[idx]]=list(sv_start_chr=seqname, sv_start=sv_start+.5*is-0.9*rl, sv_end_chr=seqname, sv_end=sv_end+.5*is+.9*rl, sv_chrname=seqname,sv_start_ci_low=sv_start_ci_low, sv_start_ci_high=sv_start_ci_high, sv_end_ci_low=sv_end_ci_low, sv_end_ci_high=sv_end_ci_high, sv_type="DUP", sv_mtd="disc", sv_val=disc_lst[bi,]$support, sv_remap=FALSE,sv_filter=NA,sv_qual=NA,sv_id="disc",sv_cvg=-1,sv_imprecise=TRUE)
      if(!is.null(swan_file)) {
	#cat("disc_region",round(disc_call[[idx]]$sv_start/stepsize),",",round(disc_call[[idx]]$sv_end/stepsize),"\n")
	#print(cvg_track[round(disc_call[[idx]]$sv_start/stepsize):round(disc_call[[idx]]$sv_end/stepsize)])
	disc_cvg=mean(cvg_track[round(disc_call[[idx]]$sv_start/stepsize):round(disc_call[[idx]]$sv_end/stepsize)],na.rm=TRUE)
	disc_call[[idx]]$sv_cvg=disc_cvg 
	#cat("disc_cvg",disc_cvg,"\n")
        if(na_less(disc_cvg,1.3*disc_opt$cvg[[seqname]])) { disc_call[[idx]]$sv_type="INV" } #
      }
      idx=idx+1
    }
  }
  disc_call
}

#FUNC: merge two tracks if peak1 in track1 is overlapping or within distance of tele_win from peak2 in track2 (ordered)
merge_track=function(call_track1,call_track2,tele_win=0){ #merging binary tracks of NA, 0, 1 within tele_win 
  #list[call_track1,call_track2,tele_win]=list(as.vector(ldl_run),as.vector(ldr_run),round(tele_thresh/stepsize))
  #list[call_track1,call_track2,tele_win]=list(call1,call2,0)
  call_track=call_track1+2*call_track2; merge_win=list(); idx=1; pos=1;
  call_run=Rle(call_track)
  call_value=runValue(call_run); call_length=runLength(call_run); 
  call_start=start(call_run); call_end=end(call_run); i=0
  #runValue(call_run)[!is.na(runValue(call_run))&runValue(call_run)>0]
  #call_value[!is.na(runValue(call_run))&runValue(call_run)>0]
  #runLength(call_run)[!is.na(runValue(call_run))&runValue(call_run)>0]
  #call_length[!is.na(runValue(call_run))&runValue(call_run)>0]
  #interestingly 201 and 231 dominates the pattern
  while(i<nrun(call_run)){ #could be NA,1,2,3
    i=i+1
    if(call_value[i]%in%c(3)){ #starting with 3|,3N,32,30
      if(i==nrun(call_run)) { #3|
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i],tele=FALSE);idx=idx+1;break
      } else if(call_value[i+1] %in% c(0,NA)) { #3N...,30...
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i],tele=FALSE);idx=idx+1;i=i+1;next
      } else if(call_value[i+1] %in% c(2)) { #32...
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i+1],tele=FALSE);idx=idx+1;i=i+1;next
      }
    } else if(call_value[i]%in%c(1)) { #starting with 1|,12|,13|,102,103,132,13N,130
      if(i==nrun(call_run)) { #1|
        break
      } else if (i==nrun(call_run)-1 & call_value[i+1] %in% c(2,3)){ #12|,13|
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i+1],tele=FALSE);idx=idx+1;break
      } else if (i<nrun(call_run)-1 & call_value[i+1] %in% c(3)){ #132...,130...,13N...
        if(call_value[i+2] %in% c(0,NA)){ #130...,13N...
          merge_win[[idx]]=list(start=call_start[i],end=call_end[i+1],tele=FALSE);idx=idx+1;i=i+2;next
        } else if (call_value[i+2] %in% c(2)){ #132...
          merge_win[[idx]]=list(start=call_start[i],end=call_end[i+2],tele=FALSE);idx=idx+1;i=i+2;next
        }
      } else if (i<nrun(call_run)-1 & call_value[i+1] %in% c(0) & call_value[i+2] %in% c(2) & call_length[i+1]<=tele_win & call_length[i+1]>0) { #102...
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i+2],tele=TRUE);idx=idx+1;i=i+2;next
      } else if (i<nrun(call_run)-1 & call_value[i+1] %in% c(2)) { #12..
        merge_win[[idx]]=list(start=call_start[i],end=call_end[i+1],tele=FALSE);idx=idx+1;i=i+1;next
      }
    }
  }
  merge_win
}

#' @title thresh_score
#'
#' @description \code{thresh_score} \cr
#' If called, \code{thresh_score} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' thresh_score(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
thresh_score=function(score,idx,score_track,swan_par,ovrd_file,method,thresh,seqname,min_signal=0.1) { #ovrd_file if present is first priority
  #idx=which(swan_opt$score==score)[1]; 
  #need to take care of track enabled or not here
  #print(length(score_track))
  start=swan_par$start[1]; end=swan_par$end[1]; stepsize=swan_par$stepsize[1]; 
  delta=max(swan_par$delta); n_wins=length(score_track)
  thresh_track=rep(0,n_wins) #assume n_wins is all the same
  ovrd_swan=NULL
  if(!is.null(ovrd_file))
    if(file.exists(ovrd_file)) 
      ovrd_swan=safe_read(ovrd_file,header=TRUE,stringsAsFactors=FALSE) 
  if(!is.null(ovrd_swan)&score%in%colnames(ovrd_swan)){ #override must have value for all regions
    ncol=dim(ovrd_swan)[2]
    set_colClass(ovrd_swan,c("character","integer","integer",rep("numeric",ncol-3)))
    #score_idx=which(colnames(ovrd_swan)=="lCd"|colnames(ovrd_swan)=="lcd")[1]
    sel=which(ovrd_swan$chrom==seqname)
    for(ri in sel){
      w_start=max(1,ceiling(((ovrd_swan$chromStart[ri]+1-start)+1)/stepsize))
      w_end=min(n_wins,ceiling(((ovrd_swan$chromEnd[ri]+1-start)+1)/stepsize))
      w_sel=w_start:w_end
      thresh_track[w_sel]=if(is.na(ovrd_swan[[score]][ri])) thresh_track[w_sel] else rep(ovrd_swan[[score]][ri],length(w_sel))
    }#end for(reg)
  } else {
    use_method=TRUE
    if(!method %in% c("value","theo","boot","robust","empr")) { 
      cat("=Warn: method not implemented, fall back to theo\n"); method="theo" }  
    thresh=thresh
    if(!grepl("level",thresh)){ # numeric
      thresh_value=as.numeric(thresh); 
      if(method=="theo"){ 
        use_method=FALSE
        thresh_track=rep(thresh_value,n_wins)
      }
    } else { #specified in significance level
      #thresh_value=as.integer(substr(thresh,nchar(thresh),nchar(thresh)+1)); 
			thresh_value=as.integer(gsub("[^0-9]","",thresh))
    }
    if(use_method) {
      #list[comment,score_track]=read_com(swan_file,comment_char="#",colClasses="numeric")
      alpha=alpha_level[thresh_value]
      if(method=="boot") {
        list[thresh_value,tmp1,tmp2]=
          noovlap_est(data=score_track,resample=1000,
                    buffer=round(delta*2/stepsize),alpha=alpha)
        thresh_track=rep(thresh_value,n_wins)
      } else if(method=="robust") {
        list[thresh_value,tmp1,tmp2]=
          robust_est(data=score_track,lvl=thresh_value)
        thresh_track=rep(thresh_value,n_wins)
      } else if(method=="empr"){
        thresh_value = max(min(thresh_value,20),2)
				cat("==Info: empr thresh_value=",thresh_value,"\n")
        K=round(100000/(2*stepsize))*2+1
        if(is.na(K)){ 
					cat("==WARNING: K is NA, stepsize=",stepsize,"\n")
          K=10000
        }
        Z=score_track
        sel=which(is.na(Z))
        Z[sel]=median(Z,na.rm=TRUE)
        Zmed = runmed(Z,k=K)
        ZAD=abs(Z-Zmed)
        Zmad = runmed(ZAD,k=K)
        thresh_track=Zmed+thresh_value*pmax(Zmad,median(Zmad))
        cat("=Info: Threshold calculated for ",score," using empirical moments, NMAD=",thresh_value,sep="")
       } else {
         if(method!="theo") cat("=Warn: method", method, "implemented, fall back to theo\n")
         list[thresh_value,tmp1,tmp2]=
        theo_thresh(score=score,alpha=alpha,scan_par=swan_par,verbose=F)
         thresh_track=rep(thresh_value,n_wins)
       }
     }#end if use_method 
   }#end if ovrd_file
  thresh_track[thresh_track<min_signal]=min_signal
  thresh_track
}

#FUNC: return NULL if reg in reg_lst1 overlap with any in reg_lst2 in their core part
reg_core = function(x,reg_lst1,reg_lst2,verbose=FALSE){ #run reg_lst1[x] against all reg in reg_lst2
  reg1=reg_lst1[[x]]; ovlap=FALSE; 
  if(verbose) cat("Processing",x,"-th region of",length(reg_lst1),"regions\n")
  for(ri in seq_len(length(reg_lst2))){
    reg2=reg_lst2[[ri]]
    #print("reg1"); print(reg1)
    #print("reg2"); print(reg2)
    if(ovlap_sv_end(reg1$sv_start_chr,reg1$sv_start,reg1$sv_start_ci_low,reg1$sv_start_ci_high,reg2$sv_start_chr,reg2$sv_start,reg2$sv_start_ci_low,reg2$sv_start_ci_high) || ovlap_sv_end(reg1$sv_end_chr,reg1$sv_end,reg1$sv_end_ci_low,reg1$sv_end_ci_high,reg2$sv_end_chr,reg2$sv_end,reg2$sv_end_ci_low,reg2$sv_end_ci_high)) { #double overlap -> remove
      #cat(reg1$sv_start+reg1$sv_start_ci_low,reg1$sv_start+reg1$sv_start_ci_high,reg2$sv_start+reg2$sv_start_ci_low,reg2$sv_start+reg2$sv_start_ci_high,"\n")
      #cat(reg1$sv_end+reg1$sv_end_ci_low,reg1$sv_end+reg1$sv_end_ci_high,reg2$sv_end+reg2$sv_end_ci_low,reg2$sv_end+reg2$sv_end_ci_high,"\n")
      ovlap=TRUE; break
    }
  }
  if(ovlap) NULL else reg1
}

#FUNC: return regs in reg_lst1 don't overlap with any in reg_lst2
reg_diff=function(reg_lst1,reg_lst2){ #remove region in reg_lst1 that also in reg_lst2 
  #save(list=c("reg_lst1","reg_lst2"),file=paste("t",".RData",sep=""))
  if(length(reg_lst1)==0) return(reg_lst1)
  reg_diff=lapply(seq_len(length(reg_lst1)),reg_core,reg_lst1,reg_lst2)
  reg_diff[sapply(reg_diff, is.null)] <- NULL
  reg_diff
}

#FUNC: format region as sv_event structure
prep_reg=function(ri,reg_lst) { #standardizing reg fields for loading bam file and confirmation
  reg=reg_lst[[ri]]
  if(is.null(reg$sv_start_chr) | is.null(reg$sv_end_chr)) { reg$sv_start_chr=reg$sv_chrname; reg$sv_end_chr=reg$sv_chrname }
  reg$sv_start_chr=as.character(reg$sv_start_chr)
  reg$sv_end_chr=as.character(reg$sv_end_chr)
  if(is.null(reg$sv_start) | is.null(reg$sv_end) | is.null(reg$sv_chrname)) { print(reg); return(NULL) }
  if(is.null(reg$imprecise)) reg$imprecise=TRUE
  if(is.null(reg$sv_val)) reg$sv_val=0
  if(is.null(reg$sv_stat1)) reg$sv_stat1=-1
  if(is.null(reg$sv_stat2)) reg$sv_stat2=-1
  if(is.null(reg$sv_cvg)) reg$sv_cvg=-1
  if( (reg$sv_start_chr==reg$sv_end_chr && reg$sv_start>reg$sv_end) 
         || reg$sv_start_chr > reg$sv_end_chr ) {  # by lexical order
    #either (chrE>chrS) or (chrE==chrS & posE>posS)
    #cat(paste("=Warn:",ri,"sv_start is to the right of sv_end, switching\n"))
    tmp=reg$sv_start; reg$sv_start=reg$sv_end; reg$sv_end=tmp;
    tmp=reg$sv_start_chr; reg$sv_start_chr=reg$sv_end_chr; reg$sv_end_chr=tmp;
    tmp=reg$sv_start_ci_low; reg$sv_start_ci_low=reg$sv_end_ci_low; reg$sv_end_ci_low=tmp;
    tmp=reg$sv_start_ci_high; reg$sv_start_ci_high=reg$sv_end_ci_high; reg$sv_end_ci_high=tmp;
  }
  if(reg$sv_start_ci_low>0) { reg$sv_start_ci_low=-100; cat(paste("=Warn:",ri,"sv_start_ci_low >0, set -100\n")) }
  if(reg$sv_start_ci_high<0) { reg$sv_start_ci_high=100; cat(paste("=Warn:",ri,"sv_start_ci_high <0, set 100\n")) }
  if(reg$sv_end_ci_low>0) { reg$sv_end_ci_low=-100; cat(paste("=Warn:",ri,"sv_end_ci_low >0, set -100\n")) }
  if(reg$sv_end_ci_high<0) { reg$sv_end_ci_high=100; cat(paste("=Warn:",ri,"sv_end_ci_high <0, set 100\n")) }
  reg$conf_mtd=rep(FALSE,4); names(reg$conf_mtd)=c("CV","HA","SC","ST")
  reg$conf_val=rep(NA,4); reg$conf_p=rep(NA,4); reg$conf_type=rep("",4); reg$conf_freq=rep(NA,4)
  reg$conf_st=rep(NA,4);reg$conf_ed=rep(NA,4); reg$sv_conf=FALSE; reg$sv_dep=NA; reg$sv_freq=NA; 
  reg
}
#INFO: structure of a region:
#sv_chrname, sv_start, sv_end, sv_start_ci, sv_end_ci 
#sv_start_chr, sv_end_chr
#conf_mtd, conf_val, conf_p, conf_type, conf_freq, conf_st, conf_ed
#sv_conf, sv_cvg, sv_dep, sv_freq

#FUNC: check region format, if not good, format region as sv_event structure
chk_reg=function(ri,reg_lst) { #standardizing reg fields for loading bam file and confirmation
  reg=reg_lst[[ri]]
  if( (reg$sv_start_chr==reg$sv_end_chr && reg$sv_start>reg$sv_end) || reg$sv_start_chr>reg$sv_end_chr) {
  #if(reg$sv_start>reg$sv_end && reg$sv_start_chr==reg$sv_end_chr) {
    cat(paste("=Warn:",ri,"sv_start is to the right of sv_end, switching\n"))
    tmp=reg$sv_start; reg$sv_start=reg$sv_end; reg$sv_end=tmp;
    tmp=reg$sv_start_chr; reg$sv_start_chr=reg$sv_end_chr; reg$sv_end_chr=tmp;
    tmp=reg$sv_start_ci_low; reg$sv_start_ci_low=reg$sv_end_ci_low; reg$sv_end_ci_low=tmp;
    tmp=reg$sv_start_ci_high; reg$sv_start_ci_high=reg$sv_end_ci_high; reg$sv_end_ci_high=tmp;
  }
  reg
}

#FUNC: print reg as readable
preg=function(i,reg_lst){
  reg=reg_lst[[i]]
  cat(paste(i,reg$sv_start_chr, reg$sv_start+reg$sv_start_ci_low, reg$sv_start+reg$sv_start_ci_high, 
            reg$sv_end_chr, reg$sv_end+reg$sv_end_ci_low, reg$sv_end+reg$sv_end_ci_high, reg$sv_mtd, reg$sv_type, "\n", sep="\t"))
}

#FUNC: find any self-overlapping regions within a list of regions
ovlap_reg=function(reg_first){
  #thie only overlap properly on the same chr
  reg_second=reg_first; shift=1
  ovlap=vector("list",length(reg_first))
  if(length(reg_first)==0|length(reg_second)==0) return(ovlap)
  first_start=sapply(reg_first,"[[","sv_start")+sapply(reg_first,"[[","sv_start_ci_low")
  first_end=sapply(reg_first,"[[","sv_end")+sapply(reg_first,"[[","sv_end_ci_high")
  first_chr=sapply(reg_first,"[[","sv_start_chr")
  second_start=sapply(reg_second,"[[","sv_start")+sapply(reg_second,"[[","sv_start_ci_low")
  second_end=sapply(reg_second,"[[","sv_end")+sapply(reg_second,"[[","sv_end_ci_high")
  second_chr=sapply(reg_second,"[[","sv_start_chr")
  for(i in seq_along(reg_first)) 
    if(na_unequal(reg_first[[i]]$sv_start_chr,reg_first[[i]]$sv_end_chr)){
      first_start[i]=max_chr_len+shift; first_end[i]=max_chr_len+shift #
      second_start[i]=max_chr_len+shift; second_end[i]=max_chr_len+shift #
      shift=shift+1
    }
  for(i in seq_along(reg_first)) 
    if(reg_first[[i]]$sv_mtd=="sclip"){ #use sclip for region ovlap is dangerous
      first_start[i]=max_chr_len+shift; first_end[i]=max_chr_len+shift #
      second_start[i]=max_chr_len+shift; second_end[i]=max_chr_len+shift #
      shift=shift+1
    }
  first_span=GRanges(seqnames=first_chr,IRanges(start=first_start,end=first_end))
  second_span=GRanges(seqnames=second_chr,IRanges(start=second_start,end=second_end))
  ovlap_mat=IRanges::as.matrix(IRanges::findOverlaps(first_span,second_span))
  for(i in seq(nrow(ovlap_mat)))
    ovlap[[ovlap_mat[i,1]]]=c(ovlap[[ovlap_mat[i,1]]],ovlap_mat[i,2])
  ovlap
}

#FUNC: find any both end-overlapping regions within a list of regions
ovlap_tip=function(reg_first){
  #we shall call it overlap only if start and end both overlap,
  #yet must ensure the first and second are sorted in the same order (see prep_reg)
  reg_second=reg_first
  start_ovlap=vector("list",length(reg_first))
  end_ovlap=vector("list",length(reg_first))
  if(length(reg_first)==0|length(reg_second)==0) return(ovlap)
  first_start_chr=sapply(reg_first,"[[","sv_start_chr")
  first_start_lci=sapply(reg_first,"[[","sv_start")+sapply(reg_first,"[[","sv_start_ci_low")
  first_start_hci=sapply(reg_first,"[[","sv_start")+sapply(reg_first,"[[","sv_start_ci_high")
  first_end_chr=sapply(reg_first,"[[","sv_end_chr")
  first_end_lci=sapply(reg_first,"[[","sv_end")+sapply(reg_first,"[[","sv_end_ci_low")
  first_end_hci=sapply(reg_first,"[[","sv_end")+sapply(reg_first,"[[","sv_end_ci_high")
  second_start_chr=sapply(reg_second,"[[","sv_start_chr")
  second_start_lci=sapply(reg_second,"[[","sv_start")+sapply(reg_second,"[[","sv_start_ci_low")
  second_start_hci=sapply(reg_second,"[[","sv_start")+sapply(reg_second,"[[","sv_start_ci_high")
  second_end_chr=sapply(reg_second,"[[","sv_end_chr")
  second_end_lci=sapply(reg_second,"[[","sv_end")+sapply(reg_second,"[[","sv_end_ci_low")
  second_end_hci=sapply(reg_second,"[[","sv_end")+sapply(reg_second,"[[","sv_end_ci_high")
  first_start=GRanges(seqnames=first_start_chr,IRanges(start=first_start_lci,end=first_start_hci))
  first_end=GRanges(seqnames=first_end_chr,IRanges(start=first_end_lci,end=first_end_hci))
  second_start=GRanges(seqnames=second_start_chr,IRanges(start=second_start_lci,end=second_start_hci))
  second_end=GRanges(seqnames=second_end_chr,IRanges(start=second_end_lci,end=second_end_hci))
  start_ovlap_mat=IRanges::as.matrix(IRanges::findOverlaps(first_start,second_start))
  end_ovlap_mat=IRanges::as.matrix(IRanges::findOverlaps(first_end,second_end))
  for(i in seq(nrow(start_ovlap_mat)))
    start_ovlap[[start_ovlap_mat[i,1]]]=c(start_ovlap[[start_ovlap_mat[i,1]]],start_ovlap_mat[i,2])
  for(i in seq(nrow(end_ovlap_mat)))
    end_ovlap[[end_ovlap_mat[i,1]]]=c(end_ovlap[[end_ovlap_mat[i,1]]],end_ovlap_mat[i,2])
  list(start_ovlap,end_ovlap)
}

#cat("loaded line 1500...\n")

#FUNC: test if reg size is within limit
within_limit=function(reg,len,limit){
  within=abs((reg$sv_end-reg$sv_start)-len)<limit
  if(!is.na(within)) within else FALSE
}

#val1 possibly contains na
#FUNC: return True if and only val1 is ans, return False if val1 is na
na_equal=function(val1,ans){
  chk_na=is.na(val1);
  if(!chk_na) val1==ans else FALSE
}
#FUNC: return True if val1 != ans, return False if val1 is na 
na_unequal=function(val1,ans){
  chk_na=is.na(val1); chk_type=ans!=val1; 
  if(!chk_na) chk_type else FALSE
}
#FUNC: return True if ans >= val1, return False if val1 is na 
na_less=function(val1,ans){
  chk_na=is.na(val1); chk_type=ans>=val1; 
  if(!chk_na) chk_type else FALSE
}
#FUNC: return True if val1 <= preset, return False if val1 is na 
nna_and_leq=function(val1,ans){
  chk_na=is.na(val1); chk_type=ans>=val1; 
  if(!chk_na) chk_type else FALSE
}
#FUNC: return True if val1 >= preset, return False if val1 is na 
nna_and_geq=function(val1,ans){
  chk_na=is.na(val1); chk_type=ans<=val1; 
  if(!chk_na) chk_type else FALSE
}

#TODO: further better this important code
#' @title dedup_reg
#'
#' @description \code{dedup_reg} \cr
#' If called, \code{dedup_reg} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' dedup_reg(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
dedup_reg=function(reg_lst,strict=TRUE,limit=1000,biglimit=100000) {
  if(strict) output_ldx=FALSE else output_ldx=TRUE
  if(strict) output_sclip=FALSE else output_sclip=TRUE
  list[start_ovlap,end_ovlap]=ovlap_tip(reg_lst); #print(length(start_ovlap)); print(length(end_ovlap));
  span_ovlap=ovlap_reg(reg_lst); #print(length(span_ovlap))
  reg_dedup=rep(list(NULL), length(reg_lst)); done_idx=c()
  for(idx in seq_along(reg_lst)){
    if(idx %in% done_idx) {
      next
    } else { #this search will form a closed clique eventually
      #find out which regions to combine
      tmp_todo_idx=unique(c(start_ovlap[[idx]],end_ovlap[[idx]],span_ovlap[[idx]])); tmp_done_idx=c()
      while(length(tmp_todo_idx)!=length(tmp_done_idx)){ #stop when newly add todo is already in done
        tmp_add_idx=c()
        for(jdx in tmp_todo_idx[!tmp_todo_idx %in% tmp_done_idx]){
          tmp_add_idx=unique(c(tmp_add_idx,start_ovlap[[jdx]],end_ovlap[[jdx]],span_ovlap[[jdx]]))
        }
        tmp_done_idx=tmp_todo_idx #keep old todo in done
        tmp_todo_idx=unique(c(tmp_done_idx,tmp_add_idx)) #add new todo and 
      }
      done_idx=c(done_idx,tmp_done_idx)
      #rules for how to combine
      #region clique to combine is in tmp_todo_idx
      tmp_todo_mtd=sapply(reg_lst[tmp_todo_idx],"[[","sv_mtd")
      tmp_todo_type=sapply(reg_lst[tmp_todo_idx],"[[","sv_type")
      tmp_todo_chr1=sapply(reg_lst[tmp_todo_idx],"[[","sv_start_chr")
      tmp_todo_start=sapply(reg_lst[tmp_todo_idx],"[[","sv_start")
      tmp_todo_chr2=sapply(reg_lst[tmp_todo_idx],"[[","sv_end_chr")
      tmp_todo_end=sapply(reg_lst[tmp_todo_idx],"[[","sv_end")
      #sv_mtd: lcd, ldx, hax, cvg, disc, bigd, sclip, cbs, del, ins
      #sv_type: DEL, INS, INV, DUP, TRP
      done=FALSE
      #0. if see one and only one evidence:
      #   0.1 if sclip yet type=NA, no output of this, 
      #   0.2 otherwise output this and continue
      if(length(tmp_todo_idx)==1){
        kdx = tmp_todo_idx; 
        if( reg_lst[[kdx]]$sv_mtd=="ldx" && !output_ldx ) {
          done=TRUE; if(done) next; #ignoring singleton ldx-INS cases 
        } else if( reg_lst[[kdx]]$sv_mtd=="sclip" && reg_lst[[kdx]]$sv_type=="NA" && !output_sclip) {
          done=TRUE; if(done) next;  #ignoring singleton sclip-NA cases
        } else {
          reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; if(done) next; #whatever it is output
        }
        if(done) next;
      }
      #   0.3 if only ldx in the evidence
      if( all(tmp_todo_mtd==c("ldx")) && !output_ldx ) { 
        done=TRUE; if(done) next; #ignoring multiple ldx-INS cases
      }
      #1. if we see sclip in evidence
      #    1.1 if sclip or other DEL evidence: bigd, lcd, etc, 
      #        we consider this is a deletion and output the most probable sclip, bigd, lcd (DEL)
      #    1.2 else if sclip or other INV, TRP, DUP evidence: 
      #        we consider this is a INV, TRP, DUP in order and output
      #    1.3 else if sclip and INS evidence:
      #        we consider this is a INS
      #    1.4 else if no other call mtd than ldx:
      #        we discard it
      #    1.5 else we it will be passed to following classificaitons:
      if("sclip" %in% tmp_todo_mtd){ #
        sclip_idx = tmp_todo_idx[tmp_todo_mtd=="sclip"]
        sclip_type = tmp_todo_type[tmp_todo_mtd=="sclip"]
        del_idx = tmp_todo_idx[which(tmp_todo_type=="DEL")]
        del_mtd = tmp_todo_mtd[which(tmp_todo_type=="DEL")]
        del_len = (tmp_todo_end-tmp_todo_start)[tmp_todo_type=="DEL"]
        del_type = tmp_todo_type[which(tmp_todo_type=="DEL")]
        good_type = tmp_todo_type[tmp_todo_type!="NA"]
        good_type_idx = tmp_todo_idx[tmp_todo_type!="NA"]
        if(length(del_idx)>0) { #first consider deletion
          if("sclip" %in% del_type){ #let's use sclip
            kdx=del_idx[del_mtd=="sclip"][which.max(del_len[del_mtd=="sclip"])]
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          } else { #otherwise use the one with biggest span
            kdx=del_idx[which.max(del_len)]; 
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          }
        } else if(length(good_type)>0) { # trust INV>TRP>DUP>INS
          #controlled verbose types: DEL,INS,INV,TRP,DUP
          inv_idx=good_type_idx[good_type=="INV"]
          tra_idx=good_type_idx[good_type=="TRA"]
          dup_idx=good_type_idx[good_type=="DUP"]
          ins_idx=good_type_idx[good_type=="INS"]
          if(length(inv_idx>1)) { #let's choose the first one for simplicity
            kdx=inv_idx[1];
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          } else if(length(tra_idx>1)) {
            kdx=tra_idx[1];
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          } else if(length(dup_idx>1)) {
            kdx=dup_idx[1];
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          } else if(length(ins_idx>1)) {
            kdx=ins_idx[1];
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
          }
        } else { #else combinations mingled with sclip
          non_sclip_idx=tmp_todo_idx[tmp_todo_mtd!="sclip"]
          if(length(non_sclip_idx)==0) done=TRUE; #cat("add", kdx, "\n");
        }
        if(done) next;
      }
      #2. if see more than two of one method of evidence (but no sclip)
      #   2.1 we find the first with a sv type and we output
      #   2.2 if we already reach the last one, we have to output
      if(length(unique(tmp_todo_mtd))==1){ 
        for(i in seq_along(tmp_todo_idx)){
          kdx = tmp_todo_idx[i]
          if(na_unequal(tmp_todo_mtd[i],"NA")) { 
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; } #cat("add", kdx, "\n");
          if(i == length(tmp_todo_idx)) { 
            reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; } #cat("add", kdx, "\n");
        }
        if(done) next;
      }
      #3. if we see bigd in evidence with others, trust bigd
      #    3.1 find the max span bigd and output  
      if("bigd" %in% tmp_todo_mtd){ #dedup DEL
        bigd_idx = tmp_todo_idx[which(tmp_todo_mtd=="bigd")]
        kdx = bigd_idx[1]
        bigd_len = reg_lst[[kdx]]$sv_end-reg_lst[[kdx]]$sv_start
        for(jdx in bigd_idx[-1]) 
          if(reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start>bigd_len) 
            { kdx=jdx; bigd_len=reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start }
        reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
        if(done) next;
      }
      #4. if we see lcd in evidence with others 
      #    3.1 if we additionally see disc, find the max length of (disc and lcd), output this
      #    3.2 else we use lcd and output this
      if("lcd" %in% tmp_todo_mtd){ #dedup DEL
        lcd_idx = tmp_todo_idx[which(tmp_todo_mtd=="lcd")]; kdx = lcd_idx[1]
        lcd_len = reg_lst[[kdx]]$sv_end-reg_lst[[kdx]]$sv_start
        for(jdx in lcd_idx[-1]) 
          if(reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start>lcd_len) 
            { kdx=jdx; lcd_len=reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start }
        if("disc" %in% tmp_todo_mtd) {
          disc_idx = tmp_todo_idx[which(tmp_todo_mtd=="disc")]
          for(jdx in disc_idx) 
            if(reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start>lcd_len)
              { kdx=jdx; lcd_len=reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start }
          reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
        } else { #trust lcd if we see del, ins, cvg, ldx, hax and seqcbs etc.
          reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; #cat("add", kdx, "\n");
        }
        if(done) next;
      }
      #5. if we see disc and ldx together (no bigd, no lcd, no sclip)
      #   5.1 we find the max_len(disc, ldx) and output this. If INV is the sv_type, we set it to DUP 
      if("disc" %in% tmp_todo_mtd & "ldx" %in% tmp_todo_mtd){
        misc_idx = tmp_todo_idx[which(tmp_todo_mtd=="disc" | tmp_todo_mtd=="ldx")]; kdx = misc_idx[1]
        misc_len = reg_lst[[kdx]]$sv_end-reg_lst[[kdx]]$sv_start
        for(jdx in misc_idx[-1]) #exclude first
          if(reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start>misc_len)
            { kdx=jdx; misc_len=reg_lst[[jdx]]$sv_end-reg_lst[[jdx]]$sv_start }
        reg_dedup[[kdx]]=reg_lst[[kdx]]; done=TRUE; 
        #if(reg_dedup[[kdx]]$sv_type == "INV") reg_dedup[[kdx]]$sv_type="DUP"; done = TRUE; #correct type #cat("add", kdx, "\n");
        if(done) next;
      }
      #6. if we see ldx and ins/del together (no bigd, no lcd, no sclip)
      #   6.1 If ldx and del together, we set it to DEL 
      #   6.2 If ldx and del together, we set it to INS
      if("del" %in% tmp_todo_mtd & "ldx" %in% tmp_todo_mtd){
        del_idx = tmp_todo_idx[which(tmp_todo_mtd=="del")]; reg_dedup[[del_idx]] = reg_lst[[del_idx]]; done = TRUE; if(done) next;
      }
      if("ins" %in% tmp_todo_mtd & "ldx" %in% tmp_todo_mtd){
        ins_idx = tmp_todo_idx[which(tmp_todo_mtd=="ins")]; reg_dedup[[ins_idx]] = reg_lst[[ins_idx]]; done = TRUE; if(done) next;
      }
    } #end of process new_idx not in done_index
    #but if we still reach here without being processed, we see something unexpected, let's complain and fix
    cat("dropping index:",paste(tmp_todo_idx,collapse=","),"\n")
    cat("dropping methods:",paste(tmp_todo_mtd,collapse=","),"\n")
    cat(paste("above are not considerred occasions:",paste(tmp_todo_mtd,collapse=","),". If needed, please use raw.bed file"))
  }
  reg_dedup
}

#TODO: further better this important code
#' @title dehotspot_reg
#'
#' @description \code{dehotspot_reg} \cr
#' If called, \code{dehotspot_reg} will return two list of sv_event structure
#' one for good seqcbs region and one for seqcbs region other than good
#'
#' @param seqname: chr name
#' @param seqcbs_file: main output from seqcbs scan
#' @param seqcbs_parf: parameter file for seqcbs
#' @param seqdbs_opt: passed in options for seqcbs
#' @param rg_files: bam files possibly read group separated
#'
#' @note
#'
#' @examples
#' dehotspot_reg(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
dehotspot_reg=function(chr,bed,DIST.TO.NEXT.IN.HOTSPOT=1000,HOTSPOT.CLUSTER.SIZE=5){ #assume input is a dataframe with $CHROM defined
    subidx=which(bed$CHROM==chr)
    if(length(subidx)==0) return(subidx) #no idx to return
    toNext=rep(NA,length(subidx)-1)
    for(i in seq_along(subidx[-1])){
        toNext[i] = bed$START[subidx[i]+1]-bed$START[subidx[i]]
    }
    runst=0; runend=-1
    inHotspot=rep(FALSE,length(subidx))
    for(i in seq_along(subidx[-1])){
      if(toNext[i]<DIST.TO.NEXT.IN.HOTSPOT && toNext[i]>=0){
        if(runst==0) runst=i+1;
        runend=i+1
        #cat(i+1,"toNext=",toNext[i],": part of run, st=",runst,", end=",runend,"\n")
      } else {
        if(runend-runst>HOTSPOT.CLUSTER.SIZE) inHotspot[runst:runend]=TRUE
        #cat(i+1,"toNext=",toNext[i],": not in run, st=",runst,", end=",runend,"\n")
        runst=0; runend=-1
      }
    }
  subidx[!inHotspot]
}

#FUNC: load bam file according to a list of bed regions
load_bam3=function(bi,bed_lst,bam_files,bam_adj,load_opt){ 
  bam=DataFrame()
  for(i in seq_along(bam_files)){
    bf=bam_files[i]
    seq_info=GRanges(bed_lst[[bi]]$chrom,IRanges(start=bed_lst[[bi]]$start-bam_adj[i],end=bed_lst[[bi]]$end+bam_adj[i]))
    bam=IRanges::rbind(bam,DataFrame(allFunction(seq_info,bf,what=load_opt$what)))
  }
  #print(ri);
  bam
}#load all reads from sv_start+sv_start_ci_low,sv_end+sv_end_ci_high 

#FUNC: load bam file according to ri-th of reg_lst
load_bam2=function(ri,reg_lst,bam_files,load_opt){ 
  reg=reg_lst[[ri]] # c("pos","mpos","rname","mrnm","flag")
  #print(reg_lst[[ri]]$sv_start); print(reg$sv_start)
  sv_range=c(max(1,reg$sv_start+reg$sv_start_ci_low),reg$sv_end+reg$sv_end_ci_high)
  seq_info=GRanges(load_opt$seqname,IRanges(start=sv_range[1],end=sv_range[2]))
  bam=DataFrame()
  for(bf in bam_files){
    bam=IRanges::rbind(bam,DataFrame(allFunction(seq_info,bf,what=load_opt$what)))
  }
  #print(ri);
  bam
}#load all reads from sv_start+sv_start_ci_low,sv_end+sv_end_ci_high 
if(FALSE){
  bam_file=c("test/test_multi.lib1.bam","test/test_multi.lib2.bam")
  reg_lst=list(list(sv_start=70000000,sv_start_ci_low=-1000,sv_end=70500000,sv_end_ci_high=1000))
  load_opt=list(seqname="11",what=c("rname","mrnm","flag","isize"))
  bam=load_bam2(1,reg_lst,bam_file,load_opt)
}

#FUNC: 
reg_isize_alert=function(reg_lst,bam_files,isize_min,isize_max,seqname,flags){
  #alert false positive if not any promising read exceed min and non exceeding max
  load_opt=list(seqname=seqname,what=c("rname","mrnm","flag","isize"))
  bam_lst=lapply(seq_along(reg_lst),load_bam2,reg_lst,bam_files,load_opt)
  alert_idx=c()
  for(i in seq_along(reg_lst)) {
    alert=any(bam_lst[[i]]$flag %in% flags & abs(bam_lst[[i]]$isize)>isize_min & abs(bam_lst[[i]]$isize)<isize_max & bam_lst[[i]]$mrnm==seqname)
    if(!alert) alert_idx=c(alert_idx,i)
  }
  alert_idx
}
if(FALSE){
  bam_files=c("test/test_multi.lib1.bam","test/test_multi.lib2.bam")
  reg_lst=list(list(sv_start=70000000,sv_start_ci_low=-1000,sv_end=70500000,sv_end_ci_high=1000))
  reg_isize_alert(reg_lst,bam_file,1000,100000,"11",c(conc_flags,impp_flags,disc_flags))
  bam_lst[[1]][1,4]=5000
  reg_isize_alert(reg_lst,bam_file,1000,100000,"11",c(conc_flags,impp_flags,disc_flags))
  bam_lst[[1]][1,4]=500
}

#FUNC: load bam file according to ri-th of reg_lst
load_bam=function(ri,reg_lst,bam_files,swan_par,load_opt){
  #list[ri,reg_lst,bam_files,swan_par,load_opt]=list(1,reg_lst,bam_files,swan_par,load_opt)
  #conf_cvg: spX region(pos,qwidth) and spY region(pos,qwidth) / spX 0.5*buffer_len(pos,qwidth)
  #conf_hang: spX buffer_isize(pos,isize,hang) and spY region+med_buffer(pos,isize,hang)
  #conf_soft: spX buffer_isize(pos,isize,soft) and spY buffer_sci(pos,isize)
  #conf_strad: spX buffer_isize(pos,isize)+buffer_eci(pos,isize) and spY buffer_sci(pos,isize)+buffer_eci(pos,isize)
  n_sp=length(bam_files); reg=reg_lst[[ri]]; bam=list()
  if(!is.null(swan_par)) buf_isize=lapply(swan_par,function(x){max(x$isize)+3*max(x$isize_sdR)}) else as.list(rep(1000,n_sp))
  sv_outter=c(max(1,reg$sv_start+reg$sv_start_ci_low),reg$sv_end+reg$sv_end_ci_high); 
  sv_inner=c(reg$sv_start+reg$sv_start_ci_high,reg$sv_end+reg$sv_end_ci_low);
  st0=sv_outter[1];ed0=sv_outter[2];buf_svlen=sv_outter[2]-sv_outter[1] #buffer_ci is 2*ci
  bam$sv_outter=sv_outter; bam$sv_inner=sv_inner
  what=c("pos","qwidth"); seq_info=GRanges(reg$sv_start_chr,IRanges(start=st0,end=ed0)); 
  bam$bam_buf_svlen=buf_svlen; bam$bam_cvg=list(DataFrame())
  #load bam_cvg only if sv on one chr
  if(reg$sv_start_chr==reg$sv_end_chr){ 
    if(n_sp>1) {
      bam$bam_reg_svlen=list(seq_info); 
      for(bf in bam_files[[1]])
        bam$bam_cvg[[1]]=IRanges::rbind(bam$bam_cvg[[1]],DataFrame(allFunction(seq_info,bf,what=what)))
    } else {
      bam$bam_reg_svlen=list(GRanges(reg$sv_start_chr,IRanges(start=c(max(1,round(st0-0.5*buf_svlen)),ed0),end=c(st0,ed0+round(.5*buf_svlen)))))
      for(bf in bam_files[[1]])
        bam$bam_cvg[[1]]=IRanges::rbind(bam$bam_cvg[[1]],DataFrame(allFunction(bam$bam_reg_svlen[[1]],bf,what=what)))
    }
    bam$bam_reg_svlen[[2]]=seq_info; bam$bam_cvg[[2]]=DataFrame()
    for(bf in bam_files[[n_sp]])
      bam$bam_cvg[[2]]=IRanges::rbind(bam$bam_cvg[[2]],DataFrame(allFunction(seq_info,bf,what=what)))
  }
  #if(load_opt$verbose) cat("=Info: loaded bam_cvg\n")
  #bam_strad, bam_soft, bam_hang share the same loading region
  seq_info_st=GRanges(reg$sv_start_chr, IRanges(start=st0-buf_isize[[n_sp]],end=st0+buf_isize[[n_sp]]))
  seq_info_ed=GRanges(reg$sv_end_chr, IRanges(start=ed0-buf_isize[[n_sp]],end=ed0+buf_isize[[n_sp]]))
  bam$bam_reg_st=seq_info_st; bam$bam_reg_st_len=width(ranges(seq_info_st))
  bam$bam_reg_ed=seq_info_ed; bam$bam_reg_ed_len=width(ranges(seq_info_ed))
  what=c("pos","isize","flag"); bam$bam_strad_st=list(); bam$bam_strad_ed=list();
  for(i in seq_len(n_sp)) { bam$bam_strad_st[[i]]=DataFrame(); bam$bam_strad_ed[[i]]=DataFrame(); }
  #load bam_strad only if sv on one chr
  if(reg$sv_start_chr==reg$sv_end_chr){ 
    if(n_sp>1) {
      for(bi in seq_len(length(bam_files[[1]]))){   
        bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[1]][bi],what=what));lf=parseFlags(bam_tmp$flag)
        bam$bam_strad_st[[1]][[bi]]=bam_tmp[lf$isProperPair==1,]
        bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[1]][bi],what=what));lf=parseFlags(bam_tmp$flag)
        bam$bam_strad_ed[[1]][[bi]]=bam_tmp[lf$isProperPair==1,]
      }
    }
    for(bi in seq_len(length(bam_files[[n_sp]]))){
      bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[n_sp]][bi],what=what));lf=parseFlags(bam_tmp$flag)
      bam$bam_strad_st[[n_sp]][[bi]]=bam_tmp[lf$isProperPair==1,]
      bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[n_sp]][bi],what=what));lf=parseFlags(bam_tmp$flag)
      bam$bam_strad_ed[[n_sp]][[bi]]=bam_tmp[lf$isProperPair==1,]
    }
  }
  #if(load_opt$verbose) cat("loaded bam_strad\n")
  what=c("pos","strand","flag"); bam$bam_hang_st=list(); bam$bam_hang_ed=list()
  for(i in seq_len(n_sp)) { bam$bam_hang_st[[i]]=DataFrame(); bam$bam_hang_ed[[i]]=DataFrame(); }
  #lead bam_hang only if sv on one chr
  if(n_sp>1) {
    #bam$bam_hang_st[[1]]=DataFrame();bam$bam_hang_ed[[1]]=DataFrame();
    for(bi in seq_len(length(bam_files[[1]]))){   
      bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[1]][bi],what=what));lf=parseFlags(bam_tmp$flag)
      bam$bam_hang_st[[1]]=IRanges::rbind(bam$bam_hang_st[[1]],bam_tmp[lf$isMateUnmapped==1,])
      bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[1]][bi],what=what));lf=parseFlags(bam_tmp$flag)
      bam$bam_hang_ed[[1]]=IRanges::rbind(bam$bam_hang_ed[[1]],bam_tmp[lf$isMateUnmapped==1,])
    }
  }
  for(bi in seq_len(length(bam_files[[n_sp]]))){
    bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[n_sp]][bi],what=what));lf=parseFlags(bam_tmp$flag)
    bam$bam_hang_st[[n_sp]]=IRanges::rbind(bam$bam_hang_st[[n_sp]],bam_tmp[lf$isMateUnmapped==1,])
    bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[n_sp]][bi],what=what));lf=parseFlags(bam_tmp$flag)
    bam$bam_hang_ed[[n_sp]]=IRanges::rbind(bam$bam_hang_ed[[n_sp]],bam_tmp[lf$isMateUnmapped==1,])
  }
  #if(load_opt$verbose) cat("loaded bam_hang\n")
  #lead bam_soft only if sv on one chr
  what=c("pos","qwidth","cigar"); bam$bam_soft_st=list(); bam$bam_soft_ed=list()
  for(i in seq_len(n_sp)) { bam$bam_hang_st[[i]]=DataFrame(); bam$bam_hang_ed[[i]]=DataFrame(); }
  if(n_sp>1) {
    #bam$bam_soft_st[[1]]=DataFrame(); bam$bam_soft_ed[[1]]=DataFrame();
    for(bi in seq_len(length(bam_files[[1]]))){   
      bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[1]][bi],what=what))
      sc=grepl(complex_pattern,bam_tmp$cigar);bam$bam_soft_st[[1]]=IRanges::rbind(bam$bam_soft_st[[1]],bam_tmp[sc,])
      bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[1]][bi],what=what))
      sc=grepl(complex_pattern,bam_tmp$cigar);bam$bam_soft_ed[[1]]=IRanges::rbind(bam$bam_soft_ed[[1]],bam_tmp[sc,])
    }
  }
  #bam$bam_soft_st[[n_sp]]=DataFrame();bam$bam_soft_ed[[n_sp]]=DataFrame()
  for(bi in seq_len(length(bam_files[[n_sp]]))){
    bam_tmp=DataFrame(allFunction(seq_info_st,bam_files[[n_sp]][bi],what=what))
    sc=grepl(complex_pattern,bam_tmp$cigar);bam$bam_soft_st[[n_sp]]=IRanges::rbind(bam$bam_soft_st[[n_sp]],bam_tmp[sc,])
    bam_tmp=DataFrame(allFunction(seq_info_ed,bam_files[[n_sp]][bi],what=what))
    sc=grepl(complex_pattern,bam_tmp$cigar);bam$bam_soft_ed[[n_sp]]=IRanges::rbind(bam$bam_soft_ed[[n_sp]],bam_tmp[sc,])
  }
  #if(load_opt$verbose) cat("loaded bam_soft\n")
  bam
}

#FUNC: confirm by removing reg in hotspot 
conf_hot=function(ri,reg_lst,conf_opt){
  reg=reg_lst[[ri]];
  reg$sv_cvg=if(reg$sv_start_chr==reg$sv_end_chr & max(conf_opt[["rle_cvg"]][[reg$sv_start_chr]][reg$sv_start:reg$sv_end],na.rm=T)>=3*conf_opt[["mean_cvg"]][[reg$sv_start_chr]]) TRUE else FALSE
  reg
}

#FUNC: confirm by testing region coverage
#options: mincvg=1000, don't conf_cvg for SV<1000bp
conf_cvg=function(ri,reg_lst,bam_lst,conf_opt) {
  #list[ri,reg_lst,bam_lst,conf_opt]=list(1,reg_lst,bam_lst,conf_opt)
  reg=reg_lst[[ri]]; bam=bam_lst[[ri]]; cov=list(); max_del=conf_opt$max_del; min_dup=conf_opt$min_dup;
  r_lst=conf_opt$r_lst
  bonferroni=conf_opt$p_thresh/length(reg_lst)
  #mono cov1; matched cov1 is control, cov2 is case
  sel=which(!is.na(bam$bam_cvg[[1]]$pos))
  cov[[1]]=sum(IRanges::coverage(IRanges(start=bam$bam_cvg[[1]]$pos[sel],width=bam$bam_cvg[[1]]$qwidth[sel])))
  sel=which(!is.na(bam$bam_cvg[[2]]$pos))
  cov[[2]]=sum(IRanges::coverage(IRanges(start=bam$bam_cvg[[2]]$pos[sel],width=bam$bam_cvg[[2]]$qwidth[sel])))
  lk=lapply(r_lst,function(x){dpois(round(cov[[2]]/bam$bam_buf_svlen),(1-x)*cov[[1]]/bam$bam_buf_svlen)})
  rmax=r_lst[which.max(lk)]; pt=poisson.test(cov[[2]],bam$bam_buf_svlen,cov[[1]]/bam$bam_buf_svlen)
  if((cov[[2]]/cov[[1]])<=max_del & pt$p.value<bonferroni){
    reg$conf_type[1]="DEL"; reg$conf_mtd[1]=TRUE
    reg$conf_st[1]=NA; reg$conf_ed[1]=NA
    reg$conf_p[1]=pt$p.value; reg$conf_val[1]=cov[[2]]/bam$bam_buf_svlen
    reg$conf_freq[1]=rmax 
    #if(conf_opt$verbose) cat("=Info: conf_cvg confirmed DEL\n")
  } else if((cov[[2]]/cov[[1]])>=min_dup & pt$p.value<bonferroni){
    reg$conf_type[1]="DUP"; reg$conf_mtd[1]=TRUE
    reg$conf_st[1]=NA; reg$conf_ed[1]=NA
    reg$conf_p[1]=pt$p.value; reg$conf_val[1]=cov[[2]]/bam$bam_buf_svlen
    reg$conf_freq[1]=2-cov[[2]]/cov[[1]] #confounded k and r, r=(k-b)/(k-1), duplicate once, k=2
    #if(conf_opt$verbose) cat("=Info: conf_cvg confirmed DUP \n")
  }
  reg
}

#FUNC: convert 3d index to 1d index
major_order=function(idx,nd1,nd2,nd3){ #2,3,4  0->0,0,0; 11->1,0,0; 18->1,2,0
  offset=idx
  x = offset %% nd1
  offset = offset%/%nd1
  y = offset %% nd2
  offset = offset%/%nd2
  z = offset %% nd3
  offset = offset%/%nd3
  #cat(z,y,x,"\n")
  list(x,y,z)
}

#cat("loaded line 2000...\n")

#FUNC: return the mode of vector x
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#FUNC: confirm by testing hanging count
conf_hang=function(ri,reg_lst,bam_lst,swan_par,conf_opt) {
  #list[ri,reg_lst,bam_lst,conf_opt]=list(1,reg_lst,bam_lst,swan_par,conf_opt)
  reg=reg_lst[[ri]];bam=bam_lst[[ri]];n_sp=conf_opt$n_sp
  #print(n_sp)
  #print(names(bam))
  #print(swan_par)
  bonferroni=conf_opt$p_thresh/length(reg_lst);r_lst=conf_opt$r_lst
  p_plus=sum(swan_par[[n_sp]]$lambda*swan_par[[n_sp]]$p_left)/sum(swan_par[[n_sp]]$lambda)
  p_minus=sum(swan_par[[n_sp]]$lambda*swan_par[[n_sp]]$p_left)/sum(swan_par[[n_sp]]$lambda)
  #st: we need mate is - strand and hang reads, which means itself is + strand with isize>0
  #ed: we need mate is + strand and hang radds, which means itself is - strand with isize<0
  n_st_plus=sum(bam$bam_hang_st[[n_sp]]$strand=="+"); n_st_len=bam$bam_reg_st_len
  n_ed_minus=sum(bam$bam_hang_ed[[n_sp]]$strand=="-"); n_ed_len=bam$bam_reg_ed_len
  #cat(n_st_plus,n_st_len,p_plus,"\n")
  pt_st=poisson.test(n_st_plus,n_st_len,p_plus)
  pt_ed=poisson.test(n_ed_minus,n_ed_len,p_minus)
  lk_st=lapply(r_lst,function(x){dbinom(n_st_plus,n_st_len,(1-x)*sum(swan_par$lambda))})
  lk_ed=lapply(r_lst,function(x){dbinom(n_ed_minus,n_ed_len,(1-x)*sum(swan_par$lambda))})
  rmax_st=r_lst[which.max(lk_st)];rmax_ed=r_lst[which.max(lk_ed)];rmax=mean(rmax_st,rmax_ed)
  if(pt_st$p.value<bonferroni & pt_ed$p.value<bonferroni){ 
    reg$conf_type[2]=NA; reg$conf_mtd[2]=TRUE
    reg$conf_p[2]=max(pt_st$p.value,pt_ed$p.value)
    reg$conf_val[2]=max(n_st_plus/n_st_len,n_ed_minus/n_ed_len)
    reg$conf_freq[2]=rmax 
    reg$conf_st[2]=NA; reg$conf_ed[2]=NA
    #if(conf_opt$verbose) cat("Info: conf_hang confirmed DEL \n")
  }
  reg
}#this is a weak confirmation

#FUNC: return the value that has max frequency
tab_max=function(V){
  Tab=as.data.frame(table(V))
  Freq=Tab[,2]; M=which.max(Freq)
  c(as.numeric(as.character(Tab[,1][M])),Freq[M]) 
  #return what value is c(v_most_frequent, count)
}

#FUNC: confirm by testing softclipping count
conf_soft=function(ri,reg_lst,bam_lst,conf_opt) {
  #list[ri,reg_lst,bam_lst,conf_opt]=list(9,reg_lst,bam_lst,conf_opt)
  reg=reg_lst[[ri]]; bam=bam_lst[[ri]]; n_sp=conf_opt$n_sp; bonferroni=conf_opt$p_thresh/length(reg_lst)
  min_sc=conf_opt$min_sc; max_sc_ins=10
  #st: we need pos table of tail clipped reads
  #ed: we need pos table of head clipped reads
  st_tail_clip=sapply(str_match_all_perl(bam$bam_soft_st[[n_sp]]$cigar,tail_pattern_new), get_al_new)
  st_head_clip=sapply(str_match_all_perl(bam$bam_soft_st[[n_sp]]$cigar,head_pattern_new), get_al_new)
  #ed_head_clip=sapply(str_match_all_perl(bam$bam_soft_ed[[n_sp]]$cigar,head_pattern_new), get_al_new)
  #print(bam$bam_soft_st[[n_sp]]$cigar); print(st_head_clip); print(st_tail_clip)
  if(length(st_head_clip)>0&length(st_tail_clip)>0){
    st_sc_pos=bam$bam_soft_st[[n_sp]]$pos-st_head_clip+bam$bam_soft_st[[n_sp]]$qwidth-st_tail_clip #recover clip pos
    ed_sc_pos=bam$bam_soft_ed[[n_sp]]$pos
    st_sc=if(length(st_sc_pos)>0) tab_max(st_sc_pos) else c(reg$sv_start,0) 
    ed_sc=if(length(ed_sc_pos)>0) tab_max(ed_sc_pos) else c(reg$sv_end,0)
  } else {
    st_sc=c(reg$sv_start,0); 
    ed_sc=c(reg$sv_end,0)
  }
  #ed_sel=ed_head_clip>0; st_sel=st_tail_clip>0
  #cat(ri,"of",length(reg_lst),st_sc,ed_sc,"\n")
  if(st_sc[2]>min_sc&ed_sc[2]>min_sc){ #week confirmation by counting
    if(abs(st_sc[1]-ed_sc[1])<max_sc_ins){ #propose for insertion
      reg$conf_type[3]="INS"
      #if(conf_opt$verbose) cat("Info: conf_soft confirmed INS \n")
    } else {
      reg$conf_type[3]="DEL"
      #if(conf_opt$verbose) cat("Info: conf_soft confirmed DEL \n")      
    }
    reg$conf_mtd[3]=TRUE
    reg$conf_p[3]=NA
    reg$conf_val[3]=max(st_sc[2],ed_sc[2])
    reg$conf_freq[3]=NA
    reg$conf_st[3]=st_sc[1]
    reg$conf_ed[3]=ed_sc[1]
  }
  reg
}

#FUNC: confirm by testing stradelling size
conf_strad=function(ri,reg_lst,bam_lst,swan_par,conf_opt) {
  #list[ri,reg_lst,bam_lst,swan_par,conf_opt]=list(24,reg_lst,bam_lst,swan_par,conf_opt)
  #appearantly we need to load more than ci defined here
  reg=reg_lst[[ri]]; bam=bam_lst[[ri]];
  stopifnot(reg$sv_end>=reg$sv_start)
  stopifnot(!any(is.na(c(bam$sv_inner,bam$sv_outter)))) #these has to be proper
  max_run=conf_opt$max_run; n_sp=conf_opt$n_sp; r_lst=conf_opt$r_lst; min_mpr=conf_opt$min_mpr
  bonferroni=conf_opt$p_thresh/length(reg_lst); lrt_thresh=conf_opt$lrt_thresh; diff_thresh=conf_opt$diff_thresh
  ci_start_len=min(bam$sv_inner[1],bam$sv_inner[2])-bam$sv_outter[1]
  ci_end_len=bam$sv_outter[2]-max(bam$sv_inner[1],bam$sv_inner[2])
  step=max(1,round(sqrt(ci_start_len*ci_end_len*length(r_lst)/max_run)))
  st=seq(bam$sv_outter[1],min(bam$sv_inner[1],bam$sv_inner[2]),by=step)
  ed=seq(max(bam$sv_inner[1],bam$sv_inner[2]),bam$sv_outter[2],by=step)
  st_len=length(st); ed_len=length(ed); r_len=length(r_lst); 
  lrt_vec=rep(0,r_len*st_len*ed_len); #isize_sec_mode=c()
  #cat(reg$sv_start,reg$sv_end,r_len,st_len,ed_len,r_len*st_len*ed_len)  
  for(bi in seq_len(length(bam$bam_strad_st[[n_sp]]))){
    #bi=1, function(ti,r_len,st_len,ed_len,r_lst,st,ed,is,sd,rDF_st,rDF_ed)
    w_set=seq(min(ed)-max(st),max(ed)-min(st),by=step)+swan_par[[n_sp]]$rl[bi] #isize=x2-x1+rl
    is_set=unique(c(abs(bam$bam_strad_st[[n_sp]][[bi]]$isize),abs(bam$bam_strad_ed[[n_sp]][[bi]]$isize)))
    f_wi=lapply(w_set,function(w,is_set,is,sd){ fIw(is_set,w,is,sd) },is_set,swan_par[[n_sp]]$isize[bi],swan_par[[n_sp]]$isize_sdR[bi])
    lrt_vec=lrt_vec+sapply(seq(1,r_len*st_len*ed_len),strad_lrt,r_len,st_len,ed_len,r_lst,st,ed,bam$bam_strad_st[[n_sp]][[bi]],bam$bam_strad_ed[[n_sp]][[bi]],f_wi,is_set,swan_par[[n_sp]]$isize[bi]+3*swan_par[[n_sp]]$isize_sdR[bi])
    isize_v=c(abs(bam$bam_strad_st[[n_sp]][[bi]]$isize),abs(bam$bam_strad_ed[[n_sp]][[bi]]$isize))
    isize_p=c(abs(bam$bam_strad_st[[n_sp]][[bi]]$pos),abs(bam$bam_strad_ed[[n_sp]][[bi]]$pos))
    #isize_sec_mode=c(isize_sec_mode,Mode(abs(isize_v)[abs(isize_v)>swan_par[[n_sp]]$isize[bi]+swan_par[[n_sp]]$isize_sdR[bi]]))
  }
  
  ti=which.max(lrt_vec);list[rj,sj,ej]=major_order(ti-1,r_len,st_len,ed_len);rmax=r_lst[rj+1];smax=st[sj+1];emax=ed[ej+1]
  #print(lrt_vec); print(ti)
  #isize_diff=abs(mean(isize_sec_mode,na.rm=TRUE)-(emax-smax))
  #isize_diff=if(is.na(isize_diff)) diff_thresh else isize_diff  
  #print(emax-smax); print(isize_sec_mode)
  if(lrt_vec[ti]>lrt_thresh & sum(((isize_p<=(smax+step))&isize_v>(emax-smax-2*step))|((isize_p>=(emax-step))&abs(isize_v)>(emax-smax-2*step)))>=min_mpr){
    reg$conf_type[4]="DEL"; reg$conf_mtd[4]=TRUE 
    reg$conf_p[4]=NA; reg$conf_val[4]=lrt_vec[ti]
    reg$conf_st[4]=smax; reg$conf_ed[4]=emax; reg$conf_freq[4]=rmax
    #if(conf_opt$verbose) cat("Info: conf_strad confirmed DEL \n")
  }
  reg 
}

###############################


rDF_st_lrt = function(idx,rDF_st,rj,sj,wj,st,delta,is_set,r_lst,f_wi){
  lrt=0
  if(rDF_st$isize[idx]>0){
    if(rDF_st$pos[idx]<=st[sj] & rDF_st$pos[idx]>st[sj]-delta) { #only these as mixture others background
      isj=which(is_set==rDF_st$isize[idx])
      lrt=log((1-r_lst[rj])+r_lst[rj]*f_wi[[wj]][isj])
    }
  }
  return(lrt)
}

rDF_ed_lrt = function(idx,rDF_ed,rj,ej,wj,ed,delta,is_set,r_lst,f_wi){
  lrt=0
  if(rDF_ed$isize[idx]<0){
    if(rDF_ed$pos[idx]>=ed[ej] & rDF_ed$pos[idx]>=ed[ej]+delta) { #only these as mixture others background
      isj=which(is_set==abs(rDF_ed$isize[idx])) 
      lrt=log((1-r_lst[rj])+r_lst[rj]*f_wi[[wj]][isj])
    }
  }
  return(lrt)
}

strad_lrt=function(ti,r_len,st_len,ed_len,r_lst,st,ed,rDF_st,rDF_ed,f_wi,is_set,delta) {
  #list[ti,r_len,st_len,ed_len,r_lst,st,ed,is,sd,rDF_st,rDF_ed,f_wi,is_set]=list(1,r_len,st_len,ed_len,r_lst,st,ed,swan_par[[n_sp]]$isize[1],swan_par[[n_sp]]$isize_sdR[1],bam$bam_strad_st[[n_sp]][[1]],bam$bam_strad_ed[[n_sp]][[1]],f_wi,is_set) 
  list[rj,sj,ej]=major_order(ti-1,r_len,st_len,ed_len); rj=rj+1; sj=sj+1; ej=ej+1; lrt=0; wj=ej-sj+st_len
  lrt=sum(sapply(seq_len(nrow(rDF_st)),rDF_st_lrt,rDF_st,rj,sj,wj,st,delta,is_set,r_lst,f_wi))
  lrt=lrt+sum(sapply(seq_len(nrow(rDF_ed)),rDF_ed_lrt,rDF_ed,rj,ej,wj,ed,delta,is_set,r_lst,f_wi))
  #for(idx in seq_len(nrow(rDF_st))){ #idx=1 }
  #for(idx in seq_len(nrow(rDF_ed))){}
  #cat("bam_strad, done",ti,"of",r_len*st_len*ed_len,"\n")
  lrt
}

MIN.CVG.CHANGE.LEN=10000
MIN.ISIZE.CHANGE.LEN=500

conf_join=function(ri,reg_lst,bam_lst,swan_par,conf_opt,hot_only=TRUE) {
  reg=reg_lst[[ri]]; verbose=conf_opt$verbose; idx=0
  #cat("=Info: reg",ri,"of",length(reg_lst),", size",reg$sv_end-reg$sv_start)
  #rule of thumb, let's do easy ones first; no confirmation at all, throw away
  if(!hot_only){
    bam=bam_lst[[1]]
    if( (reg$sv_end-reg$sv_start) < MIN.ISIZE.CHANGE.LEN) {
      reg=conf_hang(1,list(reg),list(bam),swan_par,conf_opt)
      reg=conf_soft(1,list(reg),list(bam),conf_opt)
      if(reg$conf_mtd[2] & reg$conf_mtd[3]){
        reg$sv_conf=TRUE; idx=3
      }
    }
    if( (reg$sv_end-reg$sv_start) >= MIN.ISIZE.CHANGE.LEN & (reg$sv_end-reg$sv_start)<MIN.CVG.CHANGE.LEN ) {
      reg=conf_strad(1,list(reg),list(bam),swan_par,conf_opt)
      reg=conf_cvg(1,list(reg),list(bam),conf_opt)
      if(reg$conf_mtd[4] & reg$conf_mtd[1]){
        reg$sv_conf=TRUE; idx=4
      }
    }
    if( (reg$sv_end-reg$sv_start) >= MIN.CVG.CHANGE.LEN) {
      reg=conf_cvg(1,list(reg),list(bam),conf_opt)
      if(reg$conf_mtd[1]) {
        reg$sv_conf=TRUE; idx=1
      }
    }
  }
  output=(!hot_only&reg$sv_conf&!reg$sv_cvg)|(hot_only&!reg$sv_cvg) #reg that will be kept
  if(output&idx!=0){ #only carry this out if doing serious confirmation
    reg$sv_start=if(!is.na(reg$conf_st[idx])) reg$conf_st[idx] else reg$sv_start
    reg$sv_end=if(!is.na(reg$conf_ed[idx])) reg$conf_ed[idx] else reg$sv_end
    reg$sv_type=if(!is.na(reg$conf_type[idx])) reg$conf_type[idx] else reg$sv_type
    reg$sv_freq=if(!is.na(reg$conf_freq[idx])) reg$conf_freq[idx] else reg$sv_freq
    #cat("mtd",reg$sv_mtd,"conf",reg$conf_mtd,"->",reg$sv_conf,"hot",reg$sv_cvg,taggie(conf_opt$gtk),"\n")
  } else{
    #no conf change, keep original information
  }
  if(output) reg else NULL
} #basically decide if conf=T/F and sort out sv_type, add stat and pvalue in vcf_merge 

vcf_merge=function(final_set,sample_id,sample_info,contig_info,ref_seq,ref_file,verbose=TRUE){
  vcf_field=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  bed_field=c("CHROM","POS","END","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") #bed format containing the END field
  vcf_class=c("character","integer","character","character","character","character","character","character","character","character")
  bed_class=c("character","integer","integer","character","character","character","character","character","character","character","character")
  n_sv = length(final_set); #if(verbose) cat("==Info: number of final calls=", n_sv,"\n")
  #cat(length(rep(rep("",length(vcf_field)+1),n_sv)),length(vcf_field)+1,"\n")
  #cat(length(rep(rep("",length(bed_field)+1),n_sv)),length(bed_field)+1,"\n")
  #vcf_empty = base::matrix(rep(rep("",length(vcf_field)+1),n_sv), ncol=length(vcf_field)+1)
  #bed_empty = base::matrix(rep(rep("",length(bed_field)+1),n_sv), ncol=length(bed_field)+1)
  #vcf_set = data.frame(vcf_empty, stringsAsFactors=F)
  #bed_set = data.frame(bed_empty, stringsAsFactors=F)
  vcf_set = data.frame(); bed_set=data.frame()
  
  for(vi in seq_len(n_sv)){
    #if break points is not found or sketchy
    #print(final_set[[vi]])
    sv_imprecise=if(final_set[[vi]]$sv_imprecise) "IMPRECISE:" else ""
    #print("-------------------------------------")
    #print(final_set[[vi]])
    sv_filter=if(!is.na(final_set[[vi]]$sv_filter)) "." else final_set[[vi]]$sv_filter 
    sv_qual=if(is.na(final_set[[vi]]$sv_qual)) "." else final_set[[vi]]$sv_qual
    sv_type=paste(final_set[[vi]]$sv_type)
    sv_len=if(sv_type=="INS") 0 else if (sv_type %in% c("DEL","TANDEM.DUP","INV","INS.DUP")) final_set[[vi]]$sv_end-final_set[[vi]]$sv_start else NA
    sv_len=if(is.na(sv_len)) NA else if(sv_len<0) 0 else round(sv_len) 
    tryCatch( {
        ref_base=toString(subseq(ref_seq[[final_set[[vi]]$sv_chrname]],final_set[[vi]]$sv_start,final_set[[vi]]$sv_start))
      }, error=function(e) { 
        cat("===Warning: can't read ", final_set[[vi]]$sv_chrname, ":", final_set[[vi]]$sv_start, "-", final_set[[vi]]$sv_start, "\n")
        ref_base="N"
      })
    sv_format="MTD:VAL:STAT1:STAT2:CVG"
    sv_dummy=paste(final_set[[vi]]$sv_mtd,round(final_set[[vi]]$sv_val,2),round(final_set[[vi]]$sv_stat1,2),round(final_set[[vi]]$sv_stat2,2),round(final_set[[vi]]$sv_cvg,2),sep=":")
    #sv_format="TSC_FSC_SSC_PSC:TST_FST_SST_PST:THA_FHA_SHA_PHA:TCV_FCV_SCV_PCV";
    #sv_dummy=paste(final_set[[vi]]$conf_type,round(final_set[[vi]]$conf_freq,digits=2),
    #               round(final_set[[vi]]$conf_val,digits=2),
    #               round(final_set[[vi]]$conf_p,digits=2),sep="_",collapse=":")
    sv_info=gettextf("%sSVTYPE=<%s>;END=%s;ENDCHR=%s;SVLEN=%s;CIPOS=%s,%s;CIEND=%s,%s;AA=%s;DP=%s;",sv_imprecise,sv_type,final_set[[vi]]$sv_end,final_set[[vi]]$sv_end_chr,sv_len,final_set[[vi]]$sv_start_ci_low,final_set[[vi]]$sv_start_ci_high,final_set[[vi]]$sv_end_ci_low,final_set[[vi]]$sv_end_ci_high,paste(final_set[[vi]]$sv_freq),paste(final_set[[vi]]$sv_dep))
    #if(length(sv_info)==0) cat(sv_imprecise,sv_type,final_set[[vi]]$sv_end,sv_len,final_set[[vi]]$sv_start_ci_low,final_set[[vi]]$sv_start_ci_high,final_set[[vi]]$sv_end_ci_low,final_set[[vi]]$sv_end_ci_high,final_set[[vi]]$sv_freq,final_set[[vi]]$sv_dep,"\n")
    #print(sv_info)
    #we might need better ideas here...
    vcf_tmp=as.data.frame(list(final_set[[vi]]$sv_chrname,final_set[[vi]]$sv_start,final_set[[vi]]$sv_id,ref_base,sv_type,sv_qual,sv_filter,sv_info,sv_format,sv_dummy),stringsAsFactors=FALSE)
    colnames(vcf_tmp) = c(vcf_field,sample_id); 
    #print(vcf_tmp)
    vcf_set=rbind(vcf_set,vcf_tmp)
    bed_tmp=as.data.frame(list(final_set[[vi]]$sv_chrname,final_set[[vi]]$sv_start,final_set[[vi]]$sv_end,final_set[[vi]]$sv_id,ref_base,sv_type,sv_qual,sv_filter,sv_info,sv_format,sv_dummy),stringsAsFactors=FALSE)
    #print(bed_tmp)
    colnames(bed_tmp) = c(bed_field,sample_id); bed_set=rbind(bed_set,bed_tmp)
  }
  vcf_meta=paste(call_meta(
    c(meta_keys, sv_keys), c(meta_values, sv_values), ref_file, "swan_join.R"),sample_meta(sample_info[,1], sample_info[,2], sample_info[,3], sample_info[,4]),contig_meta(contig_info[,1], contig_info[,2], contig_info[,3], contig_info[,4]),gettextf("#%s", paste(colnames(vcf_set), collapse='\t')), sep='\n')
  vcf_set=set_colClass(vcf_set,vcf_class)
  bed_set=set_colClass(bed_set,bed_class)
  return(list(vcf_set=vcf_set,vcf_meta=vcf_meta,bed_set=bed_set))
}

#then we call tracks based on this rle
#we need to read in the tracks
#we then call peaks based on thresh rle
#we then merge paks within tele distance

merge_rclu=function(lstart,lend,rstart,rend,support){
  #lstart, lend, rstart, rend, support
  #print(lstart);print(lend);print(rstart);print(rend);print(support)
  list[lCluster,lMap]=tryCatch({
    lCluster=IRanges::reduce(IRanges(start=lstart,end=lend),with.mapping=T)
    lMap=mcols(lCluster)$mapping
    list(lCluster=lCluster,lMap=lMap)
  }, error=function(e){
    lCluster=IRanges::reduce(IRanges(start=lstart,end=lend),with.revmap=T)
    lMap=mcols(lCluster)$revmap
    list(lCluster=lCluster,lMap=lMap)
  })
  #cat("lCluster\n"); print(lCluster); cat("mapping\n"); print(lMap)
  read_lC = hash::hash()
  for(i in seq_len(length(lMap)))
    for(j in lMap[[i]])
      read_lC[j]=i
  list[rCluster,rMap]=tryCatch({
    rCluster=IRanges::reduce(IRanges(start=rstart,end=rend),with.mapping=T)
    rMap=mcols(rCluster)$mapping
    list(rCluster=rCluster,rMap=rMap)
  }, error=function(e){
    rCluster=IRanges::reduce(IRanges(start=rstart,end=rend),with.revmap=T)
    rMap=mcols(rCluster)$revmap
    list(rCluster=rCluster,rMap=rMap)
  })
  read_rC = hash::hash()
  for(i in seq_len(length(rMap)))
    for(j in rMap[[i]])
      read_rC[j]=i
  #cat("rCluster\n"); print(rCluster); cat("mapping\n"); print(rMap)
  pm=new(PairMap)
  for(i in seq_len(length(lstart))){
    pm$add(read_lC[[as.character(i)]],read_rC[[as.character(i)]],support[i])
  }
  lCb = pm$read()
  as.data.frame(list(lstart=start(lCluster)[as.vector(lCb$x)],
                     lend=end(lCluster)[as.vector(lCb$x)],
                     rstart=start(rCluster)[as.vector(lCb$y)],
                     rend=end(rCluster)[as.vector(lCb$y)],
                     support=as.vector(lCb$n)))
}

set_isize=function(isize_text,bam_file,seq_info,max_isize=NA){
  if(verbose) cat("==Info: setting isize by",isize_text,"\n")
  if(grepl("^[0-9]",isize_text)){ 
    list[isize,isize_sdR,isize_sdL]=as.numeric(strsplit(isize_text,split="_")[[1]])
  } else if(isize_text=="learn"){
    #print(seq_info)
    #print(bam_file)
    isize_V=allFunction(seq_info,bam_file,what=c("isize"),index=bam_file)$isize
    if(length(isize_V)==0) {
      isize=-1; isize_sdR=-1; isize_sdL=-1; good_isize=NA;
      cat(paste("-Error: no reads loaded from",bam_file,", check .bam and .bai\n"))
    } else {
      pars=isizePar(isize_V[isize_V>0&!is.na(isize_V)],MAX.ISIZE=NA,doplot=F,method=2)
      isize=pars$target.isize; isize_sdR=round(pars$sdR); isize_sdL=round(pars$sdL); 
      good_isize=pars$good.for.insertion
    }
  } else {
    isize=-1; isize_sdR=-1; isize_sdL=-1; good_isize=NA
  } 
  return(list(isize,isize_sdR,isize_sdL,good_isize))  
}

set_coverage=function(coverage_text){ 
  if(grepl("^[0-9]",coverage_text)) return(as.numeric(coverage_text)) else return(-1)
}
set_RL=function(RL_text){ if(grepl("^[0-9]",RL_text)) return(as.numeric(RL_text)) else return(-1) }
set_Delta=function(Delta_text){ 
  if(grepl("^[0-9]",Delta_text)) return(as.numeric(Delta_text)) else return(-1) }
set_bigDel=function(bigDel_text){
  if(grepl("^[0-9]",bigDel_text)) return(as.numeric(bigDel_text)) else return(-1) }
set_smallDel=function(smallDel_text){
  if(grepl("^[0-9]",smallDel_text)) return(as.numeric(smallDel_text)) else return(-1) }
set_smallIns=function(smallIns_text){
  if(grepl("^[0-9]",smallIns_text)) return(as.numeric(smallIns_text)) else return(-1) }
set_maxInsert=function(maxInsert_text){
  if(grepl("^[0-9]",maxInsert_text)) return(as.numeric(maxInsert_text)) else return(-1) }
set_p=function(p_text){ if(grepl("^[0-9]",p_text)) return(as.numeric(p_text)) else return(-1) }
set_q=function(q_text){ if(grepl("^[0-9]",q_text)) return(as.numeric(q_text)) else return(-1) }
set_hang_clip=function(hang_clip_text) { 
  if(grepl("^[0-9]",hang_clip_text)) return(as.numeric(hang_clip_text)) else return(-1) }
set_prop_clip=function(prop_clip_text) {
  if(grepl("^[0-9]",prop_clip_text)) return(as.numeric(prop_clip_text)) else return(-1) }

### memory utilities
showMemoryUse='
sort="size"
limit=10
decreasing=TRUE
objectList <- ls()
oneKB <- 1024
oneMB <- 1048576
oneGB <- 1073741824
memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
memListing <- sapply(memoryUse, function(size) {
if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
else return(paste(size, "bytes"))
})
memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
if (sort=="alphabetical") { memListing = memListing[order(memListing$objectName,decreasing=decreasing),] 
} else { memListing = memListing[order(memoryUse,decreasing=decreasing),] }
if(!missing(limit)) memListing <- memListing[1:limit,]
print(memListing, row.names=FALSE)
'
get_rule_code=function(prules,nrules,namespace){ 
  pcode_set=c(); ncode_set=c(); 
  np_set=lapply(namespace,sets::as.set)
  np_power=2^sets::as.set(np_set)
  for(np in np_power){
    for(r in prules){
      prule=sets::as.set(lapply(r,sets::as.set))
      if(prule<=np) pcode_set=c(pcode_set,rule2code(unlist(np),namespace))
    }
  }
  pcode_set=unique(pcode_set)
  for(np in np_power){
    for(r in nrules){
      nrule=sets::as.set(lapply(r,sets::as.set))
      if(nrule<=np) ncode_set=c(ncode_set,rule2code(unlist(np),namespace))
    }
  }
  ncode_set=unique(ncode_set)
  pcode_set[!(pcode_set %in% ncode_set)]
}

load_RData=function(ixs,dprefix,trunked=TRUE){
  if(!trunked) {
    all_RData=new.env()
    load(paste(dprefix,".RData",sep=""),envir=all_RData)
    all_RData=as.list(all_RData)
  } else {
    tmp_RData=new.env() #must be indexed by character
    for(i in seq_along(ixs)){
      ci=as.character(i)
      tmp_RData[[ci]]=new.env()
      cat("loading",paste(dprefix,".trunk",ixs[i],".RData",sep=""),"\n")
      load(paste(dprefix,".trunk",ixs[i],".RData",sep=""),envir=tmp_RData[[ci]])
    }
    var_names=ls(tmp_RData[["1"]])
    var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },tmp_RData[["1"]])
    all_RData=lapply(seq_along(var_names),merge_RData,var_names,var_classes,ixs,tmp_RData)
    for(di in seq_along(all_RData)) { all_RData[[di]]=if(length(all_RData[[di]])==0) list() else all_RData[[di]][[1]] }
    names(all_RData)=var_names
  }
  return(all_RData)
}

merge_RData=function(j,var_names,var_classes,ixs,data){ #merege character indexed list, data[[ci]][[ci]]
  var_name=var_names[j]; var_class=var_classes[[j]]; var_val=NULL; var_list=list()
  #cat("var_name=",var_name,"var_class=",var_class,"\n")
  for(i in seq_along(ixs)){
    ci=as.character(i)
    switch(var_class,
           "integer"={var_val=if(is.null(var_val)) data[[ci]][[var_name]] else c(var_val,data[[ci]][[var_name]])},"numeric"={var_val=if(is.null(var_val)) data[[ci]][[var_name]] else c(var_val,data[[ci]][[var_name]])},"data.frame"={var_val=if(is.null(var_val)) data[[ci]][[var_name]] else rbind(var_val,data[[ci]][[var_name]])},"IRanges"={var_val=if(is.null(var_val)) data[[ci]][[var_name]] else c(var_val,data[[ci]][[var_name]])},"Rle"={var_val=if(is.null(var_val)) data[[ci]][[var_name]] else c(var_val,data[[ci]][[var_name]])})
  }
  var_list[[var_name]]=var_val
  return(var_list)
} #this may possibly replacable by new_merge_RData but let's keep this for now

new_merge_RData=function(j,var_names,var_classes,data){ #merge number indexed list
  var_name=var_names[j]; var_class=var_classes[[j]]; var_val=NULL; var_list=list()
  #cat("var_name=",var_name,"var_class=",var_class,"\n")
  if (var_class == "data.frame") {
    var_val = do.call(rbind, lapply(data, function (o) { o[[var_name]] }))
  } else {
    var_val = do.call(c, lapply(data, function (o) { o[[var_name]] }))
  }
  var_list[[var_name]]=var_val; return(var_list)
}

#jumpy.redroo <- function(x) .Internal(sum(2^.subset((length(x)-1):0, x)))
jumpy.roo=function(x) { sum(2^.subset((length(x)-1):0, x)) } #most left is highest bit

rule2code=function(rule,namespace){ jumpy.roo(namespace %in% rule) }

code2call=function(track,codes){ track %in% codes }

modeDecomp=function(V,CI=100,RMIN=0.001,RDIP=0.25){
  # Charlie's function go give modes and their densities.
  Freq=as.data.frame(table(V))[,2]; N=length(V)
  Val=as.numeric(as.data.frame(table(V),stringsAsFactors=FALSE)[,1])
  M=c(); R=c(); RCUR=1; MASK=rep(FALSE,length(Freq));
  #cat("\nFreq=",Freq)
  #cat("\nVal=",Val)
  while(RCUR>RMIN & !all(MASK)){
    m=which.max(ifelse(MASK,0,Freq));
    RCUR=Freq[m]/N
    MASK=ifelse(!MASK,Val>Val[m]-CI & Val<Val[m]+CI,MASK); #turn off nearby region
    l=m; while(l>1&Freq[l]>=RCUR*RDIP*N) l=l-1; #turn off extending humps
    r=m; while(r<length(Val)&Freq[r]>=RCUR*RDIP*N) r=r+1; #turn off extending humps
    MASK[l:r]=TRUE;
    if(RCUR>=RMIN) { #cat("\nfound new mode at",Val[m],"with",RCUR)
      M=c(M,Val[m]); R=c(R,RCUR); }
  }
  ord=order(M)
  M=M[ord]; R=R[ord]
  return(list(M=M,R=R))
}


findModes=function(V,MAX.V=NA,CI=20,bw="nrd0",adjust=2){
  # same input and output modeDecomp, but uses a different algorithm.
  V = V[which(V<10*median(V))] # First do a very rough thresholding, to avoid bogus sd estimates in the next step.
  if(is.na(MAX.V)) MAX.V = median(V)+6*sd(V)  # Only consider data within 6 sd from median.
  V = V[which(V<MAX.V & V>0)]
  d=density(V, bw=bw, adjust=adjust)
  maxinwin = rep(0,length(d$x))
  for(i in seq_along(d$x)){
    ix = which(abs(d$x-d$x[i])<CI & d$x!=d$x[i])
    maxinwin[i] = max(d$y[ix])
  }
  modes = which(d$y > maxinwin)
  M = round(d$x[modes])
  R = d$y[modes]
  ord=order(M)
  M=M[ord]; R=R[ord]
  return(list(M=M,R=R))
}

#cat("loaded line 2500...\n")

isizePar=function(isize,MAX.ISIZE=NA,doplot=TRUE,method=2,opre="isize_diagnosis"){
  ## Estimates library parameters used by SWAN.
  if(method==1){
    res=modeDecomp(isize)
    pdf_file=paste(opre,"char","hist","pdf",sep=".")
  } else {
    res=findModes(isize)
    pdf_file=paste(opre,"hist","pdf",sep=".")  
  }
  
  maxr = max(res$R) # this mode has the maximum mass.
  counter = length(res$R)
  while(res$R[counter] < maxr*0.1){
    counter = counter-1
  }
  # targetmode: assume that this was the "targeted" library size
  # Set it to be the right-most mode that is >= 10% of the size of the mode with largest mass.
  # Then, when there is contamination at the lower end, the contamination can not exceed 90%.
  # Also, when there  is contamination at the upper end, the proportion can not exceed 10%.
  
	cat("=Info: find rg isize modes at",res$M,"with densities",res$R,"\n")
  
  targetmode = res$M[counter]
  targetProportion= res$R[counter]
  good.for.deletion=TRUE
  if(targetProportion < maxr*0.5) good.for.deletion=FALSE
  
  # MAX.ISIZE is the maximum library size plausible.
  # It's best if a value is given based on what was chosen during the library prep stage.
  # If no value given, we will get a crude estimate.
  if(is.na(MAX.ISIZE)) MAX.ISIZE = targetmode+10*sd(isize[isize<targetmode*10])
  
  # Get the parameters for deletion detection.
  z = isize[which(isize>targetmode & isize<MAX.ISIZE)]
  sdhat = sqrt(mean((z-targetmode)^2))
  
  # Get the parameters for insertion detection.
  leftz = isize[which(isize<targetmode)]
  leftsdhat = sqrt(mean((leftz-targetmode)^2))
  # The library is not good for insertion if:
  #   (1) There are other modes to the left of targetmode, whose value is > 10% of the targetmode.
  #   (2) In the single mode case the proportion of data less than max(100,targetmode-2*leftsdhat) is greater than 0.1.
  no.dominating.left.modes=counter==1 || (counter>1 && max(res$R[seq_len(counter-1)])<0.1*res$R[counter])
  not.too.skewed=sum(leftz<max(100,targetmode-2*leftsdhat))/(2*length(leftz)) < 0.05
  good.for.insertion = no.dominating.left.modes & not.too.skewed
  
  if(doplot){
    temp=isize[which(isize<targetmode+6*sdhat)]
    pdf(pdf_file)
    h=hist(temp, 100,xlim=c(0,targetmode+6*sdhat),main=paste("Distribution of insert sizes (Insertions: ",good.for.insertion,", Deletions: ",good.for.deletion,")",sep=""), xlab="insert size", col="gray", border="gray")
    for(i in seq_along(res$M)){
      segments(res$M[i],0,res$M[i], max(h$counts)*2, col="red", lwd=1, lty=2)
    }
    segments(targetmode,0,targetmode, max(h$counts)*2, col="red", lwd=3, lty=1)
    
    segments(median(isize), 0, median(isize), max(h$counts)*2, col="green", lty=2)
    stepsize=1; x=seq(targetmode, max(isize), stepsize)
    lines(x, dnorm(x,mean=targetmode, sd=sdhat)*(2*length(z))*(h$breaks[2]-h$breaks[1]),col="cornflowerblue", lwd=3)
    stepsize=1; x=seq(0, targetmode, stepsize)
    lines(x, dnorm(x,mean=targetmode, sd=leftsdhat)*(2*length(leftz))*(h$breaks[2]-h$breaks[1]),col="orange", lwd=3)
    
    legend(x="topright", col=c("red","red","green", "cornflowerblue", "orange"), lty=c(1,2,2,1,1),lwd=c(3,1,1,2,2), legend=c("target insert size","local mode","median", "fitted normal (right)", "fitted normal (left)"))
    suppressMessages(dev.off())  
  }
  
  list(target.isize=targetmode, sdR = sdhat, sdL = leftsdhat, modes = res$M, good.for.insertion=good.for.insertion,good.for.deletion=good.for.deletion)
}

#change dataframe colClass
# Coerces data.frame columns to the specified classes
set_colClass <- function(d, colClasses) {
    colClasses <- rep(colClasses, len=length(d))
    d[] <- lapply(seq_along(d), function(i) switch(colClasses[i], 
        numeric=as.numeric(d[[i]]), 
        character=as.character(d[[i]]), 
        Date=as.Date(d[[i]], origin='1970-01-01'), 
        POSIXct=as.POSIXct(d[[i]], origin='1970-01-01'), 
        factor=as.factor(d[[i]]),
        as(d[[i]], colClasses[i]) ))
    return(d)
}

# call various method to calculate threshold

sig_thresh=function(data,scores,thresh_level,scan_par,method,adaptive=T,verbose=F){
  #assuming centromere and telomere region is minority and if existing is masked as NA
  data_mean=list(); data_sd=list(); thresh = list()
  if(verbose) cat("====Info: reading the call table using method=", method, "adaptive=", adaptive,
   "threshold=", thresh_level, "\n")
  #if(verbose) cat("*********\n")
  if(method == "boot"){
    for(i in seq_len(length(scores))) {
      s=scores[i]; alpha=alpha_level[thresh_level]
      #boot=(s %in% c("lD","lDl","lDr")) & adaptive # adaptive+lD => boot
      list[thresh[[s]],data_mean[[s]],data_sd[[s]]] =
        noovlap_est(data=data[[s]],resample=1000,
          buffer=round(scan_par$delta*2/scan_par$stepsize),alpha=alpha)
    }
  } else if(method == "robust") {
    for(i in seq_len(length(scores))) {
      s=scores[i]
      list[thresh[[s]],data_mean[[s]],data_sd[[s]]] = 
        robust_est(data=data[[s]],lvl=thresh_level)
    }
  } else if(method == "theo") {
    for(i in seq_len(length(scores))) {
      s=scores[i]; alpha=alpha_level[thresh_level]
      list[thresh[[s]],data_mean[[s]],data_sd[[s]]] =
        theo_thresh(score=s,alpha=alpha,scan_par=scan_par,verbose=verbose)
    }
  } else if(method == "mad"){
    stop("method not yet implemented, bail out...\n")
  } else if(method == "fisher") {
    stop("method not implemented, bail out...\n")
  } else if(method == "stouffer") {
    stop("method not implemented, bail out...\n")
  } else {
    stop("method not implemented, bail out...\n")
  }
  #if(verbose) cat("*********\n")
  return(list(thresh=thresh,data_mean=data_mean,data_sd=data_sd))
}

### moving average function
pair<-function(x,y){ 0.5*(x+y)*(x+y+1)+x }

unpair<-function(z){
  w= floor( (sqrt(8*z+1) - 1)/2 )
  t = w*(w+1)/2
  cbind(z-t,w-z+t)
}

ma = function(x,n=mvsum_fail){ filter(x,rep(1,n), sides=2) }

robust_est= function(data,lvl){
  data=data[!is.na(data)]
  med_val = median(data, na.rm=T)
  par_val = Sn(data)
  if(par_val==0) { cat("====Warn: estimated sd is 0!\n") }
  lvl_thresh = med_val+par_val*lvl
  return(list(thresh=lvl_thresh,data_mean=med_val,data_sd=par_val))
}

noovlap_est= function(data,resample,buffer,alpha,max_try=10, verbose=FALSE){
  i_try=1
  resampleB=floor(length(data)/buffer) #resampleB is how many blocks we can divide the region w/o overlap
  #print(resampleB)
  R=floor((resample-1)/resampleB)+1 #R is how many repeats we need to reach at least resample samples, R*resampleB>=resample
  if(verbose) cat("====Warn: need to sample ", R, "times \n")
  if(resampleB<=resample) {
    Rresample=resampleB # we need to sample R times 
  } else {
    Rresample=resample  # we just sample once
  }
  alpha_data=rep(0,length(alpha))
  while(any(alpha_data==0) & i_try<=max_try){
    sample_idx=c()
    for(i in seq_len(R)){
      blocks=sample(resampleB,Rresample); shifts=sample(buffer,Rresample,replace=T)
      sample_idx=c(sample_idx,(blocks-1)*buffer+shifts)
    }
    sample_data = sort(data[sample_idx][!is.na(data[sample_idx])])
    sample_mean = median(sample_data)
    sample_sd = mad(sample_data)
    #cat("sample_data="); print(sample_data); print(length(sample_data))
    #cat("alpha="); print(alpha)
    #cat("alpha="); print(R*Rresample*(1-alpha)) 
    alpha_data = sample_data[round(length(sample_data)*(1-alpha))]
    i_try = i_try+1
    #cat("R=",R,"Rresample=",Rresample,"alpha=",round(R*Rresample*(1-alpha)),"sample=",length(sample_data),"data=",alpha_data,"\n")
  }
  if(any(alpha_data==0)) { cat("====Warn: some estimated threshold is 0! \n") }
  return(list(thresh=alpha_data,data_mean=sample_mean,data_sd=sample_sd))
}

bsearch=function(val,vec){
  L=1; U=length(vec); M=1;
  while(!(vec[M+1]>=val & vec[M]<=val)) {
    M=floor((L+U)/2)
    if(vec[M]>val) { U=M 
    } else if(vec[M]<val) { L=M }
  }
  if(vec[M+1]==val) M=M+1
  return(M)
}

#https://stat.ethz.ch/pipermail/r-help/2011-April/274182.html
bsearch1 <-
  function(val, tab, L=1L, H=length(tab))
  {
    while (H >= L) {
      M <- L + (H - L) %/% 2L
      if (tab[M] > val) H <- M - 1L
      else if (tab[M] < val) L <- M + 1L
      else return(M)
    }
    return(L - 1L)
  }

score_std=function(data,score,scaled,centered,halfed,resample,buffer,lvl,boot=F){ 
  if(boot) { #use resampling
    #cat("score =", score, "std by boot", "resample=", resample, "buffer=", buffer, "\n")
    #cat("lvl=",lvl)
    #cat("alpha_level=",alpha_level[lvl])
    list[data_mean,data_par] = noovlap_est(data,resample,buffer,alpha_level[lvl])
    data_par=(data_par-data_mean)/lvl #lvl used in estimating emperical sd, thin it first and will convert back
  } else {
    #cat("score =", score, "std by Sn\n") #use robust sd
    list[data_mean,data_par] = robust_est(data) #no lvl considered in estimating sd
  }
  if(verbose) cat("====Info: score=",score,"data_mean=",data_mean,"data_par=",data_par,"\n")
  data[is.na(data)] = data_mean #fill NA with data_mean
  if(scaled=='y' & centered=='y') {
    data=(data-data_mean)/data_par
    data_par=1; data_mean=0
    if(halfed=='y') data=pmax(0,data)
  }else if(scaled=='y' & centered=='n') {
    data=data/data_par
    data_sd=1
  }else if(scaled=='n' & centered=='y') {
    data=data-data_mean
    data_mean=0
    if(halfed=='y') data=pmax(0,data)
  }
  return(list(data=data,data_mean=data_mean,data_par=data_par))
} #use data_mean+lvl*data_par to get threshold at the lvl

### data operating functions ###
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}  # data structure to support list[a,b,c,...] = func(x) { list(a=a, b=b, c=c, ...) }

### fix numeric issues in ggplot2 ###
double.min.exp=.Machine$double.min.exp
double.max.exp=.Machine$double.max.exp
double.xmin=-.Machine$integer.max #this is used to replace -inf in plotting
double.xmax=.Machine$integer.max  #this is used to replace inf in plotting

meta_keys=c("fileformat","fileDate","source","reference","phasing")
meta_values=c("None","None","None","None","partial")
info_values=c(
'<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
'<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
'<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
'<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
'<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
'<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
'<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
'<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">',
'<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
'<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
'<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
'<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
'<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
'<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
'<ID=SOMATIC,Number=0,Type=Flag,Description="SOMATIC structural variation">',
'<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">',
'<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
'<ID=METRANS,Number=.,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">',
'<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
'<ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">',
'<ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVAR">',
'<ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">',
'<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">',
'<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">',
'<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
'<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the length of the inserted material between breakends">',
'<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">',
'<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">',
'<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of adjacency">',
'<ID=CICN,Number=2,Type=Integer,Description="Confidence interval around copy number for the segment">',
'<ID=CICNADJ,Number=.,Type=Integer,Description="Confidence interval around copy number for the adjacency">'
)
format_values=c(
'<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',
'<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">',
'<ID=CNL,Number=.,Type=Float,Description="Copy number genotype likelihood for imprecise events">',
'<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel with respect to the genome ancestor">',
'<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">',
'<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">',
'<ID=GT,Number=1,Type=String,Description="Genotype">',
'<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
'<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
'<ID=AD,Number=1,Type=Integer,Description="Read Supporting ALT">',#TCGA
'<ID=BQ,Number=1,Type=Integer,Description="Base Quality ofr Read Supporting ALT">',#TCGA
'<ID=SS,Number=1,Type=Integer,Description="Somatic Status 0=wt,1=gm,3=so,4=LOH,5=PTM/unknown">',#TCGA
'<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
'<ID=GL,Number=G,Type=Integer,Description="Genotype Likelihood">',
'<ID=lW,Number=1,Type=Integer,Description="Max Coverage Likelihood Ratio">',
'<ID=lWc,Number=1,Type=Integer,Description="Cumulative Max Coverage Likelihood Ratio">',
'<ID=lCd,Number=1,Type=Integer,Description="Deletion Insert Size Likelihood Ratio">',
'<ID=lCi,Number=1,Type=Integer,Description="Insertion Insert Size Likelihood Ratio">',
'<ID=lDl,Number=1,Type=Integer,Description="Left/Plus Anchored Hang Read Likelihood Ratio">',
'<ID=lDr,Number=1,Type=Integer,Description="Right/Minus Anchored Hang Read Likelihood Ratio">',
'<ID=RGl,Number=1,Type=Integer,Description="Left Limit of Voted Region">',
'<ID=RGr,Number=1,Type=Integer,Description="Right Limit of Voted Region">',
'<ID=BPl,Number=1,Type=Integer,Description="Softclipped Reads Left Break Point">',
'<ID=BPr,Number=1,Type=Integer,Description="Softclipped Reads Right Break Point">',
'<ID=MTD,Number=1,Type=String,Description="Calling track">',
'<ID=VAL,Number=1,Type=Float,Description="Track value">',
'<ID=STAT,Number=1,Type=Float,Description="Z score">'
)
alt_values=c(
'<ID=DEL,Description="Deletion">',
'<ID=DEL:ME:ALU,Description="Deletion of ALU element">',
'<ID=DEL:ME:L1,Description="Deletion of L1 element">',
'<ID=DUP,Description="Duplication">',
'<ID=DUP:TANDEM,Description="Tandem Duplication">',
'<ID=INS,Description="Insertion of novel sequence">',
'<ID=INS:ME:ALU,Description="Insertion of ALU element">',
'<ID=INS:ME:L1,Description="Insertion of L1 element">',
'<ID=INV,Description="Inversion">',
'<ID=CNV,Description="Copy number variable region">'
)
filter_values=c(
'<ID=q10,Description="Quality below 10">',
'<ID=s50,Description="Less than 50% of samples have data">'
)
info_keys=rep("INFO",length(info_values)) 
format_keys=rep("FORMAT",length(format_values)) 
alt_keys=rep("ALT",length(alt_values)) 
filter_keys=rep("FILTER",length(filter_values))
sv_keys=c(info_keys,format_keys,alt_keys,filter_keys)
sv_values=c(info_values,format_values,alt_values,filter_values)

sample_meta = function(id, genomes, mixture, desc){
  template=rep('##SAMPLE=<ID=%s,Genomes=%s,Mixture=%s,Description="%s">', length(id))
  paste(gettextf(template, id, genomes, mixture, desc), collapse='\n')
}

contig_meta = function(id, len, md5, species){
  template=rep('##contig=<ID=%s,length=%s,md5=%s,species="%s">', length(id))
  paste(gettextf(template, id, len, md5, species), collapse='\n')
}

# id=c("sample1","sample2")
# genomes=c("somatic","somatic,tumor")
# mixture=c(1, 0.2)
# desc=c("pure", "mixture")
# sample_meta(id, genomes, mixture, desc)

vcf_meta=function(keys, values, program="calls2vcf.R"){
  values[keys=="fileformat"]="VCFv4.1"
  values[keys=="fileDate"]=date()
  values[keys=="source"]=program
  return(paste("##",paste(keys,values,sep="="), sep="", collapse="\n"))
} #generate meta info for writing .vcf files

parse_vcf=function(comment){
  keys=new.env(); values=new.env()
  list(keys=keys,values=values)
} #parse meta info when reading .vcf files

call_meta=function(keys, values, ref_info, cmd_info){
  values[keys=="fileformat"]="VCFv4.1"
  values[keys=="fileDate"]=date()
  values[keys=="source"]=cmd_info
  values[keys=="reference"]=ref_info
  return(paste("##",paste(keys,values,sep="="), sep="", collapse="\n"))
} #modify vcf info when generating calling files

#for each region get all reads
#detach("package:swan", unload = TRUE)
#library.dynam.unload("swan", system.file(package = "swan"))
region_break=function(region, bam_file, flank, smallDel){
  #print(region)
  chr=as.character(region[1])
  start=as.numeric(region[2])-flank+1 #convert to 1-based
  end=as.numeric(region[3])+flank+1   #convert to 1-based
  if(end-start+1>smallDel){
    br_pos=as.numeric(c(-1,-1)); br_freq=as.numeric(c(0,0))
  } else {
    if(verbose) cat("==Info: chr=", chr, "start=", start, "end=", end, "\n")
    seq_info=GRanges(chr, IRanges(start, end))
    all_reads <- allFunction(seq_info, bam_file, what=what, index=bam_file)
    all_reads <- do.call("DataFrame", all_reads); N_all=nrow(all_reads);
    head_clip=sapply(str_match_all_perl(all_reads$cigar, head_pattern_new), get_al_new)
    tail_clip=sapply(str_match_all_perl(all_reads$cigar, tail_pattern_new), get_al_new)
    al_clip=sapply(str_match_all_perl(all_reads$cigar, al_pattern), get_al)
    soft_reads <- all_reads[head_clip>0 | tail_clip>0,]
    soft_head <- head_clip[head_clip>0 | tail_clip>0] #clipped at tail
    soft_tail <- tail_clip[head_clip>0 | tail_clip>0] #clipped at head
    soft_al <- al_clip[head_clip>0 | tail_clip>0]
    #keep a record of br supports in thre future
    br_distr=unlist(mapply(read_break, soft_reads$pos, soft_head, soft_tail, soft_al))
    if(is.null(br_distr)) {
      br_pos=c(-1,-1); br_freq=c(0,0)
    } else {
      br_tab=as.data.frame(table(br_distr),stringAsFactors=F)
      if(nrow(br_tab)==1) {
        br_pos=rep(as.numeric(as.character(br_tab$br_distr[1])),2);
        br_freq=rep(as.numeric(as.character(br_tab$Freq[1])),2)
      } else {
        br_pos_all=as.numeric(as.character(br_tab$br_distr[sort(br_tab$Freq,index.return=T,decreasing=T)$ix]))
        br_freq_all=as.numeric(as.character(br_tab$Freq[sort(br_tab$Freq,index.return=T,decreasing=T)$ix]))
        br_freq_2max=br_freq_all[2];
        if(br_freq_all[1]==br_freq_2max) {
          br_pos=c(min(br_pos_all[br_freq_all==br_freq_2max]),max(br_pos_all[br_freq_all==br_freq_2max]))
          br_freq=rep(br_freq_2max,2)
        } else {
          br_pos=c(max(br_pos_all[br_freq_all==br_freq_2max]),min(br_pos_all[br_freq_all==br_freq_2max]))
          if(abs(br_pos_all[1]-br_pos[1])>abs(br_pos_all[1]-br_pos[2])){
            br_pos=c(br_pos_all[1],br_pos[2]); br_freq=c(br_freq_all[1],br_freq_2max)
          } else {
            br_pos=c(br_pos_all[1],br_pos[1]); br_freq=c(br_freq_all[1],br_freq_2max)
          }
        }
      }
    }
  }
  if(br_pos[1]>br_pos[2]) { br_pos=rev(br_pos); br_freq=rev(br_freq) }
  br_info=matrix(c(br_pos,br_freq),ncol=4)
  if(verbose) cat("==Info: br_info=",br_info,"\n")
  return(br_info)
}

#for soft clip reads, index the clipping pos
read_break=function(read_pos, soft_head, soft_tail, soft_al){
  br_pos=c()
  if(soft_head!=0) br_pos=c(br_pos,read_pos)
  if(soft_tail!=0) br_pos=c(br_pos,read_pos+soft_al-1)
  return(br_pos)
}

#rotate through the region
merge_region=function(ins_regions,del_regions){
  delCluster=IRanges(start=del_regions[,2],end=del_regions[,3])
  insCluster=IRanges(start=ins_regions[,2],end=ins_regions[,3])
  #print(delCluster)
  #print(insCluster)
  ins_ovlap_del=IRanges::as.matrix(IRanges::findOverlaps(insCluster,delCluster))
  #get non ovlap ones first
  if(nrow(ins_ovlap_del)>0){
    non_ovlap_ins=ins_regions[-unique(ins_ovlap_del[,1]),]
    #print(unique(ins_ovlap_del[,2]))
    #print(del_regions)
    non_ovlap_del=del_regions[-unique(ins_ovlap_del[,2]),]
    #keep only deletion if it conflicts with some insertion
    ovlap_del=del_regions[unique(ins_ovlap_del[,2]),]
  } else {
    non_ovlap_ins=ins_regions
    non_ovlap_del=del_regions
    ovlap_del=NULL
  }
  #print("ovlap_del")
  #print(ovlap_del)
  #print("non_ovlap_ins")
  #print(non_ovlap_ins)
  #print("non_ovlap_del")
  #print(non_ovlap_del)
  #print("merge_region")
  #print(rbind(non_ovlap_del,ovlap_del,non_ovlap_ins))
  rbind(non_ovlap_del,ovlap_del,non_ovlap_ins)
}

write_call=function(breaks,regions,data,data_mean,data_sd,sample_id,adjust=0,del_limit=10){
  #print(dim(breaks))
  n_sv = nrow(breaks); if(verbose) cat("==Info: number of calls=", n_sv,"\n")
  vcf_empty = matrix(rep(rep("",length(vcf_field)+1),n_sv), ncol=length(vcf_field)+1)
  bed_empty = matrix(rep(rep("",length(bed_field)+1),n_sv), ncol=length(bed_field)+1)
  vcf_set = data.frame(vcf_empty, stringsAsFactors=F)
  bed_set = data.frame(bed_empty, stringsAsFactors=F)
  colnames(vcf_set) = c(vcf_field,sample_id)
  colnames(bed_set) = c(bed_field,sample_id)
  for(i in seq_len(n_sv)){
    #if break points is not found or scarry
    if(breaks[i,1]==-1|breaks[i,2]==-1|(breaks[i,3]<2&breaks[i,4]<2)){
      sv_raw_pos=regions[i,2]; reg_left=regions[i,2]-adjust
      sv_raw_end=regions[i,3]; reg_right=regions[i,3]-adjust
      sv_pos=regions[i,2]-adjust; sv_cipos_left=0; sv_cipos_right=scan_par$rl;
      sv_end=regions[i,3]-adjust; sv_ciend_left=0; sv_ciend_right=scan_par$rl;
      sv_win_start=floor((regions[i,2]-scan_par$start)/scan_par$stepsize)+1
      sv_win_end=floor((regions[i,3]-scan_par$start)/scan_par$stepsize)+1
    } else { # we can assume br_pos[1]<br_pos[2]
      sv_raw_pos=breaks[i,1]; reg_left=regions[i,2]-adjust
      sv_raw_end=breaks[i,2]; reg_right=regions[i,3]-adjust
      sv_pos=breaks[i,1]-adjust; sv_end=breaks[i,2]-adjust;
      sv_cipos_left=min(0,regions[i,2]-breaks[i,1]);sv_cipos_right=scan_par$rl
      sv_ciend_left=0;sv_ciend_right=max(0,regions[i,3]-breaks[i,2]+scan_par$rl)
      sv_win_start=floor((breaks[i,1]-scan_par$start)/scan_par$stepsize)+1
      sv_win_end=floor((breaks[i,2]-scan_par$start)/scan_par$stepsize)+1
    }
    del_flag=T #how to detect deletion/insertion, default is deletion
    #break point based rules: break point 1bp apart => insertion else => deletion
    if(sv_end-sv_pos<=del_limit) { #can't distinguish 5bp lower as deletion/insertion, default to ins
      del_flag=F;
    } else {
      del_flag=T;
    }
    #               2, lCi signal (weak) => insertion
    #               3, lDr and lDl peak overlap + no lCd signal => insertion (weak)
    #if(sum(data[["lCi"]][sv_win_start:sv_win_end]>2*data_sd[["lCi"]])>2) del_flag=F
    #if(sum(data$lC[sv_win_start:sv_win_end]>2*lC_sd)>support) del_flag=T
    #peak_1st=sort(data$lD[sv_win_start:sv_win_end],index.return=T,decreasing=T)$ix[1]
    #peak_2nd=sort(data$lD[sv_win_start:sv_win_end],index.return=T,decreasing=T)$ix[2]
    #if(abs(peak_2nd-peak_1st)>support) del_flag=T
    if(del_flag){
      sv_type="<DEL>"; sv_type_short="DEL"
    } else {
      sv_type="<INS>"; sv_type_short="INS"
    } #if two peaks => deletion, else => insertion
    sv_size=sv_end-sv_pos
    #cat("IMPRECISE;SVTYPE=%s;END=%d;SVLEN=%d;CIPOS=%d,%d;CIEND=%d,%d",sv_type_short,sv_end,sv_size,sv_cipos_left,sv_cipos_right,sv_ciend_left,sv_ciend_right,"\n")
    sv_info=gettextf("IMPRECISE;SVTYPE=%s;END=%d;SVLEN=%d;CIPOS=%d,%d;CIEND=%d,%d",sv_type_short,sv_end,sv_size,sv_cipos_left,sv_cipos_right,sv_ciend_left,sv_ciend_right)
    sv_id=gettextf("%s", i)
    ref_base=toString(subseq(ref_seq[[seqname]],sv_raw_pos,sv_raw_pos))
    qual=0 #maximum of all signals, standardized
    for(s in scores){
      qual=max(qual,mean(data[[s]][sv_win_start:sv_win_end]-data_mean[[s]]/data_sd[[s]]),na.rm=T) }
    qual=round(qual)
    filter="PASS"
    sv_format=paste(c(scores,"RGl","RGr","BPl","BPr"),collapse=":")
    signals=c()
    for(s in scores) {
       signals=c(signals, round((max(c(data[[s]][sv_win_start:sv_win_end],0),na.rm=T)-data_mean[[s]])/data_sd[[s]]))
    }
    sv_dummy=paste(c(signals,reg_left,reg_right,breaks[i,3],breaks[i,4]),collapse=":")
    vcf_set[i,]=c(seqname,sv_pos,sv_id,ref_base,sv_type,qual,filter,sv_info,sv_format,sv_dummy)
    bed_set[i,]=c(seqname,sv_pos,sv_end,sv_id,ref_base,sv_type,qual,filter,sv_info,sv_format,sv_dummy)
  }
  vcf_meta=paste(call_meta(
    c(meta_keys, sv_keys), c(meta_values, sv_values), ref_file, "bed2vcf.R"),
    sample_meta(sample_info[,1], sample_info[,2], sample_info[,3], sample_info[,4]),
    contig_meta(contig_info[,1], contig_info[,2], contig_info[,3], contig_info[,4]),
    gettextf("#%s", paste(colnames(vcf_set), collapse='\t')), sep='\n')
  if(verbose) cat("==Info: writing vcf file", "length=", dim(vcf_set)[1], "\n")
  return(list(vcf_set=vcf_set,vcf_meta=vcf_meta,bed_set=bed_set))
}

#cat("loaded line 3000... \n")

inGap<-function(pos,seqname,gap){
  gap=gap[gap$chrom==seqname,]
  if(is.null(gap)) return(rep(FALSE,length(pos)))
  sapply(pos,function(x,gap){any(x>gap[,2] & x<gap[,3])},gap)
}

distanceFromGap<-function(pos,seqname,gap){
  gap=gap[gap$chrom==seqname,]
  if(is.null(gap)) return(rep(Inf,length(pos)))
  sapply(pos,function(x,gap){if(any(x>gap[,2] & x<gap[,3])) return(0) else return(min(c(abs(gap[,2]-x),abs(gap[,3]-x))))},gap)
}

inf2num=function(x){ #x is a unbounded vector need to be plotted
  x[which(is.infinite(x))]=ifelse(x[is.infinite(x)]>0,double.xmax,double.xmin)
  return(x)
} #convert +/- infnity to max +/- double numbers 

pct <- function(x, digits = 2, format = "f", ...){ #x is formatted to percentage
  paste(formatC(100 * x, format = format, digits = digits, ...), "%", sep = "")
} #format a number in xx% format

nsf <- function(x){
  format(x, scientific=F)
}

write_com <- function(x, comment, filename, gzip=T, ...){
  if(gzip){
    gzf=gzfile(filename,"w") #T
  } else {
    gzf=file(filename,"w") #F 
  }
  if(!is.null(comment)) write(comment, gzf)
  write.table(x, gzf, ...)          #write content in non-scientific format
  close(gzf)
} #write a table with comment, for example meta info in vcf file

read_com <- function(filename, comment_char='#', ...){
  tryCatch({
    if(!file.exists(filename)){
      stop("====Error: ", filename, "does not exist, bail out now\n")
    }  
    content=readLines(filename, n=10000) #max comments 10000
    if(is.na(content[length(content)][1]))
      stop("====Error: ", filename, "is not correctly formatted, bail out now\n")
    if(content[length(content)][1]==comment_char){ 
      cat("on reading file:",filename,"\n")
      cat("the final, which is the ",length(content)," -th line:",content[length[content]][1]," contains comment_char", comment_char, "\n")
      stop("====Error: maximum allowed comment lines 10000, over that truncated\n")
    }
    comment=content[substring(content,1,1)==comment_char]
    data=read.delim(filename, comment.char=comment_char, ...) #TODO: this sometimes halt without error, need find safe alternative
  },error=function(x) {stop("read_com: ",filename," error!\n")}) 
  return (list(comment=comment, data=data))
} #read a table with comment, for example meta info in vcf file

### spike-in functions ###
## get trunk out of reference, save trunk1
spike_sv<-function(ref,pos,w,sv){
  
  ref_seq=ref[[1]]
  if(w<0){ #insertion at target
    sved=DNAStringSet(c(subseq(ref_seq,1,pos), sv, 
                        subseq(ref_seq,pos+1,length(ref_seq))))  #target
    sv=DNAStringSet(sv)
  } #insertion right after pos in target
  else if(w>0){ #deletion at target
    sved=DNAStringSet(c(subseq(ref_seq,1,pos), 
                        subseq(ref_seq,pos+1+w,length(ref_seq)))) #target
    sv=DNAStringSet(subseq(ref_seq,pos+1,pos+w))
  } #deletion right after pos in target
  else{
    sved=ref_seq
    sv=DNAStringSet(subseq(ref_seq,pos+1,pos)) #"this is just placeholder"
  } #no-sv
  return(list(ref=ref,sved=sved,sv=sv,w=w))
}

ref_ends <- function(ref_seq){ 
  head = regexpr("[ACGTacgt]", ref_seq, perl=TRUE) # first match 
  tail = length(ref_seq)+1-regexpr("[ACGTacgt]", reverse(ref_seq), perl=TRUE) 
  list(head=as.numeric(head),tail=as.numeric(tail))
} # determine the non-N head and tail of ref_seq

###
countOverlaps_omp = function(win_intvl, read_intvl){
  #cat("in countOverlaps_omp ..."); ptm=proc.time(); pmm=gc();
  #cat("out countOverlaps_omp ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  db_start = start(read_intvl); db_end=end(read_intvl);
  query_start = start(win_intvl); query_end=end(win_intvl);
  #overlapCount = .Call("countOverlaps_omp", db_start, db_end, query_start, query_end, PACKAGE="swan");
  #gced=gc(); cat("out countOverlaps_omp ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(overlapCount)
}

findOverlaps_within = function(qry_intvl, sbj_intvl){ #find qry_intvl(query) complete within sbj_intvl(subject)
  result=IRanges::findOverlaps(qry_intvl,sbj_intvl,type="within")
  #IRanges::.IntervalTreeCall(subject, "overlap_all", query, query_ord) #query[query_ord] must be sorted  
  #sbj_tree=IntervalTree(sbj_intvl)
  #result = .Call("IntegerIntervalTree_overlap_all", sbj_tree@ptr, qry_intvl, seq_along(qry_intvl), PACKAGE="IRanges") 
  ##win_intvl[seq_along_win_intvl] must be sorted  
  #m <- IRanges::as.matrix(result)
  #r <- ranges(result, qry_intvl, sbj_intvl)
  #m <- m[width(qry_intvl)[m[,1L]] - width(r) <= 0, , drop=FALSE]
  ### unname() required because in case 'm' has only 1 row
  ### 'm[ , 1L]' and 'm[ , 2L]' will return a named atomic vector
  #result@queryHits <- unname(m[ , 1L])
  #result@subjectHits <- unname(m[ , 2L])
  return(result)
}

findOverlaps_within_old = function(qry_intvl, sbj_intvl){ #find qry_intvl(query) complete within sbj_intvl(subject)
  #IRanges::.IntervalTreeCall(subject, "overlap_all", query, query_ord) #query[query_ord] must be sorted  
  sbj_tree=IntervalTree(sbj_intvl)
  result = .Call("IntegerIntervalTree_overlap_all", sbj_tree@ptr, qry_intvl, seq_along(qry_intvl), PACKAGE="IRanges") 
  #win_intvl[seq_along_win_intvl] must be sorted  
  m <- IRanges::as.matrix(result)
  r <- ranges(result, qry_intvl, sbj_intvl)
  m <- m[width(qry_intvl)[m[,1L]] - width(r) <= 0, , drop=FALSE]
  ## unname() required because in case 'm' has only 1 row
  ## 'm[ , 1L]' and 'm[ , 2L]' will return a named atomic vector
  result@queryHits <- unname(m[ , 1L])
  result@subjectHits <- unname(m[ , 2L])
  return(result)
}

### bam scan functions ###
lC_R = function(winC_rPb, fy, n_wins, mixing_rate){
  #cat("in lC_R ..."); ptm=proc.time(); pmm=gc();
  lC_term = rep(0,n_wins); i=1; N=nrow(winC_rPb); if(is.null(N)) N=0;   
  while(i<=N){
    lC_term[winC_rPb[i,1]]=lC_term[winC_rPb[i,1]]+log((1-mixing_rate)+mixing_rate*fy[winC_rPb[i,2]])
    i=i+1
  }
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out lC_R ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(lC_term)
}

lC_omp <- function(winC_rPb, rPb_isize, fy, n_win, mixing_rate){
  #cat("in lC_omp ..."); ptm=proc.time(); pmm=gc();
  lC_term=.Call( "lC_omp", winC_rPb, rPb_isize, fy, n_win, mixing_rate, PACKAGE="swan")
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out lC_omp ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(lC_term)
}

lD_R = function(rS_winD, winD, rS, Fx, n_wins, r, p, left=TRUE) { #left/right = winDl/winDr
  #cat("in lD_R ..."); ptm=proc.time(); pmm=gc();
  lD_term = rep(0, n_wins); i=1; N=nrow(rS_winD); if(is.null(N)) N=0; sign=ifelse(left,1,-1)
  while(i<=N){
    lD_term[rS_winD[i,2]]=lD_term[rS_winD[i,2]]+log(1+r*(1-p)/p*Fx[sign*(winD[rS_winD[i,2]]-rS[rS_winD[i,1]])+1])
    i=i+1
  }
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out lD_R ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(lD_term)
}

lD_omp <- function(rS_winD, winD, rS, Fx, n_wins, r, p, left=TRUE){
  #cat("in lD_omp ..."); ptm=proc.time(); pmm=gc();
  lD_term=.Call( "lD_omp", rS_winD, winD, rS, Fx, n_wins, r, p, left, PACKAGE="swan")
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out lD_omp ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(lD_term)
}

aggregate_sum_omp <- function(coverage, start, end){
  #cat("in aggregate_sum_omp ..."); ptm=proc.time(); pmm=gc();
  aggregate_sum=.Call( "aggregate_sum_omp", coverage, start, end )
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out aggregate_sum_omp ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(aggregate_sum)
}

anch_clust <- function(Lfwd, Lchr, Lst, Led, Rfwd, Rchr, Rst, Red, chr, chr_size){
  anchor=.Call( "anch_clust", Lfwd, Lchr, Lst, Led, Rfwd, Rchr, Rst, Red, chr, chr_size )
  return(anchor)
}

#' This gonna do substractive clustering between spX (case) and spY (control) of read pairs
#' 
#' x is a case sample,  y is a control sample 
#' L is for left read, R is for right read
#' fwd is true if direction is forward, false otherwise
#' st is the start position of a read, ed is the end position of the read
#' sup is the number of read pairs support this anchor cluster, here is always one (because it is a rdp)
contrast_clust <- function(spX_anch_data,spY_anch_data,chr,chr_size,spX_opt,spY_opt){
  spX_sup=as.integer(spX_opt$sup)
  spY_sup=as.integer(spY_opt$sup)
  #print(spX_sup); print(spY_sup)
  xLfwd=spX_anch_data[,1]=="+"; xLchr=as.character(spX_anch_data[,2]); xLst=spX_anch_data[,3]; xLed=spX_anch_data[,4];
  xRfwd=spX_anch_data[,5]=="+"; xRchr=as.character(spX_anch_data[,6]); xRst=spX_anch_data[,7]; xRed=spX_anch_data[,8];
  yLfwd=spY_anch_data[,1]=="+"; yLchr=as.character(spY_anch_data[,2]); yLst=spY_anch_data[,3]; yLed=spY_anch_data[,4];
  yRfwd=spY_anch_data[,5]=="+"; yRchr=as.character(spY_anch_data[,6]); yRst=spY_anch_data[,7]; yRed=spY_anch_data[,8];
  #print(xLfwd); print(xLchr); print(xLst); print(xLed); print(xRfwd); print(xRchr); print(xRst); print(xRed);
  #print(yLfwd); print(yLchr); print(yLst); print(yLed); print(yRfwd); print(yRchr); print(yRst); print(yRed);
  anchor=.Call( "contrast_clust", xLfwd,  xLchr,  xLst,  xLed,  xRfwd,  xRchr,  xRst,  xRed,
                     yLfwd,  yLchr,  yLst,  yLed,  yRfwd,  yRchr,  yRst,  yRed,
                     chr,  chr_size, spX_sup, spY_sup )
  return(anchor)
}

# lDl_omp <- function(rSr_winDl, winW_start, rS_right, Fx, n_wins, r, p){
#   .Call( "lDl_omp", rSr_winDl, winW_start, rS_right, Fx, n_wins, r, p, PACKAGE="swan")
# }
# 
# lDr_omp <- function(rSl_winDr, winW_end, rS_left, Fx, n_wins, r, p){
#   .Call( "lDr_omp", rSl_winDr, winW_end, rS_left, Fx, n_wins, r, p, PACKAGE="swan")
# }

statFunction <- function(seqname, bamFile, ...){
  param <- ScanBamParam(what = what,
                        which = GRanges(seqname, IRanges(1,max_chr_len)),
                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  x <- scanBam(bamFile, ..., param = param)[[1]]
  list(cvg=coverage(IRanges(x[["pos"]], width = x[["qwidth"]])),
       rl=x[["qwidth"]],
       isize=x[["isize"]])
} # used in bam_stat.R

allFunction <- function(seqinfo, bamFile, what, ...){
  param <- ScanBamParam(what=what,
                        which=seqinfo,
                        flag=scanBamFlag(isDuplicate=F,isNotPassingQualityControls=F))
  x <- scanBam(bamFile, ..., param = param)[[1]]
  return(x)
} # used in scan.R

parFunction <- function(bamFile, param, ...){
  #ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE,
  #                     reverseComplement = FALSE, tag = character(0),
  #                             what = character(0), which)
  x <- scanBam(bamFile, ..., param = param)[[1]]
  return(x)
} # used in scan.R


### swan functions ###

fIw<-function(y,w,is,issd){     #w=wi-is 
  exp(w*(2*(y-is)-w)/(2*issd^2))
} ## fI(w) function

fI<-function(y,w,is,issd){
  exp(w*(2*(y-is)-w)/(2*issd^2))
}

FI<-function(y,is,issd){
  pnorm(y,mean=is,sd=issd,lower.tail=FALSE)
} ## FI(x) function

im3 <- function(mat1, mat2){
  stopifnot(ncol(mat1)==2, ncol(mat1)==ncol(mat2))
  s=new(HashSet)
  for(i in seq_len(nrow(mat1))) s$set(mat1[i,1], mat1[i,2])
  idx=rep(0,nrow(mat2))
  for(i in seq_len(nrow(mat2))) idx[i]=s$get(mat2[i,1],mat2[i,2])
  matrix(mat2[idx==1],ncol=2)
} # use HashSet to intersect two matrices, replacing im2 and im1
#speed testing
#size=100000
#A  <- matrix(sample(1:1000, 2*size, replace=TRUE), size, 2)
#B  <- matrix(sample(1:1000, 2*size, replace=TRUE), size, 2)
#system.time(n <- im3(A, B))
#size=10000000, time=
#size=1000000, time=25s
#size=100000, time=1.3s

z_scan <- function(x){
  scale_switch=TRUE
  if(length(unique(x$logL))!=1) x$ZlogL=scale(inf2num(x$logL)) else x$ZlogL=0
  if(length(unique(x$logLdel))!=1) x$ZlogLdel=scale(inf2num(x$logLdel)) else x$ZlogLdel=0
  if(length(unique(x$logLins))!=1) x$ZlogLins=scale(inf2num(x$logLins)) else x$ZlogLins=0
  if(length(unique(x$lW))!=1) x$ZlW=scale(inf2num(x$lW), scale=scale_switch) else x$ZlW=0
  if(length(unique(x$lC))!=1) x$ZlC=scale(inf2num(x$lC), scale=scale_switch) else x$ZlC=0
  if(length(unique(x$lDl))!=1) x$ZlDl=scale(inf2num(x$lDl), scale=scale_switch) else x$ZlDl=0
  if(length(unique(x$lDr))!=1) x$ZlDr=scale(inf2num(x$lDr), scale=scale_switch) else x$ZlDr=0
  if(length(unique(x$lD))!=1) x$ZlD=scale(inf2num(x$lD), scale=scale_switch) else x$ZlD=0
  return(x)
} # scale and center the LR scores

sv_call<-function(score_array, score_label, thresh){
  call_idx=which(score_array[[score_label]]>=thresh)
  call_wins=IRanges(start=score_array$start[call_idx], end=score_array$end[call_idx])
} 
# call SVs by reducing connecting high score windows

get_al=function(match){ return(sum(as.numeric(match[,2]))) } #if not found will return zero
# alignment length summary for results from str_match_all_perl

get_al_new=function(match){ return(
    sum(ifelse(is.na(as.numeric(match[,2])),0,as.numeric(match[,2]))) 
    +sum(ifelse(is.na(as.numeric(match[,3])),0,as.numeric(match[,3])))) 
} #return sum of size1 and size2

get_al_old=function(cigar){
  pattern="(?<size>[0-9]+)[MX=IN]"
  match=str_match_all_perl(cigar,pattern)
  return(sapply(match,get_al))
}

str_match_all_perl <- function(string,pattern){
  string[is.na(string)] = ""
  parsed <- gregexpr(pattern,string,perl=TRUE)
  #print(parsed)
  lapply(seq_along(parsed),function(i){
    r <- parsed[[i]]
    starts <- attr(r,"capture.start")
    #if(is.na(r[1])) return(matrix(nrow=0,ncol=1+ncol(starts)))
    if(r[1]==-1) return(matrix(nrow=0,ncol=1+ncol(starts)))
    names <- attr(r,"capture.names")
    lengths <- attr(r,"capture.length")
    full <- substring(string[i],r,r+attr(r,"match.length")-1)
    subs <- substring(string[i],starts,starts+lengths-1)
    m <- matrix(c(full,subs),ncol=length(names)+1)
    colnames(m) <- c("",names)
    m
  })
} 
# accelerated string matchwith multiple patterns

get_al_omp <- function(cigar, al_string){
  al=.Call( "get_al_omp", cigar, al_string, PACKAGE="swan" )
  return(as.vector(al))
}

mate_MPRs_omp <- function(h1st_idx, c1st_MPRs, c2nd_MPRs, impute=F, left=T, self=T){
  mate_MPRs=.Call( "mate_MPRs_omp", h1st_idx, 
                 c1st_MPRs$qname, c1st_MPRs$cigar, c1st_MPRs$isize, 
                 c2nd_MPRs$qname, c2nd_MPRs$cigar, c2nd_MPRs$isize, 
                 impute, left, self, PACKAGE="swan")
  #elapsed=proc.time()-ptm; gced=gc(); 
  return(mate_MPRs)
}

#write a fully cpp function mate_MPRs_cpp
mate_MPRs_hash <- function(h1st_idx, c1st_MPRs, c2nd_MPRs, impute=F, left=T, self=T){ #default: input is left, impute self
  #input=left, impute=mate | input=right, impute=self
  #cat("in mate_MPRs_hash ..."); ptm=proc.time(); pmm=gc();
  p1st=tail_pattern_new; p2nd=head_pattern_new; impute_isize=NA;
  if(impute) impute_isize=rep(NA,length(h1st_idx))
  if(left) { p1st=head_pattern_new; p2nd=tail_pattern_new }
  qname_h1st = new(HashMap) #new.env(hash=h) #the hash way
  for(idx in h1st_idx) { qname_h1st$set(c1st_MPRs$qname[idx],idx) } #hang_idx collect all reads hang and record name
  #cat("done build qname_h1st ...")
  mh1st_idx=rep(0,length(h1st_idx)); j=1 #only corret isize
  #cat("length of c2nd_MPRs", nrow(c2nd_MPRs))
  for(midx in seq_len(nrow(c2nd_MPRs))){
    idx=qname_h1st$get(c2nd_MPRs$qname[midx])
    #if(midx%%1000==0) cat(midx, "of", nrow(c2nd_MPRs), "j=", j, "\n")
    if(idx!=0){ #midx is a mate of idx
      #if(j==1) cat("midx=",midx,c2nd_MPRs$qname[midx],"idx=",idx,c1st_MPRs$qname[idx],"\n")
      if(impute){
        clip_1st = sapply(str_match_all_perl(c1st_MPRs$cigar[idx], p1st), get_al_new)
        clip_2nd = sapply(str_match_all_perl(c2nd_MPRs$cigar[midx], p2nd), get_al_new)      
        if(left & self){ impute_isize[j]=c1st_MPRs$isize[idx]+clip_1st+clip_2nd } #impute a left self
        if(left & !self){ impute_isize[j]=c2nd_MPRs$isize[midx]-clip_1st-clip_2nd } #impute a left mate
        if(!left & self){ impute_isize[j]=c1st_MPRs$isize[idx]-clip_1st-clip_2nd } #impute a right self
        if(!left & !self){ impute_isize[j]=c2nd_MPRs$isize[midx]+clip_1st+clip_2nd } #impute a right mate
      }
      mh1st_idx[j]=midx; j=j+1; #however, mh1st_idx[j] may not be the mate of h1st_idx[j]
    }
  } #there is the possibility hang_idx overlap with mhang_idx, below we only keep single hang, both hang is thrown
  stopifnot(any(!is.na(mh1st_idx))) #ensure all input has a mate
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out mate_MPRs_hash ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(list(mate_idx=mh1st_idx, impute_isize=impute_isize))
}

mate_MPRs_list <- function(h1st_idx, c1st_MPRs, c2nd_MPRs, impute=F, left=T, self=T, hash=T){ #input=left/right, impute=mate/self
  #input=left, impute=mate | input=right, impute=self
  #cat("in mate_MPRs_list ..."); ptm=proc.time(); pmm=gc();
  #head_pattern="^(?<size>[0-9]+)S"; tail_pattern="(?<size>[0-9]+)S$"; 
  p1st=tail_pattern_new; p2nd=head_pattern_new; impute_isize=NA;
  if(impute) impute_isize=rep(NA,length(h1st_idx))
  if(left) { p1st=head_pattern_new; p2nd=tail_pattern_new }
  qname_h1st = new.env(hash=hash) #the hash way
  for(idx in h1st_idx) { qname_h1st[[c1st_MPRs$qname[idx]]]=idx } #hang_idx collect all reads hang and record name
  mh1st_idx=rep(0,length(h1st_idx)); j=1 #only corret isize
  for(midx in seq_len(nrow(c2nd_MPRs))){
    idx=qname_h1st[[c2nd_MPRs$qname[midx]]]
    #if(midx%%1000==0) cat(midx, "of", nrow(c2nd_MPRs), "j=", j, "\n")
    if(!is.null(idx)){ #midx is a mate of idx
      if(impute){
        clip_1st = sapply(str_match_all_perl(c1st_MPRs$cigar[idx], p1st), get_al_new)
        clip_2nd = sapply(str_match_all_perl(c2nd_MPRs$cigar[midx], p2nd), get_al_new)      
        if(left & self){ impute_isize[j]=c1st_MPRs$isize[idx]+clip_1st+clip_2nd } #impute a left self
        if(!left & !self){ impute_isize[j]=c2nd_MPRs$isize[midx]+clip_1st+clip_2nd } #impute a left mate
        if(!left & self){ impute_isize[j]=c1st_MPRs$isize[idx]-clip_1st-clip_2nd } #impute a right self
        if(left & !self){ impute_isize[j]=c2nd_MPRs$isize[midx]-clip_1st-clip_2nd } #impute a right mate
      }
      mh1st_idx[j]=midx; j=j+1; #however, mh1st_idx[j] may not be the mate of h1st_idx[j]
    }
  } #there is the possibility hang_idx overlap with mhang_idx, below we only keep single hang, both hang is thrown
  stopifnot(any(!is.na(mh1st_idx))) #ensure all input has a mate
  #elapsed=proc.time()-ptm; gced=gc(); 
  #cat("out mate_MPRs_list ...",elapsed[3],"s", (gced[1,2]+gced[2,2])-(pmm[1,2]+pmm[2,2]), colnames(pmm)[2], "\n")
  return(list(mate_idx=mh1st_idx, impute_isize=impute_isize))
}

mate_strand=function(flags,strands){
  ifelse(flags %in% disc_flags, strands, switch_strand(strands))
}

switch_strand=function(strands){
  ifelse(strands=="+","-","+")
}

#' @title get_MPRs
#' 
#' @description 
#' get_MPRs will return anchor read pairs
#' 
#' @param seq_info
#' @param bam_file
#' @param RL
#' @param Delta
#' @param bigDel
#' @param smallDel
#' @param smallIns
#' @param hang_clip
#' @param prop_clip
#' @param wplus
#' 
#' @note \code{get_MPRs} 
#' \cr
#' If called, \code{funs} get_MPRs will return anchor read pairs
#' 
#' @examples
#' get_MPRs(seq_info=GenomicRanges("11",IRanges(start=1,end=1000000)),
#'              bam_file="NA12878.mem.bam",
#'              RL=100,
#'              Delta=500,
#'              bigDel=1200,
#'              smallDel=50,
#'              smallIns=50,
#'              hang_clip=50,
#'              prop_clip=50,
#'              wplus=10,
#'              )
#'              
#' @seealso \code{\link{scan_joint}} and \code{\link{scan_bam}}
get_MPRs <- function(seq_info,bam_file,RL,Delta,bigDel,smallDel,smallIns,maxInsert,hang_clip,prop_clip,wplus,gtk=NULL,debug=F,verbose=F,index=bam_file){ 
  # this function replace get_MPRs and reduce to only one round of scan
  # let's use rSC to store clipping point
  # clipping theory:
  #    1. read mapped to forward and clipped on left  -> rSCpL, flags=rmfw_flags
  #    2. read mapped to forward and clipped on right -> rSCpR, flags=rmfw_flags
  #   3. read mapped to reverse and clipped on left  -> rSCnL, flags=rmrv_flags
  #   4. read mapped to reverse and clipped on right -> rSCnR, flags=rmrv_flags
  #   a read can have both but should be rare
  #TODO: output rSCpL, rSCpR, rSCnL, rSCnR replacing rSCHp, rSCHn, rSCISp, rSCISn
  #   the softclipping information can help genotyping
  if(verbose) cat("---Info: in get_MPRs ...\n")
  if(verbose) cat("---Info: input options ...\n")
  if(verbose) cat(paste(c("bam_file","RL","Delta","bigDel","smallDel","smallIns","hang_clip","prop_clip","maxInsert"),sep="\t"),"\n")
  if(verbose) cat(paste(c(bam_file,RL,Delta,bigDel,smallDel,smallIns,hang_clip,prop_clip,maxInsert),sep="\t"),"\n")
  if(verbose) cat("===Info: seq_info\n")
  if(verbose) print(seq_info)
  seq_name=IRanges::as.vector(seqnames(seq_info))[1]
  all_what=c("qname","pos","mrnm","mpos","isize","cigar","strand","flag","qwidth")
  all_param=ScanBamParam(flag = scanBamFlag(isDuplicate=F, isNotPassingQualityControls=F),
                            simpleCigar=F, reverseComplement = F, tag = character(0),
                            what = all_what, which=seq_info) #read in all pass qual reads
  all_MPRs=as.data.frame(parFunction(bam_file, all_param), stringsAsFactors=F)
  #all_flags=parseFlags(all_MPRs$flag) #all_flags
  N_all=nrow(all_MPRs)
  #all initial values for varaibles
  rSl=IRanges();rSr=IRanges();rHr=IRanges();rHl=IRanges();rSDa=IRanges();rSIa=IRanges();
  rPbi=IRanges();rPbo=IRanges();rPbd=IRanges();rPc=IRanges();rMc=IRanges();
  
  rSCPbi=IRanges();rSCPbo=IRanges();
  rPn=0;rMn=0;rMHp=IRanges();rMHn=IRanges();rMDp=IRanges();rMDn=IRanges();
  rSCHp=IRanges();rSCHn=IRanges();rSCISp=IRanges();rSCISn=IRanges()
  hend_MPRs=data.frame();disc_MPRs=data.frame();
  soft_MPRs=data.frame();soft_lMPRs=data.frame();soft_rMPRs=data.frame();
  hang_rMPRs=data.frame();hang_lMPRs=data.frame();merge_MPRs=data.frame();
  anchor=data.frame() #left_strand, left_chr, left_start, left_end, right_strand, right_chr, right_start, right_end
  isize=NULL;isize_sd=NULL;

  if(verbose) cat("===Info: N_all=",N_all,"\n")
  if(N_all<min_trunk_MPRs) { 
    cat("---Warn: this trunk has too small number of MPRs for inference!\n") 
    return(list(
    rPbi=rPbi,rPbo=rPbo,rPbd=rPbd,rPc=rPc,rMc=rMc,rPn=rPn,rMn=rMn,
    rSr=rSr,rSl=rSl,rHr=rHr,rHl=rHl,
    rMHp=rMHp,rMHn=rMHn,rMDp=rMDp,rMDn=rMDn,rSCHp=rSCHp,rSCHn=rSCHn,
    rSCISp=rSCISp,rSCISn=rSCISn,rSDa=rSDa,rSIa=rSIa,rSCPbi=rSCPbi,rSCPbo=rSCPbo,
    hang_clip=NA,prop_clip=NA,
    RL=NA,isize=NA,isize_sd=NA,Delta=NA,
    smallDel=NA,smallIns=NA,bigDel=NA,maxInsert=NA,anchor=data.frame()))
  }

  all_MPRs$pos <- as.integer(all_MPRs$pos)
  all_MPRs$mpos <- as.integer(all_MPRs$mpos)
  all_MPRs$isize <- as.integer(all_MPRs$isize)
  all_MPRs$qwidth <- as.integer(all_MPRs$qwidth)
  all_MPRs$flag <- as.character(all_MPRs$flag)
  isize=median(abs(all_MPRs$isize[all_MPRs$flag %in% conc_flags & abs(all_MPRs$isize)>RL]),na.rm=T)
  isize_sd=round(mad(abs(all_MPRs$isize[all_MPRs$flag %in% conc_flags & abs(all_MPRs$isize)>RL])))
  if(Delta==-1) Delta=min(max(minDelta_fail,isize+Delta_sd*isize_sd,na.rm=T),maxDelta_fail)
  if(bigDel==-1) bigDel=min(max(minbigDel_fail,isize+bigDel_sd*isize_sd,na.rm=T),maxbigDel_fail)
  if(RL==-1) RL=median(all_MPRs$qwidth[all_MPRs$flag %in% conc_flags],na.rm=T)
  if(maxInsert==-1) maxInsert=bigDel+20000; minInsert=RL #20000 is buffer zone between bigDel and maxInsert
  #we define pos<=mpos as left and mpos<=pos as right regardless of rname and mrnm
  #disc marks if orientation is  fr or rf pairs
  #impp makrs if isize is NA or over large
  #we supllement these by concordant pairs that have bigDel embedded or these repairs actually collides
  anch_Lidx = (all_MPRs$mpos>=all_MPRs$pos) & (all_MPRs$flag %in% disc_flags | all_MPRs$flag %in% impp_flags | ((all_MPRs$flag %in% conc_flags) & (all_MPRs$isize<RL | all_MPRs$isize>bigDel))) #this is to select all anchor read pairs, select only if the read pair is marked as discordant or the read pair is marked improper or the read pair is concordant but the insert size is either too large or too small. And if the mate read has a mpos larger than the pos of current read, which means current read is to the Left of the mate read, and this anchor is registered in anch_Lidx; and similarily for anch_Ridx if its to the right of its mate.
  anch_Lidx[is.na(anch_Lidx)]=FALSE
  anch_Ridx = (all_MPRs$pos>=all_MPRs$mpos) & (all_MPRs$flag %in% disc_flags | all_MPRs$flag %in% impp_flags | ((all_MPRs$flag %in% conc_flags) & (all_MPRs$isize>-RL | all_MPRs$isize<(-bigDel))))
  anch_Ridx[is.na(anch_Ridx)]=FALSE
  field_idx = !is.na(all_MPRs$mrnm) & !is.na(all_MPRs$strand) # check if the rdp has the full iformation needed, otherwsie would discard.
  tmp_lMPRs = all_MPRs[anch_Lidx & field_idx, c("strand","pos","mrnm","mpos","qwidth","flag")] 
  tmp_rMPRs = all_MPRs[anch_Ridx & field_idx, c("strand","pos","mrnm","mpos","qwidth","flag")] 
  Lanch = cbind(as.character(tmp_lMPRs$strand),rep(seq_name,nrow(tmp_lMPRs)),tmp_lMPRs$pos,tmp_lMPRs$pos+tmp_lMPRs$qwidth,mate_strand(tmp_lMPRs$flag,as.character(tmp_lMPRs$strand)),tmp_lMPRs$mrnm,tmp_lMPRs$mpos,tmp_lMPRs$mpos+RL) #left registered anchors
  Ranch = cbind(mate_strand(tmp_rMPRs$flag,as.character(tmp_rMPRs$strand)),tmp_rMPRs$mrnm,tmp_rMPRs$mpos,tmp_rMPRs$mpos+RL,as.character(tmp_rMPRs$strand),rep(seq_name,nrow(tmp_rMPRs)),tmp_rMPRs$pos,tmp_rMPRs$pos+tmp_rMPRs$qwidth) #right registered anchors
  anchor = rbind(Lanch,Ranch)
  if(nrow(anchor) != 0)
    names(anchor) = c("Lstrand","Lchr","Lstart","Lend","Rstrand","Rchr","Rstart","Rend")
  #set_colClass(anchor,rep(c("character","character","integer","integer"),2))
  if(verbose) cat("---Info: constructed anch_MPRs, row=",nrow(anchor),"col=",ncol(anchor),"\n")
  #anch_left = L & ( disc_flag | (impp_flag | conc_flag) & (isize<RL | isize>bigDel) ) 
  #anch_right = R & ( disc_flag | (impp_flag | conc_flag) & (isize>RL | isize<-bigDel) )
  #anch_1st will be all reads that is a left read
  #we need to deduct the orientation of the other end
  #if discordant is set, that means either ff or rr
  #if discordant is not set, that means either fr or rf
  #we exclude anything that either mpos is missing or mrnm is missing
  #we need to include flag though

  if(verbose) cat("===Info: concordant MPRs",sum(all_MPRs$flag %in% conc_flags),"\n")
  if(verbose) cat("===Info: improper MPRs",sum(all_MPRs$flag %in% impp_flags),"\n")

  if(prop_clip==-1) prop_clip=as.integer(RL*prop_clip_fail) #use -1 to auto, use RL+1 to disable 
  if(hang_clip==-1) hang_clip=RL-soft_cut_fail              #use -1 to auto, use 0 to disable
  
  del_size=get_al_omp(all_MPRs$cigar,ald_pattern)
  ins_size=get_al_omp(all_MPRs$cigar,ali_pattern)
  if(smallDel<0) smallDel=min(smallDel_fail,round(.2*RL))
  if(smallIns<0) smallIns=min(smallIns_fail,round(.2*RL))
  del_aidx=which(del_size>=smallDel); ins_aidx=which(ins_size>=smallIns)
  #print(strsplit(all_MPRs$cigar[del_aidx],"D"))
  if(length(del_aidx)>0){
    del_head_cigar=sapply(strsplit(all_MPRs$cigar[del_aidx],"D"),"[",1)
    del_al_head=get_al_omp(del_head_cigar,alg_pattern)
    rSDa=IRanges(start=all_MPRs$pos[del_aidx]+del_al_head,width=del_size[del_aidx])
  }
  if(length(ins_aidx)>0){
    ins_head_cigar=sapply(strsplit(all_MPRs$cigar[ins_aidx],"I"),"[",1)
    ins_al_head=get_al_omp(ins_head_cigar,alg_pattern)
    rSIa=IRanges(start=all_MPRs$pos[ins_aidx]+ins_al_head,width=smallIns)
  }
  #cat(length(del_aidx),length(del_al_head),length(del_size),length(ins_aidx),length(ins_al_head),"\n")
  #del_tail_cigar=sapply(strsplit(all_MPRs$cigar[del_aidx],"D"),"[",-1)
  #del_al_tail=get_al_omp(del_tail_cigar,alg_patern)
  #ins_tail_cigar=sapply(strsplit(all_MPRs$cigar[ins_aidx],"I"),"[",-1)
  #ins_al_tail=get_al_omp(ins_tail_cigar,alg_patern)
  #we simply assume each read can only contain one smallDel or smallIns and we keep read (+/- insertion as evidence)
  
  #here we deal with the difficulty with *improper pairs*, see impp_flags=c(97,145,81,161):
  #i.e. too close, too far, and cross cochromosome
  #all cases for impp: two far apart, in reverse order or in two chromosome
  #we can use two far apart in conc, use reverse order in disc, and use cross chromosome in hend
  #Fidx,Ridx,Aidx: indices for Forward, Reverse and All reads
  ### 1. FR still in correct order ###
  impp_conc_Fidx=all_MPRs$flag %in% c(97,161) & all_MPRs$pos<=all_MPRs$mpos & all_MPRs$mrnm==seq_name
  impp_conc_Fidx[is.na(impp_conc_Fidx)]=FALSE
  #cat("read in impp_conc_Fidx:",any(all_MPRs[impp_conc_Fidx,]$pos==27436679),"\n")
  impp_conc_Ridx=all_MPRs$flag %in% c(81,145) & all_MPRs$pos>=all_MPRs$mpos & all_MPRs$mrnm==seq_name 
  impp_conc_Ridx[is.na(impp_conc_Ridx)]=FALSE
  ### 2. FR in wrong order, basically RF ###
  #cat("read in impp_conc_Ridx:",any(all_MPRs[impp_conc_Ridx,]$pos==27437485),"\n")
  impp_disc_Fidx=all_MPRs$flag %in% c(97,161) & all_MPRs$pos>all_MPRs$mpos & all_MPRs$mrnm==seq_name #(- +) pairs
  impp_disc_Fidx[is.na(impp_disc_Fidx)]=FALSE
  impp_disc_Ridx=all_MPRs$flag %in% c(81,145) & all_MPRs$pos<all_MPRs$mpos & all_MPRs$mrnm==seq_name
  impp_disc_Ridx[is.na(impp_disc_Ridx)]=FALSE
  ### 3. FR cross chromosome ###
  impp_hend_Aidx=all_MPRs$flag %in% c(97,161,81,145) & (all_MPRs$mrnm!=seq_name | is.na(all_MPRs$mrnm)) 
  impp_hend_Aidx[is.na(impp_hend_Aidx)]=FALSE

  #handling *concordant pairs*
  #to consider, how do we decide orientation of read, by mpos and pos or by isize?
  #conc i.e. seq_name the same, concordant must have both read mapped to the same chromosome in order
  conc_Fidx=all_MPRs$flag %in% conc_flags & all_MPRs$mpos>=all_MPRs$pos & all_MPRs$mrnm==seq_name  #this already include 0<=isize<RL
  conc_Fidx[is.na(conc_Fidx)]=FALSE
  conc_Ridx=all_MPRs$flag %in% conc_flags & all_MPRs$pos>=all_MPRs$mpos & all_MPRs$mrnm==seq_name  #this already include -RL<isize<=0
  conc_Ridx[is.na(conc_Ridx)]=FALSE

  #handling *discordant pairs*
  #but discordant pairs need to be analyzed based on its mate
  #if on the same chromosome go disc track, else merge to hend_MPRs
  #if we use in disc track, we need to take start=max(pos,mpos) and end=(mpos,pos)
  disc_Aidx=all_MPRs$flag %in% disc_flags & all_MPRs$mrnm==seq_name
  disc_Aidx[is.na(disc_Aidx)]=FALSE
  disc_hend_Aidx=all_MPRs$flag %in% disc_flags & (all_MPRs$mrnm!=seq_name | is.na(all_MPRs$mrnm))
  disc_hend_Aidx[is.na(disc_hend_Aidx)]=FALSE
  disc_MPRs=all_MPRs[disc_Aidx|impp_disc_Fidx|impp_disc_Ridx, # Fidx and Ridx should be moved to new category revs_MPRs
                      -which(names(all_MPRs) %in% c("qwidth"))]
  rMDp=IRanges(start=disc_MPRs$pos[disc_MPRs$strand=="+"&!is.na(disc_MPRs$pos)],width=RL) #right pos
  rMDn=IRanges(start=disc_MPRs$pos[disc_MPRs$strand=="-"&!is.na(disc_MPRs$pos)],width=RL) #left pos
  #rMDp and rMDn contains all discordant and improper pairs on the same chr, doulbe counted

  #handling *hanging flags*
  hend_Aidx=all_MPRs$flag %in% hang_flags #anything that is truely hang
  hend_Aidx[is.na(hend_Aidx)]=FALSE       
  hend_MPRs=all_MPRs[hend_Aidx|disc_hend_Aidx|impp_hend_Aidx, #adding in disc and impp
                      -which(names(all_MPRs) %in% c("qwidth"))]

  rPn=length(all_MPRs[all_MPRs$strand=="+",]$pos[!is.na(all_MPRs[all_MPRs$strand=="+",]$pos)])
  rPc_idx=which((!is.na(all_MPRs$pos))&(!is.na(all_MPRs$qwidth))&(all_MPRs$strand=="+"))
  rPc=IRanges(start=all_MPRs$pos[rPc_idx],width=all_MPRs$qwidth[rPc_idx])
  rMn=length(all_MPRs[all_MPRs$strand=="-",]$pos[!is.na(all_MPRs[all_MPRs$strand=="-",]$pos)])
  rMc_idx=which((!is.na(all_MPRs$pos))&(!is.na(all_MPRs$qwidth))&(all_MPRs$strand=="-"))
  rMc=IRanges(start=all_MPRs$pos[rMc_idx],width=all_MPRs$qwidth[rMc_idx])

  if(verbose) cat(
    "---Info: [1] in concordant MPRs we always have FR order",
    "\n---Info: i.e. strand==+ <=> isize>0 <=> pos<mpos ",
    "\n---Info: in principle all these numbers shall be zero ",
    "\n---Info: flag bit 1st and 2nd in pair is for read in fastq only",
    "\n===Info: conc_lMPRs with isize>0 yet strand -:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$isize>0)&(all_MPRs$strand=='-')),
    "\n===Info: conc_lMPRs with mpos<pos yet strand +:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$strand=='+')&(all_MPRs$mpos<all_MPRs$pos)),
    "\n===Info: conc_lMPRs with isize>0 yet mpos<pos:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$isize>0)&(all_MPRs$mpos<all_MPRs$pos)),
    "\n===Info: conc_lMPRs with isize<0 yet strand +:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$isize<0)&(all_MPRs$strand=='+')),
    "\n===Info: conc_lMPRs with pos<mpos yet strand -:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$strand=='-')&(all_MPRs$pos<all_MPRs$mpos)), 
    "\n===Info: conc_lMPRs with isize<0 yet pos<mpos:", 
    sum((all_MPRs$flag %in% conc_flags)&(all_MPRs$isize<0)&(all_MPRs$pos<all_MPRs$mpos)),
    "\n---Info: meaningfulness of isize is currently ensured by mpos<pos",
    "\n---Info: assuming this trans-chromosome number should be strictly zero:",
    sum(conc_Fidx & all_MPRs$mrnm!=seq_name), 
    "\n---Info: if this number is not zero and above numbers are far from zero raise concerns!",
    "\n---Info: [2] additionally we want to see FR and RF order in improper pairs",
    "\n---Info: in principle these numbers shall be equal to each other",
    "\n===Info: improper Forward Read in proper order:",
    sum(impp_conc_Fidx), "Backword Read in proper order:", sum(impp_conc_Ridx),
    "\n===Info: improper Forward Read in inproper order:",
    sum(impp_disc_Fidx), "Backword Read in inproper order:", sum(impp_disc_Ridx),
    "\n---Info: both number not equal each other raise concerns!",
    "\n"
  )
  #cat("NA in conc_lMPRs?",any(is.na(c(conc_Fidx,impp_conc_Fidx))),"\n")
  conc_lMPRs=all_MPRs[conc_Fidx|impp_conc_Fidx,
                      -which(names(all_MPRs) %in% c("qwidth"))] 
  conc_rMPRs=all_MPRs[conc_Ridx|impp_conc_Ridx,
                      -which(names(all_MPRs) %in% c("qwidth"))]
  
  #here we deal with the difficulty with softclipped reads
  #TODO: if we want to reuse softclip reads for hanging peaks, we need
  #      1. find reads that were substantially clipped
  #      2. report its mate's position and derived flag (since it is concordant
  #      3. mate pair extend to the right append to rHr (isize>0)
  #      4. mate pair extend to the left append to rHl (isize<0)
  full_match=paste(RL,"M",sep="")
  schang_idn=grep(softclip_pattern,all_MPRs$cigar,value=F) #all cigars containing I/D/S/H, numeric
  al_left=get_al_omp(all_MPRs$cigar[schang_idn],alr_pattern)-get_al_omp(all_MPRs$cigar[schang_idn],alih_pattern)-get_al_omp(all_MPRs$cigar[schang_idn],alit_pattern)
  al_right=get_al_omp(all_MPRs$cigar[schang_idn],alr_pattern)-get_al_omp(all_MPRs$cigar[schang_idn],alih_pattern)-get_al_omp(all_MPRs$cigar[schang_idn],alit_pattern)
  pstran_idx=all_MPRs$strand=="+"
  nstran_idx=all_MPRs$strand=="-"
  schang_idx=rep(F,nrow(all_MPRs)); al_left_idx=rep(F,nrow(all_MPRs)); al_right_idx=rep(F,nrow(all_MPRs));
  schang_idx[schang_idn]=T; al_left_idx[schang_idn[al_left<=hang_clip]]=T; al_right_idx[schang_idn[al_right<=hang_clip]]=T
	#mate position of soft clipped - read
  rSCHp=IRanges(start=all_MPRs$mpos[nstran_idx&al_right_idx&!is.na(all_MPRs$mpos)&(all_MPRs$mrnm==seq_name)],width=1) 
  #mate position of soft clipped + read
  rSCHn=IRanges(start=all_MPRs$mpos[pstran_idx&al_left_idx&!is.na(all_MPRs$mpos)&(all_MPRs$mrnm==seq_name)],width=1) 
  N_shang=length(rSCHp)+length(rSCHn)
  
  #part_lidn=grep(softclip_pattern,conc_lMPRs$cigar,value=F) #all cigars containing I/D/S/H, numeric
  #conc_lMPRs[part_lidn[(al_left<=hang_clip)&(al_left<=RL)]]=TRUE # if aligned length<hang_clip, mark sclip as hang
  #pre_soft_lidx=part_lidn[(al_left>=prop_clip)&(al_left<RL)] # if aligned length>prop_clip, use isize from sclip
  #part_ridn=grep(complex_pattern,conc_rMPRs$cigar,value=F) #all cigars containing I/D/S/H, numeric
  #stopifnot(all(al_right>=0))
  #shang_ridx=rep(FALSE,nrow(conc_rMPRs))
  #shang_ridx[part_ridn[(al_right<=hang_clip)&(al_right<RL)]]=TRUE # if aligned length<hang_clip, mark sclip as hang
  #pre_soft_ridx=part_ridn[(al_right>=prop_clip)&(al_right<RL)] # if aligned length>prop_clip, use isize from sclip
  #hang_lMPRs=conc_lMPRs[shang_lidx,]; hang_rMPRs=conc_rMPRs[shang_ridx,]
  #N_shang=nrow(hang_rMPRs)+nrow(hang_lMPRs)

  #for now to reduce running time and allow for read name mess ups we don't use cigar code for rl 
  #isize of softclipped pair is extracted from the softclipped read itself
  #given TLEN is calculated by CIGAR, I doubt following is necessary
  #we will just use whatever is already in conc_lMPRs
  #soft_lidx=pre_soft_lidx
  #soft_lidx_impute=conc_lMPRs$isize[soft_lidx]
  #soft_lMPRs=conc_lMPRs[soft_lidx,]
  #soft_ridx=pre_soft_ridx
  #soft_ridx_impute=conc_rMPRs$isize[soft_ridx]
  #soft_rMPRs=conc_rMPRs[soft_ridx,]; 
  #soft_MPRs=soft_rMPRs #swith soft_rMPRs to be partner presented
  #soft_MPRs$pos=soft_rMPRs$mpos
  #soft_MPRs$mpos=soft_rMPRs$pos
  #soft_MPRs$isize=abs(soft_ridx_impute)+1 #now it is converted to _lMPRs and +1 to avoid zero
  #soft_MPRs=rbind(soft_MPRs,soft_lMPRs) #we only need to include soft_rMPRs into conc_lMPRs 
  #N_soft=nrow(soft_MPRs)
  #rSCISp=IRanges(start=soft_lMPRs$pos, width=abs(soft_lMPRs$isize)+1)
  #rSCISn=IRanges(start=soft_rMPRs$mpos, width=abs(soft_rMPRs$isize)+1)
  N_soft=NA
  rSCISp=IRanges()
  rSCISn=IRanges()
  #rm(conc_rMPRs); gc();

  #now we gathering the read information
  N_prop=nrow(conc_lMPRs);  #conc_lMPRs should contain
  mout_lidx=(conc_lMPRs$mpos-conc_lMPRs$pos>Delta)        #isize>delta, in lCd and lDx
  mover_lidx=(conc_lMPRs$mpos-conc_lMPRs$pos>=bigDel)     #isize>bigDel, in bigDel
  min_lidx=(conc_lMPRs$mpos-conc_lMPRs$pos<=maxInsert)    #isize<maxInsert, in lCd and bigD
  mout_ridx=(conc_rMPRs$pos-conc_rMPRs$mpos>Delta)
  mover_ridx=(conc_rMPRs$pos-conc_rMPRs$mpos>=bigDel)  #bigDel<isize<maxInsert double used
  min_ridx=(conc_rMPRs$pos-conc_rMPRs$mpos<=maxInsert)
  #cat("NA in conc_lMPRs$mpos?",any(is.na(conc_lMPRs$mpos)),"\n")
  #cat("NA in conc_lMPRs$pos?",any(is.na(conc_lMPRs$pos)),"\n")
  #cat("NA in all_MPRs$pos[conc_Fidx]?",any(is.na(all_MPRs$pos[conc_Fidx])),"\n")
  #cat("NA in all_MPRs$pos[impp_conc_Fidx]?",any(is.na(all_MPRs$pos[impp_conc_Fidx])),"\n")
  #cat("NA in all_MPRs$mpos[conc_Fidx]?",any(is.na(all_MPRs$mpos[conc_Fidx])),"\n")
  #cat("NA in all_MPRs$mpos[impp_conc_Fidx]?",any(is.na(all_MPRs$mpos[impp_conc_Fidx])),"\n")
  #cat("NA in mout_idx?",any(is.na(mout_idx)),"\n"); print(which(is.na(mout_idx)))
  #print(conc_lMPRs$pos[which(is.na(mout_idx))])
  #print(conc_lMPRs$mpos[which(is.na(mout_idx))])
  #cat("NA in mover_idx?",any(is.na(mover_idx)),"\n"); print(which(is.na(mover_idx))) 
  #cat("NA in min_idx?",any(is.na(min_idx)),"\n"); print(which(is.na(min_idx)))
  if(verbose) cat("---Info: found",sum(c(mout_lidx,mout_ridx)),"conc_MPRs out of delta=",Delta,"counted in lDl\n")
  if(verbose) cat("---Info: found",sum(c(min_lidx,min_ridx)),"conc_MPRs in delta=",Delta,"counted in lCd \n")
  if(verbose) cat("---Info: found",sum(c(mover_lidx,mover_ridx)),"conc_MPRs embed bigdel",bigDel,"counted in bigD\n")
  if(verbose) cat("---Info: found",sum(c(conc_lMPRs$mpos-conc_lMPRs$pos<=maxInsert & conc_lMPRs$mpos-conc_lMPRs$pos>=bigDel,conc_rMPRs$pos-conc_rMPRs$mpos<=maxInsert & conc_rMPRs$pos-conc_rMPRs$mpos>=bigDel)),"got double counted in lCd and bigD\n")
  conc_isize_lidx=!is.na(conc_lMPRs$isize) & conc_lMPRs$isize>=RL  #only use consistent isize
  conc_isize_ridx=!is.na(conc_rMPRs$isize) & conc_rMPRs$isize<=-RL #only use consistent isize
  #problematic isize goes to disc 
  if(sum(c(min_lidx,min_ridx))!=0)
    rPbi=IRanges(start=c(conc_lMPRs$pos[min_lidx&conc_isize_lidx],conc_rMPRs$mpos[min_ridx&conc_isize_ridx]),
                 width=c(conc_lMPRs$isize[min_lidx&conc_isize_lidx],-conc_rMPRs$isize[min_ridx&conc_isize_ridx]))
    #rPbi=IRanges(start=c(conc_lMPRs$pos[conc_isize_lidx],conc_rMPRs$mpos[conc_isize_ridx]),
    #             width=c(conc_lMPRs$isize[conc_isize_lidx],-conc_rMPRs$isize[conc_isize_ridx]))
  if(sum(c(mover_lidx,mover_ridx))!=0) 
    rPbo=IRanges(start=c(conc_lMPRs$pos[mover_lidx&conc_isize_lidx],conc_rMPRs$mpos[mover_ridx&conc_isize_ridx]),
                 width=c(conc_lMPRs$isize[mover_lidx&conc_isize_lidx],-conc_rMPRs$isize[mover_ridx&conc_isize_ridx]))
    rPbd=IRanges(start=pmin(disc_MPRs$pos,disc_MPRs$mpos,na.rm=T),
               	 end=pmax(disc_MPRs$pos,disc_MPRs$mpos,na.rm=T))

  #soft_mout_idx=which(soft_MPRs$mpos-soft_MPRs$pos>Delta)
  #soft_mover_idx=which(soft_MPRs$mpos-soft_MPRs$pos>=bigDel)
  #soft_min_idx=which(soft_MPRs$mpos-soft_MPRs$pos<=Delta)
  #isize based for softclipped
  #if(length(soft_min_idx)!=0)
  #  rSCPbi=IRanges(start=soft_MPRs$pos[soft_min_idx],width=soft_MPRs$isize[soft_min_idx])
  #if(length(soft_mover_idx)!=0) 
  #  rSCPbo=IRanges(start=soft_MPRs$pos[soft_mover_idx],width=soft_MPRs$isize[soft_mover_idx])

  #important information based on strand and order in hang pair
  #hend_MPRs contains not only truely hang but those likely hang from disc and impp
  #hang_flags=c(73,137,121,185,105,169,89,153) 
  #impp_flags=c(97,161,81,145)
  #conc_flags=c(99,163,83,147)
  #disc_flags=c(65,129,113,177)
  #rHl should be defined as read HANG to the left of break point: 
    #case i: read is first AND read is forward (73,97,99); 
    #case ii: read is second AND read is forward (137,161,163).
  #rHr should be defined as read HANG to the left of break point: 
    #case iii: read is first AND read is reverse (89,81,83); 
    #case iv: read is second AND read is reverse (153,145,147).
  #rHl will contribute to lDl; rHr will contribute to lDr

  #if(verbose) cat("---Info: hang_lMPRs:",nrow(hang_lMPRs),"hend_lMPRs",sum(hend_MPRs$strand=='+'),"\n")
  #if(verbose) cat("---Info: hang_rMPRs:",nrow(hang_rMPRs),"hend_rMPRs",sum(hend_MPRs$strand=='-'),"\n")
  #flags could be both way c(113,129,65,177)
  rMHp=IRanges(start=sort(
    c(hend_MPRs$pos[hend_MPRs$strand=="+"&!is.na(hend_MPRs$pos)],
      conc_lMPRs$pos[mout_lidx&conc_lMPRs$flag %in% c(97,99,161,163)],
      conc_rMPRs$pos[mout_ridx&conc_rMPRs$flag %in% c(97,99,161,163)]
		)),width=1)
  rMHn=IRanges(start=sort(
		c(hend_MPRs$pos[hend_MPRs$strand=="-"&!is.na(hend_MPRs$pos)],
      conc_lMPRs$pos[mout_lidx&conc_lMPRs$flag %in% c(81,83,145,147)],
      conc_rMPRs$pos[mout_ridx&conc_rMPRs$flag %in% c(81,83,145,147)]
		)),width=1) #left pos
  #rMHp and rMHn contains all discordant and improper pairs on the same chr, double counted
  if(verbose) print(table(hend_MPRs$flag))
  if(verbose) cat("---Info: hend_MPRs", nrow(hend_MPRs),"\n")
  rHl=IRanges(start=sort( 
    c(start(rMHp),
      start(rSCHp)
      )),width=1)
  if(verbose) cat("---Info: rHl:",length(rHl),"\n")
  rHr=IRanges(start=sort( 
    c(start(rMHn),
      start(rSCHn)
      )),width=1)
  if(verbose) cat("---Info: rHr:",length(rHr),"\n")

  #this part may need to be updated in the future
  #rSl=IRanges(start=sort(hang_lMPRs$pos),width=1)
  #rSr=IRanges(start=sort(hang_rMPRs$pos),width=1) #self of softhang
  rSl=IRanges()
  rSr=IRanges()

  merge_MPRs = if(verbose) conc_lMPRs else NULL #only keep for plotting if verbose
  if(verbose) cat("---Info: read set summary\n"); 
  if(verbose) cat("===Info: rPbi:\n"); if(verbose) print(summary(width(rPbi)))
  if(verbose) cat("===Info: total rPbi:",length(rPbi),"\n")
  if(verbose) cat("===Info: rPbo:\n"); if(verbose) print(summary(width(rPbo)))
  if(verbose) cat("===Info: total rPbo:",length(rPbo),"\n")
  if(verbose) cat("===Info: rPbd:\n"); if(verbose) print(summary(width(rPbd)))
  if(verbose) cat("===Info: total rPbd:",length(rPbd),"\n")
  if(verbose) cat("===Info: total MPRs (double entry)",N_all, pct(1), "(all)\n")
  if(verbose) cat("===Info: prop MPRs (double entry)",N_prop*2, pct(2*N_prop/N_all), "(all)\n") 
  if(verbose) cat("===Info: soft MPRs (double entry)",N_soft*2, pct(2*N_soft/N_all), "(all)\n")  
  if(verbose) cat("===Info: soft_hang MPRs (double entry)",N_shang*2, pct(2*N_shang/N_all), "(all)\n") 
  if(verbose) cat("===Info: hang MPRs (single entry)",nrow(hend_MPRs),pct(nrow(hend_MPRs)/N_all),"(all)\n") 
  if(verbose) cat("===Info: disc MPRs (double entry)",nrow(disc_MPRs),pct(nrow(disc_MPRs)/N_all), "(all)\n")
  if(verbose) cat("---Info: out get_MPRs, start=",start(ranges(seq_info)),"end=",end(ranges(seq_info)),"\n")
  return(list(rPbi=rPbi,rPbo=rPbo,rPbd=rPbd,rPc=rPc,rMc=rMc,rPn=rPn,rMn=rMn,rSr=rSr,rSl=rSl,rHr=rHr,rHl=rHl,
    rMHp=rMHp,rMHn=rMHn,rMDp=rMDp,rMDn=rMDn,rSCHp=rSCHp,rSCHn=rSCHn,
    rSCISp=rSCISp,rSCISn=rSCISn,rSDa=rSDa,rSIa=rSIa,rSCPbi=rSCPbi,rSCPbo=rSCPbo,
    hang_clip=hang_clip,prop_clip=prop_clip,
    RL=RL,isize=isize,isize_sd=isize_sd,Delta=Delta,
    smallDel=smallDel,smallIns=smallIns,bigDel=bigDel,maxInsert=maxInsert,anchor=anchor))
} # get MPRs from bam_file

#' @title scan_joint
#' 
#' @description 
#' scan_joint will return SWAN score tracks
#' 
#' @param srange
#' @param width
#' @param lw_width
#' @param stepsize
#' @param block_size
#' @param Delta
#' @param RL
#' @param maxInsert
#' @param p_input
#' @param q_input
#' @param scan_cols
#' @param isize
#' @param isize_sdR
#' @param isize_sdL
#' @param coverage
#' @param mixing_rate
#' @param rPbi
#' @param rPbo
#' @param rPbd
#' @param rPc
#' @param rMc
#' @param rPn
#' @param rMn
#' @param rSr
#' @param rSl
#' @param rHr
#' @param rHl
#' 
#' @note \code{scan_joint} 
#' \cr
#' If called, \code{scan_joint} will return SWAN score tracks
#' 
#' @examples
#' scan_joint<-function(srange, width, lw_width, stepsize, block_size,
#' Delta, RL, maxInsert, p_input, q_input, scan_cols, 
#' isize, isize_sdR, isize_sdL, coverage, mixing_rate,  
#' rPbi, rPbo, rPbd, rPc, rMc, rPn, rMn, rSr, rSl, rHr, rHl,
#' debug, verbose
#' )
#'              
#' @seealso \code{\link{get_MPRs}} and \code{\link{scan_bam}}
scan_joint<-function(srange, width, lw_width, stepsize, block_size,
                     Delta, RL, maxInsert, p_input, q_input, scan_cols, 
                     isize, isize_sdR, isize_sdL, coverage, mixing_rate,  
                     rPbi, rPbo, rPbd, rPc, rMc, rPn, rMn, rSr, rSl, rHr, rHl,
                     debug, verbose
                     ){ 
  if(verbose) cat("---Info: into scan_joint ...\n")
  if(verbose) cat("---Info: options: ... \n")
  if(verbose) cat("width","lw_width","stepsize","block_size","Delta","RL","maxInsert","p_input","q_input","scan_cols","isize","isize_sdR","isize_sdL","coverage","mixing_rate","maxInsert","\n")
  if(verbose) cat(paste(width,lw_width,stepsize,block_size,Delta,RL,maxInsert,p_input,q_input,paste(scan_cols,collapse="-"),isize,isize_sdR,isize_sdL,coverage,mixing_rate,maxInsert,sep="\t"),"\n")
  wplus=max(width,0)
  spanW=wplus #effective window size, now neglect RL
  #w_start=floor((Delta+RL)/stepsize)+1   # leading windows discarded
  #w_end=length(srange)-(floor((wplus+Delta+RL)/stepsize)+1) # ending windows discarded
  if(verbose) cat("===Info: width=",width,";wplus=",wplus,";spanW=",spanW,"\n")
  winW=IRanges(start=srange,end=srange+wplus-1) #only use positive windows to save time 
  winW2=IRanges(start=srange,end=srange+lw_width-1)
  #NOTE: ideally we shall separate winW for w<0 and w>0 respectively
  winDl=IRanges(start=pmax(srange-Delta,1),width=pmin(srange,Delta)) #shrink Dl win to srange
  winDr=IRanges(start=srange+wplus,width=Delta)
  winZ=IRanges(start=srange,width=1)
  n_wins = length(winW)
  if(verbose) cat("===Info: wins: start",start(winW)[1],";end",end(winW)[n_wins],"bp;total",n_wins,"\n")
  
  #estimate p, zero rPn and rMn is handled
  q_left = length(rSl)/rPn # p estimated for + anchored read
  q_right = length(rSr)/rMn # p estimated for - anchored read
  p_left = length(rHl)/(rPn+rMn)/2.0 # p estimated for + anchored read
  p_right = length(rHr)/(rMn+rPn)/2.0 # p estimated for - anchored read
  if( p_input == -1 & p_left<p_fail & p_left>0 & p_right<p_fail & p_right>0 ) {
    #p = (length(rSl)+length(rSr))/(length(rM)+length(rP))
    if(verbose) cat("===Info: using estimated p+ =", p_left, "p- =", p_right, "\n")
  } else if (p_input>0 & p_input<p_fail) {
    p_left = p_input; p_right = p_input
    if(verbose) cat("===Info: using inputed p+ and p- =", p_input, "\n")
  } else {
    if(verbose) cat("===Info: using shaky estimated p+ =", p_left, "p- =", p_right, "\n")
    p_left = p_fail; p_right = p_fail
    if(verbose) cat("---Warn: estimated hanging probability p too high, forcing p_fail,",p_fail,"
        considerring reducing hang_clip and rerun the scan \n")
  }
  if( q_input == -1 & q_left<q_fail & q_left>0 & q_right<q_fail & q_right>0 ) {
    #p = (length(rSl)+length(rSr))/(length(rM)+length(rP))
    if(verbose) cat("===Info: using estimated q+ =", q_left, "q- =", q_right, "\n")
  } else if (q_input>0 & q_input<q_fail) {
    q_left = q_input; q_right = q_input
    if(verbose) cat("===Info: using inputed q+ and q- =", q_input, "\n")
  } else {
    if(verbose) cat("===Info: using shaky estimated q+ =", q_left, "q- =", q_right, "\n")
    q_left = q_fail; q_right = q_fail
    if(verbose) cat("---Warn: estimated softclipping probability q too high, forcing q_fail,",q_fail,"
        considerring reducing hang_clip and rerun the scan \n")
  }
  
  #estimate lambda
  #if(verbose) cat("===Info: coverage\n")
  #if(verbose) print(coverage)
  #cat("srange[1]=", srange[1], "end(winW[n_wins])=", end(winW[n_wins]), 
  #    "min(start(rPbi))=", min(start(rPbi)), "max(end(rPbi))=", max(end(rPbi)), "\n")
  coverage_rle=IRanges::coverage(c(rPc,rMc))
  if(end(winW[n_wins])-length(coverage_rle)>0) # coverage: 0000 t_start ---- t_end 0000, append to length
    coverage=Rle(c(IRanges::as.vector(coverage_rle),rep(0,end(winW[n_wins])-length(coverage_rle))))
  cvg_track=IRanges::as.vector(coverage_rle)[start(winW)]
  #may slightly over estimate coverage due to random zeros but avoids difficulty handling gaps
  if(verbose) cat("===Info: reliable range=", max(min(start(rPbi)),srange[1]), 
                  min(max(end(rPbi)),end(winW[n_wins])), "\n")
  reliable_start=max(min(start(rPbi)),srange[1]); reliable_end=min(max(end(rPbi)),end(winW[n_wins]))
  if(reliable_end-reliable_start>0) {
    reliable_range=IRanges(start=reliable_start, end=reliable_end)
  } else {
    #no reliable range, return nothing
    if(verbose) cat("==Warn: no reliable range in this trunk, skipping...\n")
    scan_result=data.frame(matrix(rep(NA,n_wins*length(scan_cols)),ncol=length(scan_cols)))
    scan_result[,1]=srange
    scan_par=list(p_left=NA,p_right=NA,q_left=NA,q_right=NA,coverage=NA,r_start=NA,r_end=NA,
                  is=NA,issdR=NA,issdL=NA,rl=NA,maxInsert=maxInsert)
    scan_bigd=data.frame(); scan_disc=data.frame()
    return(list(scan_result=scan_result,scan_par=scan_par,scan_bigd=scan_bigd,scan_disc=scan_disc))
  }
  #coverage_mean=IRanges::aggregate(coverage,start=as.integer(start(reliable_range)),width=as.integer(width(reliable_range)),FUN=IRanges::mean) 
  #or we take mean of non-zero entries ?
  coverage_mean=sum(coverage_rle/sum(coverage_rle>0)) 
  lambda_mean=coverage_mean/RL
  if(verbose) cat("===Info: estimated average coverage=",coverage_mean,"lambda=",lambda_mean,"\n")
  if(verbose) cat("===Info: if the program broke below due to sporadic extreme high coverage, please rerun with extra memory\n")

  #estimate F and f
  maxInsert=max(width(rPbi),maxInsert) #keep this if observed maxInsert is bigger than set
  
  Fx = FI(0:(Delta+1000),isize,isize_sdR)            #Fx=FI(s-x), s-x>=0, first is 0
  fy_cap = fy_cap_fail                               #this is if(fy>is_cut) use is_cut for fy
  fyR = fIw(0:(maxInsert+1),width,isize,isize_sdR)   #fIw for w=width and lCd
  fyR = ifelse(fyR>fy_cap,fy_cap,fyR)
  fyL = fIw(0:(maxInsert+1),-width,isize,isize_sdL)  #fIw for w=-width and lCi
  fyL = ifelse(fyL>fy_cap,fy_cap,fyL)

  #blockfied IRanges to reduce memory requirement; within block we maybe able to speedup later
  lW_term = rep(0,n_wins); lCd_term=rep(0,n_wins); lCi_term=rep(0,n_wins); lambda_term=rep(0,n_wins)
  lDl_term=rep(0,n_wins); lDr_term=rep(0,n_wins); lSl_term=rep(0,n_wins); lSr_term=rep(0,n_wins)
  cCd_term=rep(0,n_wins); cCi_term=rep(0,n_wins); 
  cDl_term=rep(0,n_wins); cDr_term=rep(0,n_wins);

  block_num=floor((n_wins-1)/block_size)+1;  #total number of blocks we got, always 1 if n_wins==block_size

  for(i in seq_len(block_num)){
    #gtk_now=proc.time() #ptm=proc.time(); pmm=gc();
    block_w_start=(i-1)*block_size+1; block_w_end=min(block_size*i, n_wins); 
    #if(verbose) cat("===Info: block executions",i,"-th of",block_num,",block_size=",block_w_end-block_w_start+1,"...\n")
    #cat("===Info: within-trunk block executions",i,"-th of",
    #      block_num,",block_size=",nsf(block_w_end-block_w_start+1),"...\n")
    #cat("==block_w_start=", block_w_start, "block_w_end=", block_w_end, "\n")
    block_winW=winW[block_w_start:block_w_end]; block_winW_start=start(block_winW); block_winW_end=end(block_winW);
    block_winW2=winW2[block_w_start:block_w_end];
    block_winZ=winZ[block_w_start:block_w_end];
    block_winDl=winDl[block_w_start:block_w_end]; block_winDr=winDr[block_w_start:block_w_end];
    block_start=start(winW[block_w_start]); block_end=end(winW[block_w_end]);
    #block_lambda_start=block_start
    #block_lambda_width=block_end-block_lambda_start
    #block_lambda_start=max(block_start, start(reliable_range))
    #block_lambda_width=min(block_end, end(reliable_range))-block_lambda_start
    #block_winWlam2=IRanges::aggregate(lambda, start=block_lambda_start, width=block_lambda_width, FUN=IRanges::mean)*wplus
    block_winWlam=lambda_mean*lw_width 
    block_nP_winW2=IRanges::as.matrix(countOverlaps(block_winW2,IRanges(start=start(rPc),width=1))) #positive ends in W
    block_nM_winW2=IRanges::as.matrix(countOverlaps(block_winW2,IRanges(start=start(rMc),width=1))) #negative ends in W
    #fix lW term, r increase sig decrease; w increase sig decrease; use small r and w!!! 
    block_lW_term = block_winWlam*mixing_rate + (block_nP_winW2+block_nM_winW2)*log(1-mixing_rate)
    #block_winW_rPb=IRanges::as.matrix(findOverlaps_within(block_winW,rPbi)) #crossing pair, doesn't mean corssing
    #block_winZ_rPb=IRanges::as.matrix(findOverlaps_within(block_winZ,rPbi)) #crossing pair, doesn't mean corssing
    block_winW_rPb=IRanges::as.matrix(IRanges::findOverlaps(block_winW,rPbi,type="within")) #find crossing pair
    block_winZ_rPb=IRanges::as.matrix(IRanges::findOverlaps(block_winZ,rPbi,type="within")) #find crossing pair
    #rPbl=IRanges(start=start(rPbi),width=1); 
    #rPbr=IRanges(start=pmax(1,end(rPbi)-RL),width=1) # left/right end points of rPbi
    #block_rPbl_winDl=IRanges::as.matrix(findOverlaps_within(rPbl,block_winDl))
    #block_rPbr_winDr=IRanges::as.matrix(findOverlaps_within(rPbr,block_winDr))
    #block_rPbr_winZDr=IRanges::as.matrix(findOverlaps_within(rPbr,block_winZDr))
    #block_winC_rPb=rbind(block_rPbl_winDl,block_rPbr_winDr)
    #block_winC_rPb=block_winC_rPb[duplicated(block_winC_rPb),]
    #block_winC_rPb=block_winC_rPb[,c(2,1)] #switch column to the right order
    #block_winZ_rPb=rbind(block_rPbl_winDl,block_rPbr_winZDr)
    #block_winZ_rPb=block_winZ_rPb[duplicated(block_winZ_rPb),] 
    #block_winZ_rPb=block_winZ_rPb[,c(2,1)] #switch column to the right order
    #stopifnot(nrow(unique(block_winC_rPb))==nrow(block_winC_rPb)) #That says everything duplicated is only twice
    #stopifnot(nrow(unique(block_winZ_rPb))==nrow(block_winZ_rPb)) #That says everything duplicated is only twice
    #block_rSr_winDl=IRanges::as.matrix(findOverlaps_within(rSr,block_winDl)) #sc-ends in D+, q <-> sbj
    #block_rSl_winDr=IRanges::as.matrix(findOverlaps_within(rSl,block_winDr)) #sc-ends in D-
    #block_rHr_winDr=IRanges::as.matrix(findOverlaps_within(rHr,block_winDr)) #mh-ends in D+
    #block_rHl_winDl=IRanges::as.matrix(findOverlaps_within(rHl,block_winDl)) #mh-ends in D-
    block_rSr_winDl=IRanges::as.matrix(IRanges::findOverlaps(rSr,block_winDl,type="within")) #sc-ends in D+, q <-> sbj
    block_rSl_winDr=IRanges::as.matrix(IRanges::findOverlaps(rSl,block_winDr,type="within")) #sc-ends in D-
    block_rHr_winDr=IRanges::as.matrix(IRanges::findOverlaps(rHr,block_winDr,type="within")) #mh-ends in D+
    block_rHl_winDl=IRanges::as.matrix(IRanges::findOverlaps(rHl,block_winDl,type="within")) #mh-ends in D-
    #need tryCatch here; give all zeros if error
    block_Cd_term=tryCatch({
        lC_omp(block_winW_rPb,width(rPbi),fyR,block_w_end-block_w_start+1,mixing_rate)
      }, error=function(e) { 
	if(verbose) { print(e); cat("==Warn: salvaging the scan by embedding zeros lCd and lCi terms this trunk only\n") }
        list(lC_term=rep(0,n_wins),cC_term=rep(0,n_wins))
      })
    block_Ci_term=tryCatch({
        lC_omp(block_winZ_rPb,width(rPbi),fyL,block_w_end-block_w_start+1,mixing_rate)
      }, error=function(e) { 
	if(verbose) { print(e); cat("==Warn: salvaging the scan by embedding zeros lCd and lCi terms this trunk only\n") }
        list(lC_term=rep(0,n_wins),cC_term=rep(0,n_wins))
      })
    block_Dl_term=lD_omp(block_rHl_winDl,start(block_winW),start(rHl),Fx,block_w_end-block_w_start+1,mixing_rate,p_left,left=TRUE)
    block_Dr_term=lD_omp(block_rHr_winDr,end(block_winW),start(rHr),Fx,block_w_end-block_w_start+1,mixing_rate,p_right,left=FALSE) 
    Fs=rep(1,length(Fx))
    block_Sr_term=lD_omp(block_rSr_winDl,start(block_winW),start(rSr),Fs,block_w_end-block_w_start+1,mixing_rate,q_right,left=TRUE)
    block_Sl_term=lD_omp(block_rSl_winDr,end(block_winW),start(rSl),Fs,block_w_end-block_w_start+1,mixing_rate,q_left,left=FALSE)
    lW_term[block_w_start:block_w_end]=block_lW_term
    lCd_term[block_w_start:block_w_end]=as.vector(block_Cd_term[["lC_term"]])
    cCd_term[block_w_start:block_w_end]=as.vector(block_Cd_term[["cC_term"]])
    lCi_term[block_w_start:block_w_end]=as.vector(block_Ci_term[["lC_term"]])
    cCi_term[block_w_start:block_w_end]=as.vector(block_Ci_term[["cC_term"]])
    lDl_term[block_w_start:block_w_end]=as.vector(block_Dl_term[["lD_term"]])
    cDl_term[block_w_start:block_w_end]=as.vector(block_Dl_term[["cD_term"]])
    lDr_term[block_w_start:block_w_end]=as.vector(block_Dr_term[["lD_term"]])
    cDr_term[block_w_start:block_w_end]=as.vector(block_Dr_term[["cD_term"]])
    lSl_term[block_w_start:block_w_end]=as.vector(block_Sl_term[["lD_term"]])
    lSr_term[block_w_start:block_w_end]=as.vector(block_Sr_term[["lD_term"]])
  }

  if(length(rPbo)!=0){
    scan_bigd=merge_rclu(lstart=start(rPbo),lend=start(rPbo)+RL,rstart=end(rPbo)-RL,rend=end(rPbo),support=rep(1,length(rPbo)))
  } else {
    scan_bigd=data.frame()
  }
  
  if(length(rPbd)!=0){
    scan_disc=merge_rclu(lstart=start(rPbd),lend=start(rPbd)+RL,rstart=end(rPbd)-RL,rend=end(rPbd),support=rep(1,length(rPbd)))
  } else {
    scan_disc=data.frame()
  }

  #cat("br2\n")
  #cat(class(cCd_term),length(cCd_term),"\n")
  #cat(class(lCd_term),length(lCd_term),"\n")
  scan_stat=as.data.frame(list(start=start(winW),lW=lW_term,lCd=lCd_term,lCi=lCi_term,
    lDr=lDr_term,lDl=lDl_term,lSl=lSl_term,lSr=lSr_term,cvg=cvg_track,
    cCd=cCd_term,cCi=cCi_term,cDl=cDl_term,cDr=cDr_term
    )) #NOTE: can use list to save memory space
  #cat("br3\n")
  scan_par=data.frame(list(p_left=p_left,p_right=p_right,q_left=q_left,q_right=q_right,
    is=isize,issdR=isize_sdR,issdL=isize_sdL,rl=RL,maxInsert=maxInsert,
    coverage=coverage_mean,r_start=start(reliable_range),r_end=end(reliable_range)))
  #if(verbose) 
  if(verbose) cat("---Info: trunk coverage summary:\n")
  if(verbose) print(summary(cvg_track)) 
  if(verbose) cat("---Info: out scan_join \n")
  #print(scan_stat$lCd_term[1:10])
  #print(scan_stat$cCd_term[1:10])
  return(list(scan_result=scan_stat,scan_par=scan_par,scan_bigd=scan_bigd,scan_disc=scan_disc))
}

#cat("loaded line 4000... \n")

apply_scan = function(seq_name, rg_files, rg_prefix, sp_prefix,
               what,isize_global,
               coverage_global, 
               start, end, RL_global, gap, trunk_size,
               width, lw_width, mixing_rate, stepsize, block_size,
               Delta_global, bigDel_global, smallDel_global, smallIns_global, maxInsert_global, 
               p_global, q_global,
               hang_clip_global, prop_clip_global,
               fast, usestat, debug, verbose, gtk){
  #rg_files=rg_files;rg_prefix=rg_prefix;sp_prefix=sp_prefix;what=what;isize_global=isize_global;coverage_global=coverage_global;seq_name=seq_name;start=start;end=end;RL_global=RL_global;gap=gap;trunk_size=trunk_size;width=width;mixing_rate=mixing_rate;stepsize=stepsize;block_size=block_size;Delta_global=Delta_global;bigDel_global=bigDel_global;p_global=p_global;q_global=q_global;hang_clip_global=hang_clip_global;prop_clip_global=prop_clip_global;usestat=TRUE;fast=fastSave;debug=debug;verbose=verbose;gtk=gtk
  #other_opt needs to be the same for all (unless we need otherwise, will be more complicated)
  seq_info=GRanges(seq_name, IRanges(start, end))
  scan_cols=scan_cols_preset; srange=seq(start,end,stepsize); out_cols=out_cols_preset
  list[Delta_global]=split_args(Delta_global,rg_files,1)
  if(verbose) { cat("-Info: Delta_global\n"); print(Delta_global) }
  list[RL_global]=split_args(RL_global,rg_files,1)
  if(verbose) { cat("-Info: RL_global\n"); print(RL_global) }
  list[isize_global]=split_args(isize_global,rg_files,1)
  if(verbose) { cat("-Info: isize_global\n"); print(isize_global) }
  list[coverage_global]=split_args(coverage_global,rg_files,1)
  if(verbose) { cat("-Info: coverage_global\n"); print(coverage_global) }
  list[p_global]=split_args(p_global,rg_files,1)
  if(verbose) { cat("-Info: p_global\n"); print(p_global) }
  list[q_global]=split_args(q_global,rg_files,1)
  if(verbose) { cat("-Info: q_global\n"); print(q_global) }
  list[hang_clip_global]=split_args(hang_clip_global,rg_files,1)
  if(verbose) { cat("-Info: hang_clip_global\n"); print(hang_clip_global) }
  list[prop_clip_global]=split_args(prop_clip_global,rg_files,1)
  if(verbose) { cat("-Info: prop_clip_global\n"); print(prop_clip_global) }
  list[bigDel_global]=split_args(bigDel_global,rg_files,1)
  if(verbose) { cat("-Info: bigDel_global\n"); print(bigDel_global) }
  list[smallDel_global]=set_args(smallDel_global,rg_files,1)
  if(verbose) { cat("-Info: smallDel_global\n"); print(smallDel_global) }
  list[smallIns_global]=set_args(smallIns_global,rg_files,1)
  if(verbose) { cat("-Info: smallIns_global\n"); print(smallIns_global) }
  list[maxInsert_global]=set_args(maxInsert_global,rg_files,1)
  if(verbose) { cat("-Info: maxInsert_global\n"); print(maxInsert_global) }
  n_sp=length(rg_files)
  
  for(si in seq_len(n_sp)){ #si=1
    #we need to know start, end, stepsize and cnames to construct a zero dataframe
    sp_scores=data.frame(matrix(rep(0,length(srange)*length(scan_cols)),ncol=length(scan_cols)))
    colnames(sp_scores)=scan_cols; sp_scores[["start"]]=srange
    sp_bigd=list(); sp_disc=list(); sp_anch=list(); sp_out=sp_prefix[si]; sp_par=NULL; n_sp=length(rg_files)
    n_rg = length(rg_files[[si]])
    for(gi in seq_len(n_rg)){ #gi=1
      #each read group shall have its own isize estimates
      rg_scores=data.frame(); rg_par=list(); rg_bigd=data.frame(); rg_out=rg_prefix[[si]][gi];
      cat("=Info apply_scan: group",gi,"of",n_rg,"groups within sample",si,"of",n_sp,"samples\n")
      cat("=Info apply_scan: rg_file", rg_files[[si]][gi], "\n")
      cat("=Info apply_scan: rg_prefix", rg_prefix[[si]][gi], "\n")
      cat("=Info apply_scan: rg_out", rg_out, "\n")
      bam_file=rg_files[[si]][gi]; stat_file=paste(bam_file,"stat",sep=".");  bam_stat=NULL
      if(usestat) { bam_stat=if(file.exists(stat_file)) safe_read(stat_file,header=TRUE) }
      list[isize_rg,isize_sdR_rg,isize_sdL_rg,good_lCi_rg]=
        set_isize(isize_global[[si]][gi],rg_files[[si]][gi],seq_info) 
      if(verbose) cat("=Info isize:",isize_rg,isize_sdR_rg,isize_sdL_rg,"\n")
      coverage_rg=set_coverage(coverage_global[[si]][[gi]])
      RL_rg=set_RL(RL_global[[si]][[gi]])
      Delta_rg=set_Delta(Delta_global[[si]][[gi]])
      bigDel_rg=set_bigDel(bigDel_global[[si]][[gi]])
      smallDel_rg=set_smallDel(smallDel_global[[si]][[gi]])
      smallIns_rg=set_smallIns(smallIns_global[[si]][[gi]])
      maxInsert_rg=set_maxInsert(maxInsert_global[[si]][[gi]])
      p_rg=set_p(p_global[[si]][[gi]])
      q_rg=set_q(q_global[[si]][[gi]])
      hang_clip_rg=set_hang_clip(hang_clip_global[[si]][[gi]])
      prop_clip_rg=set_prop_clip(prop_clip_global[[si]][[gi]])
      scan_par=list(bam_file=bam_file,what=what,isize_global=isize_rg,isize_sdR_global=isize_sdR_rg,
              isize_sdL_global=isize_sdL_rg,coverage_global=coverage_rg, seq_name=seq_name,  
              start=start, end=end, RL_global=RL_rg, gap=gap, trunk_size=trunk_size,
              width=width, lw_width=lw_width, mixing_rate=mixing_rate, stepsize=stepsize, block_size=block_size,
              Delta_global=Delta_rg, bigDel_global=bigDel_rg, smallDel_global=smallDel_rg, 
              smallIns_global=smallIns_rg, maxInsert_global=maxInsert_rg, p_global=p_rg, q_global=q_rg,
              hang_clip_global=hang_clip_rg, prop_clip_global=prop_clip_rg,
              debug=debug, rg_out=rg_out, verbose=verbose, fast=fastSave, gtk=gtk)
      if(debug) cat("=Info :", rg_files[[si]][gi],"start=",nsf(start),"end=",nsf(end), "gtk=",get_gtk(gtk),"s,gmk=",get_gmk(Sys.getpid()),"kb\n") 
      list[rg_scores,rg_par,rg_bigd,rg_disc,rg_anch]=do.call(scan_bam,scan_par)
      rg_par$lCd=TRUE; rg_par$lW=TRUE; rg_par$lCi=TRUE; rg_par$lDl=TRUE; 
      rg_par$lDr=TRUE; rg_par$lSl=TRUE; rg_par$lSr=TRUE
      cat("=Info apply_scan: merging rg scores to sp scores,",taggie(gtk),"\n")
      for(ci in seq(2,dim(sp_scores)[2])){
        score=scan_cols[ci]
        if(!is.null(bam_stat[[score]])) { if(!bam_stat[[score]]) { rg_par[[score]]=FALSE; next } }
        sp_scores[,ci]=sp_scores[,ci]+rg_scores[,ci]
      }
      sp_bigd=rbind(sp_bigd,rg_bigd)
      sp_disc=rbind(sp_disc,rg_disc)
      sp_anch=rbind(sp_anch,rg_anch)
      sp_par=rbind(sp_par,rg_par)
      if(debug) {
        cat("=Info apply_scan: saving rg scores,",taggie(gtk),"\n")
        bigd=paste(rg_out,"bigd","txt",sep=".")
        write.table(rg_bigd,bigd,quote=F,row.names=F,sep="\t")
        disc=paste(rg_out,"disc","txt",sep=".")
        write.table(rg_disc,disc,quote=F,row.names=F,sep="\t")
        anch=paste(rg_out,"anch","txt",sep=".")
        write.table(rg_anch,anch,quote=F,row.names=F,sep="\t")
        gzf=paste(rg_out,"swan","txt","gz",sep=".") #remove .bam from output
        comment=call_meta(meta_keys,meta_values,ref_info=ref_file,cmd_info=paste(commandArgs(),collapse=" "))
        cat("rg_names",names(rg_scores),"\n")
        write_com(round(rg_scores[,out_cols],digits=4),comment,gzf,gzip=T,quote=F,row.names=F,sep="\t")
        parf=paste(rg_out,"swan","par","txt",sep=".")
        write.table(format(rg_par,digits=2,scientific=F),parf,quote=F,row.names=F,sep="\t")
      }
    }  
    cat("=Info apply_scan: saving sp parameter and scores,",taggie(gtk),"\n")
    parf=paste(sp_out,"swan","par","txt",sep=".")
    write.table(format(sp_par,digits=2,scientific=F),parf,quote=F,row.names=F,sep="\t")
    gzf=paste(sp_out,"swan","txt","gz",sep=".")
    comment=call_meta(meta_keys,meta_values,ref_info=ref_file,cmd_info=paste(commandArgs(),collapse=" "))
    cat("sp_names",names(sp_scores),"\n")
    tmp=sp_scores[,out_cols]
    write_com(round(sp_scores[,out_cols],digits=6),comment,gzf,gzip=T,quote=F,row.names=F,sep="\t")
    sp_bigd=merge_rclu(sp_bigd$lstart,sp_bigd$lend,sp_bigd$rstart,sp_bigd$rend,sp_bigd$support)
    sp_disc=merge_rclu(sp_disc$lstart,sp_disc$lend,sp_disc$rstart,sp_disc$rend,sp_disc$support)
    cat("=Info apply_scan: saving sp_bigd (bigd) size=",nrow(sp_bigd),"\n")
    cat("=Info apply_scan: saving sp_disc (disc) size=",nrow(sp_disc),"\n")
    cat("=Info apply_scan: saving sp_anch (anch) size=",nrow(sp_anch),"\n")
    bigd=paste(sp_out,"bigd","txt",sep=".")
    write.table(sp_bigd,bigd,quote=F,col.names=F,row.names=F,sep="\t")
    disc=paste(sp_out,"disc","txt",sep=".")
    write.table(sp_disc,disc,quote=F,col.names=F,row.names=F,sep="\t")
    anch=paste(sp_out,"anch","txt",sep=".")
    write.table(sp_anch,anch,quote=F,col.names=F,row.names=F,sep="\t")
    if(debug){
      save(list=c("scan_par"),file=paste(sp_out,".RData",sep=""))
    }
  } #now we are done, move to next sample file
  return(T)
}

count_bam=function(bam_file,chr=""){
  res=system(sprintf("samtools view -c -F 4 %s %s",bam_file,chr),intern=T) #FIXME: dirty quick hack here?
  as.numeric(res)
}

apply_seqcbs = function(rg_files,sp_prefix,gap,seq_name,start,end,width,stepsize,ref_head,ref_tail,
                        debug,verbose,gtk){
  seq_info=GRanges(seq_name, IRanges(start, end));
  win_pos=seq(start,end,stepsize)
  seqcbs_what=c("pos","mapq","qwidth"); MAPQ.THRESH=0; #seqcbs_read0=IRanges()
  seqcbs_pos0=c(); seqcbs_qwidth0=c(); n_sp=length(rg_files); seqcbs_par=data.frame()
  nr=rep(0,n_sp)
  for(si in seq_len(length(rg_files))){ #si=1 is always the baseline, si>1 always compare to 1
    sp_out=sp_prefix[si]; seqcbs_pos1=c(); seqcbs_qwidth1=c(); n_rg=length(rg_files[[si]])
    for(gi in seq_len(length(rg_files[[si]]))){
      cat("=Info seqcbs: group",gi,"of",n_rg,"groups within sample",si,"of",n_sp,"samples\n")
      cat("=Info seqcbs: rg_file",rg_files[[si]][[gi]],taggie(gtk),"\n")
      bam_file=rg_files[[si]][[gi]]
      bam_full=BamFile(bam_file); bam_sum=count_bam(bam_file); nr[si]=nr[si]+bam_sum
      rg_reads=allFunction(seq_info,bam_file,what=seqcbs_what,index=bam_file)
      sel=which(rg_reads[["mapq"]]>MAPQ.THRESH&!is.null(rg_reads[["pos"]]))
      if(si==1) { #this 1st must be control, second or later is case
        seqcbs_pos0=c(seqcbs_pos0,rg_reads[["pos"]][sel])
        seqcbs_qwidth0=c(seqcbs_qwidth0,rg_reads[["qwidth"]][sel])
      } else {
        seqcbs_pos1=c(seqcbs_pos1,rg_reads[["pos"]][sel])
        seqcbs_qwidth1=c(seqcbs_qwidth1,rg_reads[["qwidth"]][sel])
      }
    }
    seqcbs_bic=rep(0,n_wins); seqcbs_stat=rep(0,n_wins); seqcbs_res=data.frame()
    if(si>1){
      if(verbose) cat("--Info: into ScanCBS",taggie(gtk),"\n")
      seqcbs_res=as.data.frame(seqCBS::ScanCBS(seqcbs_pos1,seqcbs_pos0,takeN=10,maxNCut=100,verbose=FALSE)$statHat)
      if(verbose) cat("==Info: out ScanCBS",taggie(gtk),"\n")
      in_gap=inGap(seqcbs_res[,1],seq_name,gap)&inGap(seqcbs_res[,2],seq_name,gap)
      if(verbose) cat("==Info: remove",sum(in_gap),"regions in gap\n")
      seqcbs_res=seqcbs_res[!in_gap,]
      #if(!is.null(seqcbs_res)){
        for(ri in seq_len(nrow(seqcbs_res))){
          ri_w_start=floor((seqcbs_res[ri,1]-start)/stepsize)+1 #cptL=1
          ri_w_end=floor((seqcbs_res[ri,2]-start)/stepsize)+1 #cptR=2
          seqcbs_bic[ri_w_start:ri_w_end]=seqcbs_res[ri,8] #BIC=8
          seqcbs_stat[ri_w_start:ri_w_end]=seqcbs_res[ri,5] #stat=5
        }
      #}
      sp_scores=as.data.frame(list(lCBS=seqcbs_stat))
      sp_scores=set_colClass(sp_scores,rep("numeric",1))
      gzf=paste(sp_out,"seqcbs","txt","gz",sep=".")
      comment=call_meta(meta_keys,meta_values,ref_info=ref_file,cmd_info=paste(commandArgs(),collapse=" "))
      write_com(format(sp_scores,digits=2,scientific=F),comment,gzf,gzip=T,quote=F,row.names=F,sep="\t")
      tabf=paste(sp_out,"seqcbs","txt",sep=".")
      write.table(seqcbs_res,tabf,row.names=F,quote=F,sep="\t")
      parf=paste(sp_out,"seqcbs","par","txt",sep=".")
      sp_par=as.data.frame(list(chr=seq_name,start=start,end=end,w=width,stepsize=stepsize,cvgrxy=nr[si]/nr[1]))
      write.table(format(sp_par,digits=2,scientific=F),parf,quote=F,row.names=F,sep="\t")
      if(debug) save(list=c("seqcbs_pos0","seqcbs_qwidth0","seqcbs_pos1","seqcbs_qwidth1"),file=paste(sp_out,"seqcbs","RData",sep="."))
    }
    #if(si==1) seqcbs_coverage=countOverlaps(IRanges(start=win_pos,width=width),
    #                                        IRanges(start=seqcbs_pos0,width=seqcbs_qwidth0))
    #if(si>1) seqcbs_coverage=countOverlaps(IRanges(start=win_pos,width=width),
    #                                        IRanges(start=seqcbs_pos1,width=seqcbs_qwidth1))
    #coverage=mean(seqcbs_coverage[seqcbs_coverage>0])
    #sp_scores=as.data.frame(list(lCBS=seqcbs_stat,lCVG=seqcbs_coverage))
    #sp_scores=set_colClass(sp_scores,rep("numeric",2))
  }
}

scan_bam = function(bam_file, what, isize_global, isize_sdR_global, isize_sdL_global,
                    coverage_global, seq_name,
                    start, end, RL_global, gap, trunk_size,
                    width, lw_width, mixing_rate, stepsize,
                    block_size, Delta_global, bigDel_global, 
                    smallDel_global, smallIns_global, maxInsert_global,
                    p_global, q_global,
                    hang_clip_global, prop_clip_global,
                    debug, rg_out, verbose, fast, gtk){ 
  #list[bam_file,what,isize_global,isize_sdR_global,isize_sdL_global,coverage_global,seq_name,start,end,RL_global,gap,trunk_size,width,mixing_rate,stepsize,block_size,Delta_global,bigDel_global,p_global,q_global,hang_clip_global,prop_clip_global,debug,out,verbose,fast,gtk]=scan_par
  if(verbose) cat("--Info: in scan_bam: start:",nsf(start),"end:",nsf(end),"stepsize:",nsf(stepsize), "\n")
  maxDelta=max(delta_fail,Delta_global)
  if(verbose) cat("==Info: maxDelta",nsf(maxDelta),"\n")
  wplus=max(width,0); lw_width=max(lw_width,lw_width_fail)
  Delta_trunk=c(); hang_clip_trunk=c(); prop_clip_trunk=c(); bigDel_trunk=c(); 
  RL_trunk=c(); isize_trunk=c(); isize_sdR_trunk=c(); isize_sdL_trunk=c(); smallDel_trunk=c()
  p_left_trunk=c(); p_right_trunk=c(); q_left_trunk=c(); q_right_trunk=c(); smallIns_trunk=c()
  coverage_trunk=c(); rlen_trunk=c(); maxInsert_trunk=c();
  srange=seq(start,end,stepsize); n_wins=length(srange); r_start_trunk=c(); r_end_trunk=c();
  if(verbose) cat("==Info: region chr=",seq_name,"start=",nsf(start),"end=",nsf(end),
                  "size=", nsf(end-start+1),"step=",nsf(stepsize),"n_wins=",nsf(n_wins),"\n")
  win_per_trunk=trunk_size/stepsize #this is ensured to be multiples
  n_trunk=ceiling(n_wins/win_per_trunk) #1001 windows, 100 w_p_t => 11 trunk
  if(verbose) cat("==Info: trunk_size=", nsf(trunk_size),
                  "max_win_per_trunk=",nsf(win_per_trunk),"n_trunk=",nsf(n_trunk),"\n")
  scan_cols=scan_cols_preset
  scan_result=data.frame(matrix(rep(0,n_wins*length(scan_cols)),n_wins,length(scan_cols)))
  colnames(scan_result)=scan_cols
  scan_bigd=data.frame(); scan_disc=data.frame(); scan_anch=data.frame()
  
  for(t in seq_len(n_trunk)){ #t=1
    t_start=start+(t-1)*trunk_size; t_end=min(start+t*trunk_size-1,end)
    t_win_start=(t-1)*win_per_trunk+1; t_win_end=min(n_wins,t*win_per_trunk)
    t_srange=srange[t_win_start:t_win_end]
    cat("==Info: scan_bam: t_start=",nsf(t_start),"t_end=",nsf(t_end),";", t,"-th of", n_trunk, "\n") 
    if(verbose) cat("==Info: t_win_start=",nsf(t_win_start),"t_win_end=",nsf(t_win_end),"\n") 
    if(verbose) cat("==Info: t_srange_start=",nsf(t_srange[1]),
                    "t_srange_end=",nsf(t_srange[length(t_srange)]+stepsize-1),"\n") 
    if(verbose) cat("==Info: scanning trunk from ", nsf(t_start), "to", 
                    nsf(t_end), "with stepsize", stepsize, "\n")
    if(verbose) cat("==Info: getting trunk reads from",nsf(max(t_start-maxDelta,0)),
                    "to",nsf(min(t_end,end)+maxDelta),"\n")
    seq_info=GRanges(seq_name,IRanges(max(t_start-maxDelta,0),min(t_end,end)+maxDelta))
    list[rPbi,rPbo,rPbd,rPc,rMc,rPn,rMn,rSr,rSl,rHr,rHl,
         rMHp,rMHn,rMDp,rMDn,rSCHp,rSCHn,
         rSCISp,rSCISn,rSDa,rSIa,rSCPbi,rSCPbo,
         hang_clip_est,prop_clip_est, 
         RL_est,isize_est,isize_sd_est,Delta_est,smallDel_est,smallIns_est,bigDel_est,maxInsert_est,t_anch]=
      get_MPRs(seq_info,bam_file,RL_global,Delta_global,bigDel_global,smallDel_global,smallIns_global,maxInsert_global,
                       hang_clip_global,prop_clip_global,wplus,gtk,debug,verbose)
    if(verbose) cat("---Info: out get_MPRs, start=",start(ranges(seq_info)),"end=",end(ranges(seq_info)),"\n")
    if(verbose) cat("==Info: out get_MPRs: ",t,"-th trunk of ",n_trunk," bam trunk,",taggie(gtk),"\n")
    
#     print(list(rPbi=rPbi,rPbo=rPbo,rPbd=rPbd,rPc=rPc,rMc=rMc,rPn=rPn,rMn=rMn,rSr=rSr,rSl=rSl,rHr=rHr,rHl=rHl,
#                        rMHp=rMHp,rMHn=rMHn,rMDp=rMDp,rMDn=rMDn,rSCHp=rSCHp,rSCHn=rSCHn,
#                        rSCISp=rSCISp,rSCISn=rSCISn,rSDa=rSDa,rSIa=rSIa,rSCPbi=rSCPbi,rSCPbo=rSCPbo,
#                        hang_clip_est=hang_clip_est,prop_clip_est=prop_clip_est, 
#                        RL_est=RL_est,isize_est=isize_est,isize_sd_est=isize_sd_est,
#                        Delta_est=Delta_est,smallDel_est=smallDel_est,smallIns_est=smallIns_est,bigDel_est=bigDel_est))
    #trunk parameter controls
    if(isize_global==-1) isize=isize_est else isize=isize_global
    if(isize_sdR_global==-1) isize_sdR=isize_sd_est else isize_sdR=isize_sdR_global
    if(isize_sdL_global==-1) isize_sdL=isize_sd_est else isize_sdL=isize_sdL_global
    isize_sd=max(isize_sdR,isize_sdL)
    if(RL_global==-1) RL=RL_est else RL=RL_global
    if(Delta_global==-1) Delta=Delta_est else Delta=Delta_global
    if(bigDel_global==-1) bigDel=bigDel_est else bigDel=bigDel_global
    if(smallDel_global==-1) smallDel=smallDel_est else smallDel=smallDel_global
    if(smallIns_global==-1) smallIns=smallIns_est else smallIns=smallIns_global
    if(maxInsert_global==-1) maxInsert=maxInsert_est else maxInsert=maxInsert_global
    if(hang_clip_global==-1) hang_clip=hang_clip_est else hang_clip=hang_clip_global
    if(prop_clip_global==-1) prop_clip=prop_clip_est else prop_clip=prop_clip_global

    Delta_trunk=c(Delta_trunk,Delta); hang_clip_trunk=c(hang_clip_trunk,hang_clip) 
    prop_clip_trunk=c(prop_clip_trunk,prop_clip); 
    bigDel_trunk=c(bigDel_trunk,bigDel) 
    smallDel_trunk=c(smallDel_trunk,smallDel)
    smallIns_trunk=c(smallDel_trunk,smallIns)
    maxInsert_trunk=c(maxInsert_trunk,maxInsert)
    RL_trunk=c(RL_trunk,RL); 
    isize_trunk=c(isize_trunk,isize); 
    isize_sdR_trunk=c(isize_sdR_trunk,isize_sdR)
    isize_sdL_trunk=c(isize_sdL_trunk,isize_sdL)
    
    if(fast=="super") { speed_factor=2
      } else if(fast=="fast") { speed_factor=1
      } else { speed_factor=0 }
    if(length(rPbi)!=0)
      rPbi=rPbi[width(rPbi)<=(isize-speed_factor*isize_sdL) | width(rPbi)>=(isize+speed_factor*isize_sdR)]
  
    #do the trunk scan results
    joint_par=c(list(srange=t_srange),width=width,lw_width=lw_width,stepsize=stepsize,block_size=block_size,
                Delta=Delta, RL=RL, maxInsert=maxInsert, p_input=p_global, q_input=q_global, list(scan_cols=scan_cols[1:(length(scan_cols)-4)]),
                isize=isize, isize_sdR=isize_sdR, isize_sdL=isize_sdL, coverage=coverage_global,
                mixing_rate=mixing_rate, list(rPbi=c(rPbi,rSCPbi),rPbo=c(rPbo,rSCPbo),rPbd=rPbd,
                rPc=rPc,rMc=rMc,rPn=rPn,rMn=rMn,rHr=rHr,rHl=rHl,rSr=rSr,rSl=rSl),debug=debug, verbose=verbose)
    #here maybe a little bit less isize with rSCPbi due to softclipping but not too much
    if(verbose) cat("==Info: scan pars:", "width=", width, "Delta=", Delta, "RL=", RL, 
        "isize=", isize, "isize_sdR=", isize_sdR, "isize_sdL=",isize_sdL,"mixing_rate=", mixing_rate, "\n")
    if(length(rMc)==0 | length(rPc)==0 | length(rPbi)==0){ #skip scan
      cat("==Warn: no MPRs in this trunk, skipping...\n")
      t_result=data.frame(matrix(rep(NA,(t_win_end-t_win_start+1)*length(scan_cols)),ncol=length(scan_cols)))
      t_result[,1]=t_srange
      t_par=list(p_left=NA,p_right=NA,q_left=NA,q_right=NA,coverage=NA,r_start=NA,r_end=NA,
                 is=NA,issdR=NA,issdL=NA,rl=NA,maxInsert=NA)
      t_bigd=data.frame(); t_disc=data.frame()
    } else { #in fixed order
      list[t_result, t_par, t_bigd, t_disc]=do.call(scan_joint, joint_par) #t_anch is determined
      t_result=cbind(t_result,ins=IRanges::as.vector(IRanges::coverage(rSIa))[t_srange]) #ins
      t_result=cbind(t_result,del=IRanges::as.vector(IRanges::coverage(rSDa))[t_srange]) #del
      t_result=cbind(t_result,HAF=IRanges::as.vector(IRanges::coverage(rMHp))[t_srange]) #HAF
      t_result=cbind(t_result,HAR=IRanges::as.vector(IRanges::coverage(rMHn))[t_srange]) #HAR
    }
    scan_result[t_win_start:t_win_end,]=t_result
    scan_bigd=rbind(scan_bigd,t_bigd)
    scan_disc=rbind(scan_disc,t_disc)
    scan_anch=rbind(scan_anch,t_anch)
    #additional parameters to record
    p_left_trunk=c(p_left_trunk,t_par$p_left);p_right_trunk=c(p_right_trunk,t_par$p_right) 
    q_left_trunk=c(q_left_trunk,t_par$q_left);q_right_trunk=c(q_right_trunk,t_par$q_right) 
    coverage_trunk=c(coverage_trunk,t_par$coverage)
    coverage_est=t_par$coverage
    r_start_trunk=c(r_start_trunk,t_par$r_start)
    r_end_trunk=c(r_end_trunk,t_par$r_end)
    rlen_trunk=c(rlen_trunk,t_par$r_end-t_par$r_start)
    if(debug) { 
      if(verbose) cat("--Info: saving trunk data\n"); 
      # following is to enforce trunk.RData contents are strcitly within start to end
      if(verbose) cat("--Info: trunk range:", t_start, "-", t_end, "\n")
      rPbi=rPbi[start(rPbi)>=t_start & start(rPbi)<t_end]
      rPbo=rPbo[start(rPbo)>=t_start & start(rPbo)<t_end]
      rMc=rMc[start(rMc)>=t_start & start(rMc)<t_end]
      rPc=rPc[start(rPc)>=t_start & start(rPc)<t_end]
      rSr=rSr[start(rSr)>=t_start & start(rSr)<t_end]
      rSl=rSl[start(rSl)>=t_start & start(rSl)<t_end]
      rHr=rHr[start(rHr)>=t_start & start(rHr)<t_end]
      rHl=rHl[start(rHl)>=t_start & start(rHl)<t_end]
      rMHp=rMHp[start(rMHp)>=t_start & start(rMHp)<t_end]
      rMHn=rMHn[start(rMHn)>=t_start & start(rMHn)<t_end]
      rMDp=rMDp[start(rMDp)>=t_start & start(rMDp)<t_end]
      rMDn=rMDn[start(rMDn)>=t_start & start(rMDn)<t_end]
      rSCHp=rSCHp[start(rSCHp)>=t_start & start(rSCHp)<t_end]
      rSCHn=rSCHn[start(rSCHn)>=t_start & start(rSCHn)<t_end]
      rSCISp=rSCISp[start(rSCISp)>=t_start & start(rSCISp)<t_end]
      rSCISn=rSCISn[start(rSCISn)>=t_start & start(rSCISn)<t_end]
      rSCPbo=rSCPbo[start(rSCPbo)>=t_start & start(rSCPbo)<t_end]
      rSCPbi=rSCPbi[start(rSCPbi)>=t_start & start(rSCPbi)<t_end]
      rSDa=rSDa[start(rSDa)>=t_start & start(rSDa)<t_end]
      rSIa=rSIa[start(rSIa)>=t_start & start(rSIa)<t_end]
      save(list=c("rPbi","rPbo","rMc","rPc","rSr","rSl","rHr","rHl",
          "rMHp","rMHn","rMDp","rMDn","rSCHp","rSCHn",
          "rSCISp","rSCISn","rSCPbo","rSCPbi","rSDa","rSIa",
          "hang_clip_est","prop_clip_est",
          "coverage_est","RL_est","isize_est","isize_sd_est",
          "Delta_est","smallDel_est","bigDel_est"),
          file=paste(rg_out,".trunk",t,".RData",sep=""))
    }
    if(debug) cat("==Info: done trunk scan, gtk=",get_gtk(gtk),"s,gmk=",get_gmk(Sys.getpid()),"kb\n") 
    if(verbose) cat("== out of trunk ...\n")
  }

  if(!is.null(gap)){
    gap_found=gap[which(as.character(gap$chrom)==seq_name),] #centromere region of seq_name, about 3M
    if(verbose) cat("==Info: found", nrow(gap_found), "gaps in chr", seq_name, "\n")
    if(verbose) cat("==Info: gaps:\n"); if(verbose) print(gap_found)
    if(verbose) cat("==Info: start=",start,"end=",end,"\n")
    for(i in seq_len(nrow(gap_found))){
      g_start=gap_found$chromStart[i]+1; g_end=gap_found$chromEnd[i] #start 0-based, end 1-based
      if(!(g_end<start) & !(g_start>srange[length(srange)])){ #this ensure no 0 return from bsearch
        #g_w_start=bsearch1(g_start,srange); #g_w_start=max(g_w_start,1);
        g_w_start=max(0,floor((g_start-start)/stepsize))+1
        #g_w_end=bsearch1(g_end,srange); #g_w_end=max(g_w_end,1)
        g_w_end=max(0,floor((g_end-start)/stepsize))+1
      if(verbose) cat("==Info: i=", i, "g_start=",g_start,"g_end=",g_end, "\n")
      if(verbose) cat("==Info: i=", i, "g_w_start=",g_w_start,"g_w_end=",g_w_end, "\n")
        for(s in scan_cols[seq_along(scan_cols)])
          if(!(s %in% c("start","cvg"))) scan_result[[s]][g_w_start:g_w_end]=NA
      }
    }
  }
  
  mvsum=mvsum_fail
  flank=flank_fail
  #if(verbose) cat("r_start_trunk=",r_start_trunk,"\n")
  #if(verbose) cat("r_end_trunk=",r_end_trunk,"\n")
  success=T
  if(all(is.na(r_start_trunk)) | all(is.na(r_end_trunk))){
    success=F; r_start=NA; r_end=NA; w_start=n_wins; w_end=1
  } else {
    r_start=min(r_start_trunk,na.rm=T); r_end=max(r_end_trunk,na.rm=T); #will produce inf if trunk is NA
    w_start=max(1,ceiling((r_start+flank-start)/stepsize),na.rm=T)
    w_end=min(n_wins,n_wins-ceiling((end-r_end+flank)/stepsize),na.rm=T)
  }
  for(s in scan_cols){
    if(!(s %in% c("start","cvg"))) { scan_result[[s]][1:w_start]=NA; scan_result[[s]][w_end:n_wins]=NA }
  }
  if(verbose) cat("==Info: r_start=",nsf(r_start),"start=",nsf(start),"w_start=",nsf(w_start),"\n")
  if(verbose) cat("==Info: r_end=",nsf(r_end),"end=",nsf(end),"w_end=",nsf(w_end),"\n")
  if(verbose) cat("==Info: coverage_trunk=\n")
  if(verbose) print(coverage_trunk) 
  if(verbose) cat("==Info: rlen_trunk=\n")
  if(verbose) print(rlen_trunk) 
  if(verbose) cat("==Info: maxInsert_trunk=\n")
  if(verbose) print(maxInsert_trunk) 
  samllDel=round(sum(smallDel_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T))
  maxInsert=min(max(c(maxInsert_trunk,0),na.rm=T),maxInsert_fail)
  bigDel=round(sum(bigDel_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T))
  p_left=round(sum(p_left_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T),digits=6)
  p_right=round(sum(p_right_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T),digits=6)
  q_left=round(sum(q_left_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T),digits=6)
  q_right=round(sum(q_right_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T),digits=6)
  coverage=round(sum(coverage_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T),digits=4)
  rl=round(sum(RL_trunk*(rlen_trunk/sum(rlen_trunk,na.rm=T)),na.rm=T))
  lambda=round(coverage/rl,digits=2)
  if(verbose) cat("==Info: coverage=",coverage,"lambda=",lambda,"\n") 

  swan_par=as.data.frame(list(
    delta=round(sum(Delta_trunk*(rlen_trunk/sum(rlen_trunk,na.rm=T)),na.rm=T)),
    hang_clip=round(sum(hang_clip_trunk*(rlen_trunk/sum(rlen_trunk,na.rm=T)),na.rm=T)),
    prop_clip=round(sum(prop_clip_trunk*(rlen_trunk/sum(rlen_trunk,na.rm=T)),na.rm=T)),
    rl=rl, coverage=coverage,
    isize=round(sum(isize_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T)),
    isize_sdR=round(sum(isize_sdR_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T)),
    isize_sdL=round(sum(isize_sdL_trunk*rlen_trunk,na.rm=T)/sum(rlen_trunk,na.rm=T)),
    smallDel=smallDel, smallIns=smallIns, bigDel=bigDel, maxInsert=maxInsert, p_left=p_left, p_right=p_right, q_left=q_left, q_right=q_right,
    start=start, end=end, chr=seq_name, w=width, lw_width=lw_width, r=mixing_rate,
    lambda=lambda, r_start=r_start, r_end=r_end, success=success,
    trunk_size=trunk_size, block_size=block_size, n_wins=n_wins,
    speed_factor=speed_factor, stepsize=stepsize, fy_cap=fy_cap_fail))
  
  #dedup big list
  #scan_bigd=scan_bigd[!duplicated(scan_bigd),]; 
  #scan_disc=scan_disc[!duplicated(scan_disc),];
  if(verbose) cat("--Info: out scan_bam ...  \n")
  #print(scan_result$lCd_term[1:10])
  #print(scan_result$cCd_term[1:10])
  return(list(scan_result=scan_result,scan_par=swan_par,scan_bigd=scan_bigd,scan_disc=scan_disc,scan_anch=scan_anch))
} # main function to scan bam file

#cat("loaded line 4500... \n")

### testing functions ###

rcpp_hello_world <- function(){
  .Call( "rcpp_hello_world", PACKAGE = "swan" )
}

### libsv.R ###
theo_thresh=function(score,alpha,scan_par,verbose=F){
  #R1, R2, kappa, p are variable
  #R1, R2: pairs with + strand mapping within [0, t-R1] and - strand within [t+w-R2,T] are counted
  #need to take into account whether the library is enabled or not here
  #libnum=nrow(scan_par) #only true are counted toward libnum
  r_start=scan_par$r_start[1]; r_end=scan_par$r_end[1]; stepsize=scan_par$stepsize[1] #assumed identical
  nbases=r_end-r_start; data_mean=NULL; data_sd=NULL; status=scan_par[[score]]
  R=scan_par$rl; R1=R-scan_par$hang_clip; R2=R-scan_par$hang_clip; #R must be a constant
  p_left=scan_par$p_left; p_right=scan_par$p_right;
  q_left=scan_par$q_left; q_right=scan_par$q_right;
  lambda=scan_par$coverage/scan_par$rl; N=lambda*nbases; n=round(nbases/stepsize)
  kappa=sqrt(lambda); kappa_wc=sqrt(lambda*scan_par$rl/stepsize)
  w=scan_par$w; r=scan_par$r; mu=scan_par$isize; sigmaR=scan_par$isize_sdR; sigmaL=scan_par$isize_sdL;
  D=scan_par$delta; nullmean=0; nullvar=0; nullstdev=0
  if(verbose) cat("---Info: computing thresh for score=",score,"\n")
  if(score=="lDl"|score=="lDr"){
    if(score=="lDr") cp=p_right else cp=p_left; cpbase=cp*N/nbases
      tryCatch({
        b=getThresholdHangingReadsScan(alpha=alpha,kappa=kappa,w=median(w),r=median(r),p=cp,pbase=cpbase,mu=mu,sigma=sigmaR,n=n,D=D,R=median(R),UPPER=NA,verbose=FALSE,status=status)
      }, error=function(e) { print(e); cat("---Info: forcing 20\n") }, finally={ b=20 })
    #b=getThresholdHangingReadsScan(alpha=0.01,kappa=c(0.4582576,0.3316625,0.3316625),w=20,r=0.1,p=c(0.082,0.095,0.107),pbase=c(0.01722,0.01045,0.01177),mu=c(258,368,538),sigma=c(23,29,54),n=10000009,D=c(2000,2000,2000),R=100,UPPER=NA, verbose=FALSE)
    #pvalHangingReadsScan(b=b,kappa=kappa,w=median(w),r=median(r),p=cp,pbase=cpbase,mu=mu,sigma=sigmaR,n=n,D=D,R=median(R),verbose=verbose)
  } else if(score=="lCd"){
    b=getThresholdStraddlingReadsDeletionScan(alpha=alpha,kappa=kappa,w=median(w),r=median(r),mu=mu,sigma=sigmaR,n=n,R1=median(R1),R2=median(R2),UPPER=NA,verbose=FALSE,status=status)
    #pvalStraddlingReadsDeletionScan(b=b,kappa=kappa,w=median(w),r=median(r),mu=mu,sigma=sigmaR,n=n,R1=median(R1),R2=median(R2),UPPER=NA,verbose=FALSE)
  } else {
    b=0
  }
  return(list(thresh=b,data_mean=nullmean,data_sd=nullstdev))
}

get_lWc_thresh=function(alpha=0.1,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R,n=n,verbose=F){
  nullmean=psidotbetaCoverageScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,verbose=verbose)
  nullvar=sigmabetaCoverageScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,varZt=TRUE,verbose=verbose)
  nullstdev = sqrt(nullvar)
  b=getThresholdCoverageScan(alpha,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n,UPPER=NA,verbose=verbose)
  #pval=pvalCoverageScan(b=b,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n,UPPER=NA,everbose=F)
  #if(everbose) cat("====Info: pval of the threshold ", b, " of alpha= ", alpha, " is ", pval, "\n",file=stderr())
  return(list(thresh=b,data_mean=nullmean,data_sd=nullstdev))
}

get_lD_thresh=function(alpha=0.1,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,n=n,D=D,R=R,verbose=F){
  nullmean=psidotbetaHangingReadsScan(beta=0,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,D=D,R=R)
  nullvar=sigmabetaHangingReadsScan(beta=0,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,D=D,R=R,varZt=TRUE)
  nullstdev=sqrt(nullvar)
  b=getThresholdHangingReadsScan(alpha,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,n=n,D=D,R=R,UPPER=NA,verbose=F)
  #pval=pvalHangingReadsScan(b,kappa,w,r,p,pbase,mu,sigma,n,D,R,everbose=F)
  #if(verbose) cat("====Info: pval of the threshold ", b, " of alpha= ", alpha, " is ", pval, "\n",file=stderr())
  return(list(thresh=b,data_mean=nullmean,data_sd=nullstdev))
}

get_lCd_thresh=function(alpha=0.1,kappa=kappa,w=w,r=r,p=p,mu=mu,sigma=sigma,D=D,n=n,R=R,R1=R1,R2=R2,verbose=F){
  #cat("getting lCd threshold... \n")
  nullmean=psidotbetaStraddlingReadsDeletionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2,IntLim=6)
  nullvar=sigmabetaStraddlingReadsDeletionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2, varZt=TRUE)
  nullstdev = sqrt(nullvar)
  b=getThresholdStraddlingReadsDeletionScan(alpha,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=NA, verbose=verbose, useEmpiricalDist=TRUE)
  #pval=pvalStraddlingReadsDeletionScan(b=b,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=NA,verbose=F)
  #if(verbose) cat("====Info: pval of the threshold ", b, " of alpha= ", alpha, " is ", pval, "\n",file=stderr())
  return(list(thresh=b,data_mean=nullmean,data_sd=nullstdev))
}

get_lCi_thresh=function(alpha=0.1,kappa=kappa,w=w,r=r,p=p,mu=mu,sigma=sigma,n=n,D=D,R=R,R1=R1,R2=R2,verbose=F){
  nullmean=psidotbetaStraddlingReadsInsertionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2)
  nullvar=sigmabetaStraddlingReadsInsertionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2,varZt=TRUE)
  nullstdev = sqrt(nullvar)
  b=getThresholdStraddlingReadsInsertionScan(alpha,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2,UPPER=NA,verbose=verbose)
  #pval=pvalStraddlingReadsInsertionScan(b=b,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=NA,everbose=F)
  #if(verbose) cat("====Info: pval of the threshold ", b, " of alpha= ", alpha, " is ", pval, "\n",file=stderr())
  return(list(thresh=b,data_mean=nullmean,data_sd=nullstdev))
}


### charlie's work on Coverage Scan ###
getThresholdCoverageScan<-function(alpha,kappa,w,r,mu,sigma,n,D,R,UPPER=NA,verbose=F){
  pvalminusalpha = function(b){ pvalCoverageScan(b=b,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n,D=D,R=R,UPPER=UPPER,verbose=verbose)-alpha}
  cat("\nComputing threshold for Hanging Reads Scan for FWER=",alpha,"...\n")
  
  m0=psidotbetaCoverageScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R)
  s0=sigmabetaCoverageScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R, varZt=TRUE)
  s0=sqrt(s0)
  
  if(verbose) cat("Null mean = ",m0,", standard deviation = ",s0,".\n",file=stderr())
  
  K=6
  prevroot = m0
  lower = m0+0.1*s0
  while(TRUE){
    upper = m0+K*s0
    if(verbose) cat("\nTrying to find threshold with search upper limit ", K," sd from mean: [",lower,", ",upper,"]\n",sep="",file=stderr())
    res=try(uniroot(pvalminusalpha, lower=m0+0.1*s0, upper=m0+K*s0), silent=TRUE)
    if(class(res) == "try-error"){
      if(verbose) cat("Upper limit too low, multiply K by 1.5.\n",file=stderr())
      K = K*1.5
      if(!is.finite(K)) stop("Stop Because of Infinity Loop!")
    } else {
      if(verbose) cat("\n-----> uniroot converged on b=",res$root,".\n", file=stderr())
      break
    }
  }
  cat("Threshold value: ",res$root,"\n")
  res$root
}

pvalCoverageScan<-function(b,kappa,w,r,mu,sigma,n,D,R, UPPER=NA,verbose=F){
  if(verbose) cat("\nComputing P-VALUE for Hanging Reads Scan with b=",b,"...\n",file=stderr())
  
  ### Get betab = root of psidotminusb.
  psidotminusb = function(beta){ psidotbetaCoverageScan(beta=beta,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R)-b}
  if(is.na(UPPER)){
    UPPER = 2
  }
  if(psidotminusb(0)>=0){
    # Must start with psidotminusb(0)<0.
    stop("Threshold value must be larger than null mean.\n")
  }
  prevroot=0
  while(TRUE){
    if(verbose) cat("\tTrying to compute betab with search upper limit ", UPPER,"...\n",sep="",file=stderr())
    res=try(uniroot(psidotminusb, lower=0, upper=UPPER), silent=TRUE)
    if(class(res) == "try-error"){
      if(psidotminusb(UPPER)<0){
        # Must have returned error because psidotminusb(UPPER) is also negative.
        # uniroot must have psidotminusb(0)<0 (enforced above) and psidotminusb(UPPER)>0
        UPPER = UPPER*2
        if(verbose) cat("\tUpper limit too low and psidotminusb<0, increase upper to ", UPPER,".\n",sep="",file=stderr())
      } else{
        # Must have returned error because psidotminusb(UPPER) is infinity.
        if(verbose) cat("Upper limit too high, divide by 2.\n",file=stderr())
        UPPER = UPPER/1.5            
      }
    } else {
      break
    }
  }
  
  betab=res$root
  if(verbose) cat("Final value of beta_b: ",betab,".\n",sep="",file=stderr())
  
  sigmab = sigmabetaCoverageScan(beta=betab,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R)
  deltab = deltabetaCoverageScan(beta=betab,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,D=D,R=R)
  psib = psibetaCoverageScan(beta=betab,kappa=kappa,w=w,r=r,p=p,mu=mu,sigma=sigma,D=D,R=R)
  
  pval = n* exp(-(betab*b-psib))*(1/(sqrt(2*pi*sigmab)))*deltab
  if(verbose) cat("P-value=",pval,", betab=",betab,", sigmab=",sigmab,", deltab=",deltab, ", psib=",psib,"\n",file=stderr())
  pval
}

psidotbetaCoverageScan<-function(beta,kappa,w,r,mu,sigma,D,R){
  #region 2
  f2 = function(z){ (w-z)*dnorm(z, mean=mu, sd=sigma) }
  reg2 = kappa^2*(1-r)^beta*(w-integrate(f2, lower=0, upper=w)$value)
  #region 4
  f4 = function(z){ dnorm(z, mean=mu, sd=sigma) }
  reg4 = kappa^2*(1-r)^(2*beta)*w*integrate(f4, lower=0, upper=w)$value
  #region 5
  f5 = function(z){ z*dnorm(z, mean=mu, sd=sigma) }
  reg5 = kappa^2*(1-r)^beta*(w-integrate(f5, lower=0, upper=w)$value)
  res=log(1-r)*(reg2+reg5+2*reg4)
  res
}

psibetaCoverageScan<-function(beta,kappa,w,r,mu,sigma,D,R){
  #region 1 + region 3 + region 6 - k^2T=0
  #region 2
  f2 = function(z){ (w-z)*dnorm(z, mean=mu, sd=sigma) }
  reg2 = kappa^2*(1-r)^beta*(w-integrate(f2, lower=0, upper=w)$value)
  #region 4
  f4 = function(z){ dnorm(z, mean=mu, sd=sigma) }
  reg4 = kappa^2*(1-r)^(2*beta)*w*integrate(f4, lower=0, upper=w)$value
  #region 5
  f5 = function(z){ z*dnorm(z, mean=mu, sd=sigma) }
  reg5 = kappa^2*(1-r)^beta*(w-integrate(f5, lower=0, upper=w)$value)
  res=reg2+reg4+reg5
  res
}

sigmabetaCoverageScan<-function(beta,kappa,w,r,mu,sigma,D,R,varZt=FALSE){
  #region 2
  f2 = function(z){ (w-z)*dnorm(z, mean=mu, sd=sigma) }
  reg2 = kappa^2*(1-r)^beta*(w-integrate(f2, lower=0, upper=w)$value)
  #region 4
  f4 = function(z){ dnorm(z, mean=mu, sd=sigma) }
  reg4 = kappa^2*(1-r)^(2*beta)*w*integrate(f4, lower=0, upper=w)$value
  #region 5
  f5 = function(z){ z*dnorm(z, mean=mu, sd=sigma) }
  reg5 = kappa^2*(1-r)^beta*(w-integrate(f5, lower=0, upper=w)$value)
  res=log(1-r)^2*(reg2+reg5+4*reg4)
  if(!varZt){
    sigmabeta = beta^2*res
  } else {
    sigmabeta = res
  }
  sigmabeta
}

deltabetaCoverageScan<-function(beta,kappa,w,r,p,pbase,mu,sigma,D,R){
  f1 = function(z){ dnorm(z, mean=mu, sd=sigma) }
  deltabeta = 2*kappa^2*log(1-r)*(1-r)^(beta)*(integrage(f1,lower=0,upper=w)$value+integrate(f1,lower=-w,upper=0)$value)
  deltabeta
}

# ### sclip scan routine adopted from Nancy ###
core_sclip <- function(ti,n_trunks,scan_start,scan_end,trunk_size,files,seqname,gap,MIN.READS.PER.CLUSTER=5,MIN.BASES.PER.CLUSTER=50,SC.PROPORTION.MISMATCH.THRESH=0.05,MIN.GAP.DISTANCE=10000,debug=FALSE){
   what=c("strand", "pos","seq","cigar")
   pos=c(); seqs=c(); strand=c(); cigar=c(); 
   st0 = scan_start+(ti-1)*trunk_size
   ed0 = min(scan_start+ti*trunk_size-1,scan_end)
   if(debug) cat("----- in trunk ",ti," of ",n_trunks,"from",nsf(st0),"to",nsf(ed0),taggie(gtk),"\n")
   seq_info=GRanges(seqname, IRanges(st0, ed0))
   for(bam in files){
     reads=allFunction(seq_info, bam, what=what, index=bam)
     pos=c(pos,reads$pos); strand=c(strand,reads$strand)
     cigar=c(cigar,reads$cigar)
     seqs=if(is.null(seqs)) reads$seq else c(seqs,reads$seq)
   } #we combine gi_reads from multiple bamfiles
   sc=getSoftClipClusters(pos=pos,cigar=cigar,seq=seqs,strand=strand,minpos=st0,maxpos=ed0,MIN.READS.PER.CLUSTER=MIN.READS.PER.CLUSTER,MIN.BASES.PER.CLUSTER=MIN.BASES.PER.CLUSTER)
   tabL =summarizeClusters(sc$LsoftclipClusters,includeNBases=TRUE)
   tabR =summarizeClusters(sc$RsoftclipClusters,includeNBases=TRUE)
   sel = which(tabL$ProportionMismatch<SC.PROPORTION.MISMATCH.THRESH)
   scL = sc$LsoftclipClusters[sel]
   sel = which(tabR$ProportionMismatch<SC.PROPORTION.MISMATCH.THRESH)
   scR = sc$RsoftclipClusters[sel]
   if(debug) cat("----- out of trunk ",ti,"of",n_trunks,taggie(gtk),"\n")
   if(debug) cat("\nLeft-clipped: ",length(scL),", right-clipped: ",length(scR),"\n",sep="")
   list(scL=scL,scR=scR)
} #core_sclip
#This is only softsclipping clusters by region, sc0L and sc0R not matched, same true for sc1L and sc1R

plotRegions<-function(bamfile1,chr,st0,ed0,st=NA,ed=NA,mfrow=TRUE,label="",readlen=100,doCoverage=TRUE,bamfile0=NULL,doIL=TRUE,doHanging=TRUE,doSoftClip=TRUE,covymax=NA){
    # Prepare for getting the data.
    if(mfrow) par(mfrow=c(doCoverage+doIL+doHanging+doSoftClip,1))
    header=scanBamHeader(bamfile1[1])
    chrnames=names(header[[1]]$targets)  # find out what the chromosomes are called.
    what = c("qname","rname", "flag","strand", "pos","mapq","mrnm","mpos","isize","qwidth","seq","cigar","qual")

    ### Get the data.
    which = RangesList("quack"=IRanges(st0,ed0))
    names(which) = chrnames[chr]
    param = ScanBamParam(which=which, what=what)
    bam1 = list(); for(w in what) { bam1[[w]]=NULL }
    for(bi in seq_along(bamfile1)){
      tmp = scanBam(bamfile1[bi], param=param)[[1]]
      for(w in what) { 
        if(w %in% c("seq","qual")) {
          if(is.null(bam1[[w]])) bam1[[w]]=tmp[[w]] else bam1[[w]]=c(bam1[[w]],tmp[[w]]) 
        } else if(w %in% c("strand")) {
          bam1[[w]]=c(bam1[[w]],as.character(tmp[[w]]))
        } else {
          bam1[[w]]=c(bam1[[w]],tmp[[w]])
        }
      }
    }
    if(length(bam1$pos)==0) return(NULL)
    fl=parseFlags(bam1$flag)
    if(length(bam1$pos)<10) return(bam1)

    ### Plot everything.
    res=getReadMapPositions(bam1,fl)
    if(doCoverage){
        positions = IRanges(start=c(st0:ed0),width=1)
        sel = which(!is.na(bam1$pos))
        readst=bam1$pos[sel]
        sel = which(!is.na(bam1$mpos))
        readst = c(readst,bam1$mpos[sel])
        readir = IRanges(start=readst,width=readlen)
        cov1=countOverlaps(positions,readir)
        if(is.null(bamfile0)){
            mainstr=paste("Base Coverage,", label,sep=" ")
            cat(mainstr,"\n")
            plot(positions,cov1,type="b",xlab="Position",ylab="Read Coverage",main=mainstr,cex.main=2)
            if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
            if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)   
            grid()
        } else {
            mainstr=paste("Fold Change in Coverage,",label,sep=" ")
            cat(mainstr,"\n")
            what=c("pos","mpos")
            param = ScanBamParam(which=which, what=what)
            bam0 = list(); for(w in what) { bam0[[w]]=NULL }
            for(bi in seq_along(bamfile0)){
              tmp = scanBam(bamfile0[bi], param=param)[[1]]
              for(w in what) { 
                if(w %in% c("seq","qual")) {
                  if(is.null(bam0[[w]])) bam0[[w]]=tmp[[w]] else bam0[[w]]=c(bam0[[w]],tmp[[w]]) 
                } else if(w %in% c("strand")) {
                  bam0[[w]]=c(bam0[[w]],as.character(tmp[[w]]))
                } else {
                  bam0[[w]]=c(bam0[[w]],tmp[[w]])
                }
              }
            }
            sel = which(!is.na(bam0$pos))
            readst=bam0$pos[sel]
            sel = which(!is.na(bam0$mpos))
            readst = c(readst,bam0$mpos[sel])
            readir = IRanges(start=readst,width=readlen)
            cov0=countOverlaps(positions,readir)
            fc=cov1/pmax(cov0,1)
            if(is.na(covymax)) covymax=max(fc,0,na.rm=TRUE)
            cat("data points for coverage:", length(positions), " ", length(fc),"\n")
            plot(positions,fc,type="b",xlab="Position",main=mainstr,ylim=c(0,covymax),ylab="Fold change",cex.main=2.,cex.lab=2)
            if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
            if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)   
            grid()
        }
    }
    if(doIL){
        cat("Mapped Insert Sizes\t")
        mainstr=paste("Mapped Insert Sizes")
        #cat("data points for insert:", length(res$rPstart), " ", length(res$rMstart),"\n")
        plotInsertLengths(res$rPstart,res$rMstart,minpos=st0,maxpos=ed0,Delta=1000,main=mainstr,cex.main=2,cex.lab=2)
        if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
        if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)   
        grid()
    }
    if(doHanging){
        cat(paste("Hanging Reads Count in ", 100," Windows",sep=""),"\n")
        mainstr=paste("Hanging Reads Count in ", 100," Windows",sep="")
        plotHangingPlusReads(res$rPstart,res$rMstart,minpos=st0,maxpos=ed0,Delta=1000,win=100,step=10,main=mainstr,cex.main=2,cex.lab=2)
        plotHangingMinusReads(res$rPstart,res$rMstart,minpos=st0,maxpos=ed0,Delta=1000,win=100,step=10,add=TRUE)
        if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
        if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)   
        grid()
        legend(x="topright", col=c("black","red"),pch=c(1,3),legend=c("Left-hang","Right-hang"),cex=2)
    }
    if(doSoftClip){
        mainstr=paste("Soft Clipped Reads Count")
        cat(mainstr,"\n")
        #cat("data points for hang plus:", length(res$rPstart), " ", length(res$rMstart),"\n")
        sc=getSoftClipClusters(pos=bam1$pos,cigar=bam1$cigar,seq=bam1$seq,strand=bam1$strand,minpos=st0,maxpos=ed0,MIN.READS.PER.CLUSTER=3)
        plot(sc$bcoord,sc$nclipreadsL,xlim=c(st0,ed0), ylim=c(0,max(c(sc$nclipreadsL,sc$nclipreadsR,1,na.rm=T))),xlab="Position",ylab="Number of Clippings",type="b",col="blue", lwd=1, main=mainstr,cex.main=2,cex.lab=2)
        lines(sc$bcoord,sc$nclipreadsR,lwd=1,col="red")
        #let's select top 10 and label it
        peak_scL=whichpart(sc$nclipreadsL)
        peak_scR=whichpart(sc$nclipreadsR)
        if(length(peak_scL)>0) text(sc$bcoord[peak_scL],sc$nclipreadsL[peak_scL],paste(sc$bcoord[peak_scL],sc$nclipreadsL[peak_scL],sep=":"),col="orange")
        if(length(peak_scR)>0) text(sc$bcoord[peak_scR],sc$nclipreadsR[peak_scR],paste(sc$bcoord[peak_scR],sc$nclipreadsR[peak_scR],sep=":"),col="green")
        if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
        if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)
        #grid()
        #legend(x="topright", col=c("skyblue","red"),lty=c(1,1),pch=c(1,1),legend=c("Left-clipped","Right-clipped"),cex=2)
    }
    bam1
}

whichpart <- function(x, n=10) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}


#cat("loaded line 5000... \n")
