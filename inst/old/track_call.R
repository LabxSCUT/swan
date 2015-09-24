#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(IRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(Biostrings))
suppressMessages(library(digest))
suppressMessages(library(robustbase))
suppressMessages(library(swan))
suppressMessages(library(Biobase))
suppressMessages(library(plyr))

option_list <- list(
  #complex calling control options
  make_option(c("-s", "--scores"), default="0lCBS",
    help="you can use predefined scores rule (manual) or customize [default %default]\n  e.g. lCBS:lCd,lDl+lDr"),
  make_option(c("-m", "--method"), default="0boot",
    help="you can use predefined method rule (manual) or customize [default %default]\n e.g. lCd=theo,lDr=boot:lDl=theo..."),
  make_option(c("-t", "--threshold"), default="0level3",
    help="you can use predefined thresh rule (manual) or customize [default %default]\n  e.g. lCd=level1:lCi=1..."),
  make_option(c("-f", "--filemerge"), default="0all",
    help="you can use predefined merge rule (manual)  or customize [default %default]\n  e.g. +1-2"),
  make_option(c("-n", "--namespace"), default="lCBS,lCd,lDl,lDr",
    help="you can use default namespace or customize (manual) [default %default]"),
  #merging options
  make_option(c("-x", "--gapsize"), default=20,
    help="maximum bps gap to merge windows into one SV [default %default]"),
  make_option(c("-y", "--support"), default=100,
    help="minimum bps support for windows considerred a SV [default %default]"),
  make_option(c("-z", "--bignum"), default=10, #this only merges at file level separate output
    help="minimum read# support for big isize considerred a SV [default %default]"),
  #other options
  make_option(c("-b", "--bam"), default="input",
    help="bam prefix of prefix.bam files [default %default]"),
  make_option(c("-d", "--sample"), default="DUMMY,INFO,MIX,DESCRIPTION",
    help="sample name [default %default]"),
  make_option(c("-p", "--species"), default="human",
    help="species name [default %default]"),
  #computing options
  make_option(c("-r", "--nonverified"), action="store_true", default=FALSE,
    help="if include regions not verified by softclipping analysis [default %default]"),
  make_option(c("-e", "--extra"), action="store_true", default=FALSE,
    help="if include regions considerred extra by softclipping analysis [default %default]"),
  make_option(c("-o", "--out"), default="input",
    help="prefix for output file [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
    help="save verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
    help="save debug, [default %default]")
)
parser <- OptionParser(usage = "%prog [options] refFile scoreFiles (implicitly assuming parFiles and bigFiles)", option_list=option_list)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
if(length(cmd$args)!=2){ print_help(parser); quit(); }

ptm=proc.time()
ref_file=cmd$args[1]
inputs=cmd$args[2]
text_methods=cmd$options$method
text_thresh=cmd$options$threshold
text_scores=cmd$options$scores
text_merge=cmd$options$filemerge
text_out=cmd$options$out
text_bam=cmd$options$bam
text_sample=cmd$options$sample
verbose=cmd$options$noQuiet #true to turn on
namespace=sort(strsplit(cmd$options$namespace,split=",")[[1]]) # sorted from lower to higher
adjust=cmd$options$adjust
bignum=cmd$options$bignum
support=cmd$options$support
gapsize=cmd$options$gapsize
species=cmd$options$species
nonverified=cmd$options$nonverified
extra=cmd$options$extra

#Thresholds from Nancy
MIN.READS.PER.SC.CLUSTER=5 # Minimum reads in soft clipping cluster.
SC.FROM.BK.THRESH=1000   # Maximum distance of soft clipping to breakpoint for confirmation.
SC.PROPORTION.MISMATCH.THRESH = 0.01
SCLEN.THRESH = 30  # Minium length of soft clipped sequence for remapping.
gainRelCNTHRESH = 1.2
lossRelCNTHRESH = 0.8
MIN.COV.SOFTCLIP.OVERRIDE = 100
MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE = 2
NUM.READS.FOR.CI = 100
CI.LEN = 500
REMAP.RELCN.RATIO.THRESH = 2
PRECISE.BKPOINT.OFFSET.ALLOWANCE = 20 # Should be some value that fits well within a read.
BIC.OVERSHOOT=10

score_files=strsplit(inputs,split=':')[[1]]
par_files=gsub("txt.gz","par.txt",score_files)
big_files=gsub("txt.gz","big.txt",score_files)
if(text_bam=="input"){ #one score to one bam correspondence
  bam_files=gsub("\\w+.txt.gz","bam",score_files)
} else { #one score to multiple bam corespondence
  sample_tags=strsplit(text_bam,split=':')[[1]]
  bam_files=strsplit(sample_tags,split=",")
  n_bams=lapply(bam_files,length) #number of bams in each group
}
n_files=length(score_files)
sample_id=strsplit(text_sample,split=",")[[1]][1]
sample_info=matrix(strsplit(text_sample,split=",")[[1]], ncol=4)

if(verbose) {
  cat("par_files:\n"); print(par_files)
  cat("big_files:\n"); print(big_files)
  cat("score_files:\n"); print(score_files)
  cat("bam_files:\n"); print(bam_files)
}

#let's start with all boot and replace any if specified
file_method=as.data.frame(matrix(rep("boot",length(namespace)*n_files),nrow=n_files))
if(!grepl("^[0-9]",text_methods)){
  method_text=strsplit(text_methods,split=':')[[1]]
  for(i in seq_len(method_text)){ 
    texts=strsplit(method_text[i],split=',')[[1]]
    for(t in texts){
      score=strsplit(t,split="=")[[1]][1]
      method=strsplit(t,split="=")[[1]][2]
      file_method[[score]][i]=method
    }
  }
} else {
  switch(text_methods,
    "0boot"={ file_method=as.data.frame(matrix(rep("boot",length(namespace)*n_files),nrow=n_files)) },
    "1theo"={ file_method=as.data.frame(matrix(rep("theo",length(namespace)*n_files),nrow=n_files)) },
    "2robust"={ file_method=as.data.frame(matrix(rep("robust",length(namespace)*n_files),nrow=n_files)) }
  )
}
colnames(file_method)=namespace
cat("theshold method\n"); print(file_method)

thresh_level=as.data.frame(matrix(rep(TRUE,length(namespace)*n_files),nrow=n_files))
thresh_value=as.data.frame(matrix(rep(3,length(namespace)*n_files),nrow=n_files))
colnames(thresh_level)=namespace
colnames(thresh_value)=namespace
if(!grepl("^[0-9]",text_thresh)){
  thresh_text=strsplit(text_thresh,split=':')[[1]]
  for(i in seq_len(thresh_text)){
    texts=strsplit(thresh_text[i],split=',')[[1]]
    for(t in texts){
      score=strsplit(t,split="=")[[1]][1]
      thresh=strsplit(t,split="=")[[1]][2]
      if(!grepl("level",thresh)){ # numeric
        thresh_value[[score]][i]=as.numeric(thresh)
        thresh_level[[score]][i]=FALSE
      } else {
        thresh_value[[socre]][i]=as.integer(substr(thresh,length(method)-1,1))
      }
    }
  }
} else {
  switch(text_thresh,
    "3level3"={ 
      thresh_level=as.data.frame(matrix(rep(TRUE,length(namespace)*n_files),nrow=n_files))
      thresh_value=as.data.frame(matrix(rep(3,length(namespace)*n_files),nrow=n_files))
    },
    "4level4"={ 
      thresh_level=as.data.frame(matrix(rep(TRUE,length(namespace)*n_files),nrow=n_files))
      thresh_value=as.data.frame(matrix(rep(4,length(namespace)*n_files),nrow=n_files))
    },
    "5level5"={ 
      thresh_level=as.data.frame(matrix(rep(TRUE,length(namespace)*n_files),nrow=n_files))
      thresh_value=as.data.frame(matrix(rep(5,length(namespace)*n_files),nrow=n_files))
    }
  )
}
cat("do we threshhold level?\n"); print(thresh_level); cat("what is the threshhold value\n"); print(thresh_value)

#allowed rule specification: lCBS:lCd,lDl+lDr:lCd,lDl+lDr:lCBS:lCd,lDl+lDr:lCd,lDl+lDr
file_rules=NULL
if(!grepl("^[0-9]",text_scores)){
  score_text=strsplit(text_scores,split=':')[[1]]
  file_scores=sapply(score_text,strsplit,split=",")
  file_rules=lapply(file_scores,sapply,strsplit,split="+",fixed=T)
} else {
  switch(text_scores,
    "0lCBS"={ file_rules=list(rep(n_files,"lCBS")) },
    "1lCd"={ file_rules=list(rep(n_files,"lCBS")) },
  )
}
file_tracks=lapply(lapply(lapply(file_rules,unlist),as.vector),unique)
cat("score file rules\n"); print(file_rules); cat("score file tracks\n"); print(file_tracks)
  
#allowed file merges: +1+2+3-4-5-6
merge_prules=NULL; merge_nrules=NULL
if(!grepl("^[0-9]",text_merge)){
  ppattern="\\+(?<file>[0-9]+)"
  npattern="-(?<file>[0-9]+)"
  merge_prules=as.integer(str_match_all_perl(text_merge,ppattern)[[1]][,2])
  merge_nrules=as.integer(str_match_all_perl(text_merge,npattern)[[1]][,2])
} else if (grepl("^0",text_merge)){
  merge_prules=seq(1:length(score_files))
  merge_nrules=NULL
} else {
  spattern="(?<pos>[0-9]+)split"
  split_pos=as.integer(str_match_all_perl(text_merge,spattern)[[1]][,2])
  stopifnot(split_pos<=n_files)
  merge_nrules=seq(1,split_pos-1)
  merge_prules=seq(split_pos,n_files)
}
if(length(merge_prules)==0) merge_prules=NULL; if(length(merge_nrules)==0) merge_nrules=NULL; 
cat("positive merging rules\n"); print(merge_prules); cat("negative merging rules\n"); print(merge_nrules)

#easy global options
if(text_out=="input"){
  out_prefix=gsub("\\.\\w*$","",lcPrefix(score_files),perl=T)
  out_prefix=gsub(".scores.txt.gz","",out_prefix)
} else {
  out_prefix=text_out
}
par_lst=list(); big_lst=list(); score_lst=list(); comment_lst=list(); 
stepsize=NULL; seqname=NULL; scan_start=NULL; scan_end=NULL;

### Load Scores ###
for(si in seq_along(par_files)){
  parFile=par_files[[si]]
  if(verbose) { cat("=Info:", si,"th par_file=", parFile, "\n") }
  par_lst[[si]]=read.table(parFile,header=T,sep="\t",colClass=NA,comment.char="#")
  if(verbose) { cat("=Info:", si,"th scan_par="); print(par_lst[[si]]) }
  bigFile=big_files[[si]]; big_lst[[si]]=data.frame()
  if(file.exists(bigFile))
    big_lst[[si]]=try(read.table(bigFile, header=F, sep="\t"),silent=T)
  scoreFile=score_files[[si]]
  list[comment,score_lst[[si]]]=read_com(scoreFile,comment_char="#",colClasses="numeric")
  if(class(big_lst[[si]])=='try-error') big_lst[[si]]=NULL
  if(is.null(stepsize)) stepsize=par_lst[[si]]$stepsize else stopifnot(stepsize==par_lst[[si]]$stepsize)
  if(is.null(seqname)) seqname=par_lst[[si]]$chr else stopifnot(seqname==par_lst[[si]]$chr)
  if(is.null(scan_start)) scan_start=par_lst[[si]]$start else stopifnot(scan_start==par_lst[[si]]$start)
  if(is.null(scan_end)) scan_end=par_lst[[si]]$end else stopifnot(scan_end==par_lst[[si]]$end)
}
#n_wins=ceiling(scan_end-scan_start+1)/stepsize;
win_pos=seq(scan_start,scan_end,stepsize)
n_wins=length(win_pos)
supnum=support/stepsize;gaplimit=gapsize/stepsize
ref=FaFile(ref_file)
rg=scanFa(ref, param=scanFaIndex(ref)) #this returns a DNAStringSet
ref_len=length(rg[[seqname]]) #[]->DNAStringSet, support width; [[]]->DNAString, support length
contig_info=matrix(c(seqname, ref_len, digest(toString(rg[[seqname]],algo="md5")), species), ncol=4)

### Compute Threshold ###
final_value=thresh_value #thresh value is included
tmp1=NULL; tmp2=NULL
for(si in seq_along(score_files)){
  for(ni in seq_along(namespace)){
    if(namespace[ni] %in% file_tracks[[si]]){
      if(namespace[ni]=="lCBS") {
        final_value[si,ni]=0.1 #lCBS is a dummy track with non-signal region all in zero
        next
      }
      if(thresh_level[si,ni]){ #use level
        if(file_method[si,ni]=="boot"){
          alpha=alpha_level[final_value[si,ni]]
          list[final_value[si,ni],tmp1,tmp2]=
            noovlap_est(data=score_lst[[si]][[namespace[ni]]],resample=1000,
                    buffer=round(par_lst[[si]]$delta*2/stepsize),alpha=alpha)
        } else if(file_method=="robust") {
          list[final_value[si,ni],tmp1,tmp2]=
            robust_est(data=scores_lst[[si]][[namespace[ni]]],lvl=final_value[si,ni])
        } else if(file_method=="theo") {
          list[final_value[si,ni],tmp1,tmp2]=
            theo_thresh(score=namespace[ni],alpha=alpha,scan_par=par_lst[[si]],verbose=F)
        } else {
          stop(paste("method",file_method[si,ni],"implemented"))
        }
      } #else use final_value, we don't need to do anything
    }
  }
} #all values in final_value updated
print(final_value)

### Call Tracks and Track Merge to File Calls ###
file_calls=matrix(rep(FALSE,n_wins*n_files),ncol=n_files)
if(verbose) { print(n_files); print(dim(file_calls)); print(dim(score_lst[[si]])) }
for(si in seq_along(score_files)){
  file_codes=get_rule_code(file_rules[[si]],NULL,namespace)
  tmp_tracks=matrix(rep(FALSE,n_wins*length(namespace)),ncol=length(namespace))
  dim(tmp_tracks)
  for(ni in seq_along(namespace)){
    if(namespace[ni] %in% file_tracks[[si]]) {
      tmp_tracks[,ni]=score_lst[[si]][namespace[ni]]>=final_value[si,ni]
    }
  }
  tmp_tracks=apply(tmp_tracks,MARGIN=1,FUN=jumpy.roo)
  file_calls[,si]=code2call(tmp_tracks,file_codes)
}
merge_codes=get_rule_code(merge_prules,merge_nrules,seq_len(n_files))
merge_track=apply(file_calls,MARGIN=1,FUN=jumpy.roo)
merge_calls=code2call(merge_track,merge_codes)

filter_calls=function(merge_calls, supnum=supnum, gaplimit=gaplimit) {
  merge_rle=rle(merge_calls)
  merge_rle$values[merge_rle$lengths<=gaplimit&!merge_rle$values]=TRUE #extend through small gaps
  merge_rle$values[merge_rle$lengths<=supnum&merge_rle$values]=FALSE #filter out small bumps
  inverse.rle(merge_rle)
}
merge_calls=filter_calls(merge_calls, supnum=supnum, gaplimit=gaplimit)
rle_merge_calls=Rle(merge_calls)
merge_regions=IRanges(start=scan_start+(start(rle_merge_calls)[IRanges::runValue(rle_merge_calls)]-1)*stepsize,
                      end=scan_start+(end(rle_merge_calls)[IRanges::runValue(rle_merge_calls)]-1)*stepsize)
n_sv=length(merge_regions)

offsetL =rep(NA, n_sv)  # When sc clusters are found, how far are they from the left seqCBS breakpoint?
offsetR = rep(NA,n_sv)  # When sc clusters are found, how far are they from the right seqCBS breakpoint?
toss = rep(FALSE,n_sv) # Regions that are not "validated" by downstream soft clipping will be tossed. 
scPASS = rep(FALSE,n_sv)
seqremapPASS = rep(FALSE,n_sv)

#following adapted from Nancy's dream_seqCBS_pipeline.r
#print(merge_regions)
pair_mode=TRUE
if(length(bam_files)==1) pair_mode=FALSE
if(pair_mode){
  ctr_files=bam_files[[1]]; case_files=bam_files[[2]]
} else {
  case_files=bam_files[[1]]
}
target_set=data.frame(); extra_set=data.frame()

for(rind in seq_len(n_sv)){
  st=start(merge_regions[rind]); ed=end(merge_regions[rind]); 
  sv_target=list(CHROM=seqname,POS=st,END=ed,ID=NA,SVTYPE=NA,FILTER=NA,IMPRECISE=TRUE,CIPOSL=NA,CIPOSH=NA,CIENDL=NA,CIENDH=NA)
  sv_extra=list()
  w=ed-st;w=min(max(w, 10000),50000)
  st0=st-2*w        # start of analysis region 
  ed0=ed+2*w        # end of analysis region
  what=c("qname","rname", "flag","strand", "pos","mapq","mrnm","mpos","isize","qwidth","seq","cigar","qual") 
  seq_info=GRanges(seqname, IRanges(st0, ed0))
  case_reads=DataFrame(); 
  for(bam in case_files){
    case_reads=rbind(case_reads,DataFrame(allFunction(seq_info, bam, what=what, index=bam)))
  } #we combine gi_reads from multiple bamfiles  
  if(pair_mode) { ctr_reads=DataFrame();
    for(bam in ctr_files){
      ctr_reads=rbind(ctr_reads,DataFrame(allFunction(seq_info, bam, what=what, index=bam)))
    } #we combine gi_reads from multiple bamfiles
    pos0r=IRanges(start=ctr_reads$pos[!is.na(ctr_reads$pos)],width=ctr_reads$qwidth[!is.na(ctr_reads$pos)]) 
  }
  fl=parseFlags(case_reads$flag)
  temp=unique(case_reads$qname)
  mm = match(temp,case_reads$qname)
  sel = which(case_reads$strand[mm]=="+" & as.integer(case_reads$rname[mm])==seqname)
  rPstart1 = case_reads$pos[mm[sel]]
  rMstart1 = case_reads$mpos[mm[sel]]
  isDiscordant1 = fl$QStrand[mm[sel]] == fl$MStrand[mm[sel]]
  sel = which(case_reads$strand[mm]=="-" & as.integer(case_reads$rname[mm])==seqname)
  rMstart1 = c(rMstart1, case_reads$pos[mm[sel]])
  rPstart1 = c(rPstart1, case_reads$mpos[mm[sel]])
  isDiscordant1 = c(isDiscordant1,fl$QStrand[mm[sel]] == fl$MStrand[mm[sel]])
  pos1r=IRanges(start=case_reads$pos[!is.na(case_reads$pos)],width=case_reads$qwidth[!is.na(case_reads$pos)])
  
  # -------------- Soft-clipping --------------------- #
  sc=getSoftClipClusters(pos=case_reads$pos,cigar=case_reads$cigar,seq=case_reads$seq,strand=case_reads$strand,minpos=st0,maxpos=ed0,MIN.READS.PER.CLUSTER=MIN.READS.PER.SC.CLUSTER)
  ####  Which soft-clip cluster is next to the break point?
  scL=summarizeClusters(sc$LsoftclipClusters)
  scR=summarizeClusters(sc$RsoftclipClusters)
  
  # ------------ Is this duplication, deletion, or inversion? --------------#
  # First, look for soft clip clusters near st, ed.
  # If soft-clip clusters found at both end points, then register this as true and precise.
  # If double-double clusters are found, then this is an inversion or transposition.
  if(length(sc$LsoftclipClusters)>0 && length(sc$RsoftclipClusters)>0){
    dd = intersect(scL$Position,scR$Position)
    if(length(dd)==2 && max(scL$ProportionMismatch[match(dd,scL$Position)])<SC.PROPORTION.MISMATCH.THRESH && max(scR$ProportionMismatch[match(dd,scR$Position)])<SC.PROPORTION.MISMATCH.THRESH){
      if(countOverlaps(IRanges(start=min(dd),end=max(dd)),IRanges(start=st,end=ed))>0){
        # This TRP/INV overlaps with the seqCBS detected segment.
        sv_target$POS=dd[which.min(abs(dd-st))]; sv_target$END=dd[which.min(abs(dd-ed))]
        sv_target$SVTYPE="TRPINV"; sv_target$ID=paste("TRPINV_",seqname,"_",rind,sep="")
        sv_target$IMPRECISE=FALSE; sv_target$FILTER="PASS";
        offsetL[rind] = abs(sv_target$START-st)     #
        offsetR[rind] = abs(sv_target$END-ed)       #
        #cat("\t\tFound Double Double that overlaps with original.\n",sep="")                         
      } else {
        # No overlap, this should be entered as an "extra".
        sv_extra=list(CHROM=seqname,POS=min(dd),END=max(dd),ID=".",SVTYPE="TRPINV",FILTER=".",IMPRECISE=F,CIPOSL=NA,CIPOSH=NA,CIENDL=NA,CIENDH=NA)
        #cat("\t\tTRPINV -- Registered breakpoint at chr",chrnum,", start=",min(dd)," that does not overlap with original.\n",sep="")                         
      }
      ### TODO: What if length(dd)>2?
      ### TODO: Distinguish what this is: transposition or inversion?  Where did it go?
      ### TODO: If there are more than 2 double clusters, then record those as well. 
    } 
  }
  target_set=rbind(target_set,as.data.frame(sv_target,stringsAsFactor=F))
  
  # Check relative copy number, and assign to DUP or DEL states if not already an INV state.    
  cov1 = countOverlaps(IRanges(start=st,end=ed),pos1r)
  if(pair_mode) cov0 = countOverlaps(IRanges(start=st,end=ed),pos0r) else cov0=par_lst[[1]]$coverage
  # in non-pair mode then compre with the its own average coverage, 
  # or we shall always compare to own coverage
  relCN = (cov1/length(st:ed))/(cov0/length(st:ed))
  if(relCN>1.2 && is.na(target_set$SVTYPE[rind])){
    target_set$SVTYPE[rind]="DUP"
  } else if(relCN<0.8 && is.na(target_set$SVTYPE[rind])){
    target_set$SVTYPE[rind]="DEL"
  } else if(is.na(target_set$SVTYPE[rind])){
    target_set$SVTYPE[rind]="INS"
  } 
  cat("---> This region is a ",target_set$SVTYPE[rind],".\n",sep="")
  
  # ------------ Process deletions  --------------#
  if(target_set$SVTYPE[rind]=="DEL"){
    target_set$ID[rind]=paste("DEL_",seqname,"_",rind,sep="")
    
    # Putative deletion:  Look for left-hanging-read-cluster near ed, and right-hanging-read-cluster near st.
    # If clusters found, and clipped sequences mapped successfully, then mark breakpoints as precise.  
    # Otherwise, if coverage is larger than MIN.COV.SOFTCLIP.OVERRIDE and relCN larger than MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE
    # then mark the breakpoints as imprecise.
    # If neither of the above holds, then consider this region as false positive.
    
    if(length(sc$LsoftclipClusters)>0 && length(sc$RsoftclipClusters)>0){
      Led = which.min(abs(scL$Position-ed))
      offsetR[rind] = abs(scL$Position[Led]-ed)
      Rst = which.min(abs(scR$Position-st))
      offsetL[rind] = abs(scR$Position[Rst]-st)
      if(scL$ProportionMismatch[Led]<SC.PROPORTION.MISMATCH.THRESH && scR$ProportionMismatch[Rst]<SC.PROPORTION.MISMATCH.THRESH && offsetR[rind]<SC.FROM.BK.THRESH && offsetL[rind]<SC.FROM.BK.THRESH ){ 
        scPASS[rind]=TRUE
      }
      if(scPASS[rind]){
        Lseq = sc$LsoftclipClusters[[Led]]$consmat$consensus
        Rseq = sc$RsoftclipClusters[[Rst]]$consmat$consensus
        cat("\tSoft-clip Clusters found near breakpoints.  Matching clipped sequences to genome...\n")
        localseq = getSeq(Hsapiens, names=paste("chr",chrnum,sep=""), start=st0, end=ed0)
        Lmatch=matchPattern(Lseq,localseq,min.mismatch=0,max.mismatch=1)  
        Rmatch=matchPattern(Rseq,localseq,min.mismatch=0,max.mismatch=1)  
        if(length(Lmatch)>0 && length(Rmatch)>0){
          Lbasecoord = start(Lmatch) + nchar(Lseq) + st0-1
          Rbasecoord = start(Rmatch) + st0-1
          cat("\t\tLeft sequence mapped to ",Lbasecoord,", Right sequence mapped to ",Rbasecoord,".\n",sep="")
          cat("\t\tscR Position=",scR$Position[Rst],", scL Position=",scL$Position[Led],".\n",sep="")
          if(Lbasecoord == scR$Position[Rst] && Rbasecoord==scL$Position[Led]){
            seqremapPASS[rind]=TRUE
            target_set$POS[rind]=Lbasecoord-1
            target_set$END[rind]=Rbasecoord-1
            target_set$IMPRECISE[rind]=FALSE
          }
        } else {
          cat("\t\tOne of the two ends didn't map: length(Lmatch)=",length(Lmatch),", length(Rmatch)=",length(Rmatch),".\n")
        }
      } else {
        cat("\tNo Soft-clip Clusters found.\n")
      }
    }
    if(!seqremapPASS[rind]){
      cat("\tDid not get confirmed by soft-clip remapping.\n")
      # Did not get confirmed by soft-clipped sequence remapping.
      if((cov0+cov1)< MIN.COV.SOFTCLIP.OVERRIDE || relCN> 1/MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE){
        cat("\tCoverage=",cov1,", relCN=",relCN, ", this is tossed.\n", sep="")
        toss[rind]=TRUE
      }
      target_set$POS[rind] = st
      target_set$END[rind]=ed
      target_set$CIPOSL[rind] = -CI.LEN
      target_set$CIPOSH[rind] = CI.LEN
      target_set$CIENDL[rind] = -CI.LEN
      target_set$CIENDH[rind] = CI.LEN
    }
  } # End of processing "DEL"
  
  # ------------ Process duplications  --------------#
  if(target_set$POS[rind]=="DUP"){
    target_set$ID[rind]=paste("DUP_",chrnum,"_",rind,sep="")
    
    # Putative duplication:  Look for left-hanging-read-cluster near st, and right-hanging-read-cluster near ed.
    # If clusters found, register this duplication as TRUE and precise.
    # If clipped sequences mapped successfully and left and right remap positions meet in a single genome position, then 
    # record a new insertion in that position.
    # If no clusters are found, but if coverage is larger than MIN.COV.SOFTCLIP.OVERRIDE and relCN larger than MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE
    # then mark the breakpoints as imprecise.
    # If neither of the above holds, then consider this region as false positive.
    
    if(length(sc$LsoftclipClusters)>0 && length(sc$RsoftclipClusters)>0){
      Lst = which.min(abs(scL$Position-st))
      offsetL[rind] = abs(scL$Position[Lst]-st)
      Red = which.min(abs(scR$Position-ed))
      offsetR[rind] = abs(scR$Position[Red]-ed)
      if(scL$ProportionMismatch[Lst]<SC.PROPORTION.MISMATCH.THRESH && scR$ProportionMismatch[Red]<SC.PROPORTION.MISMATCH.THRESH && offsetR[rind]<SC.FROM.BK.THRESH && offsetL[rind]<SC.FROM.BK.THRESH ){ 
        scPASS[rind]=TRUE
      }
      if(scPASS[rind]){
        target_set$IMPRECISE[rind]=FALSE
        target_set$POS[rind]=scL$Position[Lst]
        target_set$END[rind]=scR$Position[Red]
        Lseq = sc$LsoftclipClusters[[Lst]]$consmat$consensus
        Rseq = sc$RsoftclipClusters[[Red]]$consmat$consensus
        cat("\tSoft-clip Clusters found near breakpoints.")
        if(nchar(Lseq)>SCLEN.THRESH) {
          cat("\t\tMatching left clipped sequences to genome...\n")
          ptm=proc.time()
          Lmatch=vmatchPattern(Lseq,hg19template,min.mismatch=0,max.mismatch=1)  
          ptm2=proc.time(); elapsed = ptm2-ptm
          cat("\t\tLeft clip mapping took ",elapsed[3]," seconds.\n",sep="")  
          tempL= unlist(lapply(Lmatch, length))
        } else {
          cat("\tLeft clipped sequences too short to map.\n")
          tempL = rep(0,23)
        }
        if(nchar(Rseq)>SCLEN.THRESH){
          ptm2=proc.time(); 
          Rmatch=vmatchPattern(Rseq,hg19template,min.mismatch=0,max.mismatch=1)  
          elapsed=proc.time()-ptm2
          cat("\t\tRight clip mapping took ",elapsed[3]," seconds.\n",sep="")   
          tempR= unlist(lapply(Rmatch, length))
        } else {
          cat("\tRight clipped sequences too short to map.\n")
          tempR = rep(0,23)
        }
        
        cat("\t\t",sum(tempL),"matches found for left soft-clip sequence.\n",sep=" ")
        cat("\t\t",sum(tempR),"matches found for right soft-clip sequence.\n",sep=" ")
        
        if(abs(sum(tempL)/(max(1,relCN-1)))<REMAP.RELCN.RATIO.THRESH &&abs(sum(tempR)/(max(1,relCN-1)))<REMAP.RELCN.RATIO.THRESH ){
          # The number of times Lseq and Rseq are found in the genome is roughly of the same order as the
          # relative copy number, and so we register those found positions as true breakpoints.
          # Otherwise, we don't register those positions, for fear of a repetitive element
          # blowing up the number of false discoveries.
          seqremapPASS[rind]=TRUE
          for(ch in which(tempL>0)){
            sv_extra=as.data.frame(list(CHROM=rep(ch,length(Lmatch[[ch]])),POS=end(Lmatch[[ch]]),END=end(Lmatch[[ch]])+target_set$POS[rind]-target_set$END[rind],ID=rep(regionid[rind],length(Lmatch[[ch]])),SVTYPE=rep("INS",length(Lmatch[[ch]])),FILTER=".",IMPRECISE=rep(FALSE,length(Lmatch[[ch]])),CIPOSL=NA,CIPOSH=NA,CIENDL=NA,CIENDH=NA))
            cat("\t\tGain -- Registered breakpoint(s) at chr",ch,", start=",end(Lmatch[[ch]]),".\n",sep="")                         
          }
          for(ch in which(tempR>0)){
            matchpos=start(Rmatch[[ch]])
            #sel = which(abs(matchpos-startpos2[chr2==ch]) > PRECISE.BKPOINT.OFFSET.ALLOWANCE)
            sel = 1:length(matchpos)  ### TODO: Have a better way of filtering out the repeats.
            if(length(sel)>0){
              sv_extra=rbind(sv_extra,as.data.frame(list(CHROM=rep(ch,length(sel)),POS=matchpos[sel],END=matchpos[sel] + target_set$POS[rind]-target_set$END[rind],ID=rep(regionid[rind],length(sel)),SVTYPE=rep("INS",length(sel)),FILTER=".",IMPRECISE=rep(FALSE,length(sel)),CIPOSL=NA,CIPOSH=NA,CIENDL=NA,CIENDH=NA)))
              cat("\t\tGain -- Registered breakpoint(s) at chr",ch,", start=",matchpos[sel],".\n",sep="")                         
            }
          }
        }
      }
    }
    if(!scPASS[rind]){
      # Did not get confirmed by soft-clipped sequence remapping.
      if((cov0+cov1)< MIN.COV.SOFTCLIP.OVERRIDE || relCN< MIN.FOLDCHANGE.SOFTCLIP.OVERRIDE){
        cat("\tCoverage=",cov1,", relCN=",relCN, ", this is tossed.\n", sep="")
        toss[rind]=TRUE 
      }
      cat("\tDid not get confirmed by soft-clip clusters.\n")
      target_set$POS[rind] = st
      target_set$END[rind] = ed
      target_set$CIPOSL[rind] = -CI.LEN
      target_set$CIPOSH[rind] = CI.LEN
      target_set$CIENDL[rind] = -CI.LEN
      target_set$CIENDH[rind] = CI.LEN      
    }   
  }  ### End of if(alt[rind]=="DUP")
} ### End of for(rind in 1:nregions)

toss=if(nonverified) NULL else toss
extra_set=if(extra) extra_set else NULL
final_set=rbind(if(is.null(toss)) target_set else target_set[-which(toss)],extra_set)

write_set=function(final_set,scores_lst,file_tracks,sample_id){
  vcf_field=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  bed_field=c("CHROM","POS","END","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  vcf_class=c("character","numeric","character","character","character","numeric","character","character","character","character")
  bed_class=c("character","numeric","numeric","character","character","character","numeric","character","character","character","character")
  n_sv = nrow(final_set); if(verbose) cat("==Info: number of final calls=", n_sv,"\n")
  vcf_empty = matrix(rep(rep("",length(vcf_field)+1),n_sv), ncol=length(vcf_field)+1)
  bed_empty = matrix(rep(rep("",length(bed_field)+1),n_sv), ncol=length(bed_field)+1)
  vcf_set = data.frame(vcf_empty, stringsAsFactors=F)
  bed_set = data.frame(bed_empty, stringsAsFactors=F)
  colnames(vcf_set) = c(vcf_field,sample_id)
  colnames(bed_set) = c(bed_field,sample_id)
  for(vi in seq_len(n_sv)){
    #if break points is not found or scarry
    sv_imprecise=if(final_set$IMPRECISE[vi]) "IMPRECISE:" else ""
    sv_filter=if(is.na(final_set$FILTER[vi])) "." else final_set$FILTER[vi]
    sv_type=paste("<",final_set$SVTYPE[vi],">",sep="")
    sv_info=gettextf("%sSVTYPE=%s;END=%d;SVLEN=%d;CIPOS=%d,%d;CIEND=%d,%d",sv_imprecise,sv_type,final_set$END[vi],final_set$END[vi]-final_set$POS[vi],final_set$CIPOSL[vi],final_set$CIPOSH[vi],final_set$CIENDL[vi],final_set$CIENDH[vi])
    ref_base=toString(subseq(rg[[seqname]],final_set$POS[vi],final_set$POS[vi]))
    sv_win_start=floor((final_set$POS[vi]-scan_start)/stepsize)+1
    sv_win_end=floor((final_set$END[vi]-scan_start)/stepsize)+1
    sel=sv_win_start:sv_win_end
    sv_format=paste(unlist(lapply(seq_along(file_tracks),function(i,file_tracks){ paste(i,file_tracks[[i]],sep="_") },file_tracks)),collapse=":")
    sv_vals=lapply(seq_along(file_tracks),function(i,file_tracks,score_lst,sel){ apply(score_lst[[i]][sel,][file_tracks[[i]]],2,max) },file_tracks,score_lst)
    sv_qual=max(unlist(sv_vals[merge_prules]))
    sv_dummy=paste(round(unlist(sv_vals)),collapse=":")
    vcf_set[vi,]=c(seqname,final_set$POS[vi],final_set$ID[vi],ref_base,sv_type,sv_qual,sv_filter,sv_info,sv_format,sv_dummy)
    bed_set[vi,]=c(seqname,final_set$POS[vi],final_set$END[vi],final_set$ID[vi],ref_base,sv_type,sv_qual,sv_filter,sv_info,sv_format,sv_dummy)
  }
  vcf_meta=paste(call_meta(
    c(meta_keys, sv_keys), c(meta_values, sv_values), ref_file, "bed2vcf.R"),
                 sample_meta(sample_info[,1], sample_info[,2], sample_info[,3], sample_info[,4]),
                 contig_meta(contig_info[,1], contig_info[,2], contig_info[,3], contig_info[,4]),
                 gettextf("#%s", paste(colnames(vcf_set), collapse='\t')), sep='\n')
  vcf_set=set_colClass(vcf_set,vcf_class)
  bed_set=set_colClass(bed_set,bed_class)
  #print(vcf_set)
  #print(bed_set)
  return(list(vcf_set=vcf_set,vcf_meta=vcf_meta,bed_set=bed_set))
}

list[vcf_set,vcf_meta,bed_set]=write_set(final_set,scores_lst,file_tracks,sample_id)
sv_order=sort(vcf_set$POS,index.return=T)$ix
vcf_set=vcf_set[sv_order,]
bed_set=bed_set[sv_order,]
if(verbose) cat("==Info: vcf sorted\n")
write_com(vcf_set,comment=vcf_meta,filename=paste(out_prefix,".vcf",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
write_com(bed_set,comment=NULL,filename=paste(out_prefix,".bed",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)

if(verbose) cat("=Info: time=", proc.time()-ptm, "\n")
if(verbose) cat("-call safely done!\n")

#debug variable lines; use two files, use 4 types of track, use bootstrap, use alpha=0.01
#use either lCd, lCBS, or lDl+lDr, find SVs in file 1 but not in file 2
# setwd("~/Downloads/swan/real")
# inputs="synthetic.challenge.set1.normal.v2.16M-18M.scores.txt.gz:synthetic.challenge.set1.tumor.v2.16M-18M.scores.txt.gz"
# text_bam="synthetic.challenge.set1.normal.v2.16M-18M.bam:synthetic.challenge.set1.tumor.v2.16M-18M.bam"
# namespace=c("lCBS","lCd","lDl","lDr"); text_methods="1theo"; text_thresh="5level5"; text_scores="lCBS:lCBS"
# text_merge="2split"; text_out="input"; verbose=T; adjust=0; bignum=3; support=100; gapsize=20; extra=T; nonverified=T
# text_sample="DUMMY,INFO,MIX,DESCRIPTION"
# vcf_field=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
# bed_field=c("CHROM","POS","END","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
# vcf_class=c("character","numeric","character","character","character","numeric","character","character","character","character")
# bed_class=c("character","numeric","numeric","character","character","character","numeric","character","character","character","character")
# ref_file="human_g1k_v37.fasta"
#sv_field=c("CHROM","POS","END","ID","ALT","QUAL","IMPRECISE","CIPOS","CIEND","SVTYPE")
#sv_class=c("character","numeric","numeric","character","character","character","logical","character","character","character")
# par(mfrow=c(4,1))
# par(mar=rep(2,4))
# plot(file_calls[,1])
# plot(file_calls[,2])
# plot(file_calls[,2]&(!file_calls[,1]))
# plot(merge_calls)
