#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"

suppressMessages(library(optparse))

#suppressMessages(library(optparse))
#suppressMessages(library(IRanges))
#suppressMessages(library(Rsamtools))
#suppressMessages(library(Biostrings))
#suppressMessages(library(digest))
#suppressMessages(library(robustbase))
#suppressMessages(library(Biobase))
#suppressMessages(library(plyr))

#	tryCatch({
#		source("~/setup/swan/R/libswan.R")
#  	source("~/setup/swan/R/svlib.r")
#  	source("~/setup/swan/R/ThresholdFunctions.R")
#  	dyn.load("~/setup/swan/src/swan.so")
#		cat("-Info: loaded from ~/setup/swan package!\n")
#  }, error=function(e) {
#		cat("-Warn: cant load from ~/setup/swan package!\n")
#		suppressMessages(library(swan))
#		cat("-Info: loaded from installed swan package!\n")
#  }) 

#expected inputs
#spX is sample; spY is control
#if mono: expect olny swan_score and swan_opt input of spX; 
#else: expect two swan_scores and optional seqcbs and sclip input, first swan_score is control
#swan_score files is expected in [spY.swan.txt.gz:]spX.swan.txt.gz, [spY.swan.par.txt:]spX.swan.par.txt
#swan_bigd files is expected in [spY.swan.bigd.txt:]spX.swan.bigd.txt 
#seqcbs files is expected in spX.seqcbs.txt, spX.seqcbs.par.txt
#sclip files is expected in spX.sclip.RData, spX.sclip.par.txt
#bamfiles is expected in [ctrY.rg1.bam,ctrY.rg2.bam,...:]spX.rg1.bam,spX.rg2.bam,...
#parameter filenames and bai filenames are assumed and cann't be specified
#differential opt input: [spY_opt]:spX_opt
#tier1 separator (:); tier2 separator (_); tier3 separator (,)

#default1 for single input, default2 for matched input
#swan_opt_default1="track=lCd,method=empr,thresh=level7,sup=100,gap=100_track=lDr+lDl,method=theo,thresh=level3,sup=100,gap=100"
swan_opt_default1="track=lCd,method=empr,thresh=level7,sup=100,gap=100_track=lDr+lDl,method=theo,thresh=level3,sup=100,gap=100_track=ins,sup=50,cvg=5_track=del,sup=50,cvg=5"
swan_opt_default2="track=lCd,method=empr,thresh=level6,sup=50,gap=100_track=lDr+lDl,method=theo,thresh=level2,sup=50,gap=100_track=ins,sup=20,cvg=2_track=del,sup=20,cvg=2:track=lCd,method=empr,thresh=level7,sup=100,gap=100_track=lDr+lDl,method=theo,thresh=level3,sup=100,gap=100_track=ins,sup=50,cvg=5_track=del,sup=50,cvg=5"
bigd_opt_default1="minmpr=5,maxins=50000"
bigd_opt_default2="minmpr=2,maxins=50000:minmpr=5,maxins=50000"
disc_opt_default1="minmpr=5,maxins=10000"
disc_opt_default2="minmpr=2,maxins=20000:minmpr=5,maxins=20000"
seqcbs_opt_default="minstat=0,sup=1500,gap=1000,expand=2000,good=4"
sclip_opt_default=""

option_list <- list(
  make_option(c("-c", "--chrname"), default="22",
    help="chromosome name, [default: %default] \n"),
  make_option(c("-t", "--stat"), default="",
    help="stat inputs: [spY.stat:]spX.stat; 
          [spY.stat:]spX.stat implicitly assumed \n"),
  make_option(c("-i", "--swan"), default="",
    help="swan inputs: [spY.swan.txt.gz:]spX.swan.txt.gz; 
          [spY.swan.par.txt:]spX.swan.par.txt implicitly assumed \n"),
  make_option(c("-j", "--bigd"), default="",
    help="big deletion inputs: [spY.bigd.txt:]spX.bigd.txt; 
          [spY.swan.par.txt:]spX.swan.par.txt implicitly assumed \n"),
  make_option(c("-k", "--seqcbs"), default="",
    help="seqcbs inputs: spX.seqcbs.txt; 
          spX.seqcbs.par.txt implicitly assumed \n"),
  make_option(c("-l", "--sclip"), default="",
    help="sclip inputs: spX.sclip.Rdata; 
          spX.sclip.par.txt implicitly assumed \n"),
  make_option(c("-m", "--disc"), default="",
    help="discordant cluster inputs: [spY.disc.txt:]spX.disc.txt; 
          [spY.swan.par.txt:]spX.swan.par.txt implicitly assumed \n"),
  make_option(c("-u", "--swan_opt"), default="",
    help=paste("swan options: [spY_opt:]track=t1_key1=value1_key2=value2,track=t2_..., \n",
							 "default1:", swan_opt_default1, "\n",
							 "default2:", swan_opt_default2, "\n"
							)),
  make_option(c("-v", "--bigd_opt"), default="",
    help=paste("swan big deletion options: [spY_opt:]key1=value1,key2=value2,..., \n",
							 "default1:", bigd_opt_default1, "\n",
							 "default2:", bigd_opt_default2, "\n"
							 )), 
  make_option(c("-w", "--seqcbs_opt"), default=seqcbs_opt_default,
    help=paste("seqcbs options: key1=value1,key2=value2,..., default:", seqcbs_opt_default, "\n")), 
  make_option(c("-x", "--sclip_opt"), default=sclip_opt_default,
    help=paste("sclip inputs: key1=value1,key2=value2,..., default:", sclip_opt_default, "\n")),  
  make_option(c("-y", "--disc_opt"), default="",
    help=paste("swan discordent clusters options: [spY_opt:]key1=value1,key2=value2,..., \n", 
							 "default1:", disc_opt_default1, "\n",
							 "default2:", disc_opt_default2, "\n"
							 )), 
  make_option(c("-d", "--override"), default="",
    help="bed formatted with colnames, parameter overriding files for swan calling: 
          [spY.swan.ovrd.txt:]spX.swan.ovrd.txt \n"),
  make_option(c("-f", "--fineconf"), action="store_true", default=FALSE,
    help="fine conf mode and .bam is assumed for all inputs, see manual [default %default]"),
  make_option(c("-o", "--outprefix"), default="input",
    help="prefix for output file [default %default]"),
  make_option(c("-p", "--sample"), default="spX,INFO,MIX,DESCRIPTION",
    help="mannual override of spX information [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
    help="verbose mode and additional information outputs [default %default]"),
  make_option(c("-r", "--confirm"), default="dedup",
    help="use which confirmation mode? combinaitons of dedup,nostrict [default %default]"),
  make_option(c("-s", "--savevcf"), default="",
    help="whether to savevcf file (slower) and parameters, e.g. 
					species=human_sapien:other_opt=other_value
					[default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
    help="debug mode and additional .RData is assumed for all inputs, see manual [default %default]")
)
parser <- OptionParser( 
  usage="%prog [options] refFile [spY.rg1.bam,spY.rg2.bam,...:]spX.rg1.bam,spX.rg2.bam,...", 
  option_list=option_list)
cat("-Info: swan_join vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
if(length(cmd$args)!=2){ print_help(parser); quit(); }

#starting the real work after checking format
suppressMessages(library(swan))
suppressMessages(for(p in c("digest","plyr","robustbase","optparse")) myrequire(p))
suppressMessages(for(p in c("IRanges","Rsamtools","Biostrings","Biobase")) myrequire(p,repo="Bioc"))
ref_file=cmd$args[1]
text_bam=cmd$args[2]
opt_verbose=cmd$options$noQuiet #set to turn on
opt_debug=cmd$options$debug #set to turn on
verbose=opt_verbose
debug=opt_debug
if(debug){ options(warn=2) } else { options(warn=-1) }
gtk=proc.time()[3]
gmk=get_gmk(Sys.getpid())

cat("==Info: assigning options...\n")
text_stat=cmd$options$stat
text_chrname=cmd$options$chrname
text_swan=cmd$options$swan
text_bigd=cmd$options$bigd
text_disc=cmd$options$disc
text_seqcbs=cmd$options$seqcbs
text_sclip=cmd$options$sclip
text_swan_opt=cmd$options$swan_opt
text_bigd_opt=cmd$options$bigd_opt
text_disc_opt=cmd$options$disc_opt
text_seqcbs_opt=cmd$options$seqcbs_opt
text_sclip_opt=cmd$options$sclip_opt
text_override=cmd$options$override
text_savevcf=cmd$options$savevcf
text_sample=cmd$options$sample
#text_species=cmd$options$species
text_outprefix=cmd$options$outprefix
opt_verbose=cmd$options$noQuiet #set to turn on
opt_debug=cmd$options$debug #set to turn on
opt_fineconf=cmd$options$fineconf #set to turn on
opt_confirm=cmd$options$confirm #set to turn on

#debug options
#setwd("~/Downloads/swan/test")
#ref_file="~/hg/hg19/human_g1k_v37.fasta"; text_bam="test_multi.lib1.bam,test_multi.lib2.bam,test_multi.lib3.bam:test_multi.sv.lib1.bam,test_multi.sv.lib2.bam,test_multi.sv.lib3.bam"; text_chrname="11"; text_swan="test_multi.swan.txt.gz:test_multi.sv.swan.txt.gz"; text_bigd="test_multi.swan.bigd.txt:test_multi.sv.swan.bigd.txt"; text_seqcbs="test_multi.sv.seqcbs.txt"; text_sclip="test_multi.sv.sclip.RData"; text_bigd_opt="minmpr=3"; text_seqcbs_opt="minstat=0"; text_sclip_opt=""; text_sample="test_multi.sv,INFO,MIX,DESCRIPTION"; text_species="human"; text_outprefix="input"; text_override="test_multi.swan.ovrd.txt:test_multi.sv.swan.ovrd.txt"; opt_verbose=TRUE; opt_debug=TRUE; opt_confirm=FALSE; text_swan_opt="track=lCd,method=theo,thresh=level3,sup=200,gap=200:track=lDr+lDl,method=theo,thresh=level3,tele=100,sup=100,gap=100:track=ins,sup=50,cvg=5:track=del,sup=50,cvg=5";
#ref_file="~/hg/hg19/human_g1k_v37.ucsc.fasta"; text_bam="NA12878.chr22.16M-21M.bam"; text_chrname="chr22"; text_swan="NA12878.chr22.16M-21M.swan.txt.gz"; text_bigd="NA12878.chr22.16M-21M.swan.bigd.txt"; text_seqcbs=""; text_sclip="";  text_bigd_opt="minmpr=5"; text_seqcbs_opt="minstat=0"; text_sclip_opt=""; text_sample="NA12878,INFO,MIX,DESCRIPTION"; text_species="human"; text_outprefix="input"; text_override="NA12878.chr22.16M-21M.swan.ovrd.txt"; opt_verbose=TRUE; opt_debug=TRUE; opt_confirm=FALSE; text_swan_opt="track=lCd,method=theo,thresh=level3,sup=200,gap=200:track=lDr+lDl,method=theo,thresh=level3,tele=100,sup=100,gap=100:track=ins,sup=50,cvg=5:track=del,sup=20,cvg=5";

### parsing options ###
cat("==Info: parsing options...\n")
sample_tags=strsplit(text_bam,split=':')[[1]]
n_sample=length(sample_tags); rg_files=strsplit(sample_tags,split=","); n_sp=n_sample
if(text_outprefix=='input') { 
  rg_prefix=rg_files
  sp_prefix=gsub("\\.\\w+$","",sapply(rg_files,lcPrefix),perl=T)
  out_prefix=sp_prefix[n_sample]
} else {
  out_prefix=text_outprefix
}

swan_files=NULL;swan_parfs=list();disc_files=NULL;seqcbs_files=NULL;seqcbs_parf=NULL
bigd_files=NULL;sclip_parf=NULL;sclip_file=NULL;ovrd_files=NULL;swan_par=NULL
#NOTE: ovrd_file to be compatible with safe_read function
stat_files=lapply(rg_files,function(x){ gsub(".bam$",".stat",x) }) #we need stat file for coverage, rl
if(text_stat!=""){
  stat_files=as.list(strsplit(text_stat,split=':')[[1]])
}
nstat=sapply(stat_files,length)
if(max(nstat)>1){
  stop("Error Input: must provide stat file by -t if multiple lib bams \n")
}
if(text_swan!=""){ 
  swan_files=as.list(strsplit(text_swan,split=':')[[1]])
 	if(text_override!="") ovrd_files=as.list(strsplit(text_override,split=':')[[1]])
}
if(text_bigd!=""){ 
  bigd_files=as.list(strsplit(text_bigd,split=':')[[1]])
}
if(text_disc!=""){ 
  disc_files=as.list(strsplit(text_disc,split=':')[[1]])
}
if(text_seqcbs!=""){ #only available when doing one sample
  #seqcbs_files=strsplit(text_seqcbs,split=':')[[1]]
  seqcbs_files=text_seqcbs
  seqcbs_parf=gsub(".txt",".par.txt",seqcbs_files)
}
if(text_sclip!=""){ #only available when doing one sample
  #sclip_file=strsplit(text_sclip,split=':')[[1]]
  sclip_file=text_sclip
  sclip_parf=gsub(".RData",".par.txt",sclip_file)
}
savevcf=if(text_savevcf=="") FALSE else TRUE 
sample_id=strsplit(text_sample,split=",")[[1]][1]
fineconf=opt_fineconf
confirm=strsplit(opt_confirm,split=",")[[1]]
if(!fineconf & ( "dehot" %in% confirm | "bam" %in% confirm | "dream" %in% confirm )) stop("--Error: cann't use hot and bam confirmation without .bam files\n")
ref=FaFile(ref_file)
ref_seq=scanFa(ref,param=scanFaIndex(ref)) #this returns a DNAStringSet
ref_size=length(ref_seq)
ref_name=names(ref_seq)
if(text_chrname %in% c("A","a","ALL","all","All")){
  seq_name = names(ref_seq)[seq_len(min(length(ref_seq),24))]
  cat("do all chrs:",seq_name,"for",n_sp,"samples \n")
} else {
	seq_name = strsplit(text_chrname,split=",")[[1]]
}
ref_len = sapply(seq_name, function(x) { length(ref_seq[[x]]) }); 
#get coverage
cat("==Info: getting coverage...\n")
bamfile1=rg_files[[n_sp]]; mean_cvg1=list(); rl1=rep(NA,length(bamfile1)); cvg1=rep(NA,length(bamfile1));
bamfile0=NULL; mean_cvg0=NULL; 
for(ix in seq_along(ref_seq)){
  sn=names(ref_seq)[ix]; mean_cvg1[[sn]]=0; idx=1;
  for(bam_file in bamfile1){
    cat("-Info: calculating coverage...",bam_file,",",sn,",")
    #this is assumed that genome are evenly covered
    #it is user's responsibility to provide a .stat file by running swan_stat on the bam input file
    tryCatch({
      rl1[idx]=read.table(stat_files[[n_sp]], header=T)$rl[1]
      mean_cvg_tmp=sum(read.table(stat_files[[n_sp]], header=T)$cvg)
    }, error=function(e) {
      cat("-Warn: cant load from ", stat_files[[n_sp]], " ! Check your swan_stat log\n")
    })
    
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
cat("==Info: getting options..\n")
if(text_swan_opt==""){
	if(n_sp>1){
		text_swan_opt=swan_opt_default2
	} else {
		text_swan_opt=swan_opt_default1
	}
}
if(text_bigd_opt==""){
	if(n_sp>1){
		text_bigd_opt=bigd_opt_default2
	} else {
		text_bigd_opt=bigd_opt_default1
	}
}
if(text_disc_opt==""){
	if(n_sp>1){
		text_disc_opt=disc_opt_default2
	} else {
		text_disc_opt=disc_opt_default1
	}
}
swan_opt=lapply(strsplit(text_swan_opt,":")[[1]],parse_opt2)
for(i in seq_along(swan_opt)) names(swan_opt[[i]])=sapply(lapply(swan_opt[[i]],"[[",1),"[",1) #track is the keyword
bigd_opt=lapply(strsplit(text_bigd_opt,":")[[1]],parse_opt)
disc_opt=lapply(strsplit(text_disc_opt,":")[[1]],parse_opt)
disc_opt[[n_sp]]$cvg=mean_cvg1; if(n_sp>1) disc_opt[[1]]$cvg=mean_cvg0
seqcbs_opt=parse_opt(text_seqcbs_opt)
sclip_opt=parse_opt(text_sclip_opt)

cat("==Info: parsing track files..\n")
cat("==Info: bam_files:\n"); print(rg_files)
cat("==Info: swan_files:\n"); print(swan_files)
cat("==Info: swan_parfs:\n"); print(swan_parfs)
cat("==Info: bigd_files:\n"); print(bigd_files)
cat("==Info: disc_files:\n"); print(disc_files)
cat("==Info: seqcbs_files:\n"); print(seqcbs_files)
cat("==Info: sclip_file:\n"); print(sclip_file)
cat("==Info: ovrd_files:\n"); print(ovrd_files)
cat("==Info: ref_file:\n"); print(ref_file)
swan_score = NULL
if(length(seq_name)>1 || text_chrname %in% c("A","a","ALL","all","All")) { # if >1 chr, input MUST be %pf.%chr.swan.txt.gz
  if(length(swan_files)>0){
		cat("==Info: setting swan files..\n")
  	for(i in seq_len(n_sp)){
  		swan_files[[i]]=paste(swan_files[[i]],seq_name,"swan.txt.gz",sep=".")
  		swan_parfs[[i]]=gsub("swan.txt.gz","swan.par.txt",swan_files[[i]])
      swan_score[[i]]=lapply(swan_files[[i]],function(x){read_com(x,comment_char="#",colClasses="numeric")[[2]]})
  	}
  }
  if(length(bigd_files)>0){
		cat("==Info: setting bigd files..\n")
  	for(i in seq_len(n_sp)){
  		bigd_files[[i]]=paste(bigd_files[[i]],seq_name,"bigd.txt",sep=".")
  		swan_parfs[[i]]=gsub(".bigd.txt",".swan.par.txt",bigd_files[[i]])
  	}
  }
  if(length(disc_files)>0){
		cat("==Info: setting disc files..\n")
  	for(i in seq_len(n_sp)){
  		disc_files[[i]]=paste(disc_files[[i]],seq_name,"disc.txt",sep=".")
  		swan_parfs[[i]]=gsub(".disc.txt",".swan.par.txt",disc_files[[i]])
  	}
  }
} else { # only one chr
	if(length(swan_files)>0){
		cat("==Info: setting swan file..\n")
		for(i in seq_len(n_sp)){
    	#if(!grepl(swan_files[[i]],"txt")) stop("error: must specify full file name for single chr scan\n")
 			swan_parfs[[i]]=gsub("swan.txt.gz","swan.par.txt",swan_files[[i]])
		}
	}
	if(length(bigd_files)>0){
		cat("==Info: setting bigd file..\n")
		for(i in seq_len(n_sp)){
    	#if(!grepl(swan_files[[i]],"txt")) stop("error: must specify full file name for single chr scan\n")
  		swan_parfs[[i]]=gsub(".bigd.txt",".swan.par.txt",bigd_files[[i]])
		}
	}
	if(length(disc_files)>0){
		cat("==Info: setting disc file..\n")
		for(i in seq_len(n_sp)){
    	#if(!grepl(swan_files[[i]],"txt")) stop("error: must specify full file name for single chr scan\n")
  		swan_parfs[[i]]=gsub(".disc.txt",".swan.par.txt",disc_files[[i]])
		}
	}
}
cat("==Info: parsing par file(s)..\n")
if(!length(swan_parfs)==0) {
	for(i in seq_len(n_sp)){
		swan_par[[i]]=lapply(swan_parfs[[i]],read.table,header=TRUE)
	}
}
conf_opt=list(max_del=0.9,min_dup=1.1,p_thresh=0.05,max_run=100,min_sc=10,min_mpr=5,diff_thresh=100,lrt_thresh=3,n_sp=n_sp,r_lst=c(0.2,0.5,0.95),gtk=gtk,verbose=verbose,debug=debug);
cat("==Info: showing parameters..\n")
cat("==Info: swan_par:\n"); print(swan_par)
cat("==Info: swan_opt:\n"); print(swan_opt)
cat("==Info: bigd_opt:\n"); print(bigd_opt)
cat("==Info: disc_opt:\n"); print(disc_opt)
cat("==Info: seqcbs_opt:\n"); print(seqcbs_opt)
cat("==Info: sclip_opt:\n"); print(sclip_opt)
cat("==Info: conf_opt:\n"); print(conf_opt)

### check coherence of input ###
if(length(swan_files)==0&length(seqcbs_files)==0&length(sclip_file)==0&!length(bigd_files)==0&!length(disc_files)==0) { #we need at least one scan
  stop("=Error: need at least one set of input from SWAN, BigDel, SeqCBS or Sclip scan\n!") 
} else {
  if(!length(seqcbs_files)==0&length(seqcbs_files)>1) stop("=Error: can only accept upto one seqcbs file")
  if(!length(sclip_file)==0&length(sclip_file)>1) stop("=Error: can only accept upto one sclip file")
  if(!length(swan_files)==0&length(swan_files)>2) stop("=Error: can only accept upto two swan file")
  if(!length(bigd_files)==0&length(bigd_files)>2) stop("=Error: can only accept upto two bigd file")
  if(!length(disc_files)==0&length(disc_files)>2) stop("=Error: can only accept upto two disc file")
  if(n_sample==1&length(swan_files)==0&length(bigd_files)==0&length(disc_files)==0&length(sclip_file)==0) {
    stop("=Error: can only perform one sample calling with at least one SWAN,BigDel,Discordant,SCLIP file")
  } else if (n_sample!=2) {
    if(!length(seqcbs_files)==0|length(sclip_file)>=2|length(swan_files)>=2|length(bigd_files)>=2|length(disc_files)>=2) 
    stop("=Error: can only perform one sample calling with one SWAN and/or Discordant and/or BigDel file and without SeqCBS, Sclip file, missing one samples bam file?")
  }
  if(n_sample==2){
    if(!length(swan_files)==0&length(swan_files)!=n_sample)
      stop("=Error: need exact two SWAN files for perform calling with matched samples \n!")
    if(!length(bigd_files)==0&length(bigd_files)!=n_sample)
      stop("=Error: need exact two BigDel files for perform calling with matched samples \n!")
    if(!length(disc_files)==0&length(disc_files)!=n_sample)
      stop("=Error: need exact two Discordant files for perform calling with matched samples \n!")
    if(length(seqcbs_files)>1|length(sclip_file)>1)
      stop("=Error: can only accept upto to one SeqCBS/Sclip file for calling with matched samples \n!")
  }
} #what else do we need to check?

coverage1=mean(unlist(mean_cvg1),na.rm=T)
rl=median(unlist(rl1),na.rm=T)
sclip_opt=list()
sclip_opt$rl=rl1
sclip_opt$seq_name=seq_name
scalefrac=1
if(n_sp>1){
  coverage0=mean(unlist(mean_cvg0),na.rm=T)
  scalefrac=coverage1/coverage0
}
print("sclip_opt"); print(sclip_opt)

### starting join ###
reg_seqcbs=list();  reg_ldr_ldl=list(); reg_lcd=list(); reg_sclip=list(); reg_bigd=list(); reg_disc=list(); reg_cvg=list(); reg_haf_har=list();
reg_seqcbs_good=NULL; reg_del=NULL; reg_ins=NULL;
cat("scores:",names(swan_opt[[1]]),"\n")
cat("seq_name:",seq_name,"\n")

if(length(swan_files)==0|length(swan_files)>2) { #ensure swan_files is length 1 to 2
  cat("=Warn: skip SWAN call due to no inputs or incorrect inputs number\n")
} else { #standalone case
  if(!length(ovrd_files)==0){
    if(length(swan_files)!=length(ovrd_files))
      stop("override files must be the same number as swan files")
	} else {
		print(length(swan_files))
		ovrd_files=rep(list(NULL),length(swan_files))
		cat("ovrd_files\n"); print(ovrd_files)
	}
  if(length(swan_files)==1){ #only one sample

    cat("=Info: calling small deletion ...")
		reg_del=lapply(seq_along(seq_name),function(x){call_del(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],swan_opt[[1]])})
		reg_del=unlist(reg_del,recursive=F)
		reg_del[sapply(reg_del, is.null)] <- NULL
    cat("-Info: called", length(reg_del), "in", taggie(gtk), "s\n")

    cat("=Info: calling small insertion ...")
		reg_ins=lapply(seq_along(seq_name),function(x){call_ins(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],swan_opt[[1]])})
		reg_ins=unlist(reg_ins,recursive=F)
		reg_ins[sapply(reg_ins, is.null)] <- NULL
    cat("-Info: called", length(reg_ins),"in", taggie(gtk), "s\n")

    cat("=Info: calling cvg deletion / duplication ...")
		reg_cvg=lapply(seq_along(seq_name),function(x){call_cvg(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
		reg_cvg=unlist(reg_cvg,recursive=F)
		reg_cvg[sapply(reg_cvg, is.null)] <- NULL
    #reg_cvg=call_cvg(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
    cat("-Info: called", length(reg_cvg),"in", taggie(gtk), "s\n")

    if(all(sapply(swan_par[[1]],"[[",'lCd'))){
      cat("=Info: calling lcd ...")
			if(!debug){
				reg_lcd=lapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
			} else {
				reg_lcd=mclapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
			}
			reg_lcd=unlist(reg_lcd,recursive=F)
			reg_lcd[sapply(reg_lcd, is.null)] <- NULL
      cat("-Info: called", length(reg_lcd),"in", taggie(gtk), "s\n")
    }

    if(all(sapply(swan_par[[1]],"[[",'lDl'))&all(sapply(swan_par[[1]],"[[","lDr"))){
      cat("=Info: calling ldl and ldr ...")
			reg_ldr_ldl=lapply(seq_along(seq_name),function(x){call_ldr_ldl(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
			reg_ldr_ldl=unlist(reg_ldr_ldl,recursive=F)
			reg_ldr_ldl[sapply(reg_ldr_ldl, is.null)] <- NULL
      cat("-Info: called", length(reg_ldr_ldl),"in", taggie(gtk), "s\n")
    }

    cat("=Info: calling HAF and HAR ...")
		reg_haf_har=lapply(seq_along(seq_name),function(x){call_haf_har(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
		reg_haf_har=unlist(reg_haf_har,recursive=F)
		reg_haf_har[sapply(reg_haf_har, is.null)] <- NULL
    #reg_haf_har=call_haf_har(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
    cat("-Info: called", length(reg_haf_har),"in", taggie(gtk), "s\n")

  } else { #2 swan files, matched case, find reg in spX not in spY

    reg_lcd_spX=list(); reg_lcd_spY=list(); reg_ldr_ldl_spX=list(); reg_ldr_ldl_spY=list()
    reg_cvg_spX=list(); reg_cvg_spY=list(); reg_haf_har_spX=list(); reg_haf_har_spY=list()
    if(length(swan_opt)==1) swan_opt[[2]]=swan_opt[[1]] ##TODO: smart way
		reg_cvg_spX=lapply(seq_along(seq_name),function(x){call_cvg(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],ovrd_files[[2]],swan_opt[[2]],mean_cvg1[[seq_name[x]]])})
		reg_cvg_spX=unlist(reg_cvg_spX,recursive=F)
		reg_cvg_spX[sapply(reg_cvg_spX, is.null)] <- NULL
		reg_cvg_spY=lapply(seq_along(seq_name),function(x){call_cvg(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg0[[seq_name[x]]])})
		reg_cvg_spY=unlist(reg_cvg_spY,recursive=F)
		reg_cvg_spY[sapply(reg_cvg_spY, is.null)] <- NULL
    #reg_cvg_spX=call_cvg(swan_files[2],swan_par[[2]],ovrd_files[2],swan_opt[[2]],seqname)
    #reg_cvg_spY=call_cvg(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
    reg_cvg=reg_diff(reg_cvg_spX, reg_cvg_spY)
		
    if(all(sapply(swan_par[[1]],"[[",'lCd'))&all(sapply(swan_par[[2]],"[[","lCd"))){
      cat("=Info: calling lcd ...")
			if(!debug){
				reg_lcd_spX=lapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],ovrd_files[[2]],swan_opt[[2]],mean_cvg1[[seq_name[x]]])})
				reg_lcd_spY=lapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg0[[seq_name[x]]])})
			} else {
				reg_lcd_spX=mclapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],ovrd_files[[2]],swan_opt[[2]],mean_cvg1[[seq_name[x]]])})
				reg_lcd_spY=mclapply(seq_along(seq_name),function(x){call_lcd(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg0[[seq_name[x]]])})
			}
			reg_lcd_spX=unlist(reg_lcd_spX,recursive=F)
			reg_lcd_spX[sapply(reg_lcd_spX, is.null)] <- NULL
			reg_lcd_spY=unlist(reg_lcd_spY,recursive=F)
			reg_lcd_spY[sapply(reg_lcd_spY, is.null)] <- NULL
      #reg_lcd_spX=call_lcd(swan_files[2],swan_par[[2]],ovrd_files[2],swan_opt[[2]],seqname)
      #reg_lcd_spY=call_lcd(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
      reg_lcd=reg_diff(reg_lcd_spX, reg_lcd_spY)
      cat("-Info: called", length(reg_lcd),"in", taggie(gtk), "s\n")
    }

    if(all(sapply(swan_par[[1]],"[[",'lDl'))&all(sapply(swan_par[[1]],"[[","lDr"))&all(sapply(swan_par[[2]],"[[",'lDl'))&all(sapply(swan_par[[2]],"[[","lDr"))){
      cat("=Info: calling ldl and ldr ...")
			reg_ldr_ldl_spX=lapply(seq_along(seq_name),function(x){call_ldr_ldl(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],ovrd_files[[2]],swan_opt[[2]],mean_cvg1[[seq_name[x]]])})
			reg_ldr_ldl_spX=unlist(reg_ldr_ldl_spX,recursive=F)
			reg_ldr_ldl_spX[sapply(reg_ldr_ldl_spX, is.null)] <- NULL
			reg_ldr_ldl_spY=lapply(seq_along(seq_name),function(x){call_ldr_ldl(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg0[[seq_name[x]]])})
			reg_ldr_ldl_spY=unlist(reg_ldr_ldl_spY,recursive=F)
			reg_ldr_ldl_spY[sapply(reg_ldr_ldl_spY, is.null)] <- NULL
      #reg_ldr_ldl_spX=call_ldr_ldl(swan_files[2],swan_par[[2]],ovrd_files[2],swan_opt[[2]],seqname)
      #reg_ldr_ldl_spY=call_ldr_ldl(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
      reg_ldr_ldl=reg_diff(reg_ldr_ldl_spX, reg_ldr_ldl_spY)
      cat("-Info: called", length(reg_ldr_ldl),"in", taggie(gtk),"s\n")
    }

    cat("=Info: calling HAF and HAR ...")
		reg_haf_har_spX=lapply(seq_along(seq_name),function(x){call_haf_har(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],ovrd_files[[2]],swan_opt[[2]],mean_cvg1[[seq_name[x]]])})
		reg_haf_har_spX=unlist(reg_haf_har_spX,recursive=F)
		reg_haf_har_spX[sapply(reg_haf_har_spX, is.null)] <- NULL
		reg_haf_har_spY=lapply(seq_along(seq_name),function(x){call_haf_har(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],ovrd_files[[1]],swan_opt[[1]],mean_cvg1[[seq_name[x]]])})
		reg_haf_har_spY=unlist(reg_haf_har_spY,recursive=F)
		reg_haf_har_spY[sapply(reg_haf_har_spY, is.null)] <- NULL
    #reg_haf_har_spX=call_haf_har(swan_files[2],swan_par[[2]],ovrd_files[2],swan_opt[[2]],seqname)
    #reg_haf_har_spY=call_haf_har(swan_files[1],swan_par[[1]],ovrd_files[1],swan_opt[[1]],seqname)
    reg_haf_har=reg_diff(reg_haf_har_spX, reg_haf_har_spY)
    cat("-Info: called", length(reg_haf_har),"in", taggie(gtk),"s\n")

    cat("=Info: calling small deletion ...")
		reg_del_spX=lapply(seq_along(seq_name),function(x){if(debug) { cat("spX:",seq_name[x],"\n") }; call_del(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],swan_opt[[2]])})
		reg_del_spX=unlist(reg_del_spX,recursive=F)
		reg_del_spX[sapply(reg_del_spX, is.null)] <- NULL
		reg_del_spY=lapply(seq_along(seq_name),function(x){if(debug) cat("spY:",seq_name[x],"\n"); call_del(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],swan_opt[[1]])})
		reg_del_spY=unlist(reg_del_spY,recursive=F)
		reg_del_spY[sapply(reg_del_spY, is.null)] <- NULL
    reg_del=reg_diff(reg_del_spX,reg_del_spY)
    cat("-Info: called", length(reg_del),"in", taggie(gtk),"s\n")

    cat("=Info: calling small insertion ...")
		reg_ins_spX=lapply(seq_along(seq_name),function(x){if(debug) cat("spX:",seq_name[x],"\n"); call_ins(seq_name[x],swan_files[[2]][x],swan_par[[2]][[x]],swan_opt[[2]])})
		reg_ins_spX=unlist(reg_ins_spX,recursive=F)
		reg_ins_spX[sapply(reg_ins_spX, is.null)] <- NULL
		reg_ins_spY=lapply(seq_along(seq_name),function(x){if(debug) cat("spY:",seq_name[x],"\n"); call_ins(seq_name[x],swan_files[[1]][x],swan_par[[1]][[x]],swan_opt[[1]])})
		reg_ins_spY=unlist(reg_ins_spY,recursive=F)
		reg_ins_spY[sapply(reg_ins_spY, is.null)] <- NULL
    reg_ins=reg_diff(reg_ins_spX,reg_ins_spY)
    cat("-Info: called", length(reg_ins),"in", taggie(gtk),"s\n")

  }
} #end swan call

if(length(bigd_files)==0|length(bigd_files)>2) { #ensure bigd_files is length 1 to 2
  cat("=Warn: skip bigd call due to no inputs or incorrect inputs number\n")
} else { #standalone case
  cat("=Info: calling big deletion ...")
  if(length(bigd_files)==1){
		reg_bigd=lapply(seq_along(seq_name),function(x){call_bigd(seq_name[x],bigd_files[[1]][x],swan_par[[1]][[x]],bigd_opt[[1]],swan_files[[1]][x],mean_cvg1[[seq_name[x]]])})
		reg_bigd=unlist(reg_bigd,recursive=F)
		reg_bigd[sapply(reg_bigd, is.null)] <- NULL
  } else { #2 bigd files, matched case, find reg in spX not in spY
    if(length(bigd_opt)==1) bigd_opt[[2]]=bigd_opt[[1]] ##TODO: smart way
		reg_bigd_spX=lapply(seq_along(seq_name),function(x){call_bigd(seq_name[x],bigd_files[[2]][x],swan_par[[2]][[x]],bigd_opt[[2]],swan_files[[2]][x],mean_cvg1[[seq_name[x]]])})
		reg_bigd_spX=unlist(reg_bigd_spX,recursive=F)
		reg_bigd_spX[sapply(reg_bigd_spX, is.null)] <- NULL
		reg_bigd_spY=lapply(seq_along(seq_name),function(x){call_bigd(seq_name[x],bigd_files[[1]][x],swan_par[[1]][[x]],bigd_opt[[1]],swan_files[[1]][x],mean_cvg0[[seq_name[x]]])})
		reg_bigd_spY=unlist(reg_bigd_spY,recursive=F)
		reg_bigd_spY[sapply(reg_bigd_spY, is.null)] <- NULL
    cat("-Info: bigd clusters spX:",length(reg_bigd_spX),"spY:",length(reg_bigd_spY),"\n")
    reg_bigd=reg_diff(reg_bigd_spX, reg_bigd_spY)
  }
  cat("-Info: called", length(reg_bigd),"\n")
} #end bigd call

if(length(disc_files)==0|length(disc_files)>2) { #ensure disc_files is length 1 to 2
  cat("=Warn: skip disc call due to no inputs or incorrect inputs number\n")
} else { #standalone case
  cat("=Info: calling discordant...")
  if(length(disc_files)==1){
		reg_disc=lapply(seq_along(seq_name),function(x){call_disc(seq_name[x],disc_files[[1]][x],swan_par[[1]][[x]],disc_opt[[1]],swan_files[[1]][x],mean_cvg1[[seq_name[x]]])})
		reg_disc=unlist(reg_disc,recursive=F)
		reg_disc[sapply(reg_disc, is.null)] <- NULL
  } else { #2 disc files, matched case, find reg in spX not in spY
    if(length(disc_opt)==1) disc_opt[[2]]=disc_opt[[1]] ##NOTE: now have two defaults based on input number
		reg_disc_spX=lapply(seq_along(seq_name),function(x){call_disc(seq_name[x],disc_files[[2]][x],swan_par[[2]][[x]],disc_opt[[2]],swan_files[[2]][x],mean_cvg1[[seq_name[x]]])})
		reg_disc_spX=unlist(reg_disc_spX,recursive=F)
		reg_disc_spX[sapply(reg_disc_spX, is.null)] <- NULL
		reg_disc_spY=lapply(seq_along(seq_name),function(x){call_disc(seq_name[x],disc_files[[1]][x],swan_par[[1]][[x]],disc_opt[[1]],swan_files[[1]][x],mean_cvg0[[seq_name[x]]])})
		reg_disc_spY=unlist(reg_disc_spY,recursive=F)
		reg_disc_spY[sapply(reg_disc_spY, is.null)] <- NULL
    cat("-Info: disc clusters spX:",length(reg_disc_spX),"spY:",length(reg_disc_spY),"\n")
    reg_disc=reg_diff(reg_disc_spX,reg_disc_spY)
  }
  cat("-Info: called", length(reg_disc),"\n")
} #end disc call

if(length(seqcbs_files)==0|length(seqcbs_files)!=1) { #ensure seqcbs_files is length 1
  cat("=Warn: skip seqcbs call due to no inputs or incorrect inputs number\n")
} else {
	tmp_seqcbs=lapply(seq_along(seq_name),function(x){call_seqcbs(seq_name[x],seqcbs_files,seqcbs_parf,seqcbs_opt,rg_files)})
	list[reg_seqcbs,reg_seqcbs_good]=do.call(Map, c(c, tmp_seqcbs))
	reg_seqcbs[sapply(reg_seqcbs, is.null)] <- NULL
	reg_seqcbs_good[sapply(reg_seqcbs_good, is.null)] <- NULL
  #list[reg_seqcbs,reg_seqcbs_good]=call_seqcbs(seqcbs_files,seqcbs_parf,seqcbs_opt,seqname,rg_files)
} #end seqcbs call

if(length(sclip_file)==0|length(sclip_file)!=1) { #ensure sclip_file is length 1
  cat("=Warn: skip sclip call due to no inputs or incorrect inputs number\n")
} else {
  cat("=Info: calling sclip...")
  reg_sclip=call_sclip(sclip_file,sclip_opt,bamfile1,ref_seq) #sclip_file is always sclip_file
  cat("-Info: called", length(reg_sclip),"\n")
} #end sclip call
 
reg_raw=c(reg_seqcbs,reg_lcd,reg_ldr_ldl,reg_haf_har,reg_bigd,reg_disc,reg_cvg,reg_sclip,reg_seqcbs_good,reg_del,reg_ins)
reg_raw[sapply(reg_raw, is.null)] <- NULL
reg_raw=lapply(seq_len(length(reg_raw)),prep_reg,reg_raw)
reg_raw[sapply(reg_raw, is.null)] <- NULL
if(savevcf) {
	vcf_par=parse_opt(text_savevcf,fs=c(":","%"))
	species=vcf_par[["species"]]
	sample_info=matrix(strsplit(text_sample,split=",")[[1]], ncol=4)
	contig_info = cbind(seq_name,ref_len,sapply(seq_along(seq_name), function(x){ digest(toString(ref_seq[[x]]),algo="md5") }),species)
} else {
	sample_info = matrix(rep(NA,4), ncol=4)
	contig_info = cbind(seq_name,ref_len,sapply(seq_along(seq_name), function(x){ NA }),NA)
}
list[vcf_raw,vcf_meta,bed_raw]=vcf_merge(reg_raw,sample_id,sample_info,contig_info,ref_seq,ref_file)
if(nrow(bed_raw)!=0) { 
	bed_raw=bed_raw[with(bed_raw, order(CHROM, POS)), ]
	vcf_raw=vcf_raw[with(vcf_raw, order(CHROM, POS)), ] 
} else { cat("-swan_join safely done!", taggie(gtk),"\n"); 
  if(debug) { cat("=Info: warnings if any\n"); warnings() }
  cat("JOINDONE\n"); quit() }
for(vi in seq_len(nrow(vcf_raw))) vcf_raw$ID[vi]=paste(sample_id,vcf_raw$CHROM[vi],vi,vcf_raw$ALT[vi],vcf_raw$ID[vi],sep=".")
for(vi in seq_len(nrow(bed_raw))) bed_raw$ID[vi]=paste(sample_id,bed_raw$CHROM[vi],vi,bed_raw$ALT[vi],bed_raw$ID[vi],sep=".")
if(verbose) cat("==Info: raw sorted!\n")
write_com(bed_raw,comment=NULL,filename=paste(out_prefix,".raw.bed",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
if(savevcf) write_com(vcf_raw,comment=vcf_meta,filename=paste(out_prefix,".raw.vcf",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
cat("==Info: total ",nrow(bed_raw),"raw calls, ",paste(out_prefix,".raw.bed",sep=""),"written!\n")
load_opt=list(seqname=seq_name, debug=debug, verbose=verbose) #why this need seqname

cat("-Info: prepared raw calls for confirmation and genotyping\n")
reg_good=c(reg_sclip,reg_seqcbs_good,reg_del,reg_ins); cat("-Info: directly confirmed REGIONS, total:",length(reg_good),"\n")
reg_lst=c(reg_seqcbs,reg_lcd,reg_ldr_ldl,reg_haf_har,reg_bigd,reg_disc,reg_cvg)
cat("-Info: by SCLIP:",length(reg_sclip),"of",length(reg_sclip),"confirmed\n")
cat("-Info: by CIGAR:",length(reg_del),"of",length(reg_del),"deletion confirmed\n")
cat("-Info: by CIGAR:",length(reg_ins),"of",length(reg_ins),"insertion confirmed\n")
cat("-Info: by SEQCBS:",length(reg_seqcbs_good),"of",length(reg_seqcbs)+length(reg_seqcbs_good),"duplication confirmed\n")
#reg_seqcbs_bam=reg_diff(reg_seqcbs,reg_good);cat("-Info: SEQCBS:",length(reg_seqcbs_bam),"of",length(reg_seqcbs),"to be confirmed\n")
#reg_lcd_bam=reg_diff(reg_lcd,reg_good);cat("-Info: LCD:",length(reg_lcd_bam),"of",length(reg_lcd),"to be confirmed\n")
#reg_ldr_ldl_bam=reg_diff(reg_ldr_ldl,reg_good);cat("-Info: LDL+LDR",length(reg_ldr_ldl_bam),"of",length(reg_ldr_ldl),"to be confirmed\n")
#reg_haf_har_bam=reg_diff(reg_haf_har,reg_good);cat("-Info: HAF+HAR",length(reg_haf_har_bam),"of",length(reg_haf_har),"to be confirmed\n")
#reg_bigd_bam=reg_diff(reg_bigd,reg_sclip);cat("-Info: BIGD:",length(reg_bigd_bam),"of",length(reg_bigd),"to be confirmed\n")
#reg_disc_bam=reg_diff(reg_disc,reg_sclip);cat("-Info: DISC:",length(reg_disc_bam),"of",length(reg_disc),"to be confirmed\n")
#reg_cvg_bam=reg_diff(reg_cvg,reg_sclip);cat("-Info: COVERAGE:",length(reg_cvg_bam),"of",length(reg_cvg),"to be confirmed\n")
#reg_seqcbs_bam=reg_seqcbs
#reg_lcd_bam=reg_lcd
#reg_ldr_ldl_bam=reg_ldr_ldl
#reg_haf_har_bam=reg_haf_har
#reg_bigd_bam=reg_bigd
#reg_disc_bam=reg_disc
#reg_cvg_bam=reg_cvg
#0. option to check reg_seqcbs, optimized confirmation strategy for dream data
if("dream" %in% confirm){ #use seqcbs_region as dominant calling
  cat("-Info: ", length(reg_seqcbs), "unconfirmed regions before filtering dream\n")
  seqcbs_sclip_ovlap=reg_ovlap(reg_seqcbs,reg_sclip)
  seqcbs_bigd_ovlap=reg_ovlap(reg_seqcbs,reg_bigd)
  seqcbs_disc_ovlap=reg_ovlap(reg_seqcbs,reg_disc)
  exclude_idx=c(); bigd_include_idx=c(); disc_include_idx=c()
  for(i in seq_along(seqcbs_sclip_ovlap)) {  #only look at first overlap
    if(!length(seqcbs_sclip_ovlap[[i]])==0 & !(i %in% exclude_idx)) {
      exclude_idx=c(exclude_idx,i)
    }
  } 
  #cat(exclude_idx,"\n")
  for(i in seq_along(seqcbs_bigd_ovlap)) {  #only look at first overlap
    if(!length(seqcbs_bigd_ovlap[[i]])==0 & !(i %in% exclude_idx)) {
      exclude_idx=c(exclude_idx,i); bigd_include_idx=c(bigd_include_idx,seqcbs_bigd_ovlap[[i]][1])
    }
  } 
  #cat(exclude_idx,"\n")
  for(i in seq_along(seqcbs_disc_ovlap)) {  #only look at first overlap
    if(!length(seqcbs_disc_ovlap[[i]])==0 & !(i %in% exclude_idx)) {
      exclude_idx=c(exclude_idx,i); disc_include_idx=c(disc_include_idx,seqcbs_disc_ovlap[[i]][1])
    }
  } 
  #cat(exclude_idx,"\n"); 
  seqcbs_sel=!(seq_along(reg_seqcbs) %in% exclude_idx)
  #cat(bigd_include_idx,"\n"); 
  bigd_sel=(seq_along(reg_bigd) %in% bigd_include_idx)
  #cat(disc_include_idx,"\n"); 
  disc_sel=(seq_along(reg_disc) %in% disc_include_idx)
  
  reg_seqcbs_remain=reg_seqcbs[seqcbs_sel]
  alert_idx=reg_isize_alert(reg_seqcbs_remain,bamfile1,isize_min=1000,isize_max=10000,seqname,
                             c(conc_flags,impp_flags,disc_flags))
  #cat(alert_idx,"\n")
  if(!length(alert_idx)==0) reg_seqcbs_remain=reg_seqcbs_remain[!alert_idx]
  reg_lst=c(reg_seqcbs_remain,reg_bigd[bigd_sel],reg_disc[disc_sel]) 
  cat("-Info: ", length(reg_lst), "unconfirmed regions before filtering by dream mode\n")
}

#2. option to chekc hot regions (now done by providing the gapfile of hc region ahead)
if("hot" %in% confirm){
	rle_cvg1=list();
	for(ix in seq_along(ref_seq)){
	  sn=names(ref_seq)[ix]; idx=1; rle_cvg1[[sn]]=list();
	  for(bam_file in bamfile1){
	    cat("-Info: calculating coverage...",bam_file,",",sn,"\n")
	    stats=allFunction(seqinfo=GRanges(sn,IRanges(start=1,end=max_chr_len)),bamFile=bam_file,what=c("pos","qwidth"))
	    cvg_tmp=coverage(IRanges(start=stats[["pos"]][!is.na(stats[["pos"]])],width=stats[["qwidth"]][!is.na(stats[["pos"]])]))
	    cvg_tmp=IRanges::Rle(c(base::as.vector(cvg_tmp),rep(0,max(0,max_chr_len-length(base::as.vector(cvg_tmp))))))
	    mean_cvg_tmp=sum(cvg_tmp/sum(cvg_tmp>0,na.rm=T),na.rm=T)
	    if(length(rle_cvg1[[sn]])==0) rle_cvg1[[sn]]=cvg_tmp else rle_cvg1[[sn]]=rle_cvg1[[sn]]+cvg_tmp
			idx=idx+1
	    mean_cvg1[[sn]]=mean_cvg1[[sn]]+mean_cvg_tmp
	  }
	}
	if(n_sp>1){#
	  rle_cvg0=list();
		for(ix in seq_along(ref_seq)){
	    sn=names(ref_seq)[ix]; idx=1; rle_cvg0[[sn]]=list();
	    for(bam_file in bamfile0){
				cat("-Info: calculating coverage...",bam_file,",",sn,"\n")
	      stats=allFunction(seqinfo=GRanges(sn,IRanges(start=1,end=max_chr_len)),bamFile=bam_file,what=c("pos","qwidth"))
	      cvg_tmp=coverage(IRanges(start=stats[["pos"]][!is.na(stats[["pos"]])],width=stats[["qwidth"]][!is.na(stats[["pos"]])]))
	      cvg_tmp=IRanges::Rle(c(base::as.vector(cvg_tmp),rep(0,max(0,max_chr_len-length(base::as.vector(cvg_tmp))))))
	      mean_cvg_tmp=sum(cvg_tmp/sum(cvg_tmp>0,na.rm=T),na.rm=T)
	      if(length(rle_cvg0[[sn]])==0) rle_cvg0[[sn]]=cvg_tmp else rle_cvg0[[sn]]=rle_cvg0[[sn]]+cvg_tmp
				idx=idx+1
	      mean_cvg0[[sn]]=mean_cvg0[[sn]]+mean_cvg_tmp
	    }
	  }
	}
  cat("-Info: ", length(reg_lst), "regions before filtering hotspots\n")
  conf_opt=c(conf_opt, list(mean_cvg=mean_cvg1,rle_cvg=rle_cvg1))
  if(!debug){
    reg_lst=lapply(seq_len(length(reg_lst)),conf_hot,reg_lst,conf_opt) #
  } else {
    reg_lst=mclapply(seq_len(length(reg_lst)),conf_hot,reg_lst,conf_opt)
  }
  cat("-Info: ", length(reg_lst), "regions after filtering hotspots\n")
}
#3. option to chekc hot and bam regions (this is mostly outated too)
if("bam" %in% confirm | "hot" %in% confirm){
  cat("-Info: previously:",length(reg_lst),"regions\n")
  if("bam" %in% confirm){ #bam | bam,hot, exclude hot anyway
    bam_lst=lapply(seq_len(length(reg_lst)),load_bam,reg_lst,rg_files,swan_par,load_opt);cat("-Info: loaded bam_lst for conf and geno\n")
    hot_only=FALSE
  } else { #hot
    bam_lst=list()
    hot_only=TRUE
  }
  #reg_lst=lapply(seq_len(length(reg_lst)),conf_soft,reg_lst,bam_lst,conf_opt);cat("-Info: done conf_soft\n")
  #reg_lst=lapply(seq_len(length(reg_lst)),conf_hang,reg_lst,bam_lst,swan_par,conf_opt);cat("-Info: done conf_hang\n")
  #reg_lst=lapply(seq_len(length(reg_lst)),conf_cvg,reg_lst,bam_lst,conf_opt);cat("-Info: done conf_cvg\n")
  #reg_lst=lapply(seq_len(length(reg_lst)),conf_strad,reg_lst,bam_lst,swan_par,conf_opt);cat("-Info: done conf_strad\n")
  if(!debug){
    reg_lst=lapply(seq_len(length(reg_lst)),conf_join,reg_lst,bam_lst,swan_par,conf_opt,hot_only=hot_only) #
  } else {
    reg_lst=mclapply(seq_len(length(reg_lst)),conf_join,reg_lst,bam_lst,swan_par,conf_opt,hot_only=hot_only)
  }
  reg_lst[sapply(reg_lst,length)] <- NULL
  cat("-Info: confirmed:",length(reg_lst),"regions\n")
}

#output confirmed and auto-confirmed events
reg_conf=reg_lst
reg_conf=lapply(seq_len(length(reg_conf)),prep_reg,reg_conf)
reg_conf[sapply(reg_conf, is.null)] <- NULL
reg_good=lapply(seq_len(length(reg_good)),prep_reg,reg_good)
reg_good[sapply(reg_good, is.null)] <- NULL
reg_all=c(reg_conf,reg_good)
#1. option to check overlaps 
if("dedup" %in% confirm){
	strict=ifelse("nostrict" %in% confirm,FALSE,TRUE)
	cat("-Info: strict deduplication set to be", strict ,"\n")
	reg_all=dedup_reg(reg_all,strict)
	reg_all[sapply(reg_all, is.null)] <- NULL
	cat("-Info: ", length(reg_all),  "remains after deduplication\n")
}
cat("-Info: reg_all: ", length(reg_all), "\n")
list[vcf_set,vcf_meta,bed_set]=vcf_merge(reg_all,sample_id,sample_info,contig_info,ref_seq,ref_file)
if(nrow(bed_set)!=0) { 
	vcf_set=vcf_set[with(vcf_set, order(CHROM, POS)), ]
	bed_set=bed_set[with(bed_set, order(CHROM, POS)), ]
} else { cat("-swan_join safely done!", taggie(gtk),"\n"); if(debug) { cat("=Info: warnings if any\n"); warnings() }; cat("JOINDONE\n"); quit() }
#2. option to check hotspots
if("dehotspot" %in% confirm){
	dehotspot_idx=lapply(seq_name,vcf_set,dehotstpot_reg)
	dehotspot_idx=unlist(dehotspot_idx)
	if(length(dehotspot_idx)>0){
		vcf_set=vcf_set[dehotspot,]
		bed_set=bed_set[dehotspot,]
	} else { cat("-swan_join safely done!", taggie(gtk),"\n"); if(debug) { cat("=Info: warnings if any\n"); warnings() }; cat("JOINDONE\n"); quit() }
}
#we might want to replace vcf_set$ID by sample_id.seq_name.sv_idx.sv_type
for(vi in seq_len(nrow(vcf_set))) vcf_set$ID[vi]=paste(sample_id,vcf_set$CHROM[vi],vi,vcf_set$ALT[vi],vcf_set$ID[vi],sep=".")
for(vi in seq_len(nrow(bed_set))) bed_set$ID[vi]=paste(sample_id,bed_set$CHROM[vi],vi,bed_set$ALT[vi],bed_set$ID[vi],sep=".")
if(verbose) cat("==Info: vcf sorted!\n")
if(savevcf) write_com(vcf_set,comment=vcf_meta,filename=paste(out_prefix,".conf.vcf",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
write_com(bed_set,comment=NULL,filename=paste(out_prefix,".conf.bed",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
cat("==Info: total ",nrow(bed_set),"confirmed calls,",paste(out_prefix,".conf.bed",sep=""),"written!\n")
cat("==Info: swan_join safely done!", taggie(gtk),"\n")
if(debug) { cat("=Info: warnings if any\n"); warnings() }
cat("please find outputs in:",paste(out_prefix,".conf.vcf",sep=""),"and/or",paste(out_prefix,".conf.bed",sep=""),"\n",sep=" ")
cat("JOINDONE\n")
