#!/usr/bin/env Rscript
#this script currently only does germline genotyping for single sample analysis
#with homo 1/1 and hete 0/1 cases 
#all reads have to be from single library i.e. same isize
#TODO:
#1. hete 1/2 in single sample are more involved scenario which has to be worked out later
#2. somatic genotyping based on pair sample may be done indirectly in following way:
#  first, assemble normal sample of the proposed region
#  second, using the assembled normal sequence(s) as reference
#  third, do geno_reg based on newly assembled references
#  however, this is not implemented yet
version="REPLACE_WITH_COMMIT_OR_VERSION"
gtk=proc.time()[3]
suppressMessages(library(optparse))
suppressMessages(library(Biostrings))
suppressMessages(library(digest))
suppressMessages(library(parallel))
suppressMessages(library(swan))
gmk=get_gmk(Sys.getpid())
options(warn=2)

geno_default=list(sam_par="", #reserved
									fas_par="", #reserved
									fna_par="", #reserved
									bam_rl=100, #read length
									bam_is=300, #insert size mean 
									bam_issd=30, #insert size sd
									bam_cvg=50, #coverage mean
								  bam_buf=400, #buf in loaded bam 
									fna_ref="~/hg/hg19/human_g1k_v37.fasta.2bit", #2bit of reffile
									ex_ovd=50, #min end overlapping for exonerate deletion align 
									ex_ovi=100, #min end overlapping for exonerate insertion align  
									ex_clo=0.75, #min coverage reduction for exonerate heterozygous align  
									ex_chi=0.75, #min coverage increase for exonerate duplication align  
									vv_kmer=31, #kmer used in velvet
									vv_isize=300, #isize used in velvet
									vv_minctg=150, #minimum contig size used in velvet
									vv_mincvg=4) #minimum contig coverage used in velvet

option_list <- list(
  make_option(c("-o", "--outprefix"), default="input",
    help="prefix for output file [default %default]"),
	make_option(c("-s", "--sample"), default="spX,INFO,MIX,DESCRIPTION",
    help="mannual override of spX information [default %default]"),
  make_option(c("-p", "--pars"), default="",
		help="input parameters for the pipeline X=string,Y_number,... [default learned]"),
	make_option(c("-v", "--savevcf"), default="",
    help="whether to savevcf file (slower) and parameters, e.g.
          species=human_sapien:other_opt=other_value
          [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
		help="verbose mode and additional information outputs [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
    help="debug mode and additional .RData is assumed for all inputs, see manual [default %default]")
)

# move to R/libswan.R when mature
#FUNC: parse option3
parse_opt3 = function(arg_str, sep1=",", sep2=c("=","_")){
	xkeyval = function(x){ 
			key_val=list()
			if(grepl(sep2[1],x)){
				tuple=strsplit(x,sep2[1])[[1]] 
				key_val[tuple[1]]=tuple[2]
			} #return named string par
			if(grepl(sep2[2],x)){
				tuple=strsplit(x,sep2[2])[[1]] 
				key_val[tuple[1]]=as.numeric(tuple[2])
			} #return named numeric par
			return(key_val)
		}
	args=lapply(strsplit(arg_str,sep1)[[1]],xkeyval)
	args=sapply(args,"[")
	return(args)
}

parser <- OptionParser(
  usage="%prog [options] bedFile refFile bamFile",
  option_list=option_list)

cat("-Info: swan_geno.R vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE,
                 positional_arguments = TRUE)
if(length(cmd$args)!=3){ print_help(parser); quit(); }

bed_file=cmd$args[1]
ref_file=cmd$args[2]
bam_file=cmd$args[3]
print(cmd$options$pars)
debug=cmd$options$debug
par_text=cmd$options$pars
#DEBUG: bed_file="test.homo.bed"; ref_file="~/hg/hg19/human_g1k_v37.fasta"; bam_file="homo1.c50.f100.bam"
#DEBUG: bed_file="test.hete.bed"; ref_file="~/hg/hg19/human_g1k_v37.fasta"; bam_file="homo1.c50.f50.bam"
stat_file=gsub(".bam$",".stat",bam_file)
bam_stat=read.table(stat_file, header=TRUE)
learn_pars = list()
learn_pars$bam_rl = bam_stat$rl
learn_pars$bam_is = bam_stat$is
learn_pars$bam_issd = bam_stat$sdR
learn_pars$bam_cvg = bam_stat$cvg
learn_pars$bam_buf = round(ceiling((bam_stat$is + 3*bam_stat$sdR)/100.0)*100)
learn_pars$fna_ref = paste(ref_file,"2bit",sep=".")
learn_pars$ex_ovd = bam_stat$rl/2
learn_pars$ex_ovi = bam_stat$rl
learn_pars$vv_isize = bam_stat$is
learn_pars$vv_minctg = bam_stat$rl*1.5
learn_pars$vv_mincvg = max(round(bam_stat$cvg/5), 5)

#DEBUG: par_text=""
input_pars = parse_opt3(par_text)
for(nm in names(learn_pars))  if(!nm %in% names(input_pars)) input_pars[nm]=learn_pars[nm] #assign default to pars

text_sample=cmd$options$sample
text_savevcf = cmd$options$savevcf
#DEBUG: text_sample="spX,INFO,MIX,DESCRIPTION"
#DEBUG: text_savevcf=""
sample_id=strsplit(text_sample,split=",")[[1]][1]
savevcf = if(text_savevcf=="") FALSE else TRUE
ref = FaFile(ref_file)
cat("scanning reference sequences...\n")
ref_seq = scanFa(ref,param=scanFaIndex(ref)) #this returns a DNAStringSet
ref_size = length(ref_seq)
ref_name = names(ref_seq)
seq_name = names(ref_seq)[seq_len(min(length(ref_seq),24))]
ref_len = sapply(seq_name, function(x) { length(ref_seq[[x]]) })

#FUNC: use parallel sapply if not in debug mode
debug_sapply = function(debug,pars){
	nc=length(pars$X)
	if(!debug) { res=do.call("sapply",pars) } else {
		res=unlist(do.call("mclapply",pars))
		nr=length(res)/nc #before transpose, nrow is features, ncol is samples
		if(length(res)!=nc) {
			fns=names(res)[seq_len(nr)]
			print(fns)
			res=matrix(res,nrow=nr,byrow=FALSE) #filling by column
			rownames(res)=fns
		}
	}
	return(res)
}

#FUNC: use parallel lapply if not in debug mode
debug_lapply = function(debug,pars){
	if(!debug) res=do.call("lapply",pars) else res=do.call("mclapply",pars)
	return(res)
}

if(savevcf) {
  vcf_par = parse_opt(text_savevcf,fs=c(":","%"))
  species = vcf_par[["species"]]
  sample_info = matrix(strsplit(text_sample,split=",")[[1]], ncol=4)
	cat("digesting reference sequences...\n")
  contig_info = cbind(seq_name,ref_len,debug_sapply(debug,list(X=seq_along(seq_name), FUN=function(x){ digest(toString(ref_seq[[x]]),algo="md5") })), species)
} else {
  sample_info = matrix(rep(NA,4), ncol=4)
  contig_info = cbind(seq_name,ref_len,debug_sapply(debug,list(X=seq_along(seq_name), FUN=function(x){ NA })),NA)
}

#' @title geno_reg
#'
#' @description \code{geno_reg} \cr
#' If called, \code{geno_reg} will return two list of sv_event structure
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
#' geno_reg(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
#' @seealso \code{\link{call_lcd}} and \code{\link{call_ldr_ldl}}
geno_reg = function(ix, data, bams, pars){
	confident=TRUE
	chr = data$CHROM[ix]
	st = round(floor((data$POS[ix] - pars$bam_buf)/100.0)*100)
	ed = round(ceiling((data$END[ix] + pars$bam_buf)/100.0)*100)
	type = data$ALT[ix]
	cat("geno_reg", ix, ":", "[", type, "]", chr, ":", st,"-", ed, "\n")
	#DEBUG: pars=list(); bams=bam_file; chr=10; st=3999400; ed=4001200; type="DEL"; confident=T; 
	bamfn=bams[1] #only allow one bam  
	for(nm in names(geno_default))  if(!nm %in% names(pars)) pars[nm] = geno_default[nm] #assign default to pars 
	
	### TODO: sam2aln()
	#assume all bams are used to generate assembly

	#samtools view -h homo1.c50.f100.bam 1:1099380-1101300 >asm1.sam
	#bam2fastx -A -a -o asm1.fa -s asm1.sam -N
	#twoBitToFa ~/hg/hg19/human_g1k_v37.ucsc.fasta.2bit -seq=chr1 -start=1099380 -end=1101300 asm1.fa
	#velveth asm1 31 -sam -shortPaired asm1.sam
	#velvetg asm1/ -cov_cutoff 4 -min_contig_lgth 150 -ins_length 300
	#exonerate asm1/contigs.fa asm1.fa >asm1.exon
	dir.create(t<-paste(tempdir(),paste("sv",ix,sep=""),sep="/"), FALSE, TRUE, "0700")
  vvd_tmp = t
	#vvd_tmp = tempdir(tpaste(tempdir(),ix,sep="-")
	sam_tmp = paste(vvd_tmp,"swan_sam_tmp.sam",sep="/")
	fas_tmp = paste(vvd_tmp,"swan_fas_tmp.fas",sep="/")
	fna_tmp = paste(vvd_tmp,"swan_fna_tmp.fna",sep="/")
	aln_tmp = paste(vvd_tmp,"swan_aln_tmp.aln",sep="/")
	bla_tmp = paste(vvd_tmp,"swan_bla_tmp.bla",sep="/")
	vug_tmp = paste(vvd_tmp,"swan_vug_tmp.txt",sep="/")
	#sam_tmp = tempfile(pattern = "swan_sam_tmp", tmpdir = vvd_tmp, fileext = ".sam")
	#fas_tmp = tempfile(pattern = "swan_fas_tmp", tmpdir = vvd_tmp, fileext = ".fas")
	#fna_tmp = tempfile(pattern = "swan_fna_tmp", tmpdir = vvd_tmp, fileext = ".fna")
	#aln_tmp = tempfile(pattern = "swan_aln_tmp", tmpdir = vvd_tmp, fileext = ".aln")
	#bla_tmp = tempfile(pattern = "swan_bla_tmp", tmpdir = vvd_tmp, fileext = ".bla")
	#vug_tmp = tempfile(pattern = "swan_vug_tmp", tmpdir = vvd_tmp, fileext = ".txt")
	ctg_tmp = paste(vvd_tmp,"contigs.fa",sep="/")

	tryCatch({
	sam_tag = paste(pars$sam_par,bamfn,paste(chr,":",nsf(st),"-",nsf(ed),sep=""),">",sam_tmp)
	sam_cmd = sprintf("samtools view -h %s",sam_tag)
	print(sam_cmd)
	sam_run = system(sam_cmd,intern=TRUE) #check by samtools view -S $sam_tmp | wc -l
	fas_tag = paste(pars$fas_par,"-o",fas_tmp,"-s",sam_tmp)
	fas_cmd = sprintf("bam2fastx -A -a -N %s", fas_tag)
	print(fas_cmd)
	fas_run = system(fas_cmd,intern=TRUE) #check by cat $fas_tmp | wc -l
	fna_tag = paste(pars$fna_par, paste(c("-seq","-start","-end"),c(chr,nsf(st),nsf(ed)),sep="=", collapse=" "), pars$fna_ref, fna_tmp) 
	fna_cmd = sprintf("twoBitToFa %s", fna_tag)
	print(fna_cmd)
	fna_run = system(fna_cmd,intern=TRUE) #check by cat $fna_tmp
	vvh_tag = paste(vvd_tmp, pars$vv_kmer,"-reference", fna_tmp, "-shortPaired", "-sam", sam_tmp)
	vvh_cmd = sprintf("velveth %s", vvh_tag)
	print(vvh_cmd)
	vvh_run = system(vvh_cmd,intern=TRUE) #check by ls $vvd_tmp
	vvg_tag = paste(paste(vvd_tmp,"/",sep=""), "-cov_cutoff", pars$vv_mincvg, "-min_contig_lgth", pars$vv_minctg, "-ins_length", pars$vv_isize)
	vvg_cmd = sprintf("velvetg %s", vvg_tag)
	print(vvg_cmd)
	vvg_run = system(vvg_cmd,intern=TRUE) #check by cat $vvd_tmp/contigs.fa
	aln_tag = paste(paste(vvd_tmp,"/contigs.fa",sep=""), fna_tmp, ">", aln_tmp)
	aln_cmd = sprintf("exonerate %s", aln_tag)
	print(aln_cmd)
	aln_run = system(aln_cmd,intern=TRUE) #check by cat $aln_tmp
	bla_tag = paste(fna_tmp, paste(vvd_tmp,"/contigs.fa", sep=""), bla_tmp)  
	bla_cmd = sprintf("blat %s", bla_tag)
	print(bla_cmd)
	bla_run = system(bla_cmd,intern=TRUE) #check by cat $aln_tmp
	vug_tag = paste(aln_tmp, " | grep vulgar >", vug_tmp)
	vug_cmd = sprintf("cat %s",vug_tag)
	print(vug_cmd)
	vug_run = system(vug_cmd,intern=TRUE) #check by cat $vug_tmp 
	}, error=function(e){ cat("geno typing failed at STAGE1, see", vvd_tmp, "\n") } )

	### INFO: read the output and store them in some data structure
	# we need all nodes in /contigs.fa
	# for each node we need to know quality alignment from Qst-Qed on Rst-Red
	# so the data structure aln_dat is like a list of data.frame:
	### INFO: vug_tmp is a table file with following format
  #vulgar: NODE_1_length_1181_cov_22.111771 20 622 + 10:3999400-4001200 0 602 + 2974 M 602 602 
  #vulgar: Qname Qst Qed Qdir Rname Rst Red Rdir Mscore Mqual Qlen Rlen
	#DEBUG: vug_tmp="asm3/vulgar.txt"; ctg_tmp="asm3/contigs.fa"
	tryCatch( {
	if(confident==FALSE) cat("need also infer SV typing, not implemented yet!")
	vug_dat = read.table(vug_tmp)
	names(vug_dat)=c("vug","Qname","Qst","Qed","Qdir","Rname","Rst","Red","Rdir","Mscore","Mstat","Qal","Ral")
	vug_dat = vug_dat[vug_dat$Qal>pars$bam_rl,]
	ctg_dat = readDNAStringSet(ctg_tmp, format="fasta")
	cvg_dat = sapply(names(ctg_dat),function(x) { as.numeric(strsplit(x,"_")[[1]][6]) })
	len_dat = sapply(names(ctg_dat),function(x) { as.integer(strsplit(x,"_")[[1]][4]) })
	vug_dat$Qcvg = sapply(vug_dat$Qname, function(x) { round(cvg_dat[x]) })
	vug_dat$Qlen = sapply(vug_dat$Qname, function(x) { len_dat[x] })
	Q_cnt = sort(sapply(names(ctg_dat), function(x) { sum(vug_dat$Qname==x) })) #ascending
	}, error=function(e){ cat("geno typing failed at STAGE2, see", vvd_tmp, "\n") } )
	
	### TODO: aln2type()
	### RULE: always try most abundant typing first
	### here we believe: homo DEL > homo TRD > hete DEL > hete TRD
	#here are the logic to infer genotype only with confident SV type
	#assuming only high coverage homozygous and heterzygous cases
	#if DEL:
	#  if homo DEL:
	#			ideally asm produce one quality node and two quality aln for that node
	#			half aligned to the left of deleted region
	#			half aligned to the right of deleted region
	#		  let Q1st, Q1ed, R1st and R1ed be that of left most aln1 and Q2st, Q2ed, R2st, R2ed be right most aln2
	#			then the confirming criteria is: Q1ed ~= Q2st
	#			then the break points are: R1ed, R2st
	#	 if homo TRD:
	#			ideally asm produce two quality node with one and two quality aln each
	#			the three aln should spaced as Q1aln1 - Q2aln1 - Q1aln2
	#			expect second node span R1ed - R2st with full coverage
	#	 if hete:
	#			ideally asm produce three quality nodes and one quality aln for each node
	#			let it be aln1, aln2 and aln3 from the left most
	#			the confirming criteria is: Q1ed ~= Q2st, Q2ed ~= Q3st, 
	#			the cvg for Q2 about halfed for simple deletion, nonchange for translocation, multipled for translocated duplication 	
	#			then the break points are: R1ed/R2st, R2ed/R3st
	tryCatch({ 
	if(type=="DEL"){
		if(length(unique(vug_dat$Qname))==1 & nrow(vug_dat)==2){
			vug_dat = vug_dat[order(pmin(vug_dat$Qst,vug_dat$Qed)),] #order the two aln by Q
			if(abs(max(vug_dat[1,c("Qst","Qed")])-min(vug_dat[2,c("Qst","Qed")]))<=pars$ex_ovd){
				vug_dat = vug_dat[order(pmin(vug_dat$Rst,vug_dat$Red)),] #order the two aln by R
				Lbr = st+max(vug_dat[1,c("Rst","Red")]); Rbr = st+min(vug_dat[2,c("Rst","Red")]);
				geno = list(type="DEL",alt="DEL",dp=vug_dat$Qcvg[1],gt="1/1",gq="30",hq=".,.",Lbr=Lbr,Rbr=Rbr)
			}
			return(geno)
		}#INFO: if return here confirmed for homo DEL 
		if(length(unique(vug_dat$Qname))==2 & all(Q_cnt==c(1,2))){ #Q_cnt is from low to high
			vug_dat = vug_dat[order(pmin(vug_dat$Rst,vug_dat$Red)),]
			if(all(names(Q_cnt)[c(2,1,2)]==vug_dat$Qname)){ #Qname is in Q1,Q2,Q1 order
				if(abs(max(vug_dat[1,c("Rst","Red")])-min(vug_dat[2,c("Rst","Red")]))<=pars$ex_ovd){ #Q1-Q2
					if(abs(max(vug_dat[2,c("Rst","Red")])-min(vug_dat[3,c("Rst","Red")]))<=pars$ex_ovd){ #Q2-Q1
						Lbr = st+max(vug_dat[1,c("Rst","Red")]); Rbr = st+min(vug_dat[3,c("Rst","Red")]);
						geno = list(type="TRD",alt="DEL",dp=vug_dat$Qcvg[1],gt="1/1",gq="30",hq=".,.",Lbr=Lbr,Rbr=Rbr)
					}
				}
			}
			return(geno)
		}#INFO: if return here confirmed for homo TRD
		if(length(unique(vug_dat$Qname))==3 & all(Q_cnt==c(1,1,1))){ #Q_cnt is from low to high
      vug_dat = vug_dat[order(pmin(vug_dat$Rst,vug_dat$Red)),]
      if(abs(max(vug_dat[1,c("Rst","Red")])-min(vug_dat[2,c("Rst","Red")]))<=pars$ex_ovd){ #Q1-Q2
        if(abs(max(vug_dat[2,c("Rst","Red")])-min(vug_dat[3,c("Rst","Red")]))<=pars$ex_ovd){ #Q2-Q3
          Lbr = st+max(vug_dat[1,c("Rst","Red")]); Rbr = st+min(vug_dat[3,c("Rst","Red")]);
					if(vug_dat$Qcvg[2]/vug_dat$Qcvg[1]<pars$ex_clo & vug_dat$Qcvg[2]/vug_dat$Qcvg[3]<pars$ex_clo){
						geno = list(type="DEL",alt="DEL",dp=rug_dat$Qcvg[2],gt="0/1",gq="30",hq=".,.",Lbr=Lbr,Rbr=Rbr) #hete DEL
					} else if(vug_dat$Qcvg[2]/vug_dat$Qcvg[1]>pars$ex_chi & vug_dat$Qcvg[2]/vug_dat$Qcvg[3]>pars$ex_chi){
						geno = list(type="TDD",alt="DEL",dp=vug_dat$Qcvg[2],gt="0/1",gq="30",hq=".,.",Lbr=Lbr,Rbr=Rbr) #hete TDD
					} else {
						geno = list(type="TRD",alt="DEL",dp=vug_dat$Qcvg[2],gt="0/1",gq="30",hq=".,.",Lbr=Lbr,Rbr=Rbr) #hete DEL
					}
        }
      }
      return(geno)
    }#INFO: if return here confirmed for either hete DEL, TRD, TDD
	}
	}, error=function(e){ cat("geno typing DEL failed at STAGE3",vvd_tmp,"\n") } )
	#if INS:
	#	 if homo:
	#			ideally asm produce two quality node and one quality aln for each node	
	#			then the confirming is: R1ed ~= R2st
	#			then the berak points are: R1ed/R2st
	#	 if hete:
	#			ideally asm produce one quality node and one quality aln		
	#			this represents the reference allel
	#     FIXME: for now we just assume without confirming, and leave it as is
	#			FIXME: in the future we need to rely on clipping peaks within region
	tryCatch( {
	if(type=="INS"){
		if(length(unique(vug_dat$Qname))==2 & all(Q_cnt==c(1,1))){
			vug_dat = vug_dat[order(pmin(vug_dat$Rst,vug_dat$Red)),]
			if(abs(max(vug_dat[1,c("Rst","Red")])-min(vug_dat[2,c("Rst","Red")]))<=pars$ex_ovi){ #Q1-Q2
				Lbr = st+max(vug_dat[1,c("Rst","Red")]); Rbr = st+min(vug_dat[2,c("Rst","Red")]);
				geno = list(type="INS",alt="INS",dp=vug_dat$Qcvg[1],gt="1/1",gq="30",hq=".,.",Lbr=round((Lbr+Rbr)/2),Rbr=round((Lbr+Rbr)/2)) #homo INS
			}
			return(geno)
		} #INFO: if return here confirmed for homo INS
		if(length(unique(vug_dat$Qname))==1 & all(Q_cnt==c(1))){
			geno = list(type="INS",alt="INS",dp=vug_dat$Qcvg[1],gt="0/1",gq="20",hq=".,.",Lbr=data$POS[ix],Rbr=data$END[ix]) #homo INS
			return(geno)
		} #INFO: if return here confirmed for hete INS
	}
	}, error=function(e){ cat("geno typing INS failed at STAGE3",vvd_tmp,"\n") } )
	
	cat("No suitable genotyping found at STAGE4",vvd_tmp,"\n")
	file.copy(vvd_tmp, getwd(), recursive=TRUE)
	geno = list(type=type,alt=type,dp=".",gt=".",gq="0",hq=".,.",Lbr=data$POS[ix],Rbr=data$END[ix])
	return(geno)
} #INFO: if return here confirmed for nothing

list[raw_cmt,raw_bed] = read_com(bed_file, header=FALSE, stringsAsFactors=FALSE)
bed_field=c("CHROM","POS","END","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
names(raw_bed) = c(bed_field,"SAMPLE")
#DEBUG: geno_reg(1,raw_bed,bam_file,input_pars)
#DEBUG: geno_reg(2,raw_bed,bam_file,input_pars)
geno_dat=debug_sapply(debug,list(X=seq_len(nrow(raw_bed)),FUN=geno_reg,data=raw_bed,bams=bam_file,pars=input_pars))
geno_dat=as.data.frame(t(geno_dat)) #this returns matrix of strings

print(geno_dat)

geno_dat$Lbr=as.integer(geno_dat$Lbr) #convert to integer
geno_dat$Rbr=as.integer(geno_dat$Rbr) #convert to integer

#' @title geno_vcf
#'
#' @description \code{geno_vcf} \cr
#' If called, \code{geno_vcf} will return two list of sv_event structure
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
#' geno_vcf(seqname="1",
#' seqcbs_file="spX.seqcbs.txt",
#' seqcbs_parf="spX.seqcbs.par.txt",
#' seqcbs_opt="learn",
#' rg_files=c("spX.rg1.bam","spX.rg2.bam")
#' )
#'
#' @seealso \code{\link{call_lcd}} and \code{\link{call_ldr_ldl}}
geno_vcf=function(geno_dat,raw_bed,sample_id,sample_info,contig_info,ref_seq,ref_file) {
	#print(nrow(raw_bed)); print(nrow(geno_dat));
	#print(nrow(raw_bed)==nrow(geno_dat))
	stopifnot(nrow(raw_bed)==nrow(geno_dat))
  n_sv = nrow(raw_bed); bed_set=raw_bed;
	vcf_field=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")       #1-based vcf format
	bed_field=c("CHROM","POS","END","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT") #0-based bed format additionally containing the END field
	#vcf_class=c("character","integer","character","character","character","character","character","character","character","character")
	#bed_class=c("character","integer","integer","character","character","character","character","character","character","character","character")

  for(vi in seq_len(n_sv)){
		#NOTE: we need to update following of old raw_bed
		#POS -> Lbr
		#END -> Rbr
		#REF -> ref_seq(CHROM,POS)
		#INFO: END -> Rbr; SVLEN=Rbr-Lbr; DP=dp; AF=0.5 or 1.0;  
		#FORMAT -> GT:GQ:DP:HQ
		#SAMPLE -> gt,gq,dp,hq
		#only modify INFO field no introduce
		#totaly replace the FORMAT and SAMPLE filed
		#TODO: need to introduce NS and AF field in swan_join.R
		bed_set$POS[vi] = geno_dat[vi,]$Lbr	
		bed_set$END[vi] = geno_dat[vi,]$Rbr	
		info_str = bed_set$INFO[vi]
		if(substring(info_str,nchar(info_str),nchar(info_str)+1)!=";") info_str=paste(info_str,";",sep="")
		info_str = gsub(";END=[NA0-9.]+;]",paste(";END=",geno_dat[vi,]$Rbr,";",sep=""),info_str)
		info_str = gsub(";SVLEN=[NA0-9.]+;",paste(";SVLEN=",geno_dat[vi,]$Rbr-geno_dat[vi,]$Lbr,";",sep=""),info_str)
		info_str = gsub(";DP=[NA0-9.]+;",paste(";DP=",geno_dat[vi,]$dp,";",sep=""),info_str)
		if(geno_dat[vi,]$gt %in% c("0/1")) info_str = gsub(";AF=[NA0-9.]+;",paste(";AF=","0.5",";",sep=""),info_str)
		if(geno_dat[vi,]$gt %in% c("1/1")) info_str = gsub(";AF=[NA0-9.]+;",paste(";AF=","1.0",";",sep=""),info_str)
		format_str = "GT:GQ:DP"
		sample_str = paste(unlist(geno_dat[vi,][c("gt","gq","dp")]),collapse=":",sep="")
    ref_base=toString(subseq(ref_seq[[bed_set$CHROM[vi]]],geno_dat[vi,]$Lbr,geno_dat[vi,]$Lbr))
		bed_set$REF[vi] = ref_base
		bed_set$INFO[vi] = info_str
		bed_set$FORMAT[vi] = format_str
		bed_set[[names(bed_set)[ncol(bed_set)]]][vi] = sample_str
  }
	
	colnames(bed_set)[length(colnames(bed_set))] = sample_id
	vcf_set = bed_set[,-2]
  vcf_meta=paste(call_meta(
    c(meta_keys, sv_keys), c(meta_values, sv_values), ref_file, "swan_geno.R"),sample_meta(sample_info[,1], sample_info[,2], sample_info[,3], sample_info[,4]),contig_meta(contig_info[,1], contig_info[,2], contig_info[,3], contig_info[,4]),gettextf("#%s", paste(c(vcf_field,sample_id), collapse='\t')), sep='\n')
  #vcf_set=set_colClass(vcf_set,vcf_class)
  #bed_set=set_colClass(bed_set,bed_class)
  return(list(vcf_set=vcf_set,vcf_meta=vcf_meta,bed_set=bed_set))
}

list[vcf_set,vcf_meta,bed_set] = geno_vcf(geno_dat,raw_bed,sample_id,sample_info,contig_info,ref_seq,ref_file)

text_outprefix=cmd$options$outprefix
if(text_outprefix=='input') {
  out_prefix=gsub("\\.\\w+$",".geno",bed_file,perl=T)
} else {
  out_prefix=text_outprefix
}
if(savevcf) write_com(vcf_set,comment=vcf_meta,filename=paste(out_prefix,".vcf",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
write_com(bed_set,comment=NULL,filename=paste(out_prefix,".bed",sep=""),gzip=F,sep='\t',quote=F,col.names=F,row.names=F)
cat("==Info: total ",nrow(bed_set),"set calls,",paste(out_prefix,".bed",sep=""),"written!\n")
cat("==Info: swan_join safely done!", taggie(gtk),"\n")
cat("==Info: warnings:\n"); warnings()
cat("please find outputs in:",paste(out_prefix,".vcf",sep=""),"and/or",paste(out_prefix,".bed",sep=""),"\n",sep=" ")
cat("GENODONE\n")
