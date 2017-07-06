#!/usr/bin/env Rscript
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(optparse))
suppressMessages(library(robustbase))
suppressMessages(library(swan))
suppressMessages(library(hash))
suppressMessages(library(zoo))
suppressMessages(library(Cairo))
gtk=proc.time()[3]
options(warn=2)

option_list <- list(
  make_option(c("-l", "--scores"), default="lCd,lDl,lDr,cvg,lCi,lW",
              help="likelihood score list to be plotted [default %default]"),
  make_option(c("-k", "--compare"), default="cvg,isize,hang,softhang,big",
              help="compare characteristics to be plotted [default %default]"),
  make_option(c("-b", "--block"), type="integer", default=1000000,
              help="plottting block size in bp [default %default]"),
  make_option(c("-u", "--start"), type="integer", default=-1,
              help="plottting block start in bp [default %default]"),
  make_option(c("-v", "--end"), type="integer", default=-1,
              help="plottting block end in bp [default %default]"),
  make_option(c("-t", "--ticks"), default="",
              help="plottting ticks in plot x1,x2:y1,y2 ... [default %default]"),
  make_option(c("-m", "--opt"), default="track=HAF,method=value,thresh=5_track=HAR,method=value,thresh=5_track=lCd,method=empr,thresh=8_track=lDl,method=empr,thresh=8_track=lDr,method=empr,thresh=8_track=lCi,method=theo,thresh=level3_track=lW,method=theo,thresh=level3_track=cvg,method=value,thresh=30",
              help="options for score plotting: boot, robust, theo [default %default]"),
  make_option(c("-w", "--wh"), default=1.5,
              help="plotting width to height ratio [default %default]"),
  make_option(c("-c", "--chrName"), default="chr22",
              help="seqname as in spX.swan.ovrd.txt [default %default]"),
  make_option(c("-d", "--override"), default="",
              help="threshold override file, spX.swan.ovrd.txt [default %default]"),
  make_option(c("-r", "--RData"), default="input",
              help="prefix of prefix.turnk%.RData files [default %default]"),
  make_option(c("-e", "--extraRData"), default="input",
              help="prefix of prefix.RData files [default %default]"),
  make_option(c("-n", "--noTrunk"), action="store_true", default=FALSE,
              help="disable if RData not trunked [default %default]"),
  make_option(c("-p", "--platform"), default="linux",
              help="platform png utility [default %default]"),
  make_option(c("-o", "--out"), default="input",
              help="output prefix [default %default]"),
  make_option(c("-f", "--local"),action="store_true",default=FALSE,
              help="adjusting y limits by local [default %default]"),
  make_option(c("-a", "--debug"),action="store_true",default=FALSE,
              help="debug mode [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="save verbose, [default %default]")
)
parser <- OptionParser(usage = "%prog [options] scoreFiles (separated by ,)", 
                       option_list=option_list)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)

if(length(cmd$args)!=1){ print_help(parser); quit() }
text_scoreFiles=cmd$args[1]
text_override=cmd$options$override
text_spout=cmd$options$out
text_rdata=cmd$options$RData
text_edata=cmd$options$extraRData
text_scores=cmd$options$scores
text_compare=cmd$options$compare
text_ticks=cmd$options$ticks
text_opt=cmd$options$opt
verbose=cmd$options$noQuiet #true to turn on
trunked=!cmd$options$noTrunk #true to turn on
wh=cmd$options$wh
block=cmd$options$block
platform=cmd$options$platform
local=cmd$options$local
start=cmd$options$start
end=cmd$options$end
seqname=cmd$options$chrName
#debug options
#text_scoreFiles="./NA12878.chr22.16M-21M.swan.txt.gz";text_override="NA12878.chr22.16M-21M.swan.ovrd.txt";text_spout="./NA12878.chr22.16M-21M.raw.spX.chr22.5.DEL.lcd.18717741-18717941";text_rdata="./NA12878.chr22.16M-21M.bam";text_edata="input";verbose=TRUE;trunked=TRUE;wh=1.5;block=2000000;text_scores="lW,lCi,lCd,lDl,lDr";text_compare="cvg,isize,hang,softhang,big";text_ticks="18717741:18717941";platform="linux";local=FALSE;start=18715741;end=18719941;calc_method="theo";seqname="chr22";text_opt="track=lW,method=theo,thresh=level3:track=lCi,method=theo,thresh=level3:track=lCd,method=theo,thresh=level3:track=lDl,method=theo,thresh=level3:track=lDr,method=theo,thresh=level3:track=lW,method=theo,thresh=level3";
#text_scoreFiles="./NA12878.chr22.16M-21M.swan.txt.gz";text_override="NA12878.chr22.16M-21M.swan.ovrd.txt";text_spout="./NA12878.chr22.16M-21M";text_rdata="./NA12878.chr22.16M-21M.bam";text_edata="input";verbose=TRUE;trunked=TRUE;wh=1.5;block=2000000;text_scores="lW,lCi,lCd,lDl,lDr";text_compare="cvg,isize,hang,softhang,big";text_ticks="";platform="linux";local=FALSE;start=16000001;end=21000000;calc_method="theo";seqname="chr22";text_opt="track=lW,method=theo,thresh=level3:track=lCi,method=theo,thresh=level3:track=lCd,method=theo,thresh=level3:track=lDl,method=theo,thresh=level3:track=lDr,method=theo,thresh=level3";
#text_scoreFiles="VU000119.chr22.swan.txt.gz";text_override="VU000119.chr22.swan.ovrd.txt";text_spout="VU000119.chr22.raw.spX.22.177.DEL.lcd.50672591-50672821";text_rdata="VU000119.chr22.liba.bam,VU000119.chr22.libb.bam,VU000119.chr22.libc.bam,VU000119.chr22.libd.bam";text_edata="input";verbose=TRUE;trunked=TRUE;wh=1.5;block=2000000;text_scores="lW,lCi,lCd,lDl,lDr";text_compare="cvg,isize,hang,softhang,big";text_ticks="50672591:50672821";platform="linux";local=FALSE;start=50670591;end=50674821;calc_method="theo";seqname="11";text_opt="track=lW,method=theo,thresh=level3:track=lCi,method=theo,thresh=level3:track=lCd,method=theo,thresh=level3:track=lDl,method=theo,thresh=level3:track=lDr,method=theo,thresh=level3";

scoreFiles=strsplit(text_scoreFiles,split=',')[[1]] #allow multiple score file but only one par and one big file
ovrd_file=if(!text_override=="") strsplit(text_override,split=':')[[1]] else "" #always one though
thresh_level=3
pattern="(?<key>[a-z])(?<value>[-|0-9]+)"
parFile=gsub("swan.txt.gz","swan.par.txt",scoreFiles)[grepl("swan.txt.gz",scoreFiles)] #only one though
stopifnot(length(parFile)==1)
if(text_spout=="input") { prefix=gsub(".swan.par.txt","",parFile) 
} else { prefix=text_spout }
if(text_rdata=="input") { dprefix=paste(gsub(".swan.par.txt","",parFile),sep=".")
} else { dprefix=strsplit(text_rdata,split=",")[[1]] } #dprefix could be a vector
if(text_edata=="input") { eprefix=gsub(".swan.par.txt","",parFile)
} else { eprefix=strsplit(cmd$options$extraRData,split=',')[[1]] } #eprefix is only one for sample RData
scores=strsplit(text_scores,split=",")[[1]]
compares=strsplit(text_compare,split=",")[[1]]
quote=FALSE; if(quote) col=NA else col="numeric"
xticks=NULL;yticks=NULL
if(!text_ticks==""){
  xticks=as.integer(strsplit(strsplit(text_ticks,split=":")[[1]][1],split=",")[[1]])
  yticks=as.integer(strsplit(strsplit(text_ticks,split=":")[[1]][2],split=",")[[1]])
}
if(verbose) {
  cat("scoreFiles\n"); print(scoreFiles)
  cat("parFile\n"); print(parFile)
  cat("prefix\n"); print(prefix)
  cat("dprefix\n"); print(dprefix)
  cat("eprefix\n"); print(eprefix)
}

list[comment,sv_scores]=read_com(scoreFiles, comment_char="#", colClass=col)
scan_par=read.table(parFile, header=T, sep="\t", colClass=NA, comment.char="#")
scan_opt=parse_opt2(text_opt); names(scan_opt)=lapply(scan_opt,function(x){x$track})
stepsize=unique(scan_par$stepsize)[1]; width=unique(scan_par$w)[1]
trunk_size=unique(scan_par$trunk_size)[1]
scan_start=unique(scan_par$start)[1]; scan_end=unique(scan_par$end)[1]
plot_start=ifelse(start>0,start,scan_start)
plot_end=ifelse(end>0,end,scan_end)
stopifnot(plot_end>=plot_start)
win_start=seq(scan_start,scan_end,stepsize)
plot_w_start=floor((plot_start-scan_start)/stepsize)+1
plot_w_end=floor((plot_end-scan_start)/stepsize)+1
plot_t_start=floor((plot_start-scan_start)/trunk_size)+1
plot_t_end=floor((plot_end-scan_start)/trunk_size)+1
w_per_block=block/stepsize
plot_blocks=ceiling((plot_w_end-plot_w_start+1)/(w_per_block))
if(verbose) {
  cat("=Info: plot_start=",plot_start,"\n")
  cat("=Info: plot_end=",plot_end,"\n")
  cat("=Info: block=",block,"\n")
  cat("=Info: plot_w_start=",plot_w_start,"\n")
  cat("=Info: plot_w_end=",plot_w_end,"\n")
  cat("=Info: plot_t_start=",plot_t_start,"\n")
  cat("=Info: plot_t_end=",plot_t_end,"\n")
  cat("=Info: plot_blocks=",plot_blocks,"\n")
}

rg_RData=list()
for(i in seq_len(length(dprefix)))
  rg_RData[[i]]=load_RData(plot_t_start:plot_t_end,dprefix[i],trunked)
var_names=names(rg_RData[[1]])
var_classes=lapply(var_names,function(x,data) { return(class(data[[x]])) },rg_RData[[1]])
all_RData=lapply(seq_along(var_names),new_merge_RData,var_names,var_classes,rg_RData)
names(all_RData)=var_names
for(di in seq_along(all_RData)) { all_RData[[di]]=if(length(all_RData[[di]])==0) list() else all_RData[[di]][[1]] }
extra_RData=new.env()
for(efile in eprefix) load(paste(efile,"RData",sep="."),envir=extra_RData)
all_RData=c(all_RData,as.list(extra_RData))
if(verbose){
  cat("RData:\n"); print(var_names)
  cat("extraRData:\n"); print(if(length(extra_RData)==0) ls(extra_RData) else NULL)
  #print(var_classes)
  #print(length(all_RData))
  cat("=Info: scan_par="); print(scan_par)
  cat("=Info: scan_par="); print(scan_opt)
  cat("=Info: scores="); cat(scores,"\n")
  cat("=Info: compares="); cat(compares,"\n")
}
#   expected variables in all_RData
#   rPbi=rPbi,rPbo=rPbo,rM=rM,rP=rP,rMn=rMn,rPn=rPn,rSr=rSr,rSl=rSl,rHr=rHr,rHl=rHl,
#   rMHp=rMHp,rMHn=rMHn,rMDp=rMDp,rMDn=rMDn,rSCHp=rSCHp,rSCHn=rSCHn,rSCISp=rSCISp,rSCISn=rSCISn,
#   hend_MPRs=hend_MPRs,disc_MPRs=disc_MPRs,soft_lMPRs=soft_lMPRs,soft_rMPRs=soft_rMPRs,
#   hang_rMPRs=hang_rMPRs,hang_lMPRs=hang_lMPRs,
#   hang_clip=hang_clip,prop_clip=prop_clip,coverage=coverage,
#   RL=RL,isize=isize,isize_sd=isize_sd,Delta=Delta,
#   smallDel=smallDel,bigDel=bigDel

data=list(); method=list(); thresh=list(); thresh_track=list(); all_est=list() #data and thresh
for(s in scores){
  #idx=which(scan_opt$score==s)[1];
  #method[[s]]=if(is.null(scan_opt[[s]]$method) "theo" else scan_opt[[s]]$method
  thresh_track[[s]]=thresh_score(s,idx,sv_scores[[s]],scan_par,ovrd_file,scan_opt[[s]]$method,scan_opt[[s]]$thresh,seqname)
  data[[s]]=data.frame(list(x=sv_scores$start,y=inf2num(sv_scores[[s]]),z=thresh_track[[s]]))
  colnames(data[[s]]) = c("pos",s,"thresh")
	thresh[[s]]=median(thresh_track[[s]],na.rm=T)
  #thresh[[s]]=median(thresh_score(s,idx,sv_scores[[s]],scan_par,"",scan_opt[[s]]$method,scan_opt[[s]]$thresh,seqname))
  intercept=thresh[[s]];method=scan_opt[[s]][["method"]];level=thresh[[s]];type="threshold";slope=0
  all_est[[s]]=as.data.frame(list(intercept=intercept,method=method,type=type,slope=slope,level=level))
	#cat("score=",s,"\n")
	#print(all_est[[s]])
}

plot_pair = function(dat_x1,dat_y1,dat_x2,dat_y2,min_pos,max_pos,ylab,leglab,legcol,xticks,yticks){
    cat("in plot_pair ", ylab, "\n")
    stopifnot(length(dat_x1)==length(dat_y1)); #print(length(dat_x1))
    stopifnot(length(dat_x2)==length(dat_y2)); #print(length(dat_x2))
    ymax=max(c(dat_y1,dat_y2,1),na.rm=T);ymin=0;
    #if(ylab=="hang") { print(dat_x1); print(dat_y1); print(dat_x2); print(dat_y2) }
    #if(ylab=="softhang") { print(dat_x1); print(dat_y1); print(dat_x2); print(dat_y2) }
    plot(dat_x1,dat_y1,xlim=c(min_pos,max_pos),ylim=c(ymin,ceiling(ymax*1.2)),
         col=legcol[1],ylab=ylab)
    points(dat_x2,dat_y2,col=legcol[2])
    abline(v=xticks,col=legcol[1],lty=2)
    abline(v=yticks,col=legcol[2],lty=2)
    grid()
    legend(x="topright",col=legcol,pch=c(1,1),legend=leglab)
} #plot_pair(1:10,1:10,11:20,11:20,-10,30,"a",c("1","2"),c("red","blue"),10,20)

its_win_sum=function(win_start,data_x,data_y,width){
  sum(data_y[(data_x>=win_start)&(data_x<win_start+width)])
} #counting within win_start:win_start+width

for(bi in seq_len(plot_blocks)){ #bi=1   how many blocks needed to plot the whole region
	cat("plot",bi,"-th of",plot_blocks,"blocks\n")
  b_w_start=plot_w_start+(bi-1)*(w_per_block)      #start window index of this block
  b_w_end=b_w_start+(w_per_block)-1                #end window index of this block
  b_w_end=min(b_w_end,plot_w_end)
  if(plot_blocks==1){
    sig_plot=paste(prefix,"sig",sep=".")  #generate only one fig
    dat_plot=paste(prefix,"dat",sep=".")  #generate only one fig
  } else {
    sig_plot=paste(prefix,"sig",paste(nsf(sv_scores$start[b_w_start]),"_",
                   nsf(sv_scores$start[b_w_end]+stepsize-1),sep=""),sep=".")
    dat_plot=paste(prefix,"dat",paste(nsf(sv_scores$start[b_w_start]),"_",
                   nsf(sv_scores$start[b_w_end]+stepsize-1),sep=""),sep=".")
  }
  min_pos=sv_scores$start[b_w_start];max_pos=sv_scores$start[b_w_end]+stepsize-1
  tick_sel=(xticks>=min_pos)&(xticks<=max_pos)&(yticks>=min_pos)&(yticks<=max_pos)
  win_pos=seq(min_pos,max_pos,by=stepsize)
  ### plot scores fig first ###
  plots=list()
  for(l in seq(length(scores))){ #l=3
    if(verbose) cat("==Info: l=",l,"score=",scores[l],"\n")
    s=scores[l];#idx=which(scan_opt$score==s)[1];
		data[[s]][is.na(data[[s]])] = 0
		#print(dim(data[[s]][b_w_start:b_w_end,]))
    p1=ggplot(data[[s]][b_w_start:b_w_end,c(1,2)],aes_string(x="pos",y=s))+
			 geom_point()+
       geom_line(data=data[[s]][b_w_start:b_w_end,c(1,3)],aes_string(x="pos",y="thresh"),color="blue",shape=2)+ #must add data=, otherwise rais layer problem
       geom_abline(data=all_est[[s]],show_guide=T,aes(intercept=intercept,slope=slope,
                        color=as.factor(level),linetype=as.factor(method)))+ #levels
       geom_vline(xintercept=xticks[tick_sel],color="red",linetype=2)+ #sv start
       geom_vline(xintercept=yticks[tick_sel],color="blue",linetype=2)+ #sv start
       coord_cartesian(xlim=c(min_pos,max_pos))+
       theme(axis.title.y=element_text(size=rel(4),angle=90))+
       scale_y_continuous(limits=c(min(data[[s]][b_w_start:b_w_end,2:3],-1,na.rm=TRUE),max(data[[s]][b_w_start:b_w_end,2:3],thresh[[s]],na.rm=TRUE)))
    plots[[l]]=p1 #print(p1)
  }
  tryCatch({
    g=do.call("arrangeGrob",c(plots,ncol=1))
    a=ggsave(filename=paste(sig_plot,"png",sep="."),plot=g,width=14*wh,height=14,units="in")}
    ,error=function(e){
      if(verbose) { cat("==Warn: plot error, nothing in block",bi,":",min_pos,"-",max_pos,"\n"); print(e) }})
  if(!file.exists(paste(sig_plot,"png",sep="."))) next
  
  ### now we plot the RData information ###
  if(verbose) cat("==Info: min_pos", min_pos,"max_pos", max_pos, "\n")
  CairoPNG(filename=paste(dat_plot,"png",sep="."),height=14,width=14*wh,res=600,units="in")
  par(mfrow=c(length(compares),1))
  for(c in compares){ #c="isize"
    dat_x1=min_pos:max_pos;dat_x2=dat_x1;legcol=c("red","blue")
    switch(c,
      "lograt"={dat_y1=all_RData[["seqcbs_lograt"]][b_w_start:b_w_end];
            dat_x1=win_start[b_w_start:b_w_end]
            dat_y2=rep(NA,1); dat_x2=rep(NA,1);
            ylab="lograt";leglab=c("cvg log ratio","none")},
      "cvg"={dat_y1=IRanges::as.vector(IRanges::coverage(all_RData[["rPc"]]))[min_pos:max_pos];
             dat_y2=IRanges::as.vector(IRanges::coverage(all_RData[["rMc"]]))[min_pos:max_pos];
             ylab="cvg";leglab=c("plus strand","minus strand")},
      "isize"={sel=(start(all_RData[["rSCPbi"]])<=max_pos)&
                (start(all_RData[["rSCPbi"]])>=min_pos);
            dat_x1=start(all_RData[["rSCPbi"]])[sel];
            dat_y1=width(all_RData[["rSCPbi"]])[sel];
            sel=(start(all_RData[["rPbi"]])<=max_pos)&
                (start(all_RData[["rPbi"]])>=min_pos);
            dat_x2=start(all_RData[["rPbi"]])[sel];
            dat_y2=width(all_RData[["rPbi"]])[sel];
            ylab="isize";leglab=c("soft-clipped","normal")},
      "hang"={sel=(start(all_RData[["rMHp"]])<=max_pos)&
                (start(all_RData[["rMHp"]])>=min_pos);
            tab_sel=table(start(all_RData[["rMHp"]][sel]))
            dat_x1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y1=sapply(win_pos,its_win_sum,dat_x1,dat_y1,width)
            dat_x1=win_pos
            sel=(start(all_RData[["rMHn"]])<=max_pos)&
                (start(all_RData[["rMHn"]])>=min_pos);
            tab_sel=table(start(all_RData[["rMHn"]][sel]))
            dat_x2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y2=sapply(win_pos,its_win_sum,dat_x2,dat_y2,width)
            dat_x2=win_pos
            ylab="hang";leglab=c("plus strand","minus strand")},
      "big"={sel=(start(all_RData[["rPbo"]])<=max_pos)&
                (start(all_RData[["rPbo"]])>=min_pos);
            tab_sel=table(start(all_RData[["rPbo"]][sel]))
            dat_x1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y1=sapply(win_pos,its_win_sum,dat_x1,dat_y1,width)
            dat_x1=win_pos
            sel=(start(all_RData[["rSCPbo"]])<=max_pos)&
                (start(all_RData[["rSCPbo"]])>=min_pos);
            tab_sel=table(start(all_RData[["rSCPbo"]][sel]))
            dat_x2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y2=sapply(win_pos,its_win_sum,dat_x2,dat_y2,width)
            dat_x2=win_pos
            ylab="big";leglab=c("bigDel","bigDel_soft")},
      "softhang"={sel=(start(all_RData[["rSCHp"]])<=max_pos)&
                (start(all_RData[["rSCHp"]])>=min_pos);
            tab_sel=table(start(all_RData[["rSCHp"]][sel]))
            dat_x1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y1=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y1=sapply(win_pos,its_win_sum,dat_x1,dat_y1,width)
            dat_x1=win_pos
            sel=(start(all_RData[["rSCHn"]])<=max_pos)&
                (start(all_RData[["rSCHn"]])>=min_pos);
            tab_sel=table(start(all_RData[["rSCHn"]][sel]))
            dat_x2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,1])); 
            dat_y2=if(nrow(tab_sel)==0) NULL else as.integer(as.character(as.data.frame(tab_sel)[,2]));
	          dat_y2=sapply(win_pos,its_win_sum,dat_x2,dat_y2,width)
            dat_x2=win_pos
            ylab="softhang";leglab=c("plus strand","minus strand")})
    cat("subplotting...\n")
    plot_pair(dat_x1,dat_y1,dat_x2,dat_y2,min_pos,max_pos,ylab,leglab,legcol,xticks[tick_sel],yticks[tick_sel])
  }
  dev.off() #end data_plot.png
  if(verbose) cat("==Info: plotted",paste(sig_plot,"png",sep="."),"so far take ",taggie(gtk)," s\n")
  if(verbose) cat("==Info: plotted",paste(dat_plot,"png",sep="."),"so far take ",taggie(gtk)," s\n")
}
cat("-Info: plot safely done\n")
cat("warnings if any\n")
warnings(); quit();

#debug lcd false positives
sv_start=18884121;sv_end=18884461; #5
sv_start=18717741;sv_end=18717950; #5
sv_start=17048111;sv_end=17048451; #3

min_pos=sv_start-scan_par$delta;max_pos=sv_end+scan_par$delta;
srange=seq(min_pos,max_pos,stepsize);
sv_w_start=ceiling((sv_start-min_pos)/stepsize)+1
sv_w_end=ceiling((sv_end-min_pos)/stepsize)+1
sel=which((start(all_RData[["rPbi"]])<=max_pos)&(start(all_RData[["rPbi"]])>=min_pos))
rPbi=all_RData[["rPbi"]][sel]
sel=which((start(all_RData[["rSCPbi"]])<=max_pos)&(start(all_RData[["rSCPbi"]])>=min_pos))
rSCPbi=all_RData[["rSCPbi"]][sel]
winW=IRanges(start=srange,end=srange+width-1);maxInsert=scan_par$bigDel #only use positive windows to save time 
Fx = FI(0:(scan_par$delta+1000),scan_par$isize,scan_par$isize_sdR) #Fx=FI(s-x), s-x>=0, first is 0
fy_cap = fy_cap_fail                               #this is if(fy>is_cut) use is_cut for fy
fyR = fIw(0:(maxInsert+1),width,scan_par$isize,scan_par$isize_sdR)   #fIw for w=width and lCd
fyR = ifelse(fyR>fy_cap,fy_cap,fyR)
fyL = fIw(0:(maxInsert+1),-width,scan_par$isize,scan_par$isize_sdL)  #fIw for w=-width and lCi
fyL = ifelse(fyL>fy_cap,fy_cap,fyL)
n_wins=length(winW)
winW_rPb=IRanges::as.matrix(findOverlaps_within(winW,rPbi)) #crossing pair, doesn't mean corssing
table(winW_rPb[,1]); tab_max(winW_rPb[,1]) 
winW_rSCPb=IRanges::as.matrix(findOverlaps_within(winW,rSCPbi)) #crossing pair, doesn't mean corssing
table(winW_rSCPb[,1]); tab_max(winW_rSCPb[,1])

#win 212 get most hits, pos=min_pos+212*stepsize; tplt_cvg=scan_par$coverage*(scan_par$isize/(2*scan_par$rl))=
#solution expan sv_end to +scan_par$isize+3*scan_par$isize_sdR
#rational: 1st, deletion signal only if region on ref is > scan_par$isize
#  2nd, we expand to we wouldn't miss a big deletion if two reads are supporting that 
#  exp_is=which(fyR>((exp(thresh[["lCd"]]/minmpr)-(1-scan_par$r))/scan_par$r))[1]

winW_rPb[winW_rPb[,1]==205,]
rPbi[winW_rPb[winW_rPb[,1]==205,2]]
hist(width(rPbi[winW_rPb[winW_rPb[,1]==205,2]])-scan_par$isize)

