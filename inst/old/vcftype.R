#!/usr/bin/env Rscript
#to breakdown calls by type and methods
suppressMessages(library(optparse))
suppressMessages(library(swan))
option_list <- list()
parser=OptionParser(usage="%prog [options] vcf_file ", option_list=option_list)
args=commandArgs(trailingOnly=TRUE)
cmd=parse_args(parser,args,print_help_and_exit=TRUE,positional_arguments=TRUE)
if(length(cmd$args)!=1){ print_help(parser); quit(); }
vcf_file=cmd$args[1]

readVCF<-function(filename){
		tmp=strsplit(filename,".",fixed=TRUE)[[1]]
		suf=tmp[length(tmp)]
		if(suf=="vcf") { vcf=TRUE; shift=0 } else { vcf=FALSE; shift=1 } #either vcf or bed
    tab=read.table(filename,header=FALSE,sep="\t")
    temp=strsplit(paste(tab[,3+shift]),split=".",fixed=TRUE)
    method=rep("",nrow(tab));
    for(i in 1:length(temp)){
        method[i]=temp[[i]][length(temp[[i]])]
        if(substr(method[i],1,5)=="sclip") method[i]="sclip"
    }
    svfreq=rep(NA,nrow(tab))
    totreads=rep(NA,nrow(tab))
    for(i in 1:nrow(tab)){
        temp=strsplit(paste(tab[i,9+shift]),":",fixed=TRUE)
        temp=strsplit(temp[[1]][length(temp[[1]])],";",fixed=TRUE)
        temp=strsplit(temp[[1]],"=",fixed=TRUE)
        for(j in 1:length(temp)){
            if(temp[[j]][1]=="AA") svfreq[i] = as.numeric(temp[[j]][2])
            if(temp[[j]][1]=="DP") totreads[i] = as.numeric(temp[[j]][2])
        }
    }
    calls=data.frame(CHR=tab[,1],START=tab[,2],END=tab[,2+shift],TYPE=tab[,5+shift],METHOD=method,PURITY=svfreq,TOTREADS=totreads)
    calls
}

printStats<-function(tab,main=""){
    cat("#",main,"\n")
    eventLabel=c("#TYPE","all","sclip","lcd","bigd","disc","ldx","cigar(I)","cigar(D)","hax","cvg")
		cat(paste(eventLabel,sep="\t"),"\n")
		types=c("DEL","INS","INV","DUP","TRP") 
		methods=c("sclip","lcd","bigd","disc","ldx","ins","del","hax","cvg")
		cat("ALL","\t",nrow(tab)); for(m in methods) cat("\t",sum(tab$METHOD==m,na.rm=T),"(",pct(sum(tab$METHOD==m,na.rm=T)/(nrow(tab)+1e-20)),")"); cat("\n")
		for(type in types){
			typecnt=sum(tab$TYPE==type,na.rm=T)
			cat(type,"\t",typecnt); for(m in methods) cat("\t",sum(tab$TYPE==type&tab$METHOD==m,na.rm=T),"(",pct(sum(tab$TYPE==type&tab$METHOD==m,na.rm=T)/(typecnt+1e-20)),")"); cat("\n")
		}
		cat("NA","\t",sum(is.na(tab$TYPE))); for(m in methods) cat("\t",sum(is.na(tab$TYPE)&tab$METHOD==m,na.rm=T),"(",pct(sum(is.na(tab$TYPE)&tab$METHOD==m,na.rm=T)/(sum(is.na(tab$TYPE))+1e-20)),")"); cat("\n")
}

data=readVCF(vcf_file)
printStats(data,vcf_file)
