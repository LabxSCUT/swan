#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
gtk=proc.time()[3]
suppressMessages(library(optparse))
suppressMessages(library(swan))
gmk=get_gmk(Sys.getpid())
options(warn=2)

option_list <- list(
  make_option(c("-y", "--spY"), default="none",
              help="anch reads from contrast, [default: %default] \n"),
  make_option(c("-u", "--spXopt"), default="sup=3",
              help="options for spX clustering  , [default: %default] \n"),
  make_option(c("-v", "--spYopt"), default="sup=2",
              help="options for spY clustering  , [default: %default] \n")
)

parser <- OptionParser( 
  usage="%prog [options] ref_file spX.anch.txt", 
  option_list=option_list)
cat("-Info: vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
print(cmd$args)
if(length(cmd$args)!=2){ print_help(parser); quit(); }
ref_file=cmd$args[1]
spX_anch_file=cmd$args[2]
spX_opt=cmd$options$spXopt
spY_opt=cmd$options$spYopt
spY_anch_file=if(cmd$options$spY!="none") cmd$options$spY else NULL
clust_file=gsub(".anch.",".clust.",spX_anch_file)
ref=FaFile(ref_file)
rg=scanFa(ref, param=scanFaIndex(ref)) #this returns a DNAStringSe
chr_name=names(rg)[1:24]
seq_len=width(rg)[1:24]
spX_anch_data=read.table(spX_anch_file)
if(!is.null(spY_anch_file)){
  spY_anch_data=read.table(spY_anch_file)
  spX_anchor=as.data.frame(contrast_clust(spX_anch_data,spY_anch_data,chr_name,seq_len,parse_opt(spX_opt),parse_opt(spY_opt)))
} else {
  spX_anchor=as.data.frame(anch_clust(Lfwd=spX_anch_data[,1]=="+", Lchr=as.character(spX_anch_data[,2]), Lst=spX_anch_data[,3], Led=spX_anch_data[,4], 
                      Rfwd=spX_anch_data[,5]=="+", Rchr=as.character(spX_anch_data[,6]), Rst=spX_anch_data[,7], Red=spX_anch_data[,8], 
                      chr=chr_name, chr_size=seq_len))
}
#print(chr_name)
#print(seq_len)
#print(nrow(anch_data))
#print(anch_data)
spX_anchor[,1]=ifelse(spX_anchor[,1],"+","-")
spX_anchor[,5]=ifelse(spX_anchor[,5],"+","-")
#print(nsf(anchor))
write.table(nsf(spX_anchor),clust_file,quote=F,row.names=F,sep="\t")
cat("warnings if any:\n")
warnings()
