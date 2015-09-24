#!/usr/bin/env Rscript
#script to generate SVed reference with parameters
suppressMessages(library(Biostrings))
suppressMessages(library(Rsamtools))
setwd("~/Dropbox/Charlie_Nancy/analysis/test/")
debug=FALSE
seqname="11"
fa_dir = "2013Oct10_fa"
fa_name = "human_g1k_v37.fasta"
sample_size = 1
fa = FaFile(fa_name)  ## index is some.fa.gz.
rg = scanFa(fa, param=scanFaIndex(fa)) #this returns a DNAStringSet
ref_len = length(rg[[seqname]]) #[]->DNAStringSet, support width; [[]]->DNAString, support length
list[ref_head,ref_tail]=ref_ends(rg[[seqname]])
#ref_head=16050001; ref_tail=51244566; let's use 20000000, 25000000, 30000000, 35000000
ref_head=ceiling(ref_head/10000000)*10000000; block=5000000
cat("head=",ref_head,"tail=",ref_tail,"\n")
s_start=matrix(ref_head+block*seq(0,sample_size-1))
cat("s_start=", s_start, "\n")
s_label=paste("R",1:sample_size,sep="")
## start position
s_range=s_start
## SV size
w_range=c(-10000,-1000,-200,-100,-50,50,100,200,1000,10000) #-: insertion at target, +: deletion at target
## flanking size
d_range=c(50000)

#browseVignettes("Biostrings")
#orf <- readDNAStringSet(fastaFile, format="fasta")
dir.create(fa_dir, showWarnings = FALSE)
#fa = FaFile("human_g1k_v37.fasta")  ## index is some.fa.gz.

total_set=length(w_range)*length(s_label)*length(s_range[1,])*length(d_range)
print(paste("total_set=",total_set))

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

SV_set=list()
for(i in seq(1,length(w_range))){
  if(w_range[i]<0) { #insertion at target
    sv_tmp=sample(DNAString(paste(DNA_BASES, collapse="")),
                  -w_range[i],replace=TRUE)
  } else { #deletion at target
    sv_tmp=NA 
  }
  for(k in seq(1,length(s_label)))
    for(l in seq(1,length(s_range[k])))
      for(m in seq(1,length(d_range))){
          key=c("w","l","s","d")
          value=c(w_range[i],s_label[k],format(s_range[k,l],scientific=FALSE),d_range[m])
          tag=paste(paste(key,value,sep=""),collapse="_")
          cat("tag=",tag,"\n")
          ref=DNAStringSet(subseq(rg[[seqname]],s_range[k,l]-d_range[m]+1,
                        s_range[k,l]+d_range[m]),use.names=TRUE) 
          names(ref)=seqname
          sv_data=spike_sv(ref,d_range[m],w_range[i],sv_tmp) 
          #start at 1000, the 201th pos
          sv_info=paste("sv",w_range[i],format(s_range[k,l],scientific=F),sep=".")
          names(sv_data$sved)=paste(seqname,sv_info,sep="")
          names(sv_data$sv)=paste(seqname,sv_info,sep="")
          prefix=paste(file.path(fa_dir,fa_name),tag,sep=".")
          sv_save=list(w=w_range[i],g=s_label[k],
               s=s_range[k,l],d=d_range[m],ref=sv_data$ref,
               sved=sv_data$sved,sv=sv_data$sv)
          #cat("writing ",paste(prefix,".sv_save",sep="")," ...\n")
          #save(sv_save, file=paste(prefix,".sv_save",sep="")) 
          cat("writing ",paste(prefix,".fa",sep=""), "\n")
          writeXStringSet(c(sv_data$ref,sv_data$sved), 
                          file=paste(prefix,".fa",sep="")) #both
          cat("writing ",paste(prefix,".fas",sep=""), "length=", width(sv_data$sved), "\n")
          writeXStringSet(sv_data$sved, 
                          file=paste(prefix,".fas",sep="")) #sved
          cat("writing ",paste(prefix,".fna",sep=""), "length=", width(sv_data$ref), "\n")
          writeXStringSet(sv_data$ref, 
                          file=paste(prefix,".fna",sep="")) #ref
          SV_set=rbind(SV_set, sv_save)
        }
}
print("making idx ...")
system(paste("mk_idx.sh",fa_dir))
print("finished")
#substr/subseq equivalent?
#DNAStringSet [[]] get a DNAString, [] get a DNAStringSet
#deletion(insertion) marked positive(negative) in model
