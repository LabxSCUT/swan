myrequire = function(pkg, repo="CRAN", ...){
  cat("requiring package", pkg, "\n")
  tryCatch(library(pkg,character.only=T), error=function(e) {
    print(e)
    if(repo!="CRAN"){
      source("http://bioconductor.org/biocLite.R")
      biocLite(pkg,...)
    } else {
      install.packages(pkg,repo="http://cran.us.r-project.org",...)
    }
  })
  tryCatch(library(pkg,character.only=T), error=function(e) {
    print(e)
    stop(pkg," was not installed and cannot install on the fly!\n")
  })
}

for(p in c("Rsamtools","BSgenome")) myrequire(p,repo="Bioc")
#suppressMessages(library("Rsamtools"))
#suppressMessages(library("BSgenome"))
#library("BSgenome.Hsapiens.UCSC.hg19")
#library(cubature)

##################################
# General Convenience Functions
##################################

.unlist <- function (x)
{
    ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)) {
   structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

toDataFrame<-function(bam)
{
    bam <- unname(bam) # names not useful in unlisted result
    elts <- setNames(bamWhat(param), bamWhat(param))
    lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
#    DataFrame(lst)
    dat = data.frame(matrix(nrow=length(lst[[1]]),ncol=0))
    names=rep("",0)
    for(i in 1:length(lst)){
        res=try(as.vector(lst[[i]]), silent=TRUE)
        if(class(res) !="try-error"){
            dat = cbind(dat, as.vector(lst[[i]]))        
            names = c(names,names(lst)[i])
        }
    }
    names(dat) = names
    dat
}

# Convert a list of vectors to a data frame.  Non-atomic elements have to be marked in "thrownames"
# Evolved from:
# @author Jason Bryer <<jason@@bryer.org>>
# @references \url{http://stackoverflow.com/questions/4227223/r-list-to-data-frame}
as.data.frame.List <- function(x, thrownames=NULL,optional=FALSE, ...) {
    allequal <- all(unlist(lapply(x, length)) == length(x[[1]]))
    havenames <- all(unlist(lapply(x, FUN=function(x) !is.null(names(x)))))
    if(havenames) { #All the vectors in the list have names we can use
        colnames <- unique(unlist(lapply(x, names)))
        if(!is.null(thrownames)){ 
            throw = which(match(colnames, thrownames)>0)
            colnames=colnames[-throw]
        }
        df <- data.frame(matrix(unlist(lapply(x, FUN=function(x) { x[colnames] })),
                nrow=length(x), byrow=TRUE))
        names(df) <- colnames
    } else{
        stop("All elements of list must be same size and have names")    
    }
    return(df)
}


parseFlags<-function(fl){
        flags = matrix(nrow=length(fl),ncol=11)
        for(i in 1:length(fl)){
          #  if(i%% 10000 == 0) cat(i,"\n")
            x = as.integer(intToBits(as.integer(fl[i])))
            flags[i,] = x[1:11]
        }
        flags=as.data.frame(flags)
        names(flags)<-c("isPaired","isProperPair","isQueryUnmapped","isMateUnmapped","QStrand","MStrand","FirstRead","SecondRead","NotPrimary","QCFailure","Duplicate")
        flags
}


gridCount<-function(st,ed,gsize,x){
    gridst = seq(st,ed+gsize,gsize)
    xs = sort(x)
    ind = 1
    counts = rep(0,length(gridst))
    xst = findInterval(st, xs)+1
    for(i in xst:length(xs)){
        if(i %% 10000 == 0) cat("Counted",i,"out of",length(xs),"\n")
        outOfRange=FALSE
        while(xs[i]>gridst[ind+1]){
            if(ind == length(gridst)){
                outOfRange=TRUE
                break
            }
            ind = ind+1
            if(ind == length(gridst)) {
                outOfRange=TRUE
                break
            }
        }
        if(outOfRange) break
        counts[ind] = counts[ind]+1
    }
    nOutOfRange = length(xs)-i+1
    list(gridst=gridst,counts=counts,nOutOfRange=nOutOfRange)
}


getRepeatedRanges<-function(r,THRESH){
    uniq = unique(r)
    co=countOverlaps(uniq, r, type="equal")
    THRESH=5
    throw = which(co>THRESH)
    throwrds = uniq[throw]
    list(Repeated.Range.List = throwrds, Occurences=co[throw])
}


parseCigar<-function(s){
    temp=gregexpr("[MIDNSHP=X]",s)[[1]]
    stab = data.frame(matrix(nrow=length(temp),ncol=2))
    
    prevpos = 1
    for(i in 1:length(temp)){
        stab[i,1] = substr(s,temp[i],temp[i])
        stab[i,2] = as.integer(substr(s,prevpos,temp[i]-1))
        prevpos = temp[i]+1
    }
    names(stab) = c("Operation","Length")
    stab
}

plotQuality<-function(qual,...){
    plot(t(as.matrix(qual)), type="b",ylim=c(0,100), ylab="Phred Quality",xlab="Base",...)

}

InsertionFromCigar<-function(s){
    stab=parseCigar(s)
    instab = matrix(nrow=0,ncol=2)
    sel2 = which(stab$Operation=="I")
    if(length(sel2)==0) return(instab)
    
    for(i in 1:length(sel2)){
        len = stab[sel2[i],2]
        sel3 = which(stab$Operation=="M" | stab$Operation=="D")
        sel3 = sel3[sel3<sel2[i]]
        nFromPos = sum(stab$Length[sel3])
        instab = rbind(instab, c(nFromPos,len))
    }
    instab = as.data.frame(instab)
    names(instab)=c("BasesFromFirstMatch","Length")
    instab
}

DeletionFromCigar<-function(s){
    stab=parseCigar(s)
    deltab = matrix(nrow=0,ncol=2)
    sel2 = which(stab$Operation=="D")
    if(length(sel2)==0) return(deltab)
    
    for(i in 1:length(sel2)){
        len=stab[sel2[i],2]
        sel3 = which(stab$Operation=="M" | stab$Operation=="D")
        sel3 = sel3[sel3<sel2[i]]
        nFromPos = sum(stab$Length[sel3])
        deltab = rbind(deltab, c(nFromPos,len))
    }
    deltab = as.data.frame(deltab)
    names(deltab)=c("BasesFromFirstMatch","Length")
    deltab
}

SoftClipFromCigar<-function(s){
    # Returns (num_L_softclipped_bases,
    #          num_bases_between_last_L_softclipped_base_and_first_M, 
    #          num_matched_bases,
    #          num_bases_between_first_matched_base_and_first_R_softclip_pos, 
    #          num_R_softcliped_bases)
    # Some entries are "NA" if there are no soft clipped bases.
    # e.g. 69M1D29M2S -> (0,NA,98,99,2)
    #      98M2S -> (0,NA,98,98,2)
    #      68M1I29M2S -> (0,NA,97,97,2)
    
    stab=parseCigar(s)
    
    sel2 = which(stab$Operation=="M")
    if(length(sel2)==0) return(c(0L,NA,0L,NA,0L))
    nmatch = sum(stab[sel2,2])
    
    sel = which(stab$Operation == "S")    
    if(length(sel)==0) return(c(0L,NA,nmatch,NA,0L))
    if(length(sel)==1){
        # There was only one "S", is this left or right?
        if(sel==1 || (sel==2 && stab[1,1]=="H")){
            # Left clipped, no right clip.
            nLclipped = stab[sel,2]
            nRclipped = 0L
            if((sel+1)<sel2[1]){
                nfromLend = sum(stab[(sel+1):sel2[1],2])
            } else { nfromLend= 0L}
            return(c(nLclipped,nfromLend,nmatch,NA,nRclipped))
        } else {
            # No left clip, right clipped
            nLclipped = 0L
            nRclipped=stab[sel,2]
            sel3 = which(stab$Operation=="D")
            nfromRpos = sum(stab[sel3,2])+nmatch
            return(c(nLclipped,NA,nmatch,nfromRpos,nRclipped))
        }
    } 
    if(length(sel)==2){
        # There were two "S", one for left and one for right.
        nLclipped = stab[sel[1],2]
        nRclipped = stab[sel[2],2]
        if((sel[1]+1)<sel2[1]){
            nfromLend = sum(stab[(sel[1]+1):sel2[1],2])
        } else { nfromLend= 0L}
        sel3 = which(stab$Operation=="D")
        nfromRpos = sum(stab[sel3,2])+nmatch
        return(c(nLclipped,nfromLend,nmatch,nfromRpos,nRclipped))
    }
    if(length(sel)>2){
        stop("More than 2 'S' in CIGAR string!")
    }
}

##################################
# Analysis Utility Functions
##################################



kurtosis <- function(x) {  
    m4 <- mean((x-mean(x))^4) 
    kurt <- m4/(sd(x)^4)-3  
    kurt
}

skewness <-  function(x) {
    m3 <- mean((x-mean(x))^3)
    skew <- m3/(sd(x)^3)
    skew
}


simPairedEndSeqInsertion<-function(n,N,mu,sigma,R,t,w,r,p=0,pNA=1,insertType="Foreign"){

    n2 = n+w
    N2 = round(N*r)
    N1 = N-N2

    #########################
    # Simulate pairs from 
    # mutated chromosome.
    #########################
    # Simulate paired end reads from (1,...,n2).
    px = sample(n2,N2,replace=TRUE)
    islen = round(rnorm(N2, mean=mu, sd=sigma))
    py = px+islen
    while(TRUE){
        # Make sure that all (px,py) pairs are within (1,...,n2).
        id = which(py>n2)
        if(length(id)==0) break
        newpx = sample(n,length(id),replace=TRUE)
        px[id] = newpx
        py[id] = px[id]+islen[id]
    }
    # Reads within t-R,...,t+w-1 are "missing".
    ids = which(px>=(t-R) & px<t+w)
    if(insertType=="Foreign"){
        px[ids] = NA
    } else {
        px[ids] = sample(c(1:(t-R-1),(t+w):n2),length(ids),replace=TRUE)
    }
    ids = which(py>=(t-R) & py<t+w)
    if(insertType=="Foreign"){
        py[ids] = NA
    } else {
        py[ids] = sample(c(1:(t-R-1),(t+w):n2),length(ids),replace=TRUE)
    }
    # Reads within t+w,...,n2 are shifted forward by w 
    ids = which(px>=t+w)
    px[ids] = px[ids]-w
    ids = which(py>=t+w)
    py[ids] = py[ids]-w
    
    dat0=cbind(px,py)
    
    #########################
    # Simulate pairs from 
    # normal chromosome.
    #########################
    
    px = sample(n,N1,replace=TRUE)
    islen = round(rnorm(N1, mean=mu, sd=sigma))
    py = px+islen
    while(TRUE){
        # Make sure that all (px,py) pairs are within (1,...,n).
        id = which(py>n)
        if(length(id)==0) break
        newpx = sample(n,length(id),replace=TRUE)
        px[id] = newpx
        py[id] = px[id]+islen[id]
    }
    dat1=cbind(px,py)
    
    dat=rbind(dat0,dat1)
    
    if(p>0 && p<1){
        # If p>0, then p*N of the reads are improperly mapped.
        nFail = round(p*N)
        idsFail = sample(N,nFail,replace=FALSE)
        dat[idsFail,2] = sample(n,nFail, replace=TRUE)
        isFailed = rep(0,N)
        isFailed[idsFail]=1
        dat=cbind(dat,isFailed)
    }
    
    
    
    if(pNA>0 & sum(dat[,3]==1)>0){
        for(i in which(dat[,3]==1)){
            r= runif(1)
            if(r>(1-pNA/2)){
                dat[i,1] = NA
            } 
            if(r< pNA/2){
                dat[i,2] = NA
            }
        }
    }

    
    dat = cbind(dat, c(rep(1,N2),rep(0,N1)))
    dat = as.data.frame(dat)
    names(dat)<-c("PlusPosition","MinusPosition","isFailed","isFromMutated")
    
    dat
}


simPairedEndSeqDeletion<-function(n,N,mu,sigma,R,t,w,r,p=0,pNA=1){

    n2=n-w # length of mutated chromosome
    N2 = round(N*r) # number of reads from mutated chromosome
    N1 = N-N2       # number of reads from "normal" chromosome.

    #########################
    # Simulate pairs from 
    # mutated chromosome.
    #########################
    # Simulate paired end reads from (1,...,n2).
    px = sample(n2,N2,replace=TRUE)
    islen = round(rnorm(N2, mean=mu, sd=sigma))
    py = px+islen
    while(TRUE){
        # Make sure that all (px,py) pairs are within (1,...,n2).
        id = which(py>n2)
        if(length(id)==0) break
        newpx = sample(n,length(id),replace=TRUE)
        px[id] = newpx
        py[id] = px[id]+islen[id]
    }
    # Reads straddling the breakpoint (within t-R,...,t) are "missing".
    ids = which(px>=(t-R) & px<= t)
    px[ids] = NA
    ids = which(py>=(t-R) & py<= t)
    py[ids] = NA
    
    # Reads within t+1,...,n2 are shifted to the right by w 
    ids = which(px>=t+1)
    px[ids] = px[ids]+w
    ids = which(py>=t+1)
    py[ids] = py[ids]+w
    
    dat0=cbind(px,py)
    
    #########################
    # Simulate pairs from 
    # normal chromosome.
    #########################
    
    px = sample(n,N1,replace=TRUE)
    islen = round(rnorm(N1, mean=mu, sd=sigma))
    py = px+islen
    while(TRUE){
        # Make sure that all (px,py) pairs are within (1,...,n).
        id = which(py>n)
        if(length(id)==0) break
        newpx = sample(n,length(id),replace=TRUE)
        px[id] = newpx
        py[id] = px[id]+islen[id]
    }
    dat1=cbind(px,py)
    
    dat=rbind(dat0,dat1)
    
    if(p>0 && p<1){
        # If p>0, then p*N of the reads are improperly mapped.
        nFail = round(p*N)
        idsFail = sample(N,nFail,replace=FALSE)
        dat[idsFail,2] = sample(n,nFail, replace=TRUE)
        isFailed = rep(0,N)
        isFailed[idsFail]=1
        dat=cbind(dat,isFailed)
    }
    
    if(pNA>0 & sum(dat[,3]==1)>0){
        for(i in which(dat[,3]==1)){
            r= runif(1)
            if(r>(1-pNA/2)){
                dat[i,1] = NA
            } 
            if(r< pNA/2){
                dat[i,2] = NA
            }
        }
    }

    dat = cbind(dat, c(rep(1,N2),rep(0,N1)))
    dat = as.data.frame(dat)
    names(dat)<-c("PlusPosition","MinusPosition","isFailed","isFromMutated")
    
    dat
}


simPairedEndSeq<-function(n,N,mu,sigma, p=0,pNA=1){
    # Simulates a mapped paired-end read set with:
    #   n = length of genomic region mapped
    #   N = number of mapped reads
    #   mu = mean of insert length
    #   sigma = std deviation of insert length
    #   p = proportion of reads that fail to map properly.
    #   py
    
    px = sample(n,N,replace=TRUE)       # sample start positions of plus strand.
    islen = round(rnorm(N, mean=mu, sd=sigma))  # sample insert lengths.
    py = px+islen       # start position sof minus strand.

    while(TRUE){
        # Make sure that all (px,py) pairs are within (1,...,n).
        id = which(py>n)
        if(length(id)==0) break
        newpx = sample(n,length(id),replace=TRUE)
        px[id] = newpx
        py[id] = px[id]+islen[id]
    }
    
    dat=cbind(px,py)        # initial data, without improperly mapped reads.

    if(!(pNA<1 && pNA>0)){
        pNA= 1
    }    
    
    if(p>0 && p<1){
        # If p>0, then p*N of the reads are improperly mapped.
        nFail = round(p*N)
        idsFail = sample(N,nFail,replace=FALSE)
        for(i in 1:nFail){
            r = runif(1)
            if(r>(1-pNA/2)){
                dat[idsFail[i],1] = NA
            }
            if(r>0.5 && r<(1-pNA/2)){
                dat[idsFail[i],1] = sample(n,1)
            }
            if(r>pNA/2 && r<0.5){
                dat[idsFail[i],2] = sample(n,1)
            } 
            if(r<pNA/2){
                dat[idsFail[i],2] = NA
            }
        }
        isFailed = rep(0,N)
        isFailed[idsFail]=1
        dat=cbind(dat,isFailed)       
                # Last column of dat is an indicator of whether one read of that pair failed to map.
    }
    
    dat
}






# ------------------------------------------
# For debugging and exams, see rbam.r
# ------------------------------------------
DeletionScan<-function(rP,rM,minpos=NA,maxpos=NA,wins=NULL, w=NA,step=10,rs=0.5,lw,lam,verbose=TRUE){
# Parameters:
#   rP: IRanges object holding (start, end) map position of + read.
#   rM: IRanges object holding (start, end) map position of - read.  
#       rP[i] should correspond to rM[i].
#   minpos, maxpos:  Only scan within [minpos, maxpos].
#   w: window size for scan
#   step: consider start positions on the grid seq(minpos, maxpos, step).
#   rs: mixture proportion, can be a vector with multiple values
#   lw: the likelihood ratio for the insert lengths
#       for example: lw<-function(y,w,mu=0,sigma=1){exp((y-mu)*w/sigma^2 - w^2/(2*sigma^2))}
#   lam: the rate for seeing a plus or minus strand read.  Equivalent to kappa^2.

    if(is.null(wins)){
        if(is.na(w)) stop("Need to specify either wins or w.")
        if(is.na(minpos)) minpos = min(start(rP))
        if(is.na(maxpos)) maxpos = max(end(rM))
        wins = IRanges(start=seq(minpos,maxpos,step),width=w)
    } 

    if(verbose) cat("Computing deletion scores for a total of ",length(wins)," windows.\n",sep="")
    if(verbose) cat("Computing coverage within windows...\n")
    ptm=proc.time()
    nP = countOverlaps(wins,rP)
    nM = countOverlaps(wins,rM)
    elapsed=proc.time()-ptm
    rsvec = t(as.matrix(rs))
    covTerm = as.matrix(nP+nM)%*%log(1-rsvec) 
    if(verbose) cat("That took",elapsed[3],"seconds.\n\n")
    
    if(verbose) cat("Finding straddling pairs ...\n")
    ptm=proc.time()
    frags = IRanges(start=end(rP)+1,end=pmax(start(rM),end(rP)+2)-1)
    ov = findOverlaps(wins,frags,type="within")
    ov = as.matrix(ov)
    elapsed=proc.time()-ptm
    if(verbose) cat("That took",elapsed[3],"seconds.\n\n")
    
    if(verbose) cat("Computing contribution of straddling pairs...\n")
    lwY = lw(width(frags),w, mu=median(width(frags)), sigma=sd(width(frags)))
    ilTerm = matrix(nrow=length(wins), ncol=length(rs), data=0)
    ptm=proc.time()
    if(nrow(ov)>0){
        for(i in 1:nrow(ov)){
            if(i %% 1000000 == 0) cat(i,"out of",nrow(ov),"\n")
            Yterm = log((1-rs)+rs*lwY[ov[i,2]])
            ilTerm[ov[i,1],] = ilTerm[ov[i,1],]+Yterm
        }
    }
    elapsed=proc.time()-ptm  ### This seems to be the time hogger!
    if(verbose) cat("That took",elapsed[3],"seconds.\n\n")

#    cat("Computing window lambdas...\n")
#    ptm=proc.time()
#    winlam = aggregate(lam, wins, FUN=sum)
#    elapsed=proc.time()-ptm
#    cat("\nThat took",elapsed[3],"seconds.\n")
#    LLR1constant = as.matrix(winlam)%*%rsvec    # This is the non-random component for the first log-lik ratio statistic.
    
    res = list(wins=wins, covTerm=covTerm, ilTerm=ilTerm, Z=covTerm+ilTerm)
}


DeletionScanWithHangingPairs<-function(rPstart,rMstart,minpos=NA,maxpos=NA,w,R=100,step=10,r=0.5,lw,F, p=0.05,
                        Delta=NA,verbose=TRUE,doplots=FALSE,mu=NA,sigma=NA){

    if(is.na(minpos)) minpos = min(rPstart, na.rm=TRUE)
    if(is.na(maxpos)) maxpos = max(rMstart, na.rm=TRUE)
    if(is.na(Delta)) Delta = maxpos-minpos
    
    SetStar = which(is.finite(rPstart) & is.finite(rMstart) & (rMstart-rPstart<Delta)& (rMstart-rPstart>0))
    SetMinus = which(!is.finite(rMstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))
    SetPlus = which(!is.finite(rPstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))    
    
    ts = seq(minpos,maxpos,step)

    hangingMinusTerm = rep(0, length(ts))
    hangingPlusTerm = rep(0, length(ts))

    rpterm = r*(1-p)/p
    if(verbose) cat("Computing contribution of Hanging Minus Strands...\n")
    for(i in 1:length(SetMinus)){
        if(i %% 1000==0)  cat(i,"out of",length(SetMinus),"\n")
        x = rPstart[SetMinus[i]]
        tids=which(ts-x>0 & ts-x<Delta)
        hangingMinusTerm[tids] = hangingMinusTerm[tids]+ log(1+rpterm*F(ts[tids]-x))
    }
    if(verbose) cat("Computing contribution of Hanging Plus Strands...\n")
    for(i in 1:length(SetPlus)){
        if(i %% 1000==0)  cat(i,"out of",length(SetPlus),"\n")
        x = rMstart[SetPlus[i]]
        tids=which(x-(ts+w)>0 & x-(ts+w)<Delta)
        hangingPlusTerm[tids] = hangingPlusTerm[tids]+ log(1+rpterm*F(x-(ts[tids]+w)))
    }
    
    insertStar = IRanges(start=rPstart[SetStar]+R+1, end=pmax(rMstart[SetStar]-1, rPstart[SetStar]+R+1))
    tRanges = IRanges(start=ts,end=ts+w)
    straddle = as.matrix(findOverlaps(tRanges, insertStar,type="within"))
    
    if(is.na(mu) || is.na(sigma)){
        fraglens = rMstart[SetStar]-rPstart[SetStar]
        mu = median(fraglens)
        sigma = sd(fraglens)   
    }
    lwY = lw(rMstart[SetStar]-rPstart[SetStar],w, mu=mu, sigma=sigma)
        
    if(verbose) cat("Computing contribution of straddling pairs...\n")
    ilTerm = rep(0,length(ts))
    ptm=proc.time()
    if(nrow(straddle)>0){
        for(i in 1:nrow(straddle)){
            if(i %% 100000 == 0) cat(i,"out of",nrow(straddle),"\n")
            Yterm = log((1-r)+r*lwY[straddle[i,2]])
            ilTerm[straddle[i,1]] = ilTerm[straddle[i,1]]+Yterm
        }
    }
    
    
    if(doplots){
        par(mfrow=c(3,1))
        plot(ts,hangingMinusTerm)
        plot(ts,hangingPlusTerm)
        plot(ts,ilTerm)
    }
    
    res = list(ts=ts,hangingMinusTerm=hangingMinusTerm,hangingPlusTerm=hangingPlusTerm,ilTerm=ilTerm,
                mu=mu,sigma=sigma,w=w,p=p,r=r,R=R,step=step,lw=lw,F=F,Delta=Delta)

}

getPeaks<-function(Z,t,b,peaksOnly=FALSE,min.gapwidth=1){
    pass = which(Z>=b)  # note that b can be a scalar or vector.
    gridsize=t[2]-t[1]
    if(length(pass)>0){
        intervals = IRanges(start=t[pass],width=gridsize)
    } else {
        intervals = IRanges()
    }
    cat(length(pass)," positions passed threshold, reducing to intervals ...\n")
    peaks = reduce(intervals,min.gapwidth=min.gapwidth)
    if(peaksOnly){
        return(peaks)
    }
    score=rep(0,length(peaks))
    peakloc=rep(0,length(peaks))
    cat("Getting max scores and peak locs... ")

    for(i in 1:length(peaks)){
        temp = which(t>= start(peaks)[i] & t<end(peaks)[i])
        score[i] = max(Z[temp])
        peakloc[i] = t[temp[which.max(Z[temp])]]
    }
    
    # peaks is the range of the peak larger than b, score is the maximum height reached by the peak, peakloc is the location of the maximum.
    list(peaks=peaks,score=score,peakloc=peakloc) 
}


InsertionScan<-function(rPstart,rMstart,R,minpos=NA,maxpos=NA,w,step=10,r=0.5,lw,F,p=0.05,maxInsertLen=300,minInsertLen=100,verbose=TRUE,doplots=FALSE){
# Parameters:
#   rPstart: start of + read (can contain NAs)
#   rMstart: start of - read (can contain NAs)
#       rPstart[i] should correspond to rMstart[i].
#       At least one of rPstart and rMstart must be a finite value.
#   R: read length
#   minpos, maxpos:  Only scan within [minpos, maxpos].
#   w: putative size of insertion
#   step: consider insertions on the grid seq(minpos, maxpos, step).
#   rs: mixture proportion, can be a vector with multiple values
#   lw: the likelihood ratio for the insert lengths
#       for example: lw<-function(y,w,mu=0,sigma=1){exp((y-mu)*w/sigma^2 - w^2/(2*sigma^2))}
#   F: right tail of the insert length distribution.
#       for example: F<-function(x){1-pnorm(x,mean=mu,sd=sigma)}
#   p: null probability of getting an improper pair.


    SetStar = which(is.finite(rPstart) & is.finite(rMstart) & (rMstart-rPstart<maxInsertLen)& (rMstart-rPstart>0))
    SetMinus = which(!is.finite(rMstart) | (rMstart-rPstart>maxInsertLen) | (rMstart-rPstart<0))
    SetPlus = which(!is.finite(rPstart) | (rMstart-rPstart>maxInsertLen) | (rMstart-rPstart<0))    

    if(is.na(minpos)) minpos = min(rPstart, na.rm=TRUE)
    if(is.na(maxpos)) maxpos = max(rMstart, na.rm=TRUE)
    
    ts = seq(minpos,maxpos,step)
    tR = IRanges(start=ts+1,width=maxInsertLen)
    tL = IRanges(end=ts-1, width=maxInsertLen)
    tR = restrict(tR, start=minpos, end=maxpos)
    tL = restrict(tL, start=minpos, end=maxpos)
    
    hangingMinusTerm = rep(0, length(ts))
    hangingPlusTerm = rep(0, length(ts))

    rpterm = r*(1-p)/p
    if(verbose) cat("Computing contribution of Hanging Minus Strands...\n")
    for(i in 1:length(SetMinus)){
        x = rPstart[SetMinus[i]]
        tids=which(ts-x>0 & ts-x<maxInsertLen)
        hangingMinusTerm[tids] = hangingMinusTerm[tids]+ log(1+rpterm*F(ts[tids]-x))
    }
    if(verbose) cat("Computing contribution of Hanging Plus Strands...\n")
    for(i in 1:length(SetPlus)){
        x = rMstart[SetPlus[i]]
        tids=which(x-ts>0 & x-ts<maxInsertLen)
        hangingPlusTerm[tids] = hangingPlusTerm[tids]+ log(1+rpterm*F(x-ts[tids]))
    }
    
    insertStar = IRanges(start=rPstart[SetStar]+R+1, end=pmax(rMstart[SetStar]-1, rPstart[SetStar]+R))
    tRanges = IRanges(start=ts,end=ts)
    straddle = as.matrix(findOverlaps(tRanges, insertStar))
    
    SetStar2 = which(is.finite(rPstart) & is.finite(rMstart) & (rMstart-rPstart<maxInsertLen)& (rMstart-rPstart>minInsertLen))
    fraglens = rMstart[SetStar2]-rPstart[SetStar2]
    if(verbose) cat("Computing contribution of straddling pairs...\n")
    lwY = lw(rMstart[SetStar]-rPstart[SetStar],-w, mu=median(fraglens), sigma=sd(fraglens))
    ilTerm = rep(0,length(ts))
    ptm=proc.time()
    if(nrow(straddle)>0){
        for(i in 1:nrow(straddle)){
            if(i %% 1000000 == 0) cat(i,"out of",nrow(straddle),"\n")
            Yterm = log((1-r)+r*lwY[straddle[i,2]])
            ilTerm[straddle[i,1]] = ilTerm[straddle[i,1]]+Yterm
        }
    }
    
    Z = hangingMinusTerm+hangingPlusTerm+ilTerm
    
    if(doplots){
        par(mfrow=c(3,1))
        plot(ts,hangingMinusTerm)
        plot(ts,hangingPlusTerm)
        plot(ts,ilTerm)
    }
    
    res = list(ts=ts,hangingMinusTerm=hangingMinusTerm,hangingPlusTerm=hangingPlusTerm,ilTerm=ilTerm)
    res

}



#### For debugging of scanPairedEventCluster.
#T = 1e7
#nx = 5000
#ny = 5050
#w = 250
#x = sample(T,nx)
#y = sample(T, ny)
## spike in a few signals.
#nsignals = 20
#sigpos = sample(T,nsignals)
#xsize=  rpois(nsignals,5)
#ysize=  rpois(nsignals,5)
#for(i in 1:nsignals){
#    x = c(x, sigpos[i]+sample(w,xsize[i]))
#    y = c(y, sigpos[i]+sample(w,ysize[i])+w)
#}
#ord =order(sigpos)
#cbind(sigpos[ord],xsize[ord],ysize[ord])    
#mcx = 1
#mcy = 1
#hits = scanPairedEventCluster(x,y,w,mcx=mcx,mcy=mcy)
#temp = findInterval(sigpos,hits[,1])+1
#found = (hits[temp,1]-sigpos < 2*w)
#sum(found)  # verify that this = length(sigpos)
#cbind(hits[temp[found],],sigpos,xsize,ysize) # look at it.
#hits[temp[found],2] - xsize # verify that all of these are larger than 0.
#hits[temp[found],3] - ysize # verify that all of these are larger than 0.
## Seems to be working, when point process too dense, there is interference from "false" clusters".

scanPairedEventCluster<-function(x,y,w,mcx=1,mcy=1,verbose=FALSE){
# x and y are event locations, eg "+" and "-"
# locations of hanging reads.
# Look for a window of length w containing 
# at least mcx "x" events immediately followed
# by a window of length w containing at least mcy
# "y" events.  
#
# Complexity: BigO(nx+ny).

    x = sort(x)
    y = sort(y)
    nx = length(x)
    ny = length(y)
    
    wst = x[1] # Window is at [wst, wst+w-1] in x.
    xid = 1  # Index of smallest x in window.
    yid = 1 
    ydone = FALSE
    xdone = FALSE
    w2 = 2*w
    hits = matrix(nrow=0,ncol=3)  # Keep table of all clusters 
    while(TRUE){
#        if(wst>1066853) break  # For debugging
        # Find the index of smallest y larger than wst.
#        cat("Current x:",wst," Current y:",y[yid],"\n")
        while(y[yid] <= wst){
            yid = yid+1
            if(yid>ny){
                ydone = TRUE
                break
            }
        }
        if(ydone) break
        if(y[yid]-wst <= w2){
            # Found one x followed by one y 
            # within 2*w distance!
            # How many x's are in w of wst?
            if(verbose) cat("Found potential at ",wst,".\n")
            cx = 0
            while(x[xid]-wst < w){
                cx = cx+1
                xid = xid+1
                if(xid > nx){
                    xdone = TRUE
                    break
                }
            }
            # How many y's are within (x[xid]+1, wst+2*w)?
            cy = 0
            while(y[yid]-x[xid-1]<1){
                yid = yid+1
                if(yid>ny){
                    ydone = TRUE
                    break
                }
            }
            while(!ydone && y[yid]-wst<w2){
                cy = cy+1
                yid = yid+1
                if(yid>ny){
                    ydone = TRUE
                    break
                }
            }
            if(cx>=mcx && cy>=mcy) hits = rbind(hits, c(wst,cx,cy))
            wst = x[xid]  # update wst to current x.
            if(ydone) break
            if(xdone) break
        } else {
            # No y closely following this wst.  
            # Move wst to the first x >= y[yid]-2*w.
            # If this x  is smaller than y[yid], there is a cluster,
            # will skip while loop in y and find it in next round.
            # If this x is larger than y[yid], will increase yid 
            # next round.
            while(wst < y[yid]-w2){
                xid = xid+1
                if(xid> nx){
                    xdone = TRUE # none of the x's are larger than curr y-w2.
                    break
                }
                wst = x[xid]
            }
            if(xdone) break  # x ended before getting within 2*w of current y.
        }
    }
    hits= as.data.frame(hits)
    names(hits) =  c("First.X","Count.X","Count.Y")
    hits
}

getSoftClipClusters<-function(pos,cigar,seq,strand,minpos=NA,maxpos=NA,MIN.READS.PER.CLUSTER=3,MIN.BASES.PER.CLUSTER=30,computeClusters=TRUE,debug=FALSE){
    # Time complexity is O(N+n), where N is length(pos), n is maxpos-minpos.
    # Memory complexity is also O(N+n).
    # 
    # The "position of soft clipping" is defined as follows:
    # For left clippings: xxxxmmmmm, base coordinate of clipping is position of first "m"
    # For right clippings: mmmmmxxxx, base coordinate of clipping is position of fist "x".
        
    ###!#! How does strand play in here???
    if(debug) cat("Getting clusters of soft clipped reads (min cluster size = ",MIN.READS.PER.CLUSTER," reads, ", MIN.BASES.PER.CLUSTER," bases)",".  Total reads analyzed: ",length(pos),"\n",sep="")
    hasSC=grep("S", cigar, fixed=TRUE)
    pos=pos[hasSC]
    cigar=cigar[hasSC]
    seq=seq[hasSC]
    strand=strand[hasSC]
    if(debug) cat(length(hasSC)," reads have soft clipping.\n",sep="")
    
    # Parse the CIGAR string to get clipping information.
    softclip = matrix(nrow=length(pos),ncol=5)
    ptm=proc.time()
        if(length(pos)>0){
        for(i in 1:length(pos)){
            if(i %% 10000 == 0){
                elapsed = proc.time()-ptm
                if(debug) cat("\t\t",i," of ",length(pos)," analyzed, that took ",elapsed[3]," seconds.\n",sep="")
            }
            softclip[i,]=SoftClipFromCigar(cigar[i])
        }
    }
    softclip=as.data.frame(softclip)
    names(softclip) = c("nLclipped","nfromLend","nmatch","ntoRstart","nRclipped")
    
    # define the genome region.
    if(is.na(minpos) && length(pos)>0) minpos = min(pos,na.rm=TRUE)
    if(is.na(maxpos) && length(pos)>0) maxpos = max(pos,na.rm=TRUE)
    if(is.na(minpos) && length(pos)==0) minpos = 0
    if(is.na(maxpos) && length(pos)==0) maxpos = 0
    
    n = maxpos-minpos+1
    
    # For each base in region, the following keeps the number of reads clipped there, and total bases clipped.
    basecoords = c(minpos:maxpos)
    nclipreadsL = rep(0L,n)
    nclipbasesL = rep(0L,n)
    nclipreadsR = rep(0L,n)
    nclipbasesR = rep(0L,n)

    # For each entry in cigar, this keeps the exact reference base position of the clip.    
    clipposL=rep(NA, length(cigar))
    clipposR=rep(NA, length(cigar))
    
    # Get the base counts for the left clips.   
    # Position clipped is one base to the right of clipped base.
    # i.e. xxxmmmmm , position is position of first m.
    sel = which(softclip$nLclipped>0)
    if(length(sel)>0){
        for(i in 1:length(sel)){
            bcoord=pos[sel[i]]
            clipposL[sel[i]] = bcoord
            x=bcoord-minpos+1
            if(x>0 && x<=n){
                nclipreadsL[x]=nclipreadsL[x]+1L
                nclipbasesL[x]=nclipbasesL[x]+softclip$nLclipped[sel[i]]
            }
        }
    }
    
    # Get the base counts for the right clips.
    # Position clipped is defined as the position of the first clipped base.
    # i.e. mmmmmxxx , position is position of first x.
    sel = which(softclip$nRclipped>0)
    if(length(sel)>0){
        for(i in 1:length(sel)){
            bcoord=pos[sel[i]]+softclip$ntoRstart[sel[i]]
            clipposR[sel[i]] = bcoord
            x=bcoord-minpos+1
            if(x>0 && x<=n){
                nclipreadsR[x]=nclipreadsR[x]+1L
                nclipbasesR[x]=nclipbasesR[x]+softclip$nRclipped[sel[i]]
            }
        }
    }
    
    # Threshold to get L and R clusters
    Lclusterx = which(nclipreadsL>=MIN.READS.PER.CLUSTER & nclipbasesL>=MIN.BASES.PER.CLUSTER)
    LsoftclipClusters = vector("list",length(Lclusterx))
    if(debug) cat("\t\tFound ", length(Lclusterx), " left soft-clip clusters.\n",sep="")
    if(length(Lclusterx)>0 && computeClusters){
        for(i in 1:length(Lclusterx)){
            #cat("Analyzing ",i," cluster.\n")
            bcoord = Lclusterx[i]+minpos-1
            indices = which(clipposL==bcoord)
            nbases = nclipbasesL[Lclusterx[i]]
            cseqs = rep("",length(indices))
            for(j in 1:length(indices)){
                cseqs[j] = toString(subseq(seq[indices[j]], 1, softclip$nLclipped[indices[j]]))
            } 
            consmat = consensusMatrix(cseqs,startFromRight=TRUE)
            #LsoftclipClusters[[i]] = list(position=bcoord,indices=indices,nreads=length(indices),nbases=nbases,clippedseqs=cseqs,consmat=consmat)
            LsoftclipClusters[[i]] = list(position=bcoord,nreads=length(indices),nbases=nbases,consmat=consmat)
        }
    }
    Rclusterx = which(nclipreadsR>=MIN.READS.PER.CLUSTER & nclipbasesR>=MIN.BASES.PER.CLUSTER)
    RsoftclipClusters = vector("list",length(Rclusterx))
    if(debug) cat("\t\tFound ", length(Rclusterx), " right soft-clip clusters.\n",sep="")
    if(length(Rclusterx)>0 && computeClusters){
        for(i in 1:length(Rclusterx)){
            bcoord = Rclusterx[i]+minpos-1
            indices = which(clipposR==bcoord)
            nbases = nclipbasesR[Rclusterx[i]]
            cseqs = rep("",length(indices))
            for(j in 1:length(indices)){
                cseqs[j] = toString(subseq(seq[indices[j]], 
                                            width(seq[indices[j]])-softclip$nRclipped[indices[j]]+1,width(seq[indices[j]])))
            } 
            consmat = consensusMatrix(cseqs,startFromRight=FALSE)
            #RsoftclipClusters[[i]] = list(position=bcoord,indices=indices,nreads=length(indices),nbases=nbases,clippedseqs=cseqs,consmat=consmat)
            RsoftclipClusters[[i]] = list(position=bcoord,nreads=length(indices),nbases=nbases,consmat=consmat)
        }
    }
    
    list(LsoftclipClusters=LsoftclipClusters, RsoftclipClusters=RsoftclipClusters, bcoord=c(minpos:maxpos),
        nclipreadsL=nclipreadsL, nclipbasesL=nclipbasesL, 
        nclipreadsR=nclipreadsR, nclipbasesR=nclipbasesR,
        MIN.READS.PER.CLUSTER=MIN.READS.PER.CLUSTER)    
}

summarizeClusters<-function(scc,includeNBases=FALSE){
    if(includeNBases){
        sctab = matrix(nrow=length(scc),ncol=5)
    } else {
        sctab = matrix(nrow=length(scc),ncol=3)
    }
    if(length(scc)>0 && !is.null(scc[[1]]$chrom)){ 
        sctab = cbind(rep(NA,length(scc)),sctab)
    }
    if(length(scc)>0){
        for(i in 1:length(scc)){
            if(includeNBases){
                sctab[i,] = c(scc[[i]]$chrom,scc[[i]]$position,scc[[i]]$nreads,scc[[i]]$nbases,scc[[i]]$consmat$propMisMatch,nchar(scc[[i]]$consmat$consensus))
            } else {
                sctab[i,] = c(scc[[i]]$chrom,scc[[i]]$position,scc[[i]]$nreads,scc[[i]]$consmat$propMisMatch)
            }
        }
    }
    sctab=as.data.frame(sctab)
    if(includeNBases){
        if(length(scc)>0 && !is.null(scc[[1]]$chrom)){
            names(sctab) = c("Chrom","Position","nReads","nBases","ProportionMismatch","ConsensusLength")
        } else {
            names(sctab) = c("Position","nReads","nBases","ProportionMismatch","ConsensusLength")
        }
    } else {
        if(length(scc)>0 && !is.null(scc[[1]]$chrom)){
            names(sctab) = c("Chrom","Position","nReads","ProportionMismatch")
        } else {
            names(sctab) = c("Position","nReads","ProportionMismatch")
        }
    }
    if(length(scc)>0){
        sctab$Position = as.numeric(paste(sctab$Position))
        sctab$nReads=as.numeric(paste(sctab$nReads))
        if(includeNBases) sctab$nBases=as.numeric(paste(sctab$nBases))
        sctab$ProportionMismatch=as.numeric(paste(sctab$ProportionMismatch))
        if(includeNBases) sctab$ConsensusLength = as.numeric(paste(sctab$ConsensusLength))
        if(!is.null(scc[[1]]$chrom)) sctab$Chrom = paste(sctab$Chrom)
    }
    sctab
    
}

getDeletionClusters<-function(pos,cigar,strand){
    cat("Getting clusters of Deletions.  Total reads analyzed: ",length(pos),"\n")
    
    # define the genome region.
    if(is.na(minpos)) minpos = min(pos,na.rm=TRUE)
    if(is.na(maxpos)) maxpos = max(pos,na.rm=TRUE)
    n = maxpos-minpos+1

    # Parse the CIGAR string to get deletion information.
    delreads = matrix(nrow=0,ncol=3)
    for(i in 1:length(pos)){
        deltab=DeletionFromCigar(cigar[i])
        deltab[,1] = deltab[,1]+pos[i]
        delreads = rbind(delreads,cbind(deltab,rep(i,nrow(deltab))))
    }
    delreads=as.data.frame(delreads)
    names(delreads) = c("Position","Length","ReadIndex")
    
    delcount = unique(delreads[,1:2])
    count=rep(0,nrow(delcount))
    for(i in 1:nrow(delcount)){
        count[i] = sum(delreads[,1]==delcount[i,1] & delreads[,2] ==delcount[i,2])
    }
    delcount = as.data.frame(cbind(delcount,count))
    names(delcount) = c("Position","Length","nReads")
    
    list(delreads=delreads,delcount=delcount)   
}

getInsertionClusters<-function(pos,cigar,seq,strand,minpos=NA,maxpos=NA,MIN.READS.PER.CLUSTER=3){
    cat("Getting clusters of Insertions.  Total reads analyzed: ",length(pos),"\n")
    
    # define the genome region.
    if(is.na(minpos)) minpos = min(pos,na.rm=TRUE)
    if(is.na(maxpos)) maxpos = max(pos,na.rm=TRUE)
    n = maxpos-minpos+1

    # Parse the CIGAR string to get insertion information.
    insreads = matrix(nrow=0,ncol=4)
    for(i in 1:length(pos)){
        instab=InsertionFromCigar(cigar[i])
        insseq = rep("",nrow(instab))
        if(nrow(instab)>0){
            for(j in 1:nrow(instab)){
                insseq[j] = toString(subseq(seq[i], instab[j,1], instab[j,1]+instab[j,2]-1))
            }
        }
        instab[,1] = instab[,1]+pos[i]
        insreads = rbind(insreads,cbind(instab,rep(i,nrow(instab)),insseq))
    }
    insreads=as.data.frame(insreads)
    names(insreads) = c("Position","Length","ReadIndex","InsertedSequence")
    inscount = unique(insreads[,1])
    count=rep(0,length(inscount))
    for(i in 1:length(inscount)){
        count[i] = sum(insreads$Position==inscount[i])
    }
    inscount =as.matrix(cbind(inscount,count))
    names(inscount)=c("Position","nReads")
    list(insreads=insreads,inscount=inscount)
}

##################################
# Data IO functions
##################################

getRegionFromBamDir<-function(which,dirname,what=NULL,saveData=TRUE,fn="BamRegion.RData"){
# Example for which:
#   st=21675446
#   ed=25675446
#   which = RangesList("chr22"=IRanges(st,ed))
#   dirname = "../ppeichao/bamfiles_sorted/"

    if(is.null(what)) what = c("qname","rname", "flag","strand", "pos", 
    "mapq","mrnm","mpos","isize","qwidth")
    param = ScanBamParam(which=which, what=what)

    files = list.files(dirname,"*.bam.bai")
    files = unlist(strsplit(files,".bai"))
    for(i in 1:length(files)){
        bamFile = paste(dirname,files[i],sep="")
        bam = scanBam(bamFile, param=param)         # The bam file must be indexed, called ".bam.bai".
        dat0 = as.data.frame(toDataFrame(bam))
        cat("\nTotal number of reads (approx. double the number of read pairs): ", nrow(dat0),"\n\n")
        if(i==1){
            dat=dat0
        } else {
            dat = rbind(dat, dat0)
        }
    }
    ############################
    # Parse the flag field.
    ############################

    flags = matrix(nrow=nrow(dat),ncol=11)
    for(i in 1:nrow(dat)){
        x = as.integer(intToBits(as.integer(dat$flag[i])))
        flags[i,] = x[1:11]
    }
    flags=as.data.frame(flags)
    names(flags)<-c("isPaired","isProperPair","isQueryUnmapped","isMateUnmapped","QStrand","MStrand","FirstRead","SecondRead","NotPrimary","QCFailure","Duplicate")
    colSums(flags)
    
    if(saveData){
        save(dat,flags,file = fn)
    }
    list(flags=flags,dat=dat)
}



##################################
# Plotting Functions
##################################

## Got this cool little function from the IRangesOverview.pdf on bioconductor.
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
  col = "black", sep = 0.5, ...)
  {
  height <- 1
  if (is(xlim, "Ranges"))
  xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

showRegions<-function(reg,ymin, ymax,lcol=NULL,hcol=NULL){
    nr = nrow(reg)
    if(is.null(lcol)) lcol = rep("gray",nrow(reg))
    if(is.null(hcol)) hcol = rep("gray",nrow(reg))

    for(i in 1:nr){
        segments(reg[i,1],ymin,reg[i,1],ymax,col=lcol[i])
        segments(reg[i,2],ymin,reg[i,2],ymax,col=lcol[i])
        segments(reg[i,1],ymax,reg[i,2],ymax,col=hcol[i],lwd=2)
    }
}

dotplot<-function(gridst, matches, maxlen,offset,cex=0.5, horizontal=FALSE,maxSimDist=NULL,...){
    if(!horizontal){
        plot(0,0,cex=0.01,xlim=c(offset,offset+maxlen),ylim=c(offset,offset+maxlen), xlab="",ylab="",...)
        for(i in 1:length(matches)){
            points(rep(offset+gridst[i],length(matches[[i]])),offset+matches[[i]],pch=18,cex=cex)
        }
    } else {
        if(is.null(maxSimDist)) maxSimDist= maxlen
        plot(0,0,cex=0.01, xlim=c(offset,offset+maxlen), ylim=c(-maxSimDist, maxSimDist), xlab="Base 
        position", ylab="Distance to aligned region",...)
        for(i in 1:length(matches)){
            points(rep(offset+gridst[i],length(matches[[i]])), matches[[i]]-gridst[i],pch=18,cex=cex)
        }

    }
}



plotHorizontalArrow<-function(x0,x1,y,towards="right",col="red",arrheight,tipwidth,bandfr=0.8){
# plot(0,0,xlim=c(0,100),ylim=c(0,10), cex=0)
# x0=c(0,40,60,98)
# x1=c(20,50,90,100)
# y=c(4,4,4,4)
# towards=c("right","left","right","left")
# col=c("blue","green","red","purple")
# arrheight=2
# minArrowWidth= (par("xaxp")[2]-par("xaxp")[1])/25
# plotHorizontalArrow(x0,x1,y,towards,col=col,arrheight=1,tipwidth=min(min(x1-x0),minArrowWidth))
# plotHorizontalArrow(x0,x1,y=6,"right",col="red",arrheight=1,tipwidth=min(min(x1-x0),minArrowWidth))
    
    bandwidth = pmax(0.01,(x1-x0)-tipwidth)
    bandcenter=rep(0,length(x0))
    tipst = rep(0,length(x0))
    tiped = rep(0,length(x0))
    if(length(y)==1){
        y=rep(y,length(x0))
    }
    if(length(towards)==1){
        towards = rep(towards,length(x0))
    }
    if(length(col)==1){
        col = rep(col, length(x0))
    }

    ids = which(towards=="right")
    bandcenter[ids] = x0[ids]+bandwidth[ids]/2
    tipst[ids] = bandcenter[ids]+bandwidth[ids]/2
    tiped[ids]=x1[ids]

    ids = which(towards=="left")
    bandcenter[ids] = x1[ids]-bandwidth[ids]/2
    tipst[ids] = bandcenter[ids]-bandwidth[ids]/2
    tiped[ids] = x0[ids]

    for(i in 1:length(x0)){
        tip = cbind(c(tipst[i],tipst[i],tiped[i]),c(y[i]+arrheight/2, y[i]-arrheight/2,y[i]))
        filledrectangle(mid=c(bandcenter[i],y[i]),wx=bandwidth[i],wy=arrheight*bandfr, col=col[i])
        filledshape(tip, col=col[i])
    }
}


plotReadPair<-function(x0,x1,readlen,y=NULL,col="blue",height=NULL){
    # x0  is the left most position of the left read.
    # x1 is the left most position of the right read.
    # plot(0,0,xlim=c(0,100),ylim=c(0,10), cex=0)
    # x0=c(5,20,10,40,50)
    # ISL=30
    # x1=x0+rpois(length(x0),ISL)
    # plotReadPair(x0,x1,readlen=10)

    maxy = par("yaxp")[2]
    miny = par("yaxp")[1]
    n = length(x0)

    if(is.null(y)){
        gap = (maxy-miny)/n
        y = miny +  gap*c(0:(n-1))
    }

    if(is.null(height)){
        gap = (maxy-miny)/n
        height=gap*0.3
    }

    if(length(col)==1){
        col = rep(col, length(x0))
    }

    r1center = x0+readlen/2
    r2center = x1+readlen/2

    for(i in 1:length(x0)){
    filledrectangle(mid=c(r1center[i],y[i]),wx=readlen,wy=height, col=col[i])
    filledrectangle(mid=c(r2center[i],y[i]),wx=readlen,wy=height, col=col[i])
   }
    segments(r1center,y,r2center,y,col,lty=2)
}


plotCoverage<-function(co,minpos,maxpos,reg,takeSqrt){
    co2 = seqselect(co,start=minpos,end=maxpos)
    if(takeSqrt){
            plot(minpos:maxpos, sqrt(co2),
                xlab="Base position", ylab="Sqrt Coverage",pch=1, type="b")
    } else {
            plot(minpos:maxpos, sqrt(co2),
                xlab="Base position", ylab="Coverage",pch=1, type="b")
    }
}

plotReadStats<-function(dat,flags,minpos,maxpos,co,reg=NULL,
                plotMinusStrandInsertLen=FALSE, plotQuality=TRUE,
                seq=NULL,takeSqrt=FALSE, minIns=0,maxIns=NULL,maxSimDist=NULL,w=36,mm=4,lcol=NULL,hcol=NULL){
 
   nplots=1+1+plotMinusStrandInsertLen+plotQuality+(!is.null(seq))
   par(mfrow=c(nplots,1))
   if(!is.null(reg) && is.null(lcols)){
        lcols=rep("red",nrow(reg))
        hcols=rep("red",nrow(reg))
   }
   plotCoverage(co,minpos,maxpos,reg,takeSqrt)
   if(!is.null(reg)){ 
        showRegions(reg,ymin=0,ymax=par("yaxp")[2],lcol=lcols, 
            hcol=hcols)
   }
   grid()

    ####################################
    # Mapping Quality
    ####################################
    if(plotQuality){
        keep = which(dat$strand == "+"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & dat$pos>minpos)
        x = dat$pos[keep]
        y = dat$mapq[keep]
        plot(x,y,xlim=c(minpos,maxpos),xlab="read start", ylab="read quality", main="Mapping Quality")
        keep = which(dat$strand == "-"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & dat$pos>minpos)
        x = dat$pos[keep]
        y = dat$mapq[keep]
        points(x,y,col="red",pch=2)
        legend(x="topright",col=c("black","red"),pch=c(1,2),legend=c("+ strand", "- strand"))
    }
    ####################################
    # Query on + strand, insert length.
    ####################################

    keep = which(dat$strand == "+"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & 
    dat$pos>minpos)
    x = dat$pos[keep]
    if(!is.null(maxIns)){
        y = pmin(dat$isize[keep],maxIns)
    } else {
        y = dat$isize[keep]
        maxIns = max(y)
    }
    plot(x,y, xlim=c(minpos,maxpos),ylim=c(minIns,maxIns), xlab="+ read start", ylab=paste("Min(Insert 
    length,",maxIns,")",sep=""),main=paste("Insert Lengths by + Strand Read Start"))
    if(!is.null(reg)) showRegions(reg,ymin=0,ymax=maxIns,lcol=lcols[regType], hcol=lcols[regType])
    points(x,y)
    grid()

    ## Plot positions where mate is unmapped and query is on plus strand.
    keep = which(dat$strand=="+" & flags$isMateUnmapped==1 & flags$isQueryUnmapped==0 & dat$pos<maxpos & 
    dat$pos>minpos)
    points(dat$pos[keep], rep(minIns,length(keep)),col="red", pch=2)
    abline(minIns,0, col="red")

    ## Plot positions where mate is mapped but in + orientation and query is on plus strand.
    keep = which(dat$strand=="+" & flags$QStrand == flags$MStrand & flags$isMateUnmapped==0 & dat$pos<maxpos 
    & dat$pos>minpos)
    points(dat$pos[keep], rep(minIns+5,length(keep)),col="blue", pch=17)
    abline(minIns+5,0, col="blue")
    legend(x="topright", col=c("blue","red"),pch=c(2,17), legend=c("Mate unmapped", "Mate mapped in wrong 
    orientation"), cex=0.8)
    grid()

    ####################################
    # Query on - strand, insert length.
    ####################################

    if(plotMinusStrandInsertLen){
        keep = which(dat$strand == "-" & abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & 
        dat$pos>minpos)
        x = dat$pos[keep]
        y = pmin(-dat$isize[keep],maxIns)
    
        plot(x,y, xlim=c(minpos,maxpos),ylim=c(minIns,maxIns), xlab="- read start", ylab=paste("Min(Insert 
        length,",maxIns,")",sep=""),main=paste("Insert Lengths by - Strand Read Start"))
        if(!is.null(reg)) showRegions(reg,ymin=0,ymax=maxIns,lcol=lcols[regType], hcol=lcols[regType])
        points(x,y)
        grid()
    
        ## Plot positions where mate is unmapped and query is on minus strand.
        keep = which(dat$strand=="-" & flags$isMateUnmapped==1 & flags$isQueryUnmapped==0 & dat$pos<maxpos & 
        dat$pos>minpos)
        points(dat$pos[keep], rep(minIns,length(keep)),col="red", pch=2)
        abline(minIns,0, col="red")
    
        ## Plot positions where mate is mapped but in - orientation and query is on minus strand.
        keep = which(dat$strand=="-" & flags$QStrand == flags$MStrand & flags$isMateUnmapped==0 & dat$pos<maxpos 
        & dat$pos>minpos)
        points(dat$pos[keep], rep(minIns+5,length(keep)),col="blue", pch=17)
        abline(minIns+5,0, col="blue")
        legend(x="topright", col=c("blue","red"),pch=c(2,17), legend=c("Mate unmapped", "Mate mapped in wrong 
        orientation"))
    }

    ####################################
    # Plot the sequence similarity in
    # this region.
    ####################################

    if(!is.null(seq)){
        pdict=rep("",0)
        dicloc = seq(1,nchar(seq)-w,w)
        for(i in 1:length(dicloc)){
            pdict=c(pdict,substr(seq,dicloc[i],dicloc[i]+w-1))
        }
        pd = PDict(pdict, tb.start=floor(w/2)-2, tb.width=5)
        md=matchPDict(pd, DNAString(seq), max.mismatch=mm)
        matches=startIndex(md)
        dotplot(dicloc,matches,nchar(seq),offset=minpos,horizontal=TRUE,cex=1,
            maxSimDist=maxSimDist,
            main=paste("Similarity to Nearby Sequence (",w," base windows, ",mm," mismatches)",sep=""))
        grid()
    }

   
}


plotReadStats.old<-function(dat,flags,minpos,maxpos,reg=NULL,cts=NULL,gsize,recomputeCoverage=FALSE,
                plotMinusStrandInsertLen=FALSE, plotQuality=TRUE,
                seq=NULL,takeSqrt=FALSE, minIns=0,maxIns=NULL,maxSimDist=NULL,w=36,mm=4){

    nplots=1+1+plotMinusStrandInsertLen+plotQuality+(!is.null(seq))
    par(mfrow=c(nplots,1))

    if(is.null(cts)) recomputeCoverage=TRUE
    if(recomputeCoverage){
        keep=which(flags$isProperPair==1 & dat$strand=="+" & dat$pos>minpos & dat$pos<maxpos)
        cts2 = gridCount(minpos,maxpos,gsize,dat$pos[keep])
        if(takeSqrt){
            plot(cts2$gridst, sqrt(cts2$counts),
                xlab="Grid position", ylab="Sqrt Coverage",pch=1,main=paste("Chromosome",chr,"Sqrt Coverage 
                in",gsize,"Base Windows"), type="b")
            if(!is.null(reg)) showRegions(reg,ymin=0,ymax=sqrt(max(cts2$counts)),lcol=lcols[regType], 
            hcol=lcols[regType])
        } else {
            plot(cts2$gridst, cts2$counts,
                xlab="Grid position", ylab="Coverage",pch=1,main=paste("Chromosome",chr,"Coverage 
                in",gsize,"Base Windows"), type="b")
            if(!is.null(reg)) showRegions(reg,ymin=0,ymax=max(cts2$counts),lcol=lcols[regType], 
            hcol=lcols[regType])
        }
    } else {
        interval = findInterval(c(minpos,maxpos),cts$gridst)
        if(takeSqrt){
            plot(cts$gridst[interval[1]:interval[2]], sqrt(cts$counts[interval[1]:interval[2]]),
                xlab="Grid position", ylab="Sqrt Coverage",pch=1,main=paste("Chromosome",chr,"Sqrt Coverage 
                in",gsize,"Base Windows"), type="b")
            if(!is.null(reg)) 
            showRegions(reg,ymin=0,ymax=sqrt(max(cts$counts[interval[1]:interval[2]])),lcol=lcols[regType], 
            hcol=lcols[regType])
        } else {
            plot(cts$gridst[interval[1]:interval[2]], cts$counts[interval[1]:interval[2]],
                xlab="Grid position", ylab="Coverage",pch=1,main=paste("Chromosome",chr,"Coverage 
                in",gsize,"Base Windows"), type="b")
            if(!is.null(reg)) 
            showRegions(reg,ymin=0,ymax=max(cts$counts[interval[1]:interval[2]]),lcol=lcols[regType], 
            hcol=lcols[regType])
        }
    }
    grid()



    ####################################
    # Mapping Quality
    ####################################
    if(plotQuality){
        keep = which(dat$strand == "+"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & dat$pos>minpos)
        x = dat$pos[keep]
        y = dat$mapq[keep]
        plot(x,y,xlim=c(minpos,maxpos),xlab="read start", ylab="read quality", main="Mapping Quality")
        keep = which(dat$strand == "-"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & dat$pos>minpos)
        x = dat$pos[keep]
        y = dat$mapq[keep]
        points(x,y,col="red",pch=2)
        legend(x="topright",col=c("black","red"),pch=c(1,2),legend=c("+ strand", "- strand"))
    }
    ####################################
    # Query on + strand, insert length.
    ####################################

    keep = which(dat$strand == "+"& abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & 
    dat$pos>minpos)
    x = dat$pos[keep]
    if(!is.null(maxIns)){
        y = pmin(dat$isize[keep],maxIns)
    } else {
        y = dat$isize[keep]
        maxIns = max(y)
    }
    plot(x,y, ylim=c(minIns,maxIns), xlab="+ read start", ylab=paste("Min(Insert 
    length,",maxIns,")",sep=""),main=paste("Insert Lengths by + Strand Read Start"))
    if(!is.null(reg)) showRegions(reg,ymin=0,ymax=maxIns,lcol=lcols[regType], hcol=lcols[regType])
    points(x,y)
    grid()

    ## Plot positions where mate is unmapped and query is on plus strand.
    keep = which(dat$strand=="+" & flags$isMateUnmapped==1 & flags$isQueryUnmapped==0 & dat$pos<maxpos & 
    dat$pos>minpos)
    points(dat$pos[keep], rep(minIns,length(keep)),col="red", pch=2)
    abline(minIns,0, col="red")

    ## Plot positions where mate is mapped but in + orientation and query is on plus strand.
    keep = which(dat$strand=="+" & flags$QStrand == flags$MStrand & flags$isMateUnmapped==0 & dat$pos<maxpos 
    & dat$pos>minpos)
    points(dat$pos[keep], rep(minIns+5,length(keep)),col="blue", pch=17)
    abline(minIns+5,0, col="blue")
    legend(x="topright", col=c("blue","red"),pch=c(2,17), legend=c("Mate unmapped", "Mate mapped in wrong 
    orientation"), cex=0.8)
    grid()

    ####################################
    # Query on - strand, insert length.
    ####################################

    if(plotMinusStrandInsertLen){
        keep = which(dat$strand == "-" & abs(dat$isize)<100000 & (dat$mrnm == dat$rname) & dat$pos<maxpos & 
        dat$pos>minpos)
        x = dat$pos[keep]
        y = pmin(-dat$isize[keep],maxIns)
    
        plot(x,y, ylim=c(minIns,maxIns), xlab="- read start", ylab=paste("Min(Insert 
        length,",maxIns,")",sep=""),main=paste("Insert Lengths by - Strand Read Start"))
        if(!is.null(reg)) showRegions(reg,ymin=0,ymax=maxIns,lcol=lcols[regType], hcol=lcols[regType])
        points(x,y)
        grid()
    
        ## Plot positions where mate is unmapped and query is on minus strand.
        keep = which(dat$strand=="-" & flags$isMateUnmapped==1 & flags$isQueryUnmapped==0 & dat$pos<maxpos & 
        dat$pos>minpos)
        points(dat$pos[keep], rep(minIns,length(keep)),col="red", pch=2)
        abline(minIns,0, col="red")
    
        ## Plot positions where mate is mapped but in - orientation and query is on minus strand.
        keep = which(dat$strand=="-" & flags$QStrand == flags$MStrand & flags$isMateUnmapped==0 & dat$pos<maxpos 
        & dat$pos>minpos)
        points(dat$pos[keep], rep(minIns+5,length(keep)),col="blue", pch=17)
        abline(minIns+5,0, col="blue")
        legend(x="topright", col=c("blue","red"),pch=c(2,17), legend=c("Mate unmapped", "Mate mapped in wrong 
        orientation"))
    }

    ####################################
    # Plot the sequence similarity in
    # this region.
    ####################################

    if(!is.null(seq)){
        pdict=rep("",0)
        dicloc = seq(1,nchar(seq)-w,w)
        for(i in 1:length(dicloc)){
            pdict=c(pdict,substr(seq,dicloc[i],dicloc[i]+w-1))
        }
        pd = PDict(pdict, tb.start=floor(w/2)-2, tb.width=5)
        md=matchPDict(pd, DNAString(seq), max.mismatch=mm)
        matches=startIndex(md)
        dotplot(dicloc,matches,nchar(seq),offset=minpos,horizontal=TRUE,cex=1,
            maxSimDist=maxSimDist,
            main=paste("Similarity to Nearby Sequence (",w," base windows, ",mm," mismatches)",sep=""))
        grid()
    }
}

plotRegion<-function(bamfile,chr,st0,ed0,st=NA,ed=NA,mfrow=TRUE,label="",readlen=100,doCoverage=TRUE,bamfile0=NULL,doIL=TRUE,doHanging=TRUE,doSoftClip=TRUE,covymax=NA){
    # Prepare for getting the data.
    if(mfrow) par(mfrow=c(doCoverage+doIL+doHanging+doSoftClip,1))
    header=scanBamHeader(bamfile)
    chrnames=names(header[[1]]$targets)  # find out what the chromosomes are called.
    what = c("qname","rname", "flag","strand", "pos","mapq","mrnm","mpos","isize","qwidth","seq","cigar","qual")

    ### Get the data.
    which = RangesList("quack"=IRanges(st0,ed0))
    names(which) = chrnames[chr]
    param = ScanBamParam(which=which, what=what)
    bam1 = scanBam(bamfile, param=param); bam1 = bam1[[1]]
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
            bam0 = scanBam(bamfile0, param=param); bam0 = bam0[[1]]
            sel = which(!is.na(bam0$pos))
            readst=bam0$pos[sel]
            sel = which(!is.na(bam0$mpos))
            readst = c(readst,bam0$mpos[sel])
            readir = IRanges(start=readst,width=readlen)
            cov0=countOverlaps(positions,readir)
            fc=cov1/pmax(cov0,1)
            if(is.na(covymax)) covymax=max(fc,na.rm=TRUE)
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
        plot(sc$bcoord,sc$nclipreadsL,xlim=c(st0,ed0), ylim=c(0,max(c(sc$nclipreadsL,sc$nclipreadsR,1,na.rm=T))),xlab="Position",ylab="Number of Clippings",type="b",col="skyblue", lwd=2, main=mainstr,cex.main=2,cex.lab=2)
        lines(sc$bcoord,sc$nclipreadsR,col="red")
        if(!is.na(st)) segments(st,0, st, 1000, col="purple", lty=2)
        if(!is.na(ed)) segments(ed,0, ed, 1000, col="purple", lty=3)   
        grid()
        legend(x="topright", col=c("skyblue","red"),lty=c(1,1),pch=c(1,1),legend=c("Left-clipped","Right-clipped"),cex=2)
    }
    bam1
}

svlen.from.info<-function(s){
    temp=strsplit(s,split=";")[[1]]
    temp2=strsplit(temp,split="=")
    for(i in 1:length(temp2)){
        if(length(temp2[[i]])==2){
            if(temp2[[i]][1]=="SVLEN"){
                return(as.numeric(temp2[[i]][2]))
            }
        }
    }
    return(0)
}

plotInsertLengths<-function(rPstart,rMstart,minpos,maxpos,Delta=NA,CUSTOM.YLIM=TRUE,...){
    sel = which(rPstart>minpos & rPstart<maxpos)
    y = rMstart-rPstart
    if(is.na(Delta)){ Delta=max(y)} 
    y = pmin(y,Delta)
    y = pmax(y,0)
    if(CUSTOM.YLIM){
        ylim=c(0,max(y,1,na.rm=TRUE))
    } else {
        ylim=c(0,Delta)
    }
    plot(rPstart[sel],y[sel],ylim=ylim,xlab="Plus Read Start Position",ylab="Insert Length",xlim=c(minpos,maxpos),...)
}

plotHangingPlusReads<-function(rPstart,rMstart,minpos,maxpos,Delta,win=100,step=10,add=FALSE,addcol="red",addpch=1,...){
    
    sel=which(!is.finite(rPstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))    
    wins = IRanges(start=seq(minpos,maxpos,step),width=win)
    hr = IRanges(start = rMstart[sel],width=1)
    count = countOverlaps(wins,hr)
    if(add){
        points(start(wins),count, col=addcol,pch=addpch)
    } else {
        plot(start(wins),count, xlab="Window Start Position",ylab=paste("Count in",win,"base window"),xlim=c(minpos,maxpos),...)
    }
}
plotHangingMinusReads<-function(rPstart,rMstart,minpos,maxpos,Delta,win=100,step=10,add=FALSE,addcol="red",addpch=3,...){
    
    sel=which(!is.finite(rMstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))    
    wins = IRanges(start=seq(minpos,maxpos,step),width=win)
    hr = IRanges(start = rPstart[sel],width=1)
    count = countOverlaps(wins,hr)
    if(add){
        points(start(wins),count, col=addcol,pch=addpch)
    } else {
        plot(start(wins),count, xlab="Window Start Position",ylab=paste("Count in",win,"base window"),xlim=c(minpos,maxpos),...)
    }
}
plotHangingReads<-function(rPstart,rMstart,minpos,maxpos,Delta,win=100,step=10,...){
    sel=which(!is.finite(rPstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))    
    wins = IRanges(start=seq(minpos,maxpos,step),width=win)
    hr = IRanges(start = rMstart[sel],width=1)
    count = countOverlaps(wins,hr)
    plot(start(wins),count, xlab="Window Start Position",ylab=paste("Count in",win,"base window"),...)

    sel=which(!is.finite(rMstart) | (rMstart-rPstart>Delta) | (rMstart-rPstart<0))    
    wins = IRanges(start=seq(minpos,maxpos,step),width=win)
    hr = IRanges(start = rPstart[sel],width=1)
    count = countOverlaps(wins,hr)
    points(start(wins),count, col="red")
    
    legend(x="topright",col=c("black","red"),pch=c(1,1),legend=c("+ Strand Missing","- Strand Missing"))


}


addDeletions<-function(xlim,height,lwd=4,delcol="orange",verbose=FALSE){
    sel = which((delcalls[,2]>xlim[1] & delcalls[,2]<xlim[2]) | (delcalls[,3]>xlim[1] & delcalls[,3]<xlim[2]))
    for(i in 1:length(sel)){
        points(delcalls[sel[i],2], height, col=delcol, pch="|",cex=2)
        segments(delcalls[sel[i],2], height, delcalls[sel[i],3],height,lwd=lwd,col=delcol)
        points(delcalls[sel[i],3], height, col=delcol, pch=23,cex=2)

    } 
    if(verbose){
        cat("Deletions found in ",length(sel)," locations.\n")
        print(delcalls[sel,])
    }
}

addInsertions<-function(xlim,height,lwd=4,inscol="purple",verbose=FALSE){
    sel = which((inscalls[,2]>xlim[1] & inscalls[,2]<xlim[2]) | (inscalls[,3]>xlim[1] & inscalls[,3]<xlim[2]))
    for(i in 1:length(sel)){
        points(inscalls[sel[i],2], height, col=inscol, pch="|",cex=2)
        segments(inscalls[sel[i],2], height, inscalls[sel[i],3],height,lwd=lwd,col=inscol)
        points(inscalls[sel[i],3], height, col=inscol, pch=23,cex=2)

    } 
    if(verbose){
        cat("Insertions found in ",length(sel)," locations.\n")
        print(inscalls[sel,])
    }
}

consensusMatrix<-function(s,startFromRight=FALSE){
    mat = matrix(nrow = 5, ncol=max(nchar(s)), data=0)
    for(i in 1:length(s)){
        str=s[i]
        if(startFromRight) str = strReverse(str)
        c = base2num(str)
        place = c(1:length(c))
        temp=cbind(c,place)
        mat[temp] = mat[temp]+1
    }
    letters=as.factor(c("A","C","G","T","N"))
    consensus = letters[apply(mat, 2, which.max)]
    consensus=paste(consensus, collapse="")
    nseqs = apply(mat,2,sum)
    if(startFromRight) consensus = strReverse(consensus)
    propMisMatch = proportionMismatch(mat)
    mat=as.data.frame(mat, row.names=c("A","C","G","T","N"))
    #list(consensus=consensus, mat=mat, nseqs=nseqs,propMisMatch = propMisMatch)
    list(consensus=consensus, propMisMatch = propMisMatch)
}



base2num<-function(str){
    ch = strsplit(str,NULL)[[1]]
    bases = rep(0,length(ch))
    for(i in 1:length(ch)){
        if(ch[i]=="A" || ch[i]=="a"){ 
            bases[i] = 1
        } else{
            if(ch[i]=="C" || ch[i]=="c"){ 
                bases[i] = 2
            } else {
                if(ch[i]=="G" || ch[i]=="g"){
                    bases[i] = 3
                } else {
                    if(ch[i]=="T" || ch[i]=="t"){ 
                        bases[i] = 4
                    } else {
                        bases[i] = 5
                    }
                }
            }
        }        
    }
    bases

}

strReverse <- function(x)sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

entropy<-function(counts){
    p = counts/sum(counts)
    p = p[p>0]
    -sum(p*log(p))
}


proportionMismatch<-function(mat,MIN.OBS=2){
    n = ncol(mat)
    nobs = colSums(mat)
    maxchar = apply(mat,2,which.max)
    columns = c(1:n)
    maxcount = mat[cbind(maxchar,columns)]
    nmismatch = nobs-maxcount
    sel = which(nobs>=MIN.OBS)
    P = sum(nmismatch[sel])/sum(nobs[sel])
    P
}

writeVCFHeader<-function(filename){
    cat("##fileformat=VCFv4.1\n",file=filename)
    today<-Sys.Date()
    today<-format(today,"%Y%m%d")
    cat("##filedate=",today,"\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n",sep="",file=filename,append=TRUE)
    cat("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=DEL,Description=\"Deletion\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=DUP,Description=\"Duplication\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=INV,Description=\"Inversion\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=CNV,Description=\"Copy number variable region\">\n",sep="",file=filename,append=TRUE)
    cat("##ALT=<ID=TRPINV,Description=\"Transposition or Inversion\">\n",sep="",file=filename,append=TRUE)
    cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file=filename,append=TRUE)
}

writeVCF<-function(filename,chr,startpos,ref,alt,imprecise,endpos,svlen,cipos=NULL,ciend=NULL,id=NULL,qual=NULL,filter=NULL,purity=NULL,append=FALSE, somatic=TRUE){
    if(!append){
        # Write header.
        writeVCFHeader(filename)
    }
    nentries=length(chr)
    if(is.null(id)) id=rep(".",nentries)
    if(is.null(qual)) qual=rep("100",nentries)
    if(is.null(filter)) filter=rep("PASS",nentries)
    info=rep("",nentries)
    for(ei in 1:nentries){
        if(ei %%100==0) cat(ei,"\n")
        if(imprecise[ei]) info[ei] = "IMPRECISE;"
        info[ei] = paste(info[ei],"SVTYPE=",paste(alt[ei]),";",sep="")
        info[ei] = paste(info[ei],"END=",endpos[ei],";",sep="")
        info[ei] = paste(info[ei],"SVLEN=",svlen[ei],sep="")
        if(!is.null(cipos) && sum(is.na(cipos[ei,]))==0) info[ei] = paste(info[ei],";","CIPOS=",cipos[ei,1],",",cipos[ei,2],sep="")
        if(!is.null(ciend) && sum(is.na(ciend[ei,]))==0) info[ei] = paste(info[ei],";","CIEND=",ciend[ei,1],",",ciend[ei,2],sep="")    
        if(!is.null(purity)) info[ei] = paste(info[ei],";","PURITY=",purity[ei],sep="")    
        if(somatic) info[ei] = paste(info[ei],";","SOMATIC",sep="")    
    }

    for(ei in 1:nentries){
        if(ei %%100==0) cat(ei,"\n")
        cat(chr[ei],"\t",startpos[ei],"\t",id[ei],"\t",ref[ei],"\t<",paste(alt[ei]),">\t",qual[ei],"\t",filter[ei],"\t",info[ei],"\n",sep="",file=filename,append=TRUE)
    }
}

writeOneSingleMateBndVCF<-function(filename,chr1,pos1,ref1,goRight1,chr2,pos2,ref2,goRight2,ins=NULL,cipos1=NULL,cipos2=NULL,id=NULL,qual=NULL,filter=NULL,append=TRUE,somatic=TRUE){
    if(is.null(id)){
        id1 = "bnd_1"
        id2 = "bnd_2"
    } else {
        id1 = paste(id,"_1",sep="")
        id2 = paste(id,"_2",sep="")
    }
    
    if(is.null(ins)) ins=""
    
    if(goRight1 && goRight2){
        altstr1=paste("[",chr2,":",pos2,"[",ins,ref1,sep="")
        altstr2=paste("[",chr1,":",pos1,"[",ins,ref2,sep="")    
    }
    if(goRight1 && !goRight2){
        altstr1=paste("]",chr2,":",pos2,"]",ins,ref1,sep="")
        altstr2=paste(ref2,ins,"[",chr1,":",pos1,"[",sep="")
    }
    if(!goRight1 && goRight2){
        altstr1=paste(ref1,ins,"[",chr2,":",pos2,"[",sep="")
        altstr2=paste("]",chr1,":",pos1,"]",ins,ref2,sep="")
    }
    if(!goRight1 && !goRight2){
        altstr1=paste(ref1,ins,"]",chr2,":",pos2,"]",sep="")
        altstr2=paste(ref2,ins,"]",chr1,":",pos1,"]",sep="")
    }

    if(is.null(qual)){ qual1=".";qual2="."}
    if(is.null(filter)){ filter1="PASS";filter2="PASS"}

    info1=""
    if(!is.null(cipos1)) info1 = "IMPRECISE;"
    info1 = paste(info1,"SVTYPE=BND",sep="")
    info1 = paste(info1,";","MATEID=",id2,sep="")
    if(!is.null(cipos1)) info1 = paste(info1,";","CIPOS=",cipos1[1],",",cipos1[2],sep="")
    if(somatic) info1 = paste(info1,";","SOMATIC",sep="")    

    info2=""
    if(!is.null(cipos2)) info2 = "IMPRECISE;"
    info2 = paste(info2,"SVTYPE=BND",sep="")
    info2 = paste(info2,";","MATEID=",id2,sep="")
    if(!is.null(cipos2)) info2 = paste(info2,";","CIPOS=",cipos2[1],",",cipos2[2],sep="")
    if(somatic) info2 = paste(info2,";","SOMATIC",sep="")    

    cat(chr1,"\t",pos1,"\t",id1,"\t",ref1,"\t",altstr1,"\t",qual1,"\t",filter1,"\t",info1,"\n",sep="",file=filename,append=TRUE)
    cat(chr2,"\t",pos2,"\t",id2,"\t",ref2,"\t",altstr2,"\t",qual2,"\t",filter2,"\t",info2,"\n",sep="",file=filename,append=TRUE)
}



readVCF<-function(filename,header=TRUE){
    vcftab = read.table(filename,sep="\t",header=header)
    svtype=paste(vcftab[,5])
    svend = rep(NA,nrow(vcftab))
    svlen = rep(NA,nrow(vcftab))
    purity = rep(NA,nrow(vcftab))
    for(i in 1:nrow(vcftab)){
        temp=strsplit(toString(vcftab[i,8]),";")[[1]]
        temp=strsplit(temp,"=")
        for(j in 1:length(temp)){
            if(temp[[j]][1] == "END") svend[i] = as.integer(temp[[j]][2])
            if(temp[[j]][1] == "SVTYPE") svtype[i] = temp[[j]][2]
            if(temp[[j]][1] =="SVLEN") svlen[i] = as.integer(temp[[j]][2])
            if(temp[[j]][1] =="PURITY") purity[i]=as.numeric(temp[[j]][2])
        }
    }
    list(chr=vcftab[,1],start=vcftab[,2],end=svend,len=svlen,id=vcftab[,3],type=svtype,purity=purity)
}



getReadMapPositions<-function(bam,fl){
    temp=unique(bam$qname)
    mm = match(temp,bam$qname)
    sel = which(bam$strand[mm]=="+")
    rPstart1 = bam$pos[mm[sel]]
    rMstart1 = bam$mpos[mm[sel]]
    isDiscordant1 = fl$QStrand[mm[sel]] == fl$MStrand[mm[sel]]
    sel = which(bam$strand[mm]=="-")
    rMstart1 = c(rMstart1, bam$pos[mm[sel]])
    rPstart1 = c(rPstart1, bam$mpos[mm[sel]])
    isDiscordant1 = c(isDiscordant1,fl$QStrand[mm[sel]] == fl$MStrand[mm[sel]])
    sel = which(is.na(bam$strand[mm]))
    sel2=sel[which(fl$MStrand[mm[sel]]==0)]   # Mate strand is "+" but query unmapped.
    rPstart1 = c(rPstart1,bam$mpos[mm[sel2]])
    rMstart1 = c(rMstart1,rep(NA,length(sel2)))
    isDiscordant1 = c(isDiscordant1, rep(NA, length(sel2)))
    sel2 = sel[fl$MStrand[mm[sel]]==1]  # Mate strand is "-" but query unmapped.
    rPstart1 = c(rPstart1, rep(NA,length(sel2)))
    rMstart1 = c(rMstart1,bam$mpos[mm[sel2]])
    isDiscordant1 = c(isDiscordant1, rep(NA, length(sel2)))
    list(rPstart=rPstart1,rMstart=rMstart1,isDiscordant=isDiscordant1)
}



### Charlie's function to treat bimodality.

modeDecomp=function(V,CI=100,RMIN=0.001,RDIP=0.25){
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
    V = V[which(V<10*median(V))] # First do a very rough thresholding, to avoid bogus sd estimates in the next step.
    if(is.na(MAX.V)) MAX.V = median(V)+6*sd(V)  # Only consider data within 6 sd from median.
    V = V[which(V<MAX.V & V>0)]
    d=density(V, bw=bw, adjust=adjust)
    maxinwin = rep(0,length(d$x))
    for(i in 1:length(d$x)){
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

#### "nice" normal example: should get (sdL, sdR) that are about the same,  good.for.insertions should be TRUE.
#z = rnorm(200000, mean=300,sd=100)
#outliers= round(runif(2000)*1e6)    
#isize= round(c(z,outliers))
#libraryParameters(isize)
#### heavy left tail example: should get sdR that is correct,  good.for.insertions should be FALSE.
#z = rnorm(200000, mean=300,sd=100)
#outliers= round(runif(2000)*1e6)    
#left.outliers = round(runif(50000)*200)
#isize= round(c(z,outliers, left.outliers))
#libraryParameters(isize)
#### Bimodal example, should get sdR, good.for.insertions should be FALSE.
#sdtrue=30
#isize= round(c(rgamma(10000,5,10,20), rgamma(10000,8,200,10), rnorm(10000,400,sdtrue),20000,24000,15000, 2000))
#libraryParameters(isize)
#### skewed example: should get sdL, sdR, good.for.insertions should be TRUE. 
#sdR = 30
#sdL= 60
#mu = 300
#zR=rnorm(2000,mean=mu,sd=sdR); zR = zR[zR>mu]
#zL=rnorm(2000,mean=mu,sd=sdL); zL = zL[zL<mu]
#outliers= round(runif(100)*1e6)
#isize= round(c(zL,zR, outliers))
#libraryParameters(isize)

libraryParameters<-function(isize,MAX.ISIZE=NA,doplot=TRUE,method=2){

    if(method==1){
        res=modeDecomp(isize)
    } else {
        res=findModes(isize)
    }

    maxr = max(res$R) # this mode has the maximum mass.
    counter = length(res$R)
    while(res$R[counter] < maxr*0.2){
        counter = counter-1
    }
    # targetmode: assume that this was the "targeted" library size
    # Set it to be the right-most mode that is >= 20% of the size of the mode with largest mass.
    # Then, when there is contamination at the lower end, the contamination can not exceed 80%.
    # Also, when there  is contamination at the upper end, the proportion can not exceed 20%.
    targetmode = res$M[counter] 
    
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
    #   (2) In the single mode case the proportion of data less than max(100,targetmode-2*leftsdhat) is greater than 0.05.
    no.dominating.left.modes=counter==1 || (counter>1 && max(res$R[1:(counter-1)])<0.05*res$R[counter])
    not.too.skewed=sum(leftz<max(100,targetmode-2*leftsdhat))/(2*length(leftz)) < 0.05
    good.for.insertion = no.dominating.left.modes & not.too.skewed  
    
    if(doplot){
        temp=isize[which(isize<targetmode+6*sdhat)]
        h=hist(temp, 100,xlim=c(0,targetmode+6*sdhat),main=paste("Distribution of insert sizes (Useful for Insertions = ",good.for.insertion,")",sep=""), xlab="insert size", col="gray", border="gray",cex.lab=1.5)
        for(i in 1:length(res$M)){
            segments(res$M[i],0,res$M[i], max(h$counts)*2, col="red", lwd=1, lty=2)
        }
        segments(targetmode,0,targetmode, max(h$counts)*2, col="red", lwd=3, lty=1)

        segments(median(isize), 0, median(isize), max(h$counts)*2, col="green", lty=2)
        stepsize=1; x=seq(targetmode, max(isize), stepsize)
        lines(x, dnorm(x,mean=targetmode, sd=sdhat)*(2*length(z))*(h$breaks[2]-h$breaks[1]),col="cornflowerblue", lwd=3)
        stepsize=1; x=seq(0, targetmode, stepsize)
        lines(x, dnorm(x,mean=targetmode, sd=leftsdhat)*(2*length(leftz))*(h$breaks[2]-h$breaks[1]),col="orange", lwd=3)
        
       # legend(x="topright", col=c("red","red","green", "cornflowerblue", "orange"), lty=c(1,2,2,1,1),lwd=c(3,1,1,2,2), legend=c("target insert size","local mode","median", "fitted normal (right)", "fitted normal (left)"))
     legend(x="topright", col=c("red","cornflowerblue", "orange"), lty=c(1,1,1),lwd=c(3,2,2), legend=c("target insert size","fitted normal (right)", "fitted normal (left)"),cex=1.5)
    
    }
    
    list(target.isize=targetmode, sdR = sdhat, sdL = leftsdhat, modes = res$M, good.for.insertion=good.for.insertion)
}

getSubseq<-function(templateseq, rid, start, end){
  if (is(templateseq, "BWA")) {
    return(.Call("getSubseq", templateseq, rid, start, end, PACKAGE="swan"))
  } else {
    return(toString(subseq(templateseq[rid], start, end)[[1]]))
  }
}

remapSoftclipCluster<-function(cid,sc,sctab,scm,sctabm,templateseq,doLeftClip,idprefix="SC", mateidprefix="SCm", debug=FALSE){
    
    if(debug) cat("\nProcessing ",cid," out of ",length(sc)," clusters...\n")
    conseq = sc[[cid]]$consmat$consensus
    if(debug) print(sctab[cid,])
    chrnames=names(templateseq)

    # Initialize values.    
    nvcfentries=0
    nmatch=0;minDistToCluster=NA;gapseqRecipient=NA;nremap=NA;remapDistToCluster=NA;gapseqDonor=NA;
    nmatchRC=0;minDistToClusterRC=NA;gapseqRecipientRC=NA;nremapRC=NA;remapDistToClusterRC=NA;gapseqDonorRC=NA;
    bndlist=NULL;
    
    if(nchar(conseq)>SC.MIN.CONSENSUS.LENGTH){
        
        #########################################
        # Align consensus to genome.
        #########################################
     		if(debug) cat("Aligning consensus of length ",nchar(conseq)," to genome ...")
        ptm=proc.time()
        mp=vmatchPattern(conseq,templateseq,min.mismatch=0,max.mismatch=1)  
        elapsed=proc.time()-ptm; if(debug) cat("That took ",elapsed[3]," seconds.\n")
        nmatches= unlist(lapply(mp, length))
        if(debug) cat("\tNumber of matches: ",sum(nmatches),"\n")
        nmatch = sum(nmatches)
        if(sum(nmatches)>0 && sum(nmatches)<MAX.MATCHES){
            chm = which(nmatches>0)  # index of chromosomes with match.
            for(chmid in chm){  # chromosome index of match
                for(mid in 1:nmatches[chmid]){  # match number on chromosome.
                    ma = mp[[chmid]][mid]
                    ### Left clipped sequence matched to reference.
                    ### The end of mapped sequence should be 1 base position 
                    ### from the start of a right-clip cluster.  But if there is similarity
                    ### between the destination and source, the right-clip cluster may be a few
                    ### bases to the right.  The sequence between the end of mapped consensus and 
                    ### the start of the right clip cluster in destination should exactly match 
                    ### the sequence to the right of the left clip (c=left clip consensus,numeric=mapped,r=right clipped).
                    ###     source:         cccc1234567...
                    ###     destination:    cccc1234rrrrrrr
                    ### But the distance from the mapped consensus to the start of the right clip cluster
                    ### should not be longer than the length of the fragment, or else the read pair would have
                    ### mapped to the destination, not the source.
                    ### Somehow, tried getting it from the bam file fresh, didn't work (?!?)
                    selch = which(sctabm$Chrom==chrnames[chmid])
                    temp = which.min(abs(sctabm$Position[selch]-end(ma)))
                    dist = sctabm$Position[selch[temp]]-end(ma)
                    if(length(dist)>0){
                        if(is.na(minDistToCluster) || abs(dist)<abs(minDistToCluster)) minDistToCluster = dist
                        if(debug) cat("\t\tMatch on chr",chrnames[chmid],", ",start(ma),", distance to nearest cluster: ",dist,"\n",sep="")
                        if(abs(dist)<MAX.DIST.CONSENSUS.CLUSTER){
                            if(debug) cat("\t\tDestination cluster has ",sctabm$nReads[selch[temp]]," reads, ",sctabm$nBases[selch[temp]]," bases, proportion mismatch = ",sctabm$ProportionMismatch[selch[temp]],".\n",sep="")
                            if(abs(dist)>0){
                                if(end(ma)<sctabm$Position[selch[temp]]){ 
                                    #gapseq = toString(subseq(templateseq[chmid], start=end(ma)+1,end=sctabm$Position[selch[temp]])) ## gap sequence in receipient.
                                    gapseq = getSubseq(templateseq, chmid, start=end(ma)+1, end=sctabm$Position[selch[temp]]) ## gap sequence in receipient.
#                                    gapseq = getSeq(Hsapiens,names=chrnames[chmid],start=end(ma)+1,end=sctabm$Position[selch[temp]]) ## gap sequence in receipient.
                                } else {
                                    #gapseq = toString(subseq(templateseq[chmid],start=sctabm$Position[selch[temp]]+1,end=end(ma))) ## gap sequence in receipient.
                                    gapseq = getSubseq(templateseq, chmid, start=sctabm$Position[selch[temp]]+1, end=end(ma)) ## gap sequence in receipient.
#                                    gapseq = getSeq(Hsapiens,names=chrnames[chmid],start=sctabm$Position[selch[temp]]+1,end=end(ma)) ## gap sequence in receipient.
                                }
                            } else {
                                gapseq = ""
                            }
                            destclust = scm[[selch[temp]]]
                            destconseq = destclust$consmat$consensus   
                            #homeseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)[[1]]
                            homeseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)
#                            homeseq = getSeq(Hsapiens, names=chrnames[sctab$Chrom[cid]], start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)
                            remp = matchPattern(destconseq,homeseq,min.mismatch=0,max.mismatch=1)
                            nremap = length(remp)
                            if(length(remp)==1){
                                remapDistToCluster = start(remp)-REMAP.WIN
                                if(abs(start(remp)-REMAP.WIN)>0){
                                    if(start(remp)>REMAP.WIN){
                                        #gapseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid],end=sctab$Position[cid]-REMAP.WIN+start(remp)-1)
                                        gapseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid], end=sctab$Position[cid]-REMAP.WIN+start(remp)-1)
#                                        gapseq = getSeq(Hsapiens,names=chrnames[sctab$Chrom[cid]],start=sctab$Position[cid],end=sctab$Position[cid]-REMAP.WIN+start(remp)-1)
                                    } else {
                                        #gapseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid]-REMAP.WIN+start(remp)+1,end=sctab$Position[cid])
                                        gapseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid]-REMAP.WIN+start(remp)+1, end=sctab$Position[cid])
#                                        gapseq = getSeq(Hsapiens,names=chrnames[sctab$Chrom[cid]],start=sctab$Position[cid]-REMAP.WIN+start(remp)+1,end=sctab$Position[cid])
                                    }
                                } else {
                                    gapseq = ""
                                }
                                gapseqDonor=toString(gapseq)
                                if(debug) cat("\t\t\tDestination mapped to ",remapDistToCluster," from home cluster.\n",sep="")
                                if(debug) cat("\t\t\tRecipient gap:\t",gapseqRecipient, "\n",sep="")
                                if(debug) cat("\t\t\tDonor gap:\t",gapseqDonor, "\n",sep="")
                                chr1=sctab$Chrom[cid]; pos1=sctab$Position[cid]; goRight1=doLeftClip;
                                #ref1 = toString(subseq(templateseq[which(chrnames==chr1)],start=pos1,end=pos1))
                                ref1 = getSubseq(templateseq, which(chrnames==chr1), start=pos1, end=pos1)
#                                ref1 = toString(getSeq(Hsapiens,names=chrnames[chr1],start=pos1,end=pos1))
                                chr2=chrnames[chmid]; pos2=end(ma); goRight2=!doLeftClip;
                                #ref2 = toString(subseq(templateseq[which(chrnames==chr2)],start=pos2,end=pos2))
                                ref2 = getSubseq(templateseq, which(chrnames==chr2), start=pos2, end=pos2)
#                                ref2 = toString(getSeq(Hsapiens,names=chrnames[chr2],start=pos2,end=pos2))
                                nvcfentries=nvcfentries+1
                                bndlist[[nvcfentries]] = list(chr1=chr1,pos1=pos1,goRight1=goRight1,chr2=chr2,pos2=pos2,goRight2=goRight2,ins=NULL,cipos1=NULL,cipos2=NULL,cid=cid,mcid=selch[temp],id=paste(idprefix,cid,sep=""),somatic=NA,mateid=paste(mateidprefix,selch[temp],sep=""))
                                #writeOneSingleMateBndVCF(filename=vcffile,chr1=chr1,pos1=pos1,ref1=ref1,goRight1=goRight1,chr2=chr2,pos2=pos2,ref2,goRight2=goRight2,ins=NULL,cipos1=NULL,cipos2=NULL,id=paste("SC",lrstr,cid,sep=""),append=TRUE,somatic=TRUE)
                                
                            } else {
                                if(length(remp)>1){
                                    if(debug) cat("\t\tDestination consensus remapped to",length(remp),"locations in ",2*REMAP.WIN," window centered at original clipped position.\n")
                                } else {
                                    if(debug) cat("\t\tDestination consensus did not remap.\n")
                                }
                            }
                        }  ### End of if(abs(dist)<MAX.DIST.CONSENSUS.CLUSTER)
                    } ### End of if(length(dist)>0)
                }   ### End of for(mid in 1:nmatches[chmid])
            } ### End of for(chmid in chm)
        } ###  End of if(sum(nmatches)>0)
        
        #########################################
        # Align reverse complement of
        # consensus to genome.
        #########################################
        if(debug) cat("Aligning reverse complement of consensus to genome ...")
        ptm=proc.time()
        conseqrc = reverseComplement(DNAString(conseq))
        mp=vmatchPattern(conseqrc,templateseq,min.mismatch=0,max.mismatch=1)  
        elapsed=proc.time()-ptm; if(debug) cat("That took ",elapsed[3]," seconds.\n")
        nmatches= unlist(lapply(mp, length))
        if(debug) cat("\tNumber of matches: ",sum(nmatches),"\n")
        nmatchRC = sum(nmatches)
        if(sum(nmatches)>0 && sum(nmatches)<MAX.MATCHES){
            chm = which(nmatches>0)
            for(chmid in chm){
                for(mid in 1:nmatches[chmid]){
                    ma = mp[[chmid]][mid]
                    ### Reverse complement of left clipped sequence matched to reference.
                    ### The start position of the match should be a left-clip cluster.  
                    ### If there is similarity between the destination and source, the left-clip cluster 
                    ### may be a few bases to the left.  
                    selch = which(sctab$Chrom==chmid)
                    temp = which.min(abs(sctab$Position[selch]-start(ma)))
                    dist = sctab$Position[selch[temp]]-start(ma)
                    if(length(dist)>0){
                        if(is.na(minDistToClusterRC) || abs(dist)<abs(minDistToClusterRC)) minDistToClusterRC = dist
                        if(debug) cat("\t\tMatch on chr",chrnames[chmid],", ",start(ma),", distance to nearest cluster: ",dist,"\n",sep="")
                        if(abs(dist)<MAX.DIST.CONSENSUS.CLUSTER){
                            if(debug) cat("\t\tDestination cluster has ",sctab$nReads[selch[temp]]," reads, ",sctab$nBases[selch[temp]]," bases, proportion mismatch = ",sctab$ProportionMismatch[selch[temp]],".\n",sep="")
                            if(abs(dist)>0){
                                if(start(ma)>sctab$Position[selch[temp]]){ 
                                    #gapseq = toString(subseq(templateseq[chmid],start=sctab$Position[selch[temp]],end=start(ma)-1)) ## gap sequence in receipient.
                                    gapseq = getSubseq(templateseq, chmid, start=sctab$Position[selch[temp]], end=start(ma)-1) ## gap sequence in receipient.
#                                    gapseq = getSeq(Hsapiens,names=chrnames[chmid],start=sctab$Position[selch[temp]],end=start(ma)-1) ## gap sequence in receipient.
                                } else {
                                    #gapseq = toString(subseq(templateseq[chmid],start=start(ma)+1,end=sctab$Position[selch[temp]])) ## gap sequence in receipient.
                                    gapseq = getSubseq(templateseq, chmid, start=start(ma)+1, end=sctab$Position[selch[temp]]) ## gap sequence in receipient.
#                                    gapseq = getSeq(Hsapiens,names=chrnames[chmid],start=start(ma)+1,end=sctab$Position[selch[temp]]) ## gap sequence in receipient.
                                }
                            } else {
                                gapseq = ""
                            }
                            gapseqRecipientRC = toString(gapseq)
                            destclust = sc[[selch[temp]]]
                            destconseq = reverseComplement(DNAString(destclust$consmat$consensus))   
                            #homeseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)[[1]]
                            homeseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)
#                            homeseq = getSeq(Hsapiens, names=chrnames[sctab$Chrom[cid]], start=sctab$Position[cid]-REMAP.WIN, end=sctab$Position[cid]+REMAP.WIN)
                            remp = matchPattern(destconseq,homeseq,min.mismatch=0,max.mismatch=1)
                            nremapRC = length(remp)
                            if(length(remp)==1){
                                remapDistToClusterRC = start(remp)-REMAP.WIN-1
                                if(abs(start(remp)-REMAP.WIN-1)>0){
                                    if(start(remp)>(REMAP.WIN+1)){
                                        #gapseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid],end=sctab$Position[cid]-REMAP.WIN+start(remp)-2)
                                        gapseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid], end=sctab$Position[cid]-REMAP.WIN+start(remp)-2)
#                                        gapseq = getSeq(Hsapiens,names=chrnames[sctab$Chrom[cid]],start=sctab$Position[cid],end=sctab$Position[cid]-REMAP.WIN+start(remp)-2)  
                                    } else {
                                        #gapseq = subseq(templateseq[which(chrnames==sctab$Chrom[cid])],start=sctab$Position[cid]-(REMAP.WIN+1-start(remp)),end=sctab$Position[cid])
                                        gapseq = getSubseq(templateseq, which(chrnames==sctab$Chrom[cid]), start=sctab$Position[cid]-(REMAP.WIN+1-start(remp)), end=sctab$Position[cid])
#                                        gapseq = getSeq(Hsapiens,names=chrnames[sctab$Chrom[cid]],start=sctab$Position[cid]-(REMAP.WIN+1-start(remp)),end=sctab$Position[cid])
                                    }
                                } else {
                                    gapseq = ""
                                }
                                gapseqDonorRC=toString(gapseq)
                             		if(debug) cat("\t\t\tDestination mapped to ",remapDistToClusterRC," from home cluster.\n",sep="")
                                if(debug) cat("\t\t\tRecipient gap:\t",gapseqRecipientRC, "\n",sep="")
                                if(debug) cat("\t\t\tDonor gap:\t",gapseqDonorRC, "\n",sep="")
                                chr1=sctab$Chrom[cid]; pos1=sctab$Position[cid]; goRight1=doLeftClip;
                                #ref1 = toString(subseq(templateseq[which(chrnames==chr1)],start=pos1,end=pos1))
                                ref1 = getSubseq(templateseq, which(chrnames==chr1), start=pos1, end=pos1)
#                                ref1 = toString(getSeq(Hsapiens,names=chrnames[chr1],start=pos1,end=pos1))
                                chr2=chrnames[chmid]; pos2=start(ma); goRight2=doLeftClip;
                                #ref2 = toString(subseq(templateseq[which(chrnames==chr2)],start=pos2,end=pos2))
                                ref2 = getSubseq(templateseq, which(chrnames==chr2), start=pos2, end=pos2)
#                                ref2 = toString(getSeq(Hsapiens,names=chrnames[chr2],start=pos2,end=pos2))
                                nvcfentries=nvcfentries+1
                                bndlist[[nvcfentries]] = list(chr1=chr1,pos1=pos1,goRight1=goRight1,chr2=chr2,pos2=pos2,goRight2=goRight2,ins=NULL,cipos1=NULL,cipos2=NULL,cid=cid,mcid=selch[temp],id=paste(idprefix,cid,"RC",sep=""),somatic=NA, mateid=paste(idprefix,selch[temp],"RC",sep=""))
                                # writeOneSingleMateBndVCF(filename=vcffile,chr1=chr1,pos1=pos1,ref1=ref1,goRight1=goRight1,chr2=chr2,pos2=pos2,ref2,goRight2=goRight2,ins=NULL,cipos1=NULL,cipos2=NULL,id=paste("SC",lrstr,cid,"RC",sep=""),append=TRUE,somatic=TRUE)
                                
                            } else {
                                if(length(remp)>1){
                                    if(debug) cat("\t\tDestination consensus remapped to ",length(remp)," locations in ",2*REMAP.WIN," window centered at original clipped position.\n")
                                } else {
                                    if(debug) cat("\t\tDestination consensus did not remap.\n")
                                }

                            }
                        }
                    }
                }
            }             
        }
    } ### End of if(nchar(conseq)>SC.MIN.CONSENSUS.LENGTH)
    list(nmatch=nmatch,minDistToCluster=minDistToCluster,gapseqRecipient=gapseqRecipient,nremap=nremap,remapDistToCluster=remapDistToCluster,gapseqDonor=gapseqDonor,
          nmatchRC=nmatchRC,minDistToClusterRC=minDistToClusterRC,gapseqRecipientRC=gapseqRecipientRC,nremapRC=nremapRC,remapDistToClusterRC=remapDistToClusterRC,gapseqDonorRC=gapseqDonorRC,
          nvcfentries=nvcfentries,bndlist=bndlist)
}

posIncreasing<-function(chr1,pos1,chr2,pos2){
    if(paste(chr1)==paste(chr2)) return(pos2>pos1)
    else {
        return(paste(chr2)>paste(chr1))
    }
}


getCoverageFoldchange<-function(chr,st,ed, scalefac, bamfile1, bamfile0=NULL,MIN.BRACKETWIN=10000, MAX.EVENT.HALF=50000){
        what="pos"
        donorcov=0
        if(ed-st>2*MAX.EVENT.HALF) {
          event_st=c(st,ed-MAX.EVENT.HALF); event_ed=c(st+MAX.EVENT.HALF,ed)
        } else {
          event_st=c(st,st+ceiling((ed-st)/2)); event_ed=c(st+floor((ed-st)/2),ed)
        }
        event_load_size=min(2*MAX.EVENT.HALF,ed-st)
        BRACKETWIN=max(MIN.BRACKETWIN,event_load_size)
        for(fi in 1:length(bamfile1)){
            bam1=allFunction(GRanges(chr,IRanges(start=event_st,end=event_ed)),bamfile1[fi],what)$pos
            donorcov = donorcov+length(bam1)
        }
        if(is.null(bamfile0)){
            bracketrange = IRanges(start=c(st-BRACKETWIN,ed),end=c(st,ed+BRACKETWIN))
            scalefac=event_load_size/sum(width(bracketrange))
            contrastcov=1
            for(fi in 1:length(bamfile1)){
                bam0=allFunction(GRanges(chr,bracketrange),bamfile1[fi],what)$pos
                contrastcov = contrastcov+length(bam0)
            }
        } else {
            contrastcov=1
            for(fi in 1:length(bamfile0)){
                bam0=allFunction(GRanges(chr,IRanges(start=event_st,end=event_ed)),bamfile0[fi],what)$pos
                contrastcov=contrastcov+length(bam0)
            }
        }
        bamSuccess=(donorcov>0) && (contrastcov>1)
        fc = (donorcov/contrastcov)/scalefac
        if(length(fc)==0) { if(debug) cat("error fold change",paste(chr,st,ed,scalefac,donorcov,contrastcov,sep=" "),"\n"); quit() }
        list(fc=fc,contrastcov=contrastcov,bamSuccess=bamSuccess)
}

sclipPlot<-function(chr,st,ed,sctabL,sctabR,ylim=NA,...){
    # This function can be called after sclip_scan
    # Typically:
    #   sctabL=summarizeClusters(scL1)
    #   sctabR=summarizeClusters(scR1)
    
    sel = which(sctabL$Chrom==chr & sctabL$Position>st & sctabL$Position<ed)
    selR = which(sctabR$Chrom==chr & sctabR$Position>st & sctabR$Position<ed)
    
    if(is.na(ylim)){
        ymax=max(c(sctabL$nReads[sel],sctabR$nReads[selR]))*1.2
        ylim=c(0,ymax)
    } 
    
    plot(sctabL$Position[sel],sctabL$nReads[sel],xlim=c(st,ed), ylim=ylim,xlab="Position",ylab="# Reads Clipped at Position",type="p",col="skyblue", lwd=2, cex.lab=1.5, cex.main=2)
    segments(sctabL$Position[sel],ylim[1],sctabL$Position[sel],sctabL$nReads[sel],lwd=2,col="skyblue")
    points(sctabR$Position[selR],sctabR$nReads[selR],col="red")
    segments(sctabR$Position[selR],ylim[1],sctabR$Position[selR],sctabR$nReads[selR],lwd=1,col="red")
    grid()    
    legend(x="topright", col=c("skyblue","red"),lty=c(1,1),pch=c(1,1),legend=c("Left clipped","Right clipped"), cex=1.5)
}
# try:
# chr="1"
# st=20000000; ed=90000000
# par(mfrow=c(2,1))
# sclipPlot(chr,st,ed,sctabL,sctabR)
# bndPlot(chr,st,ed,bndtab,show.legend=TRUE, legend.cex=1)
bndPlot<-function(chr,st,ed,bndtab,main="Breakends Found by Soft-clipping",show.legend=TRUE,legend.cex=1.5,...){
# This function can be called after sclip_call,sclip_events has been run.
# bndtab is saved in the .RData at end of sclip_events.

    selb1 = which(bndtab$chr1==chr & bndtab$pos1>st & bndtab$pos1<ed & bndtab$chr2==chr)
    selb1n = which(bndtab$chr1==chr & bndtab$pos1>st & bndtab$pos1<ed & bndtab$chr2!=chr)
    selb2 = which(bndtab$chr2==chr & bndtab$pos2>st & bndtab$pos2<ed & bndtab$chr1==chr)
    selb2n = which(bndtab$chr2==chr & bndtab$pos2>st & bndtab$pos2<ed & bndtab$chr1!=chr)
    selb = union(selb1,selb2)
    ch = c(-8592,-8594)# (point left, point right)
    #ch = c("<",">")# (point left, point right)
    par(yaxt="n")
    ymax=length(selb)+1
    plot(0,0,cex=0,xlim=c(st,ed),ylim=c(0,ymax),ylab="",xlab="Base Position",main=main)
    textstr<-function(id,xpos,ypos){
            str=""
        if(!is.na(bndtab$EventType[id])){ 
            str=bndtab$EventType[id]
            if(!is.na(bndtab$FoldChange[id])) str=paste(str,", FC=",prettyNum(bndtab$FoldChange[id],digits=2),sep="")
        } else {
            if(!is.na(bndtab$FoldChange[id])) str=paste("FC=",prettyNum(bndtab$FoldChange[id],digits=2),sep="")
        }
        text(xpos,ypos,str)

    }
    if(length(selb)>0){
        for(i in 1:length(selb)){
            height=i
            points(bndtab$pos1[selb[i]],height,pch=ch[bndtab$goRight1[selb[i]]+1],cex=1.5, col="blue")
            points(bndtab$pos2[selb[i]],height,pch=ch[bndtab$goRight2[selb[i]]+1],cex=1.5, col="red")
            segments(bndtab$pos1[selb[i]],height,bndtab$pos2[selb[i]],height,col="black",lty=2)
            xx=0.5*(bndtab$pos1[selb[i]]+bndtab$pos2[selb[i]])
            textstr(selb[i],xx,height+0.5)
        }
    }
    if(length(selb1n)>0){
        for(i in 1:length(selb1n)){
            points(bndtab$pos1[selb1n[i]],ymax,pch=ch[bndtab$goRight1[selb1n[i]]+1],cex=1.5)    
            textstr(selb1n[i],bndtab$pos1[selb1n[i]],ymax)
        }
    }
    if(length(selb2n)>0){
        for(i in 1:length(selb2n)){
            points(bndtab$pos2[selb2n[i]],ymax,pch=ch[bndtab$goRight1[selb2n[i]]+1],cex=1.5)    
            textstr(selb2n[i],bndtab$pos2[selb2n[i]],ymax)
        }
    }
    if(show.legend) legend(x="topright", pch=ch,col=c("blue","red"),legend=c("Left Bnd of intrachromosomal pair","Right Bnd of intrachromosomal pair"), cex=legend.cex)
}



sclipEventsToVCF<-function(bndtab,scL1,scR1,SOMATIC=TRUE){

    ### Generate VCF tab from saved sclip_events RDATA for output.
    CHROM=rep("",0)
    POS=rep(NA,0)
    ID=rep("",0)
    REF=rep(NA,0)
    ALT=rep("",0)
    INFO=rep("",0)
    finished=rep(FALSE,nrow(bndtab))
    if(SOMATIC){
        somaticstr = "SOMATIC;"
    }else {
        somaticstr = ""
    }

    sel= which(!is.na(bndtab$EventType))
    for(rowi in sel){
        if(!finished[rowi]){
            if(bndtab$goRight1[rowi]){
                nsclip=scL1[[bndtab$cid[rowi]]]$nreads
            } else {
                nsclip=scR1[[bndtab$cid[rowi]]]$nreads
            }
            if(bndtab$goRight2[rowi]){
                nsclip=nsclip+scL1[[bndtab$mcid[rowi]]]$nreads
            } else {
                nsclip=nsclip+scR1[[bndtab$mcid[rowi]]]$nreads
            }
    
            if(bndtab$EventType[rowi]=="DEL"){
                CHROM=c(CHROM,bndtab$chr1[rowi])
                POS=c(POS,bndtab$pos1[rowi])
                ID=c(ID,paste("DEL",bndtab$cid[rowi],bndtab$mcid[rowi],sep="."))
                ALT=c(ALT,"<DEL>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=DEL;END=",bndtab$pos2[rowi],";SVLEN=",bndtab$Distance[rowi],";BNDNUM=",rowi,";FOLDCHANGE=",bndtab$FoldChange[rowi],";NSCLIP=",nsclip,sep=""))
            }
            if(bndtab$EventType[rowi]=="TANDEM.DUP"){
                CHROM=c(CHROM,bndtab$chr1[rowi])
                POS=c(POS,bndtab$pos1[rowi])
                ID=c(ID,paste("TANDEM.DUP",bndtab$cid[rowi],bndtab$mcid[rowi],sep="."))
                ALT=c(ALT,"<DUP>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=DUP;END=",bndtab$pos2[rowi],";SVLEN=",bndtab$Distance[rowi],";BNDNUM=",rowi,";FOLDCHANGE=",bndtab$FoldChange[rowi],";NSCLIP=",nsclip,sep=""))
            }
            if(bndtab$EventType[rowi]=="INV"){
                CHROM=c(CHROM,bndtab$chr1[rowi])
                st1 = bndtab$pos1[rowi]; st2 = bndtab$pos1[as.numeric(paste(bndtab$PartnerBnds[rowi]))]
                ed1 = bndtab$pos2[rowi]; ed2 = bndtab$pos2[as.numeric(paste(bndtab$PartnerBnds[rowi]))]
                st = round((st1+st2)/2)
                ed=round((ed1+ed2)/2)
                POS=c(POS,st)
                ID=c(ID,paste("INV",bndtab$cid[rowi],bndtab$mcid[rowi],sep="."))
                ALT=c(ALT,"<INV>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=INV;END=",ed,";SVLEN=",ed-st,";BNDNUM=",rowi,",",bndtab$PartnerBnds[rowi],";CIPOS=",min(st1,st2)-st,",",max(st1,st2)-st,";CIEND=",min(ed1,ed2)-ed,",",max(ed1,ed2)-ed,";NSCLIP=",nsclip,sep=""))
                finished[as.numeric(paste(bndtab$PartnerBnds[rowi]))]=TRUE
            }
            if(bndtab$EventType[rowi]=="INS.DUP"){
                partner.row=as.numeric(paste(bndtab$PartnerBnds[rowi]))
                CHROM=c(CHROM,bndtab$chr1[rowi],bndtab$chr1[partner.row])
                POS=c(POS,bndtab$pos1[rowi],bndtawb$pos1[partner.row])
                this.id = paste(somaticstr,"INS.DUP",bndtab$cid[rowi],bndtab$mcid[rowi],sep=".")
                partner.id=paste(somaticstr,"INS.DUP",bndtab$cid[partner.row],bndtab$mcid[partner.row],sep=".")
                ID=c(ID,this.id,partner.id)
                ALT=c(ALT,"<INS.DUP>","<INS.DUP>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=INS.DUP;END=",bndtab$pos2[rowi],";SVLEN=",bndtab$Distance[rowi],";PARTNERID=",partner.id,";BNDNUM=",rowi,";FOLDCHANGE=",bndtab$FoldChange[rowi],";NSCLIP=",nsclip,sep=""))
                INFO=c(INFO,paste(somaticstr,"SVTYPE=INS.DUP;END=",bndtab$pos2[partner.row],";SVLEN=",bndtab$Distance[partner.row],";PARTNERID=",this.id,";BNDNUM=",partner.row,";FOLDCHANGE=",bndtab$FoldChange[partner.row],";NSCLIP=",nsclip,sep=""))
                finished[as.numeric(paste(bndtab$PartnerBnds[rowi]))]=TRUE
            }
            if(bndtab$EventType[rowi]=="TRP"){
                partner.row=as.numeric(paste(bndtab$PartnerBnds[rowi]))
                CHROM=c(CHROM,bndtab$chr1[rowi],bndtab$chr1[partner.row])
                POS=c(POS,bndtab$pos1[rowi],bndtawb$pos1[partner.row])
                this.id = paste("TRP",bndtab$cid[rowi],bndtab$mcid[rowi],sep=".")
                partner.id=paste("TRP",bndtab$cid[partner.row],bndtab$mcid[partner.row],sep=".")
                ID=c(ID,this.id,partner.id)
                ALT=c(ALT,"<TRP>","<TRP>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=TRP;END=",bndtab$pos2[rowi],";SVLEN=",bndtab$Distance[rowi],";PARTNERID=",partner.id,";BNDNUM=",rowi,";FOLDCHANGE=",bndtab$FoldChange[rowi],";NSCLIP=",nsclip,sep=""))
                INFO=c(INFO,paste(somaticstr,"SVTYPE=TRP;END=",bndtab$pos2[partner.row],";SVLEN=",bndtab$Distance[partner.row],";PARTNERID=",this.id,";BNDNUM=",partner.row,";FOLDCHANGE=",bndtab$FoldChange[partner.row],";NSCLIP=",nsclip,sep=""))
                finished[as.numeric(paste(bndtab$PartnerBnds[rowi]))]=TRUE
            }
    
            if(bndtab$EventType[rowi]=="INTERCHROMOSOMAL"){
                CHROM=c(CHROM,bndtab$chr1[rowi])
                POS=c(POS,bndtab$pos1[rowi])
                ID=c(ID,paste("INTERCHROMOSOMAL",bndtab$cid[rowi],bndtab$mcid[rowi],sep="."))
                ALT=c(ALT,"<INTERCHROMOSOMAL>")
                INFO=c(INFO,paste(somaticstr,"SVTYPE=INTERCHROMOSOMAL;MATEPOS=",bndtab$pos2[rowi],";MATECHR=",bndtab$chr2[rowi],";BNDNUM=",rowi,";NSCLIP=",nsclip,sep=""))
            }
        }
    }
    
#    sel = which(is.na(bndtab$EventType))
#    
#    
#    for(rowi in sel){
#     if(bndtab$goRight1[rowi]){
#                nsclip=scL1[[bndtab$cid[rowi]]]$nreads
#            } else {
#                nsclip=scR1[[bndtab$cid[rowi]]]$nreads
#            }
#            if(bndtab$goRight2[rowi]){
#                nsclip=nsclip+scL1[[bndtab$mcid[rowi]]]$nreads
#            } else {
#                nsclip=nsclip+scR1[[bndtab$mcid[rowi]]]$nreads
#            }
#            if(bndtab$chr1[rowi]=="1"){ 
#                alt="INV"
#            } else {
#                alt="DEL"
#            }
#                CHROM=c(CHROM,bndtab$chr1[rowi])
#                POS=c(POS,bndtab$pos1[rowi])
#                ID=c(ID,paste(alt,bndtab$cid[rowi],bndtab$mcid[rowi],sep="."))
#                ALT=c(ALT,paste("<",alt,">",sep=""))
#                INFO=c(INFO,paste(somaticstr,"SVTYPE=",alt,";END=",bndtab$pos2[rowi],";SVLEN=",bndtab$Distance[rowi],";BNDNUM=",rowi,";FOLDCHANGE=",bndtab$FoldChange[rowi],";NSCLIP=",nsclip,sep=""))
#    }
    
    REF=rep(".",length(ALT))
    QUAL=rep("100",length(ALT))
    FILTER=rep("PASS",length(ALT))
    
    tab=data.frame(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO)
    tab
}

# Step 1: Hotspot filtering. Get rid of entries in sctabs that are clustered.
findHotspots<-function(x,DIST.TO.NEXT.IN.HOTSPOT=1000,HOTSPOT.CLUSTER.SIZE=5){
    if(length(x)<=2) return(rep(FALSE,0))
    toNext=rep(NA,length(x)-1)
    for(i in 1:(length(x)-1)){
        toNext[i] = x[i+1]-x[i]
    }
    runst=0; runend=-1
    inHotspot=rep(FALSE,length(x))
    for(i in 1:(length(x)-1)){
    if(toNext[i]<DIST.TO.NEXT.IN.HOTSPOT && toNext[i]>=0){
        #cat(i,"toNext=",toNext[i],": part of run, st=",runst,", end=",runend,"\n")
        if(runst==0) runst=i;
        runend=i
    } else {
        #cat(i,"toNext=",toNext[i],": not in run, st=",runst,", end=",runend,"\n")
        if(runend-runst>HOTSPOT.CLUSTER.SIZE) inHotspot[runst:runend]=TRUE
        runst=0; runend=-1
    }
  }
  inHotspot  
}

