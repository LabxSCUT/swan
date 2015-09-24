#!/usr/bin/env Rscript
####Externalities####
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(swan))
#source("ThresholdFunctions.R")
cat("-Info: invoking command:",commandArgs(),"\n")

####CommandLines####
option_list2 <- list(
  make_option(c("-b", "--binSize"), type="integer", default=10000,
              help="bin size for computing coverage, [default %default]"),
  make_option(c("-k", "--runmedk"), type="integer", default=3,
              help="smoothing parameter, [default %default]"),
  make_option(c("-a", "--alpha"), default=0.1,
              help="pvalue for computing threshold, [default %default]"),
  make_option(c("-i", "--lambdaInflation"), default=1,
              help="inflation factor for lambda, [default %default]"),
  make_option(c("-m", "--minCov"), default="learn",
              help="lower bound on coverage, [default %default]"),
  make_option(c("-o", "--out"), default="thresholds",
              help="prefix of output files, [default %default]"),
  make_option(c("-d", "--plot"), type="integer", default=1,
              help="type of diagnostic plotting (0=none,1=coverage), [default %default]")
)
parser <- OptionParser(usage = "%prog [options] bamfile parfile", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
#print(args)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
#quit if not valid cmd args
if(length(cmd$args)!=2){ 
  print_help(parser)
  quit()
}

#print(cmd$options)
bam_file=cmd$args[1]
par_file=cmd$args[2]
outfile = cmd$options$out
cat("Computing thresholds with bam_file=",bam_file,"\n")
cat("\tpar_file: ",par_file,"\n")
cat("\toutfile: ",outfile,"\n")
binSize=cmd$options$binSize
alpha=cmd$options$alpha
minCov=cmd$options$minCov
runmedk=cmd$options$runmedk
plotType=cmd$options$plot
lambdaInflation=cmd$options$lambdaInflation
cat("\tbinSize:",binSize,"\n")
cat("\talpha:",alpha,"\n")

#####################################
# Read in scan parameters.
#####################################
pars = read.table(par_file, header=TRUE, sep="\t")
w=pars$w
stepsize = pars$stepsize
R=pars$rl   
r=pars$r
mu=pars$isize
sigmaL=pars$isize_sdL
sigmaR=pars$isize_sdR
pL=pars$p_left
pR=pars$p_right
lambda=pars$lambda
goodlCi = pars$lCi
nbases=pars$r_end-pars$r_start            
D = pars$delta          
R1=R-pars$prop_clip     
R2=R-pars$prop_clip     
kappa=sqrt(lambda)
N = lambda*nbases
pbaseL= pL*N/nbases
pbaseR= pR*N/nbases
n=nbases/stepsize

### If multiple libraries, assume all libraries scanned using the same r,w,stepsize,n,R1,R2,D,R,nbases
nbases=nbases[1]
r=r[1]
w=w[1]
stepsize=stepsize[1]
n=n[1]
###these might change
#R1=R1[1]
#R2=R2[1]
#R=R[1]

if(minCov=="learn"){
    minCov=sum(pars$coverage)*0.8
}

nlib = nrow(pars)  ## Number of libraries. 
if(nlib>1){
    # There are multiple bamfiles.  In this case, bam_file stores the names of the files, separated by colon.
    # All of the bam files need to have the same target names, i.e. chromosome 1 needs to be
    # called "chr1" in all bam files, not "1" in some and "chr1" in others.  
    # Basically the targets field in the header of the bam files need to be the same.
    bam_file =strsplit(bam_file,split=",")[[1]] #change to , as : is used for separating samples
    if(length(bam_file) != nlib) stop("Number of bam files",length(bam_file),"need to equal number of rows in par file",nlib,"\n.")
    cat("\nThere are multiple libraries:\n")
    for(j in 1:length(bam_file)) cat(bam_file[j],"\n")
}

#####################################
# Get the coverages.
#####################################
what=c("pos")
header=scanBamHeader(bam_file[1])
chrnames=names(header[[1]]$targets)  # find out what the chromosomes are called.
lockappa = vector("list",length(header[[1]]$targets))
locbinStart = vector("list",length(header[[1]]$targets))
minkappa= rep(Inf,nlib)
maxkappa= rep(-Inf,nlib)
nobs = rep(0,length(header[[1]]$targets))
MINOBS=100

cat("\nGetting coverages from bam file(s)...\n")
for(chri in 1:length(header[[1]]$targets)){
    seq = chrnames[chri]
    #cat("br0\n")
    which = RangesList("quack"=IRanges(1,header[[1]]$targets[chri]))
    names(which) = seq
    param <- ScanBamParam(what = what,
                        which = which)
    pos = vector("list",nlib)
    maxpos=0; minpos = header[[1]]$targets[[chri]]
    for(libi in 1:nlib){
        x<- scanBam(bam_file[libi], param = param)[[1]]
        pos[[libi]] = x[["pos"]]
        cat(sum(!is.na(pos[[libi]]))," reads mapped to ",seq," in library ",libi,".\n",sep="")
        thismin=min(pos[[libi]],Inf,na.rm=TRUE)
        if(thismin<minpos) minpos = thismin
        thismax=max(pos[[libi]],-Inf,na.rm=TRUE)
        if(thismax>maxpos) maxpos = thismax
        nobs[chri] = nobs[chri]+length(pos[[libi]])
    }
    if(nobs[chri]>MINOBS){
        bins = IRanges(start=seq(minpos,maxpos,binSize),width=binSize)
        covrm = matrix(nrow=length(bins),ncol=nlib)
        for(libi in 1:nlib){ 
            posRanges=IRanges(start=pos[[libi]][!is.na(pos[[libi]])],width=R[libi])
            cov = countOverlaps(bins,posRanges)
            covrm[,libi] = runmed(cov,k=runmedk)
            if(plotType==1){
                pdf_file=paste(outfile,".",seq,".",libi,".cov.pdf",sep="")
                pdf(pdf_file)
                plot(start(bins),cov, main=paste("Coverage in ",binSize," bins on target ",seq,sep=""),xlab="Bin Start",ylab="Read Count")
                lines(start(bins),covrm[,libi], col="red")
                dev.off()
            }
        }
        locbinStart[[chri]] = start(bins)
        lockappa[[chri]] = sqrt(covrm/binSize)
        minkappa = pmin(minkappa,apply(lockappa[[chri]],2,min))
        maxkappa = pmax(maxkappa,apply(lockappa[[chri]],2,max))
    }
}

if(sum(nobs>MINOBS)==0) stop("Not enough mapped reads to compute thresholds.")
MAX.KAPPAS=5000
step = round(max(0.01,(prod(maxkappa-minkappa)/MAX.KAPPAS)^(1/nlib)),digits=2)
cat("Computing kappa-thresholds table with step=",step,"\n")
kappavec=vector("list",nlib)
for(libi in 1:nlib) kappavec[[libi]] = seq(minkappa[libi],maxkappa[libi]+step,step)
kappas = expand.grid(kappavec) 
tvec_lcd = rep(NA,nrow(kappas))
tvec_lci= rep(NA,nrow(kappas))
tvec_ldL = rep(NA,nrow(kappas))
tvec_ldR = rep(NA,nrow(kappas))

kappalen = unlist(lapply(kappavec,length))
if(nlib>1){
    kappafac = cumprod(kappalen)
    kappafac = c(1,kappafac[1:(nlib-1)])
} else {
    kappafac = 1
}
############

sigma=max(sigmaL,sigmaR, na.rm=TRUE)

cat("\nComputing thresholds look-up table...\n")
for(k in 1:nrow(kappas)){
    cat("\n--- Computing threshold for kappa set ",k," out of ",nrow(kappas)," ---\n")
    if(min(kappas[k,]) != 0){
        kappa=unlist(kappas[k,])
        lambda = lambdaInflation*kappa^2
        N = lambda*nbases
        #pbaseL= pL*N/nbases;pbaseR= pR*N/nbases 
        ptm=proc.time()
        tvec_lcd[k] = getThresholdStraddlingReadsDeletionScan(alpha=alpha,kappa=kappa,w=w,r=r,mu=mu,sigma=sigmaR,n=n, R1=median(R1), R2=median(R2), UPPER=NA,verbose=FALSE,status=pars[["lCd"]])
        tvec_lci[k] = getThresholdStraddlingReadsInsertionScan(alpha=alpha,kappa=kappa,w=w,r=r,mu=mu,sigma=sigmaL,n=n, R1=median(R1), R2=median(R2), UPPER=NA,verbose=FALSE,status=pars[["lCi"]])
        tvec_ldL[k] = getThresholdHangingReadsScan(alpha=alpha,kappa=kappa,w=w,r=r,p=pL,pbase=pbaseL,mu=mu,sigma=sigmaR,n=n,D=D,R=median(R),UPPER=NA,verbose=FALSE,status=pars[["lDl"]])
        tvec_ldR[k] = getThresholdHangingReadsScan(alpha=alpha,kappa=kappa,w=w,r=r,p=pR,pbase=pbaseR,mu=mu,sigma=sigmaR,n=n,D=D,R=median(R),UPPER=NA,verbose=FALSE,status=pars[["lDr"]])
        elapsed=proc.time()-ptm
    }
}

cat("\nMatching thresholds to each bin based on kappa value ...\n")
threshtab = data.frame(matrix(nrow=0,ncol=7))
for(chri in 1:length(header[[1]]$targets)){
    if(nobs[chri]>MINOBS){
        thiskappas= lockappa[[chri]]
        minkappamat = matrix(data=minkappa, nrow=nrow(thiskappas),ncol=length(minkappa),byrow=TRUE)
        inds = ceiling((thiskappas-minkappamat)/step)
        kappafacmat = matrix(data=kappafac,nrow=nrow(inds),ncol=ncol(inds),byrow=TRUE)
        ix=apply(inds*kappafacmat,1,sum) 
        ix=ix+1
        t_lcd = tvec_lcd[ix]
        t_lci = tvec_lci[ix]
        t_ldL = tvec_ldL[ix]
        t_ldR = tvec_ldR[ix]
        threshtab = rbind(threshtab,cbind(chrnames[chri],locbinStart[[chri]],locbinStart[[chri]]+binSize-1, t_lcd,t_lci,t_ldL,t_ldR))
    }
}

cat("\nWriting thresholds to file ...\n")
names(threshtab) = c("chrom","chromStart","chromEnd","lCd","lCi","lDl","lDr")
write.table(threshtab, file=paste(outfile,".txt",sep=""), quote=FALSE, sep="\t",row.names=FALSE)
cat("Done.\n")
cat("warnings if any\n")
warnings()
