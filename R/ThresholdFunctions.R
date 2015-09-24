#This file is from Nancy and maintained by Charlie
computeThresholds<-function(pars,alpha){
    # Wrapper function for computing thresholds based on SWAN's pars file.
    
    w=pars$w
    pleft=pars$p_left
    pright=pars$p_right
    p = pleft
    r=pars$r
    mu=pars$isize
    sigma=pars$isize_sd
    lambda=pars$lambda
    stepsize = pars$stepsize
    nbases=pars$end-pars$start
    R=pars$rl
    D = pars$delta
    R1=R-pars$hang_clip
    R2=R-pars$hang_clip
    kappa=sqrt(lambda)
    N = lambda*nbases
    pbase= p*N/nbases
    n=nbases/stepsize
    
    psidot0=psidotbetaStraddlingReadsDeletionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=100,R2=0,doplots=TRUE)
    nulvar=sigmabetaStraddlingReadsDeletionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=100,R2=0,varZt=TRUE,doplots=TRUE)
    m0LCd = psidot0
    s0LCd = sqrt(nulvar)
    bLCd = rep(0,length(alpha))  # thresholds corresponding to alpha.
    for(i in 1:length(alpha)) 
    bLCd[i]=getThresholdStraddlingReadsDeletionScan(alpha=alpha[i],kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n,R1=100,R2=0,verbose=FALSE,status=NULL,useEmpiricalDist=TRUE)
    
    m0HR=psidotbetaHangingReadsScan(beta=0,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,D=D,R=R)
    nullvar=sigmabetaHangingReadsScan(beta=0,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,D=D,R=R, 
    varZt=TRUE)
    s0HR = sqrt(nullvar)
    bHR =rep(0,length(alpha))
    for(i in 1:length(alpha)) 
    bHR[i]=getThresholdHangingReadsScan(alpha=alpha[i],kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,n=n,D=D,R=R,UPPER=NA, verbose=FALSE, status=NULL)
    
    m0LCi=psidotbetaStraddlingReadsInsertionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2)
    nullvar=sigmabetaStraddlingReadsInsertionScan(beta=0,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,R1=R1, R2=R2, varZt=TRUE)
    s0LCi = sqrt(nullvar)
    cat("When there is no signal, the LCi process should have mean ",m0LCi," and standard deviation ",s0LCi,".\n",sep="")
    bLCi = rep(0,length(alpha))
    for(i in 1:length(alpha)) bLCi[i]=getThresholdStraddlingReadsInsertionScan(alpha=0.1,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=NA,verbose=FALSE,status=NULL)

    list(bLCd=bLCd, bHR=bHR, bLCi=bLCi)
}


computelambda0<-function(b,pvalfun,n,n0=100000,...){
    # This function computes the expected number of false positives (lambda0)
    # for a deletion scan using the given parameters.
    # ------------------------------------------------------------------------
    # b: threshold
    # pvalfun: p-value function for computing the pvalue of b.
    #       Assumes that pvalfun takes in at least two parameters, b and n.
    # ... : additional parameters passed on to pvalfun.
    # n: total number of hypotheses.
    # n0: blocksize.
    # ------------------------------------------------------------------------    
    nblocks = floor(n/n0)
    lastn0 = n-n0*nblocks
    if(nblocks>=1){
        pblock = pvalfun(b=b,n=n0,...)
        cat("The pvalue for each block is ",pblock,"\n",sep="")
        if(pblock>0 && pblock<1){
            lamblock = -log(1-pblock)
        } else {
            lamblock = pblock
        }
    } else {
        lamblock=0
    }
    plast = pvalfun(b=b,n=lastn0,...)
    if(plast>0 && plast<1){
        lamlast = -log(1-plast)
    } else {
        lamlast = plast
    }
    lambda0 = lamblock*nblocks+lamlast
    lambda0
}

################################################################################################################
# The following pair of functions deal with scan of hanging + strand (or - strand) read pairs (i.e. lDl and lDr).
#   One gets the threshold for a given p-value (i.e. FWER)
#   The other give sthe p-value for a given threshold.
#
# The parameters are:
# 
#   b is the threshold (for pval...)
#   alpha is the pvalue (for getThreshold...)
#   n is the number of hypotheses tested (i.e. the total number of bases divided by step size).
#   r is the proportion of carriers assumed in the scan.
#   w is the window size used in the scan (currently, this parameter is not used, and so does not affect p-value).
#   D, R: Improper pairs with + strand mapping within [t-D,t-R] are counted (or with minus strand mapping within [t+R, t+D]).
#   p is the proportion of pairs meeting above criterion among all pairs with mapped + strand reads. It can be approximated by |S^-|/N, where N is the number of mapped + strand reads.
#   pbase is the base-wise rate of an improper minus-strand (or plus strand) pair.  It can be approximated by |S^-|/T, where T is the size of the reference.
#   kappa^2*T is the total number of mapped pairs where both reads are mapped, T is the effective length of the reference.
#   (mu, sigma) are mean and std dev of the insert length distribution (assumed Gaussian).
#   UPPER: technical parameter that users should not fiddle with.  It is the upper limit used in search for betab.
#   verbose: Do we want to print out status messages?
#################################################################################################################

getThresholdHangingReadsScan<-function(alpha,kappa,w,r,p,pbase,mu,sigma,n,D,R, F=NULL, UPPER=NA,verbose=TRUE,status=NULL){
    if(verbose) cat("\nComputing threshold for Hanging Reads Scan for FWER=",alpha,"...\n")
      
    nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
    if(length(mu)!= nlibs || length(sigma) != nlibs || length(p)!=nlibs || length(pbase)!=nlibs || length(D)!=nlibs) stop("The sizes of kappa, mu, sigma, p, pbase, and D need to match.\n")    
    if(!is.null(F) && length(kappa)>1) stop("Can not use empirical distribution with multiple libraries.\n") 

    pvalminusalpha = function(b){ pvalHangingReadsScan(b=b,kappa=kappa,w=w,r=r,p=p,pbase=pbase,mu=mu,sigma=sigma,n=n,D=D,R=R,F=F,UPPER=UPPER,verbose=verbose,status=status)-alpha }
 
    m0=0; for(libi in 1:nlibs){
      if(status[libi])
        m0=m0+psidotbetaHangingReadsScan(beta=0,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F)}
    s0=0; for(libi in 1:nlibs){
      if(status[libi])
        s0=s0+sigmabetaHangingReadsScan(beta=0,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F, varZt=TRUE)}
    s0=sqrt(s0)
    
    if(verbose) cat("Hanging reads score has null mean = ",m0,", standard deviation = ",s0,".\n")
    if(m0==0&s0==0) return(0)
    
    K=6
    prevroot = m0
    lower = m0+0.1*s0
    while(TRUE){
        upper = m0+K*s0
        if(verbose) cat("\nTrying to find threshold with search upper limit ", K," sd from mean: [",lower,", ",upper,"]\n",sep="")
        res=try(uniroot(pvalminusalpha, lower=m0+0.1*s0, upper=m0+K*s0), silent=TRUE)
        if(class(res) == "try-error"){
            if(verbose) cat("Upper limit too low, multiply K by 1.5.\n")
            K = K*1.5
        } else {
            if(verbose) cat("\n-----> uniroot converged on b=",res$root,".\n")
            break
        }
    }
    if(verbose) cat("Threshold value: ",res$root,"\n")
    res$root
}


pvalHangingReadsScan<-function(b,kappa,w,r,p,pbase,mu,sigma,n,D,R,F=NULL,UPPER=NA,verbose=TRUE,status=NULL){
    if(verbose) cat("\nComputing P-VALUE for Hanging Reads Scan with b=",b,"...\n")
   
   nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
   if(length(mu)!= nlibs || length(sigma) != nlibs || length(p)!=nlibs || length(pbase)!=nlibs || length(D)!=nlibs) stop("The sizes of kappa, mu, sigma, p, and pbase need to match.\n")    
   if(!is.null(F) && length(kappa)>1) stop("Can not use empirical distribution with multiple libraries.\n") 

   
    ### Get betab = root of psidotminusb.
    psidotminusb = function(beta){ 
        psidot=0
        for(libi in 1:nlibs) 
          if(status[libi])
            psidot = psidot+psidotbetaHangingReadsScan(beta=beta,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F)
        return(psidot-b)
    }
    if(is.na(UPPER)){
        UPPER = 2
    }
    if(psidotminusb(0)>=0){
        # Must start with psidotminusb(0)<0.
        stop("Threshold value must be larger than null mean.\n")
    }
    prevroot=0
    while(TRUE){
        if(verbose) cat("\tTrying to compute betab with search upper limit ", UPPER,"...\n",sep="")
        res=try(uniroot(psidotminusb, lower=0, upper=UPPER), silent=TRUE)
        if(class(res) == "try-error"){
            if(psidotminusb(UPPER)<0){
                # Must have returned error because psidotminusb(UPPER) is also negative.
                # uniroot must have psidotminusb(0)<0 (enforced above) and psidotminusb(UPPER)>0
                UPPER = UPPER*2
                if(verbose) cat("\tUpper limit too low and psidotminusb<0, increase upper to ", UPPER,".\n",sep="")
            } else{
                # Must have returned error because psidotminusb(UPPER) is infinity.
                if(verbose) cat("Upper limit too high, divide by 2.\n")
                UPPER = UPPER/1.5            
            }
        } else {
            break
        }
    }

    betab=res$root
    if(verbose) cat("Final value of beta_b: ",betab,".\n",sep="")
    
    sigmab=0; for(libi in 1:nlibs){ 
      if(status[libi])
        sigmab = sigmab+sigmabetaHangingReadsScan(beta=betab,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F)}
    deltab=0; for(libi in 1:nlibs){ 
      if(status[libi])
        deltab = deltab+deltabetaHangingReadsScan(beta=betab,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F)}
    psib = 0; for(libi in 1:nlibs){ 
      if(status[libi])
        psib = psib+ psibetaHangingReadsScan(beta=betab,kappa=kappa[libi],w=w,r=r,p=p[libi],pbase=pbase[libi],mu=mu[libi],sigma=sigma[libi],D=D[libi],R=R,F=F)}

    pval = n* exp(-(betab*b-psib))*(1/(sqrt(2*pi*sigmab)))*deltab
    if(verbose) cat("P-value=",pval,", betab=",betab,", sigmab=",sigmab,", deltab=",deltab, ", psib=",psib,"\n")
    pval
}

#######################################################################################################
# The following pair of functions deal with scan of straddling pairs by their insert length.
# They are for deletions only.
#   One gets the threshold for a given p-value (i.e. FWER)
#   The other gives the p-value for a given threshold.
#
# The parameters are:
# 
#   b is the threshold (for pval...)
#   alpha is the FWER threshold (for getThreshold...).
#   n is the number of hypotheses tested.
#   r is the proportion of carriers.
#   w is the size of deletion.
#   R1, R2: Pairs with + strand mapping within [0, t-R1] and - strand mapping within [t+w-R2,T] are counted.
#   kappa^2*T is the total number of mapped pairs where both reads are mapped, T is the effective length of the reference.
#   (mu, sigma) are mean and std dev of the insert length distribution (assumed Gaussian).
#   UPPER: technical parameter that users should not fiddle with.  It is the upper limit used in search for betab.
#   verbose: Do we want to print out status messages?
#################################################################################################################




getThresholdStraddlingReadsDeletionScan<-function(alpha,kappa,w,r,f=NULL,mu=NULL,sigma=NULL,n,R1,R2,UPPER=NA,verbose=TRUE,status=NULL,useEmpiricalDist=FALSE){
    
    #print(list(alpha,kappa,w,r,f,mu,sigma,n,R1,R2,UPPER,verbose,status,useEmpiricalDist))
    nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
    #print(status)
    if(length(mu)!= nlibs || length(sigma) != nlibs) stop("The sizes of kappa, mu, and sigma need to match.\n")    
    if((!is.null(f) || useEmpiricalDist) && length(kappa)>1) stop("Can not use empirical distribution with multiple libraries.\n") 
    pvalminusalpha = function(b){ pvalStraddlingReadsDeletionScan(b,kappa=kappa,w=w,r=r,f=f,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=UPPER,verbose=FALSE,status=status,useEmpiricalDist=useEmpiricalDist)-alpha}
    
    if(verbose) cat("\nComputing threshold for Straddling Reads Scan for FWER=",alpha,"...\n")
    
    m0=0; s0=0
    for(libi in 1:nlibs){
      if(status[libi]){
        psidot0i=psidotbetaStraddlingReadsDeletionScan(beta=0,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2,
                                            doplots=TRUE)
        nulvari=sigmabetaStraddlingReadsDeletionScan(beta=0,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2,
                                            varZt=TRUE,doplots=TRUE)
        if(length(nulvari)>1){
            if(useEmpiricalDist){
                s0 = s0+nulvari$analvarf
                m0=m0+psidot0i$analmeanf
            } else {
                s0 = s0+nulvari$analvardnorm
                m0 = m0+psidot0i$analmeandnorm
            }
        } else {
            m0 = m0+psidot0i
            s0 = s0+nulvari
        }
        if(verbose) cat("Library ",libi,": kappa=", kappa[libi], ", mu=",mu[libi],", sigma=",sigma[libi],", Null mean of score = ",psidot0i,", standard deviation of score= ",sqrt(nulvari),".\n")
      }
    }
    s0 = sqrt(s0)
    if(verbose) cat("Combined libraries: Null mean of score = ",m0,", standard deviation of score= ",s0,".\n")
    if(m0==0&s0==0) return(0)
    
    K=6
    prevroot = m0
    lower = m0+0.1*s0
    while(TRUE){
        upper = m0+K*s0
        if(upper>1e10) stop("Can not converge on a threshold value.\n")
        if(verbose) cat("\nTrying to find threshold with search upper limit ", K," sd from mean: [",lower,", ",upper,"]\n",sep="")
        res=try(uniroot(pvalminusalpha, lower=lower, upper=upper), silent=TRUE)
        if(class(res) == "try-error"){
            if(verbose) cat("Upper limit too low, multiply K by 1.5.\n")
            K = K*1.5
        } else {
            if(verbose) cat("\n-----> Threshold value coverged to b=",res$root,".\n")
            break
        }
    }
    if(verbose) cat("Threshold value: ",res$root,"\n")
    res$root
}




pvalStraddlingReadsDeletionScan<-function(b,kappa,w,r,f=NULL,mu,sigma,n, R1, R2, UPPER=NA,verbose=TRUE,status=NULL,useEmpiricalDist=FALSE){

    if(verbose) cat("\nComputing P-VALUE for Straddling Pairs Scan with b=",b,"...\n")
    nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
    if(length(mu)!= nlibs || length(sigma) != nlibs) stop("The sizes of kappa, mu, and sigma need to match.\n")
    if((!is.null(f) || useEmpiricalDist) && length(kappa)>1) stop("Can not use empirical distribution with multiple libraries.\n")
   
    psidotminusb = function(beta){ 
        psidot=0
        for(libi in 1:nlibs){
          if(status[libi]){
            psidoti=psidotbetaStraddlingReadsDeletionScan(beta=beta,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2)
            if(length(psidoti)>1){
                if(useEmpiricalDist){
                    psidot = psidot+ psidoti$analmeanf
                } else {
                    psidot = psidot + psidoti$analmeandnorm 
                }
            } else { 
                psidot=psidot + psidoti 
            }
          }
        }
        return(psidot-b)
    }
        
    if(is.na(UPPER)){
        UPPER = 2
    }
    if(psidotminusb(0)>=0){
        # Must start with psidotminusb(0)<0.
        stop("Threshold value must be larger than null mean.\n")
    }
    prevroot=0
    while(TRUE){
        if(verbose) cat("\tTrying to compute betab with search upper limit ", UPPER,"...\n",sep="")
        res=try(uniroot(psidotminusb, lower=0, upper=UPPER), silent=TRUE)
        if(class(res) == "try-error"){
            if(psidotminusb(UPPER)<0){
                # Must have returned error because psidotminusb(UPPER) is also negative.
                # uniroot must have psidotminusb(0)<0 (enforced above) and psidotminusb(UPPER)>0
                UPPER = UPPER*2
                if(verbose) cat("\tUpper limit too low and psidotminusb<0, increase upper to ", UPPER,".\n",sep="")
            } else{
                # Must have returned error because psidotminusb(UPPER) is infinity.
                if(verbose) cat("Upper limit too high, divide by 2.\n")
                UPPER = UPPER/1.5            
            }
        } else {
            break
        }
    }

    betab=res$root
    if(verbose) cat("Final value of beta_b: ",betab,".\n",sep="")

    sigmab = 0
    for(libi in 1:nlibs){
      if(status[libi]){
        sigmabi = sigmabetaStraddlingReadsDeletionScan(beta=betab,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2,doplots=FALSE)
        if(length(sigmabi)>1){ 
            sigmab = ifelse(useEmpiricalDist, sigmab+sigmabi$analvarf, sigmab+sigmabi$analvardnorm)
        } else {
            sigmab = sigmab + sigmabi
        }
      }
    }
    deltab = 0 
    for(libi in 1:nlibs){
      if(status[libi]){
        deltabi = deltabetaStraddlingReadsDeletionScan(beta=betab,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2,doplots=FALSE)
        if(length(deltabi)>1){
            deltab = ifelse(useEmpiricalDist, deltab+deltabi$deltaf, deltab+deltabi$deltadnorm)
        } else {
            deltab = deltab + deltabi
        }
      }
    }
    psib = 0
    for(libi in 1:nlibs){
      if(status[libi]){
        psibi = psibetaStraddlingReadsDeletionScan(beta=betab,kappa=kappa[libi],w=w,r=r,f=f,mu=mu[libi],sigma=sigma[libi],R1=R1,R2=R2,doplots=FALSE)
        if(length(psibi)>1){
            psib = ifelse(useEmpiricalDist, psib+psibi$psibetaf, psib+psibi$psibetadnorm)    
        } else {
            psib = psib + psibi
        }
      }
    }    
        
    pval = n* exp(-(betab*b-psib))*(1/(sqrt(2*pi*sigmab)))*deltab
    pval
}



#######################################################################################################
# The following pair of functions deal with scan of straddling pairs by their insert length.
# They are for insertions only.
#   One gets the threshold for a given p-value (i.e. FWER)
#   The other gives the p-value for a given threshold.
#
# The parameters are:
#   
#   b is the threshold (for pval...)
#   alpha is the FWER threshold (for getThreshold...).
#   n is the number of hypotheses tested.
#   r is the proportion of carriers.
#   w is the size of insertion (should be POSITIVE).  We usually don't know this, but the scan is not sensitive to this parameter.  w=10 usually works.
#   R1, R2: Pairs with + strand mapping within [0, t-R1] and - strand mapping within [t-R2,T] are counted.  
#       Note that R1 must be larter than R2!
#   kappa^2*T is the total number of mapped pairs where both reads are mapped, T is the effective length of the reference.
#   (mu, sigma) are mean and std dev of the insert length distribution (assumed Gaussian).
#   UPPER: technical parameter that users should not fiddle with.  It is the upper limit used in search for betab.
#   verbose: Do we want to print out status messages?
#################################################################################################################




getThresholdStraddlingReadsInsertionScan<-function(alpha,kappa,w,r,mu,sigma,n, R1, R2, UPPER=NA,verbose=TRUE,status=NULL){
    if(verbose) cat("\nComputing threshold for Straddling Reads Scan for FWER=",alpha,"...\n")
    #print(list(alpha,kappa,w,r,mu,sigma,n, R1, R2, UPPER,verbose,status))
    nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
    if(length(mu)!= nlibs || length(sigma) != nlibs) stop("The sizes of kappa, mu, and sigma need to match.\n")

    pvalminusalpha = function(b){ pvalStraddlingReadsInsertionScan(b=b,kappa=kappa,w=w,r=r,mu=mu,sigma=sigma,n=n, R1=R1, R2=R2, UPPER=UPPER,verbose=verbose,status=status)-alpha}
    
    
    m0=0; for(libi in 1:nlibs){ 
      if(status[libi])
        m0=m0+psidotbetaStraddlingReadsInsertionScan(beta=0,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2)}
    s0=0; for(libi in 1:nlibs){ 
      if(status[libi])
        s0=s0+sigmabetaStraddlingReadsInsertionScan(beta=0,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2, varZt=TRUE)}
    s0=sqrt(s0)
    
    cat("Straddling reads insertion scan: Null mean = ",m0,", standard deviation = ",s0,".\n")
    
    #m0=-0.493113,s0=0.2901709
    K=6
    prevroot = m0
    lower = m0+0.1*s0
    while(TRUE){
        upper = m0+K*s0
        if(verbose) cat("\nTrying to find threshold with search upper limit ", K," sd from mean: [",lower,", ",upper,"]\n",sep="")
        res=try(uniroot(pvalminusalpha, lower=m0+0.1*s0, upper=m0+K*s0), silent=TRUE)
        if(class(res) == "try-error"){
            if(verbose) cat("Upper limit too low, multiply K by 1.5.\n")
            K = K*1.5
        } else {
            if(verbose) cat("\n-----> uniroot converged on b=",res$root,".\n")
            break
        }
        if(!is.finite(upper)) { res=list(root=m0+6*s0); break } #just make safe return for now
    }
    if(verbose) cat("Threshold value: ",res$root,"\n")
    res$root
}

pvalStraddlingReadsInsertionScan<-function(b,kappa,w,r,mu,sigma,n, R1, R2, UPPER=NA,verbose=TRUE,status=NULL){

    if(verbose) cat("\nComputing P-VALUE for Straddling Pairs Scan with b=",b,"...\n")
    nlibs = length(kappa); if(is.null(status)) status=rep(TRUE,nlibs)
    if(length(mu)!= nlibs || length(sigma) != nlibs) stop("The sizes of kappa, mu, and sigma need to match.\n")

    psidotminusb = function(beta){ 
        psidot=0; for(libi in 1:nlibs){ 
          if(status[libi])
            psidot = psidot+psidotbetaStraddlingReadsInsertionScan(beta=beta,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2)}
        return(psidot-b)
    }
    if(is.na(UPPER)){
        UPPER = 2
    }
    if(psidotminusb(0)>=0){
        # Must start with psidotminusb(0)<0.
        stop("Threshold value must be larger than null mean.\n")
    }
    prevroot=0
    while(TRUE){
        if(verbose) cat("\tTrying to compute betab with search upper limit ", UPPER,"...\n",sep="")
        res=try(uniroot(psidotminusb, lower=0, upper=UPPER), silent=TRUE)
        if(class(res) == "try-error"){
            if(psidotminusb(UPPER)<0){
                # Must have returned error because psidotminusb(UPPER) is also negative.
                # uniroot must have psidotminusb(0)<0 (enforced above) and psidotminusb(UPPER)>0
                UPPER = UPPER*2
                if(verbose) cat("\tUpper limit too low and psidotminusb<0, increase upper to ", UPPER,".\n",sep="")
            } else{
                # Must have returned error because psidotminusb(UPPER) is infinity.
                if(verbose) cat("Upper limit too high, divide by 2.\n")
                UPPER = UPPER/1.5            
            }
        } else {
            break
        }
    }

    betab=res$root
    if(verbose) cat("Final value of beta_b: ",betab,".\n",sep="")
    
    sigmab = 0; for(libi in 1:nlibs){ 
      if(status[libi])
        sigmab=sigmab+sigmabetaStraddlingReadsInsertionScan(beta=betab,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2)}
    deltab = 0; for(libi in 1:nlibs){ 
      if(status[libi])
        deltab=deltab+deltabetaStraddlingReadsInsertionScan(beta=betab,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2)}
    psib = 0; for(libi in 1:nlibs){ 
      if(status[libi])
        psib=psib+psibetaStraddlingReadsInsertionScan(beta=betab,kappa=kappa[libi],w=w,r=r,mu=mu[libi],sigma=sigma[libi],R1=R1, R2=R2)}
    # cat("sigmab=",sigmab,", deltab=",deltab,", psib=",psib,sep="")
    
    pval = n* exp(-(betab*b-psib))*(1/(sqrt(2*pi*sigmab)))*deltab
    if(verbose) cat("P-value=",pval,", betab=",betab,", sigmab=",sigmab,", deltab=",deltab, ", psib=",psib,"\n")
    pval
}


#######################################################################################################
# The following function gives robust estimate of the mean, sd, and empirical distribution
# of the insert length distribution under null hypothesis
# of no structural variation.
# (rPstart,rMstart) are left most base of the plus and minus strand reads.
# UPPERBOUND is a bound on "reasonable" fragment sizes.  
#######################################################################################################

computeNullInsertLengthDistribution<-function(rPstart,rMstart,UPPERBOUND=NA, DELTA=NA){
    
    il=rMstart-rPstart
    ilmed=median(il, na.rm=TRUE)
    if(is.na(UPPERBOUND))   UPPERBOUND=10*ilmed      # Assumes that a fragment ten times the median is not possible.
    if(is.na(DELTA)) DELTA = UPPERBOUND/2
    
    p = 1-sum(il<DELTA & il>0, na.rm=TRUE)/length(il)
    pNA = sum(is.na(il))/length(il)
    il = il[is.finite(il)]
    il = pmax(pmin(il,UPPERBOUND,na.rm=TRUE),-100,na.rm=TRUE)
    il = il[which(il<DELTA & il>(-100))]
    if(length(il)==0) return(list(mu=NA,sigma=NA,f=NA))
    mu=mean(il)
    sigma=sd(il)
    
    ildens=density(il)

    list(mu=mu,sigma=sigma,f=f,ildens=ildens,pNA=pNA,p=p,DELTA=DELTA)
}



#######################################################################################################
# The following should be "internal functions".  They compute the intermediate values that go into a 
# p-value approximation: psibeta, psidotbeta, sigmabeta, deltabeta.
# The suffix show what scan they are for:  StraddlingReadsDeletionScan and HangingReadsScan.
#######################################################################################################

psibetaStraddlingReadsInsertionScan<-function(beta,kappa,w,r,mu,sigma,R1,R2,IntLim=6){
    # R1 must be larger than R2
    # w is the length of the insertion and should be POSITIVE.
    
    Delta = IntLim*sigma
    upperlim = mu+Delta
    f = function(z){ dnorm(z,mean=mu,sd=sigma) }
    wprime = R1-R2 
    lw=function(z){ f(z+w)/f(z)}
    
    integrand = function(z){(z-wprime)*((1-r+r*lw(z))^beta-1)*f(z)}
    res= integrate(integrand, lower=wprime, upper=upperlim)
    term = kappa^2*res$value
    term
    
}
psidotbetaStraddlingReadsInsertionScan<-function(beta,kappa,w,r,mu,sigma,R1,R2,IntLim=6){
        
    Delta = IntLim*sigma
    upperlim = mu+Delta
    f = function(z){ dnorm(z,mean=mu,sd=sigma) }
    wprime = R1-R2 
    lw=function(z){ f(z+w)/f(z)}
    
    integrand = function(z){(z-wprime)*(log(1-r+r*lw(z)))*(1-r+r*lw(z))^beta*f(z)}
    res= integrate(integrand, lower=wprime, upper=upperlim)
    term = kappa^2*res$value
    term
}


sigmabetaStraddlingReadsInsertionScan<-function(beta,kappa,w,r,mu,sigma,R1,R2,IntLim=6, varZt=FALSE){
    Delta = IntLim*sigma
    upperlim = mu+Delta
    f = function(z){ dnorm(z,mean=mu,sd=sigma) }
    wprime = R1-R2 
    lw=function(z){ f(z+w)/f(z)}
   
    f2 = function(z){
        (z-wprime)* (log(1-r+r*lw(z)))^2*(1-r+r*lw(z))^beta*f(z)
    }
    res = integrate(f2, lower=wprime, upper=upperlim)
    
    if(!varZt){
        sigma2 = beta^2*kappa^2*res$value
    } else {
        sigma2 = kappa^2*res$value
    }
    sigma2
}

deltabetaStraddlingReadsInsertionScan<-function(beta,kappa,w,r,mu,sigma,R1,R2,IntLim=6){
    Delta = IntLim*sigma
    upperlim = mu+Delta
    upperlim2 = mu+Delta+Delta
    wprime = R1-R2
    
    f = function(z){dnorm(z,mean=mu,sd=sigma) }
    lw = function(z){f(z+w)/f(z)}
    c = function(z){log(1-r+r*lw(z))}
       
    f2=function(z){
    (1-exp(beta*c(z)))*c(z)*f(z)
    }
    res = integrate(f2,lower=wprime,upper=upperlim2)
    B1mA1 = kappa^2*res$value
    delta = -beta*B1mA1
    delta 
}


psibetaStraddlingReadsDeletionScan<-function(beta,kappa,w,r,f=NULL,mu=NULL,sigma=NULL,R1,R2,IntLim=6,doplots=FALSE){
    if(is.null(mu) && is.null(sigma) && is.null(f)){
        stop("Must specify either the empirical IL distribution (f) or its mean (mu) and standard deviation (sigma).")
    }
    
    Delta = IntLim*sigma
    upperlim = mu+Delta+w
    wprime = w+R1-R2 
    
    x=seq(wprime,upperlim,1)
    psibetaf=NA;   psibetadnorm=NA
    if(!is.null(f)){
        lw=function(z){ f(z-w)/f(z)}
        integrand = function(z){(z-wprime)*((1-r+r*lw(z))^beta-1)*f(z)}
        val = integrand(x)
        psibetaf = kappa^2*sum(val)
    }
    if((!is.null(mu)) && (!is.null(sigma))){
        f2 = function(z){ dnorm(z,mean=mu,sd=sigma)}
        lw=function(z){ f2(z-w)/f2(z)}
        integrand2 = function(z){(z-wprime)*((1-r+r*lw(z))^beta-1)*dnorm(z,mean=mu,sd=sigma)}
        val2 = integrand2(x)
        psibetadnorm = kappa^2*sum(val2)  
    }
    if(!(is.null(f) || is.null(mu) || is.null(sigma)) && doplots){
        par(mfrow=c(2,1))
        plot(x,f(x), xlab="Insert Length", ylab="Density")
        lines(x,dnorm(x,mean=mu, sd=sigma), col="blue")
        plot(x,val, ylim=c(min(c(val,val2)),max(c(val,val2))),xlab="Insert Length", ylab="Integrand Value")
        lines(x,val2,col="blue")
        cat("The analytical psibeta computed assuming normal IL density is: ",psibetadnorm,".\n",
        "The analytical psibeta computed using empirical IL density is: ",psibetaf,".\n")
    }
    
    if(is.na(psibetaf)) return(psibetadnorm)
    if(is.na(psibetadnorm)) return(psibetaf)
    return(list(psibetaf=psibetaf,psibetadnorm=psibetadnorm))
}

psidotbetaStraddlingReadsDeletionScan<-function(beta,kappa,w,r,f=NULL,mu=NULL,sigma=NULL,R1,R2,IntLim=6,doplots=FALSE){
    if(is.null(mu) && is.null(sigma) && is.null(f)){
        stop("Must specify either the empirical IL distribution (f) or its mean (mu) and standard deviation (sigma).")
    }
    
    Delta = IntLim*sigma
    upperlim = mu+Delta+w
    wprime = w+R1-R2 

    x=seq(wprime,upperlim,1)
    analmeanf=NA;   analmeandnorm=NA
    if(!is.null(f)){
        lw=function(z){ f(z-w)/f(z)}
        integrand = function(z){(z-wprime)*(log(1-r+r*lw(z)))*(1-r+r*lw(z))^beta*f(z)}
        val = integrand(x)
        analmeanf = kappa^2*sum(val)
    }
    if((!is.null(mu)) && (!is.null(sigma))){
        f2 = function(z){ dnorm(z,mean=mu,sd=sigma)}
        lw=function(z){ f2(z-w)/f2(z)}
        integrand2 = function(z){(z-wprime)*(log(1-r+r*lw(z)))*(1-r+r*lw(z))^beta*dnorm(z,mean=mu,sd=sigma)}
        val2 = integrand2(x)
        analmeandnorm = kappa^2*sum(val2)  
    }
    
     
    if(!(is.null(f) || is.null(mu) || is.null(sigma)) && doplots){
        par(mfrow=c(2,1))
        plot(x,f(x), type="l",xlab="Insert Length", ylab="Density")
        lines(x,dnorm(x,mean=mu, sd=sigma), col="blue")
        legend(x="topright",lty=c(1,1),col=c("black","blue"),legend=c("Empirical","Normal"))
        plot(x,val, type="l",ylim=c(min(c(val,val2)),max(c(val,val2))),main=paste("Psidot(Beta=",beta,")",sep=""),xlab="Insert Length", ylab="Integrand Value")
        lines(x,val2,col="blue")
        abline(0,0,col="gray")
        legend(x="topright",lty=c(1,1),col=c("black","blue"),legend=c("Empirical","Normal"))
#        cat("The mean of empirical Z process is ",mean(resDel$ilTerm), ", its median is ",median(resDel$ilTerm),".\n",
#        "Comparing this to analytical values:\n",
        cat("The analytical mean computed assuming normal IL density is: ",analmeandnorm,".\n",
        "The analytical mean computed using empirical IL density is: ",analmeanf,".\n",sep="")
    }
    
    if(is.na(analmeanf)) return(analmeandnorm)
    if(is.na(analmeandnorm)) return(analmeanf)
    return(list(analmeanf=analmeanf,analmeandnorm=analmeandnorm))
}

sigmabetaStraddlingReadsDeletionScan<-function(beta,kappa,w,r,f=NULL,mu=NULL,sigma=NULL,R1,R2,IntLim=6,varZt=FALSE,doplots=FALSE){
    if(is.null(mu) && is.null(sigma) && is.null(f)){
        stop("Must specify either the empirical IL distribution (f) or its mean (mu) and standard deviation (sigma).")
    }
    
    Delta = IntLim*sigma
    upperlim = mu+Delta+w
    wprime = w+R1-R2 

    x=seq(wprime,upperlim,1)
    analvarf=NA;   analvardnorm=NA
    if(!is.null(f)){
        lw=function(z){ f(z-w)/f(z)}
        integrand = function(z){(z-wprime)* (log(1-r+r*lw(z)))^2*(1-r+r*lw(z))^beta*f(z)}
        val = integrand(x)
        analvarf = kappa^2*sum(val)
    }
    if((!is.null(mu)) && (!is.null(sigma))){
        f2 = function(z){ dnorm(z,mean=mu,sd=sigma)}
        lw=function(z){ f2(z-w)/f2(z)}
        integrand2 = function(z){(z-wprime)* (log(1-r+r*lw(z)))^2*(1-r+r*lw(z))^beta*dnorm(z,mean=mu,sd=sigma)}
        val2 = integrand2(x)
        analvardnorm = kappa^2*sum(val2)  
    }
    if(!(is.null(f) || is.null(mu) || is.null(sigma)) && doplots){
        par(mfrow=c(2,1))
        plot(x,f(x), xlab="Insert Length", ylab="Density")
        lines(x,dnorm(x,mean=mu, sd=sigma), col="blue")
        legend(x="topright",col=c("black","blue"),legend=c("Empirical","Normal"))
        plot(x,val, ylim=c(min(c(val,val2)),max(c(val,val2))),main=paste("Var(beta=",beta,")",sep=""),xlab="Insert Length", ylab="Integrand Value")
        lines(x,val2,col="blue")
        legend(x="topright",col=c("black","blue"),legend=c("Empirical","Normal"))        
#        cat("The var of empirical Z process is ",var(resDel$ilTerm), ".\n",
#        "Comparing this to analytical values:\n",
        cat("The analytical var computed assuming normal IL density is: ",analvardnorm,".\n",
        "The analytical var computed using empirical IL density is: ",analvarf,".\n",sep="")
    }

    if(!varZt){
        analvardnorm = beta^2*analvardnorm
        analvarf = beta^2*analvarf
    } 
        
    if(is.na(analvarf)) return(analvardnorm)
    if(is.na(analvardnorm)) return(analvarf)
    return(list(analvarf=analvarf,analvardnorm=analvardnorm))
}

deltabetaStraddlingReadsDeletionScan<-function(beta,kappa,w,r,f=NULL,mu=NULL,sigma=NULL,R1,R2,IntLim=6,doplots=FALSE){
    if(is.null(mu) && is.null(sigma) && is.null(f)){
        stop("Must specify either the empirical IL distribution (f) or its mean (mu) and standard deviation (sigma).")
    }
    Delta = IntLim*sigma
    upperlim = mu+2*Delta+w
    wprime = w+R1-R2 

    x=seq(wprime,upperlim,1)
    deltaf=NA;   deltadnorm=NA
    if(!is.null(f)){
        lw=function(z){ f(z-w)/f(z)}
        g = function(z){log(1-r+r*lw(z))}
        integrand = function(z){(1-exp(beta*g(z)))*g(z)*f(z)}
        val = integrand(x)
        deltaf = -beta*kappa^2*sum(val)
    }
    if((!is.null(mu)) && (!is.null(sigma))){
        f2 = function(z){ dnorm(z,mean=mu,sd=sigma)}
        lw=function(z){ f2(z-w)/f2(z)}
        g = function(z){log(1-r+r*lw(z))}
        integrand2 = function(z){(1-exp(beta*g(z)))*g(z)*dnorm(z,mean=mu,sd=sigma)}
        val2 = integrand2(x)
        deltadnorm = -beta*kappa^2*sum(val2)  
    }
    if(!(is.null(f) || is.null(mu) || is.null(sigma)) && doplots){
        par(mfrow=c(2,1))
        plot(x,f(x), xlab="Insert Length", ylab="Density")
        lines(x,dnorm(x,mean=mu, sd=sigma), col="blue")
        plot(x,val, ylim=c(min(c(val,val2)),max(c(val,val2))),xlab="Insert Length", ylab="Integrand Value")
        lines(x,val2,col="blue")
        cat("The analytical delta computed assuming normal IL density is: ",deltadnorm,".\n",
        "The analytical delta computed using empirical IL density is: ",deltaf,".\n")
    }
        
    if(is.null(f)) return(deltadnorm)
    if(is.null(mu) || is.null(sigma)) return(deltaf)
    return(list(deltaf=deltaf,deltadnorm=deltadnorm))
}

psidotbetaHangingReadsScan<-function(beta,kappa,w,r,p,pbase,mu,sigma,D,R,F=NULL){
    if(is.null(F)) F = function(z){1-pnorm(z, mean=mu,sd=sigma)}
    g = function(z){log(1+(r*(1-p)/p)*F(z))}
    f2 = function(x){
        g(D-x)*exp(beta*g(D-x))
    }
    res = integrate(f2, lower=0, upper=D-R)
    psidot = pbase*res$value
    psidot
}
psibetaHangingReadsScan<-function(beta,kappa,w,r,p,pbase,mu,sigma,D,R,F=NULL){
    if(is.null(F)) F = function(z){1-pnorm(z, mean=mu,sd=sigma)}
    g = function(z){log(1+(r*(1-p)/p)*F(z))}
    f2 = function(x){
        exp(beta*g(D-x))-1
    }
    res = integrate(f2, lower=0, upper=D-R)
    psi = pbase*res$value
    psi
}

sigmabetaHangingReadsScan<-function(beta,kappa,w,r,p,pbase,mu,sigma,D,R,varZt=FALSE,F=NULL){
    if(is.null(F)) F = function(z){1-pnorm(z, mean=mu,sd=sigma)}
    g = function(z){log(1+(r*(1-p)/p)*F(z))}
    f2 = function(x){
        (g(D-x)^2)*exp(beta*g(D-x))
    }
    res = integrate(f2, lower=0, upper=D-R)
    
    if(!varZt){
        sigmabeta = beta^2*pbase*res$value
    } else {
        sigmabeta = pbase*res$value
    }
    sigmabeta
}

deltabetaHangingReadsScan<-function(beta,kappa,w,r,p,pbase,mu,sigma,D,R,F=NULL){
    if(is.null(F)) F = function(z){1-pnorm(z, mean=mu,sd=sigma)}
    f = function(z){dnorm(z,mean=mu,sd=sigma)}
    g = function(z){log(1+(r*(1-p)/p)*F(z))}
    gprime = function(z){-(r*(1-p)/p)*f(z)/(1+(r*(1-p)/p)*F(z))}
    f2 = function(x){
        gprime(D-x)*exp(beta*g(D-x))
    }
    res = integrate(f2, lower=0, upper=D-R)
    deltabeta=beta*pbase*(g(D)*exp(beta*g(D)) - res$value -g(R))   
    deltabeta
}

