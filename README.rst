README
========

INTRODUCTION
--------------
  Statistical Structural Variant scAN(S2VAN, SWAN)
  Currently the package works for Linux (tested for Ubuntu and CentOS) and Mac (with Macports).
  It might also work for Windows with Cygwin (not tested).
  SWAN documentation is on wiki and it is available:
  `Wiki <http://bitbucket.org/charade/swan/wiki>`_

DEPENDENCIES
--------------

  GCC(>4.7)
        through apt-get(Ubuntu), yum(CentOS), macports(OSX), homebrew(OSX) 
  Boost C++ Library
        `boost download <http://www.boost.org>`_
  R(>=3.2)
        `R download <http://www.r-project.org>`_
  Samtools(>=1.2)
        `Samtools download <http://www.samtools.org>`_
  CRAN R Libraries
        RcppArmadillo (source), Rcpp (source);
        BH, data.table, digest, hash, methods, optparse, parallel, plyr, robustbase, sets, stringr, zoo
  BioConductor R Libraries
        Biobase, Biostrings, BSgenome, GenomeInfoDb, GenomicRanges, IRanges, Rsamtools, S4Vectors

INSTALL
-------------
  
  Assuming C++, Boost, R, devtools and Samtools are already properly installed; 
  $BOOST_ROOT properly set. 

  **Install R Package Dependencies**
  
  :: 

    # preset CRAN mirror to prevent interruption
    R> local({r <- getOption("repos"); r["CRAN"] <- "http://cran.us.r-project.org"; options(repos=r)}) 
    # Some Rcpp packages have to to installed from source, otherwise cause 'segfault'
    R> install.packages(pkgs=c("Rcpp","RcppArmadillo"),type="source") 
    # 'lgfortran' and 'lquadmath' may affect OS X, fix by:  
    sh> curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
    sh> sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
    R> install.packages(pkgs=c("BH", "Cairo", "data.table", "digest", 
    "ggplot2", "gridExtra", "hash", "methods", "optparse", "parallel", "plyr", "robustbase", 
    "seqCBS", "sets", "stringr", "zoo"))  # other CRAN packages 
    R> source("http://bioconductor.org/biocLite.R")      #Bioconductor
    R> biocLite(pkgs=c("Biobase", "Biostrings", "BSgenome", "GenomeInfoDb", 
    "GenomicRanges", "IRanges", "Rsamtools","S4Vectors"))   # other Bioconductor packages
  
  **Install SWAN**
  
  ::

    R> devtools::install_bitbucket("charade/swan",dependencies=T,noCache=T) 
  
  Note the executables will be available from $R_LIBS_USER/library/swan/bin.
  User can export $SWAN_BIN=$R_LIBS_USER/library/swan/bin and add it to $PATH.

EXECUTABLES
------------

  $SWAN_BIN/swan_stat         --  pre-scan lib-wise sequencing statistics

  $SWAN_BIN/swan_scan         --  genome-wide likelihood scan

  $SWAN_BIN/sclip_scan        --  genome-wide soft-sclip scan

  $SWAN_BIN/swan_join         --  merging evidence from multiple features


USAGE
--------
  (1) Use '-h' to read script-wise usage. 

  (2) Do a Sanity check for installation and learn single or paired sample analysis pipelines.



  ::

    # download the mock swan_test data package (approximately 603MB)
    sh> wget http://meta.usc.edu/softs/swan/swan_test.tgz
    # unzip it within this directory 
    sh> tar -zxvf swan_test.tgz
    # you will see an example directory containing necessary mock data files for successful testing
    sh> $SWAN_BIN/single.sh all
    sh> $SWAN_BIN/paired.sh all
  
WIKI
--------
  http://bitbucket.org/charade/swan/wiki/Home
  
FAQ
--------
  http://bitbucket.org/charade/swan/wiki/FAQ
  
BUG
--------
  https://bitbucket.org/charade/swan/issues

CONTACT
--------
  lixia at stanford dot edu
