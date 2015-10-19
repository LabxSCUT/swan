README
========

INTRODUCTION
--------------
  Statistical Structural Variant Analysis for NGS (S2VAN, SWAN)
  SWAN is an RcppArmadillo package. Please also read general online information to install R and Rcpp packages before proceed. 
 
  There is a detailed PPT tutorial for install SWAN on Linux/OSX: 
  `tutorial <http://bitbucket.org/charade/swan/wiki/doc/SWAN_Installation.pptx>`_

  Many of the installation questions are also answered in FAQ:
  `FAQ <http://bitbucket.org/charade/swan/wiki/FAQ>`_

  Currently the package works for Linux (tested for Ubuntu and CentOS) and Mac (with Macports and Homebrew).
  It might also work for Windows with Cygwin (not tested).
  Active SWAN documentation effort is on SWAN Wiki:
  `Wiki <http://bitbucket.org/charade/swan/wiki>`_

DEPENDENCIES
--------------

  GCC(>=3.4)
        through apt-get(Ubuntu), yum(CentOS), macports(OSX), homebrew(OSX) 
  R(>=3.1)
        `R download <http://www.r-project.org>`_
  Samtools(>=0.19)
        `Samtools download <http://www.samtools.org>`_
  CRAN R Libraries
        RcppArmadillo (source), Rcpp (source);
        BH, data.table, digest, devtools, hash, methods, optparse, parallel, plyr, robustbase, sets, stringr, zoo
  BioConductor R Libraries
        Biobase, Biostrings, BSgenome, GenomeInfoDb, GenomicRanges, IRanges, Rsamtools, S4Vectors

INSTALL
-------------
  
  Following installation process assumes: (1) C++, R, devtools and Samtools are already properly installed; (2) Correct libstdc++ is either in user's or system's $LD_LIBRARY_PATH.

  **Install R Package Dependencies**
  
  :: 

    # preset CRAN mirror to prevent interruption
    R> local({r <- getOption("repos"); r["CRAN"] <- "http://cran.us.r-project.org"; options(repos=r)}) 
    # Some Rcpp packages have to to installed from source, otherwise may cause runtime 'segfault'
    R> install.packages(pkgs=c("Rcpp","RcppArmadillo"),type="source") 
    # 'lgfortran' and 'lquadmath' may affect OS X, fix by:  
    sh> curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
    sh> sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
    R> install.packages(pkgs=c("BH", "data.table", "digest", "hash", "methods", "optparse", "parallel", "plyr", "robustbase", "sets", "stringr", "zoo"))  # other CRAN packages 
    R> source("http://bioconductor.org/biocLite.R")      #Bioconductor
    R> biocLite(pkgs=c("Biobase", "Biostrings", "BSgenome", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools","S4Vectors"))   # other Bioconductor packages
  
  **Install SWAN**
  
  ::

    R> devtools::install_bitbucket("charade/swan",dependencies=T,clean=T) 
  
  **Test SWAN**

  Note by default the SWAN executables will be available from path: $SWAN_BIN=$R_LIBS_USER/library/swan.
  However, the exact naming of the $R_LIBS_USER is system and user specific and can only be determined at install time.
  The path will show up in the final '#' surrounded banner looks like below:

  ::

    #####################################
    #
    #  Your SWAN Binaries can be found at:
    #  /Users/charlie/Library/R/3.2/library/swan/bin
    #  To use SWAN, set environment variable $SWAN_BIN to above path
    #  And add $SWAN_BIN to your $PATH evironment
    #
    #####################################
  
  In this case, to run the test scripts, the user should export $SWAN_BIN=/Users/charlie/Library/R/3.2/library/swan/bin and add this $SWAN_BIN to $PATH.
  Now, do a Sanity check for installation and learn single or paired sample analysis pipelines.

  ::
    
    export SWAN_BIN=/Users/charlie/Library/R/3.2/library/swan/bin
    $SWAN_BIN/swan_test.sh

  Afterwards, the executables can be moved to other places as the user need and the user need to update $SWAN_BIN and $PATH accordingly.

  **Use SWAN**
  
  You can use SWAN with pre-installed Ubuntu or CentOS virtual machines easily deployable to clouds. The virtual machine disk images can be found here:
  Ubuntu: http://meta.usc.edu/softs/vbox/Ubuntu_14_SWAN.vdi.gz
  CentOS: http://meta.usc.edu/softs/vbox/CentOS_7_SWAN.vdi.gz
  Oracle's free VirtualBox (https://www.virtualbox.org/) among others can be used to load the images. There are numerous how-to tutorials on Youtube, 
  for example this one (https://www.youtube.com/watch?v=fLyriYu0lU0).

EXECUTABLES
------------

  $SWAN_BIN/swan_stat         --  pre-scan lib-wise sequencing statistics

  $SWAN_BIN/swan_scan         --  genome-wide likelihood scan

  $SWAN_BIN/sclip_scan        --  genome-wide soft-sclip scan

  $SWAN_BIN/swan_join         --  merging evidence from multiple features


USAGE
--------
  Use '-h' to read script-wise usage. 
  
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