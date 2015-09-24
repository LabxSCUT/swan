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

  C++ build environment
        e.g. build-essential and libstdc++ in Ubuntu 
  Boost C++ Library
        `boost download <http://www.boost.org>`_
  R(>=3.1)
        `R download <http://www.r-project.org>`_
  Samtools(>=0.19)
        `Samtools download <http://www.samtools.org>`_
  R Libraries
        ggplot2, gridExtra, stringr, data.table, Biostrings,
        Rsamtools, optparse, GenomicRanges, multicore, BSgenome,
        Rcpp, inline, IRanges, rbenchmark, RcppArmadillo, methods, zoo, robustbase
  
  To setup the dependencies, users may refer to the author's development document:
  `Tutorial <http://dl.dropbox.com/u/35182955/development_environment.html>`_

INSTALL
-------------
  
  1. [R Package Dependencies]
  Before install SWAN please install all above dependencies required.
  Add following R packages. To spot errors easily, better install the dependencies
  one by one.
  
  :: 

    install.packages("Rcpp",type="source") #R packages
    install.packages("RcppArmadillo",type="source")
    install.packages("data.table",type="source")
    install.packages("optparse",type="source")
    install.packages("robustbase",type="source")
    install.packages("BH",type="source")
    install.packages("stringr",type="source")
    install.packages("digest",type="source")
    install.packages("gridExtra",type="source")
    install.packages("sets",type="source")
    install.packages("colorspace",type="source")
    install.packages("plyr",type="source")
    install.packages("hash",type="source")
    install.packages("ggplot2",type="source") #you can ignore if not plotting with SWAN
    install.packages("Cairo",type="source") #need install cairo-dev first in ubuntu, you can ignore if not plotting with SWAN
    
    source("http://bioconductor.org/biocLite.R") #BioCLite packages
    biocLite("BiocUpgrade")
    biocLite("BiocGenerics",type="source",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("seqCBS",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("Biobase",type="source",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("S4Vectors",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("IRanges",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("XVector",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("GenomicRanges",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("BSgenome",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("Biostrings",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
    biocLite("Rsamtools",ask=F,suppressUpdates=F,suppressAutoUpdate=F)
  
  2. [User]
  Download the latest master branch of SWAN from `master <https://bitbucket.org/charade/swan/downloads>`_.
  select branches and download the master branch in either .gz, .zip or .bz format
  unzip the archive to a folder say swan, create the ~/scripts folder.
  set BOOST_HOME to be your /boost/path/ and install the swan package.
  
  ::

    > mkdir ~/scripts
    > BOOST_HOME=/boost/path/ R CMD INSTALL swan --preclean
  
  Now the executables will be available from ~/scripts

  3. [Developer]
  Download the latest test branch of SWAN from `devel <https://bitbucket.org/charade/swan/downloads>`_.
  select branches and download the master branch in either .gz, .zip or .bz format
  unzip the archive to a folder say swan, create the ~/scripts folder.
  set BOOST_HOME to be your /boost/path/ and install the swan package
  
  ::

    > git clone https://bitbucket.org/charade/swan.git
    > git checkout testing
    > mkdir ~/scripts
    > BOOST_HOME=/boost/path/ R CMD INSTALL swan --preclean

EXECUTABLES
------------

  $SWAN_BIN/swan_stat         --  pre-scan lib-wise sequencing statistics

  $SWAN_BIN/swan_scan         --  genome-wide likelihood scan

  $SWAN_BIN/sclip_scan        --  genome-wide soft-sclip scan

  $SWAN_BIN/swan_join         --  merging evidence from multiple features


USAGE
--------
  (1) By default all above executables will be available from ~/scripts .
  Use '-h' to read script-wise usage.

  (2) Do a Sanity check for installation and learn single or paired sample analysis from toy examples.
  require installatoin of SVEngine (http://bitbucket.org/charade/svengine).

  ::

    cd test
    ./single.sh all
    ./paired.sh all
  
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
