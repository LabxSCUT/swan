.. |Logo| image:: https://bitbucket.org/charade/swan/raw/master/doc/images/swan_logo.png
   :alt: logo.png
   :height: 50px
   :width: 100px

.. |Rationale| image:: https://bitbucket.org/charade/swan/raw/master/doc/images/swan_rationale.png
   :alt: rationale.png
   :height: 450px
   :width: 540px

|Logo| SWAN - Statistical Structural Variant Analysis for NGS (S2VAN, SWAN)
==========================================================================================

QUICK LINKS
-----------

`Examples <https://bitbucket.org/charade/swan/wiki/Example>`__

`Manuals <https://bitbucket.org/charade/swan/wiki/Manual>`__

`FAQ <https://bitbucket.org/charade/swan/wiki/FAQ>`__


README
========

INTRODUCTION
--------------

  SWAN - Statistical Structural Variant Analysis for NGS (S2VAN, SWAN)

  - SWAN analyzes whole genome sequencing data sets to identify structural variants.

  - SWAN currently supports shot-gun whole genome sequencing data from Illumina platforms.

  - SWAN takes as input the BAM file output of BWA aligner or others and outputs a list of structural variants in BED/VCF format.

  - SWAN's modeling and pipeline are illustrated in Figures 1. 

  Figure 1. The principle and workflow of SWAN.

  |Rationale|

  - SWAN is an RcppArmadillo package. It works OS platform independently (through Docker). 
    It installs on Linux (tested for Ubuntu). 
    It might also install for Mac (with Macports and Homebrew) and Windows (with Cygwin). 
    However, these were not tested.

DOCKER (Platform Indepedent)
-------------------------------

  Due to the multiple R and Python dependencies involved,
  the easiest way to use SWAN is by the provided docker image build file.
  A Dockerfile is provided to build SWAN enabled docker image from a standard Ubuntu docker image.
  If you are not familiar with Docker, it is a container platform widely used in industry/academia. 
  Here is the link to the Docker community:
    `docker community <https://www.docker.com>`_ .
  If you have a docker server running, 
  just need to download the Dockerfile from: 
    `dockerfile <https://bitbucket.org/charade/swan/raw/master/Dockerfile>`_
  into $your_swan_container and run:

  ::

    docker build --no-cache $your_swan_container

DEPENDENCIES (All Platforms)
---------------------------------

  GCC(>=3.4)
        through apt-get(Ubuntu), yum(CentOS), macports(OSX), homebrew(OSX) 

  R(>=3.1)
        `R download <http://www.r-project.org>`_

  Samtools(>=0.20)
        `Samtools download <http://www.samtools.org>`_

  CRAN R Libraries
        RcppArmadillo (source), Rcpp (source),
        data.table, devtools, digest, hash, methods, optparse, parallel, plyr, 
        robustbase, sets, stringr, zoo

  BioConductor R Libraries
        Biobase, Biostrings, BSgenome, GenomeInfoDb, GenomicRanges, IRanges, Rsamtools

INSTALL (Linux/Ubuntu)
--------------------------
  
  Following installation process assumes: 
  (1) GCC(>=4.3), R(>=3.1), Samtools(>=0.20) are already properly installed and in your $PATH; 

  **Install R Package Dependencies**
  
  :: 

    # First disable slow Tk/Tcl prompts of mirrors
    R> options(menu.graphics=FALSE)
    # Some Rcpp packages have to to installed from source, otherwise may cause runtime 'segfault'
    R> install.packages(pkgs=c("Rcpp","RcppArmadillo"),type="source") 
    # If you have "-lgfortran" or "-lquadmath" not found problems from above commands, please see entry in FAQ for fix. It mostly affects Ubuntu<=12, where the libgfortran link is often broken. 
    R> install.packages(pkgs=c("data.table", "devtools", "digest", "hash", "methods", "optparse", "parallel", "plyr", "robustbase", "sets", "stringr", "zoo"))  # other CRAN packages 
    R> source("http://bioconductor.org/biocLite.R")      #Bioconductor
    R> biocLite("BiocUpgrade") #Upgrade your Bioc to latest version compatible with your R version
    # now if you have "Error in unloadNamespace(package)" after "preparing package for lazy loading", please see entry in FAQ for fix. It is most likely R sessions haven't finished updating packages, try reinstall SWAN with a new Shell and R session some time later and it will self correct.
    R> biocLite(pkgs=c("Biobase", "Biostrings", "BSgenome", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools"))   # other Bioconductor packages
    # now if you see warnings or errors during installation of any above packages, try the above two steps again and it usually self resolves.
  
  **Install SWAN**
  
  ::

    R> library(devtools)
    R> devtools::install_bitbucket("charade/swan",dependencies=T,clean=T) 
  
  **Test SWAN**

  Note by default the SWAN executables will be available from path: $SWAN_BIN=$R_LIBS_USER/library/swan.
  However, the exact naming of the $R_LIBS_USER is system and/or user specific and can only be determined at the install time.
  Your $SWAN_BIN path will show up in the final '#' surrounded banner looks like below:

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
  Alternatively, you might want to install R package through shell that you can pre-specify $SWAN_BIN before installation. For example, to install swan to SWAN_BIN=$HOME/setup/swan/inst,

  ::
  
    bash> cd $HOME/setup && git clone https://bitbucket.org/charade/swan.git
    bash> export SWAN_BIN=$HOME/setup/swan/inst && cd $HOME/setup && R CMD INSTALL swan

  After installation, please do a sanity check for and learn the usage of single or paired sample analysis pipelines.

  ::
    
    bash> export SWAN_BIN=/Users/charlie/Library/R/3.2/library/swan/bin
    bash> $SWAN_BIN/swan_test.sh $SWAN_BIN

  If the executables were moved to other places and the user has to update $SWAN_BIN and $PATH accordingly.

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

CITATIONS
-----------------

  1. Li C. Xia, Sukolsak Sakshuwong, Erik Hopmans, John Bell, Susan Grimes, David Siegmund, Hanlee Ji, Nancy Zhang. A genome-wide approach for detecting novel insertion-deletion variants of mid-range size. Nucleic Acids Research (2016) 44 (15): e126 (https://doi.org/10.1093/nar/gkw481)
  2. Nancy R Zhang, Benjamin Yakir, Li C Xia, David Siegmund. Scan statistics on poisson random fields with applications in genomics. Annals of Applied Statistics (2016) Volume 10, Number 2 (2016), 726-755 (https://doi.org/10.1214/15-AOAS892)
