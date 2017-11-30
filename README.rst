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

  Currently the package works for Linux (tested for Ubuntu).

  It might also work for Mac (with Macports and Homebrew) and Windows (with Cygwin). 
  These are not tested.

  Active SWAN documentation effort is on SWAN Wiki:
  `Wiki <http://bitbucket.org/charade/swan/wiki>`_

DEPENDENCIES
--------------

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

DOCKER
-------------

	Due to the multiple R and Python dependencies involved,
  the easiest way to use SWAN is by the provided docker image. 
  If you have a docker server running, 
  just git clone the swan package to a local dir say $SWAN_HOME.
  From $SWAN_HOME 
    
  

INSTALL
-------------
  
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
    
    sh> export SWAN_BIN=/Users/charlie/Library/R/3.2/library/swan/bin
    sh> $SWAN_BIN/swan_test.sh

  Afterwards, the executables can be moved to other places as the user need and the user need to update $SWAN_BIN and $PATH accordingly.

  **Use SWAN (Without Install)**
  
  You can use  Ubuntu or CentOS virtual machines with SWAN pre-installed - easily deployable to cloud computing. 
  The virtual machine disk images can be found here:
  Ubuntu: http://meta.usc.edu/softs/vbox/Ubuntu_14_SWAN.vdi.gz;
  CentOS: http://meta.usc.edu/softs/vbox/CentOS_7_SWAN.vdi.gz.
  Oracle's free VirtualBox (https://www.virtualbox.org/) among other softwares can be used to load the images. 
  There are numerous how-to tutorials on Youtube about how to import .vdi to VirtualBox, 
  for example this one (https://www.youtube.com/watch?v=fLyriYu0lU0). Once the virtual machine is running,
  you can login with account: **user** and password: **user** and refer to the "README.rst" file on the desktop
  to proceed. The "action.log" file also contains full commands that required to setup SWAN on the virtual machine. 
  
  Details regarding this **INSTALL** section can be also found in the PPT slides here: `tutorial <http://bitbucket.org/charade/swan/wiki/doc/SWAN_Installation.pptx>`_

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
