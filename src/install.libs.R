#put shared libs in PACKAGE/libs 
#options(warn=2) #warning will error
version=system(paste("cd",file.path(R_PACKAGE_SOURCE),"&& git log | head -n 1" ), intern=T)[1]
if(is.na(version)) version=scan(file.path(R_PACKAGE_SOURCE,"VERSION"),what="character",quiet=T)[1]
cat("Installing SWAN version=",version,"\n")

lib_files <- Sys.glob(paste("*", SHLIB_EXT, sep=''))
libarch <- if (nzchar(R_ARCH)) paste('libs',R_ARCH,sep="") else 'libs'
lib_dest <- file.path(R_PACKAGE_DIR, libarch)
#lib_build <- file.path(R_PACKAGE_SOURCE, "build")
dir.create(lib_dest, recursive = TRUE, showWarnings = FALSE)
#dir.create(lib_build, recursive = TRUE, showWarnings = FALSE)
n=length(lib_files)
message(paste(rep("Installing",n),lib_files,rep("to",n), rep(lib_dest,n),collapse="\n"))
file.copy(lib_files, lib_dest, overwrite = TRUE)
if(file.exists("symbols.rds"))
  file.copy("symbols.rds", lib_dest, overwrite = TRUE)

#put scripts in $RBIN or in ~/scripts
R_files <- Sys.glob(file.path(R_PACKAGE_SOURCE,"inst",paste("*", ".R", sep='')))
r_files <- Sys.glob(file.path(R_PACKAGE_SOURCE,"inst",paste("*", ".r", sep='')))
py_files <- Sys.glob(file.path(R_PACKAGE_SOURCE,"inst",paste("*", ".PY", sep='')))
sh_files <- Sys.glob(file.path(R_PACKAGE_SOURCE,"inst",paste("*", ".sh", sep='')))
test_scripts <- Sys.glob(file.path(R_PACKAGE_SOURCE,"test",paste("*", ".sh", sep='')))
R_builds <- unlist(sapply(R_files,function(x) {file.path(dirname(x),system(paste("basename",x,".R"),intern=T))}))
r_builds <- unlist(sapply(r_files,function(x) {file.path(dirname(x),system(paste("basename",x,".r"),intern=T))}))
py_builds <- unlist(sapply(py_files,function(x) {file.path(dirname(x),system(paste("basename",x,".PY"),intern=T))}))
sh_builds <- unlist(sapply(sh_files,function(x) {file.path(dirname(x),system(paste("basename",x,".sh"),intern=T))}))

#replacing version string
script_files <- c(R_files, r_files, sh_files, py_files)
build_files <- c(R_builds, r_builds, sh_builds, py_builds)
#print(script_files)
#print(build_files)
for( i in seq_along(script_files) ){
  x <- readLines(script_files[i])
  y <- gsub( "REPLACE_WITH_COMMIT_OR_VERSION", version, x )
  cat(y, file=build_files[i], sep="\n")
  Sys.chmod(build_files[i],mode="0755")
}
build_files <- c(build_files,test_scripts)

bin_dest=Sys.getenv('SWAN_BIN')
if(bin_dest=="") bin_dest <- dirname(file.path(Sys.getenv("R_LIBS_USER"),"swan","bin"))
dir.create(bin_dest, recursive = TRUE, showWarnings = FALSE)
n=length(build_files)
message(paste(rep("Installing",n), build_files, rep("to",n), rep(bin_dest,n), collapse="\n"))
file.copy(build_files, bin_dest, overwrite = TRUE)

cat(paste("library search paths:",.libPaths(),"\n"))
#source("http://bioconductor.org/biocLite.R")
#biocLite(pkgs=c("Biobase", "Biostrings", "BSgenome", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools"))
								#suppressUpdates=T,suppressAutoUpdate=T,ask=F)
#install.packages(pkgs=c("BH", "Cairo", "RcppArmadillo", "Rcpp", "data.table", "digest", "ggplot2", "gridExtra", "hash", "methods", "optparse", "parallel", "plyr", "robustbase", "seqCBS", "sets", "stringr", "zoo"))
								#suppressUpdates=T,suppressAutoUpdate=T,ask=F)

cat(sprintf("#####################################\n"))
cat(sprintf("#                                   \n"))
cat(sprintf("#  Your SWAN Binaries can be found at:   \n"))
cat(sprintf("#  %s   \n", bin_dest))
cat(sprintf("#  To use SWAN, set environment variable $SWAN_BIN to above path \n"))
cat(sprintf("#  And add $SWAN_BIN to your $PATH evironment \n"))
cat(sprintf("#                                   \n"))
cat(sprintf("#####################################\n"))

#quit()  #can be used to debug install process
