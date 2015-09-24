library(inline)
library(Rcpp)
library(RcppArmadillo)

code <- '
  arma::mat coeff = Rcpp::as<arma::mat>(a);
  arma::mat errors = Rcpp::as<arma::mat>(e);
  std::cout<<"a.w="<<coeff.n_rows<<"a.h="<<coeff.n_cols<<std::endl;
  std::cout<<"e.w="<<errors.n_rows<<"e.h="<<errors.n_cols<<std::endl;
  int m = errors.n_rows; int n = errors.n_cols;
  arma::mat simdata(m,n);
  simdata.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    simdata.row(row) = simdata.row(row-1)*trans(coeff)+errors.row(row);
  }
  return Rcpp::wrap(simdata);
'
rcppSim <- cxxfunction(signature(a="numeric",e="numeric"),code,plugin="RcppArmadillo")

rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
}

a <- matrix(c(0.5,0.1,0.1,0.5),nrow=2)
e <- matrix(rnorm(10000),ncol=2)
rcppData <- rcppSim(a,e)
rData <- rSim(a, e)
suppressMessages(require(compiler))
compRsim <- cmpfun(rSim)
compRData <- compRsim(a,e)
stopifnot(all.equal(rData, compRData))
stopifnot(all.equal(rData, rcppData))
suppressMessages(library(rbenchmark))
res <- benchmark(rcppSim(a,e),
                 rSim(a,e),
                 compRsim(a,e),
                 columns=c("test", "replications", "elapsed",
                           "relative", "user.self", "sys.self"),
                 order="relative")

print(res)