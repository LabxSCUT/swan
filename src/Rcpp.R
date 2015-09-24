#!/usr/bin/env Rscript

library(inline)
library(rbenchmark)

## openMPCode example from Rcpp/examples/OpenMP/ by Dirk E.

openMPCode <- '
// assign to C++ vector
std::vector<double> x = Rcpp::as<std::vector< double > >(xs);
size_t n = x.size();
#pragma omp parallel for shared(x, n)
for (size_t i=0; i<n; i++) {
x[i] = ::log(x[i]);
}
return Rcpp::wrap(x);
'

# A few other examples

inc <- '
#include <iostream>
using namespace std;
#include <omp.h>
'
helloWorld <- '
int th_id, nthreads;
#pragma omp parallel private(th_id) shared(nthreads)
{
th_id = omp_get_thread_num();
#pragma omp critical
{
cout << "Hello World from thread " << th_id << "\\n";
}
#pragma omp barrier
#pragma omp master
{
nthreads = omp_get_num_threads();
cout << "There are " << nthreads << " threads" << "\\n";
}
}
return Rcpp::wrap(0);
'
exampleOpenMP <- '
// assign to C++ vector
std::vector<double> x = Rcpp::as<std::vector< double > >(xs);
size_t n = x.size();
#pragma omp parallel for shared(x, n)
for (size_t i=0; i<n; i++) {
for (int j=0; j< 1000000; j++) {
x[i] += 1;
}
}
return Rcpp::wrap(x);
'
exampleNonMP <- '
// assign to C++ vector
std::vector<double> x = Rcpp::as<std::vector< double > >(xs);
size_t n = x.size();
for (size_t i=0; i<n; i++) {
for (int j=0; j< 1000000; j++) {
x[i] += 1;
}
}
return Rcpp::wrap(x);
'
 
## modify the plugin for Rcpp to support OpenMP
settings <- getPlugin("Rcpp")
settings$env$PKG_CXXFLAGS <- paste('-fopenmp', settings$env$PKG_CXXFLAGS)
settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)
 
funOpenMP <- cxxfunction(signature(xs="numeric"), body=openMPCode, plugin="Rcpp", settings=settings)
helloWorld <- cxxfunction(inc=inc, body=helloWorld, plugin="Rcpp", settings=settings)
exampleOpenMP <- cxxfunction(signature(xs="numeric"),inc=inc, body=exampleOpenMP, plugin="Rcpp", settings=settings)
exampleNonMP <- cxxfunction(signature(xs="numeric"),inc=inc, body=exampleNonMP, plugin="Rcpp", settings=settings)
 
funOpenMP(x=1:5)
helloWorld()
system.time(exampleNonMP(x=1:1000))
system.time(exampleOpenMP(x=1:1000)) # roughly half the time with 2 cores

res <- benchmark( exampleNonMP(x=1:1000),
                  exampleOpenMP(x=1:1000),
                  columns=c("test", "replications", "elapsed",
                            "relative", "user.self", "sys.self"),
                  order="relative",
                  replications=100)
print(res)