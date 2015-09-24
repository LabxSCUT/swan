#ifndef _LIBSWAN_SCLIP
#define _LIBSWAN_SCLIP

#include <RcppArmadillo.h>

//[[Rcpp::export]]
RcppExport SEXP core_sclip(SEXP ti, SEXP n_trunks, SEXP scan_start, SEXP scan_end, SEXP trunk_size, SEXP files, SEXP seqname, SEXP gap, SEXP MIN_READS_PER_CLUSTER, SEXP MIN_BASES_PER_CLUSTER, SEXP SC_PROPORTION_MISMATCH_THRESH, SEXP MIN_GAP_DISTANCE);

RcppExport SEXP scanFa(SEXP file);
RcppExport SEXP matchPattern(SEXP pattern, SEXP rdata);
RcppExport SEXP getSubseq(SEXP rdata, SEXP rid1, SEXP start, SEXP end);

#endif
