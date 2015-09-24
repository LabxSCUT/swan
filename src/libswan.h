#ifndef _LIBSWAN
#define _LIBSWAN

#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif /* Platform Supports Openmp, see https://cran.r-project.org/doc/manuals/r-release/R-exts.html#OpenMP-support */

// [[Rcpp::depends(BH)]]
// [[RcppArmadillo::depends(BH)]]  
#include <RcppArmadillo.h>
//#include <armadillo>
//#include <Rcpp.h> //Cause Error: Rcpp.h shouldn't be included, include RcppArmadillo.h
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
#include <unordered_map>
#include <unordered_set>
//#include <boost/unordered_map.hpp>
//#include <boost/math/common_factor.hpp>
#include <boost/unordered_set.hpp>
#include <boost/regex.hpp>
//#include <regex>  //Cause Error: smatch[] not supported currently by GCC, use boost/regex
//#include <boost/icl/interval_map.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>

//This is the right way to use Boost with RcppArmadillo
//This is the right way to use Boost with RcppArmadillo

using namespace std;
using namespace Rcpp;
RNGScope scope;
typedef unsigned long ul;
typedef unsigned long long ull;
typedef pair<ull,ull> rd;              //rd = <st,ed>
typedef pair<pair<rd,rd>,int> rdp;     //rdp = <<rd1, rd2>, count>
typedef pair<ull,pair<ull,bool> > evt; //event = <pos,<rdp_idx,in/out>>
typedef vector<ull> clust;             //clust = {rdp1,rpd2,...}

namespace std{ namespace tr1{
    template<typename T>
    struct hash< std::pair<T, T> >
    {
        typedef std::pair<T, T> argument_type;
        typedef size_t result_type;
        result_type operator()(const argument_type& a) const
        {
            hash<int> hasher;
            result_type h = 0;
            h = hasher(a.first) * 31 + hasher(a.second);
            return h;
        };
    };
} };

//typedef std::tr1::unordered_set<std::pair<int,int> > hashset;
//typedef std::tr1::unordered_map<string, int> hashmap;
typedef std::unordered_map<string, int> hashmap;
typedef boost::unordered_set<std::pair<int,int> > hashset; //use std::unordered_set cause error
typedef boost::regex regex;
typedef boost::smatch smatch;
//typedef std::regex regex; //TODO: until GCC fully implements regex
//typedef std::smatch smatch;
//typedef boost::unordered_map<string, int> hashmap;
//typedef boost::icl::interval_map interval_map;
//typedef boost::icl::interval interval;
typedef std::pair<int, int> key;
typedef std::map<key, int> pairmap;

class PairMap{
  public:
    PairMap(){
      m = pairmap();
    }
    void add(int x, int y, int s){
      key k = std::make_pair<int&, int&>(x, y);
      pit = m.find(k);  
      if(pit != m.end()){ pit->second += s; }
      else{ m[k] = s; }
    }
    SEXP read(){
      Rcpp::List ret; int C=m.size(); 
      arma::uvec X_uvec=arma::zeros<arma::uvec>(C);
      arma::uvec Y_uvec=arma::zeros<arma::uvec>(C);
      arma::uvec N_uvec=arma::zeros<arma::uvec>(C);
      int i=0;
      for(pit=m.begin();pit!=m.end();pit++){
        X_uvec(i)=(pit->first).first;
        Y_uvec(i)=(pit->first).second;
        N_uvec(i)=pit->second;
        i=i+1;
      }
      ret["x"]=Rcpp::wrap(X_uvec); 
      ret["y"]=Rcpp::wrap(Y_uvec); 
      ret["n"]=Rcpp::wrap(N_uvec);
      return ret;
    }
  private:
    pairmap m;
    pairmap::iterator pit;
};

RCPP_MODULE(pairmap_module) {

  class_<PairMap>( "PairMap")

  .constructor()

  .method( "add", &PairMap::add )
  .method( "read", &PairMap::read )
  ;

};

class HashMap{
  public:
    HashMap() {
      map = hashmap();
    }
    int set(string key, int val){ 
      hit = map.insert(map.begin(), pair<string, int>(key, val));
      return hit->second;
    }
    int get(string key) {
      hit = map.find(key);
      if(hit!=map.end()) { return hit->second; } else { return 0; };
    }
  private:
    hashmap map;
    hashmap::iterator hit;
};

RCPP_MODULE(hashmap_module) {

  class_<HashMap>( "HashMap")

    .constructor()

    .method( "set", &HashMap::set )
    .method( "get", &HashMap::get )
    ;
};

class HashSet{
  public:
    HashSet() {
      s = hashset();
    }
    int set(int x, int y){
      s.insert(pair<int, int>(x,y));
      return true; //true is 1 in R
    }
    int get(int x, int y) {
      return s.find(pair<int, int>(x,y)) != s.end(); //false is 0 in R
    }
  private:
    hashset s;
};

RCPP_MODULE(hashset_module) {

  class_<HashSet>( "HashSet")

  .constructor()

  .method( "set", &HashSet::set )
  .method( "get", &HashSet::get )
  ;
};

//[[Rcpp::export]]
RcppExport SEXP rcpp_hello_world();
RcppExport SEXP swan_unit_test();
RcppExport SEXP lC_omp(SEXP winC_rPb, SEXP rPb_isize, SEXP fy, SEXP n_win, SEXP mixing_rate);
RcppExport SEXP lD_omp(SEXP rS_winD, SEXP winD, SEXP rS, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p, SEXP left);
RcppExport SEXP lDp_omp(SEXP rPs_winDp, SEXP winW_start, SEXP rPs_start, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p);
RcppExport SEXP lDm_omp(SEXP rMs_winDm, SEXP winW_end, SEXP rMs_end, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p);
RcppExport SEXP mate_MPRs_omp(SEXP h1st_idx, SEXP c1st_qname, SEXP c1st_cigar, SEXP c1st_isize, 
                   SEXP c2nd_qname, SEXP c2nd_cigar, SEXP c2nd_isize, SEXP impute, SEXP left, SEXP self);
RcppExport SEXP get_al_omp(SEXP cigar, SEXP al_string);
RcppExport SEXP aggregate_sum_omp(SEXP coverage, SEXP start, SEXP end);
RcppExport SEXP anch_clust(SEXP Lfwd, SEXP Lchr, SEXP Lst, SEXP Led, SEXP Rfwd, SEXP Rchr, SEXP Rst, SEXP Red, SEXP chr, SEXP chr_size);
RcppExport SEXP contrast_clust(SEXP xLfwd, SEXP xLchr, SEXP xLst, SEXP xLed, SEXP xRfwd, SEXP xRchr, SEXP xRst, SEXP xRed,
                                SEXP yLfwd, SEXP yLchr, SEXP yLst, SEXP yLed, SEXP yRfwd, SEXP yRchr, SEXP yRst, SEXP yRed,
                                SEXP chr, SEXP chr_size, SEXP xSup, SEXP ySup);

/*
// create an external pointer to a Uniform object
RcppExport SEXP PairMap__new() {
  // convert inputs to appropriate C++ types
  //double min = as<double>(min_), max = as<double>(max_);
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<PairMap> ptr( new PairMap(), true );
  // return the external pointer to the R side
  return ptr;
}
*/

//RcppExport SEXP countOverlaps_omp(SEXP db_start, SEXP db_end, SEXP query_start, SEXP query_end);

#endif
