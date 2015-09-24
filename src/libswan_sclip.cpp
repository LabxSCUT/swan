#include "libswan_sclip.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "samtools/sam.h"
#include "bwa/bwamem.h"

#define NA -1  // TODO: fix

using namespace Rcpp;
using namespace std;

typedef struct {
    int *pos;
    int icnt, ncnt;
    char **seq;
    int *nLclipped, *nfromLend, *nmatch, *nfromRpos, *nRclipped;
} SCAN_BAM_DATA;

typedef struct {
  string consensus;
  vector<vector<int> > mat;
  vector<int> nseqs;
  double propMisMatch;
} ConsMatrix;

typedef struct {
  int position;
  vector<int> indices;
  int nreads;
  int nbases;
  vector<string> clippedseqs;
  ConsMatrix consmat;
} Cluster;

static void grow_sbd(SCAN_BAM_DATA* sbd) {
  int len = sbd->ncnt;
  sbd->pos = (int*)malloc(sizeof(int) * len);
  sbd->seq = (char**)malloc(sizeof(char*) * len);
  sbd->nLclipped = (int*)malloc(sizeof(int) * len);
  sbd->nfromLend = (int*)malloc(sizeof(int) * len);
  sbd->nmatch = (int*)malloc(sizeof(int) * len);
  sbd->nfromRpos = (int*)malloc(sizeof(int) * len);
  sbd->nRclipped = (int*)malloc(sizeof(int) * len);
}

static void destroy_sbd(SCAN_BAM_DATA* sbd) {
  int len = sbd->ncnt;
  free(sbd->pos);
  for (int i = 0; i < len; ++i)
    free(sbd->seq[i]);
  free(sbd->seq);
  free(sbd->nLclipped);
  free(sbd->nfromLend);
  free(sbd->nmatch);
  free(sbd->nfromRpos);
  free(sbd->nRclipped);
}

static int count_func(const bam1_t* bam, void* data) {
  uint32_t flag = bam->core.flag;
  if ((flag & BAM_FQCFAIL) || (flag & BAM_FDUP))
    return 0;
  if (flag & BAM_FUNMAP)
    return 0;

  // only count reads with 'S' in the cigar
  uint32_t n_cigar = bam->core.n_cigar;
  const uint32_t* cigar = bam1_cigar(bam);
  for (int i = 0; i < n_cigar; ++i) {
    if ((cigar[i] & BAM_CIGAR_MASK) == 4) { // 4 is 'S'
      SCAN_BAM_DATA* sbd = (SCAN_BAM_DATA*)data;
      sbd->ncnt++;
      return 0;
    }
  }

  return 0;
}

static char *_bamseq(const bam1_t * bam) {
  static const char key[] = {
    '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
  };

  const uint32_t len = bam->core.l_qseq;
  const unsigned char* seq = bam1_seq(bam);
  char* s = (char*)malloc(len + 1);
  for (uint32_t i = 0; i < len; ++i)
    s[i] = key[bam1_seqi(seq, i)];
  s[len] = '\0';
  return s;
}

static void sbd_assign(SCAN_BAM_DATA* sbd, int idx, int nLclipped, int nfromLend, int nmatch, int nfromRpos, int nRclipped) {
  sbd->nLclipped[idx] = nLclipped;
  sbd->nfromLend[idx] = nfromLend;
  sbd->nmatch[idx] = nmatch;
  sbd->nfromRpos[idx] = nfromRpos;
  sbd->nRclipped[idx] = nRclipped;
}

static bool process_cigar(const bam1_t* bam, SCAN_BAM_DATA* sbd, int idx) {
  uint32_t n_cigar = bam->core.n_cigar;
  const uint32_t* cigar = bam1_cigar(bam);
  int nLclipped, nfromLend, nmatch = 0, nfromRpos, nRclipped;

  char ops[n_cigar];
  uint32_t opns[n_cigar];

  const char lookup[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P'};
  unsigned int ns = 0;  // number of S's
  unsigned int si[2];  // indices of S's
  int mi0 = -1;  // index of the first M

  for (int i = 0; i < n_cigar; ++i) {
    char c = lookup[cigar[i] & BAM_CIGAR_MASK];
    ops[i] = c;
    opns[i] = cigar[i] >> 4;
    if (c == 'S') {
      if (ns == 2) { fprintf(stderr, "More than 2 'S' in CIGAR string!\n"); exit(1); }
      si[ns++] = i;
    } else if (c == 'M') {
      nmatch += opns[i];
      if (mi0 == -1)
        mi0 = i;
    }
  }

  if (ns == 0)
    return false;

  if (nmatch == 0) {
    sbd_assign(sbd, idx, 0, NA, 0, NA, 0);
    return true;
  }

  if (ns == 1) {
    if (ops[0] == 'S' || (ops[1] == 'S' && ops[0] == 'H')) {
      // Left clipped, no right clip.
      nLclipped = opns[si[0]];
      nRclipped = 0;
      nfromLend = 0;
      if (si[0] + 1 < mi0) { // TODO: why do we have to check this? Seems like a bug in the orig. code
        for (int j = si[0] + 1; j <= mi0; ++j) {
          nfromLend += opns[j];
        }
      }
      sbd_assign(sbd, idx, nLclipped, nfromLend, nmatch, NA, nRclipped);
      return true;
    } else {
      // No left clip, right clipped
      nLclipped = 0;
      nRclipped = opns[si[0]];
      nfromRpos = nmatch;
      for (int j = 0; j < n_cigar; ++j) {
        if (ops[j] == 'D')
          nfromRpos += opns[j];
      }
      sbd_assign(sbd, idx, nLclipped, NA, nmatch, nfromRpos, nRclipped);
      return true;
    }
  }

  // n == 2
  {
    // There were two "S", one for left and one for right.
    nLclipped = opns[si[0]];
    nRclipped = opns[si[1]];
    nfromLend = 0;
    if (si[0] + 1 < mi0) { // TODO: why do we have to check this? Seems like a bug in the orig. code
      for (int j = si[0] + 1; j <= mi0; ++j) {
        nfromLend += opns[j];
      }
    }
    nfromRpos = nmatch;
    for (int j = 0; j < n_cigar; ++j) {
      if (ops[j] == 'D')
        nfromRpos += opns[j];
    }
    sbd_assign(sbd, idx, nLclipped, nfromLend, nmatch, nfromRpos, nRclipped);
    return true;
  }
}

static int fetch_func(const bam1_t* bam, void* data) {
  uint32_t flag = bam->core.flag;
  if ((flag & BAM_FQCFAIL) || (flag & BAM_FDUP))
    return 0;
  if (flag & BAM_FUNMAP)
    return 0;

  SCAN_BAM_DATA* sbd = (SCAN_BAM_DATA*)data;
  int idx = sbd->icnt;

  if (!process_cigar(bam, sbd, idx))
    return 0;
  sbd->pos[idx] = bam->core.pos + 1;
  sbd->seq[idx] = _bamseq(bam);

  sbd->icnt++;
  return 0;
}

static double proportionMismatch(const vector<vector<int> >& mat, int MIN_OBS = 2) {
  int n = mat[0].size();
  int nobs[n];
  int maxcount[n];
  memset(nobs, 0, sizeof(int) * n);
  memset(maxcount, 0, sizeof(int) * n);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < n; ++j) {
      nobs[j] += mat[i][j];
      maxcount[j] = max(maxcount[j], mat[i][j]);
    }
  }
  int sum_nmismatch = 0;
  int sum_nobs = 0;
  for (int i = 0; i < n; ++i) {
    if (nobs[i] >= MIN_OBS) {
      sum_nmismatch += nobs[i] - maxcount[i];
      sum_nobs += nobs[i];
    }
  }
  return sum_nmismatch / (double)sum_nobs;
}

static int base2num(char c) {
  if (c == 'A') return 0;
  if (c == 'C') return 1;
  if (c == 'G') return 2;
  if (c == 'T') return 3;
  return 4;
}

static ConsMatrix consensusMatrix(const vector<string>& ss, bool startFromRight) {
  size_t max_nchar = 0;
  for (size_t i = 0; i < ss.size(); ++i)
    max_nchar = max(max_nchar, ss[i].length());

  // construct the matrix
  vector<vector<int> > mat(5);
  for (int i = 0; i < 5; ++i)
    mat[i].resize(max_nchar, 0);
  for (size_t i = 0; i < ss.size(); ++i) {
    string s = ss[i];
    if (startFromRight)
      s = string(s.rbegin(), s.rend());
    for (size_t j = 0; j < s.length(); ++j)
      mat[base2num(s[j])][j]++;
  }

  // find consensus
  string consensus(max_nchar, ' ');
  int consensus_c[max_nchar];
  memset(consensus_c, 0, sizeof(int) * max_nchar);
  static char letters[] = {'A', 'C', 'G', 'T', 'N'};
  for (int i = 0; i < 5; ++i) {
    for (size_t j = 0; j < max_nchar; ++j) {
      if (mat[i][j] > consensus_c[j]) {
        consensus[j] = letters[i];
        consensus_c[j] = mat[i][j];
      }
    }
  }
  if (startFromRight)
    consensus = string(consensus.rbegin(), consensus.rend());

  int nseqs[max_nchar];
  memset(nseqs, 0, sizeof(int) * max_nchar);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < max_nchar; ++j) {
      nseqs[j] += mat[i][j];
    }
  }

  double propMisMatch = proportionMismatch(mat);

  ConsMatrix ret;
  ret.consensus = consensus;
  ret.mat = mat;
  ret.nseqs.assign(nseqs, nseqs + max_nchar);
  ret.propMisMatch = propMisMatch;
  return ret;
}

//static int distanceFromGap(int pos, const string& seqname, const List& gap) {
//  const CharacterVector& chrom = gap[0];
//  const NumericVector& chromStart = gap[1];
//  const NumericVector& chromEnd = gap[2];
//  for (int i = 0; i < chrom.length(); ++i) {
//    if (strcmp(chrom[i].begin(), seqname.c_str()) == 0) {
//      if (pos > chromStart[i] && pos < chromEnd[i])
//        return 0;
//      return min(abs(chromStart[i] - pos), abs(chromEnd[i] - pos));
//    }
//  }
//  return INT_MAX;
//}

// TODO: too much copying
static vector<vector<Cluster> > cpp_core_sclip(int ti, int n_trunks, int scan_start, int scan_end,
  int trunk_size, const vector<string>& files, const string& seqname, const List& gap,
  int MIN_READS_PER_CLUSTER, int MIN_BASES_PER_CLUSTER, double SC_PROPORTION_MISMATCH_THRESH,
  int MIN_GAP_DISTANCE) {

  int st0 = scan_start + (ti-1) * trunk_size;
  int ed0 = min(scan_start + ti * trunk_size - 1, scan_end);

  //printf("----- in trunk %d of %d from %d to %d\n", ti, n_trunks, st0, ed0);

  //if (distanceFromGap(st0, seqname, gap) < MIN_GAP_DISTANCE ||
  //  distanceFromGap(ed0, seqname, gap) < MIN_GAP_DISTANCE) {
  //   printf("skipping gap region\n");
  //   vector<vector<Cluster> > ret(2);
  //   return ret;
  //}

  SCAN_BAM_DATA sbd;
  sbd.ncnt = 0;
  sbd.icnt = 0;

  // count the number of reads from all bam files
  for (size_t i = 0; i < files.size(); ++i) {
    const string& bam_file = files[i];
    samfile_t* bf = samopen(bam_file.c_str(), "rb", 0);
    bam_index_t* bam_index = bam_index_load(bam_file.c_str());
    int b_ref, b_start, b_end;
    bam_parse_region(bf->header, seqname.c_str(), &b_ref, &b_start, &b_end);
    bam_fetch(bf->x.bam, bam_index, b_ref, st0 - 1, ed0 - 1, &sbd, count_func);
    bam_index_destroy(bam_index);
    samclose(bf);
  }

  grow_sbd(&sbd);

  // then read the data from the bam files
  for (size_t i = 0; i < files.size(); ++i) {
    const string& bam_file = files[i];
    samfile_t* bf = samopen(bam_file.c_str(), "rb", 0);
    bam_index_t* bam_index = bam_index_load(bam_file.c_str());
    int b_ref, b_start, b_end;
    bam_parse_region(bf->header, seqname.c_str(), &b_ref, &b_start, &b_end);
    bam_fetch(bf->x.bam, bam_index, b_ref, st0 - 1, ed0 - 1, &sbd, fetch_func);
    bam_index_destroy(bam_index);
    samclose(bf);
  }

  int np = sbd.ncnt; // TODO: better name
  //printf("%d reads have soft clipping.\n", np);

  int minpos = st0;
  int maxpos = ed0;

  int n = maxpos - minpos + 1;

  // For each base in region, the following keeps the number of reads clipped there, and total bases clipped.
  int* nclipreadsL = (int*)calloc(n, sizeof(int));
  int* nclipbasesL = (int*)calloc(n, sizeof(int));
  int* nclipreadsR = (int*)calloc(n, sizeof(int));
  int* nclipbasesR = (int*)calloc(n, sizeof(int));

  // For each entry in cigar, this keeps the exact reference base position of the clip.
  int* clipposL = (int*)malloc(sizeof(int) * np); memset(clipposL, -1, sizeof(int) * np); // TODO: hacky
  int* clipposR = (int*)malloc(sizeof(int) * np); memset(clipposR, -1, sizeof(int) * np); // TODO: hacky

  // Get the base counts for the left clips.
  // Position clipped is one base to the right of clipped base.
  // i.e. xxxmmmmm , position is position of first m.
  for (int i = 0; i < np; ++i) {
    if (sbd.nLclipped[i] > 0) {
      int bcoord = sbd.pos[i];
      clipposL[i] = bcoord;
      int x = bcoord - minpos;
      if (x >= 0 && x < n){
        nclipreadsL[x]++;
        nclipbasesL[x] += sbd.nLclipped[i];
      }
    }
  }

  // Get the base counts for the right clips.
  // Position clipped is defined as the position of the first clipped base.
  // i.e. mmmmmxxx , position is position of first x.
  for (int i = 0; i < np; ++i) {
    if (sbd.nRclipped[i] > 0) {
      int bcoord = sbd.pos[i] + sbd.nfromRpos[i]; // nfromRpos -> ntoRstart
      clipposR[i] = bcoord;
      int x = bcoord - minpos;
      if (x >= 0 && x < n) {
        nclipreadsR[x]++;
        nclipbasesR[x] += sbd.nRclipped[i];
      }
    }
  }

  // Threshold to get L and R clusters
  vector<Cluster> LsoftclipClusters;
  vector<Cluster> RsoftclipClusters;

  for (int x = 0; x < n; ++x) {
    if (nclipreadsL[x] >= MIN_READS_PER_CLUSTER && nclipbasesL[x] >= MIN_BASES_PER_CLUSTER) {
      int bcoord = x + minpos;
      int nbases = nclipbasesL[x];
      vector<string> cseqs;
      vector<int> indices;
      for (int i = 0; i < np; ++i) {
        if (clipposL[i] == bcoord) {
          cseqs.push_back(string(sbd.seq[i]).substr(0, sbd.nLclipped[i]));
          indices.push_back(i + 1);
        }
      }
      const ConsMatrix& consmat = consensusMatrix(cseqs, true);
      if (consmat.propMisMatch >= SC_PROPORTION_MISMATCH_THRESH)
        continue;
      Cluster cluster;
      cluster.position = bcoord;
      cluster.indices = indices;
      cluster.nreads = indices.size();
      cluster.nbases = nbases;
      cluster.clippedseqs = cseqs;
      cluster.consmat = consmat;
      LsoftclipClusters.push_back(cluster);
    }
  }

  for (int x = 0; x < n; ++x) {
    if (nclipreadsR[x] >= MIN_READS_PER_CLUSTER && nclipbasesR[x] >= MIN_BASES_PER_CLUSTER) {
      int bcoord = x + minpos;
      int nbases = nclipbasesR[x];
      vector<string> cseqs;
      vector<int> indices;
      for (int i = 0; i < np; ++i) {
        if (clipposR[i] == bcoord) {
          cseqs.push_back(string(sbd.seq[i]).substr( strlen(sbd.seq[i]) - sbd.nRclipped[i], sbd.nRclipped[i] ));
          indices.push_back(i + 1);
        }
      }
      const ConsMatrix& consmat = consensusMatrix(cseqs, false);
      if (consmat.propMisMatch >= SC_PROPORTION_MISMATCH_THRESH)
        continue;
      Cluster cluster;
      cluster.position = bcoord;
      cluster.indices = indices;
      cluster.nreads = indices.size();
      cluster.nbases = nbases;
      cluster.clippedseqs = cseqs;
      cluster.consmat = consmat;
      RsoftclipClusters.push_back(cluster);
    }
  }

  free(nclipreadsL);
  free(nclipbasesL);
  free(nclipreadsR);
  free(nclipbasesR);

  free(clipposL);
  free(clipposR);

  destroy_sbd(&sbd);

  //printf("----- Out of trunk %d of %d\n", ti, n_trunks);
  //printf("\nLeft-clipped: %lu, right-clipped: %lu\n", LsoftclipClusters.size(), RsoftclipClusters.size());

  vector<vector<Cluster> > ret;
  ret.push_back(LsoftclipClusters);
  ret.push_back(RsoftclipClusters);
  return ret;
}

static List consmat_obj(const ConsMatrix& consmat) {
  //const vector<vector<int> >& mat = consmat.mat;
  //int n = mat[0].size();
  //List rmat(n);
  //for (int i = 0; i < n; ++i) {
  //  rmat[i] = NumericVector::create(mat[0][i], mat[1][i], mat[2][i], mat[3][i], mat[4][i]);
  //}

  //StringVector col_names(n);
  //for (int i = 0; i < n; ++i) {
  //  char name[10];
  //  sprintf(name, "V%d", i+1);
  //  col_names(i) = name;
  //}

  //static const string row_names[] = {"A", "C", "G", "T", "N"};
  //rmat.attr("names") = col_names;
  //rmat.attr("row.names") = StringVector(row_names, row_names + 5);
  //rmat.attr("class") = "data.frame";

  return List::create(
    _["consensus"] = consmat.consensus,
    //_["mat"] = rmat,
    //_["nseqs"] = NumericVector(consmat.nseqs.begin(), consmat.nseqs.end()),
    _["propMisMatch"] = consmat.propMisMatch
  );
}

static List cluster_obj(const Cluster& cluster) {
  return List::create(
    _["position"] = (double)cluster.position,
    //_["indices"] = cluster.indices,
    _["nreads"] = cluster.nreads,
    _["nbases"] = cluster.nbases,
    //_["clippedseqs"] = cluster.clippedseqs,
    _["consmat"] = consmat_obj(cluster.consmat)
  );
}

static List clusters_obj(const vector<Cluster>& clusters) {
  List list(clusters.size());
  for (size_t i = 0; i < clusters.size(); ++i)
    list[i] = cluster_obj(clusters[i]);
  return list;
}

SEXP core_sclip(SEXP ti, SEXP n_trunks, SEXP scan_start, SEXP scan_end, SEXP trunk_size, SEXP files, SEXP seqname, SEXP gap,
  SEXP MIN_READS_PER_CLUSTER, SEXP MIN_BASES_PER_CLUSTER, SEXP SC_PROPORTION_MISMATCH_THRESH, SEXP MIN_GAP_DISTANCE) {

  const vector<vector<Cluster> >& res = cpp_core_sclip(
    Rcpp::as<int>(ti),
    Rcpp::as<int>(n_trunks),
    Rcpp::as<int>(scan_start),
    Rcpp::as<int>(scan_end),
    Rcpp::as<int>(trunk_size),
    Rcpp::as<vector<string> >(files),
    Rcpp::as<string>(seqname),
    Rcpp::as<List>(gap),
    Rcpp::as<int>(MIN_READS_PER_CLUSTER),
    Rcpp::as<int>(MIN_BASES_PER_CLUSTER),
    Rcpp::as<double>(SC_PROPORTION_MISMATCH_THRESH),
    Rcpp::as<int>(MIN_GAP_DISTANCE)
  );

  return List::create(_["scL"] = clusters_obj(res[0]), _["scR"] = clusters_obj(res[1]));
}

const char* bwa_pg = ""; // to silence the compiler

typedef struct {
  bwaidx_t* idx;
  mem_opt_t* opt;
} BWA_DATA;

static void bwa_data_finalizer(SEXP ptr) {
  BWA_DATA* data = (BWA_DATA*)R_ExternalPtrAddr(ptr);
  if (!data)
    return;
  free(data->opt);
  bwa_idx_destroy(data->idx);
  free(data);
  R_ClearExternalPtr(ptr);
}

SEXP scanFa(SEXP file) {
  BWA_DATA* data = (BWA_DATA*)malloc(sizeof(BWA_DATA));

  data->idx = bwa_idx_load(Rcpp::as<string>(file).c_str(), BWA_IDX_ALL);
  if (data->idx == NULL) {
    fprintf(stderr, "Index load failed. Please use \"bwa index <ref_file>\".\n");
    exit(1);
  }

  data->opt = mem_opt_init();

  vector<int> wrapper_data;
  vector<string> wrapper_names;
  for (int i = 0; i < data->idx->bns->n_seqs; ++i) {
    wrapper_data.push_back(0);
    wrapper_names.push_back(data->idx->bns->anns[i].name);
  }

  IntegerVector wrapper(wrapper_data.begin(), wrapper_data.end());
  wrapper.attr("class") = string("BWA");
  wrapper.attr("names") = wrapper_names;

  SEXP ptr = PROTECT(R_MakeExternalPtr(data, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, bwa_data_finalizer, TRUE);
  wrapper.attr("ptr") = ptr;
  UNPROTECT(1);

  return wrapper;
}

SEXP matchPattern(SEXP pattern, SEXP rdata) {
  const int MAX_MATCHES = 10;  // FIXME: Make this optimization more explicit.

  string seq = Rcpp::as<string>(pattern);
  SEXP ptr = Rcpp::as<IntegerVector>(rdata).attr("ptr");
  BWA_DATA* data = (BWA_DATA*)R_ExternalPtrAddr(ptr);
  bwaidx_t* idx = data->idx;
  mem_opt_t* opt = data->opt;
  opt->a = 1; // match score
  opt->b = 4; // mismatch penalty
  opt->o_del = opt->o_ins = 6; // gap open penalties

  vector<vector<vector<int> > > ret(idx->bns->n_seqs);
  for (int i = 0; i < idx->bns->n_seqs; ++i)
    ret[i].resize(2); // start & width arrays

  int n_matches = 0;
  set<int64_t> starts;

  mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, seq.length(), seq.c_str());
  for (int i = 0; i < ar.n; ++i) {
    // Allow at most one mismatch. Note: the check below is not sufficient.
    if (ar.a[i].truesc != (int)seq.length() && ar.a[i].truesc != (int)seq.length() - 1 - 4) continue;
    mem_aln_t a = mem_reg2aln(opt, idx->bns, idx->pac, seq.length(), seq.c_str(), &ar.a[i]);
    if (!a.is_rev && a.n_cigar == 1 && (a.cigar[0]&0xf) == 0) { // check if the CIGAR consists only of M's
      n_matches++;
      ret[ar.a[i].rid][0].push_back(a.pos + 1);
      ret[ar.a[i].rid][1].push_back(seq.length());
      starts.insert(ar.a[i].rb);
    }
    free(a.cigar);
    if (n_matches > MAX_MATCHES)
      break;
  }
  free(ar.a);

  if (n_matches <= MAX_MATCHES) {
    // Mutate each base and find alignments for each mutation.
    // The reason we do this is that BWA-MEM doesn't return all the hits within the criteria,
    // especially when the query sequence is highly repetitive.
    // Doing this allows us to retrieve the hits that we have missed in the first round.
    static const char bases[] = "ACGT";
    for (int k = 0; k < seq.length(); ++k) {
      // TODO: check for 'N'
      string seq2 = seq;
      for (int j = 0; j < 4; ++j) {
        char base = bases[j];
        if (base == seq[k]) continue;
        seq2[k] = base;

        mem_alnreg_v ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, seq2.length(), seq2.c_str());
        for (int i = 0; i < ar.n; ++i) {
          if (ar.a[i].truesc != (int)seq.length()) continue; // Only accept exact matches.
          mem_aln_t a = mem_reg2aln(opt, idx->bns, idx->pac, seq2.length(), seq2.c_str(), &ar.a[i]);
          if (!a.is_rev) {
            if (starts.find(ar.a[i].rb) == starts.end()) {
              n_matches++;
              ret[ar.a[i].rid][0].push_back(a.pos + 1);
              ret[ar.a[i].rid][1].push_back(seq.length());
              starts.insert(ar.a[i].rb);
            }
          }
          free(a.cigar);
        }
        free(ar.a);
      }
      if (n_matches > MAX_MATCHES)
        break;
    }
  }

  // Sort IRanges
  // TODO: Is this necessary?
  for (int i = 0; i < ret.size(); ++i) {
    if (ret[i][0].size() > 1)
      sort(ret[i][0].begin(), ret[i][0].end());
  }

  return Rcpp::wrap(ret);
}

SEXP getSubseq(SEXP rdata, SEXP rid1, SEXP start, SEXP end) {
  SEXP ptr = Rcpp::as<IntegerVector>(rdata).attr("ptr");
  BWA_DATA* data = (BWA_DATA*)R_ExternalPtrAddr(ptr);
  bwaidx_t* idx = data->idx;

  int64_t rlen;
  int rid = Rcpp::as<int>(rid1) - 1;
  int64_t rb = idx->bns->anns[rid].offset + Rcpp::as<int>(start) - 1;
  int64_t re = idx->bns->anns[rid].offset + Rcpp::as<int>(end) - 1 + 1;
  uint8_t *rseq = bns_get_seq(idx->bns->l_pac, idx->pac, rb, re, &rlen);
  string ret(rlen, ' ');
  for (int64_t i = 0; i < rlen; ++i)
    ret[i] = "ACGT"[rseq[i]];
  free(rseq);

  // Add holes. We need to do this because BWA uses 2-bit encoded reference
  // and replaces N's with random bases.
  bntamb1_t* ambs = idx->bns->ambs;
  int l = 0;
  int r = idx->bns->n_holes - 1;
  while (l <= r) {
    int mid = (l + r) >> 1;
    if (ambs[mid].offset > rb) {
      r = mid - 1;
    } else if (ambs[mid].offset < rb) {
      l = mid + 1;
    } else {
      r = mid;
      break;
    }
  }
  if (r < 0)
    r = 0;
  for (int i = r; i < idx->bns->n_holes; ++i) {
    if (ambs[i].offset >= re) break;
    if (ambs[i].offset + ambs[i].len > rb) {
      int j0 = max(ambs[i].offset - rb, (int64_t)0);
      int j1 = min(ambs[i].offset + ambs[i].len - rb, rlen);
      for (int j = j0; j < j1; ++j)
        ret[j] = ambs[i].amb;
    }
  }

  return Rcpp::wrap(ret);
}
