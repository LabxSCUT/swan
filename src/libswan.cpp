#include "libswan.h"

SEXP lC_omp(SEXP winC_rPb, SEXP rPb_isize, SEXP fy, SEXP n_win, SEXP mixing_rate){
    arma::umat winC_rPb_umat = Rcpp::as<arma::umat>(winC_rPb);
    arma::uvec rPb_isize_uvec = Rcpp::as<arma::uvec>(rPb_isize);
    arma::vec fy_vec = Rcpp::as<arma::vec>(fy); 
    size_t N = winC_rPb_umat.n_rows; size_t W = Rcpp::as<int>(n_win); double r = Rcpp::as<double>(mixing_rate);
    arma::vec lC_term_vec = arma::zeros<arma::vec>(W); //score
    arma::vec cC_term_vec = arma::zeros<arma::vec>(W); //count
#pragma omp parallel for shared(lC_term_vec, rPb_isize, fy_vec, winC_rPb_umat, r, N)
    for(size_t i=0; i<N; i++) { //need convert winC_rPb 1-based index to 0-based index
      double v = log((1-r)+r*fy_vec(rPb_isize_uvec(winC_rPb_umat(i,1)-1)));
#pragma omp atomic
      lC_term_vec(winC_rPb_umat(i,0)-1) += v; //score++
      cC_term_vec(winC_rPb_umat(i,0)-1) += 1; //count++
    }
		Rcpp::List ret; ret["lC_term"]=Rcpp::wrap(lC_term_vec); ret["cC_term"]=Rcpp::wrap(cC_term_vec);
    return ret;
};

SEXP lD_omp(SEXP rS_winD, SEXP winD, SEXP rS, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p, SEXP left){
  arma::umat rS_winD_umat = Rcpp::as<arma::umat>(rS_winD);
  arma::uvec winD_uvec = Rcpp::as<arma::uvec>(winD);
  arma::uvec rS_uvec = Rcpp::as<arma::uvec>(rS);
  arma::vec Fx_vec = Rcpp::as<arma::vec>(Fx);
  size_t N = rS_winD_umat.n_rows; size_t W = Rcpp::as<int>(n_win); int sign=(Rcpp::as<bool>(left)) ? 1 : -1;
  //cout<<"sign=-1 if rSl and winDr; sign=1 if rSr and winDl, current is"<<sign<<endl;
  double r = Rcpp::as<double>(mixing_rate); double pp = Rcpp::as<double>(p);
  arma::vec lD_term_vec = arma::zeros<arma::vec>(W); //score
  arma::vec cD_term_vec = arma::zeros<arma::vec>(W); //count
#pragma omp parallel for shared(lD_term_vec, Fx_vec, rS_winD_umat, winD_uvec, rS_uvec, r, pp, N, sign)
  for(size_t i=0; i<N; i++) { //every 1-based index need to be converted to 0-based index
    //if(sign*(winD_uvec(rS_winD_umat(i,1)-1)-rS_uvec(rS_winD_umat(i,0)-1))<0 | sign*(winD_uvec(rS_winD_umat(i,1)-1)-rS_uvec(rS_winD_umat(i,0)-1))>=Fx_vec.size()){
    //  cout<<"rS="<<rS_uvec(rS_winD_umat(i,0)-1)<<"winD="<<winD_uvec(rS_winD_umat(i,1)-1)<<endl; exit(-1);
    //}
    double v = log(1+r*(1-pp)/pp*Fx_vec(sign*(winD_uvec(rS_winD_umat(i,1)-1)-rS_uvec(rS_winD_umat(i,0)-1))));
    //if(rS_uvec(rS_winD_umat(i,0)-1)>=49000 && rS_uvec(rS_winD_umat(i,0)-1)<=50000)
    //cout<<"rS="<<rS_uvec(rS_winD_umat(i,0)-1)<<"winD="<<winD_uvec(rS_winD_umat(i,1)-1)<<"v="<<v<<endl;
    //if(rS_uvec(rS_winD_umat(i,0)-1)>=60000 && rS_uvec(rS_winD_umat(i,0)-1)<=61000)
    //cout<<"rS="<<rS_uvec(rS_winD_umat(i,0)-1)<<"winD="<<winD_uvec(rS_winD_umat(i,1)-1)<<"v="<<v<<endl;
    //cout<<winD_uvec(rS_winD_umat(i,1)-1)<<"\t"<<rS_uvec(rS_winD_umat(i,0)-1)<<"\t"<<sign*(winD_uvec(rS_winD_umat(i,1)-1)-rS_uvec(rS_winD_umat(i,0)-1))<<"\t"<<Fx_vec(sign*(winD_uvec(rS_winD_umat(i,1)-1)-rS_uvec(rS_winD_umat(i,0)-1)))<<endl;
#pragma omp atomic
    lD_term_vec(rS_winD_umat(i,1)-1) += v;     //score++
    cD_term_vec(rS_winD_umat(i,1)-1) += 1;     //count++
  }
	Rcpp::List ret; ret["lD_term"]=Rcpp::wrap(lD_term_vec); ret["cD_term"]=Rcpp::wrap(cD_term_vec);
  return ret;
}

SEXP aggregate_sum_omp(SEXP coverage, SEXP start, SEXP end){
  arma::uvec coverage_uvec = Rcpp::as<arma::uvec>(coverage);
  arma::uvec start_uvec = Rcpp::as<arma::uvec>(start);
  arma::uvec end_uvec = Rcpp::as<arma::uvec>(end);
  arma::uvec sum_uvec = arma::zeros<arma::uvec>(start_uvec.size());
#pragma omp parallel for shared(coverage_uvec, start_uvec, end_uvec, sum_uvec)
  for(size_t i=0; i<start_uvec.size(); i++) {
    for(size_t j=start_uvec(i)-1; j<=end_uvec(i)-1; j++) {
      sum_uvec(i)=sum_uvec(i)+coverage_uvec(j);
    }
  }
  return Rcpp::wrap(sum_uvec);
}

//SEXP countOverlaps_omp(SEXP db_start, SEXP db_end, SEXP query_start, SEXP query_end){ //db=read, query=window
//  arma::uvec db_start_uvec = Rcpp::as<arma::uvec>(db_start);
//  arma::uvec db_end_uvec = Rcpp::as<arma::uvec>(db_end);
//  arma::uvec query_start_uvec = Rcpp::as<arma::uvec>(query_start);
//  arma::uvec query_end_uvec = Rcpp::as<arma::uvec>(query_end);
//  interval_map<size_t, arma::uvec> db;
//  arma::uvec overlap_count_uvec = arma::zeros<arma::uvec>(query_start_uvec.size()); //windowed overlap count
//  for(int i=0; i<db_start_uvec.size(); i++) {
//    db += make_pair(interval<size_t>::closed(db_start_uvec(i), db_end_uvec(i)), i+1);
//  }
//  interval_map<size_t, int>::const_iterator it=db.end();
//  interval<size_t>::type query;
//  for(int j=0; j<query_start_uvec.size(); j++){ //grow overlap as we go
//    query=interval<size_t>::closed(query_start_uvec(j),query_end_uvec(j));
//    it=db.find(query); cout<<it->first<<it->second<<endl;
//    overlap_count_uvec(j) += 0; //(it->second).size();
//  }
//  return Rcpp::wrap(overlap_count_uvec);
//}

static string named_capture(const string& text, const regex& pattern, const string& name){
  smatch what; string::const_iterator start, end;
  start=text.begin(); end=text.end();
  //cout<<"text="<<text<<",pattern="<<pattern<<endl;
  if(regex_search(start, end, what, pattern)){
    //cout<<"matched:"<<what[name]<<endl;
    return what[name];
  }
  return string("0");
}

static int named_sum(const string& text, const regex& pattern, const string& name){
  smatch what; string::const_iterator start, end;
  start=text.begin(); end=text.end(); int sum=0;
  //cout<<"text="<<text<<",pattern="<<pattern<<endl;
  while(regex_search(start, end, what, pattern)){
    start = what[0].second; //cout<<"matched:"<<what[name]<<":"<<what<<endl; 
    sum=sum+atoi(string(what[name]).c_str()); //cout<<"sum="<<sum<<endl;
  }
  return sum;
}

SEXP get_al_omp(SEXP cigar, SEXP al_string){
  //regex alr_pattern("(?<size>[0-9]+)[MX=IN]"); //aligned length of read, including insertion
  //regex alg_pattern("(?<size>[0-9]+)[MX=DN]"); //aligned length of genome, including deletion
  string al_pattern_string=Rcpp::as<std::string>(al_string);
  regex al_pattern(al_pattern_string.c_str()); string name("size");
  vector<string> cigar_svec = Rcpp::as<vector<string> >(cigar);
  arma::uvec size_uvec = arma::zeros<arma::uvec>(cigar_svec.size());
#pragma omp parallel for shared(size_uvec, cigar_svec, al_pattern, name)
  for(size_t i=0; i<cigar_svec.size(); i++) { size_uvec[i] = named_sum(cigar_svec[i], al_pattern, name); }
  return Rcpp::wrap(size_uvec);
}

SEXP mate_MPRs_omp(SEXP h1st_idx, SEXP c1st_qname, SEXP c1st_cigar, SEXP c1st_isize, 
                   SEXP c2nd_qname, SEXP c2nd_cigar, SEXP c2nd_isize, SEXP impute, SEXP left, SEXP self){
  arma::uvec h1st_idx_uvec = Rcpp::as<arma::uvec>(h1st_idx);
  vector<string> c1st_qname_svec = Rcpp::as<vector<string> >(c1st_qname);
  vector<string> c2nd_qname_svec = Rcpp::as<vector<string> >(c2nd_qname);
  arma::ivec impute_isize_ivec; vector<string> c1st_cigar_svec; arma::ivec c1st_isize_ivec;
  vector<string> c2nd_cigar_svec; arma::ivec c2nd_isize_ivec; int type_flag=0;
  regex head_pattern; regex tail_pattern; regex p1st; regex p2nd; string name;
  bool impute_flag = Rcpp::as<bool>(impute); size_t idx_size = h1st_idx_uvec.size(); size_t db_size = c2nd_qname_svec.size();
  arma::uvec mh1st_idx_uvec = arma::zeros<arma::uvec>(idx_size); 
  arma::uvec self_idx_uvec = arma::zeros<arma::uvec>(idx_size);
  if(impute_flag){   
    c1st_cigar_svec = Rcpp::as<vector<string> >(c1st_cigar);
    c1st_isize_ivec = Rcpp::as<arma::ivec>(c1st_isize);
    c2nd_cigar_svec = Rcpp::as<vector<string> >(c2nd_cigar);
    c2nd_isize_ivec = Rcpp::as<arma::ivec>(c2nd_isize);
    head_pattern = regex("^(?<size>[0-9]+)S"); tail_pattern = regex("(?<size>[0-9]+)S$");
    p1st=tail_pattern; p2nd=head_pattern; name="size";
    impute_isize_ivec.resize(idx_size);
    type_flag = (int)(Rcpp::as<bool>(left))*2 + (int)(Rcpp::as<bool>(self)); // >=2 => left; odd => self; use switch
    //cout<<"type_flag="<<type_flag<<endl;
    if(type_flag>=2) { p1st=head_pattern; p2nd=tail_pattern; }
  }  
  hashmap qname_h1st = hashmap(); hashmap::iterator idx_it;
//#pragma omp parallel for shared(qname_h1st, c1st_qname_svec, h1st_idx_uvec, idx_size)
  for(size_t i=0; i<idx_size; i++) { qname_h1st[c1st_qname_svec[h1st_idx_uvec[i]-1]]=h1st_idx_uvec[i]; } //i is 0-based, 1-based stored in hashmap
  //cout<<"length of c2nd_MPRs"<<db_size;
  
  size_t j=0;
//#pragma omp parallel for shared(j, qname_h1st, c1st_cigar_svec, c1st_isize_ivec, c2nd_qname_svec, c2nd_cigar_svec, c2nd_isize_ivec, mh1st_idx_uvec, self_idx_uvec, impute_isize_ivec, impute_flag, type_flag, p1st, p2nd, name, db_size)
  for(size_t midx=1; midx<=db_size; midx++) { //1-based midx, 0-based j
   //if(midx%1000000==0) cout<<midx<<" of "<<" j= "<<j<<endl;
    idx_it=qname_h1st.find(c2nd_qname_svec[midx-1]);  //convert to 0-based to index cpp vector
    if(idx_it!=qname_h1st.end()){ //midx is a mate of idx
      size_t idx = idx_it->second; //this is 1-based idx
      if(impute_flag){
        int clip_1st = atoi(named_capture(c1st_cigar_svec[idx-1], p1st, name).c_str());
        int clip_2nd = atoi(named_capture(c2nd_cigar_svec[midx-1], p2nd, name).c_str());      
        switch(type_flag){
          case 3: impute_isize_ivec[j]=c1st_isize_ivec[idx-1]+clip_1st+clip_2nd; break; //impute a left self
          case 2: impute_isize_ivec[j]=c2nd_isize_ivec[midx-1]-clip_1st-clip_2nd; break; //impute a right mate
          case 1: impute_isize_ivec[j]=c1st_isize_ivec[idx-1]-clip_1st-clip_2nd; break; //impute a right self
          case 0: impute_isize_ivec[j]=c2nd_isize_ivec[midx-1]+clip_1st+clip_2nd; break; //impute a left mate
        }
      }
      mh1st_idx_uvec[j]=midx; self_idx_uvec[j]=idx; j=j+1; //however, mh1st_idx[j] may not be the mate of h1st_idx[j]
    }
  } //there is the possibility hang_idx overlap with mhang_idx, below we only keep single hang, both hang is thrown
  
  Rcpp::List ret; ret["self_idx"]=Rcpp::wrap(self_idx_uvec); ret["mate_idx"]=Rcpp::wrap(mh1st_idx_uvec); 
  ret["impute_isize"]=Rcpp::wrap(impute_isize_ivec);
  return ret;
}

SEXP rcpp_hello_world(){
  cout<<"helloworld libswan"<<endl;
  CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
  NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  List z            = List::create( x, y ) ;
  return z ;
};

bool compare_as_evt (evt x, evt y) { //return true if x.first<y.first
  return x.first < y.first;
};

bool test_compare_as_evt(){
  evt x=make_pair(10,make_pair(1,true)); 
  evt y=make_pair(20,make_pair(2,false)); 
  bool cp1=compare_as_evt(x,y);
  bool cp2=!compare_as_evt(y,x);
  cout<<"passing test_compare_as_evt:"<<boolalpha<<(cp1&&cp2)<<endl;
  return(cp1&&cp2); 
};

bool overlap_both (rdp x, rdp y) { //return true if both end of rdp x and rdp y overlap
  return  ( !(x.first.first.second < y.first.first.first || y.first.first.second < x.first.first.first)
         && !(x.first.second.second < y.first.second.first || y.first.second.second < x.first.second.first) );
};

bool test_ovlap_both(){
  rdp r1 = make_pair(make_pair(make_pair(2,5),make_pair(11,15)),1);
  rdp r2 = make_pair(make_pair(make_pair(3,6),make_pair(10,14)),1);
  rdp r3 = make_pair(make_pair(make_pair(1,7),make_pair(12,13)),1);
  rdp r4 = make_pair(make_pair(make_pair(3,4),make_pair(6,7)),1);
  rdp r5 = make_pair(make_pair(make_pair(6,7),make_pair(12,13)),1);
  rdp r6 = make_pair(make_pair(make_pair(6,7),make_pair(8,9)),1);
  bool ov1 = overlap_both(r1,r2);
  //cout<<"ov1"<<boolalpha<<ov1<<endl;
  bool ov2 = overlap_both(r1,r3);
  //cout<<"ov2"<<boolalpha<<ov2<<endl;
  bool ov3 = !overlap_both(r1,r4);
  //cout<<"ov3"<<boolalpha<<ov3<<endl;
  bool ov4 = !overlap_both(r1,r5);
  //cout<<"ov4"<<boolalpha<<ov4<<endl;
  bool ov5 = !overlap_both(r1,r6);
  //cout<<"ov5"<<boolalpha<<ov5<<endl;
  cout<<"passing test_ovlap_both:"<<boolalpha<<(ov1&&ov2&&ov3&&ov4&&ov5)<<endl;
  return (ov1&&ov2&&ov3&&ov4&&ov5);
}

map<pair<bool,string>,ull> make_offset(vector<string>& chr_svec, vector<ul>& chr_size_uvec){
  map<pair<bool,string>,ull> fwdchr_to_offset = map<pair<bool,string>,ull>();
  ull n=chr_svec.size();
  ull a=0; ull t=0;
  for(ull i=0; i<n; i++) { a=a+chr_size_uvec[i]; } 
  for(ull i=0; i<n; i++) {
    fwdchr_to_offset[make_pair(true,chr_svec[i])]=t;
    fwdchr_to_offset[make_pair(false,chr_svec[i])]=a+t;
    t=t+chr_size_uvec[i];
  }
  return fwdchr_to_offset;
};

map<ull,pair<bool,string> > make_reverse(map<pair<bool,string>,ull>& fwdchr_to_offset){
  map<ull,pair<bool,string> > offset_to_fwdchr = map<ull,pair<bool,string> >();
  map<pair<bool,string>,ull>::iterator mit;
  for(mit=fwdchr_to_offset.begin(); mit!=fwdchr_to_offset.end(); mit++) {
     offset_to_fwdchr[mit->second]=mit->first; 
  }
  return offset_to_fwdchr;
};

bool test_make_reverse_offset(){
  const char *cinit[] = {"chr1", "chr2", "chr3"};
  //vector<string> chr_svec(cinit, end(cinit));
  vector<string> chr_svec(cinit,cinit+3);
  ul uinit[] = {10, 20, 30};
  vector<ul> chr_size_uvec(uinit,uinit+3);
  //vector<ul> chr_size_uvec(uinit, end(uinit));
  map<pair<bool,string>,ull> offset=make_offset(chr_svec,chr_size_uvec);
  map<ull,pair<bool,string> > reverse_offset=make_reverse(offset);
  vector<ul> offset_keys;
  for(map<ull,pair<bool,string> >::iterator it = reverse_offset.begin(); it != reverse_offset.end(); ++it) {
    offset_keys.push_back(it->first);
    //cout << it->first << "\n";
  }
  //lower_bound returns first y such that y>=x; to make it works with strict less, use x-1
  //for(vector<ul>::iterator it=offset_keys.begin(); it != offset_keys.end(); ++it) { cout<<*it<<endl; }
  ull off1=*(--lower_bound(offset_keys.begin(),offset_keys.end(),10)); //idx1=0
  ull off2=*(--lower_bound(offset_keys.begin(),offset_keys.end(),11)); //idx2=10
  ull pos1=offset[make_pair<bool,string>(true,"chr2")];
  ull pos2=offset[make_pair<bool,string>(false,"chr2")];
  bool fwd1=reverse_offset[10].first; // true
  string ch1=reverse_offset[10].second; // chr2
  bool fwd2=reverse_offset[70].first; // false
  string ch2=reverse_offset[70].second; // chr1
  //cout<<"off1="<<off1<<",off2="<<off2<<",pos1="<<pos1<<",pos2="<<pos2<<endl;
  //cout<<"ch1="<<ch1<<",fwd1="<<boolalpha<<fwd1<<",ch2="<<ch2<<",fwd2="<<boolalpha<<fwd2<<endl;
  bool t1=(off1==0);
  bool t2=(off2==10);
  bool t3=(pos1==10);
  bool t4=(pos2==70);
  bool t5=(ch1=="chr2"&&fwd1==true); 
  bool t6=(ch2=="chr2"&&fwd2==false); 
  cout<<"passing test_make_reverse_offset:"<<boolalpha<<(t1&&t2&&t3&&t4&&t5&&t6)<<endl;
  return (t1&&t2&&t3&&t4&&t5&&t6);
};

vector<rdp> anch_to_coord( vector<bool>& Lfwd_bvec, vector<string>& Lchr_svec, vector<ul>& Lst_uvec, vector<ul>& Led_uvec,
                           vector<bool>& Rfwd_bvec, vector<string>& Rchr_svec, vector<ul>& Rst_uvec, vector<ul>& Red_uvec,
                           map<pair<bool,string>,ull>& fwdchr_to_offset ) {
  ull n=Lfwd_bvec.size(); vector<rdp> rdp_vec = vector<rdp>(n);
  for(ull i=0; i<n; i++) {
    ull Loffset = fwdchr_to_offset[make_pair(Lfwd_bvec[i],Lchr_svec[i])];
    ull Lstart = Loffset+Lst_uvec[i];
    ull Lend = Loffset+Led_uvec[i];
    ull Roffset = fwdchr_to_offset[make_pair(Rfwd_bvec[i],Rchr_svec[i])];
    ull Rstart = Roffset+Rst_uvec[i];
    ull Rend = Roffset+Red_uvec[i];
    rdp_vec[i] =  (Rstart>=Lstart) ? make_pair(make_pair(make_pair(Lstart,Lend),make_pair(Rstart,Rend)),1) : make_pair(make_pair(make_pair(Rstart,Rend),make_pair(Lstart,Lend)),1); 
  }
  return rdp_vec;
};

bool test_anch_to_coord(){
  const char *chr[] = {"chr1", "chr2", "chr3"};
  //vector<string> chr_svec(chr, end(chr));
  vector<string> chr_svec(chr,chr+3);
  ul uinit[] = {10, 20, 30};
  vector<ul> chr_size_uvec(uinit,uinit+3);
  //vector<ul> chr_size_uvec(uinit, end(uinit));
  map<pair<bool,string>,ull> offset=make_offset(chr_svec,chr_size_uvec);
  bool Lbinit[] = {true,true}; bool Rbinit[] = {false,true};
  //vector<bool> Lfwd_bvec(Lbinit,end(Lbinit)); vector<bool> Rfwd_bvec(Rbinit,end(Rbinit));
  vector<bool> Lfwd_bvec(Lbinit,Lbinit+2); vector<bool> Rfwd_bvec(Rbinit,Rbinit+2);
  const char *Lcinit[] = {"chr2", "chr2"}; const char *Rcinit[] = {"chr2", "chr3"};
  //vector<string> Lchr_svec(Lcinit,end(Lcinit)); vector<string> Rchr_svec(Rcinit,end(Rcinit));
  vector<string> Lchr_svec(Lcinit,Lcinit+2); vector<string> Rchr_svec(Rcinit,Rcinit+2);
  ul LSuinit[] = {10, 15}; ul LEuinit[] = {15, 20}; ul RSuinit[] = {15, 15}; ul REuinit[] = {20, 20};
  vector<ul> Lst_uvec(LSuinit,LSuinit+2); vector<ul> Led_uvec(LEuinit,LEuinit+2);
  //vector<ul> Lst_uvec(LSuinit,end(LSuinit)); vector<ul> Led_uvec(LEuinit,end(LEuinit));
  vector<ul> Rst_uvec(RSuinit,RSuinit+2); vector<ul> Red_uvec(REuinit,REuinit+2);
  //vector<ul> Rst_uvec(RSuinit,end(RSuinit)); vector<ul> Red_uvec(REuinit,end(REuinit));
  vector<rdp> coord = anch_to_coord(Lfwd_bvec,Lchr_svec,Lst_uvec,Led_uvec,Rfwd_bvec,Rchr_svec,Rst_uvec,Red_uvec,offset);
  //cout<<coord[0].first.first.first<<coord[0].first.first.second<<coord[0].first.second.first<<coord[0].first.second.second<<endl;
  //cout<<coord[1].first.first.first<<coord[1].first.first.second<<coord[1].first.second.first<<coord[1].first.second.second<<endl;
  bool t1=(coord[0].first.first.first==20&&coord[0].first.second.second==90);
  bool t2=(coord[1].first.first.first==25&&coord[1].first.second.second==50);
  cout<<"passing test_anch_to_coord:"<<boolalpha<<(t1&&t2)<<endl;
  return (t1&&t2);
};

vector<evt> coord_to_evt( vector<rdp>& rdp_vec ){
  ull n=rdp_vec.size(); 
  vector<evt> evt_vec=vector<evt>(2*n);
  for(ull i=0; i<n; i++) {
    evt_vec[2*i]=make_pair(rdp_vec[i].first.first.first,make_pair(i,true)); // idx, in?
    evt_vec[2*i+1]=make_pair(rdp_vec[i].first.second.second,make_pair(i,false));
  }
  return evt_vec;
};

bool test_coord_to_evt(){
  const char *chr[] = {"chr1", "chr2", "chr3"};
  //vector<string> chr_svec(chr, end(chr));
  vector<string> chr_svec(chr,chr+3);
  ul uinit[] = {10, 20, 30};
  //vector<ul> chr_size_uvec(uinit, end(uinit));
  vector<ul> chr_size_uvec(uinit,uinit+3);
  map<pair<bool,string>,ull> offset=make_offset(chr_svec,chr_size_uvec);
  bool Lbinit[] = {true,true}; bool Rbinit[] = {false,true};
  //vector<bool> Lfwd_bvec(Lbinit,end(Lbinit)); vector<bool> Rfwd_bvec(Rbinit,end(Rbinit));
  vector<bool> Lfwd_bvec(Lbinit,Lbinit+2); vector<bool> Rfwd_bvec(Rbinit,Rbinit+2);
  const char *Lcinit[] = {"chr2", "chr2"}; const char *Rcinit[] = {"chr2", "chr3"};
  //vector<string> Lchr_svec(Lcinit,end(Lcinit)); vector<string> Rchr_svec(Rcinit,end(Rcinit));
  vector<string> Lchr_svec(Lcinit,Lcinit+2); vector<string> Rchr_svec(Rcinit,Rcinit+2);
  ul LSuinit[] = {10, 15}; ul LEuinit[] = {15, 20}; ul RSuinit[] = {15, 15}; ul REuinit[] = {20, 20};
  //vector<ul> Lst_uvec(LSuinit,end(LSuinit)); vector<ul> Led_uvec(LEuinit,end(LEuinit));
  //vector<ul> Rst_uvec(RSuinit,end(RSuinit)); vector<ul> Red_uvec(REuinit,end(REuinit));
  vector<ul> Lst_uvec(LSuinit,LSuinit+2); vector<ul> Led_uvec(LEuinit,LEuinit+2);
  vector<ul> Rst_uvec(RSuinit,RSuinit+2); vector<ul> Red_uvec(REuinit,REuinit+2);
  vector<rdp> coord = anch_to_coord(Lfwd_bvec,Lchr_svec,Lst_uvec,Led_uvec,Rfwd_bvec,Rchr_svec,Rst_uvec,Red_uvec,offset);
  vector<evt> evt = coord_to_evt(coord);
  //cout<<evt[0].first<<boolalpha<<evt[0].second.second<<endl;
  //cout<<evt[1].first<<boolalpha<<evt[1].second.second<<endl;
  //cout<<evt[2].first<<boolalpha<<evt[2].second.second<<endl;
  //cout<<evt[3].first<<boolalpha<<evt[3].second.second<<endl;
  bool t1=(evt[0].first==20&&evt[0].second.second);
  bool t2=(evt[1].first==90&&!evt[1].second.second);
  bool t3=(evt[2].first==25&&evt[2].second.second);
  bool t4=(evt[3].first==50&&!evt[3].second.second);
  cout<<"passing test_coord_to_event:"<<boolalpha<<(t1&&t2&&t3&&t4)<<endl;
  return (t1&&t2&&t3&&t4);
}

vector<clust> scan_clust(vector<evt>& evt_vec, vector<rdp>& rdp_vec, bool (*ov)(rdp,rdp)){
  vector<clust> clust_vec; map<rdp,ul> rdp_active; ull clust_id; 
  vector<ull> rdp_clust=vector<ull>(rdp_vec.size());
  for(ull i=0; i<evt_vec.size(); i++){
    //cout<<i<<":"<<evt_vec[i].first<<":"<<evt_vec[i].second.first<<":"<<evt_vec[i].second.second<<endl;
    if(!evt_vec[i].second.second) { //out
      rdp_active.erase(rdp_vec[evt_vec[i].second.first]); continue;
    }
    clust_id=clust_vec.size(); //assume a new cluster unless found overlap
    for(map<rdp,ul>::iterator it=rdp_active.begin(); it!=rdp_active.end(); ++it){
      if(ov(rdp_vec[evt_vec[i].second.first],(*it).first)) { 
        clust_id=rdp_clust[(*it).second]; break; 
      }
    }
    rdp_active[rdp_vec[evt_vec[i].second.first]]=evt_vec[i].second.first; 
    rdp_clust[evt_vec[i].second.first]=clust_id;
    if(clust_id==clust_vec.size()) { 
      clust_vec.push_back(clust(1,evt_vec[i].second.first)); } 
    else { 
      clust_vec[clust_id].push_back(evt_vec[i].second.first); }
  }
  return clust_vec;
};

bool test_scan_clust(){
  const char *chr[] = {"chr1", "chr2", "chr3"};
  //vector<string> chr_svec(chr, end(chr));
  vector<string> chr_svec(chr,chr+3);
  ul uinit[] = {10, 20, 30};
  //vector<ul> chr_size_uvec(uinit, end(uinit));
  vector<ul> chr_size_uvec(uinit,uinit+3);
  map<pair<bool,string>,ull> offset=make_offset(chr_svec,chr_size_uvec);
  bool Lbinit[] = {true,true,true}; bool Rbinit[] = {false,true,false};
  //vector<bool> Lfwd_bvec(Lbinit,end(Lbinit)); vector<bool> Rfwd_bvec(Rbinit,end(Rbinit));
  vector<bool> Lfwd_bvec(Lbinit,Lbinit+3); vector<bool> Rfwd_bvec(Rbinit,Rbinit+3);
  const char *Lcinit[] = {"chr2", "chr2", "chr2"}; 
  const char *Rcinit[] = {"chr2", "chr3", "chr2"};
  //vector<string> Lchr_svec(Lcinit,end(Lcinit)); vector<string> Rchr_svec(Rcinit,end(Rcinit));
  vector<string> Lchr_svec(Lcinit,Lcinit+3); vector<string> Rchr_svec(Rcinit,Rcinit+3);
  ul LSuinit[] = {10, 15, 10}; ul LEuinit[] = {15, 20, 15}; 
  ul RSuinit[] = {15, 15, 15}; ul REuinit[] = {20, 20, 20};
  //vector<ul> Lst_uvec(LSuinit,end(LSuinit)); vector<ul> Led_uvec(LEuinit,end(LEuinit));
  vector<ul> Lst_uvec(LSuinit,LSuinit+3); vector<ul> Led_uvec(LEuinit,LEuinit+3);
  //vector<ul> Rst_uvec(RSuinit,end(RSuinit)); vector<ul> Red_uvec(REuinit,end(REuinit));
  vector<ul> Rst_uvec(RSuinit,RSuinit+3); vector<ul> Red_uvec(REuinit,REuinit+3);
  vector<rdp> rdp_vec = anch_to_coord(Lfwd_bvec,Lchr_svec,Lst_uvec,Led_uvec,Rfwd_bvec,Rchr_svec,Rst_uvec,Red_uvec,offset);
  vector<evt> evt_vec = coord_to_evt(rdp_vec);
  stable_sort(evt_vec.begin(), evt_vec.end(), compare_as_evt);
  vector<clust> clust_vec = scan_clust(evt_vec, rdp_vec, overlap_both);
  //cout<<clust_vec.size()<<endl;
  //for( clust::iterator i = clust_vec[0].begin(); i != clust_vec[0].end(); ++i)
  //  cout << *i << ' ';
  //cout<<endl;
  //for( clust::iterator i = clust_vec[1].begin(); i != clust_vec[1].end(); ++i)
  //  cout << *i << ' ';
  //cout<<endl;
  bool t1=clust_vec[0][0]==0&&clust_vec[0][1]==2&&clust_vec[1][0]==1;
  cout<<"passing test_scan_clust:"<<boolalpha<<t1<<endl;
  return t1;  
};

rdp collapse_clust(clust clust_item, vector<rdp>& rdp_vec){
  ull Lst=numeric_limits<ull>::max(); ull Led=0; 
  ull Rst=Lst=numeric_limits<ull>::max(); ull Red=0; 
  for(vector<ull>::iterator it=clust_item.begin(); it!=clust_item.end(); ++it){
    if(rdp_vec[(*it)].first.first.first<Lst) Lst=rdp_vec[(*it)].first.first.first; 
    if(rdp_vec[(*it)].first.first.second>Led) Led=rdp_vec[(*it)].first.first.second; 
    if(rdp_vec[(*it)].first.second.first<Rst) Rst=rdp_vec[(*it)].first.second.first; 
    if(rdp_vec[(*it)].first.second.second>Red) Red=rdp_vec[(*it)].first.second.second; 
  }
  rdp clap=make_pair(make_pair(make_pair(Lst,Led),make_pair(Rst,Red)),clust_item.size());
  return clap;
}

rdp collapse_contrast(clust clust_item, vector<rdp>& rdp_vec, int x_sup, int y_sup){
  ull Lst=numeric_limits<ull>::max(); ull Led=0;
  ull Rst=Lst=numeric_limits<ull>::max(); ull Red=0;
  int xs=0; int ys=0; //report cluster only if combined contrast >-y_sup, combined support >=x_sup 
  for(vector<ull>::iterator it=clust_item.begin(); it!=clust_item.end(); ++it){
    if(rdp_vec[(*it)].second >0) { //test if is a positive supporting cluster
      xs += rdp_vec[(*it)].second; 
      if(rdp_vec[(*it)].first.first.first<Lst) Lst=rdp_vec[(*it)].first.first.first;
      if(rdp_vec[(*it)].first.first.second>Led) Led=rdp_vec[(*it)].first.first.second;
      if(rdp_vec[(*it)].first.second.first<Rst) Rst=rdp_vec[(*it)].first.second.first;
      if(rdp_vec[(*it)].first.second.second>Red) Red=rdp_vec[(*it)].first.second.second;
    } else {
      ys += rdp_vec[(*it)].second;
    }
  }
  rdp clap=make_pair(make_pair(make_pair(0,0),make_pair(0,0)),0); //a NULL result
  if(xs>=x_sup && !(ys<=y_sup))
    clap=make_pair(make_pair(make_pair(Lst,Led),make_pair(Rst,Red)),xs);
  return clap;
}

SEXP anch_clust(SEXP Lfwd, SEXP Lchr, SEXP Lst, SEXP Led, SEXP Rfwd, SEXP Rchr, SEXP Rst, SEXP Red, SEXP chr, SEXP chr_size){
  vector<bool> Lfwd_bvec = Rcpp::as<vector<bool> >(Lfwd);
  vector<string> Lchr_svec = Rcpp::as<vector<string> >(Lchr);
  vector<ul> Lst_uvec = Rcpp::as<vector<ul> >(Lst);
  vector<ul> Led_uvec = Rcpp::as<vector<ul> >(Led);
  vector<bool> Rfwd_bvec = Rcpp::as<vector<bool> >(Rfwd);
  vector<string> Rchr_svec = Rcpp::as<vector<string> >(Rchr);
  vector<ul> Rst_uvec = Rcpp::as<vector<ul> >(Rst);
  vector<ul> Red_uvec = Rcpp::as<vector<ul> >(Red);
  vector<string> chr_svec = Rcpp::as<vector<string> >(chr);
  vector<ul> chr_size_uvec = Rcpp::as<vector<ul> >(chr_size);
  
  /*
  for(vector<bool>::iterator it=Lfwd_bvec.begin(); it!=Lfwd_bvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=Lchr_svec.begin(); it!=Lchr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Lst_uvec.begin(); it!=Lst_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Led_uvec.begin(); it!=Led_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<bool>::iterator it=Rfwd_bvec.begin(); it!=Rfwd_bvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=Rchr_svec.begin(); it!=Rchr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Rst_uvec.begin(); it!=Rst_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Red_uvec.begin(); it!=Red_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=chr_svec.begin(); it!=chr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=chr_size_uvec.begin(); it!=chr_size_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  */

  map<pair<bool,string>,ull> fwdchr_to_offset = make_offset(chr_svec, chr_size_uvec);
  map<ull,pair<bool,string> > offset_to_fwdchr = make_reverse(fwdchr_to_offset);
  vector<rdp> rdp_vec = anch_to_coord(Lfwd_bvec,Lchr_svec,Lst_uvec,Led_uvec,Rfwd_bvec,Rchr_svec,Rst_uvec,Red_uvec,fwdchr_to_offset);
  /*
  for(ull i=0; i<rdp_vec.size(); i++){
    cout<<rdp_vec[i].first.first.first<<" "<<rdp_vec[i].first.first.second<<" ";
    cout<<rdp_vec[i].first.second.first<<" "<<rdp_vec[i].first.second.second<<endl;
  }*/
  vector<evt> evt_vec = coord_to_evt(rdp_vec);
  stable_sort(evt_vec.begin(), evt_vec.end(), compare_as_evt);
  vector<clust> clust_vec = scan_clust(evt_vec,rdp_vec,overlap_both);

  ull n=clust_vec.size();
  vector<rdp> clust_rdp_vec=vector<rdp>(n);
  vector<bool> anch_Lfwd=vector<bool>(n);
  vector<string> anch_Lchr=vector<string>(n);
  vector<ul> anch_Lst=vector<ul>(n);
  vector<ul> anch_Led=vector<ul>(n);
  vector<bool> anch_Rfwd=vector<bool>(n);
  vector<string> anch_Rchr=vector<string>(n);
  vector<ul> anch_Rst=vector<ul>(n);
  vector<ul> anch_Red=vector<ul>(n);
  vector<ul> clust_Sup=vector<ul>(n);
  vector<ul> offset_keys;

  for(map<ull,pair<bool,string> >::iterator it = offset_to_fwdchr.begin(); it != offset_to_fwdchr.end(); ++it) {
    offset_keys.push_back(it->first);
    //cout<<(it->second).second<<" "<<(it->second).first<<" "<<(it->first)<<endl;
  }

  /*
  for(vector<clust>::iterator it=clust_vec.begin(); it != clust_vec.end(); ++it) {
    for(clust::iterator cit=(*it).begin(); cit != (*it).end(); ++cit) cout<<(*cit)<<" ";
    cout<<endl;
  }
  */
  
  for(ull i=0; i<n; i++){
    clust_rdp_vec[i] = collapse_clust(clust_vec[i],rdp_vec);
    //cout<<clust_rdp_vec[i].first.first.first<<" "<<clust_rdp_vec[i].first.first.second<<" ";
    //cout<<clust_rdp_vec[i].first.second.first<<" "<<clust_rdp_vec[i].first.second.second<<endl;
  }

  for(ull i=0; i<n; i++){
    ull offset = *(--lower_bound(offset_keys.begin(),offset_keys.end(),clust_rdp_vec[i].first.first.first)); //
    anch_Lfwd[i] = offset_to_fwdchr[offset].first;
    anch_Lchr[i] = offset_to_fwdchr[offset].second;
    anch_Lst[i] = clust_rdp_vec[i].first.first.first - offset;
    anch_Led[i] = clust_rdp_vec[i].first.first.second - offset;
    offset = *(--lower_bound(offset_keys.begin(),offset_keys.end(),clust_rdp_vec[i].first.second.first)); //
    anch_Rfwd[i] = offset_to_fwdchr[offset].first;
    anch_Rchr[i] = offset_to_fwdchr[offset].second;
    anch_Rst[i] = clust_rdp_vec[i].first.second.first - offset;
    anch_Red[i] = clust_rdp_vec[i].first.second.second - offset;
    clust_Sup[i] = clust_rdp_vec[i].second;
  }

  Rcpp::List ret; ret["Lfwd"]=Rcpp::wrap(anch_Lfwd); ret["Lchr"]=Rcpp::wrap(anch_Lchr); ret["Lst"]=Rcpp::wrap(anch_Lst); ret["Led"]=Rcpp::wrap(anch_Led);
  ret["Rfwd"]=Rcpp::wrap(anch_Rfwd); ret["Rchr"]=Rcpp::wrap(anch_Rchr); ret["Rst"]=Rcpp::wrap(anch_Rst); ret["Red"]=Rcpp::wrap(anch_Red); ret["Sup"]=Rcpp::wrap(clust_Sup);
  
  return ret;
};

SEXP contrast_clust(SEXP xLfwd, SEXP xLchr, SEXP xLst, SEXP xLed, SEXP xRfwd, SEXP xRchr, SEXP xRst, SEXP xRed,
                    SEXP yLfwd, SEXP yLchr, SEXP yLst, SEXP yLed, SEXP yRfwd, SEXP yRchr, SEXP yRst, SEXP yRed,
                    SEXP chr, SEXP chr_size, SEXP xSup, SEXP ySup){
  vector<bool> xLfwd_bvec = Rcpp::as<vector<bool> >(xLfwd);
  vector<string> xLchr_svec = Rcpp::as<vector<string> >(xLchr);
  vector<ul> xLst_uvec = Rcpp::as<vector<ul> >(xLst);
  vector<ul> xLed_uvec = Rcpp::as<vector<ul> >(xLed);
  vector<bool> xRfwd_bvec = Rcpp::as<vector<bool> >(xRfwd);
  vector<string> xRchr_svec = Rcpp::as<vector<string> >(xRchr);
  vector<ul> xRst_uvec = Rcpp::as<vector<ul> >(xRst);
  vector<ul> xRed_uvec = Rcpp::as<vector<ul> >(xRed);
  vector<bool> yLfwd_bvec = Rcpp::as<vector<bool> >(yLfwd);
  vector<string> yLchr_svec = Rcpp::as<vector<string> >(yLchr);
  vector<ul> yLst_uvec = Rcpp::as<vector<ul> >(yLst);
  vector<ul> yLed_uvec = Rcpp::as<vector<ul> >(yLed);
  vector<bool> yRfwd_bvec = Rcpp::as<vector<bool> >(yRfwd);
  vector<string> yRchr_svec = Rcpp::as<vector<string> >(yRchr);
  vector<ul> yRst_uvec = Rcpp::as<vector<ul> >(yRst);
  vector<ul> yRed_uvec = Rcpp::as<vector<ul> >(yRed);
  vector<string> chr_svec = Rcpp::as<vector<string> >(chr);
  vector<ul> chr_size_uvec = Rcpp::as<vector<ul> >(chr_size);
  int x_sup=Rcpp::as<int>(xSup);
  int y_sup=Rcpp::as<int>(ySup);
  /*
  for(vector<bool>::iterator it=Lfwd_bvec.begin(); it!=Lfwd_bvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=Lchr_svec.begin(); it!=Lchr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Lst_uvec.begin(); it!=Lst_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Led_uvec.begin(); it!=Led_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<bool>::iterator it=Rfwd_bvec.begin(); it!=Rfwd_bvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=Rchr_svec.begin(); it!=Rchr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Rst_uvec.begin(); it!=Rst_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=Red_uvec.begin(); it!=Red_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<string>::iterator it=chr_svec.begin(); it!=chr_svec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  for(vector<ul>::iterator it=chr_size_uvec.begin(); it!=chr_size_uvec.end(); it++) cout<<(*it)<<" "; cout<<endl;
  */

  map<pair<bool,string>,ull> fwdchr_to_offset = make_offset(chr_svec, chr_size_uvec);
  map<ull,pair<bool,string> > offset_to_fwdchr = make_reverse(fwdchr_to_offset);
  vector<rdp> x_rdp_vec = anch_to_coord(xLfwd_bvec,xLchr_svec,xLst_uvec,xLed_uvec,xRfwd_bvec,xRchr_svec,xRst_uvec,xRed_uvec,fwdchr_to_offset);
  vector<rdp> y_rdp_vec = anch_to_coord(yLfwd_bvec,yLchr_svec,yLst_uvec,yLed_uvec,yRfwd_bvec,yRchr_svec,yRst_uvec,yRed_uvec,fwdchr_to_offset);
  /*
  for(ull i=0; i<rdp_vec.size(); i++){
    cout<<rdp_vec[i].first.first.first<<" "<<rdp_vec[i].first.first.second<<" ";
    cout<<rdp_vec[i].first.second.first<<" "<<rdp_vec[i].first.second.second<<endl;
  }*/
  vector<evt> x_evt_vec = coord_to_evt(x_rdp_vec);
  stable_sort(x_evt_vec.begin(), x_evt_vec.end(), compare_as_evt);
  vector<evt> y_evt_vec = coord_to_evt(y_rdp_vec);
  stable_sort(y_evt_vec.begin(), y_evt_vec.end(), compare_as_evt);
  vector<clust> x_clust_vec = scan_clust(x_evt_vec,x_rdp_vec,overlap_both);
  vector<clust> y_clust_vec = scan_clust(y_evt_vec,y_rdp_vec,overlap_both);
  vector<rdp> x_clust_rdp_vec(x_clust_vec.size());
  vector<rdp> y_clust_rdp_vec(y_clust_vec.size());
  for(ull i=0; i<x_clust_vec.size(); i++){
    x_clust_rdp_vec[i] = collapse_clust(x_clust_vec[i],x_rdp_vec);
    //cout<<clust_rdp_vec[i].first.first.first<<" "<<clust_rdp_vec[i].first.first.second<<" ";
    //cout<<clust_rdp_vec[i].first.second.first<<" "<<clust_rdp_vec[i].first.second.second<<endl;
  }
  for(ull i=0; i<y_clust_vec.size(); i++){ //switch support number sign: + sample, - contrast event
    y_clust_rdp_vec[i] = collapse_clust(y_clust_vec[i],y_rdp_vec);
    y_clust_rdp_vec[i].second = -y_clust_rdp_vec[i].second;
    //cout<<clust_rdp_vec[i].first.first.first<<" "<<clust_rdp_vec[i].first.first.second<<" ";
    //cout<<clust_rdp_vec[i].first.second.first<<" "<<clust_rdp_vec[i].first.second.second<<endl;
  }

  vector<rdp> m_clust_rdp_vec; m_clust_rdp_vec.reserve(x_clust_rdp_vec.size()+y_clust_rdp_vec.size());
  m_clust_rdp_vec.insert(m_clust_rdp_vec.end(), x_clust_rdp_vec.begin(), x_clust_rdp_vec.end());
  m_clust_rdp_vec.insert(m_clust_rdp_vec.end(), y_clust_rdp_vec.begin(), y_clust_rdp_vec.end());
  vector<evt> m_clust_evt_vec = coord_to_evt(m_clust_rdp_vec);
  stable_sort(m_clust_evt_vec.begin(), m_clust_evt_vec.end(), compare_as_evt);
  vector<clust> m_clust_vec=scan_clust(m_clust_evt_vec, m_clust_rdp_vec, overlap_both);
  vector<rdp> clust_rdp_vec;
  for(ull i=0; i<m_clust_vec.size(); i++){
    rdp r=collapse_contrast(m_clust_vec[i],m_clust_rdp_vec,x_sup,y_sup); 
    if(r.second>0) clust_rdp_vec.push_back(r);
  }

  int n=clust_rdp_vec.size();
  vector<bool> anch_Lfwd=vector<bool>(n);
  vector<string> anch_Lchr=vector<string>(n);
  vector<ul> anch_Lst=vector<ul>(n);
  vector<ul> anch_Led=vector<ul>(n);
  vector<bool> anch_Rfwd=vector<bool>(n);
  vector<string> anch_Rchr=vector<string>(n);
  vector<ul> anch_Rst=vector<ul>(n);
  vector<ul> anch_Red=vector<ul>(n);
  vector<ul> clust_Sup=vector<ul>(n);
  vector<ul> offset_keys;

  for(map<ull,pair<bool,string> >::iterator it = offset_to_fwdchr.begin(); it != offset_to_fwdchr.end(); ++it) {
    offset_keys.push_back(it->first);
    //cout<<(it->second).second<<" "<<(it->second).first<<" "<<(it->first)<<endl;
  }

  /*
  for(vector<clust>::iterator it=clust_vec.begin(); it != clust_vec.end(); ++it) {
    for(clust::iterator cit=(*it).begin(); cit != (*it).end(); ++cit) cout<<(*cit)<<" ";
    cout<<endl;
  }
  */

  for(ull i=0; i<n; i++){
    ull offset = *(--lower_bound(offset_keys.begin(),offset_keys.end(),clust_rdp_vec[i].first.first.first)); //
    anch_Lfwd[i] = offset_to_fwdchr[offset].first;
    anch_Lchr[i] = offset_to_fwdchr[offset].second;
    anch_Lst[i] = clust_rdp_vec[i].first.first.first - offset;
    anch_Led[i] = clust_rdp_vec[i].first.first.second - offset;
    offset = *(--lower_bound(offset_keys.begin(),offset_keys.end(),clust_rdp_vec[i].first.second.first)); //
    anch_Rfwd[i] = offset_to_fwdchr[offset].first;
    anch_Rchr[i] = offset_to_fwdchr[offset].second;
    anch_Rst[i] = clust_rdp_vec[i].first.second.first - offset;
    anch_Red[i] = clust_rdp_vec[i].first.second.second - offset;
    clust_Sup[i] = clust_rdp_vec[i].second;
  }

  Rcpp::List ret; ret["Lfwd"]=Rcpp::wrap(anch_Lfwd); ret["Lchr"]=Rcpp::wrap(anch_Lchr); ret["Lst"]=Rcpp::wrap(anch_Lst); ret["Led"]=Rcpp::wrap(anch_Led);
  ret["Rfwd"]=Rcpp::wrap(anch_Rfwd); ret["Rchr"]=Rcpp::wrap(anch_Rchr); ret["Rst"]=Rcpp::wrap(anch_Rst); ret["Red"]=Rcpp::wrap(anch_Red); ret["Sup"]=Rcpp::wrap(clust_Sup);
  
  return ret;
};

/* obsolete code
SEXP lDp_omp(SEXP rPs_winDp, SEXP winW_start, SEXP rPs_start, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p){
  arma::umat rPs_winDp_umat = Rcpp::as<arma::umat>(rPs_winDp);
  arma::uvec winW_start_uvec = Rcpp::as<arma::uvec>(winW_start);
  arma::uvec rPs_start_uvec = Rcpp::as<arma::uvec>(rPs_start);
  arma::vec Fx_vec = Rcpp::as<arma::vec>(Fx);
  int N = rPs_winDp_umat.n_rows; int W = Rcpp::as<int>(n_win);
  double r = Rcpp::as<double>(mixing_rate); double pp = Rcpp::as<double>(p);
  arma::vec lDp_term_vec = arma::zeros<arma::vec>(W);
  #pragma omp parallel for shared(lDp_term_vec, Fx_vec, rPs_winDp_umat, winW_start_uvec, rPs_start_uvec, r, pp, N)
  for(int i=0; i<N; i++) { //every 1-based index need to be converted to 0-based index
    double v = ::log(1+r*(1-pp)/pp*Fx_vec(winW_start_uvec(rPs_winDp_umat(i,1)-1)-rPs_start_uvec(rPs_winDp_umat(i,0)-1)));
    #pragma omp atomic
    lDp_term_vec(rPs_winDp_umat(i,1)-1) += v;     
  }
  return Rcpp::wrap(lDp_term_vec);
};

SEXP lDm_omp(SEXP rMs_winDm, SEXP winW_end, SEXP rMs_end, SEXP Fx, SEXP n_win, SEXP mixing_rate, SEXP p){
  arma::umat rMs_winDm_umat = Rcpp::as<arma::umat>(rMs_winDm);
  arma::uvec winW_end_uvec = Rcpp::as<arma::uvec>(winW_end);
  arma::uvec rMs_end_uvec = Rcpp::as<arma::uvec>(rMs_end);
  arma::vec Fx_vec = Rcpp::as<arma::vec>(Fx);
  int N = rMs_winDm_umat.n_rows; int W = Rcpp::as<int>(n_win); 
  double r = Rcpp::as<double>(mixing_rate); double pp = Rcpp::as<double>(p);
  arma::vec lDm_term_vec = arma::zeros<arma::vec>(W);
  #pragma omp parallel for shared(lDm_term_vec, Fx_vec, rMs_winDm_umat, winW_end_uvec, rMs_end_uvec, r, pp, N)
  for(int i=0; i<N; i++) { //every 1-based index need to be converted to 0-based index
    double v = ::log(1+r*(1-pp)/pp*Fx_vec(rMs_end_uvec(rMs_winDm_umat(i,0)-1)-winW_end_uvec(rMs_winDm_umat(i,1)-1)));
    #pragma omp atomic
    lDm_term_vec(rMs_winDm_umat(i,1)-1) += v;
  }
  return Rcpp::wrap(lDm_term_vec);
};
*/

SEXP swan_unit_test(){
  test_compare_as_evt();
  test_ovlap_both();
  test_make_reverse_offset();
  test_anch_to_coord();
  test_coord_to_evt();
  test_scan_clust();
  return Rcpp::wrap(NULL);
}
