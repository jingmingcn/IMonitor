/*
 * cross.cpp
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
 *
 * This file is part of COPE.
 *
 * COPE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * COPE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with COPE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *============================================================================
 *
 * File          : cross.cpp
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Jianying Yuan, Binghang Liu
 *
 *============================================================================
 */
 
#include "cross.h"

int real_insert_size2 = -1;
uint64_t Total_num = 0;
uint64_t Fail_num = 0;
uint64_t Have_cross_num = 0;
uint64_t Have_cross_connect_num = 0;

bool output_cross(ogzstream * cross_outfile, KmerSite & kmer1, KmerSite & kmer2, Pkmer & pkmer, string & read, string & qual)
{
  bool have_cross = 0;
  int insert_size = 0;
  int overlap_size = 0;
  string seq, qu;

  if((int(kmer1.s) == int(pkmer.f1)) == (int(kmer2.s) == int(pkmer.f2))) {
    if(pkmer.f) { //kmer1 and kmer2 not exchanged in the Pkmer
      if(int(kmer1.s) == int(pkmer.f1) && int(kmer2.s) == int(pkmer.f2)) {
        overlap_size = int(pkmer.l1) - int(pkmer.p1) + int(pkmer.p2) - int(int(kmer2.p) - int(kmer1.p));
      }

      else {
        overlap_size = int(pkmer.l1) - int(pkmer.p1) + int(pkmer.p2) - int(int(kmer1.p) - int(kmer2.p));
      }

      if(overlap_size < 0) {
        if(int(kmer2.p) < int(kmer1.p)) {
          seq = read.substr(int(kmer2.p) + int(pkmer.p2) + Kmer, -overlap_size);
          seq = reverse_complement(seq);
          qu = qual.substr(int(kmer2.p) + int(pkmer.p2) + Kmer, -overlap_size);
          qu = reverse_str(qu);
        }

        else {
          seq = read.substr(int(kmer1.p) + int(pkmer.l1) - int(pkmer.p1), -overlap_size);
          qu = qual.substr(int(kmer1.p) + int(pkmer.l1) - int(pkmer.p1), -overlap_size);
        }
      }

      else {
        seq = "-";
        qu = "-";
      }

      have_cross = 1;
    }

    else { //kmer1 and kmer2 exchanged in Pkmer.
      if(int(kmer1.s) == int(pkmer.f1) && int(kmer2.s) == int(pkmer.f2)) {
        overlap_size = int(pkmer.l2) - int(pkmer.p2) + int(pkmer.p1) - int(int(kmer1.p) - int(kmer2.p));
      }

      else {
        overlap_size = int(pkmer.l2) - int(pkmer.p2) + int(pkmer.p1) - int(int(kmer2.p) - int(kmer1.p));
      }

      if(overlap_size < 0) {
        if(int(kmer1.p) < int(kmer2.p)) {
          seq = read.substr(int(kmer1.p) + int(pkmer.p1) + Kmer, -overlap_size);
          seq = reverse_complement(seq);
          qu = qual.substr(int(kmer1.p) + int(pkmer.p1) + Kmer, -overlap_size);
          qu = reverse_str(qu);
        }

        else {
          seq = read.substr(int(kmer2.p) + int(pkmer.l2) - int(pkmer.p2), -overlap_size);
          qu = qual.substr(int(kmer2.p) + int(pkmer.l2) - int(pkmer.p2), -overlap_size);
        }
      }

      else {
        seq = "-";
        qu = "-";
      }

      have_cross = 1;
    }
  }

  insert_size = int(pkmer.l1) + int(pkmer.l2) - overlap_size;

  if(insert_size < int(pkmer.l1) || insert_size < int(pkmer.l2)) {
    have_cross = 0;
  }

  if(have_cross) {
    int i = pkmer.p % Batch_num;
    cross_outfile[i] << pkmer.p << " " << insert_size << " " << seq << " " << qu << endl;
  }

  return have_cross;
}


void check_cross_read(ogzstream * cross_outfile, string & read, string & qual, HashSet * kmerHash, Pkmer * pkmers, uint64_t * cross_count)
{
  Entity * e = new Entity;
  int first_num = 0;
  uint64_t oldmer = 0;
  vector< KmerSite > used_kmers;
  map< uint64_t , vector<int> > have_kmers;
  KmerSite tmp;

  for(int i = 0; i < int(read.size()) - Kmer + 1; i++) {
    string kmer = read.substr(i, Kmer);
    uint64_t mer = seq2bit(kmer);
    e->kmer = mer;
    uint64_t idx = get_hashset(kmerHash, e);

    if(idx != kmerHash->size) {
      uint64_t rcmer = get_rev_com_kbit(mer, Kmer);
      mer = rcmer < mer ? rcmer : mer;
      e->kmer = mer;
      idx = get_hashset(kmerHash, e);

      if(idx == kmerHash->size) {
        cerr << "Error! Can't not find recmer in hash table" << endl;
        continue;
      }

      tmp.k = mer;
      tmp.p = i;
      tmp.s = mer != rcmer;
      have_kmers[mer].push_back(used_kmers.size());
      tmp.i = 0;
      tmp.l = 0;

      if(kmerHash->array[idx].len > 0) {
        tmp.i = kmerHash->array[idx].pos;
        tmp.l = kmerHash->array[idx].len;
        used_kmers.push_back(tmp);
        first_num++;
      }

      else {
        used_kmers.push_back(tmp);
      }
    }
  }

  delete e;

  if(first_num == 0 || used_kmers.size() < 2) {
    return ;
  }

  for(int i = 0; i < used_kmers.size(); i++) {
    if(used_kmers[i].l == 0) {
      continue;
    }

    uint64_t end = used_kmers[i].i + used_kmers[i].l;

    for(int k = used_kmers[i].i; k < end; k++) {
      if(have_kmers.count(pkmers[k].k2)) {
        for(int y = 0; y < have_kmers[pkmers[k].k2].size(); y++) {
          int j = have_kmers[pkmers[k].k2][y];
          bool have_cross = output_cross(cross_outfile, used_kmers[i], used_kmers[j], pkmers[k], read, qual);

          if(have_cross) {
            int l = pkmers[k].p % Batch_num;
            cross_count[l]++;
          }
        }
      }
    }
  }

  have_kmers.clear();
  used_kmers.clear();
  //end of check cross read.
}

void fetch_info(string & id, string & pid, int & pos)
{
  vector< string > svec;
  str_split(id, 32, svec); //blank value.
  pid = svec[0];

  if(svec.size() > 1) {
    pos = atoi(svec[svec.size() - 1].c_str());
  }

  svec.clear();
}

//if(TESTING)
int get_insert2(string & id, string & pid)
{
  vector< string > svec;
  str_split(id, 32, svec); //blank value.
  pid = svec[0];

  if(svec.size() == 1) {
    return -1;
  }

  int insert = atoi(svec[4].c_str());
  return insert;
}

uint64_t count_pair_kmer(string & kfile)
{
  cerr << "Count pair kmer number..." << endl;
  igzstream fin_a;
  fin_a.open(kfile.c_str());

  if(!fin_a) {
    cerr << "please check file " << kfile << endl;
    exit(0);
  }

  uint64_t num = 0;
  string kmer;

  while(getline(fin_a, kmer, '\n')) {
    num++;
  }

  fin_a.close();
  cerr << "Count kmer pair number: " << num << endl;
  return num;
}

uint64_t count_cross_read(string & cross_file)
{
  cerr << "Count cross reads number..." << endl;
  igzstream fin_a;
  fin_a.open(cross_file.c_str());

  if(!fin_a) {
    cerr << "please check file " << cross_file << endl;
    exit(0);
  }

  uint64_t num = 0;
  string line;

  while(getline(fin_a, line, '\n')) {
    num++;
  }

  fin_a.close();
  cerr << "Count line number of cross.read.info: " << num << endl;
  return num;
}

void load_pair_kmer(string & kfile, Pkmer * pkmers)
{
  igzstream fin_a;
  fin_a.open(kfile.c_str());

  if(!fin_a) {
    cerr << "please check file " << kfile << endl;
    exit(0);
  }

  uint64_t num = 0, mer1 = 0, mer2 = 0, pcount = 0, rcmer1 = 0, rcmer2 = 0;
  int pos1 = 0, pos2 = 0, read_len1 = 0, read_len2 = 0;
  string kmer1, kmer2, line;
  Pkmer tmp;
  cerr << "starting to load pair kmer file...\n";

  while(fin_a >> pcount >> kmer1 >> read_len1 >> pos1 >> kmer2 >> read_len2 >> pos2) {
    mer1 = seq2bit(kmer1);
    mer2 = seq2bit(kmer2);
    rcmer1 = get_rev_com_kbit(mer1, Kmer);
    rcmer2 = get_rev_com_kbit(mer2, Kmer);
    mer1 = mer1 < rcmer1 ? mer1 : rcmer1;
    mer2 = mer2 < rcmer2 ? mer2 : rcmer2;

    if(mer1 > mer2) {
      tmp.k1 = mer2;
      tmp.k2 = mer1;
      tmp.f1 = rcmer2 != mer2;
      tmp.f2 = rcmer1 != mer1;
      tmp.f = 0; //0 is exchange.
      tmp.p = pcount;
      tmp.p1 = pos2;
      tmp.p2 = pos1;
      tmp.l1 = read_len2;
      tmp.l2 = read_len1;
    }

    else {
      tmp.f = 1;
      tmp.k1 = mer1;
      tmp.k2 = mer2;
      tmp.f1 = rcmer1 != mer1;
      tmp.f2 = rcmer2 != mer2;
      tmp.p = pcount;
      tmp.p1 = pos1;
      tmp.p2 = pos2;
      tmp.l1 = read_len1;
      tmp.l2 = read_len2;
    }

    pkmers[num] = tmp;
    getline(fin_a, line, '\n');
    num++;
  }

  fin_a.close();
  cerr << "loading kmer file finished! Totally load kmer pair number is: " << num << endl;
}

void load_hashset(Pkmer * pkmers, uint64_t pkcount, HashSet * kmerHash)
{
  Entity * e = new Entity;
  e->pos = 0;
  e->kmer = pkmers[0].k1;
  e->len = 1;
  e->flag = pkmers[0].f1;

  if(pkcount == 1) {
    add_hashset(kmerHash, e);
    //>>add reverse kmer start.
    e->kmer = get_rev_com_kbit(e->kmer, Kmer);
    e->flag = !e->flag;
    add_hashset(kmerHash, e);
    //>>add revese kmer end.
  }

  for(uint64_t i = 1; i < pkcount; i++) {
    if(e->kmer == pkmers[i].k1) {
      e->len ++;
      continue;
    }

    add_hashset(kmerHash, e);
    //>>add reverse kmer start.
    e->kmer = get_rev_com_kbit(e->kmer, Kmer);
    e->flag = !e->flag;
    add_hashset(kmerHash, e);
    //>>add revese kmer end.
    e->pos = i;
    e->kmer = pkmers[i].k1;
    e->len = 1;
    e->flag = pkmers[i].f1;
  }

  //keep other kmer used but not as the first kmer.
  for(uint64_t i = 0; i < pkcount; i++) {
    e->kmer = pkmers[i].k2;
    e->pos = 0;
    e->len = 0;
    e->flag = pkmers[i].f2;
    uint64_t idx = get_hashset(kmerHash, e);

    if(idx == kmerHash->size) {
      add_hashset(kmerHash, e);
      //>>add reverse kmer start.
      e->kmer = get_rev_com_kbit(e->kmer, Kmer);
      e->flag = !e->flag;
      add_hashset(kmerHash, e);
      //>>add revese kmer end.
    }
  }

  cerr << "final hash set size is " << kmerHash->size << endl;
  delete e;
}

void load_pmap(CrossInf * crossInf, uint64_t count, map<int, int> & pmap)
{
  for(uint64_t i = 0; i < count; i++) {
    if(!pmap.count(crossInf[i].p)) {
      pmap[crossInf[i].p] = i;
      //cerr<<i<<"\t"<<crossInf[i].p<<endl;
    }
  }
}

double get_ratio(int over, string & r1, string & r2)
{
  int match = 0;
  int mismatch = 0;
  int start = int(r1.size()) - over;

  if(start < 0 || over <= 0) {
    return 0;
  }

  for(int i = 0; i < over; i++) {
    if(r2[i] == r1[start + i]) {
      match++;
    }

    else {
      mismatch++;
    }
  }

  if(mismatch >= MaxMismatch) {
    return 0;
  }

  return double(match) / double(over);
}

void get_consensus(int insert_size, int freq, CrossInf * crossInf, uint64_t pos, int total_cross_num, int pid, string & seq, string & qual,
                   string & read_a, string & read_b, string & a_q, string & b_q)
{
  if(insert_size == read_a.size() + read_b.size()) {
    seq = read_a + read_b;
    qual = a_q + b_q;
    return ;
  }

  string seqs[freq];
  string quals[freq];
  int i = 0;

  for(uint64_t p = pos; p < total_cross_num; p++) {
    if(crossInf[p].p == pid) {
      if(int(crossInf[p].i) == insert_size) {
        seqs[i] = crossInf[p].seq;
        quals[i] = crossInf[p].qual;
        i++;
      }
    }

    else {
      break;
    }
  }

  seq = seqs[0];
  qual = quals[0];

  for(int j = 0; j < seq.size(); j++) {
    int base = alphabet2[seq[j]], maxb = base, max = 0;
    int type[5] = {0};

    for(int i = 0; i < freq; i++) {
      type[alphabet2[seqs[i][j]]]++;
    }

    for(int k = 0; k < 4; k++) {
      if(max < type[k]) {
        max = type[k];
        maxb = k;
      }
    }

    seq[j] = bases[ maxb ];
    max = phred + 40;

    for(int i = 0; i < freq; i++) {
      if(seqs[i][j] == seq[j]) {
        if(max > quals[i][j]) {
          max = quals[i][j];
        }
      }
    }

    qual[j] = max;
  }

  seq = read_a + seq + read_b;
  qual = a_q + qual + b_q;
}

bool get_cross_seq(string & read_a, string & read_b, string & a_q, string & rb_q, uint64_t pid, uint64_t pos, uint64_t total_cross_num, CrossInf * crossInf,
                   string & seq, string & qual)
{
  bool is_cross_connected = 0;
  map<int, int> inserts;
  int insert_size, total_cover = 0;

  for(uint64_t p = pos; p < total_cross_num; p++) {
    if(crossInf[p].p == pid) {
      if(int(crossInf[p].i) < int(read_a.size()) || int(crossInf[p].i) < int(read_b.size())) {
        continue;
      }

      inserts[int(crossInf[p].i)]++;
      total_cover++;
    }

    else {
      break;
    }
  }

  if(TESTING) {
    CROSSCERRINF << total_cover << "\t" << inserts.size() << "\t";
  }

  if(inserts.size() == 0) {
    return 0;
  }

  double scores[inserts.size()];
  string seqs[inserts.size()];
  string quals[inserts.size()];
  int sizes[inserts.size()];
  int real_i = 0;

  if(inserts.size() == 1) {
    is_cross_connected = 1;
    sizes[0] = int(read_a.size() + read_b.size()) - int(crossInf[pos].i);

    if(sizes[0] > 0 && inserts[int(crossInf[pos].i)] > 3) {
      scores[0] = get_ratio(sizes[0], read_a, read_b);

      if(scores[0] > MATCH_CUTOFF * DIF) {
        double tmp;
        connect_overlap(tmp, sizes[0], read_a, read_b, a_q, rb_q, seqs[0], quals[0]);
      }

      else {
        is_cross_connected = 0;
      }
    }

    else if(sizes[0] <= 0 && inserts[int(crossInf[pos].i)] > 3) {
      is_cross_connected = 1;
    }

    else {
      is_cross_connected = 0;
    }

    real_i = 0;

    if(TESTING) {
      CROSSCERRINF << "@ " << sizes[0] << "\t" << inserts[int(crossInf[pos].i)] << "\t" << scores[0] << "\t" << is_cross_connected << endl;
    }
  }

  else {
    if(total_cover > RepeatCoverNum && inserts.size() > MinCoverNum) { //max cover num and min cover num
      return 0;
    }

    map<int, int>::iterator iter = inserts.begin();
    int i = 0;
    int j = 0;
    int have_larger = 0, full_i = -1, full_num = 0;

    while(iter != inserts.end()) {
      sizes[i] = int(read_a.size() + read_b.size()) - iter->first;

      if(TESTING) {
        CROSSCERRINF << sizes[i] << "\t";
      }

      if(sizes[i] <= 0) {
        scores[i] = 1;

        if(iter->second >= total_cover * (1.0 / double(inserts.size()))) {
          have_larger++;
        }

        else {
          scores[i] = 0;
        }

        if(TESTING) {
          CROSSCERRINF << scores[i] << "\t" << have_larger << "\t";
        }
      }

      else {
        scores[i] = get_ratio(sizes[i], read_a, read_b);
      }

      if(TESTING) {
        CROSSCERRINF << scores[i] << "\t";
      }

      if(scores[i] < MATCH_CUTOFF * DIF || sizes[i] * (1 - scores[i]) >= MaxMismatch) {
        scores[i] = 0;

        if(TESTING) {
          CROSSCERRINF << "#" << "\t";
        }
      }

      if(scores[i] == 1) {
        full_num++;
        full_i = i;
      }

      if(sizes[i] > 0 && scores[i] >= MATCH_CUTOFF * DIF) {
        double tmp = 0;
        connect_overlap(tmp, sizes[i], read_a, read_b, a_q, rb_q, seqs[i], quals[i]);//from connect.cpp

        if(TESTING) {
          CROSSCERRINF << "con_over" << "\t";
        }
      }

      if(scores[i] >= MATCH_CUTOFF * DIF) {
        scores[i] = iter->second;

        if(TESTING) {
          CROSSCERRINF << "$:" << scores[i];
        }
      }

      if(TESTING) {
        CROSSCERRINF << endl;
      }

      i++;
      j += int(iter->second);
      iter++;
    }

    if(TESTING) {
      CROSSCERRINF << "%: " << have_larger << "\t" << inserts.size() << "\t";
    }

    if(have_larger > 0 && have_larger < inserts.size()) {
      is_cross_connected = 0;
    }

    else {
      real_i = get_max(scores, inserts.size(), 3, 0.4, 1, 0);

      if(real_i > 0) {
        is_cross_connected = 1;
      }

      else {
        is_cross_connected = 0;
      }
    }

    if(TESTING) {
      CROSSCERRINF << is_cross_connected << endl;
    }
  }

  if(is_cross_connected) {
    if(TESTING) {
      CROSSCERRINF << "^: " << sizes[real_i] << endl;
    }

    if(sizes[real_i] <= 0) {
      get_consensus(int(read_a.size() + read_b.size()) - sizes[real_i], inserts[ int(read_a.size() + read_b.size()) - sizes[real_i] ], crossInf, pos, total_cross_num, pid, seq, qual, read_a, read_b, a_q, rb_q);
    }

    else { //there are overlap.
      seq = seqs[real_i];
      qual = quals[real_i];
    }
  }

  return is_cross_connected;
}

void pair_read_connect(ofstream & fout1, ofstream & fout2, ofstream & out, string & fq1, string & fq2, map<int, int> & pmap, CrossInf * crossInf, uint64_t total, int batch_index)
{
  igzstream fin_a;
  fin_a.open(fq1.c_str());
  igzstream fin_b;
  fin_b.open(fq2.c_str());

  if(!fin_a || !fin_b) {
    cerr << "please check file " << fq1 << ". and file " << fq2 << endl;
    exit(0);
  }

  else {
    cerr << "loading file: \n" << fq1 << "\n" << fq2 << endl;
  }

  string seq, qual;
  uint64_t num = 0, no_marker = 0, have_cross = 0,  fail_num1 = 0, fail_num = 0, have_cross_connect = 0;
  string read_a, read_b, a_id, b_id, a_s, b_s, a_q, b_q, pid, kmer1, kmer2;
  uint64_t mer1, mer2;
  cerr << "Begin load and connect pairs...\n";
  uint64_t pair_count = 1;
  uint64_t id_count = 0;

  while(getline(fin_a, a_id, '\n')) {
    getline(fin_b, b_id, '\n');

    if(a_id[0] == '@' || b_id[0] == '>') {
      if(a_id[0] == '@') {
        getline(fin_a, read_a, '\n');
        getline(fin_b, read_b, '\n');
        getline(fin_a, a_s, '\n');
        getline(fin_a, a_q, '\n');
        getline(fin_b, b_s, '\n');
        getline(fin_b, b_q, '\n');
      }

      else {
        getline(fin_a, read_a, '\n');
        getline(fin_b, read_b, '\n');
      }

      a_id = a_id.substr(1);
      b_id = b_id.substr(1);
    }

    id_count++;

    if(id_count % Batch_num != batch_index) {
      pair_count++;
      continue;
    }

    int pos1 = -1, pos2 = -1;
    string rread_b = reverse_complement(read_b);
    string rb_q = reverse_str(b_q);
    string pid, real_a(read_a), real_b(read_b);
    fetch_info(b_id, pid, pos2);
    fetch_info(a_id, pid, pos1);

    if(TESTING) {
      real_insert_size2 = get_insert2(a_id, pid);
    }

    if(! pmap.count(pair_count)) {
      Total_num++;
      Fail_num++;
      pair_count++;

      if(TESTING) {
        CROSSCERRINF << a_id << "pmap not found" << endl;
      }

      if(a_q.size() > 0) {
        fout1 << "@" << a_id << endl
              << read_a << endl
              << a_s << endl
              << a_q << endl;
        fout2 << "@" << b_id << endl
              << read_b << endl
              << b_s << endl
              << b_q << endl;
      }

      else {
        fout1 << ">" << a_id << endl
              << read_a << endl;
        fout2 << ">" << b_id << endl
              << read_b << endl;
      }

      continue;
    }

    Total_num++;
    Have_cross_num++;
    string seq, qual;

    if(TESTING) {
      CROSSCERRINF << a_id << endl;
    }

    bool is_cross_connected = get_cross_seq(read_a, rread_b, a_q, rb_q, pair_count, pmap[pair_count], total, crossInf, seq, qual);

    if(! is_cross_connected) {
      if(a_q.size() > 0) {
        fout1 << "@" << a_id << " " << pair_count << endl
              << read_a << endl
              << a_s << endl
              << a_q << endl;
        fout2 << "@" << b_id << endl
              << read_b << endl
              << b_s << endl
              << b_q << endl;
      }

      else {
        fout1 << ">" << a_id << " " << pair_count << endl
              << read_a << endl;
        fout2 << ">" << b_id << endl
              << read_b << endl;
      }
    }

    else {
      Have_cross_connect_num++;

      if(a_q.size() > 0) {
        //trainning...
        if(TESTING && real_insert_size2 > 0 && real_insert_size2 != seq.size()) {
          FALSECON2 << pid << "\t" << seq.size()  << "\t" << real_insert_size2  << endl;
        }

        out << "@" << pid << "\t" << seq.size() << endl
            << seq << endl
            << "+" << endl
            << qual << endl;
      }

      else {
        out << ">" << pid << "\t" << seq.size() << endl
            << seq << endl;
      }
    }

    pair_count++;
  }

  fin_a.close();
  fin_b.close();
}

void load_one_read_file(ogzstream * cross_outfile, string & file, HashSet * kmerHash, Pkmer * pkmers, uint64_t * cross_count)
{
  igzstream fin_a;
  fin_a.open(file.c_str());

  if(!fin_a) {
    cerr << "fail open cross file: " << file << endl;
    return;
  }

  string read, id, s, q;
  uint64_t num = 0, keep_num = 0, fail_num = 0;
  cerr << "Begin load read file: " << file << endl;

  while(getline(fin_a, read, '\n')) {
    if(read[0] == '@' || read[0] == '>') {
      id = read.substr(1);
    }

    //load strand and sequence quality
    if(read[0] == '@') {              //for *.fq
      getline(fin_a, read, '\n');
      getline(fin_a, s, '\n');
      getline(fin_a, q, '\n');
    }

    else {
      getline(fin_a, read, '\n');
    }

    check_cross_read(cross_outfile, read, q, kmerHash, pkmers, cross_count);
    num++;

    if(num % 1000000 == 0) {
      cerr << "Process reads number: " << num << endl;
    }
  }

  fin_a.close();
  cerr << "finished load read file: " << file << " ! totally load " << num << " reads\n";
}

void load_connected_read(ogzstream * cross_outfile, string & file, HashSet * kmerHash, Pkmer * pkmers, uint64_t * cross_count)
{
  load_one_read_file(cross_outfile, file, kmerHash, pkmers, cross_count);
}

void load_raw_read(ogzstream * cross_outfile, string & a_file, string & b_file, HashSet * kmerHash, Pkmer * pkmers, uint64_t * cross_count)
{
  load_one_read_file(cross_outfile, a_file, kmerHash, pkmers, cross_count);
  load_one_read_file(cross_outfile, b_file, kmerHash, pkmers, cross_count);
}

void load_other_read(ogzstream * cross_outfile, string & other_file_list, HashSet * kmerHash, Pkmer * pkmers, uint64_t * cross_count)
{
  igzstream fin_lst;
  fin_lst.open(other_file_list.c_str());

  if(!fin_lst) {
    cerr << "fail open cross file list: " << other_file_list << endl;
    return;
  }

  string read_file;

  while(getline(fin_lst, read_file, '\n')) {
    if(read_file.size() == 0) {
      continue;
    }

    load_one_read_file(cross_outfile, read_file, kmerHash, pkmers, cross_count);
  }

  cerr << "finished load cross read file list!" << endl;
}

int comparePkmer(const void * a, const void * b)
{
  Pkmer * ia = (Pkmer *)a;
  Pkmer * ib = (Pkmer *)b;

  if(ia->k1 != ib->k1) {
    return (int)((int64_t)ia->k1 - (int64_t)ib->k1);
  }

  if(ia->k2 != ib->k2) {
    return (int)((int64_t)ia->k2 - (int64_t)ib->k2);
  }

  return (int)((int)ia->p - (int)ib->p);
}

int compareCrossInf(const void * a, const void * b)
{
  struct CrossInf * ia = (struct CrossInf *)a;
  struct CrossInf * ib = (struct CrossInf *)b;

  if(ia->p != ib->p) {
    return (int)(ia->p - ib->p);
  }

  return (int)(ia->i - ib->i);
}

void load_crossinf(igzstream & fin, CrossInf * crossInf)
{
  if(!fin) {
    cerr << "fail to creat cross read inf file!\n";
    exit(0);
  }

  int pair, insert, count = 0;
  string seq;
  CrossInf tmp;
  string bstr, qstr;

  while(fin >> tmp.p >> insert >> bstr >> qstr) {
    tmp.i = insert;
    int blen = bstr.size();
    int qlen = qstr.size();

    if(blen != qlen) {
      cerr << "Error: in cross read file, length of base string not equal to length of qual string! line number: " << count + 1 << endl;
      exit(-1);
    }

    tmp.seq = new char[blen + 1];
    strcpy(tmp.seq, bstr.c_str());
    tmp.qual = new char[blen + 1];
    strcpy(tmp.qual, qstr.c_str());
    crossInf[count] = tmp;
    count++;
    getline(fin, seq, '\n');
  }

  cerr << "load cross inf finished! and load " << count << " cross inf\n";
}

string int2str(int a)
{
  stringstream stream;
  string n = "";
  stream << a;
  stream >> n;
  return n;
}

void cross_connect_reads(ofstream & fout1, ofstream & fout2, ofstream & out,  string & fq1, string & fq2, string cross_info_file, uint64_t cross_read_count, int batch_index)
{
  cerr << "Begin to cross connect reads, processing file: " << cross_info_file << " ..." << endl;

  if(cross_read_count < 1) {
    cerr << "No cross reads, cross reads connect finished!\n";
    return ;
  }

  igzstream pin;
  pin.open(cross_info_file.c_str());

  if(!pin) {
    cerr << "fail open cross file: " << cross_info_file << endl;
    exit(3);
  }

  CrossInf * crossInf = new CrossInf[cross_read_count];
  map<int, int> pmap;
  cerr << "Loading cross read information!\n";
  load_crossinf(pin, crossInf);
  pin.close();
  cerr << "sort cross read information!\n";
  qsort(crossInf, cross_read_count, sizeof(CrossInf), compareCrossInf);
  cerr << "keep pair id in map!\n";
  load_pmap(crossInf, cross_read_count, pmap);
  cerr << "connect pairs using cross read information! and have cross read pair number is " << pmap.size() << "\n";
  pair_read_connect(fout1, fout2, out, fq1, fq2, pmap, crossInf, cross_read_count, batch_index);
  delete []crossInf;
  pmap.clear();
  cerr << "Finished processing file: " << cross_info_file << " !" << endl;
}

void cross_connect(string & fq1, string & fq2, string & connect_file, string & cross_read_list_file, string & kmer_file)
{
  cerr << "Use connected file to cross\n";
  string left1 = "cross_left_1.fq";
  string left2 = "cross_left_2.fq";
  string con = "cross_connect.fq";
  ofstream fout1(left1.c_str());
  ofstream fout2(left2.c_str());
  ofstream out(con.c_str());

  if(mode == 4) {
    igzstream CrossInfoList;
    CrossInfoList.open(Cross_info_list_file.c_str());

    if(!CrossInfoList) {
      cerr << "fail open cross information list file: " << Cross_info_list_file << endl;
      return;
    }

    string cross_info_file;
    int i = 0;

    while(getline(CrossInfoList, cross_info_file, '\n')) {
      if(cross_info_file.size() == 0) {
        continue;
      }

      uint64_t cross_read_count = count_cross_read(cross_info_file);
      cross_connect_reads(fout1, fout2, out, fq1, fq2, cross_info_file, cross_read_count, i);
      i++;
    }

    CrossInfoList.close();
  }

  else {
    string cross_file_tmp = "cross.read.inf";
    string * Cross_outfile = new string[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      Cross_outfile[i] = cross_file_tmp + int2str(i) + ".gz" ;
    }

    ogzstream * CrossOut = new ogzstream[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      CrossOut[i].open(Cross_outfile[i].c_str());

      if(!CrossOut[i]) {
        cerr << "fail open cross file: " << Cross_outfile[i] << endl;
        exit(3);
      }
    }

    kmerHash = init_hashset(initial_size, load_factor);
    uint64_t pkcount = count_pair_kmer(kmer_file);
    Pkmer * pkmers = new Pkmer[pkcount];
    cerr << "loading pair kmer information!\n";
    load_pair_kmer(kmer_file, pkmers);
    cerr << "sort pair kmers!\n";
    qsort(pkmers, pkcount, sizeof(Pkmer), comparePkmer);
    cerr << "load hashset!\n";
    load_hashset(pkmers, pkcount, kmerHash);
    cerr << "load all reads to get cross reads\n";
    uint64_t * Cross_read_count = new uint64_t[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      Cross_read_count[i] = 0;
    }

    load_connected_read(CrossOut, connect_file, kmerHash, pkmers, Cross_read_count); //use kmerHash.

    if(cross_read_list_file.size() > 0) {
      load_other_read(CrossOut, cross_read_list_file, kmerHash, pkmers, Cross_read_count);
    }

    uint64_t total_count = 0;

    for(int i = 0; i < Batch_num; i++) {
      total_count += Cross_read_count[i];
    }

    cerr << "Get cross read finished and totally cross " << total_count << " reads!\n";

    for(int i = 0; i < Batch_num; i++) {
      CrossOut[i].close();
    }

    free_hash(kmerHash);
    delete []pkmers;
    delete []CrossOut;

    if(Stop_connect_cross_read) {
      delete[] Cross_read_count;
      delete[] Cross_outfile;
      return;
    }

    for(int i = 0; i < Batch_num; i++) {
      cross_connect_reads(fout1, fout2, out, fq1, fq2, Cross_outfile[i], Cross_read_count[i], i);
    }

    delete[] Cross_outfile;
    delete[] Cross_read_count;
  }

  cerr << "All done!\n";
  cout << "\nCross connection table:\ntotal_pairs\tno_marker\thave_cross\thave_cross_ratio(%)\thave_cross_connect\tcross_connect_ratio(%)\n"
       << Total_num << "\t" << Fail_num << "\t" << Have_cross_num << "\t" << double(Have_cross_num) / double(Total_num) * 100.0 << "\t" << Have_cross_connect_num << "\t" << double(Have_cross_connect_num) / double(Have_cross_num) * 100.0 << endl;
  fout1.close();
  fout2.close();
  out.close();
}

void cross_connect(string & fq1, string & fq2, string & a_file, string & b_file, string & cross_read_list_file, string & kmer_file)
{
  cerr << "Use raw file to cross\n";
  string left1 = "cross_left_1.fq";
  string left2 = "cross_left_2.fq";
  string con = "cross_connect.fq";
  ofstream fout1(left1.c_str());
  ofstream fout2(left2.c_str());
  ofstream out(con.c_str());

  if(mode == 4) {
    igzstream CrossInfoList;
    CrossInfoList.open(Cross_info_list_file.c_str());

    if(!CrossInfoList) {
      cerr << "fail open cross information list file: " << Cross_info_list_file << endl;
      return;
    }

    string cross_info_file;
    int i = 0;

    while(getline(CrossInfoList, cross_info_file, '\n')) {
      if(cross_info_file.size() == 0) {
        continue;
      }

      uint64_t cross_read_count = count_cross_read(cross_info_file);
      cross_connect_reads(fout1, fout2, out, fq1, fq2, cross_info_file, cross_read_count, i);
      i++;
    }

    CrossInfoList.close();
  }

  else {
    string cross_file_tmp = "cross.read.inf";
    string * Cross_outfile = new string[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      Cross_outfile[i] = cross_file_tmp + int2str(i) + ".gz" ;
    }

    ogzstream * CrossOut = new ogzstream[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      CrossOut[i].open(Cross_outfile[i].c_str());

      if(!CrossOut[i]) {
        cerr << "fail open cross file: " << Cross_outfile[i] << endl;
        exit(3);
      }
    }

    kmerHash = init_hashset(initial_size, load_factor);
    Pkmer * pkmers;
    uint64_t pkcount = count_pair_kmer(kmer_file);
    pkmers = new Pkmer[pkcount];
    cerr << "loading pair kmer information!\n";
    load_pair_kmer(kmer_file, pkmers);
    cerr << "sort pair kmers!\n";
    qsort(pkmers, pkcount, sizeof(Pkmer), comparePkmer);
    cerr << "load hashset!\n";
    load_hashset(pkmers, pkcount, kmerHash);
    cerr << "load all reads to get cross reads\n";
    uint64_t * Cross_read_count = new uint64_t[Batch_num];

    for(int i = 0; i < Batch_num; i++) {
      Cross_read_count[i] = 0;
    }

    load_raw_read(CrossOut, a_file, b_file, kmerHash, pkmers, Cross_read_count); //use kmerHash.

    if(cross_read_list_file.size() > 0) {
      load_other_read(CrossOut, cross_read_list_file, kmerHash, pkmers, Cross_read_count);
    }

    uint64_t total_count = 0;

    for(int i = 0; i < Batch_num; i++) {
      total_count += Cross_read_count[i];
    }

    cerr << "Get cross read finished and totally cross " << total_count << " reads!\n";

    for(int i = 0; i < Batch_num; i++) {
      CrossOut[i].close();
    }

    free_hash(kmerHash);
    delete []pkmers;
    delete []CrossOut;

    if(Stop_connect_cross_read) {
      delete[] Cross_read_count;
      delete[] Cross_outfile;
      return;
    }

    for(int i = 0; i < Batch_num; i++) {
      cross_connect_reads(fout1, fout2, out, fq1, fq2, Cross_outfile[i], Cross_read_count[i], i);
    }

    delete[] Cross_outfile;
    delete[] Cross_read_count;
  }

  cerr << "All done!\n";
  cout << "\nCross connection table:\ntotal_pairs\tno_marker\thave_cross\thave_cross_ratio(%)\thave_cross_connect\tcross_connect_ratio(%)\n"
       << Total_num << "\t" << Fail_num << "\t" << Have_cross_num << "\t" << double(Have_cross_num) / double(Total_num) * 100.0 << "\t" << Have_cross_connect_num << "\t" << double(Have_cross_connect_num) / double(Have_cross_num) * 100.0 << endl;
  fout1.close();
  fout2.close();
  out.close();
}


