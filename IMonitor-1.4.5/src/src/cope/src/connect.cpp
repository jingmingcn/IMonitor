/*
 * connect.cpp
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
 * File          : connect.cpp
 * Revision      : 1.1.3
 * Revision_date : 2012/10/08
 * Author(s)     : Binghang Liu, Jianying Yuan
 * Update:
 *    add get_pair_id and change a_id if connected.
 *
 *============================================================================
 */
 
#include "connect.h"

//if(TESTING)
int real_insert_size = -1;
void str_split(string & line, int clim, vector<string> & svec)
{
  int start = 0, len = 0;

  for(int i = 0; i < line.size(); i++) {
    if(line[i] == clim) {
      if(len > 0) {
        svec.push_back(line.substr(start, len));
      }

      start = i + 1;
      len = 0;
      continue;
    }

    len++;
  }

  if(len > 0) {
    svec.push_back(line.substr(start, len));
  }
}

int get_insert(string & id, string & pid)
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

string get_pair_id(string &aid, string &bid)
{
	//if(aid.size() <=2 || bid.size() <=2)
	return aid + "_" + bid;
}

void get_mis_table(double match_ratio, int * table, int len)
{
  for(int i = 0; i <= len; i++) {
    table[i] = int(i * (1 - match_ratio));
  }
}

void inital_matrix(bool ** matrix, int jnum)
{
  int inum = (int)pow(2.0, double(jnum));

  for(int i = 0; i < inum; i++) {
    int num = i;

    for(int j = 0; j < jnum; j++) {
      if(num == 0) {
        matrix[i][j] = 0;
      }

      else {
        matrix[i][j] = num % 2;
        num = num / 2;
      }
    }
  }
}

//loading jellyfish result table.
void load_jellyfish_kmerfreq(string & kmer_freq_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high)
{
  cerr << "loading kmerfreq table counted by jellyfish...\n";
  //check the kmer freq format
  igzstream INTMP;
  INTMP.open(kmer_freq_file.c_str());

  if(!INTMP) {
    cerr << "Can't open the kmerfreq table counted by jellyfish, please check the input option -j !" << endl;
    exit(-1);
  }

  string line_tmp;
  int format = 0; //format: >freq\nkmer
  getline(INTMP, line_tmp, '\n');

  if(line_tmp[0] != '>') {
    format = 1;  //format: kmer\tfreq or kmer(blank)freq
  }

  INTMP.close();
  //read the frequence value
  igzstream inFile;
  inFile.open(kmer_freq_file.c_str());

  if(!inFile) {
    cerr << "Can't open the kmerfreq table counted by jellyfish, please check the input option -j !" << endl;
    exit(-1);
  }

  int freqy[256];

  for(int i = 0; i < 256; i++) {
    freqy[i] = 0;
  }

  uint64_t error_freq_num = 0, suspicious_freq_num = 0, num_total_kmers = 0,
           normal_freq_num = 0, repeat_freq_num = 0, num_effect_kmers = 0, hig_freq_num = 0;
  uint64_t count_num = 0;
  cerr << "start loading kmer freq: " << kmer_freq_file << endl;

  if(format == 0) {
    string freq_line, kmer_line;

    while(getline(inFile, freq_line, '\n')) {
      getline(inFile, kmer_line, '\n');
      string freq_seq = freq_line.substr(1);
      int frequency = atoi(freq_seq.c_str());

      if(frequency > 255) {
        frequency = 255;
      }

      if(frequency > 0) {
        uint64_t idx = seq2bit(kmer_line);
        num_total_kmers += frequency;
        num_effect_kmers++;
        freqy[frequency]++;

        if(frequency > freq_low) {
          if(frequency >= freq_high) { //repeat kmer
            repeat_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll3[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll3[rc_idx % 4];
          }

          else if(frequency > freq_normal) { //normal kmer
            normal_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll2[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll2[rc_idx % 4];
          }

          else { //suspicious kmer
            suspicious_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll1[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll1[rc_idx % 4];
          }
        }

        else { //error kmer
          error_freq_num++;
        }
      }
    }
  }

  else {
    int frequency = 0;
    string kmer_line, seq;

    while(inFile >> kmer_line >> frequency) {
      if(frequency > 255) {
        frequency = 255;
      }

      if(frequency > 0) {
        uint64_t idx = seq2bit(kmer_line);
        num_total_kmers += frequency;
        num_effect_kmers++;
        freqy[frequency]++;

        if(frequency > freq_low) {
          if(frequency >= freq_high) { //repeat kmer
            repeat_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll3[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll3[rc_idx % 4];
          }

          else if(frequency > freq_normal) { //normal kmer
            normal_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll2[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll2[rc_idx % 4];
          }

          else { //suspicious kmer
            suspicious_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll1[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll1[rc_idx % 4];
          }
        }

        else { //error kmer
          error_freq_num++;
        }
      }

      getline(inFile, seq, '\n');
    }
  }

  cerr << "kmer freq table: \n";
  double error_freq_ratio = (double)error_freq_num / (double)num_effect_kmers;
  double suspicious_freq_ratio = (double)suspicious_freq_num / (double)num_effect_kmers;
  double normal_freq_ratio = (double)normal_freq_num / (double)num_effect_kmers;
  double repeat_freq_ratio = (double)repeat_freq_num / (double)num_effect_kmers;
  cerr << "kmer_num\tnode_num\terror_freq_num\terror_ratio\tsuspicious_freq_num\tsuspicious_ratio\tnormal_freq_num\tnormal_ratio\trepeat_freq_num\trepeat_ratio\n"
       << num_total_kmers << "\t" << num_effect_kmers << "\t" << error_freq_num << "\t" << error_freq_ratio << "\t" << suspicious_freq_num << "\t"
       << suspicious_freq_ratio << "\t" << normal_freq_num << "\t" << normal_freq_ratio << "\t" << repeat_freq_num << "\t" << repeat_freq_ratio << endl;
  inFile.close();
}

//load pallal compress kmer freq table, need two file: *.freq.cz and *.freq.cz.len
void load_kmerfreq(string & kmer_freq_file, string & kmer_freq_len_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high)
{
  cerr << "loading kmerfreq table...\n";
  vector<uint64_t> BlockSizeVec;
  uint64_t SrcBlockSize = 8 * 1024 * 1024;
  uint64_t error_freq_num = 0, suspicious_freq_num = 0, normal_freq_num = 0, repeat_freq_num = 0, num_total_kmers = 0, num_effect_kmers = 0;
  uint8_t * buffer_compress = new uint8_t[SrcBlockSize];
  uint8_t * buffer_uncompress = new uint8_t[SrcBlockSize];
  ifstream infile;
  infile.open(kmer_freq_len_file.c_str());

  if(!infile) {
    cerr << "Error: can't open file: " << kmer_freq_len_file << ", please check the input file!" << endl;
    exit(-1);
  }

  string str;
  int freqy[256];

  for(int i = 0; i < 100; i++) {
    freqy[i] = 0;
  }

  while(getline(infile, str, '\n')) { //each line is block size.
    BlockSizeVec.push_back(atoi(str.c_str()));
  }

  infile.close();
  infile.open(kmer_freq_file.c_str(), ifstream::binary);

  if(!infile) {
    cerr << "Error: can't open file: " << kmer_freq_file << ", please check the input file!" << endl;
    exit(-1);
  }

  for(uint64_t i = 0; i < BlockSizeVec.size(); i++) {
    uint64_t BlockStartPos = i * SrcBlockSize;
    infile.read((char *)buffer_compress, BlockSizeVec[i]);
    // uint64_t uncompLen = SrcBlockSize;
    uLongf uncompLen = SrcBlockSize;
    uncompress((Bytef *)buffer_uncompress, &uncompLen, (const Bytef *)buffer_compress, BlockSizeVec[i]);

    for(uint64_t j = 0; j < uncompLen; j++) {
      uint64_t idx = BlockStartPos + j;

      if(buffer_uncompress[j] > 0) {
        num_total_kmers += buffer_uncompress[j];
        num_effect_kmers++;
        freqy[buffer_uncompress[j]]++;

        if(buffer_uncompress[j] > freq_low) {
          if(buffer_uncompress[j] >=  freq_high) { //repeat kmer
            repeat_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll3[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll3[rc_idx % 4];
          }

          else if(buffer_uncompress[j] > freq_normal) { //normal kmer
            normal_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll2[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll2[rc_idx % 4];
          }

          else { //suspicious kmer
            suspicious_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll1[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll1[rc_idx % 4];
          }
        }

        else {  //error kmer
          error_freq_num++;
        }
      }
    }
  }

  cerr << "kmer freq table:\n";
  double error_freq_ratio = (double)error_freq_num / (double)num_effect_kmers;
  double suspicious_freq_ratio = (double)suspicious_freq_num / (double)num_effect_kmers;
  double normal_freq_ratio = (double)normal_freq_num / (double)num_effect_kmers;
  double repeat_freq_ratio = (double)repeat_freq_num / (double)num_effect_kmers;
  cerr << "kmer_num\tnode_num\terror_freq_num\terror_ratio\tsuspicious_freq_num\tsuspicious_ratio\tnormal_freq_num\tnormal_ratio\trepeat_freq_num\trepeat_ratio\n"
       << num_total_kmers << "\t" << num_effect_kmers << "\t" << error_freq_num << "\t" << error_freq_ratio << "\t" << suspicious_freq_num << "\t"
       << suspicious_freq_ratio << "\t" << normal_freq_num << "\t" << normal_freq_ratio << "\t" << repeat_freq_num << "\t" << repeat_freq_ratio << endl;
  delete []buffer_compress;
  delete []buffer_uncompress;
  BlockSizeVec.clear();
  infile.close();
}

//read the kmers from 8bit gz file into memory, store frequency value in 1 bit, range 0-1, need one file *.freq.gz
void load_kmerfreq(string & kmer_freq_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high)
{
  cerr << "loading kmerfreq table...\n";
  //read the frequence value from gz file
  uint64_t bufSize = 10000000;
  uint8_t * buf = new  uint8_t[bufSize];
  gzFile inFile;
  inFile = gzopen(kmer_freq_file.c_str(), "rb");
  int freqy[256];

  for(int i = 0; i < 256; i++) {
    freqy[i] = 0;
  }

  uint64_t error_freq_num = 0, suspicious_freq_num = 0, num_total_kmers = 0, normal_freq_num = 0, repeat_freq_num = 0, num_effect_kmers = 0, hig_freq_num = 0;
  uint64_t count_num = 0;
  uint64_t i = 0;
  cerr << "start loading kmer freq: " << kmer_freq_file << endl;

  do {
    count_num = gzread(inFile, (voidp)buf, bufSize);

    for(uint64_t j = 0; j < count_num; j++) {
      uint64_t idx = i + j;

      if(buf[j] > 0) {
        num_total_kmers += buf[j];
        num_effect_kmers++;
        freqy[buf[j]]++;

        if(buf[j] > freq_low) {
          if(buf[j] >= freq_high) { //repeat kmer
            repeat_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll3[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll3[rc_idx % 4];
          }

          else if(buf[j] > freq_normal) { //normal kmer
            normal_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll2[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll2[rc_idx % 4];
          }

          else { //suspicious kmer
            suspicious_freq_num++;
            freq[(idx * 2 + 1) / 8] |= bitAll1[idx % 4];
            uint64_t rc_idx = get_rev_com_kbit(idx, Kmer);
            freq[(rc_idx * 2 + 1) / 8] |= bitAll1[rc_idx % 4];
          }
        }

        else { //error kmer
          error_freq_num++;
        }
      }
    }

    i += count_num;
  }
  while(count_num == bufSize);

  cerr << "kmer freq table: \n";
  delete []buf;
  double error_freq_ratio = (double)error_freq_num / (double)num_effect_kmers;
  double suspicious_freq_ratio = (double)suspicious_freq_num / (double)num_effect_kmers;
  double normal_freq_ratio = (double)normal_freq_num / (double)num_effect_kmers;
  double repeat_freq_ratio = (double)repeat_freq_num / (double)num_effect_kmers;
  cerr << "kmer_num\tnode_num\terror_freq_num\terror_ratio\tsuspicious_freq_num\tsuspicious_ratio\tnormal_freq_num\tnormal_ratio\trepeat_freq_num\trepeat_ratio\n"
       << num_total_kmers << "\t" << num_effect_kmers << "\t" << error_freq_num << "\t" << error_freq_ratio << "\t" << suspicious_freq_num << "\t"
       << suspicious_freq_ratio << "\t" << normal_freq_num << "\t" << normal_freq_ratio << "\t" << repeat_freq_num << "\t" << repeat_freq_ratio << endl;
  gzclose(inFile);
}

void get_scaffold_kmer(string & read1, int & pos1, string & read2, int & pos2)
{
  //get kmer1
  pos1 = -1;
  int flag1 = 0;
  int cont_num1 = 0;
  int total_kmers1 = int(read1.size()) - Kmer + 1;
  int type1[total_kmers1];
  uint64_t kbit1[total_kmers1];
  string kseq1;

  for(int i = 0; i < total_kmers1; i++) {
    kseq1 = read1.substr(i, Kmer);
    kbit1[i] = seq2bit(kseq1);
    type1[i] = get_freq2(freq, kbit1[i]); //0 1 2 3
  }

  //get kmer2
  pos2 = -1;
  int flag2 = 0;
  int cont_num2 = 0;
  int total_kmers2 = int(read2.size()) - Kmer + 1;
  int type2[total_kmers2];
  uint64_t kbit2[total_kmers2];
  string kseq2;

  for(int i = 0; i < total_kmers2; i++) {
    kseq2 = read2.substr(i, Kmer);
    kbit2[i] = seq2bit(kseq2);
    type2[i] = get_freq2(freq, kbit2[i]); //0 1 2 3
  }

  vector<int> pos1_tem;

  for(int i = total_kmers1 - 1; i >= 0; i--) {
    if(type1[i] <= 1) {
      cont_num1 = 0;
    }

    else {
      cont_num1++;
    }

    if(cont_num1 >= Kmer) { //0:error-kmer 1:suspicious-kmer 2:normal-kmer 3:repeat-kmer
      pos1_tem.push_back(i);
    }
  }

  vector<int> pos2_tem;

  for(int i = 0; i < total_kmers2; i++) {
    if(type2[i] <= 1) {
      cont_num2 = 0;
    }

    else {
      cont_num2++;
    }

    if(cont_num2 >= Kmer) { //be sure kmer1 not equal to kmer2
      pos2_tem.push_back(i);
    }
  }

  int flag = 0;

  for(int i = 0; i < pos1_tem.size(); i++) {
    if(type1[pos1_tem[i] + Kmer / 2] == 2) {
      for(int j = 0; j < pos2_tem.size(); j++) {
        if(kbit1[pos1_tem[i] + Kmer / 2] != kbit2[pos2_tem[j] - Kmer / 2]) {
          pos1 = pos1_tem[i] + Kmer / 2;
          pos2 = pos2_tem[j] - Kmer / 2;
          flag = 1;
          break;
        }
      }
    }

    else {
      for(int j = 0; j < pos2_tem.size(); j++) {
        if(type2[pos2_tem[j] - Kmer / 2] == 2 && kbit1[pos1_tem[i] + Kmer / 2] != kbit2[pos2_tem[j] - Kmer / 2]) {
          pos1 = pos1_tem[i] + Kmer / 2;
          pos2 = pos2_tem[j] - Kmer / 2;
          flag = 1;
          break;
        }
      }
    }

    if(flag) {
      break;
    }
  }
}

void itostr(char * buf, int value)
{
  sprintf(buf, "%d", value);
}


void update_id(string & id, int pos)
{
  char posStr[100];
  itostr(posStr, pos);
  id = id + " " + posStr;
}

int get_max(double * score, int num, double cutoff, double DIF, vector<int> & maxs)
{
  int max_i = 0;
  double max = score[max_i];

  for(int i = 1; i < num; i++) {
    if(score[i] >= max) {
      max_i = i;
      max = score[i];
    }
  }

  if(max < cutoff) {
    return -1;
  }

  double low_cut = max * DIF;

  for(int i = 0; i < num; i++) {
    if(i == max_i) {
      continue;
    }

    if(score[i] >= low_cut && score[i] >= cutoff) {
      maxs.push_back(i);
    }
  }

  if(max == 0) {
    max_i = -1;
  }

  return max_i;
}

//for allale combination
double get_freq_score(double score, string & seq, string & qual, bool & have_high) //for base combination.
{
  int total_kmers = seq.size() - Kmer + 1;
  int type[total_kmers];
  string kseq;
  uint64_t kbit;

  for(int i = 0; i < total_kmers; i++) {
    kseq = seq.substr(i, Kmer);
    kbit = seq2bit(kseq);
    type[i] = get_freq2(freq, kbit); //0 or 1 or 2 or 3
  }

  double ratio, s = 1.0;
  double prob_score;
  int num = 0, high_num = 0, repeat_num = 0;

  for(int i = 0; i < total_kmers; i++) {
    if(type[i] >= 2) { //normal and repeat kmer
      high_num++;

      if(type[i] > 2) {
        repeat_num++;
      }
    }

    num++;
  }

  ratio =  double(high_num) / double(num);
  have_high = 0;

  if(Repeat_num_cutoff >= 0 && repeat_num > Repeat_num_cutoff) {
    have_high = 1;
  }

  s = ratio * score;
  return s;
}

double get_freq_score(int readlen1, int over, double score, string & seq, string & qual, bool & cross_check, int & full_type)
{
  int total_kmers = seq.size() - Kmer + 1;
  int type[total_kmers];
  cross_check = 0;
  string kseq;
  uint64_t kbit;

  for(int i = 0; i < total_kmers; i++) {
    kseq = seq.substr(i, Kmer);
    kbit = seq2bit(kseq);
    type[i] = get_freq2(freq, kbit); //0 1 2 3
  }

  int start = readlen1 - Kmer;
  int end = readlen1 - over;
  bool before_have_low = 0, after_have_low = 0, cross_have_low = 0, cross_have_error = 0;

  for(int i = end - Kmer; i <= start; i++) {
    if(i >= 0 && type[i] <= 1) {
      before_have_low = 1;
      break;
    }
  }

  for(int i = end; i <= end + Kmer; i++) {
    if(i < total_kmers && type[i] <= 1) {
      after_have_low = 1;
      break;
    }
  }

  if(!before_have_low && !after_have_low) {
    full_type = 2;

    if(over < Kmer - 1) {
      cross_check = 1;

      for(int i = start; i < end; i++) {
        if(i >= total_kmers) {
          break;
        }

        if(type[i] == 0) {
          cross_have_error = 1;
          break;
        }

        if(type[i] == 1) {
          cross_have_low = 1;
        }
      }

      if(cross_have_error) {
        return 0;
      }

      else if(cross_have_low) {
        return score;
      }

      else {
        return 1;
      }
    }
  }

  else if(before_have_low && after_have_low) {
    full_type = 0;
  }

  else {
    full_type = 1;
  }

  return score;
}


//used to get best allele combination.
int get_freq_max(double * score, int num, vector<int> & maxs, int old_max_i, string * seqs, string * quals)
{
  double new_score;
  double max = 0, max2 = 0;;
  double max_i = 0, max2_i = 0;
  bool have_high_freq_region = 0;

  //get max score for sub-optimals.
  for(int i = 0; i < maxs.size(); i++) {
    new_score = get_freq_score(score[ maxs[i] ], seqs[ maxs[i] ], quals[ maxs[i] ], have_high_freq_region);

    if(have_high_freq_region) {
      return -2;
    }

    if(new_score > max) {
      max2 = max;
      max2_i = max_i;
      max_i = maxs[i];
      max = new_score;
    }

    else if(new_score > max2) {
      max2 = new_score;
      max2_i = maxs[i];
    }
  }

  //get update score for optimal
  new_score = get_freq_score(score[old_max_i], seqs[old_max_i], quals[old_max_i], have_high_freq_region);

  if(have_high_freq_region) {
    return -2;
  }

  if(new_score > max) {
    return old_max_i;
  }

  else if(new_score == max) { //changed here.
    return -1;
  }

  //here not further use sub optimals.
  return (int)max_i;
}
bool no_other_mis(double * value, vector<int> & maxs, int max_min_num, int max_i, int uneque)
{
  bool no_other = 1;
  int nownum = 0;
  int len_cutoff = max_i;

  if(uneque != -1) {
    len_cutoff = Kmer - 2;
  }

  for(int i = 0; i < maxs.size(); i++) {
    if(maxs[i] == uneque || maxs[i] == max_i) {
      continue;
    }

    nownum = int(maxs[i] * (1 - value[maxs[i]]) + 0.1);

    if(nownum <= max_min_num && maxs[i] > len_cutoff) {
      return 0;
    }
  }

  return no_other;
}

int get_freq_max2(int readlen1, double * value, double * score, int num, vector<int> & maxs, int old_max_i, string * seqs, string * quals)
{
  double new_score;
  double max = 0, max2 = 0;;
  int max_i = -1, max2_i = -1;
  int full_num = 0, full_type = 0, double_full = -1, full_len_type = 0; //full_type: 0 no all high, 1 one all high, 2 two all high.
  bool cross_check = 0;
  map<int, bool> len_check;
  map<int, int> len_type;

  for(int i = 0; i < maxs.size(); i++) {
    cross_check = 0;
    full_type = 0;
    new_score = get_freq_score(readlen1, maxs[i], score[ maxs[i] ], seqs[ maxs[i] ], quals[ maxs[i] ], cross_check, full_type);
    len_check[ maxs[i] ] = cross_check;
    len_type[ maxs[i] ] = full_type;

    if(full_type == 2) {
      full_len_type++;

      if(new_score > max2) {
        max2_i = maxs[i];
        max2 = new_score;
      }
    }

    if(new_score >= max) { //the longest.
      //max2 = max;
      //max2_i = max_i;
      max_i = maxs[i];
      max = new_score;
    }

    if(new_score == 1) {
      full_num++;

      if(full_type == 2) {
        if(double_full == -1) {
          double_full = maxs[i];
        }

        else {
          double_full = -2;
        }
      }
    }

    if(TESTING) {
      CERRINF << " " << maxs[i] << ":" << new_score << ":" << cross_check << ":" << full_type;
    }
  }

  cross_check = 0;
  full_type = 0;
  new_score = get_freq_score(readlen1, old_max_i, score[old_max_i], seqs[old_max_i], quals[old_max_i], cross_check, full_type);

  if(full_type == 2 && new_score != 0) {
    full_len_type++;
  }

  len_check[ old_max_i ] = cross_check;
  len_type[ old_max_i ] = full_type;

  if(TESTING) {
    CERRINF << " " << old_max_i << ":" << new_score << ":" << cross_check << ":" << full_type << endl;
  }

  if(new_score == 1) {
    full_num++;
  }

  else if(new_score == 0 && full_len_type > 1 &&  max2_i == max_i) {
    if(max_i < Kmer - 1) { //old max i is lost and multi full type, and now max i is short.
      int max_mis_num = int(max_i * (1 - value[max_i]) + 0.1);

      for(int i = 0; i < maxs.size(); i++) {
        if(len_type[maxs[i]] == 2 && maxs[i] >= Kmer - 1) { //another is full check and longer and with less mismatch number.
          if(int(maxs[i] * (1 - value[maxs[i]]) + 0.1) <= max_mis_num) {
            return maxs[i];
          }
        }
      }
    }
  }

  if(max2 != 0 && max != 0 && max2 < 1 && max < 1 && max2_i != max_i) { //both have mismatch and full type is does not have the highest score.
    if(int(max2_i * (1 - value[max2_i]) + 0.1) <= int(max_i  * (1 - value[max_i]) + 0.1) && max2_i > max_i) { //onely have one mismatch
      max_i = max2_i;
      max = max2;

      if(TESTING) {
        CERRINF << ":c";
      }
    }

    else {
      if(TESTING) {
        CERRINF << ":f";
      }
    }
  }

  else {
    if(TESTING) {
      CERRINF << ":s";
    }
  }

  if(full_num > 1) {
    if(full_num <= 2 && double_full == -1 && full_type == 2 && new_score == 1) {
      return old_max_i;
    }

    if(full_num <= 2 && double_full > 0 && (full_type != 2 || new_score != 1)) {
      return double_full;
    }

    return -3; //same score.
  }

  if(old_max_i < Kmer - 1 && new_score == 1) {
    if(max == 0) { //1. max have been checked and not passed, so accept new score. 2. there is no other scores.
      if(cross_check || old_max_i >= Min) { //the longest passed check, then accept.
        return old_max_i;
      }

      else {
        return -100;  //too short and not cross checked, meaning there is mismatch near the region, in high risk, so not accept.
      }
    }

    if(!cross_check) { //maybe have sequencing error.
      int max_mis_num = int(max_i * (1 - value[max_i]) + 0.1);

      if(len_check[max_i]) { //max is short too and passed check, so max = 1. two short, one passed, the other not, then choose the passed one.
        return max_i;
      }

      else if(max_i < Kmer - 1) { //not pass either, so two short and neither passed, then can not choose.
        if(old_max_i > max_i && old_max_i >= Min && value[max_i] != 1 && len_type[old_max_i] >= len_type[max_i]) { //try to accept the full match short.
          return old_max_i;
        }

        return -1;
      }

      else if(len_type[max_i] == 2 && len_type[old_max_i] <= 1 &&
              no_other_mis(value, maxs, max_mis_num, max_i, old_max_i)) { //max_i is long, the two edge is high freq now, and have low for short overlap.
        //but in order to remove the short overlap, we need to make sure there are mismatch in the overlap region, or else we cannot remove it.
        return max_i;
      }//else if(max_i*(1-value[max_i]) <=1.1) //there are still mismatch outside the longer overlap region, and the mismatch num for long over is small.

      //	return max_i;
      return -1;
    }

    else { //maybe no sequencing error.
      if(len_check[max_i] || full_len_type > 2) { //multi pass check, at least three full length type.
        return -3; //same score.
      }

      else if(max_i < Kmer - 1) { //|| len_type[max_i] < 2) //max not pass , but only when both of the two are short, then the old max i can be accept.
        return old_max_i;
      }

      else if(new_score > max && max < new_score * DIF) { //max longer
        return old_max_i;
      }

      else {
        return -1;
      }
    }
  }

  //now there is no confident evidence for the old_max_i.
  //>>>add at 3.13

  //to accept long high score.
  if(new_score == 1) { //have high score for old_max_i, old_max_i >= Kmer -1 and max < 1.
    if(len_type[old_max_i] == 2 && (max_i < Kmer - 1 || len_type[max_i] < 2)) { //old max i is confident with high qualtiy and max i is weak.
      return old_max_i;
    }

    else if(max < DIF) {
      return old_max_i;
    }

    else if(len_type[old_max_i] < 2 && len_type[max_i] == 2 && maxs.size() == 1 && max_i >= Kmer - 1 && old_max_i < Min && //very very carefull.
            value[max_i] > 0.9) { //strict the total mismatch number.
      return max_i;
    }

    else {
      return -1;  //two long over, and both full suport,just one have no mismatch, and other have accepted mismatch.
    }
  }

  //to accept long or get the wrong one.
  if(old_max_i < Kmer - 1) { //new_score < 1, have mismatch here, so old_max_i is quite weak.
    int max_mis_num = int(max_i * (1 - value[max_i]) + 0.1);

    if(len_type[old_max_i] == 1 && max_i >= Kmer - 1 && len_type[max_i] == 2 &&
        //(max_i * (1-value[max_i]) <= old_max_i * (1-value[old_max_i]) + 0.5 || len_type[max_i] == 2)) //accept long max.
        (no_other_mis(value, maxs, max_mis_num, max_i, old_max_i) || full_len_type == 1))
      //the edge have mismatch, the long max can solve more mismatch and the oldmaxi is quite possibly wrong.
    {
      return max_i;
    }

    if(new_score == 0 && max > 0) { //the old_max_i must be wrong.
      if(max_i < Min && len_type[max_i] < 2 && (max != 1 || !no_other_mis(value, maxs, max_mis_num, max_i, -1))) { //not check for the short or there are mismatch for the short.
        return -100;
      }

      if(value[old_max_i] == 1 && max_mis_num > 0 && !no_other_mis(value, maxs, max_mis_num, max_i, -1)) {
        return -1;
      }

      return max_i;
    }

    if(len_type[old_max_i] == 1 && len_check[max_i] && max == 1 && full_len_type == 1 &&
        no_other_mis(value, maxs, max_mis_num, max_i, old_max_i))  //repeat may cause the real length is neither oldmax nor max, so to accept a short one is difficult.
      //max i is short, but to accept it , we still need to make sure this region is not related with repeat, then full len type is necessary.
    {
      return max_i;
    }
  }

  else if(max_i < Kmer - 1 && max < 1) {
    //old max i is long and have mismatch. drop the possible weak max with short over and mismatch.
    if(!len_check[max_i]) { //have base error here.
      int max_mis_num = int(old_max_i * (1 - value[old_max_i]) + 0.1);

      if(len_type[max_i] == 1 && len_type[old_max_i] == 2 && full_len_type == 1 && value[old_max_i] > 0.9 //define max mis match number here, attention.
          && no_other_mis(value, maxs, max_mis_num, old_max_i, max_i)) { //&& old_max_i *(1-value[old_max_i])<= max_i * (1-value[max_i]) + 0.5 ) //max i is possibly wrong.
        return old_max_i;
      }
    }

    else if(max == 0) { //max i is wrong.
      return old_max_i;
    }
  }

  if(full_len_type == 0 && max > 0 && max != 1 && new_score > 0 && new_score != 1) { //no type and not check, so freq is useless, then the short over will not quite effect the long over.
    //when the longer overlap have fewer mismatch and there score are quite similar, and the shorter one is too small and the longer one has under considered length, then the longer one can be accepted.
    if(max_i >= Kmer - 1 && old_max_i < Min) { //the length distance between the two should be large enough.
      int max_mis_num = int(max_i * (1 - value[max_i]) + 0.1);

      if(max_i * (1 - value[max_i]) <= old_max_i * (1 - value[old_max_i]) + 0.5 && max >= new_score * DIF &&
          no_other_mis(value, maxs, max_mis_num, max_i, -1)) { //more strict cut off.
        return max_i;
      }
    }

    else if(old_max_i >= Kmer - 1 && max_i < Min) {
      int max_mis_num = int(old_max_i * (1 - value[old_max_i]) + 0.1);

      if(old_max_i * (1 - value[old_max_i]) <= max_i * (1 - value[max_i]) + 0.5 && new_score >= max * DIF &&
          no_other_mis(value, maxs, max_mis_num, old_max_i, -1)) { //more strict cut off.
        return old_max_i;
      }
    }
  }

  //<<<end at 3.13

  if(new_score == 0 && max == 0) {
    return -100;
  }

  if(new_score < 1 && max == 0 && old_max_i < Min && cross_check == 0) { //small over with mismatch and not used cross check.
    return -100;
  }

  if(new_score > max && max < new_score * DIF) {
    return old_max_i;
  }

  else if(new_score == max) { //changed here.
    return -3;  //same score.
  }

  if(max * DIF > max2 &&  max * DIF > new_score) {
    return max_i;
  }

  return -1;
}

//get max combination type for each overlap length.
bool get_max_set2(double & ratio, int over, string & r1, string & r2, string & q1, string & q2, string & seq, string & qual)
{
  string before = r1.substr(0, r1.size() - over);
  string over1 = r1.substr(r1.size() - over);
  string after = r2.substr(over);
  string over2 = r2.substr(0, over);
  string beforeq, over1q, over2q, afterq;

  if(q1.size() > 0) {
    beforeq = q1.substr(0, q1.size() - over);
    over1q = q1.substr(q1.size() - over);
    afterq = q2.substr(over);
    over2q = q2.substr(0, over);
  }

  int unsame_num = 0;

  if(TESTING) {
    CERRINF << over << " " << ratio;
  }

  for(int i = 0; i < over; i++) {
    if(over1[i] != over2[i]) {
      if(TESTING) {
        CERRINF << " " << over1[i] << "<->" << over2[i];
      }

      unsame_num++;
    }
  }

  int unsame[unsame_num];
  char base_type[2][unsame_num];
  char qual_type[2][unsame_num];
  unsame_num = 0;

  for(int i = 0; i < over; i++) {
    if(over1[i] != over2[i]) {
      unsame[unsame_num] = before.size() + i;
      base_type[0][unsame_num] = r1[before.size() + i];
      base_type[1][unsame_num] = r2[i];

      if(q1.size() > 0) {
        qual_type[0][unsame_num] = q1[before.size() + i];
        qual_type[1][unsame_num] = q2[i];
      }

      unsame_num++;
    }
  }

  //cerr<<"unsame number is "<< unsame_num<<endl;
  seq = r1 + after;

  if(q1.size() > 0) {
    qual = q1 + afterq;
  }

  int start = int(r1.size()) - over - Kmer;
  int len = over + 2 * Kmer;

  if(start < 0) {
    start = 0;
    len = r1.size() + r2.size() - over;
  }

  if(r2.size() < over + Kmer) {
    len = r1.size() + r2.size() - over - start;
  }

  int num = (int)pow(2.0, double(unsame_num));
  double * score = new double[num];
  string seqs[num];
  string quals[num];
  //set all seq and score.
  int p;

  for(int i = 0; i < num; i++) {
    double prob_score = 1;

    for(int j = 0; j < unsame_num; j++) {
      seq[ unsame[j] ] = base_type[ matrix[i][j] ][j];

      if(qual.size() > 0) {
        qual[ unsame[j] ] = qual_type[ matrix[i][j] ][j];
        double e1 = pow(10, -double((qual_type[ matrix[i][j] ][j] - phred)) / 10.0); //10^(-i/10)
        double e2 = pow(10, -double((qual_type[1 - matrix[i][j] ][j] - phred)) / 10.0);
        prob_score *= 1 - e1 * (1 - e2); //probability of right choice.
      }
    }

    seqs[i] = seq.substr(start, len);

    if(qual.size() > 0) {
      quals[i] = qual.substr(start, len);
    }

    score[i] = prob_score;

    if(TESTING) {
      CERRINF << " " << i << "->" << score[i];
    }
  }

  //get the best choice.
  double cutoff = 0;
  vector<int> maxs;
  int max_type = get_max(score, num, cutoff, DIF, maxs);

  if(TESTING) {
    CERRINF << "\t" << max_type << "\t" << maxs.size();
  }

  if(max_type >= 0) {
    if(TESTING) {
      CERRINF << "\t" << score[max_type];
    }

    if(maxs.size() > 0 && (table_file.size() > 0 || table_file_by_jellyfish.size() > 0) && (mode == 1 || mode == 2 || mode == 3)) { //when there is freq information, we recalculate the probability.
      max_type = get_freq_max(score, num, maxs, max_type, seqs, quals);

      if(max_type < 0) {
        delete []score;
        return 1;
      }
    }

    ratio = score[max_type];
    seq = seq.substr(0, start) + seqs[max_type] + seq.substr(start + len);

    if(qual.size() > 0) {
      qual = qual.substr(0, start) + quals[max_type] + qual.substr(start + len);
    }
  }

  else {
    //same score.
    delete []score;
    return 1;
  }

  delete []score;
  return 0;
}


bool connect_overlap(double & ratio, int over, string & r1, string & r2, string & q1, string & q2, string & seq, string & qual)
{
  if(q1.size() == 0 && (table_file.size() == 0 || table_file_by_jellyfish.size() > 0)) {
    seq = r1 + r2.substr(over);
    return 0;
  }

  bool same_score = get_max_set2(ratio, over, r1, r2, q1, q2, seq, qual); //high quality and high freq.

  if(same_score) { //fastq format, all high quality base.
    qual = q1.substr(0, r1.size() - over) + q2;
    string overlap_part = r2.substr(0, over);

    for(int i = 0; i < over; i++) {
      string::size_type p = r1.size() - over + i;

      if(q1[p] > q2[i]) {
        overlap_part[i] = r1[p];
        qual[p] = q1[p];
      }

      else {
        overlap_part[i] = r2[i];
      }
    }

    seq = r1.substr(0, r1.size() - over) + overlap_part + r2.substr(over);
    return 1;
  }

  else {
    return 0;
  }

  return 0;
}

int get_max(double * value, int size, double value_cut_off, double diff_cut_off, bool flag, bool use_max)
{
  double max_value = 0.00;
  int max_len = 0;
  int max2_len = 0;
  double max2_value = -1.0;
  int match = 0;

  for(int i = 0; i < size; i++) {
    if(value[i] == 0) {
      continue;
    }

    if(max_value < value[i]) {
      max2_value = max_value;
      max2_len = max_len;
      max_value = value[i];
      max_len = i;
      continue;
    }

    if(max2_value < value[i] && max_value > value[i]) {
      max2_value = value[i];
      max2_len = i;
    }

    if(max_value == value[i]) {
      max2_value = max_value;
      max2_len = max_len;
      max_value = value[i];
      max_len = i;
    }
  }

  if(use_max) {
    if(max_value == max2_value) {
      same_score_num++;
    }

    return max_len;
  }

  if(max_value >= value_cut_off) { //kgf1.2
    if(max2_value < double(max_value) * diff_cut_off) { //for six bp, accept one mismatch.
      if(flag && max2_value < value_cut_off) {
        return max_len;
      }

      else if(!flag) {
        return max_len;
      }
    }

    return -1;
  }

  return -2;
}

bool get_scores(string & pid, string & read_a, string & read_b, string & q1, string & q2, int start, double * value, double * score, string * seqs, string * quals, vector<int> & pos, int & match_num, int & full_num)
{
  int len = read_a.size();
  int sum, mismatch;
  bool fail, score_is_same = 0;

  for(int i = start; i <= Max; i++) {
    if(i > read_a.length() || i > read_b.length()) {
      break;
    }

    string sa = read_a.substr(len - i);
    string sb = read_b.substr(0, i);
    string::size_type a_len = sa.size();
    sum = 0;
    mismatch = 0;
    fail = 0;

    for(int j = 0; j < a_len; j++) {
      if(sa[j] == 'N' || sb[j] == 'N') {
        continue;
      }

      if(sa[j] == sb[j]) {
        sum++;
      }

      else {
        mismatch++;
      }

      if(mismatch > mis_array[i] || mismatch >= MaxMismatch) {
        fail = 1;
        break;
      }
    }

    //if(fail)
    if(fail || mismatch >= MaxMismatch) {
      continue;
    }

    match_num++;
    value[i] = double(sum) / double(a_len);

    if(value[i] == 1) {
      full_num++;  //keep the full match length.
    }

    pos.push_back(i);
    score[i] = value[i];
    bool same_score = connect_overlap(score[i], i, read_a, read_b, q1, q2, seqs[i], quals[i]);

    if(same_score) {
      score_is_same = 1;
    }

    if(TESTING) {
      SCOREOUT << pid << "\t" << int(read_a.size() + read_b.size()) - real_insert_size << "\t" << i << "\t" << mismatch << "\t" << score[i] << endl;
    }

    if(TESTING) {
      CERRINF << endl;
    }
  }

  return score_is_same;
}

int get_real_len(int readlen1, double * value, double * score, string * seqs, string * quals, int num, int full_num, bool have_qual)
{
  int over_len;
  double base_value, diff;

  //get cut off ratio.
  if((table_file.size() > 0 || table_file_by_jellyfish.size() > 0) && (mode == 1 ||  mode == 2 || mode == 3)) { //have freq information.
    if(have_qual) { //fq.
      base_value = 0;
    }

    else { //fa.
      base_value = 0.01;
    }

    diff = DIF;//*1.3;
  }

  else {
    if(have_qual) {
      base_value = 0.001; //only average qualtiy.
    }

    else {
      base_value = MATCH_CUTOFF;
    }

    diff = DIF;
  }

  //get the best choice.
  vector<int> maxs;
  int max_type = get_max(score, num, base_value, diff, maxs);

  if(max_type > 0) {
    if((table_file.size() > 0 || table_file_by_jellyfish.size() > 0) && (mode == 1 ||  mode == 2 || mode == 3)) { //when there is freq information, we recalculate the probability.
      max_type = get_freq_max2(readlen1, value, score, num, maxs, max_type, seqs, quals); //recalculate all the short

      if(max_type < 0) {
        return max_type;  //same score. -3
      }

      over_len = max_type;
    }

    else if(maxs.size() == 0) {  //have no full and have no other higher than diff.
      over_len = max_type;
    }

    else {
      return -1;  //diff filter.
    }
  }

  else {
    return -2;  //all low.
  }

  return over_len;
}

int align_connect(string & aid, string & read_a, string & q1, string & bid, string & read_b, string & q2, string & seq, string & qual)
{
  //keep the match value of each match type
  if(read_a.size() <= Min) {
    cerr << "read too short! " << read_a.size() << "\t" << Min << endl << read_a << endl;
    return 0;
  }

  int pos1 = -1, pos2 = -1, last_change = -1;
  string pid, raw_a(read_a), raw_b(read_b);

  //trainning...
  if(TESTING) {
    real_insert_size = get_insert(aid, pid);
  }

  pid = get_pair_id(aid, bid);
  //	get_scaffold_kmer(read_a, pos1, read_b, pos2);
  double * value = new double[Max + 1];
  double * score = new double[Max + 1];
  string * seqs = new string[Max + 1];
  string * quals = new string[Max + 1];
  vector<int> pos;

  for(int i = 0; i < Max + 1; i++) {
    value[i] = 0;
    score[i] = 0;
  }

  int start = Min, score_num = 0, match_num = 0, full_num = 0;

  if((table_file.size() > 0 || table_file_by_jellyfish.size() > 0) && (mode == 1 ||  mode == 2 || mode == 3)) {
    start = 3;
  }

  //get insert size;
  bool score_is_same = get_scores(pid, read_a, read_b, q1, q2, start, value, score, seqs, quals, pos, match_num, full_num);
  int over_len;

  if(match_num > 0 && !score_is_same) {
    over_len = get_real_len(read_a.size(), value, score, seqs, quals, Max + 1, full_num, q1.size() > 0);
  }

  else if(score_is_same) {
    over_len = -3;
  }

  else {
    over_len = -100;  //have no match.
  }

  if(TESTING) {
    CERRINF << pid << "\t" << real_insert_size << "\t" << match_num << "\t" << full_num << "\t" << over_len << "\t";
  }

  //connected.
  if(over_len >= start) {
    seq = seqs[over_len];

    if(q1.size() > 0) {
      qual = quals[over_len];
    }

    //trainning...
    if(TESTING && real_insert_size > 0 && real_insert_size != seq.size()) {
      if(real_insert_size <= (read_a.size() + read_b.size())) {
        FALSECON << pid << "\t" << over_len << "\t" << score[over_len] << "\t" << real_insert_size << "\t" << score[(read_a.size() + read_b.size()) - real_insert_size] << endl;
      }

      else {
        FALSECON << pid << "\t" << over_len << "\t" << score[over_len] << "\t" << real_insert_size << "\t0" << endl;
      }
    }

    delete []value;
    delete []score;
    delete []seqs;
    delete []quals;

    if(TESTING) {
      CERRINF << endl;
    }

	//update for 1.1.3
	aid = pid;
    return over_len;
  }

  delete []value;
  delete []score;
  delete []seqs;
  delete []quals;

  if((table_file.size() > 0 || table_file_by_jellyfish.size() > 0) && mode != 0) {
    get_scaffold_kmer(read_a, pos1, read_b, pos2);
    //unconnected.
    update_id(aid,  pos1); //for unconnect read pairs.
    update_id(bid,  pos2);
  }

  unconnect_count++;

  //output kmer.
  if(mode != 0 && pos1 > 0 && pos2 > 0) {
    kout << unconnect_count << " " << read_a.substr(pos1, Kmer) << " " << read_a.size() << " " << pos1 << " " << read_b.substr(pos2, Kmer) << " " << read_b.size() << " " << pos2 << endl;
  }

  switch(over_len) {
    case -1:
      if(TESTING) {
        CERRINF << "diff fail\n";
      }

      dif_filter_num++;
      break;

    case -2:
      if(TESTING) {
        CERRINF << "low filter\n";
      }

      low_filter_num++;
      break;

    case -3:
      if(TESTING) {
        CERRINF << "have same score\n";
      }

      same_score_num++;
      break;

    case -100:
      if(TESTING) {
        CERRINF << "no score\n";
      }

      no_score++;
      break;

    default:
      if(TESTING) {
        CERRINF << "small over\n";
      }

      small_over++;
      break;
  }

  return 0;
}

