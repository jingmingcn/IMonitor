/*
 * seqkmer.cpp
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
 * File          : seqkmer.cpp
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Wei Fan, Binghang Liu, Jianying Yuan
 *
 *============================================================================
 */

#include "seqKmer.h"

char alphabet[128] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

char alphabet2[128] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//from 0 1 2 3 4 to A C G T N
char bases[5] = {
  'A', 'C', 'G', 'T', 'N'
};
char l_bases[5] = {
  'a', 'c', 'g', 't', 'n'
};

//from 0 1 2 3 4 to T G C A N
char c_bases[5] = {
  'T', 'G', 'C', 'A', 'N'
};
char cl_bases[5] = {
  't', 'g', 'c', 'a', 'n'
};

//used for bit-operation
uint8_t bitAll[8] = {128, 64, 32, 16, 8, 4, 2, 1};

uint8_t bitAll1[4] = {64, 16, 4, 1}; //1

uint8_t bitAll2[4] = {128, 32, 8, 2}; //2

uint8_t bitAll3[4] = {192, 48, 12, 3}; //3

uint64_t bitLeft[4] = {
  0x0000000000000000,
  0x4000000000000000,
  0x8000000000000000,
  0xc000000000000000
};

//get the frequence value for a given index
int get_freq(uint8_t * freq, uint64_t idx)
{
  //	cerr<<"9"<<"\tidx:"<<idx<<endl;
  int value = (freq[idx / 8] >> (7 - idx % 8)) & 0x1u;
  return value;
}

//get the frequence value for a given index
int get_freq2(uint8_t * freq, uint64_t idx)
{
  //	cerr<<"9"<<"\tidx:"<<idx<<endl;
  int value = (freq[(idx * 2 + 1) / 8] >> ((3 - idx % 4) * 2)) & 0x3u;
  return value;
}

//convert kmer-seq to kmer-biT
uint64_t seq2bit(string & kseq)
{
  uint64_t kbit = 0;

  for(int i = 0; i < kseq.size(); i++) {
    kbit = (kbit << 2) | alphabet[kseq[i]];
  }

  return kbit;
}

//convert kmer-bit to kmer-seq
string bit2seq(uint64_t kbit, int kmerSize)
{
  string kseq;

  for(int i = 0; i < kmerSize; i++) {
    kseq.push_back(bases[(kbit >> (kmerSize - 1 - i) * 2) & 0x3]);
  }

  return kseq;
}

//check whether a sequence contain non base characters, such as "N"
int check_seq(string & seq)
{
  int is_good = 1;

  for(int i = 0; i < seq.size(); i++) {
    if(alphabet[seq[i]] == 4) {
      is_good = 0;
      break;
    }
  }

  return is_good;
}

//get the reverse and complement sequence
void reverse_complement(string & in_str, string & out_str)
{
  for(int64_t i = in_str.size() - 1; i >= 0; i--) {
    out_str.push_back(c_bases[alphabet[in_str[i]]]);
  }
}

//get the complement sequence
void complement_sequence(string & str)
{
  for(int64_t i = 0; i < str.size(); i++) {
    str[i] = c_bases[alphabet[str[i]]];
  }
}

string reverse_complement(string & seq)
{
  string nseq;

  for(int i = seq.size() - 1; i >= 0; i--) {
    nseq.push_back(c_bases[alphabet2[seq[i]]]);
  }

  return nseq;
}

string reverse_str(string & seq)
{
  if(seq.size() == 0) {
    return seq;
  }

  string rseq(seq);

  for(int i = seq.size() - 1; i >= 0; i--) {
    rseq[i] = seq[seq.size() - i - 1];
  }

  return rseq;
}

//get the reverse and complement kbit, independent of the sequence
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize)
{
  kbit = ~kbit;
  kbit = ((kbit & 0x3333333333333333LLU) <<  2) | ((kbit & 0xCCCCCCCCCCCCCCCCLLU) >>  2);
  kbit = ((kbit & 0x0F0F0F0F0F0F0F0FLLU) <<  4) | ((kbit & 0xF0F0F0F0F0F0F0F0LLU) >>  4);
  kbit = ((kbit & 0x00FF00FF00FF00FFLLU) <<  8) | ((kbit & 0xFF00FF00FF00FF00LLU) >>  8);
  kbit = ((kbit & 0x0000FFFF0000FFFFLLU) << 16) | ((kbit & 0xFFFF0000FFFF0000LLU) >> 16);
  kbit = ((kbit & 0x00000000FFFFFFFFLLU) << 32) | ((kbit & 0xFFFFFFFF00000000LLU) >> 32);
  return kbit >> (64 - (ksize << 1));
}

//read file_list into a vector
void reading_file_list(string & file_list, vector<string> & files)
{
  ifstream infile(file_list.c_str());

  if(! infile) {
    cerr << "fail to open input file" << file_list << endl;
  }

  string file_name;

  while(getline(infile, file_name, '\n')) {
    if(file_name.size()) {
      files.push_back(file_name);
    }
  }
}

//display a number in bit format
void display_num_in_bits(uint64_t num, int len)
{
  cout << num << "\t" << len << "\t";
  len -= 1;

  for(int i = 0; i <= len; i++) {
    int bit = (num >> (len - i)) & 1;
    cerr << bit;
  }

  cerr << endl;
}

//base^exponent calculation for small integers
uint64_t pow_integer(int base, int exponent)
{
  uint64_t result = 1;

  for(int i = 1; i <= exponent; i++) {
    result *= base;
  }

  return result;
}
