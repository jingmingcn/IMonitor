/*
 * seqkmer.h
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
 * File          : seqkmer.h
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Wei Fan, Binghang Liu, Jianying Yuan
 *
 *============================================================================
 */


#ifndef __SEQ_KMER_H_
#define __SEQ_KMER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

using namespace std;
typedef unsigned long long uint64_t;
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;

extern char alphabet[128];
extern char alphabet2[128];
extern char bases[5];
extern char l_bases[5];
extern char c_bases[5];
extern char cl_bases[5];
extern uint8_t bitAll[8];
extern uint8_t bitAll1[4];
extern uint8_t bitAll2[4];
extern uint8_t bitAll3[4];
extern uint64_t bitLeft[4];

//get the frequence value for a given index
int get_freq(uint8_t * freq, uint64_t idx);

//get the frequence value for a given index
int get_freq2(uint8_t * freq, uint64_t idx);

//convert kmer-seq to kmer-bit
uint64_t seq2bit(string & kseq);

//convert kmer-bit to kmer-seq
string bit2seq(uint64_t kbit, int kmerSize);

//check whether a sequence contain non base characters, such as "N"
int check_seq(string & seq);

inline char complement_base(char base)
{
  return c_bases[alphabet[base]];
}
//get the reverse and complement sequence
void reverse_complement(string & in_str, string & out_str);

//read file_list into a vector
void reading_file_list(string & file_list, vector<string> & files);

//get the reverse and complement kbit, independent of the sequence
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize);

//display a number in bit format
void display_num_in_bits(uint64_t num, int len);

//chengfang calculation for small integers
uint64_t pow_integer(int base, int exponent);

//get the complement sequence
void complement_sequence(string & str);

string reverse_complement(string & seq);
string reverse_str(string & seq);

#endif
