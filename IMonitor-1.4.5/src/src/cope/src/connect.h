/*
 * connect.h
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
 * File          : connect.h
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Binghang Liu, Jianying Yuan
 *
 *============================================================================
 */

#ifndef __CONNECT__H__
#define __CONNECT__H__

#include <string>
#include <cstdlib>
#include <zlib.h>
#include "seqKmer.h"
#include "gzstream.h"

extern uint64_t dif_filter_num;
extern uint64_t low_filter_num;
extern uint64_t connect_succeed;
extern uint64_t connect_fail;
extern uint64_t same_score_num;
extern uint64_t no_score;
extern uint64_t small_over;
extern uint64_t unconnect_count;

extern int Kmer;
extern uint8_t * freq;
extern int Max;
extern int Min;
extern float DIF;
extern float MATCH_CUTOFF;
extern int Repeat_num_cutoff;

extern int * mis_array;
extern bool ** matrix;
extern int MaxMismatch;
extern int phred;
extern int mode;

extern string table_file;
extern string table_file_by_jellyfish;
extern ogzstream kout;

//train.....
extern ofstream FALSECON;
extern ofstream CERRINF;
extern bool TESTING;
extern ofstream SCOREOUT;

void load_jellyfish_kmerfreq(string & kmer_freq_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high);

//load pallal compress kmer freq table.
void load_kmerfreq(string & kmer_freq_file, string & kmer_freq_len_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high);

//read the kmers from 8bit gz file into memory, store frequency value in 1 bit, range 0-1
void load_kmerfreq(string & kmer_freq_file, uint8_t * freq, int freq_low, int freq_normal, int freq_high);

void get_mis_table(double match_ratio, int * table, int len);

void str_split(string & line, int clim, vector<string> & svec);

int align_connect(string & aid, string & read_a, string & q1, string & bid, string & read_b, string & q2, string & seq, string & qual);

void get_max_freq_set(string & r1, string & r2, string & q1, string & q2, int insert_size, string & seq, string & qual); //get final seq and qual.

bool get_max_set(double & ratio, int over, string & r1, string & r2, string & q1, string & q2, string & seq, string & qual);

bool connect_overlap(double & ratio, int over, string & r1, string & r2, string & q1, string & q2, string & seq, string & qual);

int get_max(double * value, int size, double value_cut_off, double diff_cut_off, bool flag, bool use_max);

bool connect_pair(string & a_id, string & read_a, string & a_s, string & a_q, string & b_id, string & read_b, string & b_s, string & b_q);

void load_file(string & file1, string & file2);

#endif

