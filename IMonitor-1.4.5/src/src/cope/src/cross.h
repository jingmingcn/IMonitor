/*
 * cross.h
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
 * File          : cross.h
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Binghang Liu, Jianying Yuan
 *
 *============================================================================
 */
 
#ifndef __CROSS__H_
#define __CROSS__H_

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "connect.h"
#include "hashSet.h"
#include "seqKmer.h"
#include "gzstream.h"


using namespace std;

extern uint8_t * freq;
extern HashSet * kmerHash;
extern uint64_t initial_size;
extern double load_factor;
extern string Cross_info_list_file;
extern int Batch_num;
extern int Stop_connect_cross_read;
extern int mode;
extern int cross_read_type;
extern int Kmer;
extern int RepeatCoverNum;
extern int MinCoverNum;
//if(TESTING)
extern ofstream FALSECON2;
extern bool TESTING;
extern ofstream CROSSCERRINF;

struct CrossInf {
  uint32_t p; //pair id.
  uint32_t i; //insert size
  char * seq; //insert seq.
  char * qual;
};

struct Pkmer {
  uint64_t k1: 38; //19mer.
  uint64_t k2: 38; //19mer.
  uint64_t f1: 1;
  uint64_t f2: 1;
  uint64_t p1: 8; //256 for 150bp reads.
  uint64_t p2: 8;
  uint64_t l1: 8; //<256 for read1 length
  uint64_t l2: 8; //<256 for read2 length
  uint64_t p: 32; //pkmer position.
  uint64_t f: 1;
};

struct KmerSite {
  uint64_t k: 38; //max k=19
  uint64_t p: 8; //7->8, max position in read is 256 for 150bp reads or connected 250bp reads.
  uint64_t i: 32; //25->32
  uint64_t l: 17; //4->17, smaller than 25 in hashset structure.
  uint64_t s: 1; //strand.

};

void cross_connect(string & fq1, string & fq2, string & connect_file, string & cross_read_list_file, string & kmer_file);
void cross_connect(string & fq1, string & fq2, string & a_file, string & b_file, string & cross_read_list_file, string & kmer_file);

#endif
