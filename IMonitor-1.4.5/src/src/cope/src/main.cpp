/*
 * main.cpp
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
 * File          : main.cpp
 * Revision      : 1.1.3
 * Revision_date : 2012/10/08
 * Author(s)     : Binghang Liu, Jianying YUan
 *
 * Update:
 * 	v1.1.3, changed the format of connected id.
 *
 *============================================================================
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <cstdlib>
#include <sstream>

#include "cross.h"
#include "connect.h"
#include "gzstream.h"

using namespace std;

#ifndef VID
#define VID "v1.1.3"
#endif

const char * PROGRAM = "COPE (Connecting Overlapped Pair-End reads)";
const char * AUTHOR  = "BGI shenzhen";
const char * VERSION = VID;
const char * CONTACT = "liubinghang@genomics.org.cn; yuanjianying@genomics.org.cn";
const char * command_line;

bool TESTING = false;

//basic parameters
int phred = 64;
int Kmer = 17;
int MaxOverLen = 200;
int Min = 10;
int Max = 70;
int Mean = 0;
uint64_t Pairs = 0;

//match ratio parameters
float  MATCH_CUTOFF = 0.75;
float  DIF = 0.7;
int DIFF = 5;

//quality parameters.
int Nfilter = 1;
float  B_CUTOFF = 0.9;

//cross parameters.
int RepeatCoverNum = 50;
int MinCoverNum = 5;

//mode parameters.
int mode = 3; //0 simple connect; 1 connect with kmer frequncy; 2 connect left reads from mode 1 with cross reads; 3 full mode
int Repeat_num_cutoff = -1; //if repeat kmer number in overlap region is larger than this value, not connected.
int Stop_connect_cross_read = 0;
int cross_read_type = 1; //1 use connected reads, 0 use raw reads.
int Batch_num = 1;

//build-in parameters
uint64_t initial_size = 1024 * 1024;
double load_factor = 0.75;
int MaxMismatch = 10;
int sum_num;

int Freq_low = 3;
int Freq_normal = 10;
int Freq_high = 60;

//stat parameters
uint64_t total_pairs = 0;
uint64_t connect_pairs = 0;
uint64_t low_quality_pairs = 0;
uint64_t dif_filter_num = 0;
uint64_t low_filter_num = 0;
uint64_t connect_succeed = 0;
uint64_t connect_fail = 0;
uint64_t same_score_num = 0;
uint64_t no_score = 0;
uint64_t small_over = 0;
uint64_t unconnect_count = 0;

//used parameters
HashSet * kmerHash;
uint8_t * freq;
bool ** matrix;
int * mis_array; //the mismatch threshold for each overlap length.

//file names.
string a_file;
string b_file;
string fq1;
string fq2;
string out_connect_file;
string table_file;
string table_len_file;
string table_file_by_jellyfish;
string Cross_read_list_file;

string Pair_kmer_file = "./pair.kmer.list.gz";
string Cross_info_list_file;

ofstream fout;
ofstream fout_fq1;
ofstream fout_fq2;
ogzstream kout;//output pair kmer.

//training....
ofstream FALSECON;
ofstream FALSECON2;
ofstream CERRINF;
ofstream SCOREOUT;
ofstream CROSSCERRINF;

void usage(void)
{
  cout << "\nProgram:\t" << PROGRAM << "\n"
       << "Version:\t" << VERSION << "\n"
       << "Author:\t\tBGI-ShenZhen\n"
       << "CompileDate:\t" << __DATE__ << " time: " << __TIME__ << "\n"
       << "Contact: \t" << CONTACT << "\n"
       << "Usage:\tcope [option]" << "\n"
       << "\t-a  <str>   query read1 file, (*.fq, *.fa, *.fq.gz, *.fa.gz)\n"
       << "\t-b  <str>   query read2 file\n"
       << "\t-o  <str>   output connected file in *.fq or *fa\n"
       << "\t-2  <str>   output fail connected read1.fq\n"
       << "\t-3  <str>   output fail connected read2.fq\n"
       << "\t-l  <int>   lower bound overlap length, 8 or 10 is preferred for 100bp read with 170~180 insert size, default=" << Min << "\n" //minimal insert size
       << "\t-u  <int>   higher bound overlap length, 70 is preferred for 100 bp reads, default=" << Max << "\n" //max_insert_size
       << "\t-c  <float>   match ratio cutoff, default=" << MATCH_CUTOFF << "\n" //mean insert size
       << "\t-d  <float>   match max2_ratio/max_ratio, important for mode 0, default=" << DIF << "\n" //overlap length
       << "\t-B  <float>   ratio cut off value of Base-quality=2, default=" << B_CUTOFF << "\n"
       << "\t-N  <int>   filter N base containing pairs(1) or not(0), default=" << Nfilter << "\n"
       << "\t-T  <int>   read pair number threshold for parameter training or testing: " << Pairs << "\n"
       << "\t-s  <int>   the smallest ASCII value used to represent base quality in fastq files. 33 for the later Illumina platforms and Sanger reads, and 64 for the earlier Illumina platforms. default=" << phred << "\n"
       << "\t-m  <int>   mode type: 0 simple connect; 1 k-mer frequency assisted connection; 2 Auxiliary reads and cross connection for left reads from mode 1; 3 full mode(including 1 and 2),default=" << mode << "\n"
       << "\n\n\twhen set mode large than 0, set the following options:\n"
       << "\t-k  <int>   k-mer size, default=" << Kmer << "\n"
       << "\t-t  <str>   compressed k-mer frequency table\n"
       << "\t-f  <str>   kmer frequency table len, for parallelly compressed k-mer frequency table\n"
       << "\t-j  <str>   compressed kmer frequency table counted by jellyfish\n"
       << "\t-L  <int>   set the threshold of freq-low(Kmer with lower frequency will not be considered for base selection), default=" << Freq_low << "\n"
       << "\t-M  <int>   set the threshold of freq-normal(Kmer whith lower frequency will not be considered for spanning k-mer selection), default=" << Freq_normal << "\n"
       << "\t-H  <int>   set the threshold of freq-high(Kmer whith higher frequency will not be considered for spanning k-mer selection, 2 times the average kmer frequency is preferred), default=" << Freq_high << "\n"
       << "\t-R  <int>   set the threshold of spanning k-mers number of extremely high frequency, threshold=3 is preferred for high repeat genome, set -1 for disable this threshold, default=" << Repeat_num_cutoff << "\n"
       << "\t-K  <str>   input pair-kmer file for mode=2, default=" << Pair_kmer_file << "\n"
       << "\n\twhen set mode=2 or 3, you also need to set the following options:\n"
       << "\t-D  <int>   set the batch number of cross-information file(cross.read.inf*.gz) output, default=" << Batch_num << "\n"
       << "\t-r  <int>   use connected reads(1) or raw reads(0) to find cross reads, default=" << cross_read_type << "\n"
       << "\t-F  <str>   input extra read-files list(each line for one read file) while need to find cross reads with the help of other read-files" << "\n"
       << "\n\t-h          print help information\n\n"
       << "Example:\n"
       << " ##simple connect mode: ./cope -a read1.fq -b read2.fq -o connect.fa -2 left1.fq -3 left2.fq -m 0 >cope.log 2>cope.error\n"
       << " ##k-mer frequency assisted connection: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 1 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error\n"
       << " ##Auxiliary reads and cross connection for left reads from mode 1: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 2 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error\n"
       << " ##connect by full mode: ./cope -a read1.fq -b read2.fq -o connect.fq -2 left1.fq -3 left2.fq -m 3 -t kmer_table.cz -f kmer_table.cz.len >cope.log 2>cope.error\n\n"
       << "\n";
  exit(1);
}

int mGetOptions(int rgc, char * rgv[])
{
  //[options]
  int i;

  for(i = 1; i < rgc; i++) {
    if(rgv[i][0] != '-') {
      return i;
    }

    switch(rgv[i][1]) {
      case 'a':
        a_file = rgv[++i];
        break;

      case 'b':
        b_file = rgv[++i];
        break;

      case 'o':
        out_connect_file = rgv[++i];
        break;

      case '2':
        fq1 = rgv[++i];
        break;

      case '3':
        fq2 = rgv[++i];
        break;

      case 't':
        table_file = rgv[++i];
        break;

      case 'f':
        table_len_file = rgv[++i];
        break;

      case 'j':
        table_file_by_jellyfish = rgv[++i];
        break;

      case 'H':
        Freq_high = atoi(rgv[++i]);
        break;

      case 'R':
        Repeat_num_cutoff = atoi(rgv[++i]);
        break;

      case 'k':
        Kmer = atoi(rgv[++i]);
        break;

      case 'l':
        Min = atoi(rgv[++i]);
        break;

      case 'u':
        Max = atoi(rgv[++i]);
        break;

      case 'c':
        MATCH_CUTOFF = atof(rgv[++i]);
        break;

      case 'd':
        DIF = atof(rgv[++i]);
        break;

      case 'B':
        B_CUTOFF = atof(rgv[++i]);
        break;

      case 'T':
        Pairs = atoi(rgv[++i]);
        break;

      case 'N':
        Nfilter = atoi(rgv[++i]);
        break;

      case 'L':
        Freq_low = atoi(rgv[++i]);
        break;

      case 'M':
        Freq_normal = atoi(rgv[++i]);
        break;

      case 's':
        phred = atoi(rgv[++i]);
        break;

      case 'r':
        cross_read_type = atoi(rgv[++i]);
        break;

      case 'm':
        mode = atoi(rgv[++i]);
        break;

      case 'D':
        Batch_num = atoi(rgv[++i]);
        break;

      case 'F':
        Cross_read_list_file = rgv[++i];
        break;

      case 'S':
        Stop_connect_cross_read = atoi(rgv[++i]);
        break;

      case 'K':
        Pair_kmer_file = rgv[++i];
        break;

      case 'C':
        Cross_info_list_file = rgv[++i];
        break;

      case 'h':
        usage();
    }
  }

  return i;
}

bool connect_pair(string & a_id, string & read_a, string & a_s, string & a_q, string & b_id, string & read_b, string & b_s, string & b_q)
{
  if(Nfilter == 1) {
    if(read_a.find("N") != string::npos || read_b.find("N") != string::npos) {
      low_quality_pairs++;
      unconnect_count++;

      if(a_q.size() > 0) {
        fout_fq1 << "@" << a_id << endl;
        fout_fq1 << read_a << endl;
        fout_fq1 << a_s << endl;
        fout_fq1 << a_q << endl;
        fout_fq2 << "@" << b_id << endl;
        fout_fq2 << read_b << endl;
        fout_fq2 << b_s << endl;
        fout_fq2 << b_q << endl;
      }

      else {
        fout_fq1 << ">" << a_id << endl;
        fout_fq1 << read_a << endl;
        fout_fq2 << ">" << b_id << endl;
        fout_fq2 << read_b << endl;
      }

      return 0;
    }
  }

  //filter by the quality value of read
  int ab_value = 0, bb_value = 0;

  if(a_q.size() > 0 && B_CUTOFF > 0.0) {
    int left_len = a_q.size() - Max - 1;

    if(left_len < 0) {
      left_len = 0;
    }

    for(int i = a_q.size() - 1; i >= left_len; i--) {
      if(a_q[i] <= 2 + phred) {
        ab_value ++ ;
      }
    }

    left_len = b_q.size() - Max - 1;

    if(left_len < 0) {
      left_len = 0;
    }

    for(int i = b_q.size() - 1; i >= left_len; i--) {
      if(b_q[i] <= 2 + phred) {
        bb_value ++ ;
      }
    }

    if(B_CUTOFF < double(ab_value) / double(Max) || B_CUTOFF < double(bb_value) / double(Max)) {
      //			cout<<"Low quality Read, fail to connect\n";
      unconnect_count++;
      low_quality_pairs++;
      fout_fq1 << "@" << a_id << endl;
      fout_fq1 << read_a << endl;
      fout_fq1 << a_s << endl;
      fout_fq1 << a_q << endl;
      fout_fq2 << "@" << b_id << endl;
      fout_fq2 << read_b << endl;
      fout_fq2 << b_s << endl;
      fout_fq2 << b_q << endl;
      return 0;
    }
  }

  //start connect
  string rread_b = reverse_complement(read_b);
  string rb_q = reverse_str(b_q);
  string seq, qual;
  int connected = align_connect(a_id, read_a, a_q, b_id, rread_b, rb_q, seq, qual); //update id.

  if(!connected) {
    read_b = reverse_complement(rread_b);

    if(a_q.size() != 0) {
      fout_fq1 << "@" << a_id << endl;
      fout_fq1 << read_a << endl;
      fout_fq1 << a_s << endl;
      fout_fq1 << a_q << endl;
      fout_fq2 << "@" << b_id << endl;
      fout_fq2 << read_b << endl;
      fout_fq2 << b_s << endl;
      fout_fq2 << b_q << endl;
    }

    else {
      fout_fq1 << ">" << a_id << endl;
      fout_fq1 << read_a << endl;
      fout_fq2 << ">" << b_id << endl;
      fout_fq2 << read_b << endl;
    }

    return 0;
  }

  if(a_q.size() != 0) {
    //fout << "@cp" << connect_pairs << "\t" << a_id << "_" << b_id << "\t" << connected << "\t" << read_b.size() + read_a.size() - connected << "\n"
	  fout << "@" << a_id << "\t" << read_b.size() + read_a.size() - connected  << "\n" //a_id has been updated to pid.
		   << seq << "\n"
           << "+\n"
           << qual << endl;
  }

  else {
    //fout << ">cp" << connect_pairs << "\t" << a_id << "_" << b_id << "\t" << connected << "\t" << read_b.size() + read_a.size() - connected << "\n"
	  fout << ">" << a_id << "\t" << read_b.size() + read_a.size() - connected  << "\n" //a_id has been updated to pid.
		   << seq << "\n";
  }

  return 1;
}

void load_file(string & file1, string & file2)
{
  igzstream fin_a(file1.c_str());
  igzstream fin_b(file2.c_str());

  if(!fin_a || !fin_b) {
    cerr << "please check file " << file1 << ". and file " << file2 << endl;
    exit(0);
  }

  fout_fq1.open(fq1.c_str());
  fout_fq2.open(fq2.c_str());
  fout.open(out_connect_file.c_str());

  if(mode != 0) {
    string kmer_file = "pair.kmer.list.gz";
    kout.open(kmer_file.c_str());
  }

  if(!fout_fq1 || !fout_fq2 || !fout || !kout) {
    cerr << "fail creat file" << endl;
    exit(1);
  }

  string read_a, read_b, a_id, b_id, a_s, b_s, a_q, b_q;
  cerr << "Begin read files and connect pairs...\n";

  while(getline(fin_a, read_a, '\n')) {
    string seq;

    //load id;
    if(read_a[0] == '@' || read_a[0] == '>') {
      a_id = read_a.substr(1);
    }

    getline(fin_b, read_b, '\n');

    if(read_b[0] == '@' || read_b[0] == '>') {
      b_id = read_b.substr(1);
    }

    //load strand and sequence quality
    if(read_a[0] == '@') {              //for *.fq
      getline(fin_a, read_a, '\n');
      getline(fin_b, read_b, '\n');
      getline(fin_a, a_s, '\n');
      getline(fin_a, a_q, '\n');
      getline(fin_b, b_s, '\n');
      getline(fin_b, b_q, '\n');
    }

    else {                             //for *.fa
      getline(fin_a, read_a, '\n');
      getline(fin_b, read_b, '\n');
    }

    total_pairs++;
    bool connected = connect_pair(a_id, read_a, a_s, a_q, b_id, read_b, b_s, b_q);

    if(connected) {
      connect_pairs++;
    }

    if(total_pairs % 1000000 == 0) {
      cerr << "Process pair reads number: " << total_pairs << endl;
    }

    if(Pairs > 0 && total_pairs == Pairs) {
      cout << "Connected pair reads number " << connect_pairs << endl;
      break;
    }
  }//end of while.

  fin_a.close();
  fin_b.close();
  fout_fq1.close();
  fout_fq2.close();
  fout.close();
  kout.close();
  cerr << "Connect reads finished! " << endl;

  if(mode == 0) {
    cout << "Simple connection table:" << endl;
  }

  else {
    cout << "Kmer frequency based connection table:" << endl;
  }

  cout << "total_pairs\tconnected_pairs\tconnect_ratio(%)\tlow_quality_pairs\tlow_quality_ratio(%)\n"
       << total_pairs << "\t" << connect_pairs << "\t" << double(connect_pairs) / double(total_pairs) * 100.0 << "\t"
       << low_quality_pairs << "\t" << double(low_quality_pairs) / double(total_pairs) * 100.0 << endl;

  if(table_file.size() > 0 && mode == 3) {
    uint64_t unconnect_pair = total_pairs - connect_pairs;
    cerr << "starting to do cross connect for " << unconnect_pair << " pair reads\n";

    if(cross_read_type) {
      cross_connect(fq1, fq2, out_connect_file, Cross_read_list_file, Pair_kmer_file);
    }

    else {
      cross_connect(fq1, fq2, file1, file2, Cross_read_list_file, Pair_kmer_file);
    }
  }
}

void set_command_line(int argc, const char * const * argv)
{
  static string command_line_string;

  for(int i = 0; i < argc; i++) {
    command_line_string += argv[i];

    if(i != argc - 1) {
      command_line_string += ' ';
    }
  }

  command_line = command_line_string.c_str();
}


int main(int argc, char * argv[])
{
  time_t time_start, time_end;
  time_start = time(NULL);

  if(argc < 2) {
    usage();
  }

  mGetOptions(argc, argv);
  set_command_line(argc, argv);
  cerr << "Program start.." << "\n"
       << "Program:\t" << PROGRAM << "\n"
       << "Version:\t" << VERSION << "\n"
       << "Author:\t\tBGI-ShenZhen\n"
       << "CompileDate:\t" << __DATE__ << " time: " << __TIME__ << "\n"
       << "Current time:\t" << ctime(&time_start)
       << "Command line:\t" << command_line << endl;

  if(TESTING) {
    FALSECON.open("kmerfreq_false_connect.inf");
    FALSECON2.open("cross_false_connect.inf");
    CERRINF.open("kmerfreq_cerr.inf");
    SCOREOUT.open("length.score.list");
    CROSSCERRINF.open("cross_cerr.inf");
  }

  if(a_file.empty() || b_file.empty() || out_connect_file.empty() || fq1.empty() || fq2.empty()) {
    cerr << "Error: input file is empty, please check the set of input file!" << endl;
    exit(-1);
  }

  igzstream check;
  check.open(a_file.c_str());

  if(!check) {
    cerr << "Error: can not open file: " << a_file << ", please check the input option -a!" << endl;
    exit(-1);
  }

  check.close();
  check.open(b_file.c_str());

  if(!check) {
    cerr << "Error: can not open file: " << b_file << ", please check the input option -b!" << endl;
    exit(-1);
  }

  check.close();

  if(mode != 0 && table_file.size() == 0 && table_file_by_jellyfish.size() == 0) {
    cerr << "Error: please input kmer frequency table file when you set mode >0 " << endl;
    exit(-1);
  }

  if(mode == 2 || mode == 4) {
    igzstream check;
    check.open(out_connect_file.c_str());

    if(!check) {
      cerr << "Error: can not open file: " << out_connect_file << ", please check the input option -o!" << endl;
      exit(-1);
    }

    check.close();
    check.open(fq1.c_str());

    if(!check) {
      cerr << "Error: can not open file: " << fq1 << ", please check the input option -2!" << endl;
      exit(-1);
    }

    check.close();
    check.open(fq2.c_str());

    if(!check) {
      cerr << "Error: can not open file: " << fq2 << ", please check the input option -3!" << endl;
      exit(-1);
    }

    check.close();
  }

  if(mode == 2) {
    igzstream check;
    check.open(Pair_kmer_file.c_str());

    if(!check) {
      cerr << "Error: can not open file: " << Pair_kmer_file << ", please check the input option -K!" << endl;
      exit(-1);
    }

    check.close();
  }

  if(mode == 4 && Cross_info_list_file.size() == 0) {
    cerr << "Error: please input cross-information list file(each line for one file) when mode=4" << endl;
    exit(-1);
  }

  if(mode == 4 && Cross_info_list_file.size() > 0) {
    igzstream CrossInfoList;
    CrossInfoList.open(Cross_info_list_file.c_str());

    if(!CrossInfoList) {
      cerr << "Fail open cross information list file: " << Cross_info_list_file << endl;
      exit(-1);
    }

    string cross_info_file;
    int file_count = 0;

    while(getline(CrossInfoList, cross_info_file, '\n')) {
      if(cross_info_file.size() == 0) {
        continue;
      }

      file_count++;
    }

    CrossInfoList.close();

    if(Batch_num != file_count) {
      cerr << "Error: Batch number do not equal to the number of cross information file, please check the option -C and -D !" << endl;
      exit(-1);
    }
  }

  //cerr create mismatch table.
  mis_array = new int[MaxOverLen];
  get_mis_table(MATCH_CUTOFF, mis_array, MaxOverLen);
  sum_num = (int)pow(2.0, double(MaxMismatch));
  matrix = new bool*[sum_num];

  for(int i = 0; i < sum_num; i++) {
    matrix[i] = new bool[MaxMismatch];
  }

  for(int i = 0; i < sum_num; i++) {
    int n = i;

    for(int j = 0; j < MaxMismatch; j++) {
      if(n == 0) {
        matrix[i][j] = 0;
      }

      else {
        matrix[i][j] = n % 2;
        n = n / 2;
      }
    }
  }

  //get the theory max number kmers
  uint64_t total = 0;

  for(int i = 0; i < Kmer; i++) {
    total = (total << 2) | 0x3;
  }

  uint64_t array_size = (total / 8 + 1);

  if((table_file_by_jellyfish.size() > 0 || table_file.size() > 0) && mode > 0) {
    freq = new uint8_t[array_size * 2];

    for(uint64_t i = 0; i < array_size * 2; i++) {
      freq[i] = 0;
    }

    if(table_file_by_jellyfish.size() > 0) {
      load_jellyfish_kmerfreq(table_file_by_jellyfish, freq, Freq_low, Freq_normal, Freq_high);
    }

    else {
      if(table_len_file.size() > 0) {
        load_kmerfreq(table_file, table_len_file, freq, Freq_low, Freq_normal, Freq_high);
      }

      else {
        load_kmerfreq(table_file, freq, Freq_low, Freq_normal, Freq_high);
      }
    }

    cerr << "finish loading kmer freq table!\n";
  }

  else {
    cerr << "No kmer freq table exists, so just connect with simple mode!\n";
  }

  time_end = time(NULL);
  cerr << " Run time: " << time_end - time_start << "s." << endl;

  //only run cross connect.
  if(mode == 2 || mode == 4) {
    igzstream fin(Pair_kmer_file.c_str());

    if(fin) {
      fin.close();
      cerr << "loading kmer pair directly!\n";

      if(cross_read_type) {
        cross_connect(fq1, fq2, out_connect_file, Cross_read_list_file, Pair_kmer_file);
      }

      else {
        cross_connect(fq1, fq2, a_file, b_file, Cross_read_list_file, Pair_kmer_file);
      }
    }

    else {
      cerr << "Can't open file, Please check the file: " << Pair_kmer_file << endl;
      exit(2);
    }

    if(table_file.size() > 0 || table_file_by_jellyfish.size() > 0) {
      delete []freq;
    }

    for(int i = 0; i < sum_num; i++) {
      delete []matrix[i];
    }

    delete []matrix;
    delete []mis_array;
    time_end = time(NULL);
    cerr << " Run time: " << time_end - time_start << "s." << endl;
    return 0;
  }

  //connection main function.
  load_file(a_file, b_file);
  cerr << "All done, Thank you!\n";

  for(int i = 0; i < sum_num; i++) {
    delete []matrix[i];
  }

  delete []matrix;

  if(table_file.size() > 0 || table_file_by_jellyfish.size() > 0) {
    delete []freq;
  }

  delete []mis_array;
  time_end = time(NULL);
  cerr << " Run time: " << time_end - time_start << "s." << endl;
  return 0;
}
