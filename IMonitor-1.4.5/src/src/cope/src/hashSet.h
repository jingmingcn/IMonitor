/*
 * hashSet.h
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
 * File          : hashSet.h
 * Revision      : 1.1.2
 * Revision_date : 2012/10/08
 * Author(s)     : Zhenyu Li, Binghang Liu
 *
 *============================================================================
 */
 
#ifndef __HASHSET_H_
#define __HASHSET_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "seqKmer.h"

using namespace std;

//define the entity data structure (elements of the hashset-array)
extern int max_freq_cut;
extern int Kmer;

typedef struct {
  uint64_t kmer: 38; //max for 17mer, 34->38.
  uint64_t pos: 32; //33M kmer pair number, 25->33.
  uint64_t len: 25; //number in array, max length is 16, 4->25.
  uint64_t flag: 1; //reverse 0 or not 1.
} Entity;

//define the hashSet data structure(include data and control parameters)
typedef struct {
  uint32_t e_size; //the byte size for each entity
  Entity * array; //the array to store the body data
  uint64_t size; //the size of the array declared above
  uint64_t count; //the number of current stored entities
  uint64_t count_conflict; //the number of conflict
  uint64_t max; //the max number of entities allowed to store
  uint64_t iter_ptr; //the pointer to the current processing entity
  float load_factor; //the ratio (allowed entity number / total array size)
  uint8_t * nul_flag; //the array to store null flags, one bit for an entity
  uint8_t * del_flag; //the array to store delete flags, one bit for an entity
} HashSet;


uint64_t hash_code(Entity * e);
int hash_equal(Entity * a, Entity * b);
void free_hash(HashSet * set);

int is_prime(uint64_t num);
uint64_t find_next_prime(uint64_t num);
HashSet * init_hashset(uint64_t init_size, float load_factor);

int is_entity_null(uint8_t * nul_flag, uint64_t idx);
int set_entity_fill(uint8_t * nul_flag, uint64_t idx);
int is_entity_delete(uint8_t * del_flag, uint64_t idx);
int set_entity_delete(uint8_t * del_flag, uint64_t idx);

void enlarge_hashset(HashSet * set, uint64_t num);
int add_hashset(HashSet * set, Entity * entity);
uint64_t get_hashset(HashSet * set, Entity * entity);
int exists_hashset(HashSet * set, Entity * entity);
int delete_hashset(HashSet * set, Entity * entity);
uint64_t get_freq(HashSet * set, uint64_t * kmer);

#endif
