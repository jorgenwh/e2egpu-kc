#ifndef KERNELS_H_
#define KERNELS_H_

#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"

namespace kernels
{

void initialize_hashtable(uint64_t *table_keys, const int table_capacity,
    const uint64_t *keys, const int size);

void count_reads(uint64_t *table_keys, uint32_t *table_values, const int table_capacity,
    const char *reads, const int num_reads, const int read_length, const int kmer_size);

void lookup(uint64_t *table_keys, uint32_t *table_values, const int table_capacity, 
    uint64_t *keys, uint32_t *values, const int size);

} // namespace kernels

#endif // KERNELS_H_
