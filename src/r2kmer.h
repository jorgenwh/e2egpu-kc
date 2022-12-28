#ifndef R2KMER_H_
#define R2KMER_H_

#include <inttypes.h>

#include "common.h"

static const uint64_t bitbase_mask = 0b11;
extern uint8_t bitbase_lookup[4];

inline void process_kmer(
    const char *bases, uint64_t &kmer, const int kmer_size, const uint64_t rightshift);

int read_to_kmers(
    const char *read, const int read_length, uint64_t **kmers, const int kmer_size);

#endif // R2KMER_H_
