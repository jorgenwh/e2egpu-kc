#include <inttypes.h>
#include <bitset>
#include <iostream>

#include "common.h"
#include "r2kmer.h"

uint8_t bitbase_lookup[4] = {
  0, 1, 3, 2
};

inline void process_kmer(
    const char *bases, uint64_t &kmer, const int kmer_size, const uint64_t rightshift)
{
  uint64_t shift;
  uint64_t bitbase;

  for (int i = 0; i < kmer_size; i++)
  {
    shift = ((32 - i) - rightshift) * 2;
    bitbase = bitbase_lookup[(bases[i] >> 1) & bitbase_mask];
    kmer |= (bitbase << shift);
  }
}

int read_to_kmers(
    const char *read, const int read_length, uint64_t **kmers, const int kmer_size)
{
  const int num_kmers = read_length - kmer_size;
  const uint64_t rightshift = (32 - kmer_size)*2;
  *kmers = new uint64_t[num_kmers];

  for (int i = 0; i < read_length-kmer_size; i++)
  {
    process_kmer(&read[i], (*kmers)[i], kmer_size, rightshift);
  }

  return num_kmers;
}
