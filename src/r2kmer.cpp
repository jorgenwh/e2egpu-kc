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

int reads_to_kmers(const char *reads, const int num_reads, 
    const int read_length, const int kmer_size, uint64_t **kmers)
{
  const int kmers_per_read = read_length - (kmer_size - 1);
  const int num_kmers = kmers_per_read * num_reads;
  const uint64_t rightshift = (32 - kmer_size)*2;
  *kmers = new uint64_t[num_kmers];

  int read_index;
  int kmer_index;

  // Iterate over reads
  for (int i = 0; i < num_reads; i++)
  {
    read_index = i*read_length;
    kmer_index = i*kmers_per_read;

    // Iterate over kmers
    for (int j = 0; j < kmers_per_read; j++)
    {
      process_kmer(&reads[read_index+j], (*kmers)[kmer_index + j], kmer_size, rightshift);
    }
  }

  return num_kmers;
}
