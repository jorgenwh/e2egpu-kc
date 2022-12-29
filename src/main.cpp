#include <iostream>
#include <chrono>
#include <inttypes.h>
#include <bitset>
#include <cuda_runtime.h>

#include "common.h"
#include "io.h"
#include "r2kmer.h"
#include "hashtable.h"

void read_unique_kmers(uint64_t **buffer, int &num_kmers)
{
  FILE *file = fopen("data/uniquekmers.bin", "rb");
  const int buffer_size = 153627144;
  (*buffer) = new uint64_t[buffer_size];
  size_t elems_read = fread((*buffer), sizeof(uint64_t), buffer_size, file);
  num_kmers = buffer_size;
}

int main(int argc, char **argv)
{
  int capacity = 200000033;
  int num_kmers;
  uint64_t *unique_kmers;
  read_unique_kmers(&unique_kmers, num_kmers);

  HashTable counter(unique_kmers, false, num_kmers, capacity);

  delete[] unique_kmers;
  return EXIT_SUCCESS;
}
