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

  const int header_length = 10;
  const int read_length = 150;
  const int reads_per_chunk = 50000;
  const int kmer_size = 31;

  auto t1 = std::chrono::high_resolution_clock::now();
  counter.count("data/testreads20m.fa", header_length, read_length, reads_per_chunk, kmer_size);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

  uint32_t *counts = new uint32_t[num_kmers];
  counter.get(unique_kmers, counts, num_kmers);

  for (int i = 0; i < 50; i++)
  {
    std::cout << counts[i] << ", ";
  }
  std::cout << "\n";
  std::cout << elapsed << "\n";

  delete[] counts;
  delete[] unique_kmers;
  return EXIT_SUCCESS;
}
