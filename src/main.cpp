#include <iostream>
#include <chrono>
#include <inttypes.h>
#include <bitset>

#include "common.h"
#include "io.h"
#include "r2kmer.h"

int main(int argc, char **argv)
{
  auto t1 = std::chrono::high_resolution_clock::now();

  FastaReader reader("data/small.fa");
  const int reads_per_chunk = 57500;

  while (!reader.done())
  {
    char *read;
    reader.read_chunk(&read, reads_per_chunk);

    uint64_t *kmers;
    int num_kmers = read_to_kmers(read, READ_LENGHT, &kmers, 31);

    std::cout << std::bitset<64>(kmers[0]) << "\n";

    delete[] read;
    delete[] kmers;
    break;
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

  std::cout << elapsed << "\n";

  return EXIT_SUCCESS;
}
