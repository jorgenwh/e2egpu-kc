#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
#include "kernels.h"

class HashTable
{
public:
  HashTable() = default;
  HashTable(const uint64_t *keys, const bool keys_on_device, 
      const int size, const int capacity);
  ~HashTable();

  int size() const { return size_m; }
  int capacity() const { return capacity_m; }

  void count(const char *filename, const int header_length, const int read_length, 
      const int reads_per_chunk, const int kmer_size);
  void get(const uint64_t *keys, uint32_t *values, int size) const;

private:
  int size_m;
  int capacity_m;

  uint64_t *keys_m;
  uint32_t *values_m;
};

#endif // HASHTABLE_H_
