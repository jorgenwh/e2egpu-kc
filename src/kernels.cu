#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
#include "kernels.h"

namespace kernels 
{

__device__ __forceinline__ static uint64_t word_reverse_complement(
    const uint64_t kmer, uint8_t kmer_size) 
{
  uint64_t res = ~kmer;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  return (res >> (2 * (32 - kmer_size)));
}

__device__ __forceinline__ static uint64_t murmur_hash(uint64_t kmer) 
{
#ifdef __USE_MURMUR_HASH__
  kmer ^= kmer >> 33;
  kmer *= 0xff51afd7ed558ccd;
  kmer ^= kmer >> 33;
  kmer *= 0xc4ceb9fe1a85ec53;
  kmer ^= kmer >> 33;
#endif // __USE_MURMUR_HASH__
  return kmer;
}

__global__ void initialize_hashtable_kernel(uint64_t *table_keys, uint32_t *table_values, 
    const uint64_t *keys, const int size, const int capacity)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_id >= size) 
  {
    return;
  }

  uint64_t insert_key = keys[thread_id];
  uint64_t hash = murmur_hash(insert_key) % capacity;

  while (true) 
  {
    unsigned long long int *table_key_ptr = 
      reinterpret_cast<unsigned long long int *>(&table_keys[hash]);
    uint64_t old = atomicCAS(table_key_ptr, kEmpty, insert_key);

    const bool inserted = (old == kEmpty || old == insert_key);

    if (inserted)
    {
      return;
    }
    hash = (hash + 1) % capacity;
  }
}

void initialize_hashtable(uint64_t *table_keys, uint32_t *table_values, 
    const uint64_t *keys, const int size, const int capacity)
{
  int min_grid_size;
  int thread_block_size;
  cuda_errchk(cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &thread_block_size, 
        initialize_hashtable_kernel, 0, 0));

  int grid_size = SDIV(size, thread_block_size);
  initialize_hashtable_kernel<<<grid_size, thread_block_size>>>(
      table_keys, table_values, keys, size, capacity);
}

} // namespace kernels
