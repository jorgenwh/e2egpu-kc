#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
#include "kernels.h"

namespace kernels 
{

__device__ uint8_t bitbase_lookup[4] = {0, 1, 3, 2};
static const uint64_t bitbase_mask = 0b11;

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

__global__ void initialize_hashtable_kernel(uint64_t *table_keys, const int table_capacity, 
    const uint64_t *keys, const int size)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_id >= size) 
  {
    return;
  }

  uint64_t insert_key = keys[thread_id];
  uint64_t hash = murmur_hash(insert_key) % table_capacity;

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
    hash = (hash + 1) % table_capacity;
  }
}

void initialize_hashtable(uint64_t *table_keys, const int table_capacity,
    const uint64_t *keys, const int size)
{
  int min_grid_size;
  int thread_block_size;
  cuda_errchk(cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &thread_block_size, 
        initialize_hashtable_kernel, 0, 0));

  int grid_size = SDIV(size, thread_block_size);
  initialize_hashtable_kernel<<<grid_size, thread_block_size>>>(
      table_keys, table_capacity, keys, size);
}

__device__ __forceinline__ static uint64_t get_bitencoded_kmer(const char *read, 
    const int kmer_size, const uint64_t rightshift)
{
  uint64_t kmer = 0;
  uint64_t shift;
  uint64_t bitbase;
  for (int i = 0; i < kmer_size; i++)
  {
    shift = ((32 - i) - rightshift) * 2;
    bitbase = bitbase_lookup[(read[i] >> 1) & bitbase_mask];
    kmer |= (bitbase << shift);
  }
  return kmer;
}

__global__ void count_reads_kernel(uint64_t *table_keys, uint32_t *table_values, 
    const int table_capacity, const char *reads, const int num_reads, const int read_length, 
    const int kmer_size, const uint64_t rightshift)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_id >= num_reads)
  {
    return;
  }

  // Iterate over valid kmers in the read
  uint64_t insert_key;
  uint64_t hash;
  for (int i = 0; i < read_length - (kmer_size - 1); i++)
  {
    // Get the current kmer bit encoding
    insert_key = get_bitencoded_kmer(&reads[thread_id*read_length + i], kmer_size, rightshift);
    hash = murmur_hash(insert_key) % table_capacity;

    // Count the kmer
    while (true)
    {
      uint64_t table_key = table_keys[hash];
      if (table_key == kEmpty)
      {
        break;
      }
      if (table_key == insert_key)
      {
        atomicAdd((unsigned int *)&table_values[hash], 1);
        break;
      }
      hash = (hash + 1) % table_capacity;
    }

    //__syncthreads();
    // TODO: add revcomp support and test correctness against old cucounter
  }
}

void count_reads(uint64_t *table_keys, uint32_t *table_values, const int table_capacity,
    const char *reads, const int num_reads, const int read_length, const int kmer_size)
{
  int min_grid_size;
  int thread_block_size;
  cuda_errchk(cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &thread_block_size, 
        count_reads_kernel, 0, 0));

  const uint64_t rightshift = (32 - kmer_size)*2;

  int grid_size = SDIV(num_reads, thread_block_size);
  count_reads_kernel<<<grid_size, thread_block_size>>>(table_keys, table_values, table_capacity, 
      reads, num_reads, read_length, kmer_size, rightshift);
}

__global__ void lookup_kernel(uint64_t *table_keys, uint32_t *table_values, const int table_capacity,
    uint64_t *keys, uint32_t *values, const int size)
{
  ;
}

void lookup(uint64_t *table_keys, uint32_t *table_values, const int table_capacity, 
    uint64_t *keys, uint32_t *values, const int size)
{
  int min_grid_size;
  int thread_block_size;
  cuda_errchk(cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &thread_block_size, 
        lookup_kernel, 0, 0));

  int grid_size = SDIV(size, thread_block_size);
  lookup_kernel<<<grid_size, thread_block_size>>>(table_keys, table_values, table_capacity, 
      keys, values, size);
}

} // namespace kernels
