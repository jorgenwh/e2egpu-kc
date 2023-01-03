#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
#include "kernels.h"

namespace kernels 
{

__device__ static uint8_t bitbase_lookup[4] = {0, 1, 3, 2};
static const uint64_t bitbase_mask = 0b11;
__device__ static uint8_t rightshifts[33] = {
  64, 62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 
  32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0
};

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

__device__ __forceinline__ static uint64_t get_bitencoded_kmer(
    const char *read, const int kmer_size)
{
  uint64_t kmer = 0;
  uint64_t shift;
  uint64_t bitbase;
  for (int i = 0; i < kmer_size; i++)
  {
    shift = ((32 - i) - rightshifts[kmer_size]) * 2;
    bitbase = bitbase_lookup[(read[i] >> 1) & bitbase_mask];
    kmer |= (bitbase << shift);
  }
  return kmer;
}

__device__ __forceinline__ static void count_kmer(
    uint64_t *table_keys, uint32_t *table_values, const int table_capacity, 
    const uint64_t kmer)
{
  uint64_t hash = murmur_hash(kmer) % table_capacity;
  while (true)
  {
    uint64_t table_key = table_keys[hash];
    if (table_key == kEmpty)
    {
      break;
    }
    if (table_key == kmer)
    {
      atomicAdd((unsigned int *)&table_values[hash], 1);
      break;
    }
    hash = (hash + 1) % table_capacity;
  }
}

__global__ void count_reads_kernel(uint64_t *table_keys, uint32_t *table_values, 
    const int table_capacity, const char *reads, const int num_reads, const int read_length, 
    const int kmer_size)
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
    insert_key = get_bitencoded_kmer(&reads[thread_id*read_length + i], kmer_size);
    count_kmer(table_keys, table_values, table_capacity, insert_key);
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

  int grid_size = SDIV(num_reads, thread_block_size);
  count_reads_kernel<<<grid_size, thread_block_size>>>(table_keys, table_values, table_capacity, 
      reads, num_reads, read_length, kmer_size);
}

__global__ void count_reads_single_kernel(uint64_t *table_keys, uint32_t *table_values, 
    const int table_capacity, const char *reads, const int num_reads, const int read_length, 
    const int kmer_size, const uint64_t num_kmers)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_id >= num_kmers)
  {
    return;
  }

  int read_index = thread_id / read_length;
  int kmer_index = thread_id % (read_length - (kmer_size - 1));
  int global_index = read_index*read_length + kmer_index;
  printf("thread_id=%d, read_index=%d, kmer_index=%d, global_index=%d\n", 
      thread_id, read_index, kmer_index, global_index);

  int kmer = get_bitencoded_kmer(&reads[global_index], kmer_size); 
  count_kmer(table_keys, table_values, table_capacity, kmer);
}

void count_reads_single(uint64_t *table_keys, uint32_t *table_values, const int table_capacity,
    const char *reads, const int num_reads, const int read_length, const int kmer_size)
{
  int min_grid_size;
  int thread_block_size;
  cuda_errchk(cudaOccupancyMaxPotentialBlockSize(
        &min_grid_size, &thread_block_size, 
        count_reads_single_kernel, 0, 0));

  const uint64_t num_kmers = (read_length - (kmer_size - 1)) * num_reads;

  int grid_size = SDIV(num_kmers, thread_block_size);
  count_reads_single_kernel<<<grid_size, thread_block_size>>>(
      table_keys, table_values, table_capacity, 
      reads, num_reads, read_length, kmer_size, num_kmers);
}

__global__ void count_raw_reads_kernel(char *chunk, const int bytes,
    const int header_length, const int read_length, const int kmer_size)
{
  ;
}

void count_raw_reads(char *chunk, const int bytes, 
    const int header_length, const int read_length, const int kmer_size)
{
  ;
}

__global__ void lookup_kernel(uint64_t *table_keys, uint32_t *table_values, 
    const int table_capacity, uint64_t *keys, uint32_t *values, const int size)
{
  int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (thread_id >= size)
  {
    return;
  }

  uint64_t lookup_key = keys[thread_id];
  uint64_t hash = murmur_hash(lookup_key) % table_capacity;

  while (true) 
  {
    uint64_t table_key = table_keys[hash];
    if (table_key == lookup_key || table_key == kEmpty) 
    {
      values[thread_id] = (table_key == lookup_key) ? table_values[hash] : 0;
      return;
    }
    hash = (hash + 1) % table_capacity;
  }
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
