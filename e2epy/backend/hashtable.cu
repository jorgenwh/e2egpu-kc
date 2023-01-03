#include <iostream>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
#include "io.h"
#include "kernels.h"
#include "hashtable.h"

HashTable::HashTable(
    const uint64_t *keys, const bool keys_on_device, const int size, const int capacity)
{
  size_m = size;
  capacity_m = capacity;

  cuda_errchk(
      cudaMalloc(&keys_m, capacity*sizeof(uint64_t)));
  cuda_errchk(
      cudaMemset(keys_m, 0xFF, capacity*sizeof(uint64_t)));
  cuda_errchk(
      cudaMalloc(&values_m, capacity*sizeof(uint32_t)));
  cuda_errchk(
      cudaMemset(values_m, 0, capacity*sizeof(uint32_t)));

  uint64_t *keys_d;
  if (!keys_on_device)
  {
    cuda_errchk(
        cudaMalloc(&keys_d, size*sizeof(uint64_t)));
    cuda_errchk(
        cudaMemcpy(keys_d, keys, size*sizeof(uint64_t), cudaMemcpyHostToDevice));
  }

  // Synchronize because cudaMemset is asynchronous with respect to host
  cuda_errchk(cudaDeviceSynchronize());

  kernels::initialize_hashtable(
      keys_m, capacity, keys_on_device ? keys : keys_d, size);

  if (!keys_on_device)
  {
    cuda_errchk(cudaFree(keys_d));
  }
}

HashTable::~HashTable()
{
  cuda_errchk(cudaFree(keys_m));
  cuda_errchk(cudaFree(values_m));
}

void HashTable::count(const char *filename, const int header_length, const int read_length, 
    const int reads_per_chunk, const int kmer_size)
{
  int i = 0;
  FastaReader reader(filename);
  while (!reader.done())
  {
    char *reads;
    int num_reads = reader.read_chunk(&reads, reads_per_chunk, header_length, read_length);

    char *reads_d;
    cuda_errchk(
        cudaMalloc(&reads_d, read_length*num_reads*sizeof(char)));
    cuda_errchk(
        cudaMemcpy(reads_d, reads, read_length*num_reads*sizeof(char), cudaMemcpyHostToDevice));

    //kernels::count_reads(
        //keys_m, values_m, capacity_m, reads_d, num_reads, read_length, kmer_size);
    kernels::count_reads_single(
        keys_m, values_m, capacity_m, reads_d, num_reads, read_length, kmer_size);

    delete[] reads;
    cuda_errchk(cudaFree(reads_d));

    i++;
    std::cout << "Chunk " << i << "\r";
  }
  std::cout << "Chunk " << i << "\n";
}


void HashTable::count_fasta_chunks(const char *filename, const int chunk_size, 
    const int header_length, const int read_length, const int kmer_size)
{
  int i = 0;
  FastaReader reader(filename);
  while (!reader.done())
  {
    char *chunk;
    int bytes_read = reader.read_raw_chunk(&chunk, chunk_size);

    kernels::count_raw_reads(chunk, bytes_read, header_length, read_length, kmer_size);

    delete[] chunk;

    i++;
    std::cout << "Chunk " << i << "\r";
  }
  std::cout << "Chunk " << i << "\n";
}

void HashTable::lookup(const uint64_t *keys, uint32_t *values, int size) const
{
  uint64_t *keys_d;
  uint32_t *values_d;
  cuda_errchk(cudaMalloc(&keys_d, size*sizeof(uint64_t)));
  cuda_errchk(cudaMalloc(&values_d, size*sizeof(uint32_t)));
  cuda_errchk(cudaMemcpy(keys_d, keys, size*sizeof(uint64_t), cudaMemcpyHostToDevice));

  kernels::lookup(keys_m, values_m, capacity_m, keys_d, values_d, size);

  cuda_errchk(cudaMemcpy(values, values_d, size*sizeof(uint32_t), cudaMemcpyDeviceToHost)); 
  cuda_errchk(cudaFree(keys_d));
  cuda_errchk(cudaFree(values_d));
}
