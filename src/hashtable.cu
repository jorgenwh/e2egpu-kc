#include <iostream>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"
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
      keys_m, values_m, keys_on_device ? keys : keys_d, size, capacity);

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
