#ifndef KERNELS_H_
#define KERNELS_H_

#include <inttypes.h>
#include <cuda_runtime.h>

#include "common.h"

namespace kernels
{

void initialize_hashtable(uint64_t *table_keys, uint32_t *table_values, 
    const uint64_t *keys, const int size, const int capacity);

} // namespace kernels

#endif // KERNELS_H_
