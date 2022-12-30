import os

import numpy as np
import cupy as cp

from e2epy_backend import HashTable

class Counter(HashTable):
    def __init__(self, keys, capacity=0):
        assert isinstance(keys, (np.ndarray, cp.ndarray)), "Keys must be numpy or cupy ndarray"
        assert keys.dtype == np.uint64, "Keys' dtype must be uint64"

        if capacity == 0:
            capacity = int(keys.size * 1.7331)

        assert capacity > keys.size, "Capacity must be greater than keys.size"

        if isinstance(keys, np.ndarray):
            super().__init__(keys, capacity)
        elif isinstance(keys, cp.ndarray):
            super().__init__(keys.data.ptr, keys.size, capacity)


    def count(self, filename: str, header_length: int, read_length: int, reads_per_chunk: int, kmer_size: int):
        super().count(filename, header_length, read_length, reads_per_chunk, kmer_size)

