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
            raise NotImplementedError
            super().__init__(keys.data.ptr, keys.size, capacity)

    def count(self, 
            filename: str, 
            header_length: int, 
            read_length: int, 
            reads_per_chunk: int, 
            kmer_size: int):
        super().count(filename, header_length, read_length, reads_per_chunk, kmer_size)

    def count_fasta_chunks(self, 
            filename: str, 
            chunk_size: int,
            header_length: int, 
            read_length: int, 
            kmer_size: int):
        super().count_fasta_chunks(filename, chunk_size, header_length, read_length, kmer_size)

    def __getitem__(self, keys):
        assert isinstance(keys, (np.ndarray, cp.ndarray)), "Keys must be numpy or cupy ndarray"
        assert keys.dtype == np.uint64, "Keys' dtype must be uint64"

        if isinstance(keys, np.ndarray):
            return super().lookup(keys)
        elif isinstance(keys, cp.ndarray):
            raise NotImplementedError
            counts = cp.zeros_like(keys, dtype=np.uint32)
            super().lookup(keys.data.ptr, counts.data.ptr, keys.size)
            return counts

