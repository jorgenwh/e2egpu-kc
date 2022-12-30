import numpy as np
import cupy as cp

from e2epy_backend import HashTable

class Counter(HashTable):
    def __init__(self, keys, capacity=0):
        assert isinstance(keys, (np.ndarray, cp.ndarray)), "Keys must be numpy or cupy ndarray"
        assert keys.dtype == np.uint64, "Keys' dtype must be uint64"

        if capacity == 0:
            capacity = int(keys.size * 17331)

        assert capacity > keys.size, "Capacity must be greater than keys.size"

        if isinstance(keys, np.ndarray):
            super().__init__(keys, capacity)
        elif isinstance(keys, cp.ndarray):
            super().__init__(keys.data.ptr, keys.size, capacity)

    def count(self, filename: str):
        raise NotImplementedError

