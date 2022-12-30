import numpy as np
import cupy as cp

from e2epy import Counter

if __name__ == "__main__":
    counter_keys = np.load("../../kmer-counting-experimentation/data/npy/uniquekmersACGT.npy")
    counter = Counter(counter_keys)
