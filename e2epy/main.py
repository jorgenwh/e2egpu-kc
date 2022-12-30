import time

import numpy as np
import cupy as cp

from e2epy import Counter

HEADER_LENGTH = 10
READ_LENGTH = 150
READS_PER_CHUNK = 50000
KMER_SIZE = 31

SMALL_FASTA = "../data/small.fa"
BIG_FASTA = "../data/testreads20m.fa.fa"

if __name__ == "__main__":
    counter_keys = np.load("../../kmer-counting-experimentation/data/npy/uniquekmersACGT.npy")
    counter = Counter(counter_keys)

    t1 = time.time()
    counter.count(SMALL_FASTA, HEADER_LENGTH, READ_LENGTH, READS_PER_CHUNK, KMER_SIZE)
    elapsed = time.time() - t1
    print(f"count time: {round(elapsed, 3)} s")
