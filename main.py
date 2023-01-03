import time

import numpy as np
import cupy as cp

from e2epy import Counter

HEADER_LENGTH = 10
READ_LENGTH = 150
READS_PER_CHUNK = 50000
KMER_SIZE = 31

SMALL_FASTA = "data/small.fa"
BIG_FASTA = "data/testreads20m.fa"

if __name__ == "__main__":
    counter_keys = np.load("../kmer-counting-experimentation/data/npy/uniquekmersACGT.npy")[:30000000]
    counter = Counter(counter_keys)

    t1 = time.time()

    """counter.count(
            filename=SMALL_FASTA, 
            header_length=HEADER_LENGTH, 
            read_length=READ_LENGTH, 
            reads_per_chunk=READS_PER_CHUNK, 
            kmer_size=KMER_SIZE)"""

    counter.count_fasta_chunks(
            filename=BIG_FASTA,
            chunk_size=12500000,
            header_length=10,
            read_length=150,
            kmer_size=31)

    elapsed = time.time() - t1
    print(f"count time: {round(elapsed, 6)} s")

    exit()

    counts = counter[counter_keys]
    print(type(counts))
    print(counts.shape)
    print(counts.dtype)
    print(counts[:10])


