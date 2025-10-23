#!/usr/bin/env python3
import argparse
from collections import Counter
import math

def parse_args():
    parser = argparse.ArgumentParser(
        description="Simplified Inchworm-like k-mer assembler with entropy filtering."
    )
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("k", type=int, help="K-mer size")
    parser.add_argument("--min_entropy", type=float, default=1.5,
                        help="Minimum Shannon entropy required for seed k-mers (default=1.5)")
    parser.add_argument("--top", type=int, default=None,
                        help="Show only the top N assembled contigs")
    return parser.parse_args()

def read_fastq_sequences(filename):
    with open(filename, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # '+'
            f.readline()  # quality
            yield seq

def count_kmers(fastq_file, k):
    kmers = Counter()
    for seq in read_fastq_sequences(fastq_file):
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            kmers[seq[i:i+k]] += 1
    return kmers

def shannon_entropy(seq):
    length = len(seq)
    if length == 0:
        return 0.0
    counts = Counter(seq)
    entropy = 0.0
    for base, count in counts.items():
        p = count / length
        entropy -= p * math.log2(p)
    return round(entropy, 4)

def find_extension(kmer, kmers, k, direction='right'):
    """Find next kmer that overlaps by k-1 bases in given direction."""
    if direction == 'right':
        suffix = kmer[1:]
        for base in "ACGT":
            candidate = suffix + base
            if kmers.get(candidate, 0) > 0:
                return candidate
    else:
        prefix = kmer[:-1]
        for base in "ACGT":
            candidate = base + prefix
            if kmers.get(candidate, 0) > 0:
                return candidate
    return None

def extend_contig(seed, kmers, k):
    """Extend a contig from a seed kmer in both directions."""
    contig = seed
    kmers[seed] = 0  # consume seed

    # extend right
    current = seed
    while True:
        next_kmer = find_extension(current, kmers, k, direction='right')
        if not next_kmer:
            break
        contig += next_kmer[-1]
        kmers[next_kmer] = 0
        current = next_kmer

    # extend left
    current = seed
    prefix_seq = ""
    while True:
        prev_kmer = find_extension(current, kmers, k, direction='left')
        if not prev_kmer:
            break
        prefix_seq = prev_kmer[0] + prefix_seq
        kmers[prev_kmer] = 0
        current = prev_kmer

    return prefix_seq + contig

def inchworm_assemble(kmers, k, min_entropy):
    assembled = []
    sorted_kmers = sorted(kmers.items(), key=lambda x: -x[1])
    for kmer, count in sorted_kmers:
        if kmers[kmer] == 0:
            continue  # already used
        if shannon_entropy(kmer) < min_entropy:
            continue
        contig = extend_contig(kmer, kmers, k)
        assembled.append((contig, len(contig)))
    return assembled

def main():
    args = parse_args()
    kmers = count_kmers(args.fastq, args.k)
    contigs = inchworm_assemble(kmers, args.k, args.min_entropy)
    contigs = sorted(contigs, key=lambda x: -x[1])

    if args.top:
        contigs = contigs[:args.top]

    print("contig_id\tlength\tsequence")
    for i, (seq, length) in enumerate(contigs, 1):
        print(f"contig_{i}\t{length}\t{seq}")

if __name__ == "__main__":
    main()
