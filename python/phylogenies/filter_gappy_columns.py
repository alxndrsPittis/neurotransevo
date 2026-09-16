#!/usr/bin/env python3
"""Remove alignment columns with more than a given fraction of gaps (Methods).

    python filter_gappy_columns.py aln.fasta 0.9  > aln.gt01.fasta    # standard
    python filter_gappy_columns.py aln.fasta 0.99 > aln.gt001.fasta   # large (>1,000 seq) alignments
"""
import sys


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            names.append(line[1:]); seqs.append([])
        elif line:
            seqs[-1].append(line)
    return names, ["".join(s) for s in seqs]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    names, seqs = read_fasta(sys.argv[1])
    max_gap = float(sys.argv[2])
    n = len(seqs)
    keep = [i for i in range(len(seqs[0])) if sum(s[i] in "-." for s in seqs) / n <= max_gap]
    for name, s in zip(names, seqs):
        sys.stdout.write(f">{name}\n{''.join(s[i] for i in keep)}\n")
    print(f"kept {len(keep)}/{len(seqs[0])} columns", file=sys.stderr)


if __name__ == "__main__":
    main()
