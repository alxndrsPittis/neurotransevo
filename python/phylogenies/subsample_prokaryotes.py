#!/usr/bin/env python3
"""Taxonomy-aware subsampling of prokaryotic hmmsearch hits (Methods).

If a family has more than --min-hits prokaryotic hits, keep at most --per-phylum
sequences per bacterial or archaeal phylum (NCBI taxonomy). Sequence IDs are
expected as "<taxid>.<protein id>" (UniProt Reference Proteomes naming).

    python subsample_prokaryotes.py hits.tblout proteomes.fasta > prokaryotes.fasta
"""
import argparse
import sys

from ete3 import NCBITaxa

BACTERIA, ARCHAEA = 2, 2157


def read_fasta(path, keep):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            name = name if name in keep else None
            if name:
                seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tblout", help="hmmsearch --tblout against prokaryotic reference proteomes")
    ap.add_argument("fasta", help="prokaryotic reference proteomes (FASTA)")
    ap.add_argument("--min-hits", type=int, default=300)
    ap.add_argument("--per-phylum", type=int, default=100)
    args = ap.parse_args()

    hits = [l.split()[0] for l in open(args.tblout) if not l.startswith("#")]
    hits = list(dict.fromkeys(hits))
    if len(hits) > args.min_hits:
        ncbi = NCBITaxa()
        taxids = {int(h.split(".")[0]) for h in hits}
        lineages = ncbi.get_lineage_translator(list(taxids))
        ranks = ncbi.get_rank({t for lin in lineages.values() for t in lin})
        kept, per_phylum = [], {}
        for h in hits:
            lineage = lineages.get(int(h.split(".")[0]), [])
            if BACTERIA not in lineage and ARCHAEA not in lineage:
                continue
            phylum = next((t for t in lineage if ranks.get(t) == "phylum"), None)
            if per_phylum.get(phylum, 0) < args.per_phylum:
                per_phylum[phylum] = per_phylum.get(phylum, 0) + 1
                kept.append(h)
        hits = kept
    seqs = read_fasta(args.fasta, set(hits))
    for h in hits:
        if h in seqs:
            sys.stdout.write(f">{h}\n{seqs[h]}\n")
    print(f"{len(hits)} prokaryotic sequences kept", file=sys.stderr)


if __name__ == "__main__":
    main()
