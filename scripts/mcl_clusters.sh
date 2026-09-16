#!/usr/bin/env bash
# Similarity-based clustering of a superfamily (MFS, GPCR-A) before phylogenetics (Methods).
#   scripts/mcl_clusters.sh <family.fasta> <inflation> [outdir]     # inflation 1.4 (MFS), 6.0 (GPCR-A)
set -euo pipefail
fasta=$1; inflation=$2; out=${3:-$(dirname "$fasta")}
stem=$out/$(basename "${fasta%.*}")
threads=${THREADS:-8}

makeblastdb -in "$fasta" -dbtype prot -out "$stem.db" > /dev/null
blastp -query "$fasta" -db "$stem.db" -evalue 1e-5 -outfmt "6 qseqid sseqid evalue" \
       -num_threads "$threads" > "$stem.allVSall.tsv"
# edge weight = -log10(E-value), capped for E = 0
awk 'BEGIN{OFS="\t"} $1!=$2 {e=$3; w=(e==0)?200:-log(e)/log(10); print $1,$2,w}' "$stem.allVSall.tsv" > "$stem.abc"
mcl "$stem.abc" --abc -I "$inflation" -te "$threads" -o "$stem.I${inflation/./}.clusters"
echo "clusters: $stem.I${inflation/./}.clusters"
