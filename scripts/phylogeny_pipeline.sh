#!/usr/bin/env bash
# Gene-family phylogeny pipeline (Methods: "Phylogenetic reconstruction analysis methods").
#
#   scripts/phylogeny_pipeline.sh <name> <pfam.hmm> <eukaryotes.fasta> <prokaryotes.fasta> [outdir]
#
# Steps: hmmsearch (--cut_ga) -> prokaryote subsampling -> alignment -> gap filtering -> ML tree.
#   <= 1,000 sequences: MAFFT L-INS-i (E-INS-i with MAFFT_MODE=einsi), columns >90% gaps removed, IQ-TREE
#                       (LG model selection, ultrafast bootstrap)
#   >  1,000 sequences: Clustal Omega, columns >99% gaps removed, FastTree
# For the superfamilies MFS and GPCR-A (Rhodopsin), sequences were first clustered by
# all-vs-all BLASTP + MCL (inflation 1.4 and 6.0) and the clusters of interest selected;
# see scripts/mcl_clusters.sh. Clades of interest were re-extracted and re-aligned with the
# standard (IQ-TREE) route. Tools: environment envs/phylogenies.yml.
set -euo pipefail

name=$1; hmm=$2; euk=$3; prok=$4; out=${5:-results/phylogenies/$name}
threads=${THREADS:-8}
mkdir -p "$out"
here=$(cd "$(dirname "$0")/.." && pwd)

hmmsearch --cut_ga --cpu "$threads" --tblout "$out/$name.euk.tblout" "$hmm" "$euk" > /dev/null
hmmsearch --cut_ga --cpu "$threads" --tblout "$out/$name.prok.tblout" "$hmm" "$prok" > /dev/null

python "$here/python/phylogenies/subsample_prokaryotes.py" "$out/$name.prok.tblout" "$prok" > "$out/$name.prok.fasta"
grep -v '^#' "$out/$name.euk.tblout" | awk '{print $1}' | sort -u > "$out/$name.euk.ids"
seqkit grep -f "$out/$name.euk.ids" "$euk" > "$out/$name.euk.fasta"
cat "$out/$name.euk.fasta" "$out/$name.prok.fasta" > "$out/$name.fasta"

n=$(grep -c '^>' "$out/$name.fasta")
if [ "$n" -le 1000 ]; then
  mode=${MAFFT_MODE:-linsi}
  "$mode" --thread "$threads" "$out/$name.fasta" > "$out/$name.$mode.fasta"
  python "$here/python/phylogenies/filter_gappy_columns.py" "$out/$name.$mode.fasta" 0.9 > "$out/$name.$mode.gt01.fasta"
  iqtree -s "$out/$name.$mode.gt01.fasta" -m MFP -mset LG -bb 1000 -nt AUTO -ntmax "$threads" \
         -pre "$out/$name.$mode.gt01"
else
  clustalo --threads "$threads" -i "$out/$name.fasta" -o "$out/$name.clustalo.fasta"
  python "$here/python/phylogenies/filter_gappy_columns.py" "$out/$name.clustalo.fasta" 0.99 > "$out/$name.clustalo.gt001.fasta"
  FastTree "$out/$name.clustalo.gt001.fasta" > "$out/$name.clustalo.gt001.fasttree"
fi
echo "done: $out"
