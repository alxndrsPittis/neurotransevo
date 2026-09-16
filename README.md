# neurotransevo

Code to reproduce the analyses and figures of

> Pittis A.A., Yañez-Guerra L.A., Ruperti F., Musser J.M., Cole A.G., Marinković M., Thiel D., Huerta-Cepas J., Technau U., Jékely G., Arendt D.
> **From promiscuity to specialization: The origin of biogenic amine signaling is pre-bilaterian.**

A concise overview of every analysis, with the resulting panels, is in [`docs/index.html`](docs/index.html).

## What is here

| Analysis | Figures | Code |
|---|---|---|
| Cell-type-averaged expression matrices for six single-cell atlases (mouse, fly, *Hydra*, *Nematostella*, *Mnemiopsis*, *Spongilla*) | all expression panels | `R/matrices.R`, `scripts/01_build_matrices.R` |
| Amine-metabolism coregulons (Pearson r ≥ 0.5, BH-adjusted p ≤ 0.05, intersected with GO:0009308) | 8, S11–S14 | `R/coregulons.R`, `scripts/02_compute_coregulons.R`, `scripts/03_plot_coregulons.R` |
| Gene-family expression heatmaps | 4D–E, 5C–D, 7, S4–S10 | `R/heatmaps.R`, `scripts/04_plot_family_heatmaps.R` |
| Boltz-2 receptor–ligand binding predictions | 4B, S15 | `python/boltz/` |
| Gene-family phylogenies (HMMER, BLAST+MCL, MAFFT/Clustal Omega, IQ-TREE/FastTree) and tree plotting | 2, 3, 4A–B, 5A–B, S1–S2 | `scripts/phylogeny_pipeline.sh`, `scripts/mcl_clusters.sh`, `python/phylogenies/` |

All expression heatmaps share one style (`R/heatmaps.R`): cell types ordered by complete-linkage clustering of
the whole transcriptome (1 − Pearson r), neuronal cell types in a right-hand block, row z-scores.
SVG and PDF outputs keep text editable.

## Quick start

```bash
git clone https://github.com/alxndrsPittis/neurotransevo.git
cd neurotransevo
conda env create -f environment.yml
conda activate neurotransevo

# point the pipeline at the downloaded data (see data/README.md)
cp config.local.yml.example config.local.yml   # then edit the paths

make all          # matrices -> coregulons -> figures -> summary tables
```

Individual steps:

```bash
Rscript scripts/01_build_matrices.R hydra nematostella   # selected datasets
Rscript scripts/02_compute_coregulons.R HYDRA            # selected species
Rscript scripts/03_plot_coregulons.R
Rscript scripts/04_plot_family_heatmaps.R Fig7B_hydra_enzymes_transporters
python  python/boltz/plot_panels.py
```

Outputs are written to `results/` (`matrices/`, `coregulons/`, `figures/`, `tables/`).

The phylogenetic tools live in a separate environment (`envs/phylogenies.yml`); see the header of
`scripts/phylogeny_pipeline.sh`.

## Configuration

`config.yml` holds all paths (relative to the repository) and the coregulon thresholds.
Override paths per machine in a git-ignored `config.local.yml`, or with environment variables
`NEUROTRANSEVO_<KEY>` (e.g. `NEUROTRANSEVO_SINGLECELL=/data/singlecell`).

## Reproducibility notes

- The six expression matrices and the six coregulons produced by this code are identical (to
  floating-point precision) to those used for the manuscript figures.
- Column order, dataset and cell-type granularity are the same in every expression figure.
- Row orders that were curated by hand for publication (Fig. 5C, 5D, S4) are stored in
  `data/families/row_orders/`.

## License

MIT (see `LICENSE`). Third-party datasets keep their original licenses.
