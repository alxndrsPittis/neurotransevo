# Data

Small curated inputs are committed here. Large public datasets are downloaded by the
user; their location is set in `config.yml` / `config.local.yml`.

## Committed

| Path | Content |
|---|---|
| `families/neurotrans_family_members.melted.list` | Curated gene-family members per species (family, module, display name, gene ID, species), derived from the phylogenies. Defines the coregulon seeds and the rows of all family heatmaps. |
| `families/row_orders/` | Hand-curated row orders used in published panels (Fig. 5C, 5D, S4). |
| `families/taxon_colors.txt` | Taxon colour scheme used in the phylogenies. |
| `mappings/mouse_uniprot2ensembl.tab` | UniProt accession → Ensembl gene (mouse). |
| `mappings/mouse_ensembl2symbol.tab` | Ensembl gene → gene symbol (Cao et al. 2019 gene annotation). |
| `mappings/spongilla_gene2protein.tab` | Trinity gene → protein ID (Musser et al. 2021). |
| `boltz2/boltz_predictions_full.csv` | Boltz-2 predictions with structure-confidence metrics (Table S4); `boltz2/inputs/` holds the input YAML files. |
| `phylogenies/Table_S5.phylogenies.tsv` | Alignments and trees behind the figures (Table S5), with original and deposited file names. |
| `boltz2/boltz_affinities_v2.csv` | Boltz-2 affinity predictions: 46 receptors × 6 ligands (ACh, dopamine, serotonin, tryptamine, 2-phenylethylamine, histamine), both prediction heads and their mean. |

## To download (not redistributed)

Place each dataset under `data/singlecell/<Study>/` (or point `paths.singlecell` at an existing copy) with the file names below.

| Directory | Species | Files used | Source |
|---|---|---|---|
| `Cao2019/` | *Mus musculus* | `gene_count_cleaned.RDS`, `cell_annotate.csv`, `gene_annotate.csv` | Cao et al. 2019, *Nature* — Mouse Organogenesis Cell Atlas (GEO GSE119945) |
| `Li2021/` | *Drosophila melanogaster* | `s_fca_biohub_body_10x.loom` | Li et al. 2022, *Science* — Fly Cell Atlas, 10x body (flycellatlas.org) |
| `Siebert2019/` | *Hydra vulgaris* | `wt_expression_matrix.txt` | Siebert et al. 2019, *Science* |
| `Arnau2021/Hydra/` | *Hydra vulgaris* | `Hvul_cell_type_assignments` | Levy et al. 2021, *Cell* — github.com/sebe-pedros-lab/Stylophora_single_cell_atlas |
| `Arnau2021/Nematostella/` | *Nematostella vectensis* | `Nvec_adult_sc_UMI_counts.RDS`, `Nvec_adult_cell_type_assignments` | Levy et al. 2021, *Cell* (same repository) |
| `Arnau2018/` | *Mnemiopsis leidyi* | `GSM3021563_Mnemiopsis_leidyi_UMI_table.txt`, `..._metacell_definition.txt`, `..._metacell_assignments.txt` | Sebé-Pedrós et al. 2018, *Nat. Ecol. Evol.* (GEO GSM3021563) |
| `Musser2020/` | *Spongilla lacustris* | `GSE134912_spongilla_10x_count_matrix.txt`, `clusters_42_34_and_oddcells_reassigned.txt`, `Table1_celltype_descriptions_final.tsv` | Musser et al. 2021, *Science* (GEO GSE134912) |

Other inputs:

| Config key | Content | Source |
|---|---|---|
| `paths.eggnog` | eggNOG-mapper v2.1.7 annotation tables (`*.emapper.annotations`) of the six proteomes | Zenodo ([10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072)): `neurotransevo_singlecell_annotations.zip` (`eggnog/`) |
| `paths.proteomes` | `UP000000803_7227.fasta` (Drosophila reference proteome; used to map UniProt accessions to FlyBase symbols) | Zenodo ([10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072)): `neurotransevo_singlecell_annotations.zip` (`proteomes/`) |
| `paths.boltz_raw` | Boltz-2 inputs, MSAs, structures, affinity and confidence files | Zenodo ([10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072)): `neurotransevo_boltz2_predictions.zip` |
| — | Alignments and trees behind Figs 2–5 and S1–S2 (Table S5) | Zenodo ([10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072)): `neurotransevo_phylogenies.zip` |
| — | Cell-type-averaged matrices and coregulon tables (to skip rebuilding from the atlases) | Zenodo ([10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072)): `neurotransevo_singlecell_annotations.zip` (`matrices/`, `coregulons/`) |
