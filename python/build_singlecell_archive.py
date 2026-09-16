"""Zenodo archive with the single-cell-derived data and annotations used by the neurotransevo repository.

Run from the manuscript root with the neurotransevo environment:
    python 06.Data/zenodo/build_singlecell_archive.py
"""
import zipfile
from pathlib import Path


ROOT = Path(".")
REPO = ROOT / "10.Github/neurotransevo"
LEGACY = ROOT / "09.Legacy"
OUT = ROOT / "06.Data/zenodo/neurotransevo_singlecell_annotations.zip"

README = """Single-cell-derived data and annotations — Pittis et al., "From promiscuity to specialization:
The origin of biogenic amine signaling is pre-bilaterian".

Companion to the code at https://github.com/alxndrsPittis/neurotransevo (unpack into data/ as described
in data/README.md). The raw single-cell atlases are public and are not redistributed.

eggnog/                      eggNOG-mapper v2.1.7 annotations of the proteomes of Hydra vulgaris, Nematostella
                             vectensis, Mnemiopsis leidyi, Spongilla lacustris, Mus musculus and Drosophila
                             melanogaster (GO terms used to define the GO:0009308 candidate pool)
proteomes/UP000000803_7227.fasta
                             UniProt Drosophila reference proteome used to map accessions to FlyBase symbols
families/                    curated gene-family members per species and curated row orders of published panels
matrices/<species>.tsv.gz    cell-type-averaged expression (UMI counts summed per cell type / number of cells);
                             rows = genes, columns = cell types
coregulons/<SPECIES>.members.tsv
                             coregulon members: seed flag, best Pearson r, minimum BH-adjusted p, best seed,
                             eggNOG preferred name and description
"""


def main():
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("singlecell/README.txt", README)
        for f in sorted((LEGACY / "02.Analysis/SingleCell/Clustering/GOannotations").glob("*.emapper.annotations")):
            z.write(f, f"singlecell/eggnog/{f.name}")
        z.write(LEGACY / "01.Data/UP000000803_7227.fasta", "singlecell/proteomes/UP000000803_7227.fasta")
        z.write(REPO / "data/families/neurotrans_family_members.melted.list",
                "singlecell/families/neurotrans_family_members.melted.list")
        for f in sorted((REPO / "data/families/row_orders").glob("*.txt")):
            z.write(f, f"singlecell/families/row_orders/{f.name}")
        # matrices exported from R (row names = gene IDs) into 06.Data/zenodo/_matrices_tsv/
        for f in sorted((ROOT / "06.Data/zenodo/_matrices_tsv").glob("*.tsv.gz")):
            z.write(f, f"singlecell/matrices/{f.name}", compress_type=zipfile.ZIP_STORED)
        for f in sorted((REPO / "results/coregulons").glob("*.members.tsv")):
            z.write(f, f"singlecell/coregulons/{f.name}")
    print(OUT, f"{OUT.stat().st_size / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
