#!/usr/bin/env python3
"""Supplementary table of Boltz-2 predictions (Table S4) from the extracted predictions CSV.

    python python/boltz/make_table.py data/boltz2/boltz_predictions_full.csv results/tables/Table_S4_boltz2_predictions.xlsx
"""
import sys

import pandas as pd

LIGANDS = {"ach": "acetylcholine", "dop": "dopamine", "ser": "serotonin",
           "try": "tryptamine", "pea": "2-phenylethylamine", "his": "histamine"}
CLADE_A = ["XP_012561693.1", "t17782aep", "t21045aep", "XP_012561685.1", "XP_012561731.1"]
CLADE_B = ["t18253aep", "XP_012562020.1", "t1966aep", "XP_012554108.1", "XP_012561998.1", "XP_012561929.1", "t19057aep", "XP_012560064.1",
           "XP_012563056.1", "t25534aep", "XP_012562068.1", "XP_012562571.1", "XP_002170550.2"]
CONSERVED = {"XP_012561685.1": "ACh-rec1", "XP_012561731.1": "ACh-rec2", "t21045aep": "t21045aep (ACh-rec1/2 transcript)"}


def annotate(p):
    if p.endswith("_HUMAN"):
        species = "Homo sapiens"
    elif p.endswith("_MOUSE"):
        species = "Mus musculus"
    elif p.endswith("_DROME"):
        species = "Drosophila melanogaster"
    elif p.startswith(("XP_", "t")):
        species = "Hydra vulgaris"
    elif p == "MneR":
        species = "Mnemiopsis leidyi"
    elif p == "TriR":
        species = "Trichoplax adhaerens"
    else:
        species = ""
    if p.startswith("ACM") or p == "Q9VHW1_DROME":
        group = "muscarinic acetylcholine receptor (control)"
    elif p.startswith("HRH"):
        group = "histamine receptor (control)"
    elif p.startswith("5HT2"):
        group = "serotonin receptor (control)"
    elif p.startswith("DRD"):
        group = "dopamine receptor (control)"
    elif species == "Hydra vulgaris":
        group = "H. vulgaris mAChR-clade receptor"
    else:
        group = "non-bilaterian outgroup receptor"
    clade = "A" if p in CLADE_A else "B" if p in CLADE_B else ""
    return species, group, clade, "yes" if p in CONSERVED else "", CONSERVED.get(p, "")


def main(src, dst):
    d = pd.read_csv(src)
    ann = d["protein"].apply(annotate).apply(pd.Series)
    ann.columns = ["species", "receptor_group", "fig4B_hydra_subclade", "fully_conserved_orthosteric_site", "alias"]
    d.insert(1, "alias", ann["alias"])
    for i, c in enumerate(["species", "receptor_group", "fig4B_hydra_subclade", "fully_conserved_orthosteric_site"]):
        d.insert(2 + i, c, ann[c])
    d.insert(7, "ligand_name", d["ligand"].map(LIGANDS))
    d = d.sort_values(["receptor_group", "protein", "ligand"])
    desc = pd.DataFrame([
        ("protein", "Receptor identifier (UniProt entry name, NCBI RefSeq protein or Hydra transcriptome ID)"),
        ("alias", "Name used in the manuscript"),
        ("species", "Source species"),
        ("receptor_group", "Receptor class; bilaterian receptors of known specificity are controls"),
        ("fig4B_hydra_subclade", "Hydra mAChR subclade in Fig. 4B"),
        ("fully_conserved_orthosteric_site", "Orthosteric residues identical to bilaterian mAChRs (Fig. 4B)"),
        ("ligand / ligand_name / ligand_smiles", "Ligand code, name and SMILES used as Boltz-2 input"),
        ("protein_length", "Receptor sequence length (aa)"),
        ("affinity_pred_value", "Predicted affinity, log10(IC50 in uM); lower = stronger; mean of the two affinity heads"),
        ("affinity_probability_binary", "Predicted probability that the ligand is a binder; mean of the two heads"),
        ("affinity_*1 / affinity_*2", "Values of the individual affinity heads"),
        ("confidence_score", "Boltz-2 structure confidence (0.8 x complex pLDDT + 0.2 x ipTM)"),
        ("ptm / iptm / ligand_iptm", "Predicted TM-score of the complex, interface and ligand interface"),
        ("complex_plddt / complex_iplddt", "Mean pLDDT of the complex and of the interface"),
        ("complex_pde / complex_ipde", "Mean predicted distance error of the complex and of the interface"),
    ], columns=["column", "description"])
    with pd.ExcelWriter(dst) as xl:
        d.to_excel(xl, sheet_name="Boltz-2 predictions", index=False)
        desc.to_excel(xl, sheet_name="Column descriptions", index=False)
    print(f"{dst}: {len(d)} predictions, {d['protein'].nunique()} receptors x {d['ligand'].nunique()} ligands")


if __name__ == "__main__":
    main(*sys.argv[1:3])
