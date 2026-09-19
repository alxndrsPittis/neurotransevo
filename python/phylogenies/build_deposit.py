"""Collect the trees and alignments behind the figures under numbered names, write Table S5 and the Zenodo zip.

Run from the manuscript root:  python 06.Data/phylogenies/build_deposit.py
Original files stay in 09.Legacy/02.Analysis; copies go to 06.Data/phylogenies/deposit/.
"""
import re
import shutil
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(".")
ANALYSIS = ROOT / "09.Legacy/02.Analysis"
DEPOSIT = ROOT / "06.Data/phylogenies/deposit"
ZIP = ROOT / "06.Data/zenodo/neurotransevo_phylogenies.zip"

# (id, figure, family, step, alignment basename, tree basename, status)
ENTRIES = [
    ("01", "2B", "carnitine_acyltransferases", "family", "Carn_acyltransf.vs.fusion.linsi.gt01", "Carn_acyltransf.vs.fusion.linsi.gt01.treefile", ""),
    ("02", "2B, S2", "carnitine_acyltransferases", "ChAT-CrAT_clade", "Carn_acyltransf.vs.fusion.linsi.gt01.treefile.CholineAcetyltransf_clade.linsi.gt01", "Carn_acyltransf.vs.fusion.linsi.gt01.treefile.CholineAcetyltransf_clade.linsi.gt01.treefile", ""),
    ("03", "2C", "biopterin_hydroxylases", "family", "06.Data/phylogenies/regenerated/Biopterin_H.vs.fusion.mft", "Biopterin_H.vs.fusion.mft.fasttree", ""),
    ("04", "2C", "biopterin_hydroxylases", "AAAH_clade", "Biopterin_H.vs.fusion.mft.fasttree_selected_NoLegionella.linsi.gt01", "Biopterin_H.vs.fusion.mft.fasttree_selected_NoLegionella.linsi.gt01.treefile", ""),
    ("05", "2D", "PLP_decarboxylases", "family", "Pyridoxal_deC.vs.fusion.mft", "Pyridoxal_deC.vs.fusion.mft.fasttree", ""),
    ("06", "2D", "PLP_decarboxylases", "AADC_clade", "Pyridoxal_deC.vs.fusion.mft.fasttree_selected.linsi.gt01", "Pyridoxal_deC.vs.fusion.mft.fasttree_selected.linsi.gt01.treefile", ""),
    ("07", "2D", "PLP_decarboxylases", "GAD_clade", "Pyridoxal_deC.vs.fusion.mft.fasttree_selected_2_GABAandProk_wholeclade.linsi.gt01", "Pyridoxal_deC.vs.fusion.mft.fasttree_selected_2_GABAandProk_wholeclade.linsi.gt01.treefile", ""),
    ("08", "2E", "DOMON_monooxygenases", "family", "DOMON.vs.fusion.linsi.gt01", "DOMON.vs.fusion.linsi.gt01.treefile", ""),
    ("09", "2E", "DOMON_monooxygenases", "DBH-MOXD_clade", "DOMON.vs.fusion.linsi.gt01.treefile_selected.linsi.gt01", "DOMON.vs.fusion.linsi.gt01.treefile_selected.linsi.gt01.treefile", ""),
    ("10", "3A", "MFS_superfamily", "family", "MFS_1.PF07690.16.vs.jake_cleaned_proteome_files_69.mft.gt01", "MFS_1.PF07690.16.vs.jake_cleaned_proteome_files_69.mft.gt01.fasttree", ""),
    ("11", "3A", "MFS_superfamily", "SLC18_MCL-I1.4-cluster0005", "MCL_I14.0005.gt01.linsi", "MCL_I14.0005.gt01.linsi.treefile", ""),
    ("12", "3A", "MFS_superfamily", "SLC22_MCL-I1.4-cluster0001", "MCL_I14.0001.SLC22.clustalo.gt01", "MCL_I14.0001.SLC22.clustalo.gt01.fasttree", ""),
    ("13", "3A", "MFS_superfamily", "SLC22_cluster0001_subset", "MCL_I14.0001.SLC22.subsampled4_extended.einsi.gt01", "MCL_I14.0001.SLC22.subsampled4_extended.einsi.gt01.treefile", ""),
    ("14", "3B", "SSF_transporters", "family", "SSF.vs.fusion.mft", "SSF.vs.fusion.mft.fasttree", ""),
    ("15", "3B", "SSF_transporters", "ChT_clade", "SSF.vs.fusion.mft.fasttree_selected.linsi.gt01", "SSF.vs.fusion.mft.fasttree_selected.linsi.gt01.treefile", ""),
    ("16", "4A", "rhodopsin_GPCRs", "MCL-I6.0-cluster0002", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.mft.gt01", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.mft.gt01.fasttree", "source cluster; not drawn"),
    ("17", "4A", "rhodopsin_GPCRs", "cluster0002_4-10TM-filtered", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.filtered_4to10.clustalo.gt001", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.filtered_4to10.clustalo.gt001.fasttree", "tree from which the H-ACh-monoamine clade was cut; not drawn"),
    ("18", "4A", "rhodopsin_GPCRs", "H-ACh-monoamine_clade", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.filtered_4to10.clustalo.gt001.fasttree.clade_selection.clustalo.gt001", "7tm_1.PF00001.21.vs.jake_cleaned_proteome_files_70.allVSall.I60.0002.filtered_4to10.clustalo.gt001.fasttree.clade_selection.clustalo.gt001.treefile", ""),
    ("19", "5A", "Cys-loop_LGICs", "family", "Neur_chan_nicotinic_receptors.vs.fusion.clustalo.gt01", "Neur_chan_nicotinic_receptors.vs.fusion.clustalo.gt01.lg_fasttree", ""),
    ("20", "5B", "Cys-loop_LGICs", "cationic_clade_subset", "Neur_chan_nicotinic_receptors.vs.fusion.subsampled3.einsi.gt01", "Neur_chan_nicotinic_receptors.vs.fusion.subsampled3.einsi.gt01.treefile", ""),
    ("21", "S1A", "NNMT-PNMT-TEMT", "family", "NNMT_PNMT_TEMT.PF01234.17.fusion.linsi.gt01", "NNMT_PNMT_TEMT.PF01234.17.fusion.linsi.gt01.treefile", ""),
    ("22", "S1B", "cytochrome_P450", "family", "p450.PF00067.VS.fusion.clustalo.gt001", "p450.PF00067.VS.fusion.clustalo.gt001.fasttree", ""),
    ("23", "S1B", "cytochrome_P450", "CYP2_clade_subset", "p450.PF00067.VS.fusion.clustalo.gt001.fasttree_clade.subclade.subsampled5_extended.einsi.gt01", "p450.PF00067.VS.fusion.clustalo.gt001.fasttree_clade.subclade.subsampled5_extended.einsi.gt01.treefile", ""),
]

index = {}
for p in ANALYSIS.rglob("*"):
    if p.is_file() and "@eaDir" not in p.parts and "iTOLannotations" not in p.parts:
        index.setdefault(p.name, []).append(p)


def find(name):
    if "/" in name:
        return ROOT / name
    hits = sorted(index.get(name, []), key=lambda p: ("FCN/TakeB" not in str(p), -p.stat().st_size))
    hits = [h for h in hits if h.stat().st_size > 0]
    return hits[0] if hits else None


def fasta_dims(p):
    n, cols, cur = 0, 0, 0
    with open(p, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if n == 1:
                    cols = cur
                n, cur = n + 1, 0
            else:
                cur += len(line.strip())
    return n, (cols if n > 1 else cur)


def program(name):
    """Aligner of the final step: the last aligner token in the file name."""
    tokens = re.findall(r"\.(linsi|einsi|clustalo|mft)(?=\.|$)", name)
    return {"linsi": "MAFFT L-INS-i", "einsi": "MAFFT E-INS-i", "clustalo": "Clustal Omega",
            "mft": "MAFFT (FFT-NS-2)"}.get(tokens[-1], "") if tokens else ""


def gap_filter(name):
    """Gap filter of the final alignment: the last gt01/gt001 token (the step order may be swapped in the name)."""
    last_aln = max((m.end() for m in re.finditer(r"\.(linsi|einsi|clustalo|mft)(?=\.|$)", name)), default=0)
    tokens = re.findall(r"\.(gt001|gt01)(?=\.|$)", name)
    after = re.findall(r"\.(gt001|gt01)(?=\.|$)", name[last_aln:])
    tok = after[-1] if after else (tokens[-1] if tokens else None)
    return {"gt001": ">99% gaps", "gt01": ">90% gaps"}.get(tok, "none")


if DEPOSIT.exists():
    shutil.rmtree(DEPOSIT)   # rebuilt from the originals each run
DEPOSIT.mkdir(parents=True)
rows = []
for eid, fig, family, step, aln_name, tree_name, status in ENTRIES:
    stem = f"{eid}_Fig{fig.split(',')[0].replace(' ', '')}_{family}_{step}"
    tree = find(tree_name)
    aln = find(aln_name)
    iq = find(re.sub(r"\.treefile$", ".iqtree", tree_name)) if tree_name.endswith(".treefile") else None
    unaligned = False
    aln_out = DEPOSIT / f"{stem}.{'unaligned' if unaligned else 'alignment'}.fasta"
    tree_out = DEPOSIT / f"{stem}.tree.nwk"
    shutil.copyfile(aln, aln_out)
    shutil.copyfile(tree, tree_out)
    leaves = len(re.findall(r"[(,]\s*[^():,;\s]+\s*:", tree.read_text(errors="ignore")))
    n, cols = fasta_dims(aln)
    model = ""
    if iq:
        shutil.copyfile(iq, DEPOSIT / f"{stem}.iqtree.txt")
        m = re.search(r"Best-fit model according to BIC:\s*(\S+)", iq.read_text(errors="ignore"))
        model = m.group(1) if m else ""
    method = "IQ-TREE 1.6.12 (ModelFinder, LG models; ultrafast bootstrap)" if tree_name.endswith(".treefile") else "FastTree 2.1.11 (-lg)"
    rows.append({
        "ID": eid, "Figure": fig, "Protein family": family.replace("_", " "), "Tree": step.replace("_", " "),
        "Sequences": n if not unaligned else n, "Alignment columns": "" if unaligned else cols, "Tree leaves": leaves,
        "Alignment method": "" if unaligned else program(aln_name), "Gap filter": "" if unaligned else gap_filter(aln_name),
        "Tree method": method, "Model (BIC)": model,
        "Alignment file": aln_out.name, "Tree file": tree_out.name,
        "Original alignment": str(aln.relative_to(ANALYSIS)) if ANALYSIS in aln.parents else "FCN/TakeB/Seqs/Biopterin_H.vs.fusion.mft (re-aligned from Biopterin_H.vs.fusion.fasta, MAFFT v7.505)", "Original tree": str(tree.relative_to(ANALYSIS)),
        "Note": status,
    })
    print(f"{eid} {fig:7s} {step:32s} seqs={n:<5} cols={cols:<6} leaves={leaves:<5} {model:10s} {status}")

table = pd.DataFrame(rows)
PROVENANCE = ["Original alignment", "Original tree", "Note"]          # working paths: internal only
table[["ID", "Figure", "Protein family", "Tree"] + PROVENANCE].to_csv(
    ROOT / "06.Data/phylogenies/Table_S5.provenance_internal.tsv", sep="\t", index=False)
table = table.drop(columns=PROVENANCE)
table.to_csv(ROOT / "06.Data/phylogenies/Table_S5.phylogenies.tsv", sep="\t", index=False)
table.to_excel(ROOT / "03.Tables/supplementary/Table_S5.xlsx", sheet_name="Phylogenies", index=False)

readme = ("Trees and alignments of the gene-family phylogenies — Pittis et al.\n\n"
          "Files are named <ID>_<Figure>_<family>_<tree>. Each entry has the alignment (FASTA), the tree (Newick)\n"
          "and, for IQ-TREE runs, the IQ-TREE report.\n"
          "Table_S5.phylogenies.tsv lists sequences, alignment columns, methods, models and the deposited file names.\n"
          "Family trees were computed with FastTree (-lg) or IQ-TREE; clades of interest were extracted, re-aligned and\n"
          "recomputed with IQ-TREE. MFS and rhodopsin GPCR superfamilies were first clustered with BLASTP + MCL.\n")
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("phylogenies/README.txt", readme)
    z.write(ROOT / "06.Data/phylogenies/Table_S5.phylogenies.tsv", "phylogenies/Table_S5.phylogenies.tsv")
    for f in sorted(DEPOSIT.iterdir()):
        z.write(f, f"phylogenies/{f.name}")
print(ZIP, f"{ZIP.stat().st_size / 1e6:.0f} MB")
