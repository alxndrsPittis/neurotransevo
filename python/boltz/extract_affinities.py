#!/usr/bin/env python3
"""Collect Boltz-2 predictions (affinity + structure confidence) into one table.

Usage:
    python python/boltz/extract_affinities.py <Output dir> [--inputs <YAML dir>] > predictions.csv

Expected layout (boltz predict):
    Output/boltz_results_<protein>.vs.<ligand>/predictions/<protein>.vs.<ligand>/
        affinity_<protein>.vs.<ligand>.json
        confidence_<protein>.vs.<ligand>_model_0.json
Input YAMLs (<protein>.vs.<ligand>.yaml) provide the ligand SMILES and protein length.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

AFFINITY = ["affinity_pred_value", "affinity_probability_binary",
            "affinity_pred_value1", "affinity_probability_binary1",
            "affinity_pred_value2", "affinity_probability_binary2"]
CONFIDENCE = ["confidence_score", "ptm", "iptm", "ligand_iptm",
              "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]


def read_yaml_fields(path):
    """Minimal parser for the Boltz input YAML: protein sequence and ligand SMILES."""
    seq, smiles = None, None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("sequence:"):
            seq = s.split(":", 1)[1].strip()
        elif s.startswith("smiles:"):
            smiles = s.split(":", 1)[1].strip().strip("'\"")
    return seq, smiles


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output_dir", type=Path)
    ap.add_argument("--inputs", type=Path, help="directory with the input YAML files")
    args = ap.parse_args()

    fields = ["protein", "ligand", "ligand_smiles", "protein_length"] + AFFINITY + CONFIDENCE
    writer = csv.DictWriter(sys.stdout, fieldnames=fields)
    writer.writeheader()
    n = 0
    for run in sorted(args.output_dir.glob("boltz_results_*")):
        name = run.name.replace("boltz_results_", "")
        if ".vs." not in name:
            print(f"WARNING: cannot parse {run.name}", file=sys.stderr)
            continue
        protein, ligand = name.split(".vs.")
        pred = run / "predictions" / name
        aff_file = pred / f"affinity_{name}.json"
        if not aff_file.exists():
            print(f"WARNING: no affinity JSON for {name}", file=sys.stderr)
            continue
        row = {"protein": protein, "ligand": ligand}
        row.update({k: json.loads(aff_file.read_text()).get(k) for k in AFFINITY})
        conf_file = pred / f"confidence_{name}_model_0.json"
        if conf_file.exists():
            conf = json.loads(conf_file.read_text())
            row.update({k: conf.get(k) for k in CONFIDENCE})
        if args.inputs and (args.inputs / f"{name}.yaml").exists():
            seq, smiles = read_yaml_fields(args.inputs / f"{name}.yaml")
            row["ligand_smiles"] = smiles
            row["protein_length"] = len(seq) if seq else None
        writer.writerow(row)
        n += 1
    print(f"Extracted {n} predictions", file=sys.stderr)


if __name__ == "__main__":
    main()
