#!/usr/bin/env python3
"""
Extract all Boltz-2 affinity predictions into a single CSV.

Usage:
    python3 python/boltz/extract_affinities.py /path/to/Output/ > data/boltz2/boltz_affinities_v2.csv

The script expects the directory structure:
    Output/boltz_results_PROTEIN.vs.LIGAND/predictions/PROTEIN.vs.LIGAND/affinity_PROTEIN.vs.LIGAND.json
"""

import json
import os
import sys
import csv

def extract_affinities(output_dir):
    """Walk through all boltz_results_* directories and extract affinity JSONs."""
    results = []
    
    for dirname in sorted(os.listdir(output_dir)):
        if not dirname.startswith("boltz_results_"):
            continue
        
        # Parse protein and ligand from directory name
        name = dirname.replace("boltz_results_", "")
        parts = name.split(".vs.")
        if len(parts) != 2:
            print(f"WARNING: Cannot parse {dirname}", file=sys.stderr)
            continue
        
        protein, ligand = parts
        
        # Find the affinity JSON
        json_path = os.path.join(
            output_dir, dirname, "predictions", name, f"affinity_{name}.json"
        )
        
        if not os.path.exists(json_path):
            print(f"WARNING: No affinity JSON found for {name}", file=sys.stderr)
            continue
        
        try:
            with open(json_path) as f:
                data = json.load(f)
            
            results.append({
                "protein": protein,
                "ligand": ligand,
                "affinity_pred_value": data.get("affinity_pred_value"),
                "affinity_probability_binary": data.get("affinity_probability_binary"),
                "affinity_pred_value1": data.get("affinity_pred_value1"),
                "affinity_probability_binary1": data.get("affinity_probability_binary1"),
                "affinity_pred_value2": data.get("affinity_pred_value2"),
                "affinity_probability_binary2": data.get("affinity_probability_binary2"),
            })
        except Exception as e:
            print(f"ERROR reading {json_path}: {e}", file=sys.stderr)
    
    return results


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 python/boltz/extract_affinities.py /path/to/Output/", file=sys.stderr)
        sys.exit(1)
    
    output_dir = sys.argv[1]
    results = extract_affinities(output_dir)
    
    # Write CSV to stdout
    writer = csv.DictWriter(sys.stdout, fieldnames=[
        "protein", "ligand",
        "affinity_pred_value", "affinity_probability_binary",
        "affinity_pred_value1", "affinity_probability_binary1",
        "affinity_pred_value2", "affinity_probability_binary2",
    ])
    writer.writeheader()
    for row in results:
        writer.writerow(row)
    
    print(f"Extracted {len(results)} predictions", file=sys.stderr)
