"""Build the Zenodo archive of Boltz-2 inputs and predictions."""
import sys
import zipfile
from pathlib import Path

BASE = Path(sys.argv[1] if len(sys.argv) > 1 else "06.Data")   # manuscript data folder (usage: python build_zenodo_archive.py <06.Data>)
RAW = BASE / "boltz2/Raw/Neurotransmission"
OUT = BASE / "zenodo/neurotransevo_boltz2_predictions.zip"
OUT.parent.mkdir(exist_ok=True)

README = """Boltz-2 receptor-ligand predictions — Pittis et al., "From promiscuity to specialization:
The origin of biogenic amine signaling is pre-bilaterian".

46 receptors x 6 ligands (acetylcholine, dopamine, serotonin, tryptamine, 2-phenylethylamine,
histamine) = 276 predictions, Boltz 2.2.1, default settings (`boltz predict <yaml> --use_msa_server`),
one affinity prediction per complex.

inputs/YAMLs/              Boltz input files (receptor sequence, ligand SMILES, affinity property)
inputs/YAMLs_with_msa/     inputs of the 70 re-run jobs, pointing to the precomputed MSA of the receptor
msa/<receptor>.csv         MSA of each receptor (ColabFold MMseqs2 server; identical for all ligands)
predictions/<receptor>.vs.<ligand>/
    <name>_model_0.cif                 predicted complex structure
    affinity_<name>.json               affinity_pred_value (log10 IC50, uM), affinity_probability_binary,
                                       and the values of the two affinity heads
    confidence_<name>_model_0.json     confidence_score, pTM, ipTM, ligand ipTM, pLDDT, PDE
    pae_/pde_/plddt_*.npz              per-token confidence arrays
    pre_affinity_<name>.npz            affinity-module input
scripts/run_missing_boltz.sh, scripts/missing_jobs.txt   batch re-run of the 70 missing jobs
boltz_predictions_full.csv  all predictions in one table (also Supplementary Table S4)

Code: https://github.com/alxndrsPittis/neurotransevo (python/boltz/)
"""

n = 0
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    z.writestr("boltz2/README.txt", README)
    z.write(BASE / "boltz2/boltz_predictions_full.csv", "boltz2/boltz_predictions_full.csv")
    z.write(RAW / "run_missing_boltz.sh", "boltz2/scripts/run_missing_boltz.sh")
    z.write(RAW / "missing_jobs.txt", "boltz2/scripts/missing_jobs.txt")
    for sub in ("YAMLs", "YAMLs_with_msa"):
        for f in sorted((RAW / "Input" / sub).glob("*.yaml")):
            z.write(f, f"boltz2/inputs/{sub}/{f.name}")
    seen = set()
    for run in sorted((RAW / "Output").glob("boltz_results_*")):
        name = run.name.replace("boltz_results_", "")
        protein = name.split(".vs.")[0]
        for f in sorted((run / "predictions" / name).iterdir()):
            z.write(f, f"boltz2/predictions/{name}/{f.name}")
            n += 1
        if protein not in seen:
            msa = run / "msa" / f"{name}_0.csv"
            if msa.exists():
                z.write(msa, f"boltz2/msa/{protein}.csv")
                seen.add(protein)
print(f"{OUT}: {n} prediction files, {len(seen)} MSAs, {OUT.stat().st_size / 1e6:.0f} MB")
