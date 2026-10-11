# Zenodo deposition and code freeze

Two deposits back the paper:

| Deposit | DOI | Status |
|---|---|---|
| Data (alignments/trees, single-cell annotations + matrices + coregulons, Boltz-2 raw) | [10.5281/zenodo.22813072](https://doi.org/10.5281/zenodo.22813072) | **draft — DOI reserved, not yet published** (resolves 404 until published) |
| Code (this repository) | — | frozen at tag `v1.0.0`; Zenodo DOI minted via GitHub release (see B) |

---

## A. Publish the data deposit (required before the bioRxiv preprint goes live)

The manuscript, `README.md` and `data/README.md` already cite 10.5281/zenodo.22813072;
that link is dead until the draft is published.

1. Log in to zenodo.org → **My uploads** → open the draft with the reserved DOI
   `10.5281/zenodo.22813072`.
2. Upload the three archives (contents as described in `data/README.md`):
   - `neurotransevo_phylogenies.zip` — alignments and trees behind Figs 2–5, S1–S2
     (file names as in `data/phylogenies/Table_S5.phylogenies.tsv`, "deposited" column).
   - `neurotransevo_singlecell_annotations.zip` — `eggnog/`, `proteomes/`, `matrices/`,
     `coregulons/`.
   - `neurotransevo_boltz2_predictions.zip` — Boltz-2 inputs, MSAs, structures, affinity
     and confidence outputs (raw outputs are on the MacBook — see project README task list).
3. Metadata:
   - **Title:** `Data for: From promiscuity to specialization: The origin of cholinergic
     and monoaminergic signaling is pre-bilaterian` (mirror the final manuscript title if
     it changes).
   - **Authors:** exactly the final manuscript author list, with ORCIDs (copy from
     `CITATION.cff`). *If the author list changes before submission, fix it here, in the
     manuscript, and in `CITATION.cff` before publishing anything.*
   - **License:** CC BY 4.0. **Type:** Dataset.
   - **Related identifiers:** `https://github.com/alxndrsPittis/neurotransevo`
     (*is supplemented by*); after submission add the bioRxiv DOI (*is supplement to*).
4. Confirm the reserved DOI shown in the form is `10.5281/zenodo.22813072`
   (this is the version DOI cited in the paper — do not delete the draft and start a new
   one, or the reserved DOI is lost) → **Publish**. The DOI resolves within a few minutes.
5. Metadata can be edited after publishing without a new version (e.g. adding the
   bioRxiv DOI). Only file changes require a new version (which mints a new version DOI).

## B. Code freeze and archival

The repository is frozen for submission at the annotated tag **`v1.0.0`**
(created 2026-10-11; the analysis code was last modified in commit `2c1ab5b`,
2026-09-29). Nothing under this tag may change;
any post-freeze fix goes into a new commit and, if it affects the paper, a `v1.0.1` tag.

Recommended path (repo becomes public at preprint time — a preprint should not cite a
private repo):

1. On GitHub: **Settings → change visibility → public** (do this right before, or
   together with, the bioRxiv submission).
2. On zenodo.org: **Account → GitHub** → toggle `alxndrsPittis/neurotransevo` **on**
   *before* creating the release (the webhook only fires for releases created after
   enabling).
3. On GitHub: **Releases → Draft a new release** → choose existing tag `v1.0.0` →
   title `v1.0.0 — bioRxiv submission` → publish. Zenodo archives the release
   automatically and mints a concept DOI + version DOI (appears under My uploads →
   GitHub within minutes).
4. Add the minted code DOI to: `CITATION.cff` (`doi:` field), the README badge, and the
   manuscript's Code Availability statement; keep the GitHub URL alongside it.

Fallback (repo must stay private for now): upload the frozen archive
`neurotransevo-v1.0.0.zip` (built with `git archive v1.0.0`, staged in
`07.Submission/zenodo/`) as a normal Zenodo upload (Type: Software, license MIT,
same authors), reserve the DOI, cite it in Code Availability, publish. Switch to the
GitHub integration at journal publication if desired.

## C. Pre-submission checklist

- [ ] Author list final (Gaspar/Luis decision) and identical in: manuscript,
      `CITATION.cff`, README citation block, Zenodo data deposit, Zenodo code deposit.
- [ ] Data deposit published; `https://doi.org/10.5281/zenodo.22813072` resolves.
- [ ] Three data zips verified: unzip cleanly, file names match `Table_S5` /
      `data/README.md`.
- [ ] Repo public + release `v1.0.0` + code DOI minted (or fallback deposit published).
- [ ] Manuscript Data/Code Availability quotes both DOIs and the GitHub URL.
- [ ] After bioRxiv acceptance: add the preprint DOI to both Zenodo records.
