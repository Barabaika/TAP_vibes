# TAP_vibes

**Paired VH/VL sequences → ABodyBuilder2 → refinement → five TAP2 descriptors:** `L_tot`, `PSH`, `PPC`, `PNC`, and `SFvCSP`. Existing IMGT-numbered Fv structures can also be scored directly.

Includes regression tests against real responses from the [OPIG TAP web service](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap): **10 experimental humanized antibody pairs**, their sequences, original server PDBs and HTML, and SHA256 checksums. No cluster-specific paths, credentials, or experiment framework are required.

## Setup

Requirements: **Linux x86_64** (or WSL2 on Windows) and [Miniforge](https://github.com/conda-forge/miniforge) or [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) on `PATH`. No sudo is required. The bundled downloader uses a Linux x86_64 PSA executable.

```bash
git clone https://github.com/Barabaika/TAP_vibes.git
cd TAP_vibes

# Creates .tap-env, installs Python 3.10, HMMER, packages, and downloads assets.
bash setup.sh --device cpu
conda activate ./.tap-env

# For NVIDIA GPUs, choose the CUDA 12.1 PyTorch build instead:
# bash setup.sh --device cuda
# conda activate ./.tap-env
```

If using micromamba, activate with `micromamba activate ./.tap-env`. Alternatively, prefix any command with `conda run --prefix .tap-env` or `micromamba run --prefix .tap-env`; shell activation is not required.

`setup.sh` installs pinned scientific packages, checks dependency consistency, downloads PSA and the **four original ABodyBuilder2 weight files** (~674 MB), and verifies their SHA256. Assets are stored in `assets/` and excluded from Git. It does not submit sequences to TAP or start folding.

For scoring/testing existing structures without the folding dependencies or weights:

```bash
bash setup.sh --metrics-only
conda run --prefix .tap-env python -m pytest -q --run-structure
```

Use `--prefix` and `--assets` to select different local directories. `bash setup.sh --help` lists options. Dependencies are declared in [pyproject.toml](pyproject.toml); [prepare_assets.py](prepare_assets.py) pins the PSA commit and all weight checksums. `tap-assets --psa-only` downloads just PSA; `tap-assets` downloads both PSA and weights. Existing files with invalid checksums cause an error.

## Fold sequences and calculate all metrics

After activating the environment, run from the repository directory:

```bash
# The ten reference pairs, CPU, four threads:
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=4 OPENMM_CPU_THREADS=4 \
tap-fold --input examples/huab348_ten.csv --kind humanized --limit 10 \
  --device cpu --threads 4 --profile web-compatible \
  --weights assets/weights --psa assets/psa --output outputs/huab348

# On GPU, use --device cuda and do not set CUDA_VISIBLE_DEVICES=''.
```

Input CSV formats: `type,name,h_seq,l_seq` or `Specific,name,hseq,lseq`. For generated `sample_humanization_result.csv` files, select `--kind humanization`. Heavy/light chains always come from the same row. Only the 20 standard amino acids and complete variable domains without additional flanking residues are accepted. Duplicate pairs are removed using a hash of both chains. `--kind ''` disables type filtering; `--limit 0` removes the local calculation limit. `--pairs-json` accepts the format in `tests/fixtures/pairs.json`.

Outputs:

- `metrics.csv`: the five descriptors for every pair;
- `local_results.json`: descriptors, CDR lengths, surface charges, salt bridges, residue annotations, and provenance;
- `<id>/model.pdb`: refined structure, alongside unrefined PDBs, IMGT numbering, and predicted errors;
- `provenance.json`: package versions, settings, code/PSA/weight hashes.

Cached pairs are reused only when sequences and provenance match. Use a new output directory after changing code or protocol. Failed refinement stops the run while preserving completed pairs.

## Score an existing structure

```bash
tap-score tests/fixtures/web/000_9f7bf0ec0890/web_model.pdb \
  --psa assets/psa --profile web-compatible --output outputs/one.json
```

The PDB must contain `H` and `L` chains with valid IMGT numbering. ANARCI checks numbering and chain identity. Arbitrarily numbered structures must be renumbered first.

## Published protocol versus observed web behavior

`--profile paper` is the default. It follows the published [TAP2 description](https://www.nature.com/articles/s42003-023-05744-8): a 4.5 Å CDR vicinity radius and salt-bridge charge neutralization. It returns flags using the paper's [Table 1 thresholds](https://www.nature.com/articles/s42003-023-05744-8/tables/1).

`--profile web-compatible` reproduces behavior **observed on these ten antibodies on 2026-09-14**: a 4.0 Å radius, glycine hydrophobicity for salt-bridged residues, and unneutralized charges in PPC/PNC/SFvCSP. This is empirically validated compatibility, not official TAP source code or proof of equivalence for every antibody. The opt-in live test can detect subsequent server changes. See [protocol details and sources](docs/protocol.md).

Actual server flags are parsed from result logs. Inconsistent thresholds on the service homepage are not substituted for the flags returned with each result.

## Tests against the web service

```bash
# Formula/input/cache tests and integrity of the archived server responses:
python -m pytest -q

# Recompute ALL five metrics for ALL ten archived server PDBs:
python -m pytest -q --run-structure

# Explicit live test: ONE new request for the first of the same ten pairs.
# Takes several minutes; never runs in normal tests or GitHub Actions.
python -m pytest -q -s tests/test_live_web.py --run-live
```

Expected values come from real server responses, not from the implementation under test. Each metric uses absolute tolerance `1e-4` and relative tolerance `0`; this allows for server rounding. CDR lengths match exactly. An explicitly enabled structural test fails if PSA/HMMER is unavailable. Structural and live tests are skipped unless their respective flags are provided.

The live test compares metrics on the **same PDB** returned by the new service job; it does not assert that independent folding generates identical structures. Pytest retains the response in its temporary test directory. For up to ten sequential submissions with persistent output:

```bash
tap-web --input examples/huab348_ten.csv --kind humanized --limit 10 \
  --output outputs/new-reference --submit-web
```

Without `--submit-web`, this command only prepares `pairs.json`. An uncertain POST is not automatically retried: inspect `submission_started.json` and recover the existing job URL before proceeding. A recorded URL allows polling to resume without resubmission.

## Original validation results

On shared server PDBs, **50/50 values** matched in `web-compatible` mode within `1e-4`; maximum absolute error was `5.7398e-5`. After independent ABodyBuilder2 folding, MAE was: L_tot `0`, PSH `3.41928`, PPC `0.00115747`, PNC `0.0148167`, SFvCSP `1.23`. The second comparison has no arbitrary equality threshold: different structures/refinement can change surface descriptors.

[Summary statistics](docs/validation_summary.json) · [All comparisons](docs/validation_comparison.csv) · [Dataset and validation details](docs/validation.md).

To compare a new local folding run with the archive:

```bash
python compare.py --web tests/fixtures/web --local outputs/huab348 \
  --psa assets/psa --output outputs/comparison
```

This benchmark expects the ten reference pairs. It saves both protocols on both structure sets, server flags, and MAE/RMSE/max errors. For other datasets use `tap-fold` and `tap-score`.

## CPU refinement compatibility

ImmuneBuilder 1.2 passes a set `{'Threads', str(n_threads)}` instead of a dictionary to OpenMM in its CPU strained-bond repair branch. [tap2_compat.py](tap2_compat.py) corrects exactly this literal in a private in-memory module. The installed dependency and geometry checks remain intact. The correction and upstream module SHA256 are recorded in provenance.

PSA is provided by OPIG TNP under BSD-3-Clause; its license is downloaded with the executable. Model weights and dependencies retain their own licenses.
