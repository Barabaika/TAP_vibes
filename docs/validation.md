# Validation on ten real antibody pairs

The dataset contains the first ten unique `type=humanized` rows from `Antibody-AIR/data/antibody_eval_data/HuAb348_data/humanization_pair_data_filter.csv`. Original row numbers, heavy/light pairing, and SHA256 hashes are retained. Machine-specific source paths have been replaced with a relative dataset identifier.

Names: h3A3-5, AB240, 31HZ, 56HZ, 74HZ, h11E6.1, unnamed HuAb348 row 15, NuHuTac, hSC73.38, hSC73.39. These are humanized pairs from the experimental reference, not new model-generated candidates. Exactly ten unique antibodies are used.

## Archived service responses

TAP responses were retrieved on 2026-09-14. `tests/fixtures/web/web_results.json` records real job URLs, retrieval times, server metrics, and PDB SHA256 hashes. Each pair has an original `result.html` and `web_model.pdb`; Git preserves their bytes. `tests/fixtures/pairs.json` and `examples/huab348_ten.csv` provide the paired inputs.

Tests check sequence/PDB integrity, agreement of expectations with original HTML, and all five recomputed descriptors on each PDB. Absolute tolerance is 1e-4. Expectations are never regenerated automatically when the implementation changes.

## Two comparison levels

Shared PDBs, `web-compatible`: 50/50 values within tolerance; maximum error 5.7397933e-5. L_tot matches exactly; SFvCSP agrees to machine precision.

Independent folding MAE: L_tot=0, PSH=3.41928, PPC=0.00115747, PNC=0.0148167, SFvCSP=1.23. Maximum errors: PSH=9.13898 and SFvCSP=4.1. These describe the original CPU run, not deterministic expectations for new refinement. The [summary](validation_summary.json) and [200 comparisons](validation_comparison.csv) include both protocols on both structure sets.

Original configuration: Python 3.10, PyTorch 2.1.0+cu121 in CPU mode, ImmuneBuilder 1.2, OpenMM 8.2.0, PDBFixer 1.12.0, HMMER 3.4, four CPU threads. Nine structures completed before the ImmuneBuilder CPU typo was detected; the tenth completed after the targeted correction. Strained-bond checks were retained.

The live test requires `--run-live`, submits the first of the same ten pairs, and compares five metrics on the newly returned server PDB. It does not run in CI or introduce an eleventh unique antibody.

## Repository verification

The repository was installed into a new isolated environment using `bash setup.sh --metrics-only`. Dependency consistency passed; `python -m pytest -q --run-structure` reported **31 passed, 1 skipped**. The skipped test was the explicitly opt-in live request.

A separate live test passed against [a new TAP result](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap_results/20260914_0454734): **1 passed**, with all five descriptors checked on the newly returned PDB. This request reused h3A3-5 from the same ten-antibody set. The original archived expectations were preserved.

The complete CPU setup also passed: it installed the folding dependencies, downloaded all four weights with matching SHA256 hashes, and reported no broken requirements. A subsequent `tap-fold` smoke test on the first reference pair completed refinement and all five metrics in 15.7 seconds in that isolated environment (excluding setup and model initialization).
