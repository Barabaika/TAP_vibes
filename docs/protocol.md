# Protocol and sources

- [TAP2: Raybould et al., Communications Biology 2024](https://www.nature.com/articles/s42003-023-05744-8), [Table 1](https://www.nature.com/articles/s42003-023-05744-8/tables/1).
- [Original TAP: Raybould et al., PNAS 2019](https://doi.org/10.1073/pnas.1810576116).
- [ImmuneBuilder / ABodyBuilder2](https://github.com/oxpig/ImmuneBuilder), [original weights](https://zenodo.org/records/7258553).
- [ANARCI](https://github.com/oxpig/ANARCI).
- [PSA from OPIG TNP](https://github.com/oxpig/TNP/tree/29dcac72f1380e8538e8870f45a699d3c6156162), commit `29dcac72f1380e8538e8870f45a699d3c6156162`, BSD-3-Clause. TNP's nanobody metric implementation is not used.

## Definitions

ABodyBuilder2 predicts a paired Fv with the four original weights. Models are ranked and refined with strained-bond checks. Metrics use the final PDB. Hydrogens are excluded before SASA. PSA calculates side-chain accessibility on the combined H/L structure, including inter-chain occlusion.

A surface residue has relative side-chain accessibility ≥7.5% of the Ala-X-Ala reference. IMGT CDR ranges are 27–38, 56–65, and 105–117 on each chain, including insertions. L_tot sums all six CDR lengths.

CDR vicinity starts with exposed CDR residues and two flanking IMGT anchor positions on either side. Exposed neighbors are added when their minimum heavy-atom distance is below the profile's radius. Expansion is a single step, not transitive.

PSH sums `H_i H_j / d_ij²` over **ordered** pairs of different vicinity residues with `d_ij < 7.5 Å`. `H_i = 1 + (KD_i + 4.5)/9` rescales Kyte–Doolittle hydrophobicity to [1,2]. PPC/PNC use analogous sums of positive/negative charges. K/R=+1, D/E=−1, H=+0.1, others=0. SFvCSP is the **product** of the total VH and VL surface charges.

A salt bridge connects surface K NZ or R NH1/NH2 with D OD1/OD2 or E OE1/OE2 at ≤3.2 Å. Both residues receive glycine hydrophobicity. Charges are neutralized in `paper` and retained in `web-compatible`.

## Separate profiles

The published TAP2 specification states 4.5 Å, implemented by `paper`. Reproducing all ten archived service results required 4.0 Å and retaining salt-bridge charges in PPC/PNC/SFvCSP. This is an observation on a fixed benchmark, not a universal reconstruction of the server's unavailable source.

Paper thresholds are explicit in `tap2_profiles.PAPER_THRESHOLDS`. Actual server flags are read from HTML rather than replaced with paper thresholds.

## Reproducibility

Weight, PSA, and code SHA256 hashes and package versions are recorded. Seeds are fixed, but refinement and platform differences can alter coordinates. Strict metric regression uses immutable server PDBs; independent folding is assessed separately.
