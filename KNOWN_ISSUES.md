# Known issues and limitations

What to expect from `saxs-icp` 1.0.0. The package reproduces the published
benchmark numbers (`benchmark/final_lam02.csv`,
`benchmark/benchmark_cifsup_lam02.csv`) exactly; the items below are
limitations of the method and of the analyses behind those numbers, not
differences between the package and the paper.

## Main limitation: the search

* **The number of initial orientations is the real bottleneck.** ICP is a
  local optimiser, and the default of 3 starts per hand is often not enough
  to find the correct pose. In a synthetic test with an exactly known true
  pose, the pose was recovered (C-alpha RMSD < 5 Å) in 36% of cases with 3
  starts, 73% with 10 and 95% with 30.
  Test: each of the 50 benchmark models, envelope built from the model
  itself (beads on a 4 Å grid within 5 Å of a C-alpha), 6 random
  rotations + translations per model, published parameters otherwise
  (`scripts/teste_busca.py`, results in `benchmark/teste_busca.csv`).
  Success rates varied strongly with the random scramble (8–80% across
  the 6 seeds with 3 starts).
* **Most failures come from the search, not from the objective function.**
  In 268 of the 289 failed runs (93%), the pose found had a higher cost than
  the true pose: the cost function ranked the correct answer better, but the
  search did not reach it. In the remaining 21 (5 of the 50 models), the
  pose found had an equal or lower cost: there, an alternative orientation
  fits the envelope as well as the true one, and more starts cannot fix it.
* **Until this is improved,** use more starts (`--restarts 30`); run time
  grows roughly in proportion. The published benchmark numbers were
  obtained with 3 starts; increasing the number of starts lowers the mean
  Chamfer distance on the same 50 entries. The size of this effect is
  deliberately not reported here until a significance test has been done.
* **Planned:** initialisation along the principal axes of inertia, plus 30
  to 100 starts.

## Method

* **Not better than ATSAS `cifsup` in NSD mode.** On the 50-entry benchmark,
  `cifsup --method=NSD` gives a lower Chamfer distance (4.997 vs 5.170 Å,
  p < 0.001) and higher coverage (p = 0.001). The two are statistically
  indistinguishable in fraction outside and centroid separation. ICP-SAXS is
  an open-source alternative that improves on unconstrained ICP, not a
  replacement that outperforms the NSD mode.
* **Chirality is not determined.** SAXS envelopes do not define handedness.
  When the mirror image fits better (30 of 50 benchmark entries), the
  envelope, not the protein, is reflected and written as
  `envelope_mirrored.pdb`. Downstream analyses must use that file together
  with the aligned model.
* **λ = 0.2 was calibrated on one benchmark.** It was chosen by
  leave-one-out cross-validation on 50 SASBDB entries with the C-alpha
  representation. It has not been validated on other envelope types
  (e.g. very low-resolution or multi-phase models) or with all-atom
  representations.
* **Large models are subsampled.** Models with more than `--max-points`
  (default 3000) C-alpha atoms are randomly subsampled. Results then depend
  on `--seed`.

## Metrics

* **Fraction outside uses the convex hull** of the envelope, so it
  underestimates violation for concave shapes. Use it to compare methods on
  the same envelope, not as an absolute measure.
* **NSD is computed by the published formula** (Kozin & Svergun, 2001) in
  `scripts/calcular_nsd.py`, which may differ slightly from the ATSAS
  internal implementation.
* **No validation against the scattering curve.** The fit is purely
  geometric (model vs. envelope). χ² against the experimental data (e.g.
  CRYSOL) is not computed.

## Multi-copy packing (`saxs-icp pack`)

* **Exploratory mode.** The published packing values were produced by an
  earlier version of the code (before the recalibration to λ = 0.2) and
  have not yet been regenerated with the current pipeline. Running
  `saxs-icp pack` today will not necessarily reproduce them. The mode should
  be treated as exploratory until the ensembles are validated against the
  experimental scattering curve I(q) by χ².
* **A representation of conformational occupancy, not a physical assembly.**
  Copies overlap freely: for GRB2 the envelope holds ~1.1 protein volumes,
  so 30 copies exceed the available volume ~27-fold. The output shows which
  conformations fit the envelope, consistent with a conformational ensemble,
  not an oligomer.
* **Not reproducible without `--seed`.** Random initial orientations are
  drawn per structure; pass `--seed N` for repeatable runs.
* Uses plain (unweighted) Procrustes and no enantiomorph search, unlike
  `saxs-icp align`.

## Pose-recovery test (`scripts/ground_truth.py`)

The synthetic ground-truth test has known defects; a corrected version is
planned. Its results should be read with these in mind:

* **The rotation-error angle includes improper reflections.** The angle is
  computed from the trace of the relative matrix even when that matrix is a
  reflection (det = -1), where the formula does not give a valid rotation
  angle.
* **Seeds are not reproducible.** The per-entry seed uses Python's
  `hash()` of the accession code, which changes between runs unless
  `PYTHONHASHSEED` is fixed.
* **`cifsup` is not run in this test,** so it does not provide a comparison
  with ATSAS on known poses.

## Input files

* **Multi-model files are aligned as one rigid body.** All MODEL records
  (e.g. an NMR ensemble) are read together; a warning is printed. Split the
  file first to align each model separately.
* **Records after `ENDMDL` are read.** Some SASBDB model files contain a
  second copy of the atoms after `ENDMDL` without a new `MODEL` record;
  `saxs-icp` reads all of them (as the published scripts did), while some
  other programs stop at `ENDMDL`.
* **Calcium ions named `CA` count as C-alpha atoms** when selecting the
  C-alpha trace from PDB files. The effect is negligible for proteins but
  is kept for consistency with the published results.
* **ANISOU records are dropped** from the aligned PDB output, since the
  anisotropic factors would no longer match the rotated frame.
* **mmCIF rows spanning several lines** (semicolon text fields inside
  `_atom_site`) are not supported. If the aligned mmCIF cannot be written,
  a C-alpha-only PDB is written instead and a warning is printed.

## Paper scripts (`scripts/`)

* **The file readers in `scripts/` can misread some mmCIF files.** They try
  the PDB fixed-column format first; for some mmCIF files whose `ATOM` lines
  happen to parse as PDB columns, they return wrong coordinates without
  error. A full scan of our local SASBDB copy found 35 such files, none of
  them among the files used for the 50 benchmark entries (every file read
  for the paper was checked against Biopython). The `saxs-icp` package detects the format from
  the content and is not affected. Use the package, not the scripts, for new
  data.
* `scripts/icp_saxs_v2.py` defaults to `--mode bidir --max-points 300`,
  which is not the published configuration; see the README for the exact
  options. `saxs-icp` uses the published configuration by default.

## Platforms

* Tested on Linux with Python 3.10–3.13, NumPy 1.21–2.5 and SciPy 1.7–1.18
  (bit-identical results). The continuous-integration workflow runs the test
  suite on Linux, macOS and Windows with Python 3.9, 3.11 and 3.13.
