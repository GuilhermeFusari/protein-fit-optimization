# Protein Fit Optimization (ICP-SAXS)

Rigid-body alignment and multi-copy packing of atomic structures into SAXS
envelopes, using ICP with a bidirectional volumetric penalty and Chamfer
Distance in Angstroms as the primary score.

Unlike NSD-based superposition, which is dimensionless, envelope occupancy
here is an explicit term in the cost function, and the score is reported in
physical units comparable to RMSD.

## Benchmark results

50 SASBDB entries with both a dummy-atom envelope and an atomic model
deposited. Identical preprocessing and random seed across all methods.
Penalty weight λ = 0.2 (calibrated with full Cα representation and
leave-one-out cross-validation; see "Choosing λ").

Effect of each component, relative to unconstrained ICP:

| Configuration | Chamfer (Å) | Fraction outside | Coverage | Centroid sep. (Å) |
|---|---|---|---|---|
| Unconstrained ICP (baseline) | 5.912 | 0.169 | 0.493 | 7.853 |
| + bidirectional penalty (λ = 0.2) | 5.435 | 0.131 | 0.509 | 4.545 |
| + enantiomorph search (final) | **5.170** | **0.124** | **0.528** | 4.733 |

All improvements over baseline are significant (paired Wilcoxon, p < 1e-6).

Comparison against ATSAS `cifsup` 3.2.1 on the same entries, all metrics
recomputed from the aligned coordinates with the same code:

| Method | Chamfer (Å) | Fraction outside | Coverage | Centroid sep. (Å) |
|---|---|---|---|---|
| cifsup (NSD) | 4.997 | 0.123 | 0.538 | 4.197 |
| cifsup (ICP) | 5.331 | 0.136 | 0.518 | 4.853 |
| ICP-SAXS | 5.170 | 0.124 | 0.528 | 4.733 |

Paired Wilcoxon, ICP-SAXS vs. cifsup NSD: cifsup NSD is better in Chamfer
(p < 0.001) and coverage (p = 0.001); the two are statistically
indistinguishable in fraction outside (p = 0.08) and centroid separation
(p = 0.13). ICP-SAXS outperforms `cifsup --method=ICP` (the unconstrained
mode) in all four metrics. **The pipeline does not surpass the NSD mode of
ATSAS; it is a comparable open-source alternative that improves substantially
over unconstrained ICP and reports its metric in Angstroms.**

Both methods also cover both metrics: cifsup optimises NSD, this pipeline
optimises Chamfer. The two are reported side by side in
`benchmark/tabela_chamfer_nsd.csv`.

## Installation

Requires Python 3.9 or later. From a terminal (Linux, macOS or Windows):

```bash
pip install git+https://github.com/GuilhermeFusari/protein-fit-optimization.git
```

This installs the `saxs-icp` command and its dependencies (numpy, scipy,
biopython). Using a virtual environment is recommended:

```bash
python3 -m venv saxs-env
source saxs-env/bin/activate        # Windows: saxs-env\Scripts\activate
pip install git+https://github.com/GuilhermeFusari/protein-fit-optimization.git
```

or with conda:

```bash
conda create -n saxs-icp python=3.11 -y && conda activate saxs-icp
pip install git+https://github.com/GuilhermeFusari/protein-fit-optimization.git
```

Check the installation with `saxs-icp --version`. If the command is not found,
`python -m saxs_icp` works the same way.

## Quick start

Fit an atomic model into a SAXS envelope:

```bash
saxs-icp align model.pdb envelope.cif -o result/
```

* **model**: atomic structure, PDB or mmCIF (e.g. from the PDB or AlphaFold).
* **envelope**: dummy-atom model from DAMMIF, DAMMIN, GASBOR, DAMAVER, etc.,
  PDB or mmCIF. The format is detected from the content, not the extension.

Output in `result/`:

| File | Content |
|---|---|
| `model_aligned.pdb` (or `.cif`) | the full model, all atoms, moved into the envelope frame |
| `envelope_mirrored.pdb` | only when the mirror image fits better (see below) |
| `report.json` | metrics, parameters and the 4×4 transformation matrix |

The defaults are the published method (λ = 0.2, `author` mode, 3 restarts,
enantiomorph search), so no options are needed.

**Mirrored envelopes.** SAXS does not determine chirality, so the envelope may
have the wrong hand. Both hands are tested. When the mirror image fits better,
the protein is *not* mirrored (that would turn it into a D-amino-acid
structure); instead the envelope is reflected and saved as
`envelope_mirrored.pdb`. Open the aligned model together with that file.

Try it with the bundled example:

```bash
git clone https://github.com/GuilhermeFusari/protein-fit-optimization.git
cd protein-fit-optimization
saxs-icp align examples/SASDTS8_model.pdb examples/SASDTS8_envelope.cif -o test_out
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--penalty` | 0.2 | Weight λ of the leak term. See "Choosing λ". 0 = plain ICP |
| `--mode` | `author` | `author` (fill + λ·leak, published), `bidir` (thresholded), `uni` (leak only) |
| `--restarts` | 3 | Initial orientations per hand |
| `--max-iter` | 50 | ICP iterations |
| `--max-points` | 3000 | C-alpha atoms used; larger models are subsampled |
| `--sample-env` | 5000 | Envelope points used (0 = all) |
| `--seed` | 42 | Random seed; results are fully reproducible |
| `--no-enantiomorphs` | off | Do not test the mirror image |

`saxs-icp align --help` lists them all.

### Multi-copy packing

Ranks every `.pdb` in a folder (e.g. frames of a simulation) against one
envelope and writes the best copies aligned:

```bash
saxs-icp pack frames/ envelope.cif -o packing_out/ --copies 20
```

Writes `RANK_01_<name>.pdb` … , `ALL_TOP_ALIGNED.pdb` (all copies as models)
and `report.txt`. Add `--seed N` for reproducible runs.

### Batch benchmark

Runs over a SASBDB-style folder (`<base>/data/sasbdb/<ACC>_envelope.cif` and
`<ACC>_model.pdb`) and writes one CSV row per entry:

```bash
saxs-icp benchmark --base <folder containing data/sasbdb/> --out results.csv
```

With the default options this reproduces `benchmark/final_lam02.csv` exactly.

### Use from Python

```python
import saxs_icp

model = saxs_icp.read_coords("model.pdb", ca_only=True)
env = saxs_icp.read_coords("envelope.cif")
r = saxs_icp.align(model, env)
print(r.metrics["chamfer"], r.mirrored)
T_model, M_env = saxs_icp.proper_transform(r.transform)
```

## Reproducing the paper

The scripts in `scripts/` regenerate every table and figure. They need a few
extra packages:

```bash
pip install -r requirements.txt
```

```bash
# alignment (same result as `saxs-icp benchmark`)
python3 scripts/icp_saxs_v2.py \
  --base <folder containing data/sasbdb/> \
  --mode author --penalty 0.2 --max-points 3000 \
  --out results.csv

# comparison against ATSAS cifsup (needs cifsup on the PATH, ATSAS 3.2+)
python3 scripts/benchmark_cifsup.py \
  --base <folder containing data/sasbdb/> \
  --ours results.csv --max-points 3000 \
  --out benchmark_cifsup.csv

# NSD cross-metric table and lambda cross-validation
python3 scripts/calcular_nsd.py \
  --base <folder> --ours poses_ours --cifsup poses_cifsup \
  --out tabela_chamfer_nsd.csv
python3 scripts/validacao_lambda.py \
  --sweep "sweep_lambda/lam_*.csv" --criterio chamfer --out loo_lambda.csv

# figures
python3 scripts/gerar_figuras_v2.py       # lambda sweep, factorial ablation, cifsup comparison
python3 scripts/gerar_figura_packing.py   # packing with random control
```

Note that `scripts/icp_saxs_v2.py` defaults to `--mode bidir --max-points 300`;
pass the options above to reproduce the published configuration. The
`saxs-icp` command uses the published configuration by default.

## Method notes

### Choosing λ

In the cost function, both terms are means of distances on the same scale:

```
J = mean d(envelope -> protein)  +  λ · mean d(protein -> envelope)
        (empty space)                    (leaking)
```

A large λ does **not** express a strong preference against leaking. It
suppresses the occupancy term, and the cost reverts to leak-only behaviour.
Sweeping λ from 0.1 to 100 with the full Cα representation shows monotonic
degradation for λ above the optimum. The best value is **λ = 0.2**, selected
by leave-one-out cross-validation (the naive optimum and the cross-validated
optimum agree to within 0.1%, so there is no test-set bias). An earlier
calibration reported λ = 2, but it was done with 300-point downsampling — the
degenerate regime — and never tested λ < 1.

### The penalty must enter the Procrustes step

An unweighted Procrustes update produces exactly the same movement with or
without a penalty. A penalty applied only to the reported score leaves the
alignment unchanged while appearing to be active. The weights must enter the
rotation and translation computation itself.

### The occupancy term is what anchors a single copy

A leak-only penalty is satisfied by any pose adjacent to the envelope,
including poses entirely outside it. The occupancy term supplies the missing
constraint. A factorial ablation (penalty × representation) shows the two
components act independently and additively, with negligible interaction:
the penalty improves the result by ~0.4 Å regardless of representation.

### Enantiomorphs

SAXS envelopes do not define chirality. Ab initio reconstruction determines
shape only up to reflection, so testing a single hand is wrong roughly half
the time. In this benchmark the mirrored hand won in **30 of 50** entries,
close to the expected 50 percent.

### File parsing

SASBDB envelopes carry a `.cif` extension but not always mmCIF content, since
some are PDB-formatted. The scripts detect format from **content**, not
extension. A parser that decides by extension returns empty coordinates
silently — this caused two entries (SASDXB2, SASDXS4) to be wrongly reported
as failures; both align normally. The same class of bug affected the ATSAS
comparison: passing an mmCIF envelope to `cifsup` with a `.pdb` extension made
it abort on 27 of 50 entries. `benchmark_cifsup.py` now converts coordinates
before calling cifsup.

## Layout

```
├─ src/saxs_icp/           installable package (the `saxs-icp` command)
│  ├─ core.py                      alignment, 3 modes plus enantiomorphs
│  ├─ io.py                        PDB / mmCIF reading and writing
│  ├─ packing.py                   multi-copy packing
│  ├─ benchmark.py                 batch run over a SASBDB folder
│  └─ cli.py                       command-line interface
├─ tests/                  automated tests (pytest)
├─ examples/               two SASBDB entries for a quick test
├─ scripts/                scripts used for the paper
│  ├─ icp_saxs_v2.py               alignment, 3 modes plus enantiomorphs
│  ├─ benchmark_cifsup.py          comparison against ATSAS cifsup
│  ├─ calcular_nsd.py              Chamfer × NSD cross-metric table
│  ├─ validacao_lambda.py          leave-one-out λ cross-validation
│  ├─ gerar_figuras_v2.py          revised figures (λ sweep, factorial, cifsup)
│  ├─ gerar_figura_packing.py      packing random control
│  ├─ extrai_reports.py            batch extraction of report.txt scores
│  └─ icp_lib_ORIGINAL_bidirecional.py   preserved historical version
├─ benchmark/              CSVs for every experiment (sweep_lambda/, fatorial/)
│  └─ _arquivo_pre_revisao/  superseded λ=2 / pre-bugfix results
├─ figuras_artigo/         current figures
└─ data/                   not versioned (envelopes and PDBs)
```

## Reproducibility

Every number in the tables above comes from the CSVs in `benchmark/` and can
be regenerated with the commands on this page. Results were verified on two
independent installations. The installable package gives bit-identical
results to `scripts/icp_saxs_v2.py` (all 50 entries, every metric), and the
test suite (`pip install ".[test]" && pytest`) checks two entries against the
published values on every change.

## Limitations

* Fraction outside uses the convex hull of the envelope, so it underestimates
  violation for concave shapes. It is useful for comparing methods on the same
  envelope, but not as an absolute measure.
* ICP is a local optimiser. Restarts and enantiomorph search mitigate this but
  do not guarantee the global optimum.
* The multi-copy packing is a representation of conformational occupancy, not
  a physical assembly: for GRB2 the envelope holds ~1.1 protein volumes, so 30
  copies exceed the available volume ~27-fold. It shows which conformations a
  flexible molecule samples, consistent with a conformational ensemble rather
  than oligomerisation.
* Direct validation of the packing against the experimental scattering curve
  (χ² via CRYSOL on the ensemble) is not yet included and is left as future
  work.
* The NSD is computed by the published formula (Kozin & Svergun, 2001), which
  may differ slightly from the ATSAS internal implementation; the Chamfer
  comparison, in this pipeline's own metric, is unambiguous.

## Licence

MIT. See `LICENSE`.
