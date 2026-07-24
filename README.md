# Protein Fit Optimization (ICP-SAXS)

Rigid-body alignment and multi-copy packing of atomic structures into SAXS
envelopes, using ICP with a bidirectional volumetric penalty and Chamfer
Distance in Angstroms as the primary score.

Unlike NSD-based superposition (dimensionless), envelope occupancy here is an
explicit term in the cost function, and the score is reported in physical
units comparable to RMSD.

---

## Benchmark results

50 SASBDB entries with both a dummy-atom envelope and an atomic model
deposited. Identical preprocessing and random seed across all methods.

| Configuration | Chamfer (Å) | Fraction outside | Coverage | p vs. baseline |
|---|---|---|---|---|
| Unconstrained ICP (baseline) | 5.912 | 0.169 | 0.493 | — |
| + bidirectional penalty (λ = 2) | 6.114 | 0.138 | 0.385 | 9.5e-02 |
| + full Cα representation | 5.508 | 0.145 | 0.502 | 8.4e-08 |
| + enantiomorph search (final) | **5.258** | **0.133** | **0.520** | **1.2e-11** |

Comparison against ATSAS `cifsup` 3.2.1 on the same entries:

| Method | Chamfer (Å) | Fraction outside | Coverage | Centroid sep. (Å) |
|---|---|---|---|---|
| cifsup (NSD) | 5.436 | 0.151 | 0.554 | 5.244 |
| cifsup (ICP) | 5.693 | 0.158 | 0.540 | 5.685 |
| ICP-SAXS | 5.258 | 0.133 | 0.520 | 4.970 |

Paired Wilcoxon against `cifsup --method=ICP`: p = 0.80 (Chamfer),
0.63 (fraction outside), 0.43 (centroid separation). **No difference is
statistically significant** — performance is equivalent to ATSAS, not
superior.

---

## Installation

```bash
conda create -n icp python=3.11 numpy scipy pandas matplotlib -y
conda activate icp
```

Or with venv:

```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt matplotlib
```

Dependencies: numpy, scipy, pandas (matplotlib only for figures).
The ATSAS comparison requires `cifsup` on the PATH (ATSAS 3.2+, free
academic licence at biosaxs.com).

---

## Usage

### Benchmark alignment

```bash
python3 scripts/icp_saxs_v2.py \
  --base <folder containing data/sasbdb/> \
  --mode author --penalty 2 --max-points 3000 \
  --out results.csv
```

Relevant parameters:

| Flag | Default | Description |
|---|---|---|
| `--mode` | `bidir` | `author` (bidirectional, plain means — recommended), `bidir` (thresholded), `uni` (leak only) |
| `--penalty` | 50 | Weight λ. **Use 2** — see "Choosing λ" below |
| `--max-points` | 300 | Model downsampling. Use 3000 or more to disable |
| `--no-enantiomorphs` | — | Disables mirror-image testing |
| `--restarts` | 3 | Random initial rotations |
| `--seed` | 42 | Seed, for reproducibility |

### Comparison against ATSAS

```bash
python3 scripts/benchmark_cifsup.py \
  --base <folder containing data/sasbdb/> \
  --ours results.csv --max-points 3000 \
  --out benchmark_cifsup.csv
```

### Paper figures and tables

```bash
python3 scripts/gerar_figuras_artigo.py    # ablation, comparison, enantiomorphs
python3 scripts/gerar_figura_packing.py    # packing with random control
```

### Packing mode (original code)

```bash
python3 src/main.py packing -i data/pdbs/ -e data/envelope.cif -o output/
```

---

## Method notes

### Choosing λ

In the cost function, both terms are means of distances on the same scale:

```
J = mean d(envelope → protein)  +  λ · mean d(protein → envelope)
        (empty space)                    (leaking)
```

A large λ does **not** express a strong preference against leaking — it
suppresses the occupancy term. With that term suppressed, the cost reverts
to leak-only behaviour and the constraint loses its effect. Sweeping λ from
1 to 200 shows monotonic degradation; λ = 2 is the best balance.

### The penalty must enter the Procrustes step

An unweighted Procrustes update produces exactly the same movement with or
without a penalty. A penalty applied only to the reported score leaves the
alignment unchanged while appearing to be active. The weights must enter
the rotation and translation computation itself.

### The occupancy term is what anchors a single copy

A leak-only penalty is satisfied by any pose adjacent to the envelope,
including poses entirely outside it. The occupancy term supplies the
missing constraint.

### Enantiomorphs

SAXS envelopes do not define chirality: ab initio reconstruction determines
shape only up to reflection. Testing a single hand is wrong roughly half the
time. In this benchmark the mirrored hand won in **26 of 50** entries —
close to the expected 50%.

### File parsing

SASBDB envelopes carry a `.cif` extension but not always mmCIF content; some
are PDB-formatted. The scripts detect format from **content**, not extension.
A parser that decides by extension returns empty coordinates silently — this
was the cause of the two "failures" originally reported in the benchmark
(SASDXB2 and SASDXS4), which in fact align normally.

---

## Layout

```
├─ src/                    original pipeline (single, packing)
├─ scripts/
│  ├─ icp_saxs_v2.py               alignment, 3 modes + enantiomorphs
│  ├─ benchmark_cifsup.py          comparison against ATSAS cifsup
│  ├─ gerar_figuras_artigo.py      figures 1-4 and tables
│  ├─ gerar_figura_packing.py      packing random control
│  ├─ extrai_reports.py            batch extraction of report.txt scores
│  └─ icp_lib_ORIGINAL_bidirecional.py   preserved historical version
├─ benchmark/              CSVs for every experiment
├─ figures_paper/          manuscript figures and tables
└─ data/                   not versioned (envelopes and PDBs)
```

---

## Reproducibility

Every number in the tables above comes from the CSVs in `benchmark/` and can
be regenerated with the commands on this page. Results were verified on two
independent installations.

`gerar_figura_packing.py` checks the recomputed Chamfer against the values in
the original `report.txt` files and warns on divergence — SASBDB entries have
several envelopes (damaver, damfilt, damstart, pddf) and only one matches the
coordinates of a given alignment.

---

## Limitations

- Fraction outside uses the convex hull of the envelope, so it underestimates
  violation for concave shapes. Useful for comparing methods on the same
  envelope, not as an absolute measure.
- ICP is a local optimiser; restarts and enantiomorph search mitigate but do
  not guarantee the global optimum.
- Packing shows the envelope accommodates more than one rigid conformer,
  which is consistent with both oligomerisation and a conformational
  ensemble. Distinguishing them requires molecular weight from I(0).
