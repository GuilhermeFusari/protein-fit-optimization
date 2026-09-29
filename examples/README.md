# Examples

Two entries from the 50-entry benchmark, taken from
[SASBDB](https://www.sasbdb.org) (accession codes SASDTS8 and SASDXB2).

| Entry | Envelope file | Model file | Best hand |
|---|---|---|---|
| SASDTS8 | PDB content (`.cif` name) | PDB | original |
| SASDXB2 | mmCIF (GASBOR) | mmCIF (AlphaFold, `.pdb` name) | mirrored |

```bash
saxs-icp align SASDTS8_model.pdb SASDTS8_envelope.cif -o out_ts8
saxs-icp align SASDXB2_model.pdb SASDXB2_envelope.cif -o out_xb2
```

With `--seed 44` (SASDTS8) and `--seed 83` (SASDXB2), the seeds these entries
receive in the batch benchmark, the metrics match `benchmark/final_lam02.csv`
exactly (Chamfer 3.0231 Å and 4.1722 Å).
