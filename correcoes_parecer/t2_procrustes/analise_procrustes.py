#!/usr/bin/env python3
"""
analise_procrustes.py  -  tabelas do relatorio da tarefa 2 a partir de
pesos_procrustes.csv e benchmark/benchmark_cifsup_lam02.csv (so leitura).

Uso:
    python3 analise_procrustes.py pesos_procrustes.csv <repo>
"""
import os
import sys

import pandas as pd
from scipy import stats

d = pd.read_csv(sys.argv[1])
cif = pd.read_csv(os.path.join(sys.argv[2], "benchmark",
                               "benchmark_cifsup_lam02.csv")).set_index("entry")


def fmt(x):
    return f"{x:.2g}" if x < 0.001 else f"{x:.3f}"


for R in sorted(d.restarts.unique()):
    print(f"\n## {R} inicios\n")
    print("| Variante | Chamfer (A) | custo J (Eq. 1) | Fracao fora | Cobertura "
          "| Sep. centroide | Chamfer: p vs atual | J: p vs atual "
          "| entradas com J menor que atual | Chamfer: p vs cifsup NSD |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    a = d[(d.restarts == R) & (d.variante == "atual")].set_index("entry").loc[cif.index]
    for v in ("atual", "eq1_ls", "eq1_irls"):
        g = d[(d.restarts == R) & (d.variante == v)].set_index("entry").loc[cif.index]
        if v == "atual":
            pa = pj = "—"
            nj = "—"
        else:
            pa = fmt(stats.wilcoxon(g.chamfer, a.chamfer).pvalue)
            pj = fmt(stats.wilcoxon(g.custo, a.custo).pvalue)
            nj = f"{int((g.custo < a.custo).sum())}/50"
        pn = fmt(stats.wilcoxon(g.chamfer, cif.nsd_chamfer).pvalue)
        print(f"| {v} | {g.chamfer.mean():.3f} | {g.custo.mean():.4f} | "
              f"{g.frac_fora.mean():.3f} | {g.cobertura.mean():.3f} | "
              f"{g.sep_centroide.mean():.3f} | {pa} | {pj} | {nj} | {pn} |")
print(f"\ncifsup NSD: Chamfer medio {cif.nsd_chamfer.mean():.3f}")
