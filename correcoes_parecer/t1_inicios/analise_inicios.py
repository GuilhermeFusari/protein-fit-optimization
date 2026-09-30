#!/usr/bin/env python3
"""
analise_inicios.py  -  tabelas do relatorio da tarefa 1 a partir de
efeito_inicios.csv e de benchmark/benchmark_cifsup_lam02.csv (so leitura).

Wilcoxon pareado bilateral (scipy.stats.wilcoxon, default), como em
scripts/gerar_figuras_v2.py, sobre as mesmas 50 entradas.

Uso:
    python3 analise_inicios.py efeito_inicios.csv <repo>
"""
import os
import sys

import pandas as pd
from scipy import stats

d = pd.read_csv(sys.argv[1])
cif = pd.read_csv(os.path.join(sys.argv[2], "benchmark", "benchmark_cifsup_lam02.csv"))
cif = cif.set_index("entry")
MET = [("chamfer", "nsd_chamfer", "Chamfer (A)"),
       ("frac_fora", "nsd_frac_fora", "Fracao fora"),
       ("cobertura", "nsd_cobertura", "Cobertura"),
       ("sep_centroide", "nsd_sep_centroide", "Sep. centroide (A)")]


def p(a, b):
    return stats.wilcoxon(a, b).pvalue


def fmt(x):
    return f"{x:.2g}" if x < 0.001 else f"{x:.3f}"


print(f"Envelopes com mais de 5000 pontos (Chamfer amostrado != completo): "
      f"{int((d[(d.seed_base == 42) & (d.restarts == 3)].n_env_full > 5000).sum())}\n")

print("## Semente publicada (42): medias e Wilcoxon contra cifsup NSD\n")
print("| Inicios | Chamfer | p | Chamfer env. completo | p | Fracao fora | p "
      "| Cobertura | p | Sep. centroide | p |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
g42 = d[d.seed_base == 42]
for R, g in g42.groupby("restarts"):
    g = g.set_index("entry").loc[cif.index]
    cells = [str(R)]
    cells += [f"{g.chamfer.mean():.3f}", fmt(p(g.chamfer, cif.nsd_chamfer))]
    cells += [f"{g.chamfer_envfull.mean():.3f}", fmt(p(g.chamfer_envfull, cif.nsd_chamfer))]
    for ours, nsd, _ in MET[1:]:
        cells += [f"{g[ours].mean():.3f}", fmt(p(g[ours], cif[nsd]))]
    print("| " + " | ".join(cells) + " |")
print("| cifsup NSD | " + " | ".join(
    [f"{cif.nsd_chamfer.mean():.3f}", "", f"{cif.nsd_chamfer.mean():.3f}", ""]
    + sum([[f"{cif[n].mean():.3f}", ""] for _, n, _ in MET[1:]], [])) + " |")

print("\n## Direcao do efeito (semente 42): entradas em que ICP-SAXS < cifsup NSD em Chamfer\n")
print("| Inicios | ICP-SAXS melhor | cifsup melhor | diferenca media (A) |")
print("|---|---|---|---|")
for R, g in g42.groupby("restarts"):
    g = g.set_index("entry").loc[cif.index]
    diff = g.chamfer - cif.nsd_chamfer
    print(f"| {R} | {int((diff < 0).sum())} | {int((diff > 0).sum())} | {diff.mean():+.3f} |")

print("\n## Robustez: Chamfer medio e p contra cifsup NSD por semente-base\n")
print("| Inicios | " + " | ".join(f"semente {s}" for s in sorted(d.seed_base.unique())) + " |")
print("|---|" + "---|" * d.seed_base.nunique())
for R in sorted(d.restarts.unique()):
    cells = [str(R)]
    for s in sorted(d.seed_base.unique()):
        g = d[(d.restarts == R) & (d.seed_base == s)].set_index("entry").loc[cif.index]
        cells.append(f"{g.chamfer.mean():.3f} (p={fmt(p(g.chamfer, cif.nsd_chamfer))})")
    print("| " + " | ".join(cells) + " |")

print("\n## Ganho sobre 3 inicios (mesma semente), Wilcoxon pareado\n")
print("| Inicios | " + " | ".join(f"semente {s}" for s in sorted(d.seed_base.unique())) + " |")
print("|---|" + "---|" * d.seed_base.nunique())
for R in sorted(d.restarts.unique()):
    if R == 3:
        continue
    cells = [str(R)]
    for s in sorted(d.seed_base.unique()):
        a = d[(d.restarts == 3) & (d.seed_base == s)].set_index("entry").chamfer
        b = d[(d.restarts == R) & (d.seed_base == s)].set_index("entry").chamfer.loc[a.index]
        cells.append(f"{(b - a).mean():+.3f} A (p={fmt(p(b, a))})")
    print("| " + " | ".join(cells) + " |")

print("\n## Custo (Eq. 1) medio: a busca melhora o proprio objetivo?\n")
print("| Inicios | " + " | ".join(f"semente {s}" for s in sorted(d.seed_base.unique())) + " |")
print("|---|" + "---|" * d.seed_base.nunique())
for R in sorted(d.restarts.unique()):
    cells = [str(R)] + [f"{d[(d.restarts == R) & (d.seed_base == s)].custo.mean():.4f}"
                        for s in sorted(d.seed_base.unique())]
    print("| " + " | ".join(cells) + " |")
