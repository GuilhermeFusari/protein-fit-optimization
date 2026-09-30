#!/usr/bin/env python3
"""
selecao_vs_aleatorio.py  -  o packing escolhe conformacoes que reproduzem a
curva experimental melhor que o acaso?

Pool: os N membros cujos perfis CRYSOL ja foram calculados por
ensemble_chi2.py (perfis_membros.dat). Para cada k:
  * roda `saxs-icp pack` (codigo atual, lambda 0.2) escolhendo k copias, com
    10 sementes; le do report.txt quais membros foram escolhidos
  * chi^2 (so escala) do perfil medio da selecao do packing
  * chi^2 de TODOS os C(N, k) subconjuntos de k membros (o "acaso")
  * percentil da selecao do packing nessa distribuicao (menor = melhor)

A selecao do packing usa so a geometria (encaixe no envelope); a curva
experimental nao entra nela. O chi^2 e uma validacao independente.

Uso:
    python3 selecao_vs_aleatorio.py <perfis_membros.dat> <pool/> <envelope> <saida/>
"""
import itertools
import json
import os
import re
import subprocess
import sys

import numpy as np

KS = (5, 10)
SEEDS = range(1, 11)


def load(path):
    names = open(path).readline().lstrip("#").split()[3:]
    d = np.loadtxt(path)
    return d[:, 1], d[:, 2], d[:, 3:].T, [re.sub(r"^RANK_\d+_", "", n) for n in names]


def chi2_many(I, s, P):
    """chi^2 reduzido (so escala) para cada linha de P (perfis candidatos)."""
    w = 1.0 / s ** 2
    c = (P * I * w).sum(1) / (P ** 2 * w).sum(1)
    return (((I - c[:, None] * P) / s) ** 2).sum(1) / (len(I) - 1)


def main():
    prof_path, pool, env, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    I, s, prof, names = load(prof_path)
    idx = {n: i for i, n in enumerate(names)}
    N = len(names)
    res = dict(n_pool=N, chi2_pool_completo=float(chi2_many(I, s, prof.mean(0)[None])[0]),
               por_k={})
    for k in KS:
        combos = np.array(list(itertools.combinations(range(N), k)))
        allchi = np.concatenate([chi2_many(I, s, prof[c].mean(1))
                                 for c in np.array_split(combos, max(1, len(combos) // 20000))])
        picks = []
        for seed in SEEDS:
            o = os.path.join(out, f"pack_k{k}_s{seed}")
            subprocess.run(["saxs-icp", "pack", pool, env, "-o", o, "--copies", str(k),
                            "--seed", str(seed), "--workers", "8"],
                           check=True, capture_output=True)
            chosen = re.findall(r"\|\s*(\S+\.pdb)\s*$", open(os.path.join(o, "report.txt")).read(),
                                re.M)
            sel = sorted(idx[c] for c in chosen)
            chi = float(chi2_many(I, s, prof[sel].mean(0)[None])[0])
            pct = float((allchi < chi).mean() * 100)
            picks.append(dict(seed=seed, membros=[names[i] for i in sel],
                              chi2=round(chi, 3), percentil=round(pct, 1)))
        pc = [p["percentil"] for p in picks]
        res["por_k"][str(k)] = dict(
            n_subconjuntos=len(combos),
            chi2_aleatorio_mediana=round(float(np.median(allchi)), 3),
            chi2_aleatorio_p5=round(float(np.percentile(allchi, 5)), 3),
            chi2_aleatorio_p95=round(float(np.percentile(allchi, 95)), 3),
            chi2_melhor_subconjunto=round(float(allchi.min()), 3),
            chi2_packing_mediana=round(float(np.median([p["chi2"] for p in picks])), 3),
            percentil_packing_mediana=round(float(np.median(pc)), 1),
            percentil_packing_min_max=[min(pc), max(pc)],
            selecoes=picks)
        print(f"k={k}: acaso mediana {np.median(allchi):.2f} (p5 {np.percentile(allchi, 5):.2f}, "
              f"p95 {np.percentile(allchi, 95):.2f}, melhor {allchi.min():.2f}) | packing "
              f"mediana {np.median([p['chi2'] for p in picks]):.2f}, percentis {sorted(pc)}",
              flush=True)
    json.dump(res, open(os.path.join(out, "selecao_vs_aleatorio.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
