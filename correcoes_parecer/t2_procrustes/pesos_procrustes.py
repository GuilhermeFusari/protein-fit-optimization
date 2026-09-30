#!/usr/bin/env python3
"""
pesos_procrustes.py  -  compara os pesos do passo de Procrustes do metodo
publicado com pesos coerentes com a Eq. 1. NAO altera o pacote.

Eq. 1 (README, "Choosing lambda"):
    J = mean_e d(e, P) + lambda * mean_p d(p, E)
        (preenchimento)          (vazamento)

Variantes do passo de Procrustes (tudo o mais identico ao metodo publicado:
custo J para escolher a melhor iteracao, 50 iteracoes, restarts, busca de
enantiomeros, amostragem, sementes):

  atual     codigo publicado (saxs_icp.core.icp, modo author):
              vazamento  w_i = 1 + lambda * d_i / mean(d)
              preench.   w_j = d_j / mean(d_env)
              todos normalizados juntos -> razao entre os blocos
              = N_P (1 + lambda) / N_E, e nao lambda
  eq1_ls    blocos com os pesos da Eq. 1, uniformes dentro de cada bloco
              vazamento  w_i = lambda / N_P
              preench.   w_j = 1 / N_E
            (minimiza a versao quadratica da Eq. 1 com correspondencias fixas)
  eq1_irls  minimos quadrados reponderados para a Eq. 1 propriamente dita
            (media de distancias, nao de quadrados):
              vazamento  w_i = lambda / (N_P * max(d_i, eps))
              preench.   w_j = 1 / (N_E * max(d_j, eps)),  eps = 0.1 A
            com correspondencias fixas, cada passo nao aumenta J
            (majorizacao do tipo Weiszfeld)

Checagem embutida: `atual` com restarts=3 reproduz final_lam02.csv.

Uso:
    python3 pesos_procrustes.py <pasta com data/sasbdb/> <repo> saida.csv
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree

from saxs_icp import core
from saxs_icp.benchmark import find_jobs
from saxs_icp.io import read_coords

LAMBDA = 0.2
EPS = 0.1
VARIANTS = ("atual", "eq1_ls", "eq1_irls")
RESTARTS = (3, 30)


def icp_eq1(src_pts, env_pts, env_tree, max_iter, penalty, threshold,
            mode="author", w_leak=1.0, w_cover=1.0, variant="eq1_ls"):
    """Mesmo laco de core.icp (modo author); muda so os pesos do Procrustes."""
    src = src_pts.copy()
    T_total = np.eye(4)
    best_cost, best_T = float("inf"), np.eye(4)
    nP, nE = len(src), len(env_pts)

    for _ in range(max_iter):
        dist, idx = env_tree.query(src, k=1)
        closest = env_pts[idx]
        fill, leak = core.author_cost(src, env_pts, env_tree)
        cost = fill + penalty * leak
        if cost < best_cost:
            best_cost, best_T = cost, T_total.copy()

        d_env, i_env = cKDTree(src).query(env_pts, k=1)
        if variant == "eq1_ls":
            w = np.full(nP, penalty / nP)
            w_extra = np.full(nE, 1.0 / nE)
        else:
            w = penalty / (nP * np.maximum(dist, EPS))
            w_extra = 1.0 / (nE * np.maximum(d_env, EPS))
        srcA = np.vstack([src, src[i_env]])
        tgtA = np.vstack([closest, env_pts])
        wA = np.concatenate([w, w_extra])

        wA = wA / wA.sum()
        sm = (wA[:, None] * srcA).sum(axis=0)
        tm = (wA[:, None] * tgtA).sum(axis=0)
        H = (srcA - sm).T @ (wA[:, None] * (tgtA - tm))
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        t = tm - R @ sm

        src = (R @ src.T).T + t
        T = np.eye(4)
        T[:3, :3], T[:3, 3] = R, t
        T_total = T @ T_total

    fill, leak = core.author_cost(src, env_pts, env_tree)
    cost = fill + penalty * leak
    if cost < best_cost:
        best_cost, best_T = cost, T_total.copy()
    return best_T, best_cost


def job(args):
    acc, env_path, mod_path, seed, restarts, variant = args
    t0 = time.time()
    env_full = read_coords(env_path)
    mod = read_coords(mod_path, ca_only=True)
    if variant == "atual":
        r = core.align(mod, env_full, restarts=restarts, seed=seed)
    else:
        orig = core.icp
        core.icp = lambda *a, **k: icp_eq1(*a, **k, variant=variant)
        try:
            r = core.align(mod, env_full, restarts=restarts, seed=seed)
        finally:
            core.icp = orig
    m = r.metrics
    return dict(entry=acc, variante=variant, restarts=restarts, seed=seed,
                custo=round(r.cost, 4), chamfer=round(m["chamfer"], 4),
                chamfer_envfull=round(core.chamfer(r.fitted, env_full), 4),
                frac_fora=round(m["frac_outside"], 4),
                cobertura=round(m["coverage"], 4),
                sep_centroide=round(m["centroid_sep"], 3),
                espelhado=bool(r.mirrored), segundos=round(time.time() - t0, 2))


def main():
    base, repo, out = sys.argv[1:4]
    import pandas as pd

    jobs0 = find_jobs(base, seed=42)
    rows = []
    for R in RESTARTS:
        for v in VARIANTS:
            t0 = time.time()
            with Pool() as pool:
                res = pool.map(job, [(j[0], j[1], j[2], j[8], R, v) for j in jobs0])
            rows += res
            print(f"restarts={R:2d} {v:9s}: chamfer medio "
                  f"{np.mean([r['chamfer'] for r in res]):.4f}  custo medio "
                  f"{np.mean([r['custo'] for r in res]):.4f} ({time.time() - t0:.0f}s)",
                  flush=True)
    df = pd.DataFrame(rows).sort_values(["restarts", "variante", "entry"])
    df.to_csv(out, index=False)

    ref = pd.read_csv(os.path.join(repo, "benchmark", "final_lam02.csv")).set_index("entry")
    a = df[(df.variante == "atual") & (df.restarts == 3)].set_index("entry").loc[ref.index]
    c = ["custo", "chamfer", "frac_fora", "cobertura", "sep_centroide", "espelhado"]
    ndiff = int((a[c] != ref[c]).sum().sum())
    print(f"checagem atual/restarts=3 vs final_lam02.csv: {ndiff} diferencas")
    if ndiff:
        sys.exit("ERRO: nao reproduz o CSV publicado")


if __name__ == "__main__":
    main()
