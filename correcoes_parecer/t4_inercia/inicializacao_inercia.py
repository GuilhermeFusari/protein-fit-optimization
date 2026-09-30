#!/usr/bin/env python3
"""
inicializacao_inercia.py  -  inicializacao por eixos de inercia (PCA) no
teste sintetico de recuperacao de pose. NAO altera o pacote.

Protocolo identico a scripts/teste_busca.py (mesmas 50 entradas, mesmo
envelope com pose exata, mesmas 6 sementes de embaralhamento, mesmo
criterio RMSD < 5 A, mesmos parametros publicados). Muda so a escolha das
orientacoes iniciais de cada mao (original e espelhada):

  atual  (metodo publicado)  1 = pose de entrada + (k - 1) rotacoes aleatorias
  pca                        4 = eixos principais do modelo alinhados aos do
                             envelope (as 4 combinacoes de sinal com det +1,
                             centroide no centroide do envelope)
                             + (k - 4) rotacoes aleatorias

Orcamentos k (inicios por mao): atual 3, 4, 10, 30; pca 4, 10, 30.
atual 3 = default publicado; atual 4 = mesmo orcamento que pca 4.
atual 3/10/30 tem de reproduzir benchmark/teste_busca.csv (checagem).

Uso:
    python3 inicializacao_inercia.py <pasta com data/sasbdb/> <repo> saida.csv
"""
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from saxs_icp import core
from saxs_icp.io import apply_T, read_coords

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, os.pardir, "scripts"))
import teste_busca as tb  # noqa: E402  (envelope_from, cost, mesmo protocolo)

BUDGETS = [("atual", 3), ("atual", 4), ("atual", 10), ("atual", 30),
           ("pca", 4), ("pca", 10), ("pca", 30)]
SIGNS = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]


def principal_axes(pts):
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
    R = vt.copy()
    if np.linalg.det(R) < 0:
        R[-1] *= -1
    return c, R


def align_pca(model_ca, env_full, restarts, seed, mode="author", penalty=0.2,
              max_iter=50, max_points=3000, sample_env=5000):
    """core.align com inicios por eixos de inercia; resto identico."""
    mod = np.asarray(model_ca, float)
    rng = np.random.default_rng(seed)
    if len(mod) > max_points:
        mod = mod[np.sort(rng.choice(len(mod), max_points, replace=False))]
    env = env_full
    if sample_env and len(env) > sample_env:
        env = env[np.sort(rng.choice(len(env), sample_env, replace=False))]
    thr = max(3.0, core.grid_spacing(env))
    env_tree = cKDTree(env)
    center = mod.mean(axis=0)
    ce, Re = principal_axes(env)

    best, best_mirror = None, False
    for M in (np.eye(3), core.MIRROR_Z.copy()):
        mirrored = not np.allclose(M, np.eye(3))
        if mirrored:
            base = (M @ (mod - center).T).T + center
            T_mir = np.eye(4)
            T_mir[:3, :3], T_mir[:3, 3] = M, center - M @ center
        else:
            base, T_mir = mod, np.eye(4)
        cb, Rb = principal_axes(base)

        starts = []
        for sg in SIGNS[:restarts]:
            Rc = Re.T @ np.diag(sg) @ Rb
            T = np.eye(4)
            T[:3, :3], T[:3, 3] = Rc, ce - Rc @ cb
            starts.append(T)
        for _ in range(restarts - len(starts)):
            Rm = core.random_rotation(int(rng.integers(1 << 31)))
            T = np.eye(4)
            T[:3, :3], T[:3, 3] = Rm, center - Rm @ center
            starts.append(T)

        for T_start in starts:
            T_init = T_start @ T_mir
            start = apply_T(mod, T_init)
            T_icp, cost = core.icp(start, env, env_tree, max_iter, penalty, thr,
                                   mode=mode)
            if best is None or cost < best[0]:
                best = (cost, T_icp @ T_init)
                best_mirror = mirrored
    return best[1], best_mirror


def job(args):
    base, acc, seed = args
    ca = read_coords(os.path.join(base, "data", "sasbdb", f"{acc}_model.pdb"), ca_only=True)
    env = tb.envelope_from(ca)
    rng = np.random.default_rng(1000 + seed)
    R = Rotation.from_quat(rng.normal(size=4)).as_matrix()
    t = rng.normal(0, 15, 3)
    c = ca.mean(0)
    mob = (R @ (ca - c).T).T + c + t
    c_true = tb.cost(ca, env)
    rows = []
    for variant, k in BUDGETS:
        if variant == "atual":
            r = core.align(mob, env, restarts=k, seed=seed)
            T, mirrored = r.transform, r.mirrored
        else:
            T, mirrored = align_pca(mob, env, restarts=k, seed=seed)
        rec = apply_T(mob, T)
        rmsd = float(np.sqrt(((rec - ca) ** 2).sum(1).mean()))
        rows.append(dict(entry=acc, seed=seed, inicializacao=variant, restarts=k,
                         rmsd=round(rmsd, 3), sucesso=int(rmsd < tb.SUCCESS),
                         espelhado=int(mirrored),
                         custo_achado=round(tb.cost(rec, env), 5),
                         custo_verdadeiro=round(c_true, 5)))
    return rows


def main():
    import glob
    import pandas as pd
    base, repo, out = sys.argv[1:4]
    accs = sorted(os.path.basename(p)[:-10] for p in
                  glob.glob(os.path.join(base, "data", "sasbdb", "*_model.pdb")))
    jobs = [(base, a, s) for a in accs for s in tb.SEEDS]
    rows = []
    with Pool() as pool:
        for i, rr in enumerate(pool.imap_unordered(job, jobs), 1):
            rows += rr
            if i % 30 == 0:
                print(f"{i}/{len(jobs)}", flush=True)
    df = pd.DataFrame(rows).sort_values(["entry", "seed", "inicializacao", "restarts"])
    df.to_csv(out, index=False)

    ref = pd.read_csv(os.path.join(repo, "benchmark", "teste_busca.csv"))
    a = df[df.inicializacao == "atual"].merge(ref, on=["entry", "seed", "restarts"],
                                              suffixes=("", "_ref"))
    nd = int((a.rmsd != a.rmsd_ref).sum() + (a.sucesso != a.sucesso_ref).sum())
    print(f"checagem atual 3/10/30 vs benchmark/teste_busca.csv: {len(a)} linhas, "
          f"{nd} diferencas")
    if nd or len(a) != len(ref):
        sys.exit("ERRO: nao reproduz teste_busca.csv")


if __name__ == "__main__":
    main()
