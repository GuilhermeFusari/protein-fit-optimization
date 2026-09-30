#!/usr/bin/env python3
"""
inercia_dados_reais.py  -  complemento da tarefa 4: a inicializacao por
eixos de inercia em envelopes REAIS, onde os eixos do envelope e do modelo
nao coincidem por construcao (no teste sintetico coincidem).

Parte A - benchmark das 50 entradas (envelopes do SASBDB): Chamfer contra o
  envelope completo (como o cifsup em benchmark_cifsup_lam02.csv), com
  inicializacao pca e 4/10/30 inicios por mao, semente publicada (42 + i).
  A comparacao com a inicializacao atual usa chamfer_envfull de
  ../t1_inicios/efeito_inicios.csv (mesmas sementes).

Parte B - envelopes do DAMMIF da tarefa 3 (../t3_ground_truth/envelopes_v2),
  mesmos embaralhamentos e sementes de ground_truth_v2.py, sucesso = RMSD
  < 5 A da "pose verdadeira" por PCA. Compara core.align (atual, 3 e 30
  inicios) com pca (4 e 30 inicios). Atencao: a "pose verdadeira" desse
  teste e definida por PCA, o que favorece por construcao uma
  inicializacao por PCA (ver RELATORIO.md).

Uso:
    python3 inercia_dados_reais.py <pasta com data/sasbdb/> <repo>
"""
import os
import sys
import zlib
from multiprocessing import Pool

import numpy as np
from scipy.spatial.transform import Rotation

from saxs_icp import core
from saxs_icp.benchmark import find_jobs
from saxs_icp.io import apply_T, read_coords

from inicializacao_inercia import align_pca

HERE = os.path.dirname(os.path.abspath(__file__))
GT_ENTRIES = ["SASDPT5", "SASDTS8", "SASDUD6", "SASDUN8",
              "SASDUX5", "SASDV62", "SASDVK8", "SASDWG8"]


def job_bench(args):
    acc, env_path, mod_path, seed, k = args
    env_full = read_coords(env_path)
    mod = read_coords(mod_path, ca_only=True)
    T, mirrored = align_pca(mod, env_full, restarts=k, seed=seed)
    fit = apply_T(mod, T)
    return dict(parte="A_benchmark", entry=acc, inicializacao="pca", restarts=k,
                seed=seed, chamfer_envfull=round(core.chamfer(fit, env_full), 4),
                espelhado=int(mirrored))


def read_ca_pdb(path):
    # mesmo leitor do ground_truth_v2 (C-alpha, colunas fixas)
    return read_coords(path, ca_only=True)


def job_gt(args):
    base, acc, variant, k = args
    ref = read_ca_pdb(os.path.join(base, "data", "sasbdb", f"{acc}_model.pdb"))
    env = read_coords(os.path.join(HERE, os.pardir, "t3_ground_truth",
                                   "envelopes_v2", f"{acc}_env.pdb"))
    rng = np.random.default_rng(42 + zlib.crc32(acc.encode()))
    rows = []
    for tr in range(4):
        Rtrue = Rotation.from_quat(rng.normal(size=4)).as_matrix()
        ttrue = rng.normal(0, 15, 3)
        c = ref.mean(0)
        mob = (Rtrue @ (ref - c).T).T + c + ttrue
        fit_seed = int(rng.integers(1 << 31))
        if variant == "atual":
            r = core.align(mob, env, restarts=k, seed=fit_seed)
            T, mirrored = r.transform, r.mirrored
        else:
            T, mirrored = align_pca(mob, env, restarts=k, seed=fit_seed)
        rec = apply_T(mob, T)
        rmsd = float(np.sqrt(((rec - ref) ** 2).sum(1).mean()))
        rows.append(dict(parte="B_ground_truth", entry=acc, inicializacao=variant,
                         restarts=k, seed=fit_seed, trial=tr, rmsd=round(rmsd, 3),
                         sucesso=int(rmsd < 5.0), espelhado=int(mirrored)))
    return rows


def main():
    import pandas as pd
    base, repo = sys.argv[1:3]
    from scipy import stats

    # parte A
    jobs0 = find_jobs(base, seed=42)
    rows = []
    with Pool() as pool:
        for k in (4, 10, 30):
            rows += pool.map(job_bench, [(j[0], j[1], j[2], j[8], k) for j in jobs0])
            print(f"A: pca {k} feito", flush=True)
        gt = pool.map(job_gt, [(base, a, v, k) for a in GT_ENTRIES
                               for v, k in (("atual", 3), ("atual", 30),
                                            ("pca", 4), ("pca", 30))])
    for g in gt:
        rows += g
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "inercia_dados_reais.csv"), index=False)

    t1 = pd.read_csv(os.path.join(HERE, os.pardir, "t1_inicios", "efeito_inicios.csv"))
    t1 = t1[t1.seed_base == 42]
    cif = pd.read_csv(os.path.join(repo, "benchmark", "benchmark_cifsup_lam02.csv")).set_index("entry")
    A = df[df.parte == "A_benchmark"]
    print("\nParte A: Chamfer (envelope completo), 50 entradas, semente 42")
    print(f"  cifsup NSD: {cif.nsd_chamfer.mean():.3f}")
    for R in (3, 10, 30, 100):
        g = t1[t1.restarts == R].set_index("entry").loc[cif.index]
        print(f"  atual {R:3d}: {g.chamfer_envfull.mean():.3f}  p vs NSD "
              f"{stats.wilcoxon(g.chamfer_envfull, cif.nsd_chamfer).pvalue:.3g}")
    for k in (4, 10, 30):
        g = A[A.restarts == k].set_index("entry").loc[cif.index]
        a3 = t1[t1.restarts == 3].set_index("entry").loc[cif.index]
        print(f"  pca   {k:3d}: {g.chamfer_envfull.mean():.3f}  p vs NSD "
              f"{stats.wilcoxon(g.chamfer_envfull, cif.nsd_chamfer).pvalue:.3g}  "
              f"p vs atual 3 {stats.wilcoxon(g.chamfer_envfull, a3.chamfer_envfull).pvalue:.3g}"
              f"  melhor que atual 3 em {int((g.chamfer_envfull < a3.chamfer_envfull - 1e-4).sum())}/50")
    B = df[df.parte == "B_ground_truth"]
    print("\nParte B: envelopes DAMMIF (8 entradas x 4 tentativas)")
    for (v, k), g in B.groupby(["inicializacao", "restarts"]):
        print(f"  {v:5s} {k:2d}: sucesso {g.sucesso.sum()}/32  RMSD mediano {g.rmsd.median():.1f}"
              f"  espelhado {g.espelhado.sum()}")


if __name__ == "__main__":
    main()
