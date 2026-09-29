#!/usr/bin/env python3
"""
teste_busca.py  -  recuperacao de pose com pose verdadeira EXATA, em funcao
do numero de inicios (restarts). Gera benchmark/teste_busca.csv, citado em
KNOWN_ISSUES.md.

Protocolo:
  entradas   os 50 modelos do benchmark (<base>/data/sasbdb/<ACC>_model.pdb),
             representacao C-alpha
  envelope   construido a partir do proprio modelo na pose de referencia:
             esferas numa grade de 4 A dentro de 5 A de algum C-alpha. A pose
             verdadeira e, por construcao, a identidade (sem DAMMIF, sem PCA)
  sementes   6 por entrada (0..5): cada uma aplica ao modelo uma rotacao
             aleatoria uniforme + translacao N(0, 15 A) (gerador
             default_rng(1000 + semente))
  ajuste     saxs_icp.align com os parametros publicados (modo author,
             lambda 0.2, 50 iteracoes, 3000 pontos, busca de enantiomeros,
             seed = semente), variando apenas restarts = 3, 10, 30
  sucesso    RMSD pareado dos C-alpha < 5 A em relacao a pose verdadeira
  custo      fill + 0.2 * leak da pose achada e da pose verdadeira, ambos no
             envelope completo e em todos os C-alpha (custo_achado >
             custo_verdadeiro numa falha = falha da busca, nao da funcao
             objetivo)

Resultado (benchmark/teste_busca.csv): acerto de 36% / 73% / 95% com
3 / 10 / 30 inicios; 268 das 289 falhas com custo maior que o da pose
verdadeira.

Uso (requer o pacote instalado: pip install .):
    python3 scripts/teste_busca.py <pasta com data/sasbdb/> teste_busca.csv

Leva ~15 min em 8 nucleos.
"""
import csv
import glob
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from saxs_icp import core
from saxs_icp.io import read_coords

BASE = sys.argv[1]
OUT = sys.argv[2]
RESTARTS = (3, 10, 30)
SEEDS = range(6)
GRID, RADIUS, SUCCESS = 4.0, 5.0, 5.0


def envelope_from(ca):
    lo = np.floor((ca.min(0) - RADIUS) / GRID) * GRID
    hi = ca.max(0) + RADIUS
    axes = [np.arange(lo[i], hi[i] + GRID, GRID) for i in range(3)]
    g = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    d, _ = cKDTree(ca).query(g)
    return g[d <= RADIUS]


def cost(pts, env):
    fill, leak = core.author_cost(pts, env, cKDTree(env))
    return fill + 0.2 * leak


def job(args):
    acc, seed = args
    ca = read_coords(os.path.join(BASE, "data", "sasbdb", f"{acc}_model.pdb"), ca_only=True)
    env = envelope_from(ca)
    rng = np.random.default_rng(1000 + seed)
    R = Rotation.from_quat(rng.normal(size=4)).as_matrix()
    t = rng.normal(0, 15, 3)
    c = ca.mean(0)
    mob = (R @ (ca - c).T).T + c + t
    c_true = cost(ca, env)
    rows = []
    for k in RESTARTS:
        r = core.align(mob, env, restarts=k, seed=seed)
        rec = core.apply_T(mob, r.transform)
        rmsd = float(np.sqrt(((rec - ca) ** 2).sum(1).mean()))
        rows.append(dict(entry=acc, seed=seed, restarts=k, rmsd=round(rmsd, 3),
                         sucesso=int(rmsd < SUCCESS), espelhado=int(r.mirrored),
                         custo_achado=round(cost(rec, env), 5),
                         custo_verdadeiro=round(c_true, 5),
                         n_ca=len(ca), n_env=len(env)))
    return rows


if __name__ == "__main__":
    accs = sorted(os.path.basename(p)[:-10] for p in
                  glob.glob(os.path.join(BASE, "data", "sasbdb", "*_model.pdb")))
    jobs = [(a, s) for a in accs for s in SEEDS]
    rows = []
    with Pool() as pool:
        for i, rr in enumerate(pool.imap_unordered(job, jobs), 1):
            rows += rr
            if i % 30 == 0:
                print(f"{i}/{len(jobs)}", flush=True)
    rows.sort(key=lambda r: (r["entry"], r["seed"], r["restarts"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(accs)} entradas x {len(SEEDS)} sementes -> {OUT}")
