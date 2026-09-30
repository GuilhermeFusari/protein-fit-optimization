#!/usr/bin/env python3
"""
investigar_referencial.py  -  modelo e envelope depositados ja estao no
mesmo sistema de coordenadas? (so leitura: nao move nem alinha nada)

Para cada par (modelo, envelope), na pose DEPOSITADA:
  sep_centroide    distancia entre os centroides (C-alpha x envelope)
  chamfer_dep      Chamfer na pose depositada
  frac_fora_dep    fracao dos C-alpha fora do casco convexo do envelope
  pct_orientacao   % de 1000 rotacoes aleatorias em torno do centroide do
                   modelo (mesma translacao) com Chamfer MENOR que o da pose
                   depositada. Perto de 0% = a orientacao depositada e
                   especial (ja encaixa); ~50% = orientacao arbitraria.
                   Centroides coincidentes sozinhos nao provam nada: envelopes
                   do DAMMIF e muitos modelos vem centrados na origem.
  pct_espelho      o mesmo para a imagem especular da pose depositada (se
                   a pose depositada for de outra mao, o espelho encaixa)

Uso:
    python3 investigar_referencial.py <pasta com data/sasbdb/> saida.csv
"""
import glob
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.spatial.transform import Rotation

from saxs_icp.io import read_coords

N_ROT = 1000


def chamfer(a, b):
    return float((cKDTree(b).query(a)[0].mean() + cKDTree(a).query(b)[0].mean()) / 2)


def frac_out(p, env):
    try:
        return float((Delaunay(env).find_simplex(p) < 0).mean())
    except Exception:
        return float("nan")


def job(args):
    acc, env_path, mod_path = args
    env = read_coords(env_path)
    ca = read_coords(mod_path, ca_only=True)
    if len(env) < 10 or len(ca) < 4:
        return dict(entry=acc, status="coordenadas_vazias")
    if len(env) > 5000:
        env = env[np.random.default_rng(0).choice(len(env), 5000, replace=False)]
    c = ca.mean(0)
    ch = chamfer(ca, env)
    rots = Rotation.random(N_ROT, random_state=1).as_matrix()
    rand = np.array([chamfer((R @ (ca - c).T).T + c, env) for R in rots])
    mir = ca.copy()
    mir[:, 2] = 2 * c[2] - mir[:, 2]
    ch_m = chamfer(mir, env)
    return dict(entry=acc, status="ok", n_ca=len(ca), n_env=len(env),
                sep_centroide=round(float(np.linalg.norm(c - env.mean(0))), 2),
                centroide_modelo=np.round(c, 1).tolist(),
                centroide_envelope=np.round(env.mean(0), 1).tolist(),
                chamfer_dep=round(ch, 3), chamfer_espelho=round(ch_m, 3),
                frac_fora_dep=round(frac_out(ca, env), 3),
                chamfer_aleat_mediana=round(float(np.median(rand)), 3),
                chamfer_aleat_min=round(float(rand.min()), 3),
                pct_orientacao=round(float((rand < ch).mean() * 100), 1),
                pct_espelho=round(float((rand < ch_m).mean() * 100), 1))


def main():
    import pandas as pd
    base, out = sys.argv[1:3]
    sas = os.path.join(base, "data", "sasbdb")
    jobs = []
    for e in sorted(glob.glob(os.path.join(sas, "*_envelope.cif"))):
        acc = os.path.basename(e)[:-13]
        m = os.path.join(sas, f"{acc}_model.pdb")
        if os.path.exists(m):
            jobs.append((acc, e, m))
    with Pool() as pool:
        rows = pool.map(job, jobs)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
