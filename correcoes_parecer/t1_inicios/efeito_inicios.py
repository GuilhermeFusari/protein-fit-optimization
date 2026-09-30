#!/usr/bin/env python3
"""
efeito_inicios.py  -  efeito do numero de inicios (--restarts) no benchmark
das 50 entradas, comparado com o cifsup NSD.

Nao altera nada no pacote nem nos CSVs publicados: roda saxs_icp.align com
os parametros publicados, variando apenas `restarts`, e le (sem modificar)
benchmark/benchmark_cifsup_lam02.csv para os valores do cifsup NSD.

Para cada (restarts, semente-base) grava uma linha por entrada com:
  chamfer          como no artigo: contra o envelope amostrado a 5000 pts
                   (identico a coluna chamfer de final_lam02.csv)
  chamfer_envfull  contra o envelope completo, como o benchmark_cifsup.py
                   calcula o Chamfer do cifsup (difere so em envelopes com
                   mais de 5000 pontos)
  frac_fora, cobertura, sep_centroide, custo, espelhado

Semente-base 42 = a publicada (entrada i recebe 42 + i). As sementes 1042 e
2042 medem quanto o resultado depende do sorteio das orientacoes.

Checagem embutida: restarts=3, semente 42 tem de reproduzir final_lam02.csv.

Uso:
    python3 efeito_inicios.py <pasta com data/sasbdb/> <repo> saida.csv
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

RESTARTS = (3, 10, 30, 100)
SEEDS = (42, 1042, 2042)


def job(args):
    acc, env_path, mod_path, seed, restarts = args
    t0 = time.time()
    env_full = read_coords(env_path)
    mod = read_coords(mod_path, ca_only=True)
    r = core.align(mod, env_full, restarts=restarts, seed=seed)
    fit, m = r.fitted, r.metrics
    return dict(
        entry=acc, restarts=restarts, seed=seed, custo=round(r.cost, 4),
        chamfer=round(m["chamfer"], 4),
        chamfer_envfull=round(core.chamfer(fit, env_full), 4),
        frac_fora=round(m["frac_outside"], 4),
        cobertura=round(m["coverage"], 4),
        sep_centroide=round(m["centroid_sep"], 3),
        espelhado=bool(r.mirrored), n_env_full=len(env_full), n_env=r.n_env,
        segundos=round(time.time() - t0, 2))


def main():
    base, repo, out = sys.argv[1:4]
    import pandas as pd

    rows = []
    for seed_base in SEEDS:
        jobs0 = find_jobs(base, seed=seed_base)
        for R in RESTARTS:
            jobs = [(j[0], j[1], j[2], j[8], R) for j in jobs0]
            t0 = time.time()
            with Pool() as pool:
                res = pool.map(job, jobs)
            for r in res:
                r["seed_base"] = seed_base
            rows += res
            ch = np.mean([r["chamfer"] for r in res])
            print(f"seed_base={seed_base} restarts={R:3d}: chamfer medio {ch:.4f} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows).sort_values(["seed_base", "restarts", "entry"])
    cols = ["entry", "seed_base", "seed", "restarts", "custo", "chamfer",
            "chamfer_envfull", "frac_fora", "cobertura", "sep_centroide",
            "espelhado", "n_env_full", "n_env", "segundos"]
    df[cols].to_csv(out, index=False)

    # checagem: restarts=3, semente 42 == final_lam02.csv
    ref = pd.read_csv(os.path.join(repo, "benchmark", "final_lam02.csv")).set_index("entry")
    a = df[(df.seed_base == 42) & (df.restarts == 3)].set_index("entry").loc[ref.index]
    c = ["custo", "chamfer", "frac_fora", "cobertura", "sep_centroide", "espelhado"]
    ndiff = int((a[c] != ref[c]).sum().sum())
    print(f"checagem restarts=3 semente 42 vs final_lam02.csv: {ndiff} diferencas")
    if ndiff:
        sys.exit("ERRO: nao reproduz o CSV publicado")


if __name__ == "__main__":
    main()
