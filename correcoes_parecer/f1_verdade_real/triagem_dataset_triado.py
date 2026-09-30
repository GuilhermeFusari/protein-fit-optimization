#!/usr/bin/env python3
"""
triagem_dataset_triado.py  -  mesma investigacao de investigar_referencial.py
(so leitura) em todos os pares (envelope, modelo de fit) de
Dataset_Triado/3_Fila_Dammin, com 200 rotacoes aleatorias por par (triagem).

Tambem registra, a partir do texto do envelope, o software de superposicao
anotado (SUPCOMB, CIFSUP, ALPRAXIN, ...), porque uma pose superposta por
SUPCOMB/cifsup ou por eixos principais (ALPRAXIN) nao serve de verdade
independente.

Uso:
    python3 triagem_dataset_triado.py <Dataset_Triado/3_Fila_Dammin> saida.csv
"""
import glob
import os
import re
import sys
from multiprocessing import Pool

import investigar_referencial as ir

ir.N_ROT = 200
SUP = re.compile(r"SUPCOMB|CIFSUP|ALPRAXIN|DAMSUP|SUPALM|superimpos", re.I)


def software(env_path):
    try:
        txt = open(env_path, errors="replace").read(200000)
    except OSError:
        return ""
    return ";".join(sorted({m.group(0).upper() for m in SUP.finditer(txt)}))


def job(args):
    acc, env, mod = args
    try:
        r = ir.job((acc, env, mod))
    except Exception as e:
        r = dict(entry=acc, status=f"erro: {type(e).__name__}: {e}")
    r["modelo"] = os.path.basename(mod)
    r["registro_superposicao"] = software(env)
    return r


def main():
    import pandas as pd
    root, out = sys.argv[1:3]
    jobs = []
    for env in sorted(glob.glob(os.path.join(root, "*", "*_envelope.cif"))):
        acc = os.path.basename(env)[:-13]
        for mod in sorted(glob.glob(os.path.join(root, acc, acc, f"{acc}_fit*_model*"))):
            jobs.append((acc, env, mod))
    print(f"{len(jobs)} pares", flush=True)
    with Pool() as pool:
        rows = pool.map(job, jobs, chunksize=4)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
