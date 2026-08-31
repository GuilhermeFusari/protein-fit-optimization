#!/usr/bin/env python3
"""
validacao_lambda.py

Escolha do peso lambda por validacao cruzada leave-one-out (LOO), a
partir de uma varredura ja calculada (um CSV por lambda, uma linha por
entrada).

Motivacao (parecer, ponto A.4). A calibracao anterior tinha tres falhas:
  1. foi feita com downsampling de 300 pontos, o regime degenerado;
  2. so testou lambda >= 1, com o otimo aparente colado no limite;
  3. escolheu lambda e mediu desempenho nas MESMAS 50 entradas, o que
     e ajuste no proprio conjunto de teste.

Aqui:
  - a varredura de entrada usa a representacao completa (max-points alto);
  - a grade cobre lambda < 1 (escala log);
  - LOO: para cada entrada, o melhor lambda e escolhido olhando as OUTRAS
    n-1 entradas, e o resultado e medido na entrada deixada de fora. Assim
    nenhuma entrada participa da escolha do lambda que a avalia.

O criterio de selecao (qual metrica minimizar) e explicito: por padrao a
distancia de Chamfer, que e a metrica primaria reportada. Pode ser trocado
com --criterio.

Uso:
    # 1) varredura, um lambda por arquivo (rodar o alinhador antes):
    #    for L in 0.1 0.2 0.5 1 2 5 10 20 50 100; do
    #      python3 icp_saxs_v2.py --base ... --mode author --penalty $L \\
    #          --max-points 3000 --out sweep/lam_$L.csv
    #    done
    #
    # 2) validacao cruzada sobre os CSVs da varredura:
    python3 validacao_lambda.py --sweep "sweep/lam_*.csv" \\
        --criterio chamfer --out loo_lambda.csv

Cada CSV precisa das colunas: entry, e a metrica do criterio.
O lambda de cada arquivo e lido do nome (lam_<valor>.csv) ou de uma
coluna 'penalty'/'lam' se existir.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd


def lambda_do_nome(path):
    m = re.search(r"lam[_-]?([0-9]*\.?[0-9]+)", os.path.basename(path))
    return float(m.group(1)) if m else None


def carregar(sweep_glob, criterio):
    arquivos = sorted(glob.glob(sweep_glob))
    if not arquivos:
        print(f"nenhum arquivo em {sweep_glob}")
        sys.exit(1)

    tabela = {}          # lambda -> {entry: valor_do_criterio}
    for path in arquivos:
        df = pd.read_csv(path)
        if "status" in df.columns:
            df = df[df.status == "ok"]
        if criterio not in df.columns:
            print(f"{path}: sem coluna '{criterio}'")
            sys.exit(1)

        lam = None
        for c in ("penalty", "lam", "lambda"):
            if c in df.columns and df[c].nunique() == 1:
                lam = float(df[c].iloc[0])
                break
        if lam is None:
            lam = lambda_do_nome(path)
        if lam is None:
            print(f"{path}: nao consegui determinar o lambda")
            sys.exit(1)

        tabela[lam] = dict(zip(df["entry"], df[criterio]))

    lambdas = sorted(tabela)
    # entradas presentes em TODOS os lambdas (comparacao pareada)
    comuns = set.intersection(*(set(tabela[l]) for l in lambdas))
    comuns = sorted(comuns)
    if not comuns:
        print("nenhuma entrada presente em todos os lambdas")
        sys.exit(1)

    # matriz [entradas x lambdas] do valor do criterio
    M = np.array([[tabela[l][e] for l in lambdas] for e in comuns], dtype=float)
    return np.array(lambdas), np.array(comuns), M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True,
                    help="glob dos CSVs da varredura, ex: 'sweep/lam_*.csv'")
    ap.add_argument("--criterio", default="chamfer",
                    help="metrica a minimizar (default: chamfer)")
    ap.add_argument("--out", default="loo_lambda.csv")
    args = ap.parse_args()

    lambdas, entradas, M = carregar(args.sweep, args.criterio)
    n = len(entradas)

    # --- 1. otimo "ingenuo": escolhe lambda com todas as entradas ---
    media_por_lambda = M.mean(axis=0)
    i_naive = int(np.argmin(media_por_lambda))
    lam_naive = lambdas[i_naive]
    valor_naive = media_por_lambda[i_naive]

    # --- 2. leave-one-out ---
    # para cada entrada i: escolhe o lambda que minimiza a media nas OUTRAS
    # entradas; mede o criterio na entrada i com esse lambda.
    escolhidos = np.empty(n)
    medidos = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        media_sem_i = M[mask].mean(axis=0)
        j = int(np.argmin(media_sem_i))
        escolhidos[i] = lambdas[j]
        medidos[i] = M[i, j]

    loo_media = float(medidos.mean())

    # --- salva resultados por entrada ---
    df = pd.DataFrame({
        "entry": entradas,
        "lambda_escolhido_loo": escolhidos,
        f"{args.criterio}_loo": medidos,
        f"{args.criterio}_lambda_naive": M[:, i_naive],
    })
    df.to_csv(args.out, index=False)

    # --- relatorio ---
    print(f"criterio minimizado: {args.criterio}")
    print(f"entradas: {n} | lambdas: {list(lambdas)}\n")

    print("media do criterio por lambda (todas as entradas):")
    for l, v in zip(lambdas, media_por_lambda):
        marca = "  <- minimo global" if l == lam_naive else ""
        print(f"  lambda={l:6.2f}   {args.criterio}={v:.4f}{marca}")

    print(f"\nEscolha ingenua (calibra e mede em tudo):")
    print(f"  lambda={lam_naive:g}  {args.criterio}={valor_naive:.4f}")
    print(f"  -> este e o numero 'otimista', com vies de conjunto de teste.")

    vals, cont = np.unique(escolhidos, return_counts=True)
    print(f"\nLeave-one-out:")
    print(f"  lambda escolhido por entrada (frequencia):")
    for v, c in sorted(zip(vals, cont), key=lambda t: -t[1]):
        print(f"    lambda={v:6.2f}  escolhido {c}/{n} vezes")
    print(f"  {args.criterio} medio honesto (LOO): {loo_media:.4f}")

    dif = loo_media - valor_naive
    print(f"\n  diferenca LOO - ingenuo: {dif:+.4f} "
          f"({100*dif/valor_naive:+.1f}%)")
    if abs(dif) < 0.02 * valor_naive:
        print("  -> desprezivel: o lambda era robusto, nao havia vies "
              "relevante.")
    else:
        print("  -> relevante: parte do desempenho reportado vinha da "
              "escolha do lambda no proprio conjunto.")

    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
