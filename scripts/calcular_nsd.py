#!/usr/bin/env python3
"""
calcular_nsd.py

Tabela cruzada Chamfer x NSD para as tres solucoes (ICP-SAXS, cifsup NSD,
cifsup ICP) sobre as MESMAS poses ja alinhadas, com a MESMA funcao.

Motivacao (parecer, ponto A.1): a comparacao anterior usava o Chamfer,
que e a propria funcao objetivo do pipeline, o que favorece o pipeline
por construcao. Aqui reportamos as duas metricas lado a lado, cada uma
calculada identicamente para os tres metodos. Se o cifsup vencer no NSD
(a metrica dele) e empatar no Chamfer, a leitura honesta e "cada metodo
se sai melhor na propria metrica, e sao equivalentes na do outro".

NSD (Kozin & Svergun, 2001), forma simetrica:

    rho(A->B) = sqrt( (1/N_A) * sum_i min_j |a_i - b_j|^2 / d_B^2 )

    NSD(A,B)  = sqrt( ( rho(A->B)^2 + rho(B->A)^2 ) / 2 )

onde d_A e d_B sao os espacamentos tipicos entre pontos vizinhos de cada
conjunto (mediana da distancia ao vizinho mais proximo). A normalizacao
por d torna o NSD adimensional: ~1 para formas iguais na resolucao dos
modelos, >1 para formas sistematicamente diferentes.

Uso:
    python3 calcular_nsd.py \\
        --base ICP_SAXS_Final_50/ICP_SAXS_Final_50 \\
        --ours poses_ours --cifsup poses_cifsup \\
        --out tabela_chamfer_nsd.csv

Espera:
    poses_ours/<ENTRY>_final.pdb
    poses_cifsup/<ENTRY>_nsd.pdb
    poses_cifsup/<ENTRY>_icp.pdb
    <base>/data/sasbdb/<ENTRY>_envelope.cif
"""

import argparse
import glob
import os
import sys

import numpy as np
from scipy.spatial import cKDTree


# --------------------------------------------------------------------------
# leitura (PDB ou mmCIF, detectado pelo conteudo, nao pela extensao)
# --------------------------------------------------------------------------

def _pdb(path):
    pts = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    pts.append((float(line[30:38]), float(line[38:46]),
                                float(line[46:54])))
                except ValueError:
                    pass
    return pts


def _cif(path):
    cols, rows, in_loop = {}, [], False
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols[s.split(".", 1)[1].split()[0]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if s.startswith("#") or s.startswith("loop_") or not s:
                    if rows:
                        break
                    continue
                if s.startswith("_"):
                    break
                rows.append(s.split())
    need = ("Cartn_x", "Cartn_y", "Cartn_z")
    if not all(k in cols for k in need):
        return []
    ix, iy, iz = (cols[k] for k in need)
    pts = []
    for p in rows:
        if len(p) > max(ix, iy, iz):
            try:
                pts.append((float(p[ix]), float(p[iy]), float(p[iz])))
            except ValueError:
                pass
    return pts


def read_coords(path):
    pts = _pdb(path) or _cif(path)
    return np.asarray(pts, dtype=float)


# --------------------------------------------------------------------------
# metricas
# --------------------------------------------------------------------------

def spacing(pts):
    """Espacamento tipico: mediana da distancia ao vizinho mais proximo."""
    if len(pts) < 2:
        return 1.0
    d, _ = cKDTree(pts).query(pts, k=2)
    s = float(np.median(d[:, 1]))
    return s if s > 1e-6 else 1.0


def nsd(A, B):
    """NSD simetrico de Kozin & Svergun (2001)."""
    dA, dB = spacing(A), spacing(B)
    dAB, _ = cKDTree(B).query(A, k=1)     # A -> B
    dBA, _ = cKDTree(A).query(B, k=1)     # B -> A
    rho_AB = np.sqrt(np.mean((dAB / dB) ** 2))
    rho_BA = np.sqrt(np.mean((dBA / dA) ** 2))
    return float(np.sqrt((rho_AB ** 2 + rho_BA ** 2) / 2.0))


def chamfer(A, B):
    return float((cKDTree(B).query(A)[0].mean() + cKDTree(A).query(B)[0].mean()) / 2.0)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--ours", default="poses_ours")
    ap.add_argument("--cifsup", default="poses_cifsup")
    ap.add_argument("--out", default="tabela_chamfer_nsd.csv")
    args = ap.parse_args()

    sasbdb = os.path.join(os.path.expanduser(args.base), "data", "sasbdb")
    if not os.path.isdir(sasbdb):
        print(f"nao encontrei {sasbdb}")
        return 1

    entries = sorted(os.path.basename(p).replace("_final.pdb", "")
                     for p in glob.glob(os.path.join(args.ours, "*_final.pdb")))
    if not entries:
        print(f"nenhuma pose em {args.ours}/*_final.pdb")
        return 1

    rows = []
    for acc in entries:
        env_path = os.path.join(sasbdb, f"{acc}_envelope.cif")
        p_ours = os.path.join(args.ours, f"{acc}_final.pdb")
        p_nsd = os.path.join(args.cifsup, f"{acc}_nsd.pdb")
        p_icp = os.path.join(args.cifsup, f"{acc}_icp.pdb")

        if not all(os.path.exists(p) for p in (env_path, p_ours, p_nsd, p_icp)):
            rows.append(dict(entry=acc, status="pose_faltando"))
            continue

        env = read_coords(env_path)
        poses = {"ours": read_coords(p_ours),
                 "nsd": read_coords(p_nsd),
                 "icp": read_coords(p_icp)}
        if len(env) < 10 or any(len(v) < 4 for v in poses.values()):
            rows.append(dict(entry=acc, status="coordenadas_vazias"))
            continue

        r = dict(entry=acc, status="ok", n_env=len(env))
        for k, pts in poses.items():
            r[f"{k}_chamfer"] = round(chamfer(pts, env), 4)
            r[f"{k}_nsd"] = round(nsd(pts, env), 4)
        rows.append(r)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df[df.status == "ok"]
    print(f"{len(ok)}/{len(df)} entradas com as tres poses -> {args.out}\n")
    if not len(ok):
        return 0

    # tabela cruzada
    print("=" * 58)
    print(f"{'':16s}{'ICP-SAXS':>12s}{'cifsup NSD':>14s}{'cifsup ICP':>14s}")
    print("-" * 58)
    labels = {"chamfer": "Chamfer (A)", "nsd": "NSD (adim.)"}
    for met in ("chamfer", "nsd"):
        line = f"{labels[met]:16s}"
        for pre in ("ours", "nsd", "icp"):
            line += f"{ok[f'{pre}_{met}'].mean():>{12 if pre=='ours' else 14}.3f}"
        print(line)
    print("=" * 58)

    # Wilcoxon pareado: ICP-SAXS vs cada modo do cifsup, nas duas metricas
    try:
        from scipy import stats
        print("\nWilcoxon pareado (ICP-SAXS vs cifsup):")
        for met in ("chamfer", "nsd"):
            for pre, nome in (("nsd", "cifsup NSD"), ("icp", "cifsup ICP")):
                a = ok[f"ours_{met}"]
                b = ok[f"{pre}_{met}"]
                _, p = stats.wilcoxon(a, b)
                vencedor = "ICP-SAXS" if a.mean() < b.mean() else nome
                print(f"  {met:8s} vs {nome:11s} p={p:.4f}  ({vencedor} menor)")
    except Exception as e:
        print("Wilcoxon indisponivel:", e)

    print("\nLeitura: cada metodo tende a se sair melhor na metrica que")
    print("otimiza (ICP-SAXS no Chamfer, cifsup NSD no NSD). O que importa")
    print("e se as diferencas sao significativas nas DUAS metricas.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
