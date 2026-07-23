#!/usr/bin/env python3
"""
benchmark_cifsup.py

Compara o pipeline ICP-SAXS (penalidade bidirecional) com o cifsup do
ATSAS nas mesmas entradas e com o mesmo pre-processamento.

Para cada entrada roda:
    cifsup --method=NSD    equivalente moderno do SUPCOMB
    cifsup --method=ICP    ICP do ATSAS, sem restricao volumetrica

e le o resultado do pipeline proprio a partir do CSV gerado por
icp_saxs_bidirectional.py.

Todas as metricas comparaveis (chamfer, frac_fora, sep_centroide) sao
RECALCULADAS aqui com o mesmo codigo para as tres solucoes. Os scores
nativos nao sao comparaveis entre si: o NSD e adimensional e o score do
metodo ICP do cifsup esta em A.

Uso:
    python3 benchmark_cifsup.py \\
        --base ICP_SAXS_Final_50/ICP_SAXS_Final_50 \\
        --ours resultados_author_lam2.csv \\
        --out benchmark_cifsup.csv

Requer cifsup no PATH (ATSAS 3.2+).

Detalhes do cifsup 3.2.1 descobertos empiricamente:
  - decide o parser pela EXTENSAO do arquivo; os envelopes do SASBDB tem
    extensao .cif mas conteudo PDB, e falham com exit=5 se passados
    diretamente. O script copia para .pdb antes de chamar.
  - nao imprime nada em stdout; o score sai no cabecalho do PDB de saida,
    na linha "REMARK 265 score".
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.spatial import Delaunay, cKDTree

SCORE_RE = re.compile(r"REMARK 265 score\s*:\s*([-+]?\d+(?:\.\d+)?)", re.I)
TIMEOUT = 600


# --------------------------------------------------------------------------
# leitura (PDB e mmCIF, detectados pelo conteudo, nao pela extensao)
# --------------------------------------------------------------------------

def _pdb_fmt(path, ca_only):
    pts = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if ca_only and line[12:16].strip() != "CA":
                continue
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return pts


def _cif_fmt(path, ca_only):
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
    ia = cols.get("label_atom_id", cols.get("auth_atom_id"))
    pts = []
    for p in rows:
        if len(p) <= max(ix, iy, iz):
            continue
        if ca_only and ia is not None and len(p) > ia:
            if p[ia].strip('"').strip("'") != "CA":
                continue
        try:
            pts.append((float(p[ix]), float(p[iy]), float(p[iz])))
        except ValueError:
            continue
    return pts


def read_coords(path, ca_only=False):
    pts = _pdb_fmt(path, ca_only) or _cif_fmt(path, ca_only)
    if not pts and ca_only:
        pts = _pdb_fmt(path, False) or _cif_fmt(path, False)
    return np.asarray(pts, dtype=float)


def write_pdb(pts, path, resname="ALA"):
    with open(path, "w") as f:
        for i, (x, y, z) in enumerate(pts, 1):
            f.write(f"ATOM  {i:5d}  CA  {resname} A{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("END\n")


# --------------------------------------------------------------------------
# metricas (identicas as do icp_saxs_bidirectional.py)
# --------------------------------------------------------------------------

def chamfer(a, b):
    return float((cKDTree(b).query(a)[0].mean() + cKDTree(a).query(b)[0].mean()) / 2.0)


def frac_outside(pts, env):
    try:
        return float((Delaunay(env).find_simplex(pts) < 0).mean())
    except Exception:
        return float("nan")


def coverage(fit, env, thr):
    return float((cKDTree(fit).query(env)[0] <= thr).mean())


def grid_spacing(env):
    if len(env) < 2:
        return 3.0
    d, _ = cKDTree(env).query(env, k=2)
    return float(np.median(d[:, 1]))


# --------------------------------------------------------------------------

def run_cifsup(method, template, movable, out_pdb):
    cmd = ["cifsup", f"--method={method}", "-o", str(out_pdb), str(template), str(movable)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return None, TIMEOUT, "timeout"
    except FileNotFoundError:
        print("ERRO: 'cifsup' nao encontrado no PATH.", file=sys.stderr)
        sys.exit(1)
    dt = time.time() - t0

    if p.returncode != 0 or not os.path.exists(out_pdb):
        return None, dt, f"exit={p.returncode}"

    score = None
    with open(out_pdb, "r", errors="replace") as f:
        for line in f:
            if not line.startswith("REMARK"):
                break
            m = SCORE_RE.search(line)
            if m:
                score = float(m.group(1))
                break
    return score, dt, "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="pasta com data/sasbdb/")
    ap.add_argument("--ours", help="CSV do icp_saxs_bidirectional.py (opcional)")
    ap.add_argument("--out", default="benchmark_cifsup.csv")
    ap.add_argument("--max-points", type=int, default=300)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    sasbdb = os.path.join(os.path.expanduser(args.base), "data", "sasbdb")
    if not os.path.isdir(sasbdb):
        print(f"nao encontrei {sasbdb}")
        return 1

    entries = sorted(os.path.basename(p).replace("_envelope.cif", "")
                     for p in glob.glob(os.path.join(sasbdb, "*_envelope.cif")))
    if args.limit:
        entries = entries[:args.limit]

    ours = {}
    if args.ours and os.path.exists(args.ours):
        import pandas as pd
        for _, r in pd.read_csv(args.ours).iterrows():
            ours[r["entry"]] = r

    print(f"{len(entries)} entradas\n")
    rows = []

    for i, acc in enumerate(entries, 1):
        env_src = os.path.join(sasbdb, f"{acc}_envelope.cif")
        mod_src = os.path.join(sasbdb, f"{acc}_model.pdb")
        if not os.path.exists(mod_src):
            continue

        env = read_coords(env_src)
        mod = read_coords(mod_src, ca_only=True)
        if len(mod) == 0:
            mod = read_coords(mod_src)
        if len(env) < 10 or len(mod) < 4:
            rows.append(dict(entry=acc, status="coordenadas_vazias"))
            print(f"[{i}/{len(entries)}] {acc}: coordenadas vazias")
            continue

        # mesmo downsampling do pipeline proprio, mesma semente
        rng = np.random.default_rng(args.seed + i - 1)
        if len(mod) > args.max_points:
            mod = mod[np.sort(rng.choice(len(mod), args.max_points, replace=False))]

        thr = max(3.0, grid_spacing(env))
        r = dict(entry=acc, status="ok", n_env=len(env), n_model=len(mod))

        tmp = tempfile.mkdtemp(prefix=f"cs_{acc}_")
        try:
            env_pdb = os.path.join(tmp, "envelope.pdb")
            mod_pdb = os.path.join(tmp, "model.pdb")
            shutil.copy(env_src, env_pdb)
            write_pdb(mod, mod_pdb)

            for method in ("NSD", "ICP"):
                out = os.path.join(tmp, f"fit_{method}.pdb")
                score, dt, st = run_cifsup(method, env_pdb, mod_pdb, out)
                k = method.lower()
                r[f"{k}_score"] = score
                r[f"{k}_secs"] = round(dt, 2)
                r[f"{k}_status"] = st
                if st == "ok":
                    fit = read_coords(out, ca_only=True)
                    if len(fit) == 0:
                        fit = read_coords(out)
                    if len(fit):
                        r[f"{k}_chamfer"] = round(chamfer(fit, env), 3)
                        r[f"{k}_frac_fora"] = round(frac_outside(fit, env), 4)
                        r[f"{k}_cobertura"] = round(coverage(fit, env, thr), 4)
                        r[f"{k}_sep_centroide"] = round(
                            float(np.linalg.norm(fit.mean(0) - env.mean(0))), 3)

            if acc in ours:
                o = ours[acc]
                for c in ("chamfer", "frac_fora", "cobertura", "sep_centroide"):
                    if c in o:
                        r[f"ours_{c}"] = o[c]

            print(f"[{i}/{len(entries)}] {acc:9s} "
                  f"NSD={r.get('nsd_score')}  "
                  f"cham: nsd={r.get('nsd_chamfer')} icp={r.get('icp_chamfer')} "
                  f"ours={r.get('ours_chamfer')}")
            rows.append(r)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df[df.get("status") == "ok"] if "status" in df else df

    print(f"\n{len(ok)}/{len(df)} concluidas -> {args.out}")
    if len(ok):
        print("\n" + "=" * 66)
        print(f"{'metrica':18s}{'cifsup NSD':>14s}{'cifsup ICP':>14s}{'ICP-SAXS':>14s}")
        print("-" * 66)
        for met in ("chamfer", "frac_fora", "cobertura", "sep_centroide"):
            vals = []
            for pre in ("nsd", "icp", "ours"):
                c = f"{pre}_{met}"
                vals.append(f"{ok[c].mean():14.3f}" if c in ok and ok[c].notna().any()
                            else f"{'n/d':>14s}")
            print(f"{met:18s}" + "".join(vals))
        print("\nO NSD e adimensional e nao entra nesta tabela; as colunas acima")
        print("sao todas recalculadas com o mesmo codigo, logo comparaveis.")

        try:
            from scipy import stats
            print("\nWilcoxon pareado (ICP-SAXS vs cifsup ICP):")
            for met in ("chamfer", "frac_fora", "sep_centroide"):
                a, b = f"ours_{met}", f"icp_{met}"
                if a in ok and b in ok:
                    m = ok[[a, b]].dropna()
                    if len(m) >= 5:
                        _, p = stats.wilcoxon(m[a], m[b])
                        melhor = "ICP-SAXS" if m[a].mean() < m[b].mean() else "cifsup"
                        print(f"  {met:16s} p={p:.4f}  ({melhor} melhor)")
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
