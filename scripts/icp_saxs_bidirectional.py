#!/usr/bin/env python3
"""
rerun_single.py

Re-executa o alinhamento das entradas do benchmark implementando a
metodologia DESCRITA NO ARTIGO, que o modo single do codigo antigo nao
implementava:

  - representacao por atomos CA com downsampling para no maximo 300 pontos
  - penalidade volumetrica assimetrica (lambda configuravel)
  - N reinicializacoes aleatorias, ficando com a de menor custo
  - envelope amostrado para 5000 pontos

Roda em dois bracos para permitir o ablation:
    --penalty 0    -> ICP puro (reproduz o comportamento antigo)
    --penalty 50   -> ICP com restricao volumetrica

Metricas reportadas por entrada (todas recalculadas com o mesmo codigo,
independentemente do custo usado na otimizacao):
    chamfer, hausdorff, dice, frac_fora, sep_centroide, d_med

Uso:
    python3 rerun_single.py --base ICP_SAXS_Final_50/ICP_SAXS_Final_50 \\
        --penalty 50 --out resultados_lambda50.csv

Reaproveita a logica de icp_lib.py (penalidade assimetrica + rotacao
aleatoria inicial); nao depende do SAXS_Protein_Aligner instalado.
"""

import argparse
import glob
import os
import sys
import time
from multiprocessing import Pool, cpu_count

import numpy as np
from scipy.spatial import cKDTree, Delaunay
from scipy.spatial.transform import Rotation


# --------------------------------------------------------------------------
# leitura
# --------------------------------------------------------------------------

def _read_pdb_format(path, ca_only):
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


def _read_cif_format(path, ca_only):
    """
    mmCIF: le o loop _atom_site pelas posicoes declaradas no cabecalho,
    em vez de assumir uma ordem fixa de colunas.
    """
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
                    in_loop = bool(cols)
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
    for parts in rows:
        if len(parts) <= max(ix, iy, iz):
            continue
        if ca_only and ia is not None and len(parts) > ia:
            if parts[ia].strip('"').strip("'") != "CA":
                continue
        try:
            pts.append((float(parts[ix]), float(parts[iy]), float(parts[iz])))
        except ValueError:
            continue
    return pts


def read_pdb(path, ca_only=False):
    """
    Le coordenadas de PDB ou mmCIF.

    Os arquivos do SASBDB tem extensao .cif mas parte deles traz conteudo
    PDB (dummy atoms) e parte traz mmCIF real. O formato e detectado pelo
    conteudo, nao pela extensao: ler o segundo com um parser de PDB devolve
    zero coordenadas silenciosamente.
    """
    pts = _read_pdb_format(path, ca_only)
    if not pts:
        pts = _read_cif_format(path, ca_only)
    if not pts and ca_only:          # sem CA identificavel: usa tudo
        pts = _read_pdb_format(path, False) or _read_cif_format(path, False)
    return np.asarray(pts, dtype=float)


def grid_spacing(env):
    """Espacamento tipico entre dummy atoms vizinhos."""
    if len(env) < 2:
        return 3.0
    d, _ = cKDTree(env).query(env, k=2)
    return float(np.median(d[:, 1]))


# --------------------------------------------------------------------------
# custo e ICP
# --------------------------------------------------------------------------

def bidirectional_cost(src, env_tree, src_tree_pts, env_pts, penalty, threshold,
                       w_leak=1.0, w_cover=1.0):
    """
    Custo bidirecional, como descrito pelo autor:

      termo de vazamento : atomos longe do envelope (proteina saindo)
      termo de cobertura : pontos do envelope longe de qualquer atomo
                           (espaco vazio dentro do envelope)

    O termo de cobertura e o que falta na weighted_asymmetric_score
    original: sem ele nada empurra a estrutura para PREENCHER o volume,
    so para nao sair dele. Com uma copia so, o termo de vazamento
    sozinho e satisfeito por qualquer pose adjacente ao envelope.
    """
    d_src, _ = env_tree.query(src, k=1)
    base = float(d_src.mean())
    if penalty <= 0:
        return base, base, 0.0, 0.0

    over = d_src[d_src > threshold]
    leak = float((over - threshold).sum()) / len(src) if len(over) else 0.0

    src_tree = cKDTree(src)
    d_env, _ = src_tree.query(env_pts, k=1)
    unc = d_env[d_env > threshold]
    cover = float((unc - threshold).sum()) / len(env_pts) if len(unc) else 0.0

    total = base + penalty * (w_leak * leak + w_cover * cover)
    return total, base, leak, cover


def author_cost(src, env_pts, env_tree):
    """
    Formulacao original do autor (commit 6b50d10):

        score = media(dist envelope->proteina)          # preenchimento
              + penalty * media(dist proteina->envelope) # vazamento

    Medias simples, sem limiar: todo atomo contribui, e o vazamento
    custa `penalty` vezes mais que deixar espaco vazio. Mais simples que
    a versao por excesso-acima-do-limiar e sem parametro de limiar para
    calibrar.
    """
    d_fill, _ = cKDTree(src).query(env_pts, k=1)
    d_leak, _ = env_tree.query(src, k=1)
    return float(d_fill.mean()), float(d_leak.mean())


def asymmetric_cost(src, env_tree, penalty, threshold):
    """Versao original, unidirecional: apenas vazamento."""
    d, _ = env_tree.query(src, k=1)
    base = float(d.mean())
    if penalty <= 0:
        return base
    over = d[d > threshold]
    if len(over) == 0:
        return base
    return base + penalty * float((over - threshold).sum()) / len(src)


def icp(src_pts, env_pts, env_tree, max_iter, penalty, threshold,
        mode="bidir", w_leak=1.0, w_cover=1.0):
    """
    ICP com custo bidirecional opcional.

    Os pesos entram no passo de Procrustes, nao apenas na selecao da
    melhor iteracao: um Procrustes nao ponderado produz exatamente o
    mesmo movimento com ou sem penalidade, que e por que a versao
    original nao alterava o alinhamento.

    - atomos que vazam recebem peso maior -> puxados para dentro
    - regioes vazias do envelope entram como alvos extras -> a estrutura
      e atraida para preencher o volume
    """
    src = src_pts.copy()
    T_total = np.eye(4)
    best_cost, best_T = float("inf"), np.eye(4)

    for _ in range(max_iter):
        dist, idx = env_tree.query(src, k=1)
        closest = env_pts[idx]

        if mode == "author":
            fill, leak = author_cost(src, env_pts, env_tree)
            cost = fill + penalty * leak
        elif mode == "bidir":
            cost, _, _, _ = bidirectional_cost(src, env_tree, None, env_pts,
                                               penalty, threshold, w_leak, w_cover)
        else:
            cost = asymmetric_cost(src, env_tree, penalty, threshold)
        if cost < best_cost:
            best_cost, best_T = cost, T_total.copy()

        if mode == "author" and penalty > 0:
            # sem limiar: o peso cresce com a propria distancia
            w = 1.0 + penalty * dist / max(float(dist.mean()), 1e-9)
        elif penalty > 0:
            w = 1.0 + penalty * np.maximum(dist - threshold, 0.0) / max(threshold, 1e-9)
        else:
            w = np.ones(len(src))

        srcA, tgtA, wA = src, closest, w

        if mode == "author" and penalty > 0:
            # termo de preenchimento: TODO ponto do envelope e alvo,
            # ponderado pela distancia ate a proteina mais proxima
            src_tree = cKDTree(src)
            d_env, i_env = src_tree.query(env_pts, k=1)
            srcA = np.vstack([src, src[i_env]])
            tgtA = np.vstack([closest, env_pts])
            w_extra = d_env / max(float(d_env.mean()), 1e-9)
            wA = np.concatenate([w, w_extra])
        elif mode == "bidir" and penalty > 0 and w_cover > 0:
            # pontos do envelope sem atomo por perto: alvos extras que
            # puxam a estrutura para as regioes vazias
            src_tree = cKDTree(src)
            d_env, i_env = src_tree.query(env_pts, k=1)
            unc = d_env > threshold
            if unc.any():
                src_extra = src[i_env[unc]]
                tgt_extra = env_pts[unc]
                w_extra = (w_cover * penalty *
                           (d_env[unc] - threshold) / max(threshold, 1e-9))
                srcA = np.vstack([src, src_extra])
                tgtA = np.vstack([closest, tgt_extra])
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

    if mode == "author":
        fill, leak = author_cost(src, env_pts, env_tree)
        cost = fill + penalty * leak
    elif mode == "bidir":
        cost, _, _, _ = bidirectional_cost(src, env_tree, None, env_pts,
                                           penalty, threshold, w_leak, w_cover)
    else:
        cost = asymmetric_cost(src, env_tree, penalty, threshold)
    if cost < best_cost:
        best_cost, best_T = cost, T_total.copy()
    return best_T, best_cost


def apply_T(pts, T):
    return (T[:3, :3] @ pts.T).T + T[:3, 3]


# --------------------------------------------------------------------------
# metricas de avaliacao
# --------------------------------------------------------------------------

def chamfer(a, b):
    return float((cKDTree(b).query(a)[0].mean() + cKDTree(a).query(b)[0].mean()) / 2.0)


def hausdorff(a, b):
    return float(max(cKDTree(b).query(a)[0].max(), cKDTree(a).query(b)[0].max()))


def dice(a, b, voxel=1.0):
    """Sobreposicao volumetrica por voxelizacao numa grade comum."""
    allp = np.vstack([a, b])
    origin = allp.min(axis=0)
    ka = set(map(tuple, np.floor((a - origin) / voxel).astype(int)))
    kb = set(map(tuple, np.floor((b - origin) / voxel).astype(int)))
    if not ka or not kb:
        return 0.0
    return 2.0 * len(ka & kb) / (len(ka) + len(kb))


def frac_outside(pts, env):
    """Fracao fora do casco convexo do envelope. Limite inferior da
    violacao real, ja que o casco e maior que um envelope concavo."""
    try:
        return float((Delaunay(env).find_simplex(pts) < 0).mean())
    except Exception:
        return float("nan")


# --------------------------------------------------------------------------
# uma entrada
# --------------------------------------------------------------------------

def process(job):
    (acc, env_path, mod_path, penalty, restarts, max_iter,
     max_pts, sample_env, seed, mode, w_leak, w_cover) = job
    t0 = time.time()
    try:
        env_full = read_pdb(env_path)
        mod = read_pdb(mod_path, ca_only=True)
        if len(mod) == 0:                      # sem CA: usa todos os atomos
            mod = read_pdb(mod_path)
        if len(env_full) < 10 or len(mod) < 4:
            return dict(entry=acc, status="coordenadas_vazias")

        rng = np.random.default_rng(seed)

        # downsampling do modelo para 300 pontos, como descrito no artigo
        n_model_full = len(mod)
        if len(mod) > max_pts:
            mod = mod[np.sort(rng.choice(len(mod), max_pts, replace=False))]

        # amostragem do envelope para 5000 pontos
        env = env_full
        if sample_env and len(env) > sample_env:
            env = env[np.sort(rng.choice(len(env), sample_env, replace=False))]

        thr = max(3.0, grid_spacing(env))
        env_tree = cKDTree(env)
        center = mod.mean(axis=0)

        best = None
        for k in range(restarts):
            if k == 0:
                start, T_init = mod, np.eye(4)   # orientacao original
            else:
                Rm = Rotation.random(random_state=int(rng.integers(1 << 31))).as_matrix()
                start = (Rm @ (mod - center).T).T + center
                T_init = np.eye(4)
                T_init[:3, :3] = Rm
                T_init[:3, 3] = center - Rm @ center

            T_icp, cost = icp(start, env, env_tree, max_iter, penalty, thr,
                              mode=mode, w_leak=w_leak, w_cover=w_cover)
            if best is None or cost < best[0]:
                best = (cost, T_icp @ T_init)

        cost, T = best
        fit = apply_T(mod, T)

        return dict(
            entry=acc, status="ok",
            custo=round(cost, 4),
            chamfer=round(chamfer(fit, env), 4),
            hausdorff=round(hausdorff(fit, env), 4),
            dice=round(dice(fit, env), 5),
            frac_fora=round(frac_outside(fit, env_full), 4),
            cobertura=round(float((cKDTree(fit).query(env_full)[0] <= thr).mean()), 4),
            sep_centroide=round(float(np.linalg.norm(fit.mean(0) - env_full.mean(0))), 3),
            d_med=round(float(np.median(cKDTree(env_full).query(fit)[0])), 3),
            n_model=len(mod), n_model_full=n_model_full,
            n_env=len(env), grid=round(thr, 2),
            segundos=round(time.time() - t0, 2),
        )
    except Exception as e:
        return dict(entry=acc, status=f"erro: {type(e).__name__}: {e}")


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="pasta com data/sasbdb/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--penalty", type=float, default=50.0,
                    help="lambda; use 0 para o braco de controle")
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--max-iter", type=int, default=50)
    ap.add_argument("--max-points", type=int, default=300)
    ap.add_argument("--sample-env", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--mode", default="bidir",
                    choices=["uni", "bidir", "author"],
                    help="uni=so vazamento; bidir=limiar nos dois termos; "
                         "author=formulacao original (medias simples)")
    ap.add_argument("--w-leak", type=float, default=1.0)
    ap.add_argument("--w-cover", type=float, default=1.0)
    args = ap.parse_args()

    sasbdb = os.path.join(os.path.expanduser(args.base), "data", "sasbdb")
    if not os.path.isdir(sasbdb):
        print(f"nao encontrei {sasbdb}")
        return 1

    entries = sorted(os.path.basename(p).replace("_envelope.cif", "")
                     for p in glob.glob(os.path.join(sasbdb, "*_envelope.cif")))
    if args.limit:
        entries = entries[:args.limit]

    jobs = []
    for i, acc in enumerate(entries):
        env = os.path.join(sasbdb, f"{acc}_envelope.cif")
        mod = os.path.join(sasbdb, f"{acc}_model.pdb")
        if os.path.exists(mod):
            jobs.append((acc, env, mod, args.penalty, args.restarts,
                         args.max_iter, args.max_points, args.sample_env,
                         args.seed + i, args.mode,
                         args.w_leak, args.w_cover))

    workers = args.workers or min(cpu_count(), len(jobs))
    print(f"{len(jobs)} entradas | modo={args.mode} | lambda={args.penalty} | "
          f"restarts={args.restarts} | {workers} processos\n")

    rows = []
    with Pool(processes=workers) as pool:
        for i, r in enumerate(pool.imap_unordered(process, jobs), 1):
            rows.append(r)
            if r.get("status") == "ok":
                print(f"[{i}/{len(jobs)}] {r['entry']:9s} "
                      f"chamfer={r['chamfer']:7.3f}  fora={r['frac_fora']:.3f}  "
                      f"sep={r['sep_centroide']:6.2f}")
            else:
                print(f"[{i}/{len(jobs)}] {r['entry']:9s} {r.get('status')}")

    import pandas as pd
    df = pd.DataFrame(rows).sort_values("entry")
    df.to_csv(args.out, index=False)

    ok = df[df.status == "ok"]
    print(f"\n{len(ok)}/{len(df)} concluidas -> {args.out}")
    if len(ok):
        print(f"\nchamfer        media {ok.chamfer.mean():6.2f}  mediana {ok.chamfer.median():6.2f}")
        print(f"frac_fora      media {ok.frac_fora.mean():6.3f}  mediana {ok.frac_fora.median():6.3f}")
        if "cobertura" in ok:
            print(f"cobertura      media {ok.cobertura.mean():6.3f}  mediana {ok.cobertura.median():6.3f}")
        print(f"sep_centroide  media {ok.sep_centroide.mean():6.2f}  mediana {ok.sep_centroide.median():6.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
