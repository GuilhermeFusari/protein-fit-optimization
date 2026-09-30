#!/usr/bin/env python3
"""
ground_truth_v2.py  -  copia corrigida de scripts/ground_truth.py
(o original fica intacto). Experimento de recuperacao de pose conhecida.

Correcoes em relacao ao original (defeitos apontados no parecer):
  (a) erro angular so para rotacoes PROPRIAS. Se a transformacao total
      (ajuste x embaralhamento) tem det = -1, o metodo escolheu a imagem
      especular: erro_ang = NaN e espelhado = 1, e essas tentativas saem da
      media angular (continuam contando para RMSD e sucesso).
  (b) sementes reprodutiveis: a semente de cada entrada e
      seed + crc32(codigo) (o original usava hash(), que muda a cada
      execucao do Python). O DAMMIF tambem recebe --seed fixo, e os
      envelopes gerados podem ser guardados (--keep-envelopes) e reusados
      (--reuse-envelopes) para repetir o teste sem rodar o DAMMIF.
  (c) o cifsup (ATSAS) roda no MESMO teste, na mesma estrutura embaralhada
      e no mesmo envelope, nos metodos NSD e ICP (busca de enantiomeros
      ligada, o default do cifsup, como no nosso pipeline).

Nada mais foi mudado: mesma cadeia CRYSOL -> GNOM -> DAMMIF, mesmo
alinhamento do envelope por PCA para definir a "pose verdadeira", mesmo
pipeline proprio (fit_pipeline: author, lambda 0.2, 3 inicios, enantiomeros),
mesmo criterio de sucesso (RMSD < 5 A).

Opcao extra (nao muda o metodo): --restarts N roda o ICP-SAXS com N inicios
por mao nos MESMOS envelopes (--reuse-envelopes) e embaralhamentos, para
separar falha da busca de falha do envelope.

Colunas novas: metodo (icp_saxs | cifsup_nsd | cifsup_icp), espelhado,
pca_rmsd_propria / pca_rmsd_impropria (diagnostico: se o alinhamento por
PCA com reflexao encaixa melhor, o envelope do DAMMIF saiu com a mao
trocada e a "pose verdadeira" por rotacao propria e so aproximada).

Uso:
    python3 ground_truth_v2.py --base <pasta com data/sasbdb/> \
        --entries SASDPT5 SASDTS8 ... --trials 4 --out gt_v2.csv \
        --keep-envelopes envelopes_v2

Requer crysol, gnom, dammif e cifsup no PATH (ATSAS 3.2+).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time
import zlib

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


# leitura
def read_ca(path):
    pts = []
    with open(path, "r", errors="replace") as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "CA":
                try:
                    pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    pass
    if not pts:  # fallback: todos os atomos
        with open(path, "r", errors="replace") as f:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    try:
                        pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                    except ValueError:
                        pass
    return np.asarray(pts, dtype=float)


def read_pdb_any(path):
    pts = []
    with open(path, "r", errors="replace") as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                try:
                    pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    pass
    return np.asarray(pts, dtype=float)


def write_ca_pdb(pts, path):
    with open(path, "w") as f:
        for i, (x, y, z) in enumerate(pts, 1):
            f.write(f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("END\n")


def write_dummy_cif(pts, path):
    """Grava pontos como mmCIF de dummy atoms, formato que o Chimera
    renderiza como superficie de envelope."""
    with open(path, "w") as f:
        f.write("data_envelope\n#\nloop_\n")
        f.write("_atom_site.group_PDB\n_atom_site.id\n")
        f.write("_atom_site.type_symbol\n_atom_site.label_atom_id\n")
        f.write("_atom_site.label_comp_id\n_atom_site.label_asym_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n")
        for i, (x, y, z) in enumerate(pts, 1):
            f.write(f"ATOM {i} O DUM DUM A {i} {x:.3f} {y:.3f} {z:.3f} 1.00 0.00\n")
        f.write("#\n")


# geometria
def principal_axes(pts):
    c = pts.mean(axis=0)
    u, s, vt = np.linalg.svd(pts - c, full_matrices=False)
    R = vt                       # linhas = eixos principais
    if np.linalg.det(R) < 0:
        R[-1] *= -1
    return c, R


def align_by_pca(src, ref):
    """
    Traz src ao referencial de ref por eixos principais, testando as 4
    combinacoes de sinal com determinante +1; retorna a de menor RMSD.
    """
    cs, Rs = principal_axes(src)
    cr, Rr = principal_axes(ref)
    src0 = (src - cs)
    best = None
    signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]  # det +1
    for sg in signs:
        D = np.diag(sg)
        R = Rr.T @ D @ Rs
        if np.linalg.det(R) < 0:
            continue
        moved = (R @ src0.T).T + cr
        r = rmsd_min_count(moved, ref)
        if best is None or r < best[0]:
            best = (r, moved)
    return best[1]


def rmsd_min_count(a, b):
    """RMSD entre conjuntos de tamanhos possivelmente diferentes: casa cada
    ponto de a ao mais proximo de b (proxy quando nao ha correspondencia 1:1)."""
    d, _ = cKDTree(b).query(a, k=1)
    return float(np.sqrt((d ** 2).mean()))


def rmsd_paired(a, b):
    """RMSD com correspondencia 1:1 (mesma ordem de atomos)."""
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def chamfer_env(pts, env):
    """Chamfer bidirecional entre a estrutura e o envelope."""
    return float((cKDTree(env).query(pts)[0].mean() +
                  cKDTree(pts).query(env)[0].mean()) / 2.0)


def anisometry(pts):
    """Razao entre o maior e o menor eixo principal. ~1 = globular,
    grande = alongada. Indica quanto o envelope define a orientacao."""
    c = pts - pts.mean(axis=0)
    _, s, _ = np.linalg.svd(c, full_matrices=False)
    return float(s[0] / s[-1]) if s[-1] > 1e-9 else float("inf")


def angle_proper(R):
    """(a) Angulo (graus) de uma rotacao PROPRIA; NaN se det(R) < 0
    (reflexao: a formula do traco nao da um angulo de rotacao valido)."""
    if np.linalg.det(R) < 0:
        return float("nan")
    cos = (np.trace(R) - 1) / 2
    return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))


def fit_rigid(a, b):
    """Transformacao ortogonal (reflexao permitida) que leva a em b,
    pontos pareados. Usada para recuperar a matriz aplicada pelo cifsup."""
    ca, cb = a.mean(0), b.mean(0)
    U, _, Vt = np.linalg.svd((a - ca).T @ (b - cb))
    R = Vt.T @ U.T
    return R, cb - R @ ca


def pca_hand_diagnostic(src, ref):
    """RMSD (sem pareamento) do melhor alinhamento por PCA com det +1 e
    com det -1."""
    cs, Rs = principal_axes(src)
    cr, Rr = principal_axes(ref)
    out = {}
    for name, dets in (("propria", 1), ("impropria", -1)):
        best = None
        for sg in [(a, b, c) for a in (1, -1) for b in (1, -1) for c in (1, -1)]:
            R = Rr.T @ np.diag(sg) @ Rs
            if np.sign(np.linalg.det(R)) != dets:
                continue
            r = rmsd_min_count((R @ (src - cs).T).T + cr, ref)
            best = r if best is None else min(best, r)
        out[name] = best
    return out


def angle_between(R1, R2):
    """Angulo (graus) entre duas rotacoes."""
    Rrel = R1 @ R2.T
    cos = (np.trace(Rrel) - 1) / 2
    return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))


# ATSAS
def run(cmd, cwd, timeout=1800):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return -1, "timeout"
    except FileNotFoundError:
        print(f"ERRO: '{cmd[0]}' nao encontrado no PATH", file=sys.stderr)
        sys.exit(1)


def _crysol_curve(model_path, workdir):
    """Roda CRYSOL; retorna (caminho .int, diametro do envelope em A) ou (None, None)."""
    rc, log = run(["crysol", model_path], workdir)
    ints = glob.glob(os.path.join(workdir, "*.int"))
    if not ints:
        return None, None
    # extrai o diametro do envelope do log (linha "Envelope Diameter")
    diam = None
    for line in log.splitlines():
        if "Envelope Diameter" in line:
            try:
                diam = float(line.split(":")[-1])
            except ValueError:
                pass
    return ints[0], diam


def _convert_curve(int_path, out_dat):
    """Converte o .int do CRYSOL para .dat de 3 colunas que o GNOM aceita:
    q>0, intensidade total (coluna 2) normalizada, erro sintetico 1%."""
    q, I = [], []
    for l in open(int_path, errors="replace"):
        p = l.split()
        if len(p) >= 2:
            try:
                qi, Ii = float(p[0]), float(p[1])
                if qi > 1e-6:
                    q.append(qi); I.append(Ii)
            except ValueError:
                continue
    if not q:
        return False
    q = np.array(q); I = np.array(I)
    I = I / I[0] * 100.0
    err = 0.01 * I
    with open(out_dat, "w") as f:
        for a, b, c in zip(q, I, err):
            f.write(f"{a:.6f}\t{b:.6f}\t{c:.6f}\n")
    return True


def make_envelope(model_path, workdir, dammif_seed=None):
    """
    Cadeia CRYSOL -> GNOM -> DAMMIF (ATSAS 3.2).
    Retorna (pts do envelope, status). O DAMMIF 3.2 grava o modelo em .cif.
    """
    int_path, diam = _crysol_curve(model_path, workdir)
    if not int_path:
        return None, "crysol falhou", None
    if not diam or diam <= 0:
        diam = 100.0                       # fallback conservador

    dat = os.path.join(workdir, "curve.dat")
    if not _convert_curve(int_path, dat):
        return None, "conversao da curva falhou", None

    gout = os.path.join(workdir, "gnom.out")
    rc, log = run(["gnom", dat, f"--rmax={diam:.1f}", f"--output={gout}"], workdir)
    if not os.path.exists(gout):
        return None, f"gnom falhou (rc={rc})", None

    cmd = ["dammif", "--mode=fast", "--prefix=dam", gout]
    if dammif_seed is not None:
        cmd.insert(1, f"--seed={dammif_seed}")
    rc, log = run(cmd, workdir)
    envs = glob.glob(os.path.join(workdir, "dam*.cif")) + \
           glob.glob(os.path.join(workdir, "dam*.pdb"))
    if not envs:
        return None, f"dammif falhou (rc={rc})", None
    pts = read_pdb_any(envs[0])
    if len(pts) < 10:                      # .cif: le pelo parser de conteudo
        pts = _read_cif_atoms(envs[0])
    if len(pts) < 10:
        return None, "envelope vazio", None
    return pts, "ok", envs[0]              # devolve tambem o caminho do .cif


def _read_cif_atoms(path):
    """Le coordenadas de um mmCIF (dummy atoms do DAMMIF)."""
    cols, rows, in_loop = {}, [], False
    for l in open(path, errors="replace"):
        s = l.strip()
        if s.startswith("_atom_site."):
            cols[s.split(".", 1)[1].split()[0]] = len(cols); in_loop = True; continue
        if in_loop:
            if s.startswith("#") or s.startswith("loop_") or not s:
                if rows: break
                continue
            if s.startswith("_"): break
            rows.append(s.split())
    need = ("Cartn_x", "Cartn_y", "Cartn_z")
    if not all(k in cols for k in need):
        return np.empty((0, 3))
    ix, iy, iz = (cols[k] for k in need)
    pts = []
    for p in rows:
        if len(p) > max(ix, iy, iz):
            try: pts.append((float(p[ix]), float(p[iy]), float(p[iz])))
            except ValueError: pass
    return np.asarray(pts, dtype=float)


# pipeline proprio (author, lambda 0.2, com enantiomorfos)
def author_cost(src, env_pts, env_tree):
    d_fill, _ = cKDTree(src).query(env_pts, k=1)
    d_leak, _ = env_tree.query(src, k=1)
    return float(d_fill.mean()), float(d_leak.mean())


def icp_author(src_pts, env_pts, env_tree, penalty=0.2, max_iter=50):
    src = src_pts.copy()
    T = np.eye(4)
    best_cost, best_T = float("inf"), np.eye(4)
    for _ in range(max_iter):
        dist, idx = env_tree.query(src, k=1)
        closest = env_pts[idx]
        fill, leak = author_cost(src, env_pts, env_tree)
        cost = fill + penalty * leak
        if cost < best_cost:
            best_cost, best_T = cost, T.copy()
        w = 1.0 + penalty * dist / max(float(dist.mean()), 1e-9)
        st = cKDTree(src)
        d_env, i_env = st.query(env_pts, k=1)
        srcA = np.vstack([src, src[i_env]])
        tgtA = np.vstack([closest, env_pts])
        wA = np.concatenate([w, d_env / max(float(d_env.mean()), 1e-9)])
        wA = wA / wA.sum()
        sm = (wA[:, None] * srcA).sum(0)
        tm = (wA[:, None] * tgtA).sum(0)
        H = (srcA - sm).T @ (wA[:, None] * (tgtA - tm))
        U, _, Vt = np.linalg.svd(H)
        Rm = Vt.T @ U.T
        if np.linalg.det(Rm) < 0:
            Vt[-1] *= -1; Rm = Vt.T @ U.T
        t = tm - Rm @ sm
        src = (Rm @ src.T).T + t
        Tn = np.eye(4); Tn[:3, :3] = Rm; Tn[:3, 3] = t
        T = Tn @ T
    return best_T


def fit_pipeline(mob, env, restarts=3, seed=0):
    """Alinha mob (embaralhado) ao env; retorna a transformacao 4x4 aplicada."""
    env_tree = cKDTree(env)
    c = mob.mean(0)
    rng = np.random.default_rng(seed)
    best = None
    mirrors = [np.eye(3), np.diag([1, 1, -1.0])]
    for M in mirrors:
        base = (M @ (mob - c).T).T + c if not np.allclose(M, np.eye(3)) else mob
        Tm = np.eye(4)
        if not np.allclose(M, np.eye(3)):
            Tm[:3, :3] = M; Tm[:3, 3] = c - M @ c
        for k in range(restarts):
            if k == 0:
                start, Ti = base, Tm
            else:
                Rr = Rotation.random(random_state=int(rng.integers(1 << 31))).as_matrix()
                start = (Rr @ (base - c).T).T + c
                Tr = np.eye(4); Tr[:3, :3] = Rr; Tr[:3, 3] = c - Rr @ c
                Ti = Tr @ Tm
            Ticp = icp_author(start, env, env_tree)
            T = Ticp @ Ti
            moved = (T[:3, :3] @ mob.T).T + T[:3, 3]
            _, leak = author_cost(moved, env, env_tree)
            fill, _ = author_cost(moved, env, env_tree)
            cost = fill + 0.2 * leak
            if best is None or cost < best[0]:
                best = (cost, T)
    return best[1]


def run_cifsup(method, template, movable, out_pdb, timeout=600):
    """(c) cifsup do ATSAS; devolve (status, segundos)."""
    t0 = time.time()
    try:
        p = subprocess.run(["cifsup", f"--method={method}", "-o", out_pdb,
                            template, movable],
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "timeout", time.time() - t0
    except FileNotFoundError:
        print("ERRO: 'cifsup' nao encontrado no PATH", file=sys.stderr)
        sys.exit(1)
    ok = p.returncode == 0 and os.path.exists(out_pdb)
    return ("ok" if ok else f"exit={p.returncode}"), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--entries", nargs="+", required=True)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--success-A", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="ground_truth_v2.csv")
    ap.add_argument("--keep-envelopes", default="",
                    help="pasta para guardar os envelopes sinteticos (ja alinhados)")
    ap.add_argument("--reuse-envelopes", default="",
                    help="pasta com envelopes guardados: nao roda CRYSOL/DAMMIF")
    ap.add_argument("--no-cifsup", action="store_true")
    ap.add_argument("--restarts", type=int, default=3,
                    help="inicios por mao do ICP-SAXS (3 = publicado); os "
                         "embaralhamentos nao dependem deste valor")
    args = ap.parse_args()

    sasbdb = os.path.join(os.path.expanduser(args.base), "data", "sasbdb")
    rows = []

    for acc in args.entries:
        model = os.path.abspath(os.path.join(sasbdb, f"{acc}_model.pdb"))
        if not os.path.exists(model):
            print(f"{acc}: modelo nao encontrado"); continue

        ref = read_ca(model)
        if len(ref) < 4:
            print(f"{acc}: poucas coordenadas"); continue

        # (b) semente reprodutivel por entrada
        entry_seed = args.seed + zlib.crc32(acc.encode())
        work = tempfile.mkdtemp(prefix=f"gt2_{acc}_")
        try:
            t0 = time.time()
            saved = os.path.join(args.reuse_envelopes, f"{acc}_env.pdb") \
                if args.reuse_envelopes else ""
            if saved and os.path.exists(saved):
                env = read_pdb_any(saved)
                diag = pca_hand_diagnostic(env, ref)
                print(f"\n[{acc}] envelope reusado: {len(env)} pts")
            else:
                print(f"\n[{acc}] gerando envelope sintetico (crysol + gnom + dammif)...")
                env_raw, st, _ = make_envelope(model, work,
                                               dammif_seed=entry_seed % (2**31 - 1))
                if env_raw is None:
                    print(f"  {st}")
                    rows.append(dict(entry=acc, status=st)); continue
                diag = pca_hand_diagnostic(env_raw, ref)
                env = align_by_pca(env_raw, ref)
                print(f"  envelope: {len(env)} pts | {time.time()-t0:.0f}s")
            if args.keep_envelopes:
                os.makedirs(args.keep_envelopes, exist_ok=True)
                write_ca_pdb(env, os.path.join(args.keep_envelopes, f"{acc}_env.pdb"))

            env_pdb = os.path.join(work, "envelope.pdb")
            write_ca_pdb(env, env_pdb)
            cham_true = chamfer_env(ref, env)
            aniso = anisometry(ref)

            rng = np.random.default_rng(entry_seed)
            for tr in range(args.trials):
                Rtrue = Rotation.from_quat(rng.normal(size=4)).as_matrix()
                ttrue = rng.normal(0, 15, 3)
                c = ref.mean(0)
                mob = (Rtrue @ (ref - c).T).T + c + ttrue
                fit_seed = int(rng.integers(1 << 31))

                results = {}
                t_fit = time.time()
                T = fit_pipeline(mob, env, restarts=args.restarts, seed=fit_seed)
                name = "icp_saxs" if args.restarts == 3 else f"icp_saxs_r{args.restarts}"
                results[name] = ("ok", (T[:3, :3] @ mob.T).T + T[:3, 3],
                                 time.time() - t_fit)

                if not args.no_cifsup:
                    mob_pdb = os.path.join(work, f"mob_{tr}.pdb")
                    write_ca_pdb(mob, mob_pdb)
                    for method in ("NSD", "ICP"):
                        out = os.path.join(work, f"cs_{method}_{tr}.pdb")
                        st, dt = run_cifsup(method, env_pdb, mob_pdb, out)
                        rec = read_ca(out) if st == "ok" else None
                        if rec is not None and len(rec) != len(mob):
                            st, rec = f"atomos {len(rec)} != {len(mob)}", None
                        results[f"cifsup_{method.lower()}"] = (st, rec, dt)

                for name, (st, rec, dt) in results.items():
                    row = dict(entry=acc, metodo=name, trial=tr, status=st,
                               anisometria=round(aniso, 2), n_env=len(env),
                               n_model=len(ref),
                               pca_rmsd_propria=round(diag["propria"], 3),
                               pca_rmsd_impropria=round(diag["impropria"], 3),
                               chamfer_verdadeira=round(cham_true, 3),
                               segundos=round(dt, 2))
                    if rec is not None:
                        R, _ = fit_rigid(mob, rec)
                        Rtot = R @ Rtrue
                        r = rmsd_paired(rec, ref)
                        cham_rec = chamfer_env(rec, env)
                        row.update(rmsd=round(r, 3),
                                   erro_ang=round(angle_proper(Rtot), 1),
                                   espelhado=int(np.linalg.det(Rtot) < 0),
                                   chamfer_recuperada=round(cham_rec, 3),
                                   delta_chamfer=round(cham_rec - cham_true, 3),
                                   sucesso=int(r < args.success_A))
                    rows.append(row)
                    if rec is not None:
                        print(f"  trial {tr} {name:11s}: RMSD={row['rmsd']:6.2f} A  "
                              f"ang={row['erro_ang']:5.1f}deg  "
                              f"{'espelhado ' if row['espelhado'] else ''}"
                              f"{'OK' if row['sucesso'] else 'falhou'}")
                    else:
                        print(f"  trial {tr} {name:11s}: {st}")
        finally:
            shutil.rmtree(work, ignore_errors=True)

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df[df.status == "ok"]
    print(f"\n{'='*60}")
    for name, g in ok.groupby("metodo"):
        prop = g[g.espelhado == 0]
        print(f"{name:11s} n={len(g):3d}  sucesso={g.sucesso.mean()*100:5.1f}%  "
              f"RMSD medio={g.rmsd.mean():6.2f} A  espelhado={int(g.espelhado.sum())}  "
              f"erro ang (so proprias, n={len(prop)})={prop.erro_ang.mean():5.1f} deg  "
              f"delta_cham={g.delta_chamfer.mean():+.3f}")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
