"""
Batch alignment over a SASBDB-style folder (<base>/data/sasbdb/<ACC>_envelope.cif
and <ACC>_model.pdb). Reproduces scripts/icp_saxs_v2.py: same per-entry seeds,
same CSV columns.
"""
import csv
import glob
import os
import time
from multiprocessing import Pool, cpu_count

from . import core
from .io import read_coords, write_points_pdb


def process(job):
    (acc, env_path, mod_path, penalty, restarts, max_iter,
     max_pts, sample_env, seed, mode, w_leak, w_cover, enantiomorphs, save_dir) = job
    t0 = time.time()
    try:
        env_full = read_coords(env_path)
        mod = read_coords(mod_path, ca_only=True)
        if len(env_full) < 10 or len(mod) < 4:
            return dict(entry=acc, status="coordenadas_vazias")

        r = core.align(mod, env_full, mode=mode, penalty=penalty, restarts=restarts,
                       max_iter=max_iter, max_points=max_pts, sample_env=sample_env,
                       seed=seed, enantiomorphs=enantiomorphs,
                       w_leak=w_leak, w_cover=w_cover)

        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            write_points_pdb(r.fitted, os.path.join(save_dir, f"{acc}_final.pdb"))

        m = r.metrics
        return dict(
            entry=acc, status="ok",
            custo=round(r.cost, 4),
            chamfer=round(m["chamfer"], 4),
            hausdorff=round(m["hausdorff"], 4),
            dice=round(m["dice"], 5),
            frac_fora=round(m["frac_outside"], 4),
            cobertura=round(m["coverage"], 4),
            sep_centroide=round(m["centroid_sep"], 3),
            d_med=round(m["median_dist"], 3),
            espelhado=bool(r.mirrored),
            n_model=r.n_model, n_model_full=r.n_model_full,
            n_env=r.n_env, grid=round(r.grid, 2),
            segundos=round(time.time() - t0, 2),
        )
    except Exception as e:
        return dict(entry=acc, status=f"erro: {type(e).__name__}: {e}")


def find_jobs(base, penalty=0.2, restarts=3, max_iter=50, max_points=3000,
              sample_env=5000, seed=42, mode="author", w_leak=1.0, w_cover=1.0,
              enantiomorphs=True, save_aligned="", limit=0):
    sasbdb = os.path.join(os.path.expanduser(base), "data", "sasbdb")
    if not os.path.isdir(sasbdb):
        raise FileNotFoundError(f"folder not found: {sasbdb}")

    entries = sorted(os.path.basename(p).replace("_envelope.cif", "")
                     for p in glob.glob(os.path.join(sasbdb, "*_envelope.cif")))
    if limit:
        entries = entries[:limit]

    jobs = []
    for i, acc in enumerate(entries):
        env = os.path.join(sasbdb, f"{acc}_envelope.cif")
        mod = os.path.join(sasbdb, f"{acc}_model.pdb")
        if os.path.exists(mod):
            jobs.append((acc, env, mod, penalty, restarts, max_iter, max_points,
                         sample_env, seed + i, mode, w_leak, w_cover,
                         enantiomorphs, save_aligned))
    return jobs


def write_csv(rows, path):
    rows = sorted(rows, key=lambda r: r["entry"])
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if isinstance(v, float) and v != v else v)
                        for k, v in r.items()})


def run(base, out, workers=0, log=print, **kw):
    jobs = find_jobs(base, **kw)
    if not jobs:
        raise FileNotFoundError(
            f"no <ACC>_envelope.cif + <ACC>_model.pdb pairs in {base}/data/sasbdb")

    workers = max(1, min(workers or cpu_count(), len(jobs)))
    log(f"{len(jobs)} entries | mode={kw.get('mode', 'author')} | "
        f"lambda={kw.get('penalty', 0.2)} | restarts={kw.get('restarts', 3)} | "
        f"{workers} processes\n")

    rows = []
    with Pool(processes=workers) as pool:
        for i, r in enumerate(pool.imap_unordered(process, jobs), 1):
            rows.append(r)
            if r.get("status") == "ok":
                log(f"[{i}/{len(jobs)}] {r['entry']:9s} "
                    f"chamfer={r['chamfer']:7.3f}  outside={r['frac_fora']:.3f}  "
                    f"sep={r['sep_centroide']:6.2f}")
            else:
                log(f"[{i}/{len(jobs)}] {r['entry']:9s} {r.get('status')}")

    write_csv(rows, out)

    ok = [r for r in rows if r.get("status") == "ok"]
    log(f"\n{len(ok)}/{len(rows)} done -> {out}")
    if ok:
        def mean(k):
            return sum(r[k] for r in ok) / len(ok)
        log(f"chamfer        mean {mean('chamfer'):6.3f}")
        log(f"frac_outside   mean {mean('frac_fora'):6.3f}")
        log(f"coverage       mean {mean('cobertura'):6.3f}")
        log(f"centroid_sep   mean {mean('sep_centroide'):6.3f}")
        log(f"mirrored       {sum(r['espelhado'] for r in ok)}/{len(ok)} entries")
    return rows
