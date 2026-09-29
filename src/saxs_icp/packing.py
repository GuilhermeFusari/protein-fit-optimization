"""
Multi-structure packing: ranks a folder of conformations by penalized ICP
against one envelope and writes the top copies aligned.

Alignment uses the bidirectional score (fill + lambda * leak); the reported
score is the plain Chamfer distance, for objective reporting.
"""
import os
import time
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as R_scipy

from .core import random_rotation
from .io import read_coords


def _parser():
    from Bio.PDB import PDBParser
    return PDBParser(QUIET=True)


def extract_coords(pdb_path):
    struct = _parser().get_structure("X", pdb_path)
    pts = np.array([atom.coord for atom in struct.get_atoms()], dtype=float)
    if len(pts) < 4:
        raise ValueError(f"{len(pts)} atoms read (need at least 4)")
    return pts


def weighted_asymmetric_score(src, env_tree, penalty, env_pts=None):
    """
    Bidirectional cost (published formulation):
        score = mean(dist envelope->protein) + penalty * mean(dist protein->envelope)
    """
    dists_leak, _ = env_tree.query(src, k=1)
    score_leaking = np.mean(dists_leak)

    if env_pts is None:
        env_pts = env_tree.data
    prot_tree = cKDTree(src)
    dists_fill, _ = prot_tree.query(env_pts, k=1)
    score_filling = np.mean(dists_fill)

    return score_filling + penalty * score_leaking


def icp_align_with_penalty(src_pts, env_pts, max_iter=30, penalty=0.2):
    src = src_pts.copy()
    env_tree = cKDTree(env_pts)
    final_T = np.eye(4)
    best_score = float('inf')

    for _ in range(max_iter):
        _, indices = env_tree.query(src, k=1)
        closest = env_pts[indices]

        score = weighted_asymmetric_score(src, env_tree, penalty=penalty, env_pts=env_pts)
        best_score = min(best_score, score)

        src_mean, tgt_mean = src.mean(axis=0), closest.mean(axis=0)
        H = (src - src_mean).T @ (closest - tgt_mean)
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        t = tgt_mean - R @ src_mean
        src = (R @ src.T).T + t

        T = np.eye(4)
        T[:3, :3], T[:3, 3] = R, t
        final_T = T @ final_T

    return final_T, best_score


def evaluate_pdb_packing(args):
    pdb_path, env_coords, max_iter, sample_size, penalty, seed = args
    try:
        pts = extract_coords(pdb_path)
        rs = np.random if seed is None else np.random.RandomState(seed)

        if sample_size and len(env_coords) > sample_size:
            idx = rs.choice(len(env_coords), sample_size, replace=False)
            env_sampled = env_coords[idx]
        else:
            env_sampled = env_coords

        # random initial rotation
        center = pts.mean(axis=0)
        if seed is None:
            rot_matrix = R_scipy.random().as_matrix()
        else:
            rot_matrix = random_rotation(int(rs.randint(1 << 31)))
        pts_rotated = (rot_matrix @ (pts - center).T).T + center

        init_T = np.eye(4)
        init_T[:3, :3] = rot_matrix
        init_T[:3, 3] = center - rot_matrix @ center

        T_iter, score = icp_align_with_penalty(pts_rotated, env_sampled, max_iter, penalty)

        final_T = T_iter @ init_T

        return True, score, pdb_path, final_T, ""
    except Exception as e:
        return False, float('inf'), pdb_path, np.eye(4), f"{type(e).__name__}: {e}"


def get_best_candidates(input_folder, env_coords, top_n=20, penalty=0.2, max_iter=30,
                        sample_env=2000, workers=0, seed=None, log=print):
    files = sorted(os.path.join(input_folder, f) for f in os.listdir(input_folder)
                   if f.lower().endswith('.pdb'))
    if not files:
        raise FileNotFoundError(f"no .pdb files in {input_folder}")

    work_args = [(p, env_coords, max_iter, sample_env, penalty,
                  None if seed is None else seed + i) for i, p in enumerate(files)]
    results, failed = [], []

    workers = max(1, min(workers or cpu_count(), len(files)))
    with Pool(processes=workers) as pool:
        for success, score, p_path, T, err in pool.imap_unordered(evaluate_pdb_packing,
                                                                  work_args):
            if success:
                results.append((score, p_path, T))
            else:
                failed.append((p_path, err))

    if failed:
        log(f"Warning: {len(failed)} of {len(files)} files could not be read:")
        for p, err in failed[:10]:
            log(f"  {Path(p).name}: {err}")
        if len(failed) > 10:
            log(f"  ... and {len(failed) - 10} more")

    results.sort(key=lambda x: (x[0], x[1]))
    return results[:top_n]


def calculate_chamfer(setA, setB):
    if len(setA) == 0 or len(setB) == 0:
        return float('inf')

    treeA, treeB = cKDTree(setA), cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1)
    dists_B_to_A, _ = treeA.query(setB, k=1)

    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0


def run_packing(input_folder, envelope_path, output_folder, copies=20, penalty=0.2,
                max_iter=30, sample_env=2000, workers=0, seed=None, log=print):
    from Bio.PDB import PDBIO, Structure

    start_time = time.time()
    if not os.path.isdir(input_folder):
        raise NotADirectoryError(f"input folder not found: {input_folder}")
    env_coords = read_coords(envelope_path)
    if len(env_coords) < 10:
        raise ValueError(f"envelope has {len(env_coords)} points (need at least 10)")
    os.makedirs(output_folder, exist_ok=True)

    log(f"Packing started (top {copies}, penalty {penalty})")

    candidates = get_best_candidates(input_folder, env_coords, top_n=copies,
                                     penalty=penalty, max_iter=max_iter,
                                     sample_env=sample_env, workers=workers,
                                     seed=seed, log=log)

    if not candidates:
        with open(os.path.join(output_folder, "report.txt"), "w", encoding="utf-8") as f:
            f.write("ALIGNMENT REPORT\nStatus: FAILED (No candidates)\n")
        raise RuntimeError("no structure could be aligned")

    final_structs, summary_data, all_atoms_combined = [], [], []
    real_score_top1 = None

    for rank, (penalized_score, path, T) in enumerate(candidates):
        try:
            filename = Path(path).name
            struct = _parser().get_structure("X", path)
            atoms = list(struct.get_atoms())

            R, t = T[:3, :3], T[:3, 3]
            struct_coords = []
            for atom in atoms:
                new_coord = np.dot(R, atom.coord) + t
                atom.set_coord(new_coord)
                struct_coords.append(new_coord)

            all_atoms_combined.extend(struct_coords)

            current_atoms_np = np.array(struct_coords)
            real_dist = calculate_chamfer(current_atoms_np, env_coords)

            if rank == 0:
                real_score_top1 = real_dist

            indiv_name = f"RANK_{rank+1:02d}_{filename}"
            io = PDBIO()
            io.set_structure(struct)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                io.save(os.path.join(output_folder, indiv_name))

            final_structs.append(struct)
            summary_data.append({"rank": rank+1, "name": filename, "score": real_dist})

        except Exception as e:
            log(f"Error processing {path}: {e}")

    all_atoms_np = np.array(all_atoms_combined)
    global_score = calculate_chamfer(all_atoms_np, env_coords)

    if final_structs:
        try:
            io = PDBIO()
            c = Structure.Structure("Combined")
            for i, m in enumerate(final_structs):
                model = next(iter(m))
                model.detach_parent()
                model.id = i
                model.serial_num = i + 1
                c.add(model)
            io.set_structure(c)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                io.save(os.path.join(output_folder, "ALL_TOP_ALIGNED.pdb"))
        except Exception as e:
            log(f"Warning: could not save combined file: {e}")

    real_score_top1 = real_score_top1 or float('inf')
    improvement = "YES" if global_score < real_score_top1 else "NO"
    elapsed_time = time.time() - start_time

    report_path = os.path.join(output_folder, "report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("ALIGNMENT REPORT\n================\n\n")
        f.write(f"Top 1 Error (Individual): {real_score_top1:.6f} A\n")
        f.write(f"Global Error (Combined):  {global_score:.6f} A\n")
        f.write(f"Improvement:              {improvement}\n")
        f.write(f"Number of Copies:         {len(final_structs)}\n")
        f.write(f"Penalty Weight Used:      {penalty}\n")
        f.write(f"Execution Time:           {elapsed_time:.2f} seconds\n")
        f.write(f"Envelope Size (points):   {len(env_coords)}\n\n")

        f.write("-" * 50 + "\n")
        f.write(f"{'RANK':<5} | {'CHAMFER (A)':<12} | {'FILE'}\n")
        for item in summary_data:
            f.write(f"#{item['rank']:<4} | {item['score']:.4f}       | {item['name']}\n")

    log(f"Report saved: {report_path} | Global: {global_score:.4f} A | "
        f"Time: {elapsed_time:.2f}s")
    return dict(top1=real_score_top1, global_score=global_score,
                copies=len(final_structs), report=report_path)
