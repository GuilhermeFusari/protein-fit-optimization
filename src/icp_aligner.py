#!/usr/bin/env python3
"""
Baseline ICP alignment (pure Chamfer distance, no penalties).
Used for evaluating single structure fits against the envelope.
"""
import os
import time
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
from multiprocessing import Pool

def extract_structure_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PROT", pdb_path)
    coords = np.array([atom.coord for atom in structure.get_atoms()], dtype=float)
    return coords, structure

def coords_from_cif_dict(cif_dict):
    xs = np.array([float(x) for x in cif_dict['_atom_site.Cartn_x']], dtype=float)
    ys = np.array([float(y) for y in cif_dict['_atom_site.Cartn_y']], dtype=float)
    zs = np.array([float(z) for z in cif_dict['_atom_site.Cartn_z']], dtype=float)
    return np.stack([xs, ys, zs], axis=1)

def chamfer_distance(setA, setB):
    treeA = cKDTree(setA)
    treeB = cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1)
    dists_B_to_A, _ = treeA.query(setB, k=1)
    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0

def icp_align_with_chamfer(src_pts, tgt_pts, max_iter=50, tol=1e-6):
    src, tgt = src_pts.copy(), tgt_pts.copy()
    final_T = np.eye(4)
    best_score = float('inf')
    tree_tgt = cKDTree(tgt)

    for _ in range(max_iter):
        _, indices = tree_tgt.query(src, k=1)
        closest = tgt[indices]
        
        current_score = chamfer_distance(src, tgt)
        best_score = min(best_score, current_score)

        # Procrustes analysis
        src_mean, tgt_mean = src.mean(axis=0), closest.mean(axis=0)
        H = (src - src_mean).T @ (closest - tgt_mean)
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Handle reflection
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
            
        t = tgt_mean - R @ src_mean
        src = (R @ src.T).T + t

        T = np.eye(4)
        T[:3, :3], T[:3, 3] = R, t
        final_T = T @ final_T

        if current_score < tol:
            break

    return final_T, chamfer_distance(src, tgt)

def apply_transform_to_structure(structure, T4):
    R, t = T4[:3, :3], T4[:3, 3]
    for atom in structure.get_atoms():
        atom.set_coord(R @ atom.coord + t)

def evaluate_pdb_alignment(args_tuple):
    pdb_path, env_coords, align_what, max_iter, sample_env = args_tuple
    try:
        prot_coords, _ = extract_structure_from_pdb(pdb_path)
        
        # Downsample envelope for speed if requested
        env_for_icp = env_coords
        if sample_env and sample_env > 0 and env_coords.shape[0] > sample_env:
            idx = np.random.choice(env_coords.shape[0], size=sample_env, replace=False)
            env_for_icp = env_coords[idx]

        if align_what == 'protein':
            _, score = icp_align_with_chamfer(prot_coords, env_for_icp, max_iter=max_iter)
        else:
            _, score = icp_align_with_chamfer(env_for_icp, prot_coords, max_iter=max_iter)

        return True, score, os.path.basename(pdb_path), pdb_path
    except Exception as e:
        return False, float('inf'), f"{os.path.basename(pdb_path)}: {str(e)}", pdb_path

def find_best_pdb(args):
    os.makedirs(args.output, exist_ok=True)
    start_time = time.time()

    cif_dict = MMCIF2Dict(open(args.envelope, 'r', encoding='utf-8', errors='ignore'))
    env_coords = coords_from_cif_dict(cif_dict)
    
    input_path = Path(args.input)
    pdb_files = [input_path] if input_path.is_file() else list(input_path.glob("*.pdb"))
    
    if not pdb_files:
        print("No PDBs found.")
        return

    work_args = [(str(p), env_coords, args.align_what, args.max_iter, args.sample_env) for p in pdb_files]
    
    best_score, best_pdb_path, best_pdb_name = float('inf'), None, None
    processed_count = 0

    with Pool(processes=args.workers) as pool:
        for success, score, pdb_name, pdb_path in pool.imap_unordered(evaluate_pdb_alignment, work_args):
            if success:
                processed_count += 1
                if score < best_score:
                    best_score, best_pdb_path, best_pdb_name = score, pdb_path, pdb_name

    # Re-apply transformation to the best candidate and save
    if best_pdb_path:
        prot_coords, prot_struct = extract_structure_from_pdb(best_pdb_path)
        env_for_icp = env_coords
        if args.sample_env and args.sample_env > 0 and env_coords.shape[0] > args.sample_env:
            idx = np.random.choice(env_coords.shape[0], size=args.sample_env, replace=False)
            env_for_icp = env_coords[idx]

        T, _ = icp_align_with_chamfer(prot_coords, env_for_icp, max_iter=args.max_iter)
        apply_transform_to_structure(prot_struct, T)

        out_path = os.path.join(args.output, f"BEST_ALIGNMENT_{best_pdb_name}")
        io = PDBIO()
        io.set_structure(prot_struct)
        io.save(out_path)
    
    elapsed_time = time.time() - start_time

    # Generate standardized metrics report for data extraction
    report_path = os.path.join(args.output, "report.txt")
    with open(report_path, "w") as f:
        f.write("ALIGNMENT REPORT (SINGLE MODE)\n==============================\n\n")
        if best_pdb_path:
            f.write(f"Top 1 Error (Individual): {best_score:.6f} A\n")
            f.write(f"Best PDB File:            {best_pdb_name}\n")
        else:
            f.write("Status: FAILED\n")
            
        f.write(f"Execution Time:           {elapsed_time:.2f} seconds\n")
        f.write(f"Envelope Size (points):   {env_coords.shape[0]}\n")
        f.write(f"PDBs Evaluated:           {processed_count}\n")
    
    print(f"Report saved: {report_path} | Best Score: {best_score:.4f} A")