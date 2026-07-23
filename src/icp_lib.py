#!/usr/bin/env python3
"""
Modified ICP core library.
Implements asymmetric penalization for volumetric boundaries,
with Random Initial Rotation for robust statistical sampling.
"""
import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as R_scipy
from multiprocessing import Pool, cpu_count
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def extract_coords(pdb_path):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("X", pdb_path)
    return np.array([atom.coord for atom in struct.get_atoms()], dtype=float)

def get_env_coords(cif_path):
    d = MMCIF2Dict(open(cif_path, 'r', encoding='utf-8', errors='ignore'))
    xs = np.array(d['_atom_site.Cartn_x'], dtype=float)
    ys = np.array(d['_atom_site.Cartn_y'], dtype=float)
    zs = np.array(d['_atom_site.Cartn_z'], dtype=float)
    return np.stack([xs, ys, zs], axis=1)

def weighted_asymmetric_score(src, env_tree, penalty):
    dists, _ = env_tree.query(src, k=1)
    base_dist = np.mean(dists)
    
    threshold = 3.0
    outliers = dists[dists > threshold]
    
    if len(outliers) > 0:
        penalty_term = penalty * np.sum(outliers - threshold) / len(src)
        return base_dist + penalty_term
    return base_dist

def icp_align_with_penalty(src_pts, env_pts, max_iter=30, penalty=50.0):
    src = src_pts.copy()
    env_tree = cKDTree(env_pts)
    final_T = np.eye(4)
    best_score = float('inf')

    for _ in range(max_iter):
        _, indices = env_tree.query(src, k=1)
        closest = env_pts[indices]
        
        score = weighted_asymmetric_score(src, env_tree, penalty=penalty)
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
    pdb_path, env_coords, max_iter, sample_size, penalty = args
    try:
        pts = extract_coords(pdb_path)
        
        if sample_size and len(env_coords) > sample_size:
            idx = np.random.choice(len(env_coords), sample_size, replace=False)
            env_sampled = env_coords[idx]
        else:
            env_sampled = env_coords

        # ========================================================
        # ALEATORIEDADE: Gira a proteína em um ângulo X, Y, Z aleatório
        # ========================================================
        center = pts.mean(axis=0)
        rot_matrix = R_scipy.random().as_matrix()
        pts_rotated = (rot_matrix @ (pts - center).T).T + center
        
        init_T = np.eye(4)
        init_T[:3, :3] = rot_matrix
        init_T[:3, 3] = center - rot_matrix @ center
        
        # Roda o ICP a partir desse ângulo novo
        T_iter, score = icp_align_with_penalty(pts_rotated, env_sampled, max_iter, penalty)
        
        # Junta o giro inicial com a rota do ICP
        final_T = T_iter @ init_T

        return True, score, pdb_path, final_T
    except Exception as e:
        return False, float('inf'), pdb_path, np.eye(4)

def get_best_candidates(input_folder, envelope_path, top_n=5, limit=None, penalty=50.0):
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.pdb')]
    if limit: files = files[:limit]
    
    env_coords = get_env_coords(envelope_path)
    
    work_args = [(p, env_coords, 30, 2000, penalty) for p in files]
    results = []

    with Pool(processes=cpu_count()) as pool:
        for success, score, p_path, T in pool.imap_unordered(evaluate_pdb_packing, work_args):
            if success:
                results.append((score, p_path, T))

    results.sort(key=lambda x: x[0])
    return results[:top_n], env_coords
