#!/usr/bin/env python3
"""
icp_lib.py
Baseado 100% no código original do usuário.
Função: Calcular alinhamento ICP e retornar a Matriz de Transformação (T).
"""

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count

# --- SEU CÓDIGO ORIGINAL (Funções de Alinhamento) ---

def extract_structure_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("PROT", pdb_path)
    except:
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
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
    """
    Sua implementação original do ICP.
    """
    src = src_pts.copy()
    tgt = tgt_pts.copy()
    final_T = np.eye(4)
    best_score = float('inf')

    tree_tgt = cKDTree(tgt)

    for _ in range(max_iter):
        _, indices = tree_tgt.query(src, k=1)
        closest = tgt[indices]

        current_score = chamfer_distance(src, tgt)
        best_score = min(best_score, current_score)

        src_mean = src.mean(axis=0)
        tgt_mean = closest.mean(axis=0)
        src_c = src - src_mean
        tgt_c = closest - tgt_mean
        H = src_c.T @ tgt_c
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        t = tgt_mean - R @ src_mean
        
        # Atualiza coords
        src = (R @ src.T).T + t

        # Acumula Transformação
        T = np.eye(4)
        T[:3, :3] = R
        T[:3,  3] = t
        final_T = T @ final_T

        if current_score < tol:
            break

    final_score = chamfer_distance(src, tgt)
    return final_T, final_score

# --- WORKER MODIFICADO (Apenas para retornar a Matriz T) ---

def evaluate_pdb_alignment(args_tuple):
    pdb_path, env_coords, max_iter, sample_env = args_tuple
    try:
        prot_coords, _ = extract_structure_from_pdb(pdb_path)

        # Amostragem (igual ao seu original)
        if (sample_env is not None) and (sample_env > 0) and (env_coords.shape[0] > sample_env):
            idx = np.random.choice(env_coords.shape[0], size=sample_env, replace=False)
            env_for_icp = env_coords[idx]
        else:
            env_for_icp = env_coords

        # Roda o SEU ICP original
        # Note: passamos prot_coords cruas, igual você fazia
        T_matrix, score = icp_align_with_chamfer(prot_coords, env_for_icp, max_iter=max_iter)

        return True, score, pdb_path, T_matrix
    except Exception as e:
        return False, float('inf'), pdb_path, None

# --- ORQUESTRADOR PARA O PACKING ---

def get_best_candidates(input_folder, envelope_path, top_n=50, limit=None):
    # Carrega envelope
    cif = MMCIF2Dict(open(envelope_path, 'r', encoding='utf-8', errors='ignore'))
    env_coords = coords_from_cif_dict(cif)
    
    files = sorted(list(Path(input_folder).glob("*.pdb")))
    if limit: files = files[:limit]
    
    print(f"⚡ [ICP Original] Screening {len(files)} arquivos...")
    
    # Prepara args igual ao seu script
    work_args = [(str(p), env_coords, 30, 2000) for p in files]
    
    results = []
    processed = 0
    
    with Pool(processes=cpu_count()) as pool:
        for res in pool.imap_unordered(evaluate_pdb_alignment, work_args, chunksize=10):
            success, score, path, T = res
            if success and score < 1000:
                results.append((score, path, T))
            
            processed += 1
            if processed % 100 == 0:
                print(f"   Analilsado: {processed}/{len(files)}", end='\r')
                
    results.sort(key=lambda x: x[0])
    return results[:top_n], env_coords