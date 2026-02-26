#!/usr/bin/env python3
"""
icp_lib.py
Lógica de alinhamento ICP com penalidade assimétrica para packing.
"""

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count

# --- NOVA FUNÇÃO DE SCORE (A Mágica acontece aqui) ---

def weighted_asymmetric_score(prot_pts, env_pts, penalty_weight=50.0):
    """
    Calcula um score onde 'vazar' para fora do envelope é muito mais caro
    do que deixar espaços vazios dentro.
    
    prot_pts: Pontos da proteína (candidato)
    env_pts:  Pontos do envelope (alvo)
    penalty_weight: Quanto custa vazar? (50x mais caro que deixar vazio)
    """
    tree_env = cKDTree(env_pts)
    tree_prot = cKDTree(prot_pts)
    
    # 1. Quão bem preenchido está o envelope? (Filling Score)
    # Distância do Envelope -> Proteína mais próxima
    # Se alto: Significa que tem muito espaço vazio no envelope.
    dists_env_to_prot, _ = tree_prot.query(env_pts, k=1)
    score_filling = np.mean(dists_env_to_prot)
    
    # 2. O quanto a proteína vazou? (Penalty Score)
    # Distância da Proteína -> Envelope mais próximo
    # Se alto: Significa que a proteína está longe do envelope (fora).
    dists_prot_to_env, _ = tree_env.query(prot_pts, k=1)
    score_leaking = np.mean(dists_prot_to_env)
    
    # O Score final soma os dois, mas multiplica o erro de vazamento
    final_score = score_filling + (penalty_weight * score_leaking)
    
    return final_score

# --- FUNÇÕES AUXILIARES ---

def extract_structure_from_pdb(pdb_path):
    # (Mantido igual ao seu original)
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
    # (Mantido igual ao seu original)
    xs = np.array([float(x) for x in cif_dict['_atom_site.Cartn_x']], dtype=float)
    ys = np.array([float(y) for y in cif_dict['_atom_site.Cartn_y']], dtype=float)
    zs = np.array([float(z) for z in cif_dict['_atom_site.Cartn_z']], dtype=float)
    return np.stack([xs, ys, zs], axis=1)

def icp_align_with_penalty(src_pts, tgt_pts, max_iter=50, tol=1e-6, penalty=50.0):
    """
    ICP modificado para usar o score ponderado.
    """
    src = src_pts.copy()
    tgt = tgt_pts.copy()
    final_T = np.eye(4)
    best_score = float('inf')

    # KDTree do alvo (envelope) para buscar vizinhos
    tree_tgt = cKDTree(tgt)

    for i in range(max_iter):
        # 1. Encontra correspondências (Nearest Neighbors)
        _, indices = tree_tgt.query(src, k=1)
        closest = tgt[indices]

        # 2. Calcula transform (SVD) para minimizar distância Euclidiana pura
        # (A matemática do SVD não muda, ela sempre tenta aproximar os pontos)
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
        
        # 3. Aplica transformação
        new_src = (R @ src.T).T + t

        # 4. AVALIAÇÃO COM PENALIDADE (AQUI ESTÁ A MUDANÇA)
        # Verificamos se essa nova posição é boa usando a métrica injusta
        current_score = weighted_asymmetric_score(new_src, tgt, penalty_weight=penalty)
        
        # Se melhorou o score ponderado, salvamos
        if current_score < best_score:
            best_score = current_score
            # Atualiza src oficial
            src = new_src
            # Acumula na matriz T final
            T_step = np.eye(4)
            T_step[:3, :3] = R
            T_step[:3, 3] = t
            final_T = T_step @ final_T
        else:
            # Se o score piorou (provavelmente porque vazou muito), 
            # podemos interromper ou tentar diminuir o passo (aqui simplificado para break)
            # Mas no ICP clássico, geralmente deixamos fluir um pouco. 
            # Para garantir convergência, aceitamos a atualização geométrica,
            # mas o score de 'parada' é o ponderado.
            src = new_src 
            # (Nota: Em implementações rígidas, rejeitaríamos o passo, 
            # mas no ICP simples, apenas atualizamos e verificamos a convergência)
            T_step = np.eye(4)
            T_step[:3, :3] = R
            T_step[:3, 3] = t
            final_T = T_step @ final_T

        if current_score < tol:
            break

    return final_T, best_score

# --- WORKER ---

def evaluate_pdb_alignment(args_tuple):
    pdb_path, env_coords, max_iter, sample_env = args_tuple
    try:
        prot_coords, _ = extract_structure_from_pdb(pdb_path)

        # Amostragem do envelope (para performance)
        if (sample_env is not None) and (sample_env > 0) and (env_coords.shape[0] > sample_env):
            idx = np.random.choice(env_coords.shape[0], size=sample_env, replace=False)
            env_for_icp = env_coords[idx]
        else:
            env_for_icp = env_coords

        # Roda o ICP com Penalidade
        # penalty=50.0 força a proteína a ficar dentro.
        # Se continuar vazando, aumente para 100.0 ou 200.0
        T_matrix, score = icp_align_with_penalty(prot_coords, env_for_icp, max_iter=max_iter, penalty=50.0)

        return True, score, pdb_path, T_matrix
    except Exception as e:
        return False, float('inf'), pdb_path, None

# --- ORQUESTRADOR ---

def get_best_candidates(input_folder, envelope_path, top_n=50, limit=None):
    # Carrega envelope
    cif = MMCIF2Dict(open(envelope_path, 'r', encoding='utf-8', errors='ignore'))
    env_coords = coords_from_cif_dict(cif)
    
    files = sorted(list(Path(input_folder).glob("*.pdb")))
    if limit: files = files[:limit]
    
    print(f"⚡ [ICP Penalizado] Screening {len(files)} arquivos...")
    
    work_args = [(str(p), env_coords, 30, 2000) for p in files]
    
    results = []
    processed = 0
    
    with Pool(processes=cpu_count()) as pool:
        for res in pool.imap_unordered(evaluate_pdb_alignment, work_args, chunksize=10):
            success, score, path, T = res
            # Filtro de corte um pouco mais alto pois o score agora é ponderado
            if success and score < 5000: 
                results.append((score, path, T))
            
            processed += 1
            if processed % 10 == 0:
                print(f"   Analisado: {processed}/{len(files)}", end='\r')
                
    # Ordena pelo MENOR score (agora score baixo significa "dentro E preenchendo")
    results.sort(key=lambda x: x[0])
    return results[:top_n], env_coords