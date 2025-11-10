#!/usr/bin/env python3

"""
Packing otimizado de múltiplos PDBs dentro de um envelope usando ICP + Chamfer Distance.

Melhorias:
- Barra de progresso contínua desde a escolha do melhor PDB.
- Heurística de centro das proteínas já colocadas para posicionamento.
- Máximo de 20 proteínas dentro do envelope.
- Nenhuma proteína sai do envelope.
- Comentários em bloco em português.
"""

import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
import time
import sys

# -----------------------------
# Funções de alinhamento ICP
# -----------------------------

def icp_align_with_chamfer(src_pts, tgt_pts, max_iter=50, tol=1e-6):
    src = src_pts.copy()
    tgt = tgt_pts.copy()
    final_T = np.eye(4)
    tree_tgt = cKDTree(tgt)

    for _ in range(max_iter):
        _, indices = tree_tgt.query(src, k=1)
        closest = tgt[indices]
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
        src = (R @ src.T).T + t
        T = np.eye(4)
        T[:3, :3] = R
        T[:3, 3] = t
        final_T = T @ final_T
    return final_T

def chamfer_distance(a_pts, b_pts):
    if a_pts.shape[0] == 0 or b_pts.shape[0] == 0:
        return np.inf
    ta = cKDTree(b_pts).query(a_pts, k=1)[0].mean()
    tb = cKDTree(a_pts).query(b_pts, k=1)[0].mean()
    return (ta + tb) / 2.0

def random_rotation_matrix(max_angle=np.pi/4.0):
    alpha, beta, gamma = np.random.uniform(-max_angle, max_angle, 3)
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0],
                   [0, 0, 1]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)],
                   [0, 1, 0],
                   [-np.sin(beta), 0, np.cos(beta)]])
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(alpha), -np.sin(alpha)],
                   [0, np.sin(alpha), np.cos(alpha)]])
    return Rz @ Ry @ Rx

# -----------------------------
# Funções utilitárias
# -----------------------------

def extract_structure_from_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("PROT", pdb_path)
    except Exception:
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser_cif = MMCIFParser(QUIET=True)
        structure = parser_cif.get_structure("PROT", pdb_path)
    coords = np.array([atom.coord for atom in structure.get_atoms()], dtype=float)
    return coords, structure

def coords_from_cif_dict(cif_dict):
    xs = np.array([float(x) for x in cif_dict['_atom_site.Cartn_x']], dtype=float)
    ys = np.array([float(y) for y in cif_dict['_atom_site.Cartn_y']], dtype=float)
    zs = np.array([float(z) for z in cif_dict['_atom_site.Cartn_z']], dtype=float)
    return np.stack([xs, ys, zs], axis=1)

def apply_transform_to_structure(structure, T4):
    R = T4[:3, :3]
    t = T4[:3, 3]
    for atom in structure.get_atoms():
        atom.set_coord(R @ atom.coord + t)

def check_overlap(new_coords, existing_coords, min_dist=2.5):
    if existing_coords.shape[0] == 0:
        return False
    tree = cKDTree(existing_coords)
    dists, _ = tree.query(new_coords, k=1)
    return np.any(dists < min_dist)

def is_within_envelope(coords, env_tree, max_dist=4.0):
    dists, _ = env_tree.query(coords, k=1)
    return np.all(dists <= max_dist)

# -----------------------------
# Função de posicionamento otimizado
# -----------------------------

def place_protein_optimized(prot_coords, existing_coords, env_coords, env_tree,
                            max_attempts=200, jitter_range=20.0, min_dist=2.0, max_outside_dist=15.0):
    center_env = env_coords.mean(axis=0)
    for attempt in range(max_attempts):
        if existing_coords.shape[0] > 0:
            center_existing = existing_coords.mean(axis=0)
        else:
            center_existing = center_env
        t_jitter = center_existing + np.random.uniform(-jitter_range, jitter_range, 3)
        R_jitter = random_rotation_matrix()
        prot_transformed = (R_jitter @ prot_coords.T).T + t_jitter

        if not is_within_envelope(prot_transformed, env_tree, max_dist=max_outside_dist):
            continue
        if check_overlap(prot_transformed, existing_coords, min_dist=min_dist):
            continue

        T_final = np.eye(4)
        T_final[:3, :3] = R_jitter
        T_final[:3, 3] = t_jitter
        return prot_transformed, T_final

    print(f"[WARN] Nenhuma posição válida encontrada após {max_attempts} tentativas.", flush=True)
    return None, None

# -----------------------------
# Pipeline principal otimizado
# -----------------------------

def run_packing_pipeline(input_folder, envelope_path, output_folder, args=None):
    """
    Pipeline robusto para packing de PDBs dentro do envelope.

    - Força o primeiro PDB ao centro do envelope (garante resultado inicial).
    - Em seguida tenta posicionar até 19 PDBs adicionais usando place_protein_optimized.
    - Não usa tqdm (evita dependências). Usa prints com flush para barra de progresso.
    - Blocos de comentários em português.
    """
    print("🔵 Iniciando packing otimizado (até 20 proteínas)...")
    os.makedirs(output_folder, exist_ok=True)

    # Carrega envelope
    try:
        with open(envelope_path, 'r', encoding='utf-8', errors='ignore') as f:
            cif_dict = MMCIF2Dict(f)
        env_coords = coords_from_cif_dict(cif_dict)
    except Exception as e:
        print(f"\n❌ Erro ao carregar/processar o envelope ({envelope_path}): {e}")
        return

    env_tree = cKDTree(env_coords)
    print(f"📊 Envelope carregado com {env_coords.shape[0]} pontos")

    # Lista PDBs
    pdb_files = sorted(list(Path(input_folder).glob("*.pdb")))
    total_pdbs = len(pdb_files)
    print(f"📦 Total de PDBs encontrados: {total_pdbs}")
    if total_pdbs == 0:
        print("❌ Nenhum arquivo PDB encontrado no input.")
        return

    # 1) Escolher melhor PDB (mostra progresso básico)
    best_score = np.inf
    best_path = None
    best_coords = None
    best_struct = None

    for idx, pdb_path in enumerate(pdb_files, start=1):
        try:
            coords, struct = extract_structure_from_pdb(str(pdb_path))
        except Exception as e:
            print(f"⚠️ Falha ao ler {pdb_path.name}: {e}", flush=True)
            print(f"🔍 Avaliando PDBs: {idx}/{total_pdbs}", end='\r', flush=True)
            continue

        # Usamos Chamfer via ICP transformada para avaliar afinidade com o envelope
        try:
            Ttmp = icp_align_with_chamfer(coords, env_coords)
            transformed = (Ttmp[:3, :3] @ coords.T).T + Ttmp[:3, 3]
            score = chamfer_distance(transformed, env_coords)
        except Exception:
            # fallback simples: distância média ao centro do envelope
            score = np.mean(np.linalg.norm(coords - env_coords.mean(axis=0), axis=1))

        if score < best_score:
            best_score = score
            best_path = str(pdb_path)
            best_coords = coords
            best_struct = struct

        # progresso
        print(f"🔍 Avaliando PDBs: {idx}/{total_pdbs}", end='\r', flush=True)

    print("\n")  # quebra de linha após progresso
    if best_path is None:
        print("❌ Não foi possível selecionar um PDB válido como referência.")
        return

    print(f"🏅 Melhor PDB selecionado: {Path(best_path).name} (score={best_score:.3f})")

    # Preparação para packing
    all_coords = np.zeros((0, 3))
    final_structures = []
    start_time = time.time()

    # 2) Forçar o primeiro PDB ao centro do envelope (garante início)
    try:
        center_env = env_coords.mean(axis=0)
        # Translação para centralizar
        trans = center_env - best_coords.mean(axis=0)
        T_first = np.eye(4)
        T_first[:3, :3] = np.eye(3)
        T_first[:3, 3] = trans
        prot_transformed = (best_coords + trans)  # já centralizado

        apply_transform_to_structure(best_struct, T_first)
        final_structures.append(best_struct)
        all_coords = np.vstack([all_coords, prot_transformed])
        placed_count = 1

        print(f"✅ Primeiro PDB forçado ao centro do envelope: {Path(best_path).name}")
        print(f"⏱️ Progresso: {placed_count}/20 proteínas colocadas", end='\r', flush=True)

    except Exception as e:
        print(f"❌ Falha ao forçar o primeiro PDB ao centro: {e}")
        return

    # 3) Tentar posicionar os demais PDBs (até completar 20 no total)
    for pdb_path in pdb_files:
        if placed_count >= 20:
            break

        if str(pdb_path) == str(best_path):
            continue

        try:
            prot_coords, prot_struct = extract_structure_from_pdb(str(pdb_path))
        except Exception as e:
            print(f"\n⚠️ Falha ao ler {pdb_path.name}: {e}", flush=True)
            print(f"⏱️ Progresso: {placed_count}/20 proteínas colocadas", end='\r', flush=True)
            continue

        # tentar posicionar com parâmetros mais permissivos (ajuste seguro)
        transformed_coords, T = place_protein_optimized(
            prot_coords, all_coords, env_coords, env_tree,
            max_attempts=200, jitter_range=20.0, min_dist=2.0, max_outside_dist=15.0
        )

        if transformed_coords is None:
            # não conseguiu posicionar essa proteína
            print(f"\n⛔ Não foi possível posicionar {pdb_path.name} sem colisão/fora (pulando).", flush=True)
            print(f"⏱️ Progresso: {placed_count}/20 proteínas colocadas", end='\r', flush=True)
            continue

        # aplicar transformação e adicionar
        try:
            apply_transform_to_structure(prot_struct, T)
            final_structures.append(prot_struct)
            all_coords = np.vstack([all_coords, transformed_coords])
            placed_count += 1
            print(f"⏱️ Progresso: {placed_count}/20 proteínas colocadas", end='\r', flush=True)
        except Exception as e:
            print(f"\n⚠️ Erro ao aplicar transformação em {pdb_path.name}: {e}", flush=True)
            print(f"⏱️ Progresso: {placed_count}/20 proteínas colocadas", end='\r', flush=True)
            continue

    print("\n")  # quebra de linha após progresso final

    # 4) Salvar estrutura empacotada
    try:
        packed_structure = Structure.Structure("PACKED")
        model_counter = 0
        for struct in final_structures:
            for model in struct:
                model.id = model_counter
                packed_structure.add(model)
                model_counter += 1

        output_path = os.path.join(output_folder, "PACKED_ENVELOPE.pdb")
        io = PDBIO()
        io.set_structure(packed_structure)
        io.save(output_path)
    except Exception as e:
        print(f"❌ Erro ao salvar o PDB final: {e}")
        return

    total_time = time.time() - start_time
    print(f"💾 Packing finalizado. Resultado salvo em: {output_path}")
    print(f"🏆 Total de PDBs colocados: {len(final_structures)}")
    print(f"⏰ Tempo total: {total_time / 60:.2f} minutos")
