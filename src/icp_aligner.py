#!/usr/bin/env python3
"""
Script otimizado para encontrar o melhor alinhamento entre proteínas PDB e um envelope estrutural
usando ICP (Iterative Closest Point) com Chamfer Distance + multiprocessing.
"""

import os
import argparse
import numpy as np
import time
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count


def extract_structure_from_pdb(pdb_path):
    """
    Carrega um arquivo PDB usando Biopython e extrai as coordenadas de todos os átomos.
    Retorna uma matriz Nx3 de coordenadas e o objeto de estrutura completo.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PROT", pdb_path)
    coords = np.array([atom.coord for atom in structure.get_atoms()], dtype=float)
    return coords, structure


def coords_from_cif_dict(cif_dict):
    """
    Converte os dados de um dicionário MMCIF em um array Nx3 de coordenadas atômicas.
    """
    xs = np.array([float(x) for x in cif_dict['_atom_site.Cartn_x']], dtype=float)
    ys = np.array([float(y) for y in cif_dict['_atom_site.Cartn_y']], dtype=float)
    zs = np.array([float(z) for z in cif_dict['_atom_site.Cartn_z']], dtype=float)
    return np.stack([xs, ys, zs], axis=1)


def chamfer_distance(setA, setB):
    """
    Calcula a Chamfer Distance entre dois conjuntos de pontos.
    Mede a similaridade entre a forma da proteína e do envelope.
    """
    treeA = cKDTree(setA)
    treeB = cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1)
    dists_B_to_A, _ = treeA.query(setB, k=1)
    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0


def icp_align_with_chamfer(src_pts, tgt_pts, max_iter=50, tol=1e-6):
    """
    Executa o algoritmo ICP (Iterative Closest Point) rigidamente para alinhar src_pts a tgt_pts.
    Atualiza uma matriz de transformação 4x4 iterativamente até convergência.
    Retorna a transformação final e o score de Chamfer Distance.
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
        src = (R @ src.T).T + t

        T = np.eye(4)
        T[:3, :3] = R
        T[:3,  3] = t
        final_T = T @ final_T

        if current_score < tol:
            break

    final_score = chamfer_distance(src, tgt)
    return final_T, final_score


def apply_transform_to_structure(structure, T4):
    """
    Aplica uma transformação 4x4 (R+t) diretamente nos átomos da estrutura Biopython.
    """
    R = T4[:3, :3]; t = T4[:3, 3]
    for atom in structure.get_atoms():
        atom.set_coord(R @ atom.coord + t)


def evaluate_pdb_alignment(args_tuple):
    """
    Função auxiliar para processamento paralelo.
    Recebe um PDB, realiza o alinhamento ICP com o envelope e retorna o score de Chamfer.
    """
    pdb_path, env_coords, align_what, max_iter, sample_env = args_tuple
    try:
        prot_coords, _ = extract_structure_from_pdb(pdb_path)

        if (sample_env is not None) and (sample_env > 0) and (env_coords.shape[0] > sample_env):
            idx = np.random.choice(env_coords.shape[0], size=sample_env, replace=False)
            env_for_icp = env_coords[idx]
        else:
            env_for_icp = env_coords

        if align_what == 'protein':
            _, score = icp_align_with_chamfer(prot_coords, env_for_icp, max_iter=max_iter)
        else:
            _, score = icp_align_with_chamfer(env_for_icp, prot_coords, max_iter=max_iter)

        return True, score, os.path.basename(pdb_path), pdb_path
    except Exception as e:
        return False, float('inf'), f"{os.path.basename(pdb_path)}: {str(e)}", pdb_path


def find_best_pdb(args):
    """
    Função principal que gerencia todo o fluxo do alinhamento:
    - Cria a pasta de saída.
    - Carrega o envelope e os arquivos PDB.
    - Avalia todos os PDBs em paralelo usando multiprocessing.
    - Seleciona o melhor PDB pelo menor Chamfer Distance.
    - Salva o PDB alinhado e imprime relatório detalhado.
    """
    os.makedirs(args.output, exist_ok=True)

    print("🔵 Procurando melhor PDB para o envelope...")
    print(f"📁 PDBs: {args.input}")
    print(f"📁 Envelope: {args.envelope}")
    print(f"📁 Saída: {args.output}")
    print(f"🎯 Alinhamento: {args.align_what}")
    print(f"📊 Métrica: Chamfer Distance")
    print(f"⚡ Usando {args.workers} processos em paralelo")

    print("📖 Carregando envelope...")
    cif_dict = MMCIF2Dict(open(args.envelope, 'r', encoding='utf-8', errors='ignore'))
    env_coords = coords_from_cif_dict(cif_dict)
    print(f"📊 Envelope com {env_coords.shape[0]} pontos")

    input_path = Path(args.input)
        
    if input_path.is_file():
        # Se o usuário passou um arquivo específico, usamos apenas ele
        pdb_files = [input_path]
    elif input_path.is_dir():
        # Se passou uma pasta, procuramos todos os .pdb lá dentro
        pdb_files = list(input_path.glob("*.pdb"))
    else:
        print(f"❌ Erro: O caminho '{args.input}' não existe ou não é válido.")
        return

    print(f"📊 Encontrados {len(pdb_files)} arquivos PDB")

    if not pdb_files:
        print("❌ Nenhum PDB encontrado! Verifique o caminho ou a extensão dos arquivos.")
        return

    # --- estimativa de tempo com 1 PDB ---
    t0 = time.time()
    evaluate_pdb_alignment((str(pdb_files[0]), env_coords, args.align_what, args.max_iter, args.sample_env))
    single_time = time.time() - t0
    est_total = (single_time * len(pdb_files)) / args.workers
    print(f"⏱️  Estimativa de tempo total: ~{est_total/60:.1f} minutos ({est_total:.1f} segundos)")

    start_time = time.time()

    # --- processamento paralelo ---
    work_args = [(str(p), env_coords, args.align_what, args.max_iter, args.sample_env) for p in pdb_files]

    best_score = float('inf')
    best_pdb_path = None
    best_pdb_name = None
    processed_count = 0
    error_count = 0
    errors = []

    with Pool(processes=args.workers) as pool:
        for success, score, pdb_name, pdb_path in pool.imap_unordered(evaluate_pdb_alignment, work_args, chunksize=10):
            if success:
                processed_count += 1
                if score < best_score:
                    best_score = score
                    best_pdb_path = pdb_path
                    best_pdb_name = pdb_name
                    print(f"🎯 NOVO MELHOR: {pdb_name} (Chamfer: {score:.4f}Å)")
            else:
                error_count += 1
                errors.append(pdb_name)

    # salvar melhor resultado
    if best_pdb_path:
        print(f"\n🏆 MELHOR PDB ENCONTRADO: {best_pdb_name}")
        print(f"📊 Chamfer Distance: {best_score:.6f}Å")

        prot_coords, prot_struct = extract_structure_from_pdb(best_pdb_path)
        if (args.sample_env is not None) and (args.sample_env > 0) and (env_coords.shape[0] > args.sample_env):
            idx = np.random.choice(env_coords.shape[0], size=args.sample_env, replace=False)
            env_for_icp = env_coords[idx]
        else:
            env_for_icp = env_coords

        if args.align_what == 'protein':
            T, final_score = icp_align_with_chamfer(prot_coords, env_for_icp, max_iter=args.max_iter)
            apply_transform_to_structure(prot_struct, T)

        output_name = f"BEST_ALIGNMENT_CHAMFER_{os.path.basename(best_pdb_path)}"
        output_path = os.path.join(args.output, output_name)

        io = PDBIO()
        io.set_structure(prot_struct)
        io.save(output_path)

        print(f"💾 Proteína alinhada salva em: {output_path}")
        print(f"🎯 Score final: {final_score:.6f}Å")
    else:
        print("❌ Nenhum PDB válido encontrado!")

    total_time = time.time() - start_time
    print("\n" + "="*50)
    print("📊 RELATÓRIO FINAL - CHAMFER DISTANCE")
    print("="*50)
    print(f"🏆 Melhor PDB: {best_pdb_name}")
    print(f"📊 Melhor Chamfer Distance: {best_score:.6f}Å")
    print(f"✅ PDBs avaliados: {processed_count}")
    print(f"❌ Erros: {error_count}")
    print(f"⏰ Tempo total: {total_time/60:.2f} minutos")
    print(f"⏱️  Tempo médio por PDB: {total_time/processed_count:.4f} segundos")
