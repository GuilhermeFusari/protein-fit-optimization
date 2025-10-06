#!/usr/bin/env python3
"""
pdb_packing.py — Empacotamento otimizado baseado em Chamfer Distance
Seleciona automaticamente os melhores PDBs que se encaixam bem no envelope.
"""

import time
import numpy as np
from pathlib import Path
from icp_aligner import find_best_pdb, load_pdb_as_points, save_points_as_pdb
from scipy.spatial import cKDTree


def chamfer_distance(A, B):
    """
    Calcula a Chamfer Distance entre dois conjuntos de pontos (A e B).
    Essa métrica quantifica o quão bem os pontos de um conjunto se aproximam do outro,
    sendo usada para avaliar o alinhamento de proteínas com o envelope.
    """
    tree_A = cKDTree(A)
    tree_B = cKDTree(B)
    dist_AB, _ = tree_A.query(B)
    dist_BA, _ = tree_B.query(A)
    return np.mean(dist_AB**2) + np.mean(dist_BA**2)


def run_packing_pipeline(input_path, envelope_path, output_path, args):
    """
    Pipeline principal do modo 'packing':
    - Carrega uma lista limitada de PDBs.
    - Alinha cada PDB ao envelope individualmente usando ICP.
    - Calcula a Chamfer Distance de cada alinhamento.
    - Seleciona os melhores PDBs com base em um threshold automático.
    - Combina os pontos dos PDBs selecionados e salva em um único arquivo PDB final.
    """
    start = time.time()
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Iniciando modo PACKING (seleção via Chamfer Distance)...")

    pdb_files = list(input_path.glob("*.pdb"))
    if not pdb_files:
        print("❌ Nenhum PDB encontrado.")
        return

    # Limita ao número máximo de estruturas definido pelo usuário
    pdb_files = pdb_files[:args.max_structures]
    print(f"📦 Total de {len(pdb_files)} candidatos")
    envelope_points = load_pdb_as_points(str(envelope_path))

    resultados = []  # Lista de tuplas: (pdb_name, aligned_points, chamfer_value)

    # Itera sobre cada PDB, alinhando e avaliando sua qualidade
    for i, pdb_file in enumerate(pdb_files, 1):
        print(f"\n📄 [{i}/{len(pdb_files)}] Testando {pdb_file.name}")

        try:
            class Args:
                input = str(pdb_file)
                envelope = str(envelope_path)
                output = str(output_path)
                align_what = 'protein'
                max_iter = args.max_iter
                sample_env = 5000
                workers = 1

            aligned_points, transform = find_best_pdb(Args, return_transform=True)
            cd = chamfer_distance(aligned_points, envelope_points)
            print(f"🔹 Chamfer Distance = {cd:.6f}")
            resultados.append((pdb_file.name, aligned_points, cd))

        except Exception as e:
            print(f"⚠️ Erro ao processar {pdb_file.name}: {e}")

    if not resultados:
        print("\n❌ Nenhum resultado obtido.")
        return

    # Ordena os resultados pelo menor Chamfer Distance
    resultados.sort(key=lambda x: x[2])

    # Define threshold automático para selecionar os melhores PDBs (top 30%)
    valores_cd = [r[2] for r in resultados]
    threshold = np.percentile(valores_cd, 30)
    print(f"\n📊 Threshold automático definido em {threshold:.6f}")

    selecionados = [r for r in resultados if r[2] <= threshold]

    if not selecionados:
        print("\n❌ Nenhum PDB dentro do limite.")
        return

    print(f"\n🏁 {len(selecionados)} PDBs selecionados para empacotamento final")

    # Combina todos os PDBs selecionados em uma única nuvem de pontos
    combined_points = np.vstack([p[1] for p in selecionados])
    save_points_as_pdb(combined_points, str(output_path / "packing_result.pdb"))

    print(f"\n✅ Resultado final salvo em: {output_path / 'packing_result.pdb'}")
    print(f"⏱️ Tempo total: {(time.time() - start):.1f}s")
