#!/usr/bin/env python3
"""
pdb_packing.py
Pipeline de empacotamento que separa a métrica de otimização (com penalidade)
da métrica de relatório (distância real em Angstroms).
"""

import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure
from pathlib import Path
from scipy.spatial import cKDTree
import icp_lib 

"""
BLOCO 1: FUNÇÕES AUXILIARES E DEPENDÊNCIAS
A função calculate_chamfer implementa a métrica de distância 'Chamfer' pura.
Isso é utilizado aqui para quantificar a distância REAL (física) em Angstroms
para o relatório final, ignorando as penalidades usadas durante a otimização.
"""
def calculate_chamfer(setA, setB):
    if len(setA) == 0 or len(setB) == 0:
        return float('inf')
        
    treeA = cKDTree(setA)
    treeB = cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1) 
    dists_B_to_A, _ = treeA.query(setB, k=1) 
    
    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0

def run_packing_pipeline(input_folder, envelope_path, output_folder, args=None):
    """
    BLOCO 2: CONFIGURAÇÃO E SCREENING
    """
    os.makedirs(output_folder, exist_ok=True)
    
    top_n = args.max_structures if args else 5
    limit = args.input_limit if args and hasattr(args, 'input_limit') else None
    
    print(f"Starting alignment pipeline (Top {top_n})...")

    # Aqui o icp_lib retorna os candidatos baseados no Score com Penalidade
    candidates, env_coords = icp_lib.get_best_candidates(
        input_folder, envelope_path, top_n=top_n, limit=limit
    )
    
    if not candidates:
        print("No candidates found.")
        return

    """
    BLOCO 3: TRANSFORMAÇÃO E ANÁLISE GLOBAL
    Aqui acontece a 'Limpeza' do score para o relatório.
    """
    final_structs = []
    summary_data = []
    all_atoms_combined = []
    
    # Variável para guardar o score REAL (Angstroms) do melhor candidato
    real_score_top1 = None
    
    print("Processing files...")
    
    for rank, (penalized_score, path, T) in enumerate(candidates):
        try:
            filename = Path(path).name
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure("X", path)
            atoms = list(struct.get_atoms())
            
            R = T[:3, :3]; t = T[:3, 3]
            
            struct_coords = []
            for atom in atoms:
                coord = atom.coord
                new_coord = np.dot(R, coord) + t
                atom.set_coord(new_coord)
                struct_coords.append(new_coord)
            
            all_atoms_combined.extend(struct_coords)
            
            # --- O PULO DO GATO ---
            # Recalculamos o Chamfer PURO aqui.
            # Ignoramos o 'penalized_score' que veio do icp_lib para fins de relatório.
            current_atoms_np = np.array(struct_coords)
            real_distance_angstrom = calculate_chamfer(current_atoms_np, env_coords)
            
            # Se for o Rank 1, guardamos esse valor real para comparar no final
            if rank == 0:
                real_score_top1 = real_distance_angstrom

            individual_name = f"RANK_{rank+1:02d}_{filename}"
            io = PDBIO()
            io.set_structure(struct)
            io.save(os.path.join(output_folder, individual_name))
            
            final_structs.append(struct)
            
            # Salvamos no resumo a Distância Real, não a penalizada
            summary_data.append({"rank": rank+1, "name": filename, "score": real_distance_angstrom})
            
            print(f"Rank {rank+1}: {filename} (Real Dist: {real_distance_angstrom:.4f} A | Opt Score: {penalized_score:.1f})")
            
        except Exception as e:
            print(f"Error processing {path}: {e}")

    print("Calculating Global Volumetric Score...")
    
    all_atoms_np = np.array(all_atoms_combined)
    global_score = calculate_chamfer(all_atoms_np, env_coords)

    if final_structs:
        try:
            io = PDBIO()
            c = Structure.Structure("Combined")
            for i, m in enumerate(final_structs):
                model = next(iter(m)) 
                model.id = i
                c.add(model)
            io.set_structure(c)
            io.save(os.path.join(output_folder, "ALL_TOP_ALIGNED.pdb"))
        except Exception as e:
            print(f"Warning: Could not save combined file: {e}")

    # Proteção caso a lista venha vazia
    if real_score_top1 is None: 
        real_score_top1 = float('inf')

    improvement = "YES" if global_score < real_score_top1 else "NO"
    
    report_path = os.path.join(output_folder, "report.txt")
    with open(report_path, "w") as f:
        f.write("ALIGNMENT REPORT\n")
        f.write("================\n\n")
        
        f.write(f"Top 1 Error (Individual): {real_score_top1:.6f} A\n")
        f.write(f"Global Error (Combined):  {global_score:.6f} A\n")
        f.write(f"Improvement: {improvement}\n\n")
        
        f.write("-" * 50 + "\n")
        f.write(f"{'RANK':<5} | {'CHAMFER (A)'}\n")
        for item in summary_data:
            f.write(f"#{item['rank']:<4} | {item['score']:.4f}\n")
            
    print(f"Report saved to: {report_path}")
    print("-" * 30)
    print(f"Top 1 Error:    {real_score_top1:.4f} A")
    print(f"Global Error:   {global_score:.4f} A")
    print(f"Improvement:    {improvement}")

if __name__ == "__main__":
    pass