#!/usr/bin/env python3
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure
from pathlib import Path
from scipy.spatial import cKDTree
import icp_lib 

"""
BLOCO 1: FUNÇÕES AUXILIARES E DEPENDÊNCIAS
Este script utiliza BioPython para manipulação de estruturas PDB e SciPy para cálculos espaciais.
A função calculate_chamfer implementa a métrica de distância 'Chamfer', que calcula a média
das distâncias mínimas entre dois conjuntos de pontos (A->B e B->A).
Isso é utilizado para quantificar o quão bem o conjunto de proteínas preenche o envelope.
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
    Aqui o pipeline é inicializado. Criamos a pasta de saída e definimos o número de estruturas (top_n).
    A função chama a biblioteca externa (icp_lib) para realizar o alinhamento pesado via multiprocessing.
    O retorno 'candidates' contém as N melhores estruturas com seus respectivos Scores e Matrizes de Transformação (T).
    """
    os.makedirs(output_folder, exist_ok=True)
    
    top_n = args.max_structures if args else 5
    limit = args.input_limit if args and hasattr(args, 'input_limit') else None
    
    print(f"Starting alignment pipeline (Top {top_n})...")

    candidates, env_coords = icp_lib.get_best_candidates(
        input_folder, envelope_path, top_n=top_n, limit=limit
    )
    
    if not candidates:
        print("No candidates found.")
        return

    """
    BLOCO 3: TRANSFORMAÇÃO E ANÁLISE GLOBAL
    Iteramos sobre os candidatos vencedores. Para cada estrutura:
    1. Extraímos a rotação (R) e translação (t) da matriz T.
    2. Aplicamos essa transformação átomo a átomo.
    3. Salvamos o arquivo individual e acumulamos os átomos para o arquivo combinado.
    
    Ao final, calculamos o Score Global (todos os átomos combinados contra o envelope)
    para comparar se o conjunto preenche melhor o volume do que a melhor proteína individualmente.
    """
    final_structs = []
    summary_data = []
    all_atoms_combined = []
    
    score_top1 = candidates[0][0]
    
    print("Processing files...")
    
    for rank, (score, path, T) in enumerate(candidates):
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
            
            individual_name = f"RANK_{rank+1:02d}_{filename}"
            io = PDBIO()
            io.set_structure(struct)
            io.save(os.path.join(output_folder, individual_name))
            
            final_structs.append(struct)
            summary_data.append({"rank": rank+1, "name": filename, "score": score})
            
            print(f"Rank {rank+1}: {filename} (Error: {score:.4f})")
            
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

    improvement = "YES" if global_score < score_top1 else "NO"
    
    report_path = os.path.join(output_folder, "report.txt")
    with open(report_path, "w") as f:
        f.write("ALIGNMENT REPORT\n")
        f.write("================\n\n")
        
        f.write(f"Top 1 Error (Individual): {score_top1:.6f} A\n")
        f.write(f"Global Error (Combined):  {global_score:.6f} A\n")
        f.write(f"Improvement: {improvement}\n\n")
        
        f.write("-" * 50 + "\n")
        f.write(f"{'RANK':<5} | {'SCORE'}\n")
        for item in summary_data:
            f.write(f"#{item['rank']:<4} | {item['score']:.4f}\n")
            
    print(f"Report saved to: {report_path}")
    print("-" * 30)
    print(f"Top 1 Error:    {score_top1:.4f}")
    print(f"Global Error:   {global_score:.4f}")
    print(f"Improvement:    {improvement}")

if __name__ == "__main__":
    pass