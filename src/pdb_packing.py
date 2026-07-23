#!/usr/bin/env python3
"""
Multi-structure packing pipeline.
Uses penalized ICP for spatial distribution, but calculates 
pure Chamfer distance for objective scientific reporting.
"""
import os
import time
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure
from pathlib import Path
from scipy.spatial import cKDTree
import icp_lib 

def calculate_chamfer(setA, setB):
    if len(setA) == 0 or len(setB) == 0:
        return float('inf')
        
    treeA, treeB = cKDTree(setA), cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1) 
    dists_B_to_A, _ = treeA.query(setB, k=1) 
    
    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0

def run_packing_pipeline(input_folder, envelope_path, output_folder, args):
    start_time = time.time()
    os.makedirs(output_folder, exist_ok=True)
    
    top_n = args.max_structures
    penalty_weight = getattr(args, 'penalty', 50.0)
    
    print(f"Packing pipeline started (Top {top_n}, Penalty: {penalty_weight})")

    # Delegate non-colliding candidate search to icp_lib
    candidates, env_coords = icp_lib.get_best_candidates(
        input_folder, envelope_path, top_n=top_n, penalty=penalty_weight
    )
    
    if not candidates:
        print("No candidates found.")
        with open(os.path.join(output_folder, "report.txt"), "w") as f:
            f.write("ALIGNMENT REPORT\nStatus: FAILED (No candidates)\n")
        return

    final_structs, summary_data, all_atoms_combined = [], [], []
    real_score_top1 = None
    
    for rank, (penalized_score, path, T) in enumerate(candidates):
        try:
            filename = Path(path).name
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure("X", path)
            atoms = list(struct.get_atoms())
            
            # Apply transformation matrix
            R, t = T[:3, :3], T[:3, 3]
            struct_coords = []
            for atom in atoms:
                new_coord = np.dot(R, atom.coord) + t
                atom.set_coord(new_coord)
                struct_coords.append(new_coord)
            
            all_atoms_combined.extend(struct_coords)
            
            # Recalculate pure chamfer for unbiased tracking
            current_atoms_np = np.array(struct_coords)
            real_dist = calculate_chamfer(current_atoms_np, env_coords)
            
            if rank == 0:
                real_score_top1 = real_dist

            indiv_name = f"RANK_{rank+1:02d}_{filename}"
            io = PDBIO()
            io.set_structure(struct)
            io.save(os.path.join(output_folder, indiv_name))
            
            final_structs.append(struct)
            summary_data.append({"rank": rank+1, "name": filename, "score": real_dist})
            
        except Exception as e:
            print(f"Error processing {path}: {e}")

    # Calculate combined global volumetric score
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

    real_score_top1 = real_score_top1 or float('inf')
    improvement = "YES" if global_score < real_score_top1 else "NO"
    elapsed_time = time.time() - start_time
    env_size = len(env_coords) if env_coords is not None else 0
    
    # Generate final metrics report for extraction script
    report_path = os.path.join(output_folder, "report.txt")
    with open(report_path, "w") as f:
        f.write("ALIGNMENT REPORT\n================\n\n")
        f.write(f"Top 1 Error (Individual): {real_score_top1:.6f} A\n")
        f.write(f"Global Error (Combined):  {global_score:.6f} A\n")
        f.write(f"Improvement:              {improvement}\n")
        f.write(f"Number of Copies:         {len(final_structs)}\n")
        f.write(f"Penalty Weight Used:      {penalty_weight}\n")
        f.write(f"Execution Time:           {elapsed_time:.2f} seconds\n")
        f.write(f"Envelope Size (points):   {env_size}\n\n")
        
        f.write("-" * 50 + "\n")
        f.write(f"{'RANK':<5} | {'CHAMFER (A)':<12} | {'FILE'}\n")
        for item in summary_data:
            f.write(f"#{item['rank']:<4} | {item['score']:.4f}       | {item['name']}\n")
            
    print(f"Report saved: {report_path} | Global: {global_score:.4f} A | Time: {elapsed_time:.2f}s")