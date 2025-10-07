from multiprocessing import Pool, cpu_count
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pathlib import Path
from scipy.spatial import cKDTree
import time

def extract_structure_from_pdb(pdb_path):
    """
    Carrega um arquivo PDB usando Biopython e extrai as coordenadas atômicas.

    Retorna:
    - coords: matriz Nx3 com coordenadas de todos os átomos.
    - structure: objeto de estrutura Biopython completo.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PROT", pdb_path)
    coords = np.array([atom.coord for atom in structure.get_atoms()], dtype=float)
    return coords, structure

def coords_from_cif_dict(cif_dict):
    """
    Converte os dados de um dicionário MMCIF em um array Nx3 de coordenadas atômicas.

    Retorna:
    - coords: matriz Nx3 com coordenadas do envelope.
    """
    xs = np.array([float(x) for x in cif_dict['_atom_site.Cartn_x']], dtype=float)
    ys = np.array([float(y) for y in cif_dict['_atom_site.Cartn_y']], dtype=float)
    zs = np.array([float(z) for z in cif_dict['_atom_site.Cartn_z']], dtype=float)
    return np.stack([xs, ys, zs], axis=1)

def rotation_matrix(axis, theta):
    """
    Gera uma matriz de rotação 3x3 para um eixo e ângulo fornecidos.

    Args:
    - axis: vetor de 3 dimensões representando o eixo de rotação.
    - theta: ângulo de rotação em radianos.

    Retorna:
    - Matriz 3x3 de rotação.
    """
    axis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2)
    b,c,d = -axis*np.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def apply_transform(structure, T):
    """
    Aplica uma transformação 4x4 (R+t) diretamente nos átomos de uma estrutura Biopython.

    Args:
    - structure: objeto Biopython da proteína.
    - T: matriz 4x4 contendo rotação (3x3) e translação (3x1).
    """
    R = T[:3,:3]; t = T[:3,3]
    for atom in structure.get_atoms():
        atom.set_coord(R @ atom.coord + t)

def chamfer_distance(setA, setB):
    """
    Calcula a Chamfer Distance entre dois conjuntos de pontos.

    Retorna:
    - score médio da Chamfer Distance, indicando a similaridade entre formas.
    """
    treeA = cKDTree(setA)
    treeB = cKDTree(setB)
    dists_A_to_B, _ = treeB.query(setA, k=1)
    dists_B_to_A, _ = treeA.query(setB, k=1)
    return (np.mean(dists_A_to_B) + np.mean(dists_B_to_A)) / 2.0

def check_overlap(new_coords, existing_coords, min_dist=1.5):
    """
    Verifica se as coordenadas de um novo PDB sobrepõem-se a PDBs já posicionados.

    Retorna:
    - True se houver sobreposição (distância mínima menor que min_dist).
    - False caso contrário.
    """
    if existing_coords.shape[0] == 0:
        return False
    tree = cKDTree(existing_coords)
    dists, _ = tree.query(new_coords, k=1)
    return np.any(dists < min_dist)

def place_pdb(args_tuple):
    """
    Tenta posicionar um PDB dentro do envelope sem sobreposição.

    Estratégia:
    - Aplica rotações aleatórias.
    - Centraliza o PDB no envelope.
    - Verifica sobreposição com PDBs já posicionados.
    - Calcula score de preenchimento baseado na Chamfer Distance.

    Retorna:
    - pdb_path: caminho do PDB processado.
    - best_transform: melhor transformação 4x4 encontrada.
    - pdb_struct: objeto Biopython da estrutura transformada.
    - best_score: score associado à melhor posição.
    """
    pdb_path, envelope_coords, existing_coords, n_rotations = args_tuple
    pdb_coords, pdb_struct = extract_structure_from_pdb(pdb_path)
    best_score = -1
    best_transform = None
    for _ in range(n_rotations):
        angle = np.random.rand() * 2*np.pi
        axis = np.random.rand(3) - 0.5
        axis /= np.linalg.norm(axis)
        R = rotation_matrix(axis, angle)
        transformed = (R @ pdb_coords.T).T
        centroid = transformed.mean(axis=0)
        env_centroid = envelope_coords.mean(axis=0)
        transformed += (env_centroid - centroid)
        if check_overlap(transformed, existing_coords):
            continue
        score = -chamfer_distance(transformed, envelope_coords)
        if score > best_score:
            best_score = score
            T = np.eye(4)
            T[:3,:3] = R
            T[:3,3] = env_centroid - centroid
            best_transform = T
    return pdb_path, best_transform, pdb_struct, best_score

def run_packing_pipeline(input_folder, envelope_path, output_folder, args):
    """
    Função principal que executa o packing de múltiplos PDBs:

    Fluxo:
    1. Carrega envelope CIF e lista de arquivos PDB.
    2. Processa todos os PDBs em paralelo, aplicando rotações e posicionando dentro do envelope.
    3. Evita sobreposição usando Chamfer Distance.
    4. Mostra progresso resumido e estimativa de tempo restante.
    5. Salva todas as estruturas posicionadas em um único arquivo de saída.

    Args:
    - input_folder: pasta com arquivos PDB.
    - envelope_path: arquivo CIF do envelope.
    - output_folder: pasta onde salvar o resultado final.
    - args: argumentos adicionais, como número de rotações.
    """
    os.makedirs(output_folder, exist_ok=True)
    cif_dict = MMCIF2Dict(open(envelope_path,'r',encoding='utf-8',errors='ignore'))
    envelope_coords = coords_from_cif_dict(cif_dict)
    pdb_files = list(Path(input_folder).glob("*.pdb"))

    total_pdbs = len(pdb_files)
    workers = cpu_count()  # utiliza todos os núcleos disponíveis
    print(f"📦 Total de PDBs: {total_pdbs}")
    print(f"⚡ Usando {workers} cores para o packing\n")

    all_coords = np.zeros((0,3))
    final_structures = []
    start_time = time.time()

    # Prepara argumentos para multiprocessing
    work_args = [(p, envelope_coords, all_coords, args.n_rotations) for p in pdb_files]

    # Processamento paralelo com estimativa de progresso
    with Pool(processes=workers) as pool:
        for i, (pdb_path, T, pdb_struct, score) in enumerate(pool.imap_unordered(place_pdb, work_args, chunksize=10), 1):
            if T is not None:
                apply_transform(pdb_struct, T)
                final_structures.append(pdb_struct)
                all_coords = np.vstack([all_coords, np.array([atom.coord for atom in pdb_struct.get_atoms()])])
            # Mostra progresso resumido a cada 1% completado
            if i % max(1, total_pdbs//100) == 0 or i == total_pdbs:
                elapsed = time.time() - start_time
                est_total = elapsed / i * total_pdbs
                print(f"⏱️  Processados: {i}/{total_pdbs} | Tempo estimado restante: {est_total - elapsed:.1f}s", end='\r')

    # Salva resultado final
    output_path = os.path.join(output_folder, "PACKED_ENVELOPE.pdb")
    io = PDBIO()
    class MultiStructure:
        def get_atoms(self):
            for struct in final_structures:
                yield from struct.get_atoms()
    io.set_structure(MultiStructure())
    io.save(output_path)
    total_time = time.time() - start_time
    print(f"\n💾 Resultado final salvo em: {output_path}")
    print(f"⏱️  Tempo total: {total_time:.1f}s")
