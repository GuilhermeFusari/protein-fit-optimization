import os
import zipfile
import shutil
from pathlib import Path

PASTA_ZIPS = Path("Zips_Proteinas")
PASTA_BASE_ORGANIZADA = Path("Dataset_Triado")

PASTA_PRONTAS = PASTA_BASE_ORGANIZADA / "1_Prontas_Para_Uso"
PASTA_ALPHAFOLD = PASTA_BASE_ORGANIZADA / "2_Fila_AlphaFold" # Falta PDB
PASTA_DAMMIN = PASTA_BASE_ORGANIZADA / "3_Fila_Dammin"       # Falta Envelope (CIF)
PASTA_LIXO = PASTA_BASE_ORGANIZADA / "4_Incompletos_Outros"  # Falta tudo

def criar_pastas():
    for pasta in [PASTA_PRONTAS, PASTA_ALPHAFOLD, PASTA_DAMMIN, PASTA_LIXO]:
        pasta.mkdir(parents=True, exist_ok=True)

def organizar_datasets():
    criar_pastas()

    arquivos_zip = list(PASTA_ZIPS.glob("*.zip"))
    print(f"Encontrados {len(arquivos_zip)} arquivos .zip para triagem.\n")

    for arquivo_zip in arquivos_zip:
        nome_proteina = arquivo_zip.stem
        pasta_temp = PASTA_BASE_ORGANIZADA / "temp" / nome_proteina

        try:
            with zipfile.ZipFile(arquivo_zip, 'r') as zip_ref:
                zip_ref.extractall(pasta_temp)
        except zipfile.BadZipFile:
            print(f"[{nome_proteina}] Arquivo ZIP corrompido ou inválido. Pulando.")
            continue

        pdbs = list(pasta_temp.rglob("*.pdb"))
        cifs = list(pasta_temp.rglob("*.cif"))

        tem_pdb = len(pdbs) > 0
        tem_cif = len(cifs) > 0

        if tem_pdb and tem_cif:
            destino = PASTA_PRONTAS / nome_proteina
            status = "PRONTA"
        elif tem_cif and not tem_pdb:
            destino = PASTA_ALPHAFOLD / nome_proteina
            status = "ALPHAFOLD (Falta PDB)"
        elif tem_pdb and not tem_cif:
            destino = PASTA_DAMMIN / nome_proteina
            status = "DAMMIN (Falta Envelope CIF)"
        else:
            destino = PASTA_LIXO / nome_proteina
            status = "INCOMPLETO (Falta PDB e CIF)"

        if destino.exists():
            shutil.rmtree(destino)
        shutil.move(str(pasta_temp), str(destino))

        print(f"[{nome_proteina}] -> {status}")

    pasta_temp_root = PASTA_BASE_ORGANIZADA / "temp"
    if pasta_temp_root.exists():
        shutil.rmtree(pasta_temp_root)

    print("\nTriagem finalizada com sucesso!")

if __name__ == "__main__":
    organizar_datasets()
