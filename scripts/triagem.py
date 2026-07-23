import os
import zipfile
import shutil
from pathlib import Path

# --- CONFIGURAÇÕES DE PASTAS ---
# Caminho onde estão os seus zips (ajuste se necessário, por exemplo: r"C:\caminho\Zips_Proteinas")
PASTA_ZIPS = Path("Zips_Proteinas") 
PASTA_BASE_ORGANIZADA = Path("Dataset_Triado")

# Subpastas para a triagem
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
    print(f"🔍 Encontrados {len(arquivos_zip)} arquivos .zip para triagem.\n")

    for arquivo_zip in arquivos_zip:
        nome_proteina = arquivo_zip.stem # Pega o nome sem o .zip (ex: SASDFQ9)
        pasta_temp = PASTA_BASE_ORGANIZADA / "temp" / nome_proteina
        
        # 1. Extrai o conteúdo do Zip para uma pasta temporária
        try:
            with zipfile.ZipFile(arquivo_zip, 'r') as zip_ref:
                zip_ref.extractall(pasta_temp)
        except zipfile.BadZipFile:
            print(f"❌ [{nome_proteina}] Arquivo ZIP corrompido ou inválido. Pulando.")
            continue
            
        # 2. Conta os arquivos de interesse
        pdbs = list(pasta_temp.rglob("*.pdb"))
        cifs = list(pasta_temp.rglob("*.cif"))
        
        tem_pdb = len(pdbs) > 0
        tem_cif = len(cifs) > 0

        # 3. Lógica de Roteamento (Triagem)
        if tem_pdb and tem_cif:
            destino = PASTA_PRONTAS / nome_proteina
            status = "✅ PRONTA"
        elif tem_cif and not tem_pdb:
            destino = PASTA_ALPHAFOLD / nome_proteina
            status = "🤖 ALPHAFOLD (Falta PDB)"
        elif tem_pdb and not tem_cif:
            destino = PASTA_DAMMIN / nome_proteina
            status = "🧱 DAMMIN (Falta Envelope CIF)"
        else:
            destino = PASTA_LIXO / nome_proteina
            status = "❌ INCOMPLETO (Falta PDB e CIF)"

        # 4. Move a pasta temporária para o destino final correto
        if destino.exists():
            shutil.rmtree(destino) # Evita conflito se você precisar rodar o script duas vezes
        shutil.move(str(pasta_temp), str(destino))
        
        print(f"[{nome_proteina}] -> {status}")

    # Limpa a pasta temp vazia
    pasta_temp_root = PASTA_BASE_ORGANIZADA / "temp"
    if pasta_temp_root.exists():
        shutil.rmtree(pasta_temp_root)
        
    print("\n🚀 Triagem finalizada com sucesso!")

if __name__ == "__main__":
    organizar_datasets()
