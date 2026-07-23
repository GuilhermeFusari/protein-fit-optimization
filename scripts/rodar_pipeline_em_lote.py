import os
import subprocess
import time
import datetime
from pathlib import Path

def rodar_pipeline_em_lote(diretorio_base, diretorio_saida_base, caminho_main):
    caminho_base = Path(diretorio_base)
    caminho_saida = Path(diretorio_saida_base)
    
    if not caminho_base.exists() or not caminho_base.is_dir():
        print(f"❌ O diretório base '{diretorio_base}' não foi encontrado.")
        return

    caminho_saida.mkdir(parents=True, exist_ok=True)
    pastas = [p for p in caminho_base.iterdir() if p.is_dir()]
    total_pastas = len(pastas)
    print(f"🚀 Iniciando processamento de {total_pastas} pastas em {caminho_base.name} (3 Rodadas cada)...\n")

    # Inicia o cronômetro global do lote
    tempo_inicio = time.time()
    itens_processados = 0

    for pasta_proteina in pastas:
        arquivos_pdb = list(pasta_proteina.rglob("*.pdb"))
        arquivos_cif = list(pasta_proteina.rglob("*.cif"))
        
        if len(arquivos_pdb) == 0 or len(arquivos_cif) != 1:
            print(f"⏭️ Pulando '{pasta_proteina.name}': estrutura incorreta.")
            itens_processados += 1
            continue
        
        arquivo_cif = arquivos_cif[0]
        modo = "single" if len(arquivos_pdb) == 1 else "packing"
        pasta_saida_especifica = caminho_saida / pasta_proteina.name / f"{modo}_result"
        
        input_args = str(arquivos_pdb[0]) if modo == "single" else str(arquivos_pdb[0].parent)

        print(f"🔄 [{modo.upper()}] '{pasta_proteina.name}'")
        
        # O LOOP DA CIÊNCIA: Roda 3 vezes para a mesma proteína
        for rodada in range(1, 4):
            pasta_saida_rodada = pasta_saida_especifica / f"rodada_{rodada}"
            
            if pasta_saida_rodada.exists():
                print(f"   ⏭️ Rodada {rodada} já existe, pulando...")
                continue
                
            pasta_saida_rodada.mkdir(parents=True, exist_ok=True)
            
            comando = [
                "python3", caminho_main, 
                modo, 
                "-i", input_args, 
                "-e", str(arquivo_cif), 
                "-o", str(pasta_saida_rodada)
            ]
            
            try:
                subprocess.run(comando, check=True, stdout=subprocess.DEVNULL)
                print(f"   ✅ Rodada {rodada}/3 finalizada!")
            except subprocess.CalledProcessError as e:
                print(f"   ❌ Erro na rodada {rodada}. Código: {e.returncode}")

        # ==== CÁLCULO DE TEMPO E ETA ====
        itens_processados += 1
        tempo_decorrido = time.time() - tempo_inicio
        tempo_medio = tempo_decorrido / itens_processados
        tempo_restante = tempo_medio * (total_pastas - itens_processados)

        # Formatação bonita em HH:MM:SS
        eta_formatado = str(datetime.timedelta(seconds=int(tempo_restante)))
        decorrido_formatado = str(datetime.timedelta(seconds=int(tempo_decorrido)))

        print(f"   ⏳ Progresso: {itens_processados}/{total_pastas} | Decorrido: {decorrido_formatado} | Faltam aprox: {eta_formatado}\n")

if __name__ == "__main__":
    BASE_DIR = "/home/rvr.dias/Área de Trabalho/Guilherme_IC"
    CAMINHO_MAIN_PY = f"{BASE_DIR}/SAXS_Protein_Aligner/src/main.py"
    
    PASTA_PROTEINAS_42 = f"{BASE_DIR}/Dataset_Triado/1_Prontas_Para_Uso"
    PASTA_RESULTADOS_42 = f"{BASE_DIR}/resultados_finais_42"
    
    PASTA_PROTEINAS_DAMMIN = f"{BASE_DIR}/Dataset_Triado/3_Fila_Dammin"
    PASTA_RESULTADOS_DAMMIN = f"{BASE_DIR}/resultados_finais_dammin"
    
    print(f"🌟 INICIANDO PIPELINE COM TESTES DE ROBUSTEZ E CRONÔMETRO 🌟")
    print("="*50)
    rodar_pipeline_em_lote(PASTA_PROTEINAS_42, PASTA_RESULTADOS_42, CAMINHO_MAIN_PY)
    print("="*50)
    rodar_pipeline_em_lote(PASTA_PROTEINAS_DAMMIN, PASTA_RESULTADOS_DAMMIN, CAMINHO_MAIN_PY)
    print("🎉 PIPELINE TOTALMENTE CONCLUÍDO!")
