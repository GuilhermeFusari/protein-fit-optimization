import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from pathlib import Path

# Caminhos absolutos exatos
PASTA_42 = Path("/home/rvr.dias/Área de Trabalho/Guilherme_IC/resultados_finais_42")
PASTA_DAMMIN = Path("/home/rvr.dias/Área de Trabalho/Guilherme_IC/resultados_finais_dammin")
PASTA_GRAFICOS = Path("/home/rvr.dias/Área de Trabalho/Guilherme_IC/graficos_estatisticos")

def extrair_dados_estatisticos(pasta_base, label_dataset):
    dados = []
    if not pasta_base.exists():
        print(f"⚠️ Pasta não encontrada: {pasta_base}")
        return dados

    pastas_proteinas = [p for p in pasta_base.iterdir() if p.is_dir()]
    print(f"🔍 Vasculhando {len(pastas_proteinas)} proteínas em '{label_dataset}'...")

    for pasta_prot in pastas_proteinas:
        nome_proteina = pasta_prot.name
        
        pasta_result = pasta_prot / "packing_result"
        modo = "Packing"
        if not pasta_result.exists():
            pasta_result = pasta_prot / "single_result"
            modo = "Single"
            if not pasta_result.exists():
                continue

        chamfers_finais = []
        melhorias = []

        for rodada in range(1, 4):
            report_path = pasta_result / f"rodada_{rodada}" / "report.txt"
            if report_path.exists():
                with open(report_path, "r", encoding="utf-8") as f:
                    texto = f.read()
                    
                top1 = re.search(r"Top 1 (?:Error|Chamfer) \(Individual\):\s*([0-9.]+)", texto)
                global_err = re.search(r"Global (?:Error|Chamfer) \(Combined\):\s*([0-9.]+)", texto)
                
                t1_val = float(top1.group(1)) if top1 else None
                g_val = float(global_err.group(1)) if global_err else None
                
                c_final = g_val if modo == "Packing" and g_val else t1_val
                
                if c_final and c_final < 999:
                    chamfers_finais.append(c_final)
                    
                if modo == "Packing" and t1_val and g_val and t1_val > 0:
                    melhorias.append(((t1_val - g_val) / t1_val) * 100)

        if not chamfers_finais:
            continue

        media_chamfer = np.mean(chamfers_finais)
        std_chamfer = np.std(chamfers_finais) if len(chamfers_finais) > 1 else 0.0
        media_melhoria = np.mean(melhorias) if melhorias else None

        tamanho_aminoacidos = 0
        pdbs = list(pasta_result.rglob("*.pdb"))
        if pdbs:
            with open(pdbs[0], 'r') as f:
                for linha in f:
                    if linha.startswith("ATOM") and len(linha) > 16 and linha[12:16].strip() == "CA":
                        tamanho_aminoacidos += 1
            if tamanho_aminoacidos == 0:
                with open(pdbs[0], 'r') as f:
                    total_atoms = sum(1 for linha in f if linha.startswith("ATOM"))
                    tamanho_aminoacidos = max(1, total_atoms // 8)

        if tamanho_aminoacidos > 0:
            dados.append({
                "Proteina": nome_proteina,
                "Dataset": label_dataset,
                "Modo": modo,
                "Chamfer_Media": media_chamfer,
                "Chamfer_Std": std_chamfer,
                "Melhoria_Pct_Media": media_melhoria,
                "Tamanho_Aminoacidos": tamanho_aminoacidos,
                "Rodadas_Validas": len(chamfers_finais)
            })
            
    return dados

def gerar_relatorio_txt(df):
    caminho_txt = PASTA_GRAFICOS / "resumo_estatistico.txt"
    with open(caminho_txt, "w", encoding="utf-8") as f:
        f.write("="*60 + "\n")
        f.write("RELATÓRIO ESTATÍSTICO COMPLETO - ALINHAMENTO ICP\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"Total de Proteínas Analisadas com Sucesso: {len(df)}\n")
        f.write("-" * 60 + "\n\n")
        
        resumo = df.groupby(["Dataset", "Modo"])["Chamfer_Media"].describe()
        f.write(">>> ESTATÍSTICAS DA DISTÂNCIA DE CHAMFER (Å) <<<\n\n")
        f.write(resumo.to_string())
        f.write("\n\n" + "-" * 60 + "\n\n")
        
        df_packing = df[df["Modo"] == "Packing"]
        if not df_packing.empty:
            resumo_melhoria = df_packing.groupby("Dataset")["Melhoria_Pct_Media"].describe()
            f.write(">>> REDUÇÃO MÉDIA DA DISTÂNCIA DE CHAMFER NO MODO PACKING (%) <<<\n\n")
            f.write(resumo_melhoria.to_string())
            f.write("\n\n" + "-" * 60 + "\n")
            
    print(f"📄 Relatório de texto com os números exatos salvo em: {caminho_txt}")

def plotar_comparativo(df):
    if df.empty: return
        
    PASTA_GRAFICOS.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    
    label_validado = "Envelopes Validados (SASBDB)"
    label_abinitio = "Reconstrução ab initio (DAMMIF)"
    cores = {label_validado: "#1f77b4", label_abinitio: "#ff7f0e"}
    
    df_dammin = df[df["Dataset"] == label_abinitio]
    df_validado = df[df["Dataset"] == label_validado]
    
    max_x = df["Tamanho_Aminoacidos"].quantile(0.99)

    # 1. Dispersão: Chamfer vs Tamanho (Visão de Overlap dos Datasets)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.scatterplot(x="Tamanho_Aminoacidos", y="Chamfer_Media", data=df_dammin, 
                    color=cores[label_abinitio], s=30, alpha=0.3, label=label_abinitio, edgecolor=None, ax=ax)
    sns.scatterplot(x="Tamanho_Aminoacidos", y="Chamfer_Media", data=df_validado, 
                    color=cores[label_validado], s=100, alpha=1.0, label=label_validado, edgecolor="black", linewidth=1.5, ax=ax)
    ax.set_title("Precisão Média do Algoritmo vs Tamanho da Proteína")
    ax.set_xlabel("Tamanho da Proteína (Nº de Aminoácidos)")
    ax.set_ylabel("Média da Distância de Chamfer (Å)")
    ax.legend(title="Origem do Envelope")
    ax.set_xlim(0, max_x)
    plt.tight_layout()
    plt.savefig(PASTA_GRAFICOS / "1_Chamfer_Media_vs_Tamanho.png", dpi=300)
    plt.close()

    # 2. Dispersão: Melhoria vs Tamanho (Análise de Punição Assimétrica)
    df_packing = df[df["Modo"] == "Packing"]
    if not df_packing.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(x="Tamanho_Aminoacidos", y="Melhoria_Pct_Media", data=df_packing[df_packing["Dataset"] == label_abinitio], 
                        color=cores[label_abinitio], s=30, alpha=0.3, label=label_abinitio, edgecolor=None, ax=ax)
        sns.scatterplot(x="Tamanho_Aminoacidos", y="Melhoria_Pct_Media", data=df_packing[df_packing["Dataset"] == label_validado], 
                        color=cores[label_validado], s=100, alpha=1.0, label=label_validado, edgecolor="black", linewidth=1.5, ax=ax)
        ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label="Sem Melhoria (0%)")
        ax.set_title("Benefício Médio do Empacotamento Múltiplo vs Tamanho")
        ax.set_xlabel("Tamanho da Proteína (Nº de Aminoácidos)")
        ax.set_ylabel("Redução Média da Chamfer Distance (%)")
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        ax.set_xlim(0, max_x)
        plt.tight_layout()
        plt.savefig(PASTA_GRAFICOS / "2_Melhoria_Media_vs_Tamanho.png", dpi=300)
        plt.close()

    # 3. Histograma: Perfil Demográfico (Para os Métodos)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data=df, x="Tamanho_Aminoacidos", bins=40, kde=True, color="purple", alpha=0.6, ax=ax)
    ax.set_title(f"Perfil da Amostra: Frequência de Tamanhos (Total: {len(df)} Proteínas)")
    ax.set_xlabel("Tamanho da Proteína (Nº de Aminoácidos)")
    ax.set_ylabel("Quantidade de Proteínas")
    ax.set_xlim(0, max_x)
    plt.tight_layout()
    plt.savefig(PASTA_GRAFICOS / "3_Distribuicao_Tamanhos_Amostra.png", dpi=300)
    plt.close()

    # 4. Boxplot: Validação de Escalabilidade por Faixa de Peso
    fig, ax = plt.subplots(figsize=(10, 6))
    limites = [0, 200, 500, float('inf')]
    nomes_categorias = ['Pequenas\n(< 200 resíduos)', 'Médias\n(200 - 500 resíduos)', 'Grandes\n(> 500 resíduos)']
    df['Categoria_Tamanho'] = pd.cut(df['Tamanho_Aminoacidos'], bins=limites, labels=nomes_categorias)
    sns.boxplot(x="Categoria_Tamanho", y="Chamfer_Media", data=df, palette="viridis", width=0.5, ax=ax)
    ax.set_title(f"Impacto do Tamanho na Precisão do Encaixe")
    ax.set_xlabel("Categoria de Tamanho")
    ax.set_ylabel("Média da Distância de Chamfer Final (Å)")
    plt.tight_layout()
    plt.savefig(PASTA_GRAFICOS / "4_Desempenho_Categorias_Tamanho.png", dpi=300)
    plt.close()

    # 5. Violino: A Prova Definitiva (Single vs Packing)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(x="Modo", y="Chamfer_Media", hue="Dataset", data=df, split=True, 
                   inner="quartile", palette=cores, ax=ax, alpha=0.8)
    ax.set_title("Comparação de Desempenho: Single vs Packing (Geral)")
    ax.set_xlabel("Modo de Alinhamento")
    ax.set_ylabel("Média da Distância de Chamfer Final (Å)")
    plt.tight_layout()
    plt.savefig(PASTA_GRAFICOS / "5_Comparativo_Single_vs_Packing.png", dpi=300)
    plt.close()

    print(f"✅ Todos os 5 Gráficos estatísticos gerados com sucesso na pasta: {PASTA_GRAFICOS}")

if __name__ == "__main__":
    print("🚀 Iniciando extração dos testes de robustez (3 Rodadas)...")
    dados_42 = extrair_dados_estatisticos(PASTA_42, "Envelopes Validados (SASBDB)")
    dados_dammin = extrair_dados_estatisticos(PASTA_DAMMIN, "Reconstrução ab initio (DAMMIF)")
    
    df_completo = pd.DataFrame(dados_42 + dados_dammin)
    
    if not df_completo.empty:
        PASTA_GRAFICOS.mkdir(parents=True, exist_ok=True)
        caminho_csv = PASTA_GRAFICOS / "resultados_consolidados_estatisticos.csv"
        df_completo.to_csv(caminho_csv, index=False)
        print(f"📊 Planilha estatística consolidada salva em: {caminho_csv}")
        
        gerar_relatorio_txt(df_completo)
        plotar_comparativo(df_completo)
