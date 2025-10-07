# Protein Fit Optimization

![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

Um sistema de alta performance para alinhamento e empacotamento de múltiplas estruturas de proteínas (PDB) dentro de envelopes de densidade eletrônica (CIF/MAP), utilizando o algoritmo *Iterative Closest Point (ICP)* otimizado com a métrica *Chamfer Distance*.

## 🎯 Sobre o Projeto

Em biologia estrutural, é comum ter um "envelope" de baixa resolução de uma macromolécula (obtido, por exemplo, por Microscopia Eletrônica ou SAXS) e um conjunto de estruturas de alta resolução de seus componentes. Este projeto automatiza o processo de encontrar a melhor combinação e orientação (fit) dessas estruturas menores dentro do envelope maior.

A ferramenta é ideal para:
* Modelagem de complexos proteicos.
* Validação de estruturas em dados de baixa resolução.
* Análise de encaixe molecular (docking) baseado em forma.

## 🚀 Começando

Siga estas instruções para obter uma cópia funcional do projeto em sua máquina local.

### Pré-requisitos

* Python 3.8 ou superior
* `pip` e `venv`

### Instalação

1.  **Clone o repositório:**
    ```sh
    git clone [https://github.com/seu-usuario/protein-packing.git](https://github.com/seu-usuario/protein-packing.git)
    cd protein-packing
    ```

2.  **Crie e ative um ambiente virtual:**
    ```sh
    python3 -m venv venv
    source venv/bin/activate
    # No Windows, use: venv\Scripts\activate
    ```

3.  **Instale as dependências:**
    ```sh
    pip install -r requirements.txt
    ```
    As principais dependências incluem `numpy`, `scipy` e `biopython`.

## ⚙️ Uso

O script principal `main.py` opera em dois modos: `single` para alinhamento de uma única proteína e `packing` para otimizar o encaixe de múltiplas proteínas.

### Modo 1: Alinhamento de Proteína Única (`single`)

Alinha uma única estrutura PDB a um envelope de referência. Útil para testes rápidos e alinhamentos individuais.

**Comando:**
```sh
python3 src/main.py single -i data/protein.pdb -e data/envelope.cif -o results/single_output
```

### Modo 2: Empacotamento de Múltiplas Estruturas (`packing`)

Testa um conjunto de PDBs, seleciona as que melhor se encaixam no envelope sem sobreposição e gera um modelo combinado.

**Comando:**
```sh
python3 src/main.py packing -i data/pdbs/ -e data/envelope.cif -o results/packing_output
```
O resultado final será salvo como `results/packing_output/packing_result.pdb`.

### Parâmetros Opcionais

Os parâmetros abaixo podem ser usados em ambos os modos para ajustar a precisão e o desempenho do algoritmo.

| Parâmetro            | Descrição                                                                      | Padrão |
| -------------------- | ------------------------------------------------------------------------------ | ------ |
| `--max-iter`         | Número máximo de iterações para o algoritmo ICP.                               | `50`   |
| `--max-structures`   | (Apenas `packing`) Número máximo de PDBs a serem testados do diretório de entrada. | `20`   |
| `--n-rotations`      | Número de rotações iniciais a testar para encontrar o melhor alinhamento global. | `12`   |

## 🛠️ Metodologia

O núcleo do alinhamento é baseado em dois conceitos principais:

* **Iterative Closest Point (ICP):** Um algoritmo que alinha duas nuvens de pontos (neste caso, os átomos da proteína e os pontos do envelope) minimizando iterativamente a distância entre eles.
* **Chamfer Distance:** Uma métrica de dissimilaridade entre duas nuvens de pontos. Ela mede a distância média de cada ponto em um conjunto ao seu vizinho mais próximo no outro conjunto. Foi escolhida por ser robusta a ruído e eficaz para comparar formas complexas.

## 📁 Estrutura do Projeto

```
protein-packing/
├─ data/             # (Não versionado) Armazena arquivos PDB de entrada e envelopes.
├─ src/              # Código fonte Python.
│  ├─ main.py        # Orquestrador principal (executa os modos 'single' e 'packing').
│  ├─ icp_aligner.py # Implementação do alinhamento ICP + Chamfer Distance.
│  └─ pdb_packing.py # Lógica para seleção e empacotamento de múltiplas proteínas.
├─ results/          # (Não versionado) Diretório para salvar os resultados.
├─ .gitignore
└─ README.md
```

Veja os [issues abertos](https://github.com/seu-usuario/protein-packing/issues) para uma lista completa de funcionalidades propostas (e problemas conhecidos).

## 🤝 Contribuição

Contribuições são o que tornam a comunidade de código aberto um lugar incrível para aprender, inspirar e criar. Qualquer contribuição que você fizer será **muito apreciada**.

1.  Faça um *Fork* do projeto
2.  Crie sua *Feature Branch* (`git checkout -b feature/AmazingFeature`)
3.  Faça o *Commit* de suas alterações (`git commit -m 'Add some AmazingFeature'`)
4.  Faça o *Push* para a *Branch* (`git push origin feature/AmazingFeature`)
5.  Abra um *Pull Request*

## 📜 Licença

Distribuído sob a licença MIT. Veja `LICENSE.txt` para mais informações.

---
