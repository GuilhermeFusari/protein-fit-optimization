# Protein Packing Optimization

Este projeto tem como objetivo **encaixar múltiplos arquivos PDB de proteínas dentro de um envelope estrutural** utilizando o algoritmo **Iterative Closest Point (ICP)** otimizado com **Chamfer Distance**, sem depender de RMSD.

## Estrutura do Projeto

protein-packing/
├─ data/ # Arquivos PDB e envelopes (não versionados)
├─ src/ # Scripts Python, ex: pdb_packing.py
├─ results/ # Resultados gerados (opcionalmente ignorados)
├─ .gitignore
└─ README.md

> Pastas `data/` e `results/` estão no `.gitignore` para não subir arquivos grandes.

## Funcionalidades

- Carrega arquivos PDB e o envelope estrutural.
- Converte PDBs em **nuvens de pontos**.
- Alinha proteínas dentro do envelope usando **ICP**.
- Calcula score de encaixe via **Chamfer Distance**.
- Seleciona automaticamente os melhores PDBs para o envelope.

## Requisitos

- Python 3.9+
- Bibliotecas: `numpy`, `scipy`, `open3d` (ou outra para point clouds, se usar)
- Outros pacotes padrão conforme `requirements.txt`

```bash
pip install -r requirements.txt

## Uso

Coloque seus arquivos PDB na pasta data/.

Ajuste parâmetros no script pdb_packing.py.

Execute: 
```bash
python scripts/pdb_packing.py


Os resultados serão salvos em results/.

