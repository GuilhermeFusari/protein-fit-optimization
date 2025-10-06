# Protein Fit Optimization

Sistema para alinhamento e empacotamento de proteínas PDB dentro de envelopes estruturais usando ICP com Chamfer Distance.

## Estrutura do Projeto
```
protein-packing/
├─ data/ (Arquivos PDB e envelopes, não versionados)
├─ src/ (Scripts Python)
│ ├─ main.py (Controlador principal)
│ ├─ icp_aligner.py (Alinhamento individual usando ICP + Chamfer)
│ └─ pdb_packing.py (Pipeline de seleção e empacotamento)
├─ results/ (Resultados gerados, ignorados pelo git)
├─ .gitignore
└─ README.md
```
## Instalação

Recomenda-se criar um ambiente virtual Python e instalar as dependências:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
Dependências principais: numpy, scipy, biopython

## Uso

O script principal é main.py. Ele suporta dois modos de execução:

1. Modo single (alinhamento de uma proteína)

Alinha um único PDB a um envelope CIF:
```bash
python3 src/main.py single -i data/protein.pdb -e data/envelope.cif -o results/single_output
```
Parâmetros opcionais:

--max-iter: máximo de iterações ICP (default 50)

--n-rotations: rotações iniciais a testar (modo packing)

2. Modo packing (seleção e empacotamento de múltiplos PDBs)

Seleciona os melhores PDBs que se encaixam no envelope e gera um PDB combinado:
```bash
python3 src/main.py packing -i data/pdbs/ -e data/envelope.cif -o results/packing_output
```
Parâmetros opcionais:

--max-iter: máximo de iterações ICP (default 50)

--max-structures: número máximo de PDBs a testar (default 20)

--n-rotations: rotações iniciais a testar (default 12)

O resultado final será salvo em:
```bash
results/packing_output/packing_result.pdb
```
Métrica

O alinhamento é avaliado usando Chamfer Distance, que mede a similaridade entre a forma da proteína e do envelope.

Observações

A pasta data/ não é versionada para evitar subir arquivos grandes.

Use sempre Python 3.8+ e um ambiente virtual para manter dependências isoladas.
