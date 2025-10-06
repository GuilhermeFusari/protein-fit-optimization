#!/usr/bin/env python3
"""
MAIN - Controlador dos modos 'single' e 'packing'
"""

"""
Este script `main.py` funciona como o controlador principal do programa de alinhamento e empacotamento de PDBs.

 Argumentos e configuração do programa
- Usa `argparse` para definir e ler os argumentos da linha de comando.
- O programa pode ser executado em dois modos:
    1. 'single': alinha um único arquivo PDB a um envelope CIF.
    2. 'packing': seleciona e empacota múltiplos PDBs no envelope.
- Parâmetros adicionais incluem número máximo de iterações ICP, quantidade máxima de estruturas
  para o modo packing e número de rotações iniciais a testar.

Execução dos modos
- Modo 'single': chama o script `icp_aligner.py` via `subprocess`, passando os caminhos de entrada,
  envelope e saída. Todo o processamento é feito pelo script `icp_aligner.py`.
- Modo 'packing': chama a função `run_packing_pipeline` do módulo `pdb_packing.py`, que lida com
  múltiplos PDBs, calcula Chamfer Distance e seleciona os melhores para empacotamento.
- Ambos os modos têm tratamento de exceções para capturar erros durante a execução e fornecer
  mensagens claras de falha, encerrando o programa com `sys.exit(1)`.
"""

import argparse
import subprocess
import sys
import traceback
from pathlib import Path

from pdb_packing import run_packing_pipeline


def parse_arguments():
    parser = argparse.ArgumentParser(description='Programa de alinhamento e empacotamento de PDBs')
    parser.add_argument('mode', choices=['single', 'packing'],
                        help="Modo de execução: 'single' (um PDB e um envelope) ou 'packing' (vários PDBs).")
    parser.add_argument('-i', '--input', required=True, help='Pasta com os arquivos PDB ou o arquivo único.')
    parser.add_argument('-e', '--envelope', required=True, help='Arquivo CIF do envelope.')
    parser.add_argument('-o', '--output', required=True, help='Pasta de saída dos resultados.')
    parser.add_argument('--max-iter', type=int, default=50, help='Número máximo de iterações ICP.')
    parser.add_argument('--max-structures', type=int, default=20, help='Número máximo de estruturas (packing).')
    parser.add_argument('--n-rotations', type=int, default=12, help='Número de rotações iniciais (packing).')
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.mode == 'single':
        try:
            subprocess.run(
                ['python3', 'icp_aligner.py',
                 '-i', args.input,
                 '-e', args.envelope,
                 '-o', args.output],
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"\n💥 Erro ao executar modo single: {e}")
            sys.exit(1)

    elif args.mode == 'packing':
        try:
            run_packing_pipeline(args.input, args.envelope, args.output, args)
        except Exception as e:
            print(f"\n💥 Erro no modo packing: {e}")
            traceback.print_exc()
            sys.exit(1)


if __name__ == "__main__":
    main()
