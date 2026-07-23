#!/usr/bin/env python3
"""
extrai_reports.py

Varre resultados_finais_dammin/ e extrai os scores de todos os report.txt.

Uso:
    python3 extrai_reports.py ~/"Area de Trabalho/Guilherme_IC/resultados_finais_dammin"

Saida:
    reports_extraidos.csv   uma linha por report.txt encontrado

Formatos reconhecidos:
    SINGLE  -> Top 1 Error, Best PDB File, Execution Time, Envelope Size, PDBs Evaluated
    PACKING -> Top 1 Error, Global Error, Improvement, Number of Copies,
               Penalty Weight, Execution Time, Envelope Size, tabela RANK
"""

import csv
import re
import sys
from pathlib import Path

# "Rotulo:  <numero>"  -> captura o numero, ignorando unidade
NUM = r"([-+]?\d+(?:\.\d+)?)"
PATTERNS = {
    "top1_chamfer":   re.compile(r"Top\s*1\s*Error[^:]*:\s*" + NUM, re.I),
    "global_chamfer": re.compile(r"Global\s*Error[^:]*:\s*" + NUM, re.I),
    "n_copias":       re.compile(r"Number\s*of\s*Copies\s*:\s*" + NUM, re.I),
    "lambda":         re.compile(r"Penalty\s*Weight[^:]*:\s*" + NUM, re.I),
    "tempo_s":        re.compile(r"Execution\s*Time\s*:\s*" + NUM, re.I),
    "envelope_pts":   re.compile(r"Envelope\s*Size[^:]*:\s*" + NUM, re.I),
    "pdbs_avaliados": re.compile(r"PDBs\s*Evaluated\s*:\s*" + NUM, re.I),
}
IMPROVE = re.compile(r"Improvement\s*:\s*(\w+)", re.I)
BESTPDB = re.compile(r"Best\s*PDB\s*File\s*:\s*(\S+)", re.I)
# linhas "#1    | 17.5435       | arquivo.pdb"
RANK = re.compile(r"^#(\d+)\s*\|\s*" + NUM + r"\s*\|\s*(\S+)", re.M)


def parse(path):
    txt = path.read_text(errors="replace")
    r = {}
    for key, pat in PATTERNS.items():
        m = pat.search(txt)
        r[key] = float(m.group(1)) if m else None

    m = IMPROVE.search(txt)
    r["improvement"] = m.group(1).upper() if m else None
    m = BESTPDB.search(txt)
    r["best_pdb"] = m.group(1) if m else None

    ranks = RANK.findall(txt)
    r["n_ranks"] = len(ranks) if ranks else None
    if ranks:
        # melhor chamfer da tabela (menor valor)
        r["melhor_rank_chamfer"] = min(float(v) for _, v, _ in ranks)
    else:
        r["melhor_rank_chamfer"] = None

    # modo: pelo cabecalho, com fallback no caminho
    if re.search(r"SINGLE\s*MODE", txt, re.I):
        r["modo"] = "single"
    elif re.search(r"PACKING", txt, re.I) or r["global_chamfer"] is not None:
        r["modo"] = "packing"
    else:
        r["modo"] = "single" if "single_result" in str(path) else "packing"
    return r


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1

    root = Path(sys.argv[1]).expanduser()
    if not root.is_dir():
        print(f"diretorio nao encontrado: {root}")
        return 1

    reports = sorted(root.rglob("report.txt"))
    print(f"{len(reports)} arquivos report.txt encontrados")
    if not reports:
        return 1

    rows, erros = [], 0
    for i, p in enumerate(reports, 1):
        if i % 500 == 0:
            print(f"  {i}/{len(reports)}")
        try:
            rel = p.relative_to(root)
            partes = rel.parts
            row = {
                "entry": partes[0] if partes else "?",
                "rodada": next((x for x in partes if x.startswith("rodada")), ""),
                "caminho": str(rel),
            }
            row.update(parse(p))
            rows.append(row)
        except Exception as e:
            erros += 1
            if erros <= 5:
                print(f"  [!] falha em {p}: {e}")

    cols = ["entry", "modo", "rodada", "top1_chamfer", "global_chamfer",
            "melhor_rank_chamfer", "improvement", "n_copias", "lambda",
            "tempo_s", "envelope_pts", "pdbs_avaliados", "n_ranks",
            "best_pdb", "caminho"]
    out = "reports_extraidos.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    entradas = len({r["entry"] for r in rows})
    n_single = sum(1 for r in rows if r["modo"] == "single")
    n_pack = sum(1 for r in rows if r["modo"] == "packing")
    com_score = sum(1 for r in rows if r["top1_chamfer"] is not None)

    print(f"\n{len(rows)} reports | {entradas} entradas unicas")
    print(f"  single: {n_single}   packing: {n_pack}")
    print(f"  com chamfer extraido: {com_score}/{len(rows)}")
    if erros:
        print(f"  falhas de leitura: {erros}")
    print(f"\nsalvo: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
