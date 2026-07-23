#!/usr/bin/env python3
"""
gerar_figura_packing.py

Controle do experimento de packing e figura correspondente.

Motivacao: um Chamfer menor com mais copias e esperado por construcao,
ja que mais pontos cobrem melhor o envelope. Sem um controle, a reducao
observada nao distingue "o alinhamento funcionou" de "acumulamos pontos".

O controle usa a MESMA estrutura e o MESMO numero de copias, mas com
orientacoes aleatorias e centros sorteados dentro do envelope. Se o
resultado real for muito melhor que o controle, o ganho vem do
alinhamento; se forem parecidos, vem da contagem de pontos.

Uso (a partir de ~/Area de Trabalho/Guilherme_IC):
    python3 gerar_figura_packing.py

Saidas em figures_paper/:
    figure_packing_control.png
    table_packing.csv / .tex
"""

import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

OUT = "figures_paper"
ROOT = "SAXS_Protein_Aligner"
N_CONTROL = 20          # repeticoes do controle aleatorio
SEED = 0

C_REAL = "#2C5F8A"
C_CTRL = "#B0B0B0"
C_SINGLE = "#7BA7C7"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
})

# O envelope de cada caso e dado EXPLICITAMENTE. Nao busque por padrao de
# nome: existem varios .cif por entrada no disco (damaver, damfilt,
# damstart, pddf) e apenas um corresponde ao sistema de coordenadas em que
# os RANK_*.pdb foram alinhados. Usar o errado produz numeros sem sentido
# que passam despercebidos, porque nada falha -- so o resultado fica ruim.
CASES = [
    ("GRB2 (30 copies)", "resultado_packing_30",
     "data/envelope.cif"),
    ("GRB2 (5 copies)", "resultado_packing",
     "data/envelope.cif"),
    ("RNA SASDFQ9 (20 models)", "resultado_packing_SASDFQ9",
     "data/3 TESTES/SASDFQ9/SASDFQ9/pddf/SASDFQ9-1.cif"),
]

# Chamfer combinado esperado, extraido dos report.txt originais. Serve de
# verificacao: se o recalculo divergir, o envelope usado esta errado.
EXPECTED = {
    "resultado_packing_30": 1.965,
    "resultado_packing": 2.161,
    "resultado_packing_SASDFQ9": 3.061,
}


def read_coords(path):
    pts = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    pts.append((float(line[30:38]), float(line[38:46]),
                                float(line[46:54])))
                except ValueError:
                    pass
    if pts:
        return np.asarray(pts, dtype=float)
    # fallback mmCIF
    cols, rows, in_loop = {}, [], False
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols[s.split(".", 1)[1].split()[0]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if s.startswith("#") or not s:
                    if rows:
                        break
                    continue
                if s.startswith("_"):
                    break
                rows.append(s.split())
    need = ("Cartn_x", "Cartn_y", "Cartn_z")
    if not all(k in cols for k in need):
        return np.empty((0, 3))
    ix, iy, iz = (cols[k] for k in need)
    for p in rows:
        if len(p) > max(ix, iy, iz):
            try:
                pts.append((float(p[ix]), float(p[iy]), float(p[iz])))
            except ValueError:
                pass
    return np.asarray(pts, dtype=float)


def chamfer(a, b):
    return float((cKDTree(b).query(a)[0].mean() +
                  cKDTree(a).query(b)[0].mean()) / 2.0)


def find_envelope(case_dir, hint):
    """
    Resolve o envelope apenas pelo caminho explicito do caso. Nao ha
    fallback por busca: escolher um .cif qualquer com nome parecido foi
    a origem de um resultado errado (o RNA foi comparado contra o
    envelope da GRB2 e depois contra outra reconstrucao, dando 10.99 A
    em vez de 3.06 A sem que nada acusasse erro).
    """
    if not hint:
        return None
    path = os.path.join(ROOT, hint)
    return path if os.path.exists(path) else None


def random_control(one_copy, env, n_copies, n_rep=N_CONTROL, seed=SEED):
    """
    Controle: n_copies da MESMA estrutura, cada uma com rotacao aleatoria
    e centro sorteado sobre um ponto do envelope. Repetido n_rep vezes.
    """
    rng = np.random.default_rng(seed)
    center = one_copy.mean(axis=0)
    scores = []
    for _ in range(n_rep):
        stack = []
        for _ in range(n_copies):
            M = Rotation.random(random_state=int(rng.integers(1 << 31))).as_matrix()
            pos = env[rng.integers(len(env))]
            stack.append((M @ (one_copy - center).T).T + pos)
        scores.append(chamfer(np.vstack(stack), env))
    return np.asarray(scores)


def main():
    if not os.path.isdir(ROOT):
        print(f"rode a partir da pasta que contem {ROOT}/")
        return 1
    os.makedirs(OUT, exist_ok=True)

    results = []
    for label, case_dir, env_hint in CASES:
        d = os.path.join(ROOT, case_dir)
        if not os.path.isdir(d):
            print(f"  [pulado] {case_dir} nao encontrado")
            continue
        files = sorted(glob.glob(os.path.join(d, "RANK_*.pdb")))
        if not files:
            print(f"  [pulado] sem RANK_*.pdb em {case_dir}")
            continue

        env_path = find_envelope(case_dir, env_hint)
        if not env_path:
            print(f"  [pulado] envelope nao encontrado para {case_dir}")
            print(f"           esperado em: {ROOT}/{env_hint}")
            continue
        env = read_coords(env_path)
        if len(env) < 10:
            print(f"  [pulado] envelope ilegivel: {env_path}")
            continue

        copies = [read_coords(f) for f in files]
        combined = np.vstack(copies)
        single = min(chamfer(c, env) for c in copies)
        real = chamfer(combined, env)
        ctrl = random_control(copies[0], env, len(files))

        exp = EXPECTED.get(case_dir)
        check = ""
        if exp is not None:
            if abs(real - exp) <= 0.05:
                check = f"  [confere com o report: {exp}]"
            else:
                check = (f"  [DIVERGE do report ({exp}) -- o envelope usado "
                         f"provavelmente nao e o do alinhamento]")

        results.append(dict(
            case=label, n_copies=len(files), n_env=len(env),
            single_best=round(single, 3),
            combined=round(real, 3),
            control_mean=round(float(ctrl.mean()), 3),
            control_std=round(float(ctrl.std(ddof=1)), 3),
            ratio=round(float(ctrl.mean()) / real, 2),
            report_value=exp,
        ))
        print(f"  {label}: single {single:.3f} | combined {real:.3f} | "
              f"control {ctrl.mean():.3f} +/- {ctrl.std(ddof=1):.3f}{check}")

    if not results:
        print("nenhum caso processado")
        return 1

    df = pd.DataFrame(results)
    df.to_csv(f"{OUT}/table_packing.csv", index=False)
    try:
        with open(f"{OUT}/table_packing.tex", "w") as f:
            f.write(df.to_latex(index=False, escape=True))
    except Exception:
        pass

    # ---- figura ----
    fig, ax = plt.subplots(figsize=(8.5, 4))
    x = np.arange(len(df))
    w = 0.26
    ax.bar(x - w, df.single_best, w, label="Best single copy",
           color=C_SINGLE, edgecolor="black", linewidth=0.6)
    ax.bar(x, df.combined, w, label="Multi-copy packing",
           color=C_REAL, edgecolor="black", linewidth=0.6)
    ax.bar(x + w, df.control_mean, w, yerr=df.control_std, capsize=3,
           label=f"Random control (n={N_CONTROL})",
           color=C_CTRL, edgecolor="black", linewidth=0.6)

    for i, r in df.iterrows():
        ax.text(i + w, r.control_mean + r.control_std + 0.15,
                f"{r.ratio:.1f}x", ha="center", fontsize=8, color="#444444")

    ax.set_xticks(x)
    ax.set_xticklabels(df.case, fontsize=9)
    ax.set_ylabel("Chamfer distance (Å)")
    ax.set_title("Multi-copy packing against a random-placement control\n"
                 "(lower is better; the control isolates the effect of point count)",
                 fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_packing_control.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"\n  {OUT}/figure_packing_control.png")
    print(f"  {OUT}/table_packing.csv / .tex")
    print("\nInterpretacao: se 'combined' for muito menor que 'control', o ganho")
    print("vem do alinhamento e nao da contagem de pontos. Note que isso mostra")
    print("que o envelope acomoda mais de um confomero rigido, o que e")
    print("consistente TANTO com oligomerizacao QUANTO com um ensemble")
    print("conformacional flexivel; distinguir os dois exige a massa molecular")
    print("estimada a partir de I(0).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
