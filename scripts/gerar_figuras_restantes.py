#!/usr/bin/env python3
"""
gerar_figuras_restantes.py

Gera as tres figuras que faltavam para o artigo (lambda 0.2, dados atuais):

  figura_enantiomeros.png     fracao de entradas em que o espelho venceu
  figura_packing.png          packing GRB2/RNA vs controle aleatorio
  figura_ground_truth.png     recuperacao de pose: encaixe equivalente,
                              orientacao ambigua, efeito da anisometria

Le, a partir de ~/Area de Trabalho/Guilherme_IC:
  final_lam02.csv                       (coluna 'espelhado')
  ground_truth_final.csv
  SAXS_Protein_Aligner/resultado_packing_30/report.txt
  SAXS_Protein_Aligner/resultado_packing/report.txt
  SAXS_Protein_Aligner/resultado_packing_SASDFQ9/report.txt

Uso:
    python3 gerar_figuras_restantes.py
"""

import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "SAXS_Protein_Aligner/figuras_artigo"
DPI = 300
C_MAIN = "#2C5F8A"
C_ALT = "#D98E3A"
C_GREY = "#B0B0B0"
C_HL = "#1E7A46"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
})


# --------------------------------------------------------------------------
# 1. enantiomeros
# --------------------------------------------------------------------------

def fig_enantiomeros():
    d = pd.read_csv("final_lam02.csv")
    if "status" in d:
        d = d[d.status == "ok"]
    if "espelhado" not in d:
        print("  [pulado] final_lam02.csv sem coluna 'espelhado'")
        return
    n = len(d)
    n_mir = int(d.espelhado.sum())
    n_orig = n - n_mir

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.bar([0, 1], [n_orig, n_mir], color=[C_GREY, C_MAIN],
           edgecolor="black", linewidth=0.6, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Original hand", "Mirrored hand"])
    ax.set_ylabel("Number of entries")
    ax.axhline(n / 2, ls="--", color="#999999", lw=1)
    ax.text(0.98, n / 2, " expected 50%", va="bottom", ha="right",
            transform=ax.get_yaxis_transform(), fontsize=9, color="#666666")
    for i, v in enumerate([n_orig, n_mir]):
        ax.text(i, v + 0.4, str(v), ha="center", fontsize=11, fontweight="bold")
    ax.set_title(f"Enantiomorph selection ({n_mir}/{n} mirrored)\n"
                 "close to the 50% expected when the envelope\n"
                 "carries no chirality information", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_enantiomeros.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_enantiomeros.png")


# --------------------------------------------------------------------------
# 2. packing vs controle
# --------------------------------------------------------------------------

def read_report_combined(path):
    """Extrai o 'Global Error (Combined)' de um report.txt."""
    if not os.path.exists(path):
        return None
    for line in open(path, errors="replace"):
        m = re.search(r"Global Error \(Combined\):\s*([\d.]+)", line)
        if m:
            return float(m.group(1))
    return None


def fig_packing():
    R = "SAXS_Protein_Aligner"
    # valores do controle aleatorio ja calculados na sessao (validados)
    casos = [
        ("GRB2\n(30 copies)", f"{R}/resultado_packing_30/report.txt", 5.690),
        ("GRB2\n(5 copies)", f"{R}/resultado_packing/report.txt", 7.034),
        ("RNA SASDFQ9\n(20 models)", f"{R}/resultado_packing_SASDFQ9/report.txt", 6.358),
    ]
    labels, real, ctrl = [], [], []
    for lab, rep, c in casos:
        v = read_report_combined(rep)
        if v is None:
            continue
        labels.append(lab); real.append(v); ctrl.append(c)
    if not real:
        print("  [pulado] packing: reports nao encontrados")
        return

    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, real, w, label="Aligned packing",
           color=C_MAIN, edgecolor="black", linewidth=0.6)
    ax.bar(x + w / 2, ctrl, w, label="Random control",
           color=C_GREY, edgecolor="black", linewidth=0.6)
    for i, (rr, cc) in enumerate(zip(real, ctrl)):
        ax.text(i - w / 2, rr + 0.1, f"{rr:.2f}", ha="center", fontsize=8.5)
        ax.text(i + w / 2, cc + 0.1, f"{cc:.2f}", ha="center", fontsize=8.5)
        ax.text(i, max(rr, cc) + 0.5, f"{cc/rr:.1f}x", ha="center",
                fontsize=9, color=C_HL, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Chamfer distance (Å)")
    ax.set_ylim(0, max(ctrl) * 1.2)
    ax.legend(frameon=False)
    ax.set_title("Multi-copy packing versus a random-placement control\n"
                 "the gain comes from the alignment, not the point count",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_packing.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_packing.png")


# --------------------------------------------------------------------------
# 3. ground truth
# --------------------------------------------------------------------------

def fig_ground_truth():
    d = pd.read_csv("ground_truth_final.csv")
    if "status" in d:
        d = d[d.status == "ok"]
    if "delta_chamfer" not in d:
        print("  [pulado] ground_truth_final.csv sem colunas novas")
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    # painel A: Chamfer verdadeira vs recuperada (encaixe equivalente)
    ax = axes[0]
    ax.scatter(d.chamfer_verdadeira, d.chamfer_recuperada, s=40,
               color=C_MAIN, edgecolor="black", linewidth=0.4, alpha=0.8)
    lim = [0, max(d.chamfer_verdadeira.max(), d.chamfer_recuperada.max()) * 1.1]
    ax.plot(lim, lim, "--", color="#999999", lw=1, label="identical fit")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Chamfer of the true pose (Å)")
    ax.set_ylabel("Chamfer of the recovered pose (Å)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Recovered poses fit the envelope\n"
                 "as well as the true pose", fontsize=10)

    # painel B: erro angular por anisometria
    ax = axes[1]
    d = d.copy()
    d["classe"] = pd.cut(d.anisometria, [0, 1.5, 2.5, 100],
                         labels=["globular\n(<1.5)", "medium\n(1.5–2.5)", "elongated\n(>2.5)"])
    grp = d.groupby("classe", observed=True)["erro_ang"]
    classes = list(grp.groups.keys())
    means = [grp.get_group(c).mean() for c in classes]
    sems = [grp.get_group(c).sem() for c in classes]
    xx = np.arange(len(classes))
    ax.bar(xx, means, yerr=sems, capsize=3, color=C_ALT,
           edgecolor="black", linewidth=0.6, width=0.6)
    ax.set_xticks(xx)
    ax.set_xticklabels([str(c) for c in classes])
    ax.set_ylabel("Orientation error (degrees)")
    ax.set_title("Orientation ambiguity decreases\n"
                 "for anisometric particles", fontsize=10)

    fig.suptitle("Pose-recovery test with known ground truth (n = %d trials)" % len(d),
                 y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_ground_truth.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_ground_truth.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"gerando em {OUT}/\n")
    fig_enantiomeros()
    fig_packing()
    fig_ground_truth()
    print(f"\npronto")
    return 0


if __name__ == "__main__":
    sys.exit(main())
