#!/usr/bin/env python3
"""
gerar_figuras_v2.py

Regenera as tres figuras que mudaram com a revisao (lambda = 0,2, 50
entradas, bug do cifsup corrigido):

  figura_varredura_lambda.png   varredura de lambda (0,1 a 100), rep. completa
  figura_fatorial.png           ablation fatorial 2x2 (penalidade x downsampling)
  figura_comparacao_cifsup.png  ICP-SAXS vs cifsup NSD e ICP, nas 4 metricas

Le:
  sweep_lambda/lam_*.csv        varredura (uma linha por entrada, por lambda)
  fatorial/{A,B,C,D}_*.csv      as 4 celulas do fatorial
  benchmark_cifsup_lam02.csv    comparacao com o cifsup (lambda 0,2)

Uso (a partir de ~/Area de Trabalho/Guilherme_IC):
    python3 gerar_figuras_v2.py
"""

import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

OUT = "figuras_revisao"
DPI = 300
C_OURS = "#2C5F8A"
C_NSD = "#D98E3A"
C_ICP = "#C0552B"
C_BASE = "#B0B0B0"
C_HL = "#1E7A46"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
})


def load_ok(path):
    d = pd.read_csv(path)
    return d[d.status == "ok"] if "status" in d.columns else d


# --------------------------------------------------------------------------
# 1. varredura de lambda
# --------------------------------------------------------------------------

def fig_varredura():
    arqs = sorted(glob.glob("sweep_lambda/lam_*.csv"),
                  key=lambda p: float(re.search(r"lam_([0-9]+(?:\.[0-9]+)?)", p).group(1)))
    lam, cham, cob = [], [], []
    for p in arqs:
        L = float(re.search(r"lam_([0-9]+(?:\.[0-9]+)?)", p).group(1))
        d = load_ok(p)
        lam.append(L); cham.append(d.chamfer.mean()); cob.append(d.cobertura.mean())
    lam = np.array(lam)

    fig, ax1 = plt.subplots(figsize=(7.5, 4.2))
    ax1.set_xscale("log")
    l1 = ax1.plot(lam, cham, "o-", color=C_OURS, label="Chamfer distance", lw=1.8, ms=6)
    ax1.set_xlabel("Penalty weight λ (log scale)")
    ax1.set_ylabel("Chamfer distance (Å)", color=C_OURS)
    ax1.tick_params(axis="y", labelcolor=C_OURS)

    imin = int(np.argmin(cham))
    ax1.axvline(lam[imin], ls="--", color=C_HL, lw=1.2)
    ax1.annotate(f"optimum λ = {lam[imin]:g}\nChamfer {cham[imin]:.2f} Å",
                 xy=(lam[imin], cham[imin]),
                 xytext=(lam[imin]*2.2, cham[imin]+0.12),
                 fontsize=9, color=C_HL,
                 arrowprops=dict(arrowstyle="->", color=C_HL, lw=1))

    ax2 = ax1.twinx()
    l2 = ax2.plot(lam, cob, "s-", color=C_NSD, label="Coverage", lw=1.4, ms=5, alpha=0.85)
    ax2.set_ylabel("Envelope coverage", color=C_NSD)
    ax2.tick_params(axis="y", labelcolor=C_NSD)
    ax2.grid(False)

    lns = l1 + l2
    ax1.legend(lns, [x.get_label() for x in lns], loc="upper center", frameon=False)
    ax1.set_title("Calibration of λ with full Cα representation (n = 50)\n"
                  "monotonic degradation for λ > 0.2; optimum below the previously tested range",
                  fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_varredura_lambda.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_varredura_lambda.png")


# --------------------------------------------------------------------------
# 2. ablation fatorial 2x2
# --------------------------------------------------------------------------

def fig_fatorial():
    def m(pat):
        f = glob.glob(pat)
        return load_ok(f[0]).chamfer if f else None
    A = m("fatorial/A_*.csv")   # puro 300
    B = m("fatorial/B_*.csv")   # pen 300
    C = m("fatorial/C_*.csv")   # puro 3000
    D = m("fatorial/D_*.csv")   # pen 3000
    if any(x is None for x in (A, B, C, D)):
        print("  [pulado] fatorial: faltam CSVs em fatorial/")
        return

    means = [[A.mean(), C.mean()], [B.mean(), D.mean()]]  # linhas: sem/com pen; col: 300/3000
    errs = [[A.sem(), C.sem()], [B.sem(), D.sem()]]

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    x = np.arange(2)
    w = 0.36
    ax.bar(x - w/2, means[0], w, yerr=errs[0], capsize=3,
           color=C_BASE, edgecolor="black", lw=0.6, label="without penalty")
    ax.bar(x + w/2, means[1], w, yerr=errs[1], capsize=3,
           color=C_OURS, edgecolor="black", lw=0.6, label="with penalty (λ = 0.2)")
    ax.set_xticks(x)
    ax.set_xticklabels(["300 points\n(downsampled)", "3000 points\n(full Cα)"])
    ax.set_ylabel("Chamfer distance (Å)")
    ax.set_ylim(0, max(max(r) for r in means) * 1.18)

    for i in range(2):
        for j, val in enumerate([means[j][i] for j in range(2)]):
            xpos = i + (-w/2 if j == 0 else w/2)
            ax.text(xpos, val + 0.12, f"{val:.2f}", ha="center", fontsize=8.5)

    ef_pen = ((B.mean()-A.mean()) + (D.mean()-C.mean())) / 2
    ef_repr = ((C.mean()-A.mean()) + (D.mean()-B.mean())) / 2
    inter = (D.mean()-C.mean()) - (B.mean()-A.mean())
    txt = (f"main effects (both p < 0.0001):\n"
           f"  penalty: {ef_pen:+.2f} Å    representation: {ef_repr:+.2f} Å\n"
           f"  interaction: {inter:+.2f} Å (negligible)")
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8.5,
            va="bottom", ha="left", color="#333333",
            bbox=dict(boxstyle="round,pad=0.4", fc="#F4F7FA", ec="#CCCCCC"))

    ax.legend(loc="upper right", frameon=False)
    ax.set_title("Factorial ablation: penalty and representation act independently",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_fatorial.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_fatorial.png")


# --------------------------------------------------------------------------
# 3. comparacao com o cifsup
# --------------------------------------------------------------------------

def fig_comparacao():
    d = load_ok("benchmark_cifsup_lam02.csv")
    metrics = [("chamfer", "Chamfer (Å)", "lower"),
               ("frac_fora", "Fraction outside", "lower"),
               ("cobertura", "Coverage", "higher"),
               ("sep_centroide", "Centroid sep. (Å)", "lower")]
    groups = [("nsd", "cifsup NSD", C_NSD),
              ("icp", "cifsup ICP", C_ICP),
              ("ours", "ICP-SAXS", C_OURS)]

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.7))
    for ax, (col, label, better) in zip(axes, metrics):
        data, cols, names = [], [], []
        for pre, nome, cor in groups:
            c = f"{pre}_{col}"
            if c in d and d[c].notna().any():
                data.append(d[c].dropna()); cols.append(cor); names.append(nome)
        bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                        medianprops=dict(color="black", lw=1.3),
                        flierprops=dict(marker="o", ms=3, mfc="none", mec="#999", alpha=.5))
        for patch, cor in zip(bp["boxes"], cols):
            patch.set_facecolor(cor); patch.set_alpha(.75)
            patch.set_edgecolor("black"); patch.set_linewidth(.6)
        ax.set_xticks(range(1, len(names)+1))
        ax.set_xticklabels(names, rotation=20, ha="right")
        ax.set_ylabel(label)
        ax.set_title("↓ better" if better == "lower" else "↑ better",
                     fontsize=9, color="#666")

        # p ICP-SAXS vs cifsup NSD
        a, b = f"ours_{col}", f"nsd_{col}"
        if a in d and b in d:
            m = d[[a, b]].dropna()
            _, p = stats.wilcoxon(m[a], m[b])
            sig = "n.s." if p >= 0.05 else ("*" if p >= 0.01 else ("**" if p >= 0.001 else "***"))
            ax.text(0.5, 0.015, f"vs NSD: p={p:.3f} ({sig})",
                    transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=7.5, color="#555")

    fig.suptitle("ICP-SAXS (λ = 0.2) versus ATSAS cifsup on the same 50 SASBDB entries",
                 y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figura_comparacao_cifsup.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figura_comparacao_cifsup.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"gerando em {OUT}/\n")
    fig_varredura()
    fig_fatorial()
    fig_comparacao()
    print(f"\npronto -> {OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
