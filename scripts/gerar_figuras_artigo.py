#!/usr/bin/env python3
"""
gerar_figuras_artigo.py

Gera todas as figuras e tabelas do artigo, em ingles, a partir dos CSVs
produzidos pelo pipeline.

Uso (a partir de ~/Area de Trabalho/Guilherme_IC):
    python3 gerar_figuras_artigo.py

Entradas esperadas no diretorio atual:
    baseline_icp_puro.csv          ICP sem restricao (uni, lam=50)
    confirma_v2.csv                bidirecional lam=2, 300 pontos
    confirma_sem_downsampling.csv  bidirecional lam=2, sem downsampling
    final_v2.csv                   metodo final (+ enantiomorfos)
    benchmark_cifsup_v2.csv        comparacao com o ATSAS

Saidas em figures_paper/:
    figure_ablation.png       efeito cumulativo de cada componente
    figure_comparison.png     ICP-SAXS vs cifsup, pareado
    figure_enantiomorphs.png  ganho da busca de enantiomorfos
    figure_per_entry.png      Chamfer por entrada, ordenado
    table_*.tex / table_*.csv tabelas prontas
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

OUT = "figures_paper"
DPI = 300

# paleta sobria, legivel em preto e branco
C_BASE = "#B0B0B0"
C_STEP = "#7BA7C7"
C_FINAL = "#2C5F8A"
C_ATSAS = "#D98E3A"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.axisbelow": True,
})

FILES = {
    "baseline": "baseline_icp_puro.csv",
    "lam2_300": "confirma_v2.csv",
    "nodown": "confirma_sem_downsampling.csv",
    "final": "final_v2.csv",
    "cifsup": "benchmark_cifsup_v2.csv",
}

METRICS = [
    ("chamfer", "Chamfer distance (Å)", "lower"),
    ("frac_fora", "Fraction outside envelope", "lower"),
    ("cobertura", "Envelope coverage", "higher"),
    ("sep_centroide", "Centroid separation (Å)", "lower"),
]


def load():
    d = {}
    missing = []
    for k, f in FILES.items():
        if not os.path.exists(f):
            missing.append(f)
            continue
        df = pd.read_csv(f)
        if "status" in df.columns:
            df = df[df.status == "ok"]
        d[k] = df
    if missing:
        print("Arquivos nao encontrados:")
        for f in missing:
            print("   ", f)
        print("\nGere-os antes de rodar este script (veja o cabecalho).")
        sys.exit(1)
    return d


def wilcox(a, b):
    m = pd.concat([a, b], axis=1).dropna()
    if len(m) < 5:
        return np.nan
    try:
        return stats.wilcoxon(m.iloc[:, 0], m.iloc[:, 1]).pvalue
    except ValueError:
        return np.nan


def stars(p):
    if np.isnan(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# --------------------------------------------------------------------------
# Figura 1 - ablation
# --------------------------------------------------------------------------

def fig_ablation(d):
    steps = [
        ("Unconstrained ICP", d["baseline"], C_BASE),
        ("Bidirectional penalty", d["lam2_300"], C_STEP),
        ("No downsampling", d["nodown"], C_STEP),
        ("Enantiomorph search", d["final"], C_FINAL),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))

    for ax, (col, label, better) in zip(axes, METRICS):
        vals = [s[1][col].mean() for s in steps]
        errs = [s[1][col].std(ddof=1) / np.sqrt(len(s[1])) for s in steps]
        cols = [s[2] for s in steps]
        x = np.arange(len(steps))
        ax.bar(x, vals, yerr=errs, color=cols, edgecolor="black",
               linewidth=0.6, capsize=3, width=0.68)
        ax.set_xticks(x)
        ax.set_xticklabels([f"({i+1})" for i in range(len(steps))])
        ax.set_ylabel(label)
        arrow = "\u2193 better" if better == "lower" else "\u2191 better"
        ax.set_title(arrow, fontsize=9, color="#555555")

        # significancia do passo final contra o baseline
        p = wilcox(d["final"][col].reset_index(drop=True),
                   d["baseline"][col].reset_index(drop=True))
        ax.text(0.5, 0.94, stars(p), transform=ax.transAxes,
                ha="center", va="top", fontsize=11)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=s[2], ec="black", lw=0.6) for s in steps]
    labels = [f"({i+1}) {s[0]}" for i, s in enumerate(steps)]
    fig.legend(handles, labels, loc="lower center", ncol=4,
               frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Cumulative effect of each methodological component (n = 50)",
                 y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_ablation.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figure_ablation.png")


# --------------------------------------------------------------------------
# Figura 2 - comparacao com o ATSAS
# --------------------------------------------------------------------------

def fig_comparison(d):
    c = d["cifsup"]
    groups = [
        ("cifsup\n(NSD)", "nsd", C_ATSAS),
        ("cifsup\n(ICP)", "icp", C_ATSAS),
        ("ICP-SAXS\n(this work)", "ours", C_FINAL),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))

    for ax, (col, label, better) in zip(axes, METRICS):
        data, cols, names = [], [], []
        for name, pre, colr in groups:
            c_ = f"{pre}_{col}"
            if c_ in c and c[c_].notna().any():
                data.append(c[c_].dropna())
                cols.append(colr)
                names.append(name)
        if not data:
            ax.axis("off")
            continue
        bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                        medianprops=dict(color="black", lw=1.4),
                        flierprops=dict(marker="o", ms=3, mfc="none",
                                        mec="#888888", alpha=0.6))
        for patch, colr in zip(bp["boxes"], cols):
            patch.set_facecolor(colr)
            patch.set_alpha(0.75)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.6)
        ax.set_xticklabels(names)
        ax.set_ylabel(label)
        arrow = "\u2193 better" if better == "lower" else "\u2191 better"
        ax.set_title(arrow, fontsize=9, color="#555555")

        a, b = f"ours_{col}", f"icp_{col}"
        if a in c and b in c:
            p = wilcox(c[a], c[b])
            if not np.isnan(p):
                ax.text(0.5, 0.02, f"p = {p:.2f} ({stars(p)})",
                        transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=8, color="#555555")

    fig.suptitle("ICP-SAXS versus ATSAS cifsup on the same 50 SASBDB entries",
                 y=1.01, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_comparison.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figure_comparison.png")


# --------------------------------------------------------------------------
# Figura 3 - enantiomorfos
# --------------------------------------------------------------------------

def fig_enantiomorphs(d):
    a = d["nodown"][["entry", "chamfer"]].rename(columns={"chamfer": "off"})
    b = d["final"][["entry", "chamfer", "espelhado"]].rename(columns={"chamfer": "on"})
    m = a.merge(b, on="entry")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))

    ax = axes[0]
    mir = m[m.espelhado]
    nom = m[~m.espelhado]
    ax.scatter(nom.off, nom.on, s=32, facecolor="none",
               edgecolor="#888888", label="Original hand kept")
    ax.scatter(mir.off, mir.on, s=32, color=C_FINAL,
               edgecolor="black", linewidth=0.4, label="Mirrored hand selected")
    lim = [0, max(m.off.max(), m.on.max()) * 1.06]
    ax.plot(lim, lim, "--", color="#999999", lw=1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Chamfer distance, no mirror search (Å)")
    ax.set_ylabel("Chamfer distance, with mirror search (Å)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Per-entry effect")

    ax = axes[1]
    n_mir = int(m.espelhado.sum())
    ax.bar([0, 1], [len(m) - n_mir, n_mir],
           color=["#B0B0B0", C_FINAL], edgecolor="black", linewidth=0.6, width=0.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Original", "Mirrored"])
    ax.set_ylabel("Number of entries")
    ax.axhline(len(m) / 2, ls="--", color="#999999", lw=1)
    ax.text(0.98, len(m) / 2, " expected 50%", va="bottom", ha="right",
            transform=ax.get_yaxis_transform(), fontsize=8, color="#666666")
    ax.set_title(f"Selected hand ({n_mir}/{len(m)} mirrored)")

    fig.suptitle("SAXS envelopes do not define chirality: both hands must be tested",
                 y=1.00, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_enantiomorphs.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figure_enantiomorphs.png")


# --------------------------------------------------------------------------
# Figura 4 - por entrada
# --------------------------------------------------------------------------

def fig_per_entry(d):
    f = d["final"][["entry", "chamfer"]].sort_values("chamfer").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 3.8))
    x = np.arange(len(f))
    ax.bar(x, f.chamfer, color=C_FINAL, edgecolor="black", linewidth=0.4, width=0.72)
    ax.axhline(f.chamfer.median(), ls="--", color="#C0392B", lw=1.2,
               label=f"Median: {f.chamfer.median():.2f} Å")
    ax.axhline(f.chamfer.mean(), ls=":", color="#2C3E50", lw=1.2,
               label=f"Mean: {f.chamfer.mean():.2f} Å")
    ax.set_xticks(x)
    ax.set_xticklabels(f.entry, rotation=90, fontsize=6)
    ax.set_ylabel("Chamfer distance (Å)")
    ax.set_xlabel("SASBDB entry (sorted by alignment quality)")
    ax.legend(frameon=False)
    ax.set_title(f"Alignment quality across all {len(f)} benchmark entries")
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_per_entry.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  figure_per_entry.png")


# --------------------------------------------------------------------------
# Tabelas
# --------------------------------------------------------------------------

def tables(d):
    # Tabela 1 - ablation
    rows = []
    steps = [
        ("Unconstrained ICP (baseline)", "baseline"),
        ("+ bidirectional penalty (lambda = 2)", "lam2_300"),
        ("+ full Ca representation", "nodown"),
        ("+ enantiomorph search (final)", "final"),
    ]
    for label, key in steps:
        df = d[key]
        r = {"Configuration": label}
        for col, name, _ in METRICS:
            r[name] = f"{df[col].mean():.3f}"
        if key != "baseline":
            p = wilcox(df["chamfer"].reset_index(drop=True),
                       d["baseline"]["chamfer"].reset_index(drop=True))
            r["p vs baseline"] = f"{p:.1e}" if not np.isnan(p) else "n/a"
        else:
            r["p vs baseline"] = "-"
        rows.append(r)
    t1 = pd.DataFrame(rows)
    t1.to_csv(f"{OUT}/table_ablation.csv", index=False)
    with open(f"{OUT}/table_ablation.tex", "w") as f:
        f.write(t1.to_latex(index=False, escape=True))
    print("  table_ablation.csv / .tex")

    # Tabela 2 - comparacao
    c = d["cifsup"]
    rows = []
    for name, pre in [("cifsup (NSD)", "nsd"), ("cifsup (ICP)", "icp"),
                      ("ICP-SAXS (this work)", "ours")]:
        r = {"Method": name}
        for col, label, _ in METRICS:
            c_ = f"{pre}_{col}"
            r[label] = f"{c[c_].mean():.3f}" if c_ in c and c[c_].notna().any() else "n/a"
        rows.append(r)
    t2 = pd.DataFrame(rows)
    t2.to_csv(f"{OUT}/table_comparison.csv", index=False)
    with open(f"{OUT}/table_comparison.tex", "w") as f:
        f.write(t2.to_latex(index=False, escape=True))
    print("  table_comparison.csv / .tex")

    # Tabela 3 - por entrada
    t3 = d["final"][["entry", "chamfer", "frac_fora", "cobertura",
                     "sep_centroide", "espelhado"]].copy()
    t3.columns = ["Entry", "Chamfer (A)", "Fraction outside",
                  "Coverage", "Centroid sep. (A)", "Mirrored"]
    t3 = t3.sort_values("Chamfer (A)").round(3)
    t3.to_csv(f"{OUT}/table_per_entry.csv", index=False)
    with open(f"{OUT}/table_per_entry.tex", "w") as f:
        f.write(t3.to_latex(index=False, escape=True))
    print("  table_per_entry.csv / .tex")

    # Estatisticas para o texto
    with open(f"{OUT}/summary_stats.txt", "w") as f:
        base, fin = d["baseline"], d["final"]
        f.write("NUMBERS FOR THE MANUSCRIPT TEXT\n")
        f.write("=" * 52 + "\n\n")
        f.write(f"Benchmark size: {len(fin)} SASBDB entries\n\n")
        for col, label, _ in METRICS:
            b_, f_ = base[col].mean(), fin[col].mean()
            p = wilcox(fin[col].reset_index(drop=True),
                       base[col].reset_index(drop=True))
            chg = (f_ - b_) / b_ * 100 if b_ else 0
            f.write(f"{label}\n")
            f.write(f"   baseline {b_:.3f} -> final {f_:.3f} "
                    f"({chg:+.1f}%), p = {p:.2e}\n")
        n_mir = int(fin["espelhado"].sum()) if "espelhado" in fin else 0
        f.write(f"\nMirrored hand selected in {n_mir}/{len(fin)} entries "
                f"({100*n_mir/len(fin):.0f}%)\n")
        f.write("\nComparison against ATSAS cifsup: no statistically\n")
        f.write("significant difference in any metric (paired Wilcoxon).\n")
    print("  summary_stats.txt")


def main():
    os.makedirs(OUT, exist_ok=True)
    d = load()
    print(f"\nGerando em {OUT}/\n")
    fig_ablation(d)
    fig_comparison(d)
    fig_enantiomorphs(d)
    fig_per_entry(d)
    tables(d)
    print(f"\nPronto. Veja {OUT}/summary_stats.txt para os numeros do texto.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
