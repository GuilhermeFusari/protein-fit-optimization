#!/usr/bin/env python3
"""
gerar_todas_figuras.py

Gera todas as 7 figuras do artigo, em Inglês e Português.
Corrige sobreposições nas barras de erro (ablation) e no p-valor (cifsup).

Uso:
    python3 gerar_todas_figuras.py
"""

import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DPI = 300
C_MAIN = "#2C5F8A"
C_ALT = "#D98E3A"
C_GREY = "#B0B0B0"
C_HL = "#1E7A46"
C_STEP = "#7BA7C7"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
})

STRINGS = {
    "en": {
        "orig_hand": "Original hand", "mirr_hand": "Mirrored hand",
        "n_entries": "Number of entries", "exp_50": " expected 50%",
        "title_enantio": "Enantiomorph selection ({n_mir}/{n} mirrored)\nclose to the 50% expected when the envelope\ncarries no chirality information",
        "al_packing": "Aligned packing", "rand_ctrl": "Random control",
        "chamfer_dist": "Chamfer distance (Å)",
        "title_packing": "Multi-copy packing versus a random-placement control\nthe gain comes from the alignment, not the point count",
        "id_fit": "identical fit",
        "chamfer_true": "Chamfer of the true pose (Å)",
        "chamfer_rec": "Chamfer of the recovered pose (Å)",
        "title_gt_a": "Recovered poses fit the envelope\nas well as the true pose",
        "globular": "globular\n(<1.5)", "medium": "medium\n(1.5–2.5)", "elongated": "elongated\n(>2.5)",
        "ori_error": "Orientation error (degrees)",
        "title_gt_b": "Orientation ambiguity decreases\nfor anisometric particles",
        "title_gt_sup": "Pose-recovery test with known ground truth (n = {n} trials)",
        "unc_icp": "Unconstrained ICP", "bidir_pen": "Bidirectional penalty",
        "no_down": "No downsampling", "enant_search": "Enantiomorph search",
        "frac_out": "Fraction outside envelope", "env_cov": "Envelope coverage",
        "cent_sep": "Centroid separation (Å)",
        "better_down": "↓ better", "better_up": "↑ better",
        "title_ablation_sup": "Cumulative effect of each methodological component (n = {n})",
        "cif_nsd": "cifsup\n(NSD)", "cif_icp": "cifsup\n(ICP)", "ours": "ICP-SAXS\n(ours)",
        "title_comp_sup": "ICP-SAXS versus ATSAS cifsup on the same 50 SASBDB entries",
        "orig_kept": "Original hand kept", "mirr_sel": "Mirrored hand selected",
        "ch_no_mirr": "Chamfer distance, no mirror search (Å)", "ch_mirr": "Chamfer distance, with mirror search (Å)",
        "per_entry": "Per-entry effect", "sel_hand": "Selected hand",
        "title_mirr_sup": "SAXS envelopes do not define chirality: both hands must be tested",
        "median": "Median", "mean": "Mean",
        "entry_x": "SASBDB entry (sorted by alignment quality)",
        "title_entry": "Alignment quality across all {n} benchmark entries"
    },
    "pt": {
        "orig_hand": "Mão original", "mirr_hand": "Mão espelhada",
        "n_entries": "Número de entradas", "exp_50": " esperado 50%",
        "title_enantio": "Seleção de enantiomorfos ({n_mir}/{n} espelhados)\npróximo aos 50% esperados quando o envelope\nnão carrega informação de quiralidade",
        "al_packing": "Empacotamento alinhado", "rand_ctrl": "Controle aleatório",
        "chamfer_dist": "Distância Chamfer (Å)",
        "title_packing": "Empacotamento de múltiplas cópias vs controle aleatório\no ganho vem do alinhamento, não da contagem de pontos",
        "id_fit": "encaixe idêntico",
        "chamfer_true": "Chamfer da pose verdadeira (Å)",
        "chamfer_rec": "Chamfer da pose recuperada (Å)",
        "title_gt_a": "Poses recuperadas se encaixam no envelope\ntão bem quanto a pose verdadeira",
        "globular": "globular\n(<1.5)", "medium": "médio\n(1.5–2.5)", "elongated": "alongado\n(>2.5)",
        "ori_error": "Erro de orientação (graus)",
        "title_gt_b": "A ambiguidade de orientação diminui\npara partículas anisométricas",
        "title_gt_sup": "Teste de recuperação de pose com ground truth conhecido (n = {n})",
        "unc_icp": "ICP sem restrição", "bidir_pen": "Penalidade bidirecional",
        "no_down": "Sem downsampling", "enant_search": "Busca de enantiomorfos",
        "frac_out": "Fração fora do envelope", "env_cov": "Cobertura do envelope",
        "cent_sep": "Separação de centroides (Å)",
        "better_down": "↓ melhor", "better_up": "↑ melhor",
        "title_ablation_sup": "Efeito cumulativo de cada componente metodológico (n = {n})",
        "cif_nsd": "cifsup\n(NSD)", "cif_icp": "cifsup\n(ICP)", "ours": "ICP-SAXS\n(nosso)",
        "title_comp_sup": "ICP-SAXS versus ATSAS cifsup nas mesmas 50 entradas do SASBDB",
        "orig_kept": "Mão original mantida", "mirr_sel": "Mão espelhada selecionada",
        "ch_no_mirr": "Distância Chamfer, sem busca de espelho (Å)", "ch_mirr": "Distância Chamfer, com busca de espelho (Å)",
        "per_entry": "Efeito por entrada", "sel_hand": "Mão selecionada",
        "title_mirr_sup": "Envelopes SAXS não definem quiralidade: ambas as mãos devem ser testadas",
        "median": "Mediana", "mean": "Média",
        "entry_x": "Entrada SASBDB (ordenada por qualidade de alinhamento)",
        "title_entry": "Qualidade do alinhamento em todas as {n} entradas do benchmark"
    }
}

def t(key, lang, **kwargs):
    text = STRINGS[lang].get(key, key)
    if kwargs:
        return text.format(**kwargs)
    return text

FILES = {
    "baseline": "baseline_icp_puro.csv",
    "lam2_300": "confirma_v2.csv",
    "nodown": "confirma_sem_downsampling.csv",
    "final": "final_v2.csv",
    "cifsup": "benchmark_cifsup_v2.csv",
    "final_lam02": "final_lam02.csv",
    "gt": "ground_truth_final.csv"
}

def get_metrics(lang):
    return [
        ("chamfer", t("chamfer_dist", lang), "lower"),
        ("frac_fora", t("frac_out", lang), "lower"),
        ("cobertura", t("env_cov", lang), "higher"),
        ("sep_centroide", t("cent_sep", lang), "lower"),
    ]

def load():
    d = {}
    for k, f in FILES.items():
        if os.path.exists(f):
            df = pd.read_csv(f)
            if "status" in df.columns:
                df = df[df.status == "ok"]
            d[k] = df
    return d

def wilcox(a, b):
    m = pd.concat([a, b], axis=1).dropna()
    if len(m) < 5: return np.nan
    try: return stats.wilcoxon(m.iloc[:, 0], m.iloc[:, 1]).pvalue
    except ValueError: return np.nan

def stars(p):
    if np.isnan(p): return "n.s."
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."

def read_report_combined(path):
    if not os.path.exists(path): return None
    for line in open(path, errors="replace"):
        m = re.search(r"Global Error \(Combined\):\s*([\d.]+)", line)
        if m: return float(m.group(1))
    return None

def fig_enantiomeros(d, out, lang):
    if "final_lam02" not in d or "espelhado" not in d["final_lam02"]: return
    df = d["final_lam02"]
    n = len(df)
    n_mir = int(df.espelhado.sum())
    n_orig = n - n_mir

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.bar([0, 1], [n_orig, n_mir], color=[C_GREY, C_MAIN], edgecolor="black", lw=0.6, width=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([t("orig_hand", lang), t("mirr_hand", lang)])
    ax.set_ylabel(t("n_entries", lang))
    ax.axhline(n / 2, ls="--", color="#999999", lw=1)
    ax.text(0.98, n / 2, t("exp_50", lang), va="bottom", ha="right", transform=ax.get_yaxis_transform(), fontsize=9, color="#666666")
    for i, v in enumerate([n_orig, n_mir]):
        ax.text(i, v + 0.4, str(v), ha="center", fontsize=11, fontweight="bold")
    ax.set_title(t("title_enantio", lang, n_mir=n_mir, n=n), fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{out}/figura_enantiomeros.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_packing(out, lang):
    R = "SAXS_Protein_Aligner"
    casos = [
        ("GRB2\n(30 copies)", f"{R}/resultado_packing_30/report.txt", 5.690),
        ("GRB2\n(5 copies)", f"{R}/resultado_packing/report.txt", 7.034),
        ("RNA SASDFQ9\n(20 models)", f"{R}/resultado_packing_SASDFQ9/report.txt", 6.358),
    ]
    labels, real, ctrl = [], [], []
    for lab, rep, c in casos:
        v = read_report_combined(rep)
        if v is not None:
            labels.append(lab); real.append(v); ctrl.append(c)
    if not real: return

    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(len(labels))
    w = 0.38
    ax.bar(x - w / 2, real, w, label=t("al_packing", lang), color=C_MAIN, edgecolor="black", lw=0.6)
    ax.bar(x + w / 2, ctrl, w, label=t("rand_ctrl", lang), color=C_GREY, edgecolor="black", lw=0.6)
    for i, (rr, cc) in enumerate(zip(real, ctrl)):
        ax.text(i - w / 2, rr + 0.1, f"{rr:.2f}", ha="center", fontsize=8.5)
        ax.text(i + w / 2, cc + 0.1, f"{cc:.2f}", ha="center", fontsize=8.5)
        ax.text(i, max(rr, cc) + 0.5, f"{cc/rr:.1f}x", ha="center", fontsize=9, color=C_HL, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(t("chamfer_dist", lang))
    ax.set_ylim(0, max(ctrl) * 1.2)
    ax.legend(frameon=False)
    ax.set_title(t("title_packing", lang), fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{out}/figura_packing.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_ground_truth(d, out, lang):
    if "gt" not in d or "delta_chamfer" not in d["gt"]: return
    df = d["gt"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.scatter(df.chamfer_verdadeira, df.chamfer_recuperada, s=40, color=C_MAIN, edgecolor="black", lw=0.4, alpha=0.8)
    lim = [0, max(df.chamfer_verdadeira.max(), df.chamfer_recuperada.max()) * 1.1]
    ax.plot(lim, lim, "--", color="#999999", lw=1, label=t("id_fit", lang))
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(t("chamfer_true", lang))
    ax.set_ylabel(t("chamfer_rec", lang))
    ax.legend(frameon=False, loc="upper left")
    ax.set_title(t("title_gt_a", lang), fontsize=10)

    ax = axes[1]
    df = df.copy()
    c_labels = [t("globular", lang), t("medium", lang), t("elongated", lang)]
    df["classe"] = pd.cut(df.anisometria, [0, 1.5, 2.5, 100], labels=c_labels)
    grp = df.groupby("classe", observed=True)["erro_ang"]
    classes = list(grp.groups.keys())
    means = [grp.get_group(c).mean() for c in classes]
    sems = [grp.get_group(c).sem() for c in classes]
    xx = np.arange(len(classes))
    ax.bar(xx, means, yerr=sems, capsize=3, color=C_ALT, edgecolor="black", lw=0.6, width=0.6)
    ax.set_xticks(xx)
    ax.set_xticklabels([str(c) for c in classes])
    ax.set_ylabel(t("ori_error", lang))
    ax.set_title(t("title_gt_b", lang), fontsize=10)

    fig.suptitle(t("title_gt_sup", lang, n=len(df)), y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}/figura_ground_truth.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_ablation(d, out, lang):
    if not all(k in d for k in ["baseline", "lam2_300", "nodown", "final"]): return
    steps = [
        (t("unc_icp", lang), d["baseline"], C_GREY),
        (t("bidir_pen", lang), d["lam2_300"], C_STEP),
        (t("no_down", lang), d["nodown"], C_STEP),
        (t("enant_search", lang), d["final"], C_MAIN),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))
    metrics = get_metrics(lang)

    for ax, (col, label, better) in zip(axes, metrics):
        vals = [s[1][col].mean() for s in steps]
        errs = [s[1][col].std(ddof=1) / np.sqrt(len(s[1])) for s in steps]
        cols = [s[2] for s in steps]
        x = np.arange(len(steps))

        ax.bar(x, vals, yerr=errs, color=cols, edgecolor="black", lw=0.6, capsize=3, width=0.68)

        # valor acima da barra de erro
        for i, (v, e) in enumerate(zip(vals, errs)):
            offset = max(vals) * 0.03
            ax.text(i, v + e + offset, f"{v:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

        ax.set_ylim(0, max([v + e for v, e in zip(vals, errs)]) * 1.15)

        ax.set_xticks(x)
        ax.set_xticklabels([f"({i+1})" for i in range(len(steps))])
        ax.set_ylabel(label)
        arrow = t("better_down", lang) if better == "lower" else t("better_up", lang)
        ax.set_title(arrow, fontsize=9, color="#555555")

        p = wilcox(d["final"][col].reset_index(drop=True), d["baseline"][col].reset_index(drop=True))
        ax.text(0.5, 0.94, stars(p), transform=ax.transAxes, ha="center", va="top", fontsize=11)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=s[2], ec="black", lw=0.6) for s in steps]
    labels = [f"({i+1}) {s[0]}" for i, s in enumerate(steps)]
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle(t("title_ablation_sup", lang, n=len(d['final'])), y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}/figure_ablation.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_comparison(d, out, lang):
    if "cifsup" not in d: return
    c = d["cifsup"]
    groups = [
        (t("cif_nsd", lang), "nsd", C_ALT),
        (t("cif_icp", lang), "icp", C_ALT),
        (t("ours", lang), "ours", C_MAIN),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6))
    metrics = get_metrics(lang)

    for ax, (col, label, better) in zip(axes, metrics):
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
                        flierprops=dict(marker="o", ms=3, mfc="none", mec="#888888", alpha=0.6))
        for patch, colr in zip(bp["boxes"], cols):
            patch.set_facecolor(colr)
            patch.set_alpha(0.75)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.6)

        ax.set_xticklabels(names)
        ax.set_ylabel(label)
        arrow = t("better_down", lang) if better == "lower" else t("better_up", lang)
        ax.set_title(arrow, fontsize=9, color="#555555")

        a, b = f"ours_{col}", f"icp_{col}"
        if a in c and b in c:
            p = wilcox(c[a], c[b])
            if not np.isnan(p):
                # p-valor no topo do eixo
                ymin, ymax = ax.get_ylim()
                ax.set_ylim(ymin, ymax + (ymax - ymin) * 0.18)
                ax.text(0.5, 0.96, f"vs NSD: p = {p:.3f} ({stars(p)})",
                        transform=ax.transAxes, ha="center", va="top",
                        fontsize=8.5, color="#333333", fontweight="bold")

    fig.suptitle(t("title_comp_sup", lang), y=1.01, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}/figure_comparison.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_enantiomorphs(d, out, lang):
    if not all(k in d for k in ["nodown", "final"]): return
    a = d["nodown"][["entry", "chamfer"]].rename(columns={"chamfer": "off"})
    b = d["final"][["entry", "chamfer", "espelhado"]].rename(columns={"chamfer": "on"})
    m = a.merge(b, on="entry")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    mir = m[m.espelhado]
    nom = m[~m.espelhado]
    ax.scatter(nom.off, nom.on, s=32, facecolor="none", edgecolor="#888888", label=t("orig_kept", lang))
    ax.scatter(mir.off, mir.on, s=32, color=C_MAIN, edgecolor="black", lw=0.4, label=t("mirr_sel", lang))
    lim = [0, max(m.off.max(), m.on.max()) * 1.06]
    ax.plot(lim, lim, "--", color="#999999", lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(t("ch_no_mirr", lang))
    ax.set_ylabel(t("ch_mirr", lang))
    ax.legend(frameon=False, loc="upper left")
    ax.set_title(t("per_entry", lang))

    ax = axes[1]
    n_mir = int(m.espelhado.sum())
    ax.bar([0, 1], [len(m) - n_mir, n_mir], color=["#B0B0B0", C_MAIN], edgecolor="black", lw=0.6, width=0.55)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([t("orig_hand", lang), t("mirr_hand", lang)])
    ax.set_ylabel(t("n_entries", lang))
    ax.axhline(len(m) / 2, ls="--", color="#999999", lw=1)
    ax.text(0.98, len(m) / 2, t("exp_50", lang), va="bottom", ha="right", transform=ax.get_yaxis_transform(), fontsize=8, color="#666666")
    ax.set_title(f"{t('sel_hand', lang)} ({n_mir}/{len(m)})")

    fig.suptitle(t("title_mirr_sup", lang), y=1.00, fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}/figure_enantiomorphs.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def fig_per_entry(d, out, lang):
    if "final" not in d: return
    f = d["final"][["entry", "chamfer"]].sort_values("chamfer").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 3.8))
    x = np.arange(len(f))
    ax.bar(x, f.chamfer, color=C_MAIN, edgecolor="black", lw=0.4, width=0.72)
    ax.axhline(f.chamfer.median(), ls="--", color="#C0392B", lw=1.2, label=f"{t('median', lang)}: {f.chamfer.median():.2f} Å")
    ax.axhline(f.chamfer.mean(), ls=":", color="#2C3E50", lw=1.2, label=f"{t('mean', lang)}: {f.chamfer.mean():.2f} Å")
    ax.set_xticks(x)
    ax.set_xticklabels(f.entry, rotation=90, fontsize=6)
    ax.set_ylabel(t("chamfer_dist", lang))
    ax.set_xlabel(t("entry_x", lang))
    ax.legend(frameon=False)
    ax.set_title(t("title_entry", lang, n=len(f)))
    fig.tight_layout()
    fig.savefig(f"{out}/figure_per_entry.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)

def main():
    d = load()
    if not d:
        print("Nenhum arquivo CSV encontrado! Verifique o diretório.")
        return 1

    for lang, folder in [("en", "ingles"), ("pt", "portugues")]:
        out_dir = f"imagens final/{folder}"
        os.makedirs(out_dir, exist_ok=True)
        print(f"\nGerando figuras em {out_dir}/ ...")

        fig_enantiomeros(d, out_dir, lang)
        fig_packing(out_dir, lang)
        fig_ground_truth(d, out_dir, lang)
        fig_ablation(d, out_dir, lang)
        fig_comparison(d, out_dir, lang)
        fig_enantiomorphs(d, out_dir, lang)
        fig_per_entry(d, out_dir, lang)

    print("\nPronto! Todas as figuras (inglês e português) foram geradas.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
