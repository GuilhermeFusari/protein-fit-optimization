#!/usr/bin/env python3
"""
ensemble_chi2.py  -  valida um ensemble empacotado contra a curva de
espalhamento experimental por chi^2.

Fisica: moleculas em solucao espalham de forma independente (solucao
diluida), entao o perfil de um ensemble conformacional e a MEDIA dos perfis
dos membros. As posicoes das copias dentro do envelope NAO entram: o
chi^2 avalia quais conformacoes foram escolhidas, nao onde foram colocadas.
Calcular o CRYSOL no arquivo combinado (ALL_TOP_ALIGNED.pdb) estaria
errado: trataria as copias como um unico complexo.

Pipeline:
  1. para cada membro (arquivos .pdb da pasta, exceto ALL_TOP_ALIGNED.pdb):
     copia temporaria com formato corrigido para o CRYSOL (residuos rA/rU/rG/rC
     -> A/U/G/C; ocupancia vazia -> 1.00; coordenadas e atomos inalterados)
     e CRYSOL 3.2 com --explicit-hydrogens=yes (so os atomos do arquivo; os
     modelos nao tem hidrogenios, entao eles nao entram) e parametros de
     hidratacao default (nao ajustados)
  2. perfil de cada membro (coluna I total do .int) interpolado no q
     experimental
  3. perfil do ensemble = media simples (pesos iguais)
  4. chi^2 reduzido contra I(q) experimental com um fator de escala c
     ajustado analiticamente: chi2 = sum(((I - c*Icalc)/sigma)^2) / (N - 1)
     (e, como variante, com constante aditiva: N - 2)
  5. o mesmo para cada membro isolado (controle)

Uso:
    python3 ensemble_chi2.py <pasta do ensemble> <experimental.dat> <saida/>
        [--nome ROTULO]
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

RNA = {"rA": "A", "rU": "U", "rG": "G", "rC": "C"}


def read_dat(path):
    q, i, s = [], [], []
    for line in open(path, errors="replace"):
        p = line.split()
        if len(p) < 3:
            continue
        try:
            a, b, c = float(p[0]), float(p[1]), float(p[2])
        except ValueError:
            continue
        if c > 0:
            q.append(a); i.append(b); s.append(c)
    return np.array(q), np.array(i), np.array(s)


def sanitize(src, dst):
    out = []
    for line in open(src, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            l = line.rstrip("\n").ljust(80)
            res = l[17:20].strip()
            res = RNA.get(res, res)
            occ = l[54:60] if l[54:60].strip() else "  1.00"
            out.append((l[:17] + f"{res:>3s}" + l[20:54] + occ + l[60:]).rstrip() + "\n")
        elif line.startswith(("END", "TER")):
            out.append(line)
    open(dst, "w").writelines(out)


def crysol_profile(pdb, workdir, smax, ns):
    name = os.path.splitext(os.path.basename(pdb))[0]
    fixed = os.path.join(workdir, f"{name}.pdb")
    sanitize(pdb, fixed)
    p = subprocess.run(["crysol", fixed, "--explicit-hydrogens=yes",
                        f"--smax={smax}", f"--ns={ns}"],
                       cwd=workdir, capture_output=True, text=True)
    intf = os.path.join(workdir, f"{name}.int")
    log = open(os.path.join(workdir, f"{name}.log"), errors="replace").read() \
        if os.path.exists(os.path.join(workdir, f"{name}.log")) else ""
    if p.returncode != 0 or not os.path.exists(intf):
        raise RuntimeError(f"crysol falhou em {pdb}: {p.stdout[-300:]} {p.stderr[-300:]}")
    unknown = [l for l in log.splitlines() if "unknown atoms" in l.lower()]
    natoms = [l for l in log.splitlines() if "Total number of atoms read" in l]
    d = np.loadtxt(intf, skiprows=1)
    return d[:, 0], d[:, 1], (natoms[0].split(":")[-1].strip() if natoms else "?"), \
        (unknown[0].split(":")[-1].strip() if unknown else "?")


def chi2(I, s, calc, constant=False):
    w = 1.0 / s ** 2
    if constant:
        A = np.stack([calc, np.ones_like(calc)], 1)
        coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], I * np.sqrt(w), rcond=None)
        fit = A @ coef
        dof = len(I) - 2
    else:
        c = (w * I * calc).sum() / (w * calc ** 2).sum()
        fit, coef, dof = c * calc, np.array([c]), len(I) - 1
    return float((((I - fit) / s) ** 2).sum() / dof), fit, coef


def plot(q, I, s, ens_fit, best_fit, best_name, chi_e, chi_b, title, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    BLUE, ORANGE = "#2a78d6", "#eb6834"
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.2), sharex=True,
                                  gridspec_kw=dict(height_ratios=[3, 1.2], hspace=0.08))
    fig.patch.set_facecolor(SURF)
    for a in (ax, axr):
        a.set_facecolor(SURF)
        a.grid(True, color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a.spines[sp].set_color(INK2)
        a.tick_params(colors=INK2, labelsize=9)
    ax.errorbar(q, I, yerr=s, fmt="o", ms=2.5, color=INK2, ecolor=GRID, elinewidth=0.8,
                label="Experimental I(q)", zorder=2)
    ax.plot(q, best_fit, color=ORANGE, lw=2, ls=(0, (4, 2)), zorder=3,
            label=f"Best single member ({best_name}), χ² = {chi_b:.2f}")
    ax.plot(q, ens_fit, color=BLUE, lw=2, zorder=4,
            label=f"Ensemble average, χ² = {chi_e:.2f}")
    ax.set_yscale("log")
    ax.set_ylabel("I(q) (arb. units)", color=INK, fontsize=10)
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK)
    axr.axhline(0, color=INK2, lw=0.8)
    axr.plot(q, (I - ens_fit) / s, color=BLUE, lw=1.2)
    axr.set_ylabel("(I − fit)/σ\nensemble", color=INK, fontsize=9)
    axr.set_xlabel("q (Å⁻¹)", color=INK, fontsize=10)
    fig.savefig(out_png, dpi=160, bbox_inches="tight", facecolor=SURF)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ensemble")
    ap.add_argument("dat")
    ap.add_argument("out")
    ap.add_argument("--nome", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    q, I, s = read_dat(a.dat)
    smax = float(np.ceil((q.max() + 0.005) * 100) / 100)
    ns = int(round(smax / 0.001)) + 1
    members = sorted(f for f in glob.glob(os.path.join(a.ensemble, "*.pdb"))
                     if os.path.basename(f) != "ALL_TOP_ALIGNED.pdb")
    if not members:
        raise SystemExit(f"nenhum .pdb em {a.ensemble}")

    work = tempfile.mkdtemp(prefix="ens_chi2_")
    try:
        prof, info = [], []
        for m in members:
            qc, Ic, nat, unk = crysol_profile(m, work, smax, ns)
            prof.append(np.interp(q, qc, Ic))
            info.append(dict(membro=os.path.basename(m), atomos_lidos=nat,
                             atomos_desconhecidos=unk))
            print(f"  {os.path.basename(m)}: {nat} atomos, {unk} desconhecidos", flush=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    prof = np.array(prof)
    ens = prof.mean(0)

    chi_e, fit_e, _ = chi2(I, s, ens)
    chi_ec, _, _ = chi2(I, s, ens, constant=True)
    rows = []
    for k, p in enumerate(prof):
        c1, f1, _ = chi2(I, s, p)
        c2, _, _ = chi2(I, s, p, constant=True)
        rows.append(dict(info[k], chi2=round(c1, 3), chi2_com_constante=round(c2, 3)))
        info[k]["_fit"] = f1
    best = int(np.argmin([r["chi2"] for r in rows]))

    import csv
    with open(os.path.join(a.out, "chi2_membros.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    np.savetxt(os.path.join(a.out, "perfis_membros.dat"),
               np.column_stack([q, I, s, prof.T]), fmt="%.6g",
               header="q I_exp sigma " + " ".join(r["membro"] for r in rows))
    np.savetxt(os.path.join(a.out, "perfil_ensemble.dat"),
               np.stack([q, I, s, fit_e], 1), fmt="%.6g",
               header="q I_exp sigma I_ensemble_escalado")
    res = dict(nome=a.nome, ensemble=os.path.abspath(a.ensemble),
               experimental=os.path.abspath(a.dat), n_membros=len(members),
               n_pontos=len(q), q_min=float(q.min()), q_max=float(q.max()),
               chi2_ensemble=round(chi_e, 3), chi_ensemble=round(float(np.sqrt(chi_e)), 3),
               chi2_ensemble_com_constante=round(chi_ec, 3),
               chi2_melhor_membro=rows[best]["chi2"], melhor_membro=rows[best]["membro"],
               chi2_membros_mediana=round(float(np.median([r["chi2"] for r in rows])), 3),
               chi2_membros_pior=round(float(max(r["chi2"] for r in rows)), 3))
    json.dump(res, open(os.path.join(a.out, "resultado.json"), "w"), indent=2)
    plot(q, I, s, fit_e, info[best]["_fit"], rows[best]["membro"].replace(".pdb", ""),
         chi_e, rows[best]["chi2"],
         f"{a.nome}: ensemble of {len(members)} members vs experimental SAXS",
         os.path.join(a.out, "perfil_ensemble_vs_experimental.png"))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
