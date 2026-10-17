#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | The figure of the series ranking (scripts/46): for each target, how long a dosing interval
     each of the mimic series allows, over the whole scan of the repression strength rho, with
     the renders of the molecules that the target belongs to above the plots.

     Row 1 places the PyMOL renders of scripts/39 (BACE1 and the Abeta42 fibril for the
     Alzheimer's column, the alpha-synuclein fibril for the two Parkinson's columns); they are
     illustrations and no value is read from them. Row 2 draws, for every series and every
     value of rho on the grid of scripts/40, the longest dosing interval at which the repressor
     peak needed to hold the cycle-averaged protein down by the target stays within 10 times the
     baseline pool (the interpolation of scripts/40 itself, `interval_for_peak`); the series of
     first place in the primary reading of scripts/46 is drawn thicker and shaded between its 3x
     and 30x limits, and every other series is drawn beside it with the same rule. Row 3 gives
     the three anchors of rho that scripts/46 ranks at, and the 3x to 30x range at the midpoint.
     Nothing is computed here that scripts/40 and scripts/46 do not already hold.

PT | A figura da ordenacao das series (scripts/46): para cada alvo, que intervalo de dose cada
     serie de mimetico permite, ao longo de toda a varredura da forca de repressao rho, com os
     renders das moleculas a que o alvo pertence acima dos graficos.

     A linha 1 coloca os renders do PyMOL do scripts/39 (BACE1 e a fibrila de Abeta42 na coluna
     da Alzheimer, a fibrila de alfa-sinucleina nas duas colunas da Parkinson); sao ilustracao
     e nenhum valor e lido deles. A linha 2 desenha, para toda serie e todo valor de rho da grade
     do scripts/40, o maior intervalo de dose em que o pico do repressor necessario para manter a
     proteina media do ciclo abaixo pelo alvo fica ate 10 vezes o reservatorio basal (a propria
     interpolacao do scripts/40, `interval_for_peak`); a serie do primeiro lugar na leitura
     primaria do scripts/46 e desenhada mais grossa e sombreada entre seus limites de 3x e 30x, e
     toda outra serie e desenhada ao lado pela mesma regra. A linha 3 da as tres ancoras de rho em
     que o scripts/46 ordena, e a faixa de 3x a 30x no ponto medio. Nada e calculado aqui que o
     scripts/40 e o scripts/46 nao tenham.

    python scripts/47_mimic_series_figure.py [--preview]     # needs scripts/39, 40 and 46
"""

import importlib.util
import json
import math
import os
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m38 = load("m38", "38_mechanistic_figures.py")
m40 = load("m40", "40_mimic_design_space.py")
import _viz_style as vs
import _viz_mech as vm
from _bilingual import LANGS, t as tr, fig_path

FIG, TAB, MM, DPI = m38.FIG, m38.TAB, m38.MM, m38.DPI
STEM = "mech_series_ranking"
FW, FH = 174.0, 168.0
THR = (3, 10, 30)
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]
STYLE = {   # species -> (line style, label EN, label PT)
    "miR-7": ("-", "miR-7 (1.7 h)", "miR-7 (1,7 h)"),
    "miR-7 upper bound": ((0, (5, 2)), "miR-7, upper bound (5 h)", "miR-7, limite superior (5 h)"),
    "miR-29b": ("-", "miR-29b (7 h)", "miR-29b (7 h)"),
    "miR-29c": ((0, (5, 2)), "miR-29c (10.6 h)", "miR-29c (10,6 h)"),
    "median miRNA": ((0, (1, 1.5)), "median miRNA (34 h), reference", "miRNA mediano (34 h), referência"),
}


def series_matrix(d, axis, species, kd):
    ds = d["ds"]
    Ts, rhos = list(ds["interval_grid_h"]), list(ds["rho_grid"])
    P = [[None] * len(Ts) for _ in rhos]
    for r in d["ds_rows"]:
        if (r["axis"], r["species"], round(float(r["target_knockdown"]), 4)) != (axis, species, round(kd, 4)):
            continue
        ti = Ts.index(float(r["interval_h"]))
        ri = int(np.argmin(np.abs(np.log10(rhos) - math.log10(float(r["rho"])))))
        if r["reachable"] == "True":
            P[ri][ti] = float(r["peak_pool_over_baseline"])
    return Ts, rhos, P


def frontier(Ts, rhos, P, level):
    out = []
    for ri in range(len(rhos)):
        h, flag = m40.interval_for_peak(Ts, P[ri], level)
        out.append((h, flag))
    return out


def colour_of(axis, species):
    return vs.INK2 if species == "median miRNA" else (vs.AD if axis == "AD" else vs.PD)


def text_colour(axis, species):
    return vs.INK2 if species == "median miRNA" else vs.DISEASE_TEXT[axis]


def draw(d, rank, lang):
    ds = d["ds"]
    cases = ds["interval_within_threshold"]
    kap_ad = d["ode"]["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
    kap_a, kap_b = m38.kv(d, "K058"), m38.kv(d, "K063")
    cols = [("AD", kap_ad, tr(lang, f"Alzheimer's disease · fall required {100 * kap_ad:.1f}%", f"Alzheimer · queda exigida {100 * kap_ad:.1f}%"), ["miR-29b", "miR-29c", "median miRNA"]),
            ("PD", kap_a, tr(lang, f"Parkinson's disease · fall {100 * kap_a:.0f}% (K058)", f"Parkinson · queda {100 * kap_a:.0f}% (K058)"), ["miR-7", "miR-7 upper bound", "median miRNA"]),
            ("PD", kap_b, tr(lang, f"Parkinson's disease · fall {100 * kap_b:.0f}% (K063)", f"Parkinson · queda {100 * kap_b:.0f}% (K063)"), ["miR-7", "miR-7 upper bound", "median miRNA"])]
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    gs_f = fig.add_gridspec(1, 3, left=0.082, right=0.985, top=0.745, bottom=0.470, wspace=0.30)
    gs_d = fig.add_gridspec(1, 3, left=0.082, right=0.985, top=0.352, bottom=0.150, wspace=0.30)
    xt = [0.01, 1, 100, 10000]
    for k, (axis, kd, title, order) in enumerate(cols):
        feat = rank["featured_by_axis"][axis]["species"][0]
        anchors = {a: cases[f"{axis}|{feat}|{kd:.4f}|{a}"]["rho"] for a in ("scan low", "scan midpoint", "scan high")}
        floor_key = f"{axis}|{feat}|{kd:.4f}|K065 floor"
        # ---------------------------------------------------------------- row 2: frontier over rho
        ax = fig.add_subplot(gs_f[0, k])
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(0.01, 1e4); ax.set_ylim(1.4, 260)
        ax.set_xticks(xt); ax.set_xticklabels(["0.01", "1", "100", "10⁴"])
        ax.set_yticks([2, 4, 8, 24, 48, 168]); ax.set_yticklabels(["2", "4", "8", "24", "48", "168"])
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        m38.style_axes(ax)
        for h_, lab in ((24, tr(lang, "daily", "diária")), (168, tr(lang, "weekly", "semanal"))):
            ax.axhline(h_, color=vs.RULE, lw=0.6, ls=(0, (3, 2)), zorder=1)
            if k == 0:
                ax.text(0.0115, h_ * 1.07, lab, fontsize=6.5, color=vs.NOTE, va="bottom", ha="left")
        for a, rho in anchors.items():
            ax.axvline(rho, color=vs.RULE, lw=0.7, ls=(0, (1.5, 2.0)), zorder=1)
        if floor_key in cases:
            ax.axvline(cases[floor_key]["rho"], color=vs.INK, lw=0.9, zorder=1)
            ax.text(cases[floor_key]["rho"] * 1.12, 1.55, tr(lang, "K065 floor", "piso K065"), fontsize=6.5, color=vs.INK2, va="bottom", ha="left", rotation=90)
        ends = []
        for sp in order:
            Ts, rhos, P = series_matrix(d, axis, sp, kd)
            f10 = frontier(Ts, rhos, P, 10)
            ys = np.array([h if h > 0 else np.nan for h, _ in f10])
            col, ls = colour_of(axis, sp), STYLE[sp][0]
            if sp == feat:
                lo = np.array([h if h > 0 else np.nan for h, _ in frontier(Ts, rhos, P, 3)])
                hi = np.array([h if h > 0 else np.nan for h, _ in frontier(Ts, rhos, P, 30)])
                ax.fill_between(rhos, np.where(np.isnan(lo), 1.4, lo), hi, color=col, alpha=0.16, lw=0, zorder=2)
                ax.plot(rhos, ys, color=col, lw=2.4, ls="-", zorder=5, solid_capstyle="round")
            else:
                ax.plot(rhos, ys, color=col, lw=1.15, ls=ls, zorder=3)
            none_x = [r for r, (h, fl) in zip(rhos, f10) if fl == "none"]
            if none_x:
                ax.scatter(none_x, [1.55] * len(none_x), s=6, facecolor="white", edgecolor=col, lw=0.7, zorder=6, clip_on=False)
            last = ys[~np.isnan(ys)][-1] if np.any(~np.isnan(ys)) else 1.6
            ends.append((last, sp))
        ends.sort(reverse=True)      # EN/PT: the highest line is labelled above, the others below their lines | a mais alta acima, as outras abaixo
        placed = []
        for n_, (y, sp) in enumerate(ends):
            above = n_ == 0
            yy = y * (1.10 if above else 0.90)
            while any(abs(math.log(yy) - math.log(p_)) < 0.17 for p_ in placed):
                yy *= 1.10 if above else 0.90
            placed.append(yy)
            lab = (STYLE[sp][1] if lang == "en" else STYLE[sp][2]).replace(", reference", "").replace(", referência", "")
            ax.text(9000, yy, lab, fontsize=6.3, color=text_colour(axis, sp), va="bottom" if above else "top", ha="right",
                    fontweight="bold" if sp == feat else "normal", path_effects=HALO, zorder=8)
        ax.set_title(title, fontsize=7.4, fontweight="bold", color=vs.DISEASE_TEXT[axis], loc="left", pad=6)
        ax.set_xlabel(tr(lang, "ρ, miRNA-driven over intrinsic mRNA decay", "ρ, decaimento por miRNA sobre o intrínseco"), fontsize=6.8, color=vs.INK2, labelpad=2)
        if k == 0:
            ax.set_ylabel(tr(lang, "Longest dosing interval within 10× (h)", "Maior intervalo de dose dentro de 10× (h)"), fontsize=6.8, color=vs.INK2, labelpad=2)
        m38.label_panel(ax, "abc"[k], x=-0.20 if k == 0 else -0.15, y=1.10)
        # ---------------------------------------------------------------- row 3: anchors and thresholds
        bx = fig.add_subplot(gs_d[0, k])
        bx.set_xscale("log"); bx.set_xlim(1.0, 260); bx.set_ylim(-0.7, len(order) - 0.3)
        bx.set_xticks([1, 2, 4, 8, 24, 48, 168]); bx.set_xticklabels(["0", "2", "4", "8", "24", "48", "168"])
        bx.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        bx.set_yticks([])
        m38.style_axes(bx)
        bx.grid(axis="y", visible=False)
        for j, sp in enumerate(order):
            y = len(order) - 1 - j
            col = colour_of(axis, sp)
            def val(a, t_):
                c = cases[f"{axis}|{sp}|{kd:.4f}|{a}"]
                return max(c[f"interval_within_{t_}x_h"] or 0.0, 1.0)
            lo, hi = val("scan midpoint", 3), val("scan midpoint", 30)
            bx.plot([lo, hi], [y, y], color=col, lw=4.6 if sp == feat else 3.2, alpha=0.30, solid_capstyle="butt", zorder=2)
            for a, mk in (("scan low", "o"), ("scan midpoint", "D"), ("scan high", "s")):
                v = val(a, 10)
                bx.scatter([v], [y], s=24 if a == "scan midpoint" else 14, marker=mk,
                           facecolor=col if v > 1.0 else "white", edgecolor=col, lw=0.9, zorder=4 + (a == "scan midpoint"))
            mid = cases[f"{axis}|{sp}|{kd:.4f}|scan midpoint"]["interval_within_10x_h"] or 0.0
            lab = STYLE[sp][1] if lang == "en" else STYLE[sp][2]
            val_txt = f"{mid:.0f} h" if mid < 167.5 else "≥168 h"
            bx.text(1.06, y + 0.30, f"{lab}: {val_txt}", fontsize=6.4, color=text_colour(axis, sp), va="bottom", ha="left",
                    fontweight="bold" if sp == feat else "normal")
        bx.set_xlabel(tr(lang, "Longest dosing interval within 10× (h)", "Maior intervalo de dose dentro de 10× (h)"), fontsize=6.8, color=vs.INK2, labelpad=2)
        m38.label_panel(bx, "def"[k], x=-0.20 if k == 0 else -0.15, y=1.06)
    # ---------------------------------------------------------------- row 1: renders
    items = [("a2_bace1_inhibitor", 0.150, 27.0, 25.0, "BACE1 · 4D8C"),
             ("a2_abeta42_fibril_axis", 0.290, 27.0, 25.0, "Aβ42 fibril · 5OQV"),
             ("a2_asyn_fibril_axis", 0.725, 44.0, 25.0, "α-synuclein fibril · 6CU7")]
    for stem, cx, bw, bh, cap in items:
        vm.place_image_fit(fig, f"{FIG}/structures/{stem}.png", cx, 0.905, bw, bh, FW, FH)
        fig.text(cx, 0.832, cap, fontsize=6.5, color=vs.NOTE, ha="center", va="center", zorder=10)
    # ---------------------------------------------------------------- keys
    ky = 0.070
    fig.add_artist(Line2D([0.082, 0.112], [ky, ky], color=vs.INK2, lw=2.4, transform=fig.transFigure))
    fig.text(0.119, ky, tr(lang, "first place, primary reading; band: 3× to 30×", "primeiro lugar, leitura primária; faixa: 3× a 30×"), fontsize=6.5, color=vs.INK2, va="center")
    fig.add_artist(Line2D([0.082, 0.112], [ky - 0.022, ky - 0.022], color=vs.INK2, lw=1.15, ls="-", transform=fig.transFigure))
    fig.text(0.119, ky - 0.022, tr(lang, "thin line: other measured series, same rule; dotted: reference of median stability (34 h)", "linha fina: outras séries medidas, mesma regra; pontilhada: referência de estabilidade mediana (34 h)"), fontsize=6.5, color=vs.INK2, va="center")
    # EN/PT: no text may leave the page, or the tight box grows past 174 mm | nenhum texto pode sair da pagina
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W_px, H_px = fig.get_size_inches() * fig.dpi
    out_of_page = [tx.get_text()[:40] for tx in fig.findobj(matplotlib.text.Text)
                   if tx.get_text().strip() and tx.get_visible()
                   and (tx.get_window_extent(rend).x1 > W_px + 1 or tx.get_window_extent(rend).x0 < -1)]
    if out_of_page:
        sys.exit(f"EN/PT: text outside the page in {lang}: {out_of_page}")
    fig.savefig(fig_path(FIG, STEM, lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    d = m38.read_inputs()
    rp = f"{TAB}/mimic_series_ranking.json"
    if not os.path.exists(rp):
        sys.exit("EN/PT: missing results/tables/mimic_series_ranking.json; run scripts/46 first")
    rank = json.load(open(rp, encoding="utf-8"))
    os.makedirs(FIG, exist_ok=True)
    for lang in LANGS:
        draw(d, rank, lang)
    print("=" * 78)
    print("EN | Figure of the series ranking | PT | Figura da ordenação das séries")
    print("=" * 78)
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG, STEM, lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
