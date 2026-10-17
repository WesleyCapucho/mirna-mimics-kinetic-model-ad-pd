#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | The seven main-text figures of the mechanistic study (Article 2), drawn in the shared
     style of scripts/_viz_style.py and fitted to the journal's artwork limits by
     scripts/_journal_fit.py. Every number on a figure is read from the result tables of
     scripts/12, 17, 36 and 40 and from the kinetic parameter table at draw time, or follows
     from a closed form written next to the line it draws; nothing is typed in as a result.
     Where a drawing carries data in its geometry (the tanks of Fig. 3, the ribbons of Fig. 2)
     the sizes are computed from the table the figure is about. The molecular images are the
     renders of scripts/39 and are illustration, not a result. The one simulation, the time
     course of Fig. 6a, re-runs the cascade of scripts/36 at the dose that the scenario table
     of that script reports.

         mech_model_schematic   Fig. 1  four deposited structures, and the two cascades with their constants
         mech_provenance        Fig. 2  from 15 primary sources, through the system measured, to the model
         mech_alzheimer         Fig. 3  production, clearance and the pool; required knockdown; rates at the pool
         mech_parkinson         Fig. 4  saturation, transmission, mimic factor against rho
         mech_turnover          Fig. 5  time scales; washout, peak to average, stabilisation needed
         mech_target_response   Fig. 6  what a repeated dose does at the protein
         mech_design_space      Fig. 7  the peak a regimen needs, over repression strength and dosing interval

PT | As sete figuras do texto principal do estudo mecanistico (Artigo 2), no estilo
     compartilhado do scripts/_viz_style.py e ajustadas aos limites de arte do periodico
     pelo scripts/_journal_fit.py. Todo numero numa figura e lido, na hora de desenhar, das
     tabelas de resultado dos scripts 12, 17, 36 e 40 e da tabela de parametros cineticos, ou
     decorre de uma forma fechada escrita ao lado da linha que desenha; nada e digitado como
     resultado. Onde um desenho carrega dados na geometria (os tanques da Fig. 3, as fitas da
     Fig. 2) os tamanhos sao calculados da tabela de que a figura trata. As imagens
     moleculares sao os renders do scripts/39 e sao ilustracao, nao resultado. A unica
     simulacao, a curva no tempo da Fig. 6a, roda de novo a cascata do scripts/36 na dose que
     a tabela de cenarios daquele script reporta.

    python scripts/38_mechanistic_figures.py [--only=mech_alzheimer,mech_turnover] [--preview]
"""

import csv
import importlib.util
import json
import math
import os
import sys
import textwrap
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _viz_style as vs
import _viz_mech as vm
import _provenance as pc
from _bilingual import LANGS, t as tr, fig_path

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = "results/figures"
STRUCT_FIG = f"{FIG}/structures"
TAB = "results/tables"
KIN = "data/extracted/kinetic_parameters.csv"
LN2 = math.log(2.0)
MM = 1.0 / 25.4
ONLY = [a.split("=", 1)[1] for a in sys.argv if a.startswith("--only=")]
ONLY = [x for x in (ONLY[0].split(",") if ONLY else []) if x]
DPI = 110 if "--preview" in sys.argv else 300

SP_COLOUR = {   # EN/PT: species -> (colour, line style, label EN, label PT)
    "miR-7": (vs.PD, "-", "miR-7 (1.7 h)", "miR-7 (1.7 h)"),
    "miR-7 upper bound": (vs.PD, "--", "miR-7, upper bound (5 h)", "miR-7, limite superior (5 h)"),
    "miR-29b": (vs.AD, "-", "miR-29b (7 h)", "miR-29b (7 h)"),
    "miR-29c": (vs.AD, "--", "miR-29c (10.6 h)", "miR-29c (10.6 h)"),
    "median miRNA": (vs.INK2, "-", "median miRNA (34 h)", "miRNA mediano (34 h)"),
}


def load_module(name, filename):
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        sys.exit(f"EN/PT: missing {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_inputs():
    for f in (f"{TAB}/ode_calibrated_results.json", f"{TAB}/mimic_dosing_feasibility.json",
              f"{TAB}/mimic_target_response.json", f"{TAB}/mimic_target_response.csv",
              f"{TAB}/mimic_target_rho_sweep.csv", f"{TAB}/mimic_design_space.csv",
              f"{TAB}/mimic_design_space.json", f"{TAB}/structure_figure_provenance.json",
              f"{TAB}/structure_renders_article2.json", KIN):
        if not os.path.exists(f):
            sys.exit(f"EN/PT: missing {f}; run scripts/12, 13, 17, 36, 39 and 40 first")
    d = OrderedDict()
    d["kin"] = {r["param_id"]: r for r in csv.DictReader(open(KIN, encoding="utf-8"))}
    d["kin_rows"] = list(csv.DictReader(open(KIN, encoding="utf-8")))
    d["ode"] = json.load(open(f"{TAB}/ode_calibrated_results.json", encoding="utf-8"))
    d["dos"] = json.load(open(f"{TAB}/mimic_dosing_feasibility.json", encoding="utf-8"))
    d["tgt"] = json.load(open(f"{TAB}/mimic_target_response.json", encoding="utf-8"))
    d["rows"] = list(csv.DictReader(open(f"{TAB}/mimic_target_response.csv", encoding="utf-8")))
    d["sweep"] = list(csv.DictReader(open(f"{TAB}/mimic_target_rho_sweep.csv", encoding="utf-8")))
    d["ds"] = json.load(open(f"{TAB}/mimic_design_space.json", encoding="utf-8"))
    d["ds_rows"] = list(csv.DictReader(open(f"{TAB}/mimic_design_space.csv", encoding="utf-8")))
    d["prov"] = json.load(open(f"{TAB}/structure_figure_provenance.json", encoding="utf-8"))
    d["renders"] = json.load(open(f"{TAB}/structure_renders_article2.json", encoding="utf-8"))["renders"]
    return d


def kv(d, pid):
    return float(d["kin"][pid]["value_si"])


def label_panel(ax, letter, x=-0.15, y=1.08):
    """EN/PT: bold panel letter above and left of the axes, clear of the title | letra do painel fora do titulo."""
    vs.panel_letter(ax, letter, x=x, y=y, size=10)


def style_axes(ax, grid=True):
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(True, color=vs.GRID, lw=0.5)
        ax.set_axisbelow(True)


def fold_for(kappa, rho):
    """EN/PT: mimic factor that a knockdown kappa requires at strength rho | fator de mimetico."""
    return 1.0 + kappa / (1.0 - kappa) * (1.0 + rho) / rho


def knockdown_at(f, rho):
    """EN/PT: knockdown produced by a mimic factor f at strength rho | queda produzida."""
    return 1.0 - (1.0 + rho) / (1.0 + rho * f)


def method_short(prov, lang):
    m = prov["method"].upper()
    return tr(lang, "X-ray", "raios X") if "X-RAY" in m else "cryo-EM"


def structure_line(p, lang, extra=""):
    """EN/PT: 'PDB id, method, resolution' read from the provenance table | linha lida da tabela de proveniencia."""
    res = f"{p['resolution_angstrom']:g}"
    res += "" if "." in res else ".0"
    s = f"{p['pdb_id']} · {method_short(p, lang)}, {res} Å"
    if p.get("resolution_method") == "FSC 0.5 CUT-OFF":
        s += " (FSC 0.5)"
    return s + extra


# ----------------------------------------------------------------------------- Fig. 3
def fig_alzheimer(d, lang):
    ode, tgt = d["ode"], d["tgt"]
    clr = ode["ad_clearance_vs_production"]
    dose = ode["mimic_dose_vs_measured_knockdown"]
    kappa = dose["required_BACE1_knockdown"]
    measured = dose["measured_BACE1_knockdown_K047"]
    rs = tgt["repression_strength"]["AD"]
    nuc = tgt["nucleation_scaling"]
    ratio = clr["abeta_ratio_AD_over_control"]
    FW, FH = 174.0, 112.0
    TOP = 0.505                      # share of the figure taken by the tanks

    fig = plt.figure(figsize=(FW * MM, FH * MM))

    # ------------------------------------------------------------------ (a) tanks
    ax = fig.add_axes([0.0, 1.0 - TOP, 1.0, TOP])
    top_h = 100.0 * (FH * TOP) / FW
    ax.set_xlim(0, 100)
    ax.set_ylim(0, top_h)
    ax.axis("off")
    kp_c, kp_a = kv(d, "K001"), kv(d, "K002")
    kc_c, kc_a = kv(d, "K003"), kv(d, "K004")
    pool0 = kp_c / kc_c
    tanks = [(tr(lang, "Controls", "Controles"), kp_c, kc_c, vs.AD),
             (tr(lang, "Alzheimer's disease", "Doença de Alzheimer"), kp_a, kc_a, vs.AD),
             (tr(lang, f"AD, production\nlowered by {100 * kappa:.1f}%", f"DA, produção\nreduzida em {100 * kappa:.1f}%"),
              kp_a * (1.0 - kappa), kc_a, vs.TIER_RAMP[0])]
    TW, TH, base_y = 11.0, 15.0, 3.6
    unit_h = 8.8                    # height of the controls' pool; every other pool is scaled by the measured ratio
    w_rate = 0.30                   # arrow width per % per hour
    centres = [8.0, 29.0, 50.0]
    for (name, kp, kc, col), cx in zip(tanks, centres):
        lvl = (kp / kc) / pool0
        x0 = cx - TW / 2
        body = Rectangle((x0, base_y), TW, TH, fc="none", ec="none")
        ax.add_patch(body)
        vm.gradient_rect(ax, x0, x0 + TW, base_y, base_y + lvl * unit_h, "#174a8f", "#2f6fc0", z=2, clip=body)   # EN/PT: white lettering stays at 4.5:1 or more over the whole fill | texto branco fica a 4,5:1 ou mais
        ax.plot([x0, x0, x0 + TW, x0 + TW], [base_y + TH, base_y, base_y, base_y + TH], color=vs.INK2, lw=1.3,
                zorder=5, solid_capstyle="round")
        ax.text(cx, base_y + lvl * unit_h * 0.58, f"{lvl:.2f}", ha="center", va="center", fontsize=11,
                fontweight="bold", color="white", zorder=6)
        ax.text(cx, base_y + lvl * unit_h * 0.58 - 2.2, "pool", ha="center", va="center", fontsize=6.6, color="white", zorder=6)
        if cx != centres[0]:
            ax.plot([x0 - 1.0, x0 + TW + 1.0], [base_y + unit_h] * 2, color=vs.INK2, lw=0.8, ls=(0, (3, 2)), zorder=6)
        wp = 100.0 * kp * w_rate
        y_tip = base_y + TH - 0.4
        vm.block_arrow(ax, x0 + 3.0, y_tip + 7.2, x0 + 3.0, y_tip, wp, col, head=1.8, z=7)
        t = ax.text(x0 + 5.6, y_tip + 5.6, tr(lang, f"production {100 * kp:.1f} %/h", f"produção {100 * kp:.1f} %/h"),
                    fontsize=6.6, color=vs.INK, va="center", ha="left")
        vm.halo(t)
        wc = 100.0 * kc * w_rate
        vm.block_arrow(ax, x0 + TW, base_y + 1.7, x0 + TW + 6.2, base_y + 1.7, wc, vs.INK2, head=1.8, z=7)
        t = ax.text(x0 + TW + 0.8, base_y + 5.7, tr(lang, f"clearance\n{100 * kc:.1f} %/h", f"depuração\n{100 * kc:.1f} %/h"),
                    fontsize=6.6, color=vs.INK, va="center", ha="left", linespacing=1.1)
        vm.halo(t)
        ax.text(cx, top_h - 0.3, name, ha="center", va="top", fontsize=7.4, fontweight="bold", color=vs.INK)
    ax.text(centres[2], 0.0, tr(lang, "what a production-side mimic would restore", "o que um mimético do lado da produção restauraria"),
            ha="center", va="bottom", fontsize=6.6, color=vs.TIER_RAMP[0], style="italic")
    ax.text(centres[1] + TW / 2 + 1.2, base_y + unit_h + 0.6, tr(lang, "control\nlevel", "nível do\ncontrole"),
            fontsize=6.6, color=vs.INK2, ha="left", va="bottom", linespacing=1.1)
    # the fibril the pool feeds: a render of the deposited structure, as illustration
    FIB_X = 0.840
    p5 = d["prov"]["5OQV"]
    vm.place_image_fit(fig, f"{STRUCT_FIG}/a2_abeta42_fibril_end.png", FIB_X, 1.0 - TOP * 0.50, 54.0, 36.0, FW, FH)
    ax.text(100 * FIB_X, top_h - 0.3, tr(lang, "Aβ42 fibril", "Fibrila de Aβ42"), ha="center", va="top", fontsize=7.4,
            fontweight="bold", color=vs.INK)
    ax.text(100 * FIB_X, top_h - 3.1, structure_line(p5, lang), ha="center", va="top", fontsize=6.6, color=vs.INK2)
    ax.text(100 * FIB_X, 0.0, tr(lang, "the pool feeds aggregation:\nnucleation scales with the\npool squared",
                                 "o pool alimenta a agregação:\na nucleação cresce com o\npool ao quadrado"),
            ha="center", va="bottom", fontsize=6.6, color=vs.INK2, linespacing=1.15)
    vs.panel_letter(ax, "a", x=0.002, y=0.93, size=10)

    # ------------------------------------------------------------------ (b) knockdown against mimic factor
    ax = fig.add_axes([0.075, 0.095, 0.40, 0.335])
    f = np.linspace(1.0, 3.2, 300)
    for rho, lw, ls, lab in ((rs["rho_scan_high"], 1.6, "-", tr(lang, "ρ = {:.0f} (scan, strong)", "ρ = {:.0f} (varredura, forte)")),
                             (rs["rho_scan_mid"], 1.2, "-", tr(lang, "ρ = {:.0f} (midpoint)", "ρ = {:.0f} (ponto médio)")),
                             (rs["rho_scan_low"], 1.2, "--", tr(lang, "ρ = {:.1f} (scan, weak)", "ρ = {:.1f} (varredura, fraca)")),
                             (1.0, 1.2, ":", tr(lang, "ρ = {:.0f}", "ρ = {:.0f}"))):
        ax.plot(f, 100 * knockdown_at(f, rho), color=vs.AD, lw=lw, ls=ls, label=lab.format(rho))
    ax.fill_between(f, 100 * kappa, 100 * measured, color=vm.blend(vs.AD_SOFT, 0.65), lw=0, zorder=0.4)
    ax.axhline(100 * kappa, color=vs.TIER_RAMP[1], lw=1.0)
    ax.axhline(100 * measured, color=vs.INK2, lw=1.0, ls="--")
    ax.text(3.17, 100 * kappa + 1.3, tr(lang, f"required {100 * kappa:.1f}%", f"exigida {100 * kappa:.1f}%"),
            ha="right", va="bottom", fontsize=6.8, color=vs.TIER_RAMP[1])
    ax.text(1.04, 100 * measured + 1.3, tr(lang, f"measured in cells ≈{100 * measured:.0f}%", f"medida em células ≈{100 * measured:.0f}%"),
            ha="left", va="bottom", fontsize=6.8, color=vs.INK2)
    ax.set_xlim(1.0, 3.2)
    ax.set_ylim(0, 70)
    ax.set_xlabel(tr(lang, "miR-29 pool, multiple of baseline", "Pool de miR-29, múltiplo do basal"))
    ax.set_ylabel(tr(lang, "BACE1 knockdown (%)", "Queda de BACE1 (%)"))
    ax.set_title(tr(lang, "Mimic needed for the required fall", "Mimético para a queda exigida"), fontsize=8, loc="left")
    ax.legend(frameon=False, fontsize=6.6, loc="lower right", handlelength=1.8, borderaxespad=0.2)
    style_axes(ax)
    label_panel(ax, "b", x=-0.17)

    # ------------------------------------------------------------------ (c) rate ratios, from the measured reaction orders
    ax = fig.add_axes([0.62, 0.095, 0.365, 0.335])
    n_ = nuc["order_secondary_nucleation"]
    items = [(tr(lang, "Monomer\npool", "Pool de\nmonômero"), ratio, "1"),
             (tr(lang, "Elongation", "Alongamento"), nuc["elongation_rate_ratio_AD_over_control"], "1"),
             (tr(lang, "Nucleation", "Nucleação"), nuc["secondary_nucleation_rate_ratio_AD_over_control"], f"{n_:.0f}")]
    ax.bar(range(3), [v - 1 for _, v, _ in items], bottom=1.0, color=[vs.TIER_RAMP[1], vs.AD, vs.AD], width=0.62)
    ax.axhline(1.0, color=vs.INK2, lw=0.8)
    for i, (_, v, order) in enumerate(items):
        ax.text(i, v + 0.04, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5, color=vs.INK)
        ax.text(i, 1.0 - 0.06, tr(lang, f"order {order}", f"ordem {order}") if i else "", ha="center", va="top",
                fontsize=6.8, color=vs.INK2)
    ax.set_xticks(range(3))
    ax.set_xticklabels([x[0] for x in items], fontsize=6.8)
    ax.set_ylim(0.8, 2.3)
    ax.set_ylabel(tr(lang, "AD relative to control", "DA em relação ao controle"))
    ax.set_title(tr(lang, "Rates at the same pool", "Taxas no mesmo pool"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "c", x=-0.17)

    fig.savefig(fig_path(FIG, "mech_alzheimer", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 4
def fig_parkinson(d, lang):
    ode, tgt = d["ode"], d["tgt"]
    sat = ode["alpha_synuclein_saturation"]
    chain = ode["mir7_knockdown_to_elongation"]
    rs = tgt["repression_strength"]
    C0 = sat["alpha_synuclein_uM_in_bouton"]
    mhs = sorted(p["m_half_uM"] for p in sat["points"])
    kd_a, kd_b = kv(d, "K058"), kv(d, "K063")
    FW, FH = 174.0, 118.0

    fig = plt.figure(figsize=(FW * MM, FH * MM))

    # ------------------------------------------------------------------ (a) saturation curve
    AX_B, AX_H = 0.585, 0.355
    ax = fig.add_axes([0.075, AX_B, 0.40, AX_H])
    S = np.linspace(0, 120, 400)
    cols = (vs.PD, "#c04a1c")
    for mh, c in zip(mhs, cols):
        ax.plot(S, S / (mh + S), color=c, lw=1.4, label=f"m½ {mh:.1f} µM")
    ax.fill_between(S[S <= C0], 0, (S / (mhs[0] + S))[S <= C0], color=vm.blend(vs.PD_SOFT, 0.8), lw=0, zorder=0.4)
    for mh, c in zip(mhs, cols):
        for kd, mk in ((kd_a, "o"), (kd_b, "s")):
            s1 = C0 * (1 - kd)
            ax.plot([s1], [s1 / (mh + s1)], mk, ms=4.2, color=c, mec="white", mew=0.5, zorder=4)
    ax.axvline(C0, color=vs.INK2, lw=0.9, ls="--")
    t = ax.text(C0 + 2, 0.76, tr(lang, f"{C0:.1f} µM in a bouton", f"{C0:.1f} µM em um botão"), fontsize=6.8, color=vs.INK2,
                va="top", ha="left")
    vm.halo(t)
    for kd, y_txt in ((kd_a, 0.20), (kd_b, 0.115)):
        s1 = C0 * (1 - kd)
        ax.annotate(f"−{100 * kd:.0f}%", (s1, s1 / (mhs[-1] + s1)), xytext=(27, y_txt), fontsize=6.6, color=vs.INK2,
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color=vs.INK3, lw=0.6, shrinkA=1, shrinkB=2))
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel(tr(lang, "α-synuclein monomer (µM)", "Monômero de α-sinucleína (µM)"))
    ax.set_ylabel(tr(lang, "Elongation rate, fraction of maximum", "Taxa de alongamento, fração do máximo"))
    ax.legend(frameon=False, fontsize=6.8, loc="lower right", handlelength=1.6)
    ax.set_title(tr(lang, "Where the protein sits", "Onde a proteína está"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "a", x=-0.17)

    # the fibril the monomers join: render of the deposited structure, with the next layer drawn translucent
    p6 = d["prov"]["6CU7"]
    BX0, BW = 0.52, 0.46
    tax = fig.add_axes([BX0, AX_B, BW, AX_H], zorder=10)
    tax.patch.set_alpha(0.0)
    tax.set_xlim(0, 1)
    tax.set_ylim(0, 1)
    tax.axis("off")
    img_cy = AX_B + AX_H * 0.385
    vm.place_image_fit(fig, f"{STRUCT_FIG}/a2_asyn_fibril_end.png", BX0 + BW / 2, img_cy, 76.0, 25.0, FW, FH, zorder=5)
    tax.text(0.5, 0.995, tr(lang, "Elongation: monomers join the growing end", "Alongamento: monômeros entram na ponta que cresce"),
             ha="center", va="top", fontsize=8, fontweight="bold", color=vs.INK)
    tax.text(0.5, 0.905, tr(lang, f"α-synuclein fibril, {structure_line(p6, lang)}", f"Fibrila de α-sinucleína, {structure_line(p6, lang)}"),
             ha="center", va="top", fontsize=6.6, color=vs.INK2)
    tax.annotate("", xy=(0.50, 0.705), xytext=(0.50, 0.845), arrowprops=dict(arrowstyle="-|>", color=vs.PD, lw=1.8, mutation_scale=11,
                                                                          shrinkA=0, shrinkB=0))
    tax.text(0.525, 0.775, tr(lang, "S, the free monomer", "S, o monômero livre"), ha="left", va="center", fontsize=6.8, color=vs.PD_TEXT,
             fontweight="bold")
    tax.text(0.5, 0.01, tr(lang, "translucent layer: where the next one of the same fibril sits",
                           "camada translúcida: onde fica a próxima da mesma fibrila"),
             ha="center", va="bottom", fontsize=6.6, color=vs.INK2, style="italic")

    # ------------------------------------------------------------------ (b) transmission
    ax = fig.add_axes([0.075, 0.09, 0.40, 0.38])
    x = np.linspace(0, 0.6, 200)
    ax.plot(100 * x, 100 * x, color=vs.INK3, lw=0.9, ls="--")
    ax.text(57, 59, tr(lang, "proportional", "proporcional"), ha="right", va="bottom", fontsize=6.8, color=vs.NOTE)
    lo_curve = [1 - ((C0 * (1 - k)) / (mhs[0] + C0 * (1 - k))) / (C0 / (mhs[0] + C0)) for k in x]
    hi_curve = [1 - ((C0 * (1 - k)) / (mhs[-1] + C0 * (1 - k))) / (C0 / (mhs[-1] + C0)) for k in x]
    ax.fill_between(100 * x, 100 * np.array(lo_curve), 100 * np.array(hi_curve), color=vm.blend(vs.PD, 0.35), lw=0, zorder=0.4)
    ax.plot(100 * x, 100 * np.array(hi_curve), color=vs.PD, lw=1.3)
    for lab, kd in (("a", kd_a), ("b", kd_b)):
        cs = [c["elongation_rate_reduction"] for c in chain["cases"] if abs(c["knockdown"] - kd) < 1e-12]
        ax.plot([100 * kd] * 2, [100 * min(cs)] * 2, "o", color=vs.PD, ms=4, mec="white", mew=0.5, zorder=4)
        ax.annotate(f"{100 * min(cs):.1f}–{100 * max(cs):.1f}%", (100 * kd, 100 * min(cs)),
                    xytext=(100 * kd + 3, 100 * min(cs) - 7), fontsize=6.8, color=vs.INK,
                    arrowprops=dict(arrowstyle="-", color=vs.INK3, lw=0.6))
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 62)
    ax.set_xlabel(tr(lang, "Fall in α-synuclein (%)", "Queda de α-sinucleína (%)"))
    ax.set_ylabel(tr(lang, "Fall in elongation rate (%)", "Queda da taxa de alongamento (%)"))
    ax.set_title(tr(lang, "Less than proportional", "Menos que proporcional"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "b", x=-0.17)

    # ------------------------------------------------------------------ (c) mimic factor against rho, both axes
    ax = fig.add_axes([0.60, 0.09, 0.385, 0.38])
    rho = np.geomspace(0.01, 1000, 400)
    pd_ = rs["PD"]
    ad_ = rs["AD"]
    kappa_ad = d["ode"]["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
    # the two declared scans of rho as opaque bands; where they overlap, the colour that the two translucent bands would make
    edges = sorted({pd_["rho_scan_low"], pd_["rho_scan_high"], ad_["rho_scan_low"], ad_["rho_scan_high"]})
    pd_col, ad_col = vm.blend(vs.PD_SOFT, 0.7), vm.blend(vs.AD_SOFT, 0.7)
    both_col = vm.blend(vs.AD_SOFT, 0.7, over=pd_col)
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        mid_ = math.sqrt(lo_ * hi_)
        in_pd = pd_["rho_scan_low"] <= mid_ <= pd_["rho_scan_high"]
        in_ad = ad_["rho_scan_low"] <= mid_ <= ad_["rho_scan_high"]
        if in_pd or in_ad:
            ax.axvspan(lo_, hi_, color=both_col if (in_pd and in_ad) else (pd_col if in_pd else ad_col), lw=0, zorder=0.4)
    ax.plot(rho, fold_for(kd_a, rho), color=vs.PD, lw=1.4, label=tr(lang, f"PD, {100 * kd_a:.0f}%", f"DP, {100 * kd_a:.0f}%"))
    ax.plot(rho, fold_for(kd_b, rho), color=vs.PD, lw=1.4, ls="--", label=tr(lang, f"PD, {100 * kd_b:.0f}%", f"DP, {100 * kd_b:.0f}%"))
    ax.plot(rho, fold_for(kappa_ad, rho), color=vs.AD, lw=1.4, label=tr(lang, f"AD, {100 * kappa_ad:.1f}%", f"DA, {100 * kappa_ad:.1f}%"))
    for kap, c in ((kd_a, vs.PD), (kd_b, vs.PD), (kappa_ad, vs.AD)):
        ax.axhline(1.0 / (1.0 - kap), color=c, lw=0.6, ls=":")
    ax.axvline(pd_["rho_floor_K065"], color=vs.INK2, lw=0.9, ls=":")
    ax.text(pd_["rho_floor_K065"] * 0.9, 70, tr(lang, "floor from K065", "piso do K065"), ha="right", va="top", fontsize=6.8,
            color=vs.INK2, rotation=90)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 1000)
    ax.set_ylim(1.0, 100)
    ax.set_yticks([1, 2, 5, 10, 20, 50, 100])
    ax.set_yticklabels(["1", "2", "5", "10", "20", "50", "100"])
    ax.set_xlabel(tr(lang, "ρ, miRNA-driven over intrinsic mRNA decay", "ρ, decaimento do mRNA pelo miRNA sobre o intrínseco"))
    ax.set_ylabel(tr(lang, "Mimic factor required", "Fator de mimético exigido"))
    ax.legend(frameon=False, fontsize=6.8, loc="lower left", handlelength=1.8, ncol=1, borderaxespad=0.3)
    ax.set_title(tr(lang, "How much mimic, by repression strength", "Quanto mimético, por força de repressão"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "c", x=-0.17)

    fig.savefig(fig_path(FIG, "mech_parkinson", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 5
def fig_turnover(d, lang):
    dos = d["dos"]
    dec = dos["decay_constants"]
    dd = dos["daily_dosing_comparison"]
    order = ["miR-7", "miR-7 upper bound", "miR-29b", "miR-29c", "median miRNA"]
    FW, FH = 174.0, 118.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))

    # ------------------------------------------------------------------ (a) the ladder of time scales
    hl = lambda pid: LN2 / kv(d, pid)
    sp_hl = {k: LN2 / dec[k]["decay_per_hour"] for k in order}
    rows = [  # group, label EN, label PT, half-life h, colour, filled
        ("m", "miR-7", "miR-7", sp_hl["miR-7"], vs.PD, True),
        ("m", "miR-7, upper bound", "miR-7, limite superior", sp_hl["miR-7 upper bound"], vs.PD, False),
        ("m", "miR-29b", "miR-29b", sp_hl["miR-29b"], vs.AD, True),
        ("m", "miR-29c", "miR-29c", sp_hl["miR-29c"], vs.AD, False),
        ("m", "median miRNA", "miRNA mediano", sp_hl["median miRNA"], vs.INK2, True),
        ("r", "SNCA mRNA", "mRNA de SNCA", hl("K066"), vs.PD, True),
        ("r", "BACE1 mRNA", "mRNA de BACE1", hl("K068"), vs.AD, True),
        ("p", "α-synuclein", "α-sinucleína", hl("K026"), vs.PD, True),
        ("p", "BACE1", "BACE1", hl("K025"), vs.AD, True),
        ("p", "brain proteome, median", "proteoma do cérebro, mediana", hl("K029"), vs.INK2, True),
        ("a", "Aβ42 clearance, controls", "depuração de Aβ42, controles", hl("K003"), vs.AD, False),
        ("a", "Aβ42 clearance, AD", "depuração de Aβ42, DA", hl("K004"), vs.AD, True),
    ]
    group_name = {"m": tr(lang, "miRNA", "miRNA"), "r": "mRNA", "p": tr(lang, "Protein", "Proteína"),
                  "a": "Aβ42"}
    ax = fig.add_axes([0.235, 0.575, 0.585, 0.385])
    ys, y, last = [], 0.0, None
    for g, *_ in rows:
        if last is not None and g != last:
            y += 0.6
        ys.append(y)
        y += 1.0
        last = g
    ys = np.array(ys)
    groups = OrderedDict()
    for yv, r in zip(ys, rows):
        groups.setdefault(r[0], []).append(yv)
    for gi, (g, yy) in enumerate(groups.items()):
        ax.axhspan(min(yy) - 0.5, max(yy) + 0.5, color=vs.BAND if gi % 2 == 0 else "white", lw=0, zorder=0)
        ax.text(930, (min(yy) + max(yy)) / 2, group_name[g], ha="right", va="center", fontsize=7.4, fontweight="bold", color=vs.NOTE, zorder=1)
    for T in (8, 12, 24, 48, 168):
        ax.axvline(T, color=vs.INK3, lw=0.7, ls=(0, (3, 2)), zorder=1)
    tx = ax.get_xaxis_transform()
    for T, lab in ((8, "8 h"), (12, "12 h"), (24, tr(lang, "24 h\ndaily", "24 h\ndiária")), (48, "48 h"), (168, tr(lang, "168 h\nweekly", "168 h\nsemanal"))):
        ax.plot([T, T], [1.0, 1.025], transform=tx, color=vs.INK2, lw=0.7, clip_on=False, zorder=2)
        ax.text(T, 1.03, lab, transform=tx, ha="center", va="bottom", fontsize=6.6, color=vs.INK2, linespacing=1.1)
    for yv, (g, en, pt, h, c, filled) in zip(ys, rows):
        ax.plot([h], [yv], "o", ms=6.2, mfc=c if filled else "white", mec=c, mew=1.2, zorder=3)
        t = ax.text(h * 1.17, yv, f"{h:.1f} h" if h < 100 else f"{h:.0f} h", va="center", ha="left", fontsize=6.8, color=vs.INK, zorder=4)
        vm.halo(t, 2.4)
    ax.set_yticks(ys)
    ax.set_yticklabels([tr(lang, r[1], r[2]) for r in rows], fontsize=6.8)
    ax.set_xscale("log")
    ax.set_xlim(1, 1000)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1 h", "10 h", "100 h", "1000 h"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(ys.max() + 0.7, -0.7)
    ax.set_title(tr(lang, "Half-lives of the steps, against the dosing intervals", "Meias-vidas dos passos, frente aos intervalos de dose"),
                 fontsize=8, loc="left", pad=22)
    label_panel(ax, "a", x=-0.43, y=1.2)
    ax2 = fig.add_axes([0.835, 0.575, 0.15, 0.385], sharey=ax)
    ax2.set_xlim(0, 1)
    ax2.axis("off")
    ax2.text(0.0, -0.7, tr(lang, "half-lives per 24 h", "meias-vidas em 24 h"), fontsize=6.6, color=vs.INK2, ha="left", va="bottom",
             fontweight="bold")
    for yv, (g, en, pt, h, c, filled) in zip(ys, rows):
        ax2.text(0.0, yv, f"{24.0 / h:.1f}", va="center", ha="left", fontsize=6.8, color=vs.INK)

    # ------------------------------------------------------------------ (b) washout
    ax = fig.add_axes([0.075, 0.085, 0.265, 0.34])
    t = np.linspace(0, 120, 400)
    for sp in order:
        c, ls, en, pt = SP_COLOUR[sp]
        ax.plot(t, 1 + 9 * np.exp(-dec[sp]["decay_per_hour"] * t), color=c, ls=ls, lw=1.3, label=tr(lang, en, pt))
    ax.axhline(1.1, color=vs.INK3, lw=0.8, ls="--")
    ax.text(119, 1.13, tr(lang, "within 10% of baseline", "a 10% do basal"), ha="right", va="bottom", fontsize=6.6, color=vs.INK2)
    ax.set_yscale("log")
    ax.set_ylim(1.0, 11)
    ax.set_yticks([1, 2, 5, 10])
    ax.set_yticklabels(["1", "2", "5", "10"])
    ax.set_xlim(0, 120)
    ax.set_xlabel(tr(lang, "Time after a 10-fold bolus (h)", "Tempo após um bolus de 10x (h)"))
    ax.set_ylabel(tr(lang, "Level, multiple of baseline", "Nível, múltiplo do basal"))
    ax.set_title(tr(lang, "Washout", "Eliminação"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "b", x=-0.27)

    # ------------------------------------------------------------------ (c) peak to average
    ax = fig.add_axes([0.415, 0.085, 0.265, 0.34])
    T = np.geomspace(0.5, 72, 300)
    for sp in order:
        c, ls, en, pt = SP_COLOUR[sp]
        x = dec[sp]["decay_per_hour"] * T
        ax.plot(T, x / (1 - np.exp(-x)), color=c, ls=ls, lw=1.3, label=tr(lang, en, pt))
    ax.axvline(24, color=vs.INK3, lw=0.8, ls="--")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", handlelength=1.4, borderaxespad=0.12, labelspacing=0.3)
    ax.annotate(tr(lang, "daily", "diário"), xy=(24, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 1.5),
                textcoords="offset points", ha="center", va="bottom", fontsize=6.8, color=vs.NOTE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 72)
    ax.set_ylim(1, 40)
    ax.set_xticks([1, 2, 4, 8, 24, 48])
    ax.set_xticklabels(["1", "2", "4", "8", "24", "48"])
    ax.set_yticks([1, 2, 5, 10, 20])
    ax.set_yticklabels(["1", "2", "5", "10", "20"])
    ax.set_xlabel(tr(lang, "Dosing interval (h)", "Intervalo entre doses (h)"))
    ax.set_ylabel(tr(lang, "Peak / average", "Pico / média"))
    ax.set_title(tr(lang, "Cost of instability", "Custo da instabilidade"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "c", x=-0.27)

    # ------------------------------------------------------------------ (d) stabilisation needed
    ax = fig.add_axes([0.775, 0.085, 0.21, 0.34])
    sps = ["miR-7", "miR-7 upper bound", "miR-29b", "miR-29c"]
    vals = [dd[s]["fold_stabilisation_required"] for s in sps]
    cols = [SP_COLOUR[s][0] for s in sps]
    y = np.arange(len(sps))
    bars = ax.barh(y, vals, color=cols, height=0.6)
    for b, s in zip(bars, sps):
        if SP_COLOUR[s][1] == "--":
            b.set_facecolor(vm.blend(b.get_facecolor()[:3], 0.6))
    for yi, v in zip(y, vals):
        ax.text(v + 0.4, yi, f"{v:.1f}×", va="center", fontsize=7, color=vs.INK)
    ax.set_yticks(y)
    short = {"miR-7": "miR-7", "miR-7 upper bound": "miR-7 (5 h)", "miR-29b": "miR-29b", "miR-29c": "miR-29c"}
    ax.set_yticklabels([short[s_] for s_ in sps], fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlim(0, 26)
    ax.set_xlabel(tr(lang, "Fold increase in half-life", "Aumento da meia-vida (vezes)"))
    ax.set_title(tr(lang, "Stabilisation for daily dosing", "Estabilização para dose diária"), fontsize=8, loc="left")
    style_axes(ax, grid=False)
    label_panel(ax, "d", x=-0.45)

    fig.savefig(fig_path(FIG, "mech_turnover", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 6
def trajectory(cas, dose, T, dx, cycles=2, per=240):
    """EN/PT: E, mu, pi over `cycles` periodic cycles of repeated dosing | E, mu, pi ao longo de ciclos periodicos."""
    e_peak = dose / (1.0 - math.exp(-dx * T))
    x0 = list(cas.periodic_state(e_peak, T))
    ts, ys = [], []
    for c in range(cycles):
        sol = solve_ivp(cas.rhs, (0.0, T), [e_peak] + x0 + [0.0] * cas.n_int, method="LSODA",
                        rtol=1e-10, atol=1e-13, t_eval=np.linspace(0.0, T, per))
        ts.append(sol.t + c * T)
        ys.append(sol.y[:3 + cas.L])
        x0 = list(sol.y[1:1 + cas.n_lin, -1])
    return np.concatenate(ts), np.concatenate(ys, axis=1)


def scen(rows, axis, rho_case, sp, T, kap, prot=("BACE1 in", "alpha-syn")):
    out = [r for r in rows if r["axis"] == axis and r["rho_case"] == rho_case and r["species"] == sp
           and float(r["interval_h"]) == T and abs(float(r["target_knockdown"]) - kap) < 1e-9
           and r["protein_case"].startswith(prot)]
    if len(out) != 1:
        sys.exit(f"EN/PT: expected one scenario for {axis}/{rho_case}/{sp}/{T}, found {len(out)}")
    return out[0]


def fig_target(d, lang, t36):
    rows = d["rows"]
    meas, table = t36.load_module("ode12", "12_ode_models_calibrated.py").load_measured()
    tgt = d["tgt"]
    kap_a = kv(d, "K058")
    kap_ad = d["ode"]["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
    T24 = 24.0
    mid = tgt["repression_strength"]["PD"]["rho_scan_mid"]
    cap = tgt["dose_cap_in_baseline_pools"]
    intervals = [8.0, 12.0, 24.0, 48.0, 168.0]
    FW, FH = 174.0, 142.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    outer = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.46, left=0.075, right=0.985, top=0.955, bottom=0.075)

    # ------------------------------------------------------------------ (a) time courses
    inner = outer[0, 0].subgridspec(3, 1, hspace=0.18)
    dp = meas["d_aSyn"]
    d_mrna = meas["d_mRNA_SNCA"]
    series = []
    for sp in ("miR-7", "median miRNA"):
        r = scen(rows, "PD", "scan midpoint", sp, T24, kap_a)
        dx = LN2 / float(r["half_life_h"])
        cas = t36.Cascade(d_mrna, mid, dp, dx)
        tt, yy = trajectory(cas, float(r["dose"]), T24, dx)
        series.append((sp, tt, yy))
    names = [tr(lang, "Repressor pool\n(× baseline)", "Pool do repressor\n(× basal)"),
             tr(lang, "mRNA\n(× baseline)", "mRNA\n(× basal)"), tr(lang, "Protein\n(× baseline)", "Proteína\n(× basal)")]
    axs = []
    for k in range(3):
        ax = fig.add_subplot(inner[k])
        for sp, tt, yy in series:
            c, ls, en, pt = SP_COLOUR[sp]
            y = 1 + yy[0] if k == 0 else yy[k]
            ax.plot(tt, y, color=c, ls="-", lw=1.3, label=tr(lang, en, pt))
            if k == 2 and sp == "miR-7":
                ax.fill_between(tt, y, 1.0, color=vm.blend(c, 0.14), lw=0, zorder=0.4)
        if k == 0:
            ax.set_yscale("log")
            ax.set_ylim(1, 100)
            ax.set_yticks([1, 3, 10, 30])
            ax.set_yticklabels(["1", "3", "10", "30"])
        elif k == 1:
            ax.set_ylim(0, 1.15)
        else:
            ax.set_ylim(0.55, 1.05)
        ax.set_ylabel(names[k], fontsize=6.8)
        style_axes(ax)
        if k == 2:
            ax.axhline(1.0 - kap_a, color=vs.INK3, lw=0.8, ls="--")
            ax.text(47.5, 0.565, tr(lang, f"target: {100 * kap_a:.0f}% below baseline, as a cycle mean",
                                                   f"meta: {100 * kap_a:.0f}% abaixo do basal, na média do ciclo"),
                    ha="right", va="bottom", fontsize=6.6, color=vs.INK2)
        if k < 2:
            ax.set_xticklabels([])
        axs.append(ax)
    axs[2].set_xlabel(tr(lang, "Time (h), daily dosing", "Tempo (h), dose diária"))
    axs[0].legend(frameon=False, fontsize=6.6, loc="upper right", handlelength=1.6, borderaxespad=0.1)
    axs[0].set_title(tr(lang, f"PD, {100 * kap_a:.0f}% target, ρ = {mid:.1f}", f"DP, meta de {100 * kap_a:.0f}%, ρ = {mid:.1f}"),
                     fontsize=8, loc="left")
    label_panel(axs[0], "a")

    # ------------------------------------------------------------------ (b) and (c): peak against interval, ribbons over the scan
    def peak_panel(ax, axis, items, kap, title, letter, floor_case=None, ymax=None):
        for sp, colour, ls in items:
            lo, mid_, hi = [], [], []
            for T in intervals:
                rl = scen(rows, axis, "scan low", sp, T, kap)
                rm = scen(rows, axis, "scan midpoint", sp, T, kap)
                rh = scen(rows, axis, "scan high", sp, T, kap)
                f_ = lambda r: float(r["peak_pool_over_baseline"]) if r["reachable"] == "True" else None
                lo.append(f_(rl)); mid_.append(f_(rm)); hi.append(f_(rh))
            ok = [i for i in range(len(intervals)) if lo[i] is not None and hi[i] is not None]
            if ok:
                ax.fill_between([intervals[i] for i in ok], [hi[i] for i in ok], [lo[i] for i in ok], color=vm.blend(colour, 0.16), lw=0, zorder=0.4)
            okm = [i for i in range(len(intervals)) if mid_[i] is not None]
            ax.plot([intervals[i] for i in okm], [mid_[i] for i in okm], color=colour, ls=ls, lw=1.4, marker="o", ms=3.4,
                    mec="white", mew=0.4, zorder=3)
            if floor_case and sp == items[0][0]:
                fl = [scen(rows, axis, floor_case, sp, T, kap) for T in intervals]
                okf = [(T, float(r["peak_pool_over_baseline"])) for T, r in zip(intervals, fl) if r["reachable"] == "True"]
                ax.plot([a for a, _ in okf], [b for _, b in okf], color=colour, ls=":", lw=1.2, zorder=2)
            for i, T in enumerate(intervals):
                if mid_[i] is None:
                    ax.plot([T], [cap], marker="o", ms=4.2, mfc="white", mec=colour, mew=0.9, ls="none", zorder=3)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xticks(intervals)
        ax.set_xticklabels(["8", "12", "24", "48", "168"])
        ax.set_ylim(1, ymax or cap * 2)
        if ymax:
            ax.set_yticks([1, 10, 100])
            ax.set_yticklabels(["1", "10", "100"])
        else:
            ax.set_yticks([1, 10, 100, 1000, 10000])
            ax.set_yticklabels(["1", "10", "100", "1000", "10⁴"])
        ax.set_xlabel(tr(lang, "Dosing interval (h)", "Intervalo entre doses (h)"))
        ax.set_ylabel(tr(lang, "Peak of the repressor pool (× baseline)", "Pico do pool do repressor (× basal)"))
        ax.set_title(title, fontsize=8, loc="left")
        style_axes(ax)
        label_panel(ax, letter)

    ax = fig.add_subplot(outer[0, 1])
    peak_panel(ax, "PD", [("miR-7", vs.PD, "-"), ("miR-7 upper bound", vs.PD, "--"), ("median miRNA", vs.INK2, "-")], kap_a,
               tr(lang, f"PD, {100 * kap_a:.0f}% fall: scan band, midpoint, floor", f"DP, queda de {100 * kap_a:.0f}%: faixa da varredura, ponto médio, piso"),
               "b", floor_case="K065 lower bound")
    ax.text(165, cap * 1.25, tr(lang, "not reached", "não alcançada"), fontsize=6.8, color=vs.INK2, ha="right", va="bottom")
    ax.legend([plt.Line2D([], [], color=SP_COLOUR[s][0] if s != "median miRNA" else vs.INK2, ls=SP_COLOUR[s][1], lw=1.4)
               for s in ("miR-7", "miR-7 upper bound", "median miRNA")],
              [tr(lang, *SP_COLOUR[s][2:]) for s in ("miR-7", "miR-7 upper bound", "median miRNA")],
              frameon=False, fontsize=6.6, loc="upper left", handlelength=1.8, borderaxespad=0.2)
    ax = fig.add_subplot(outer[1, 0])
    peak_panel(ax, "AD", [("miR-29b", vs.AD, "-"), ("miR-29c", vs.AD, "--"), ("median miRNA", vs.TIER_RAMP[0], "-")], kap_ad,
               tr(lang, f"AD, {100 * kap_ad:.1f}% fall: scan band and midpoint", f"DA, queda de {100 * kap_ad:.1f}%: faixa da varredura e ponto médio"),
               "c", ymax=300)
    ax.legend([plt.Line2D([], [], color=c, ls=ls, lw=1.4) for c, ls in ((vs.AD, "-"), (vs.AD, "--"), (vs.TIER_RAMP[0], "-"))],
              [tr(lang, *SP_COLOUR[s][2:]) for s in ("miR-29b", "miR-29c", "median miRNA")],
              frameon=False, fontsize=6.6, loc="upper left", handlelength=1.8, borderaxespad=0.2)

    # ------------------------------------------------------------------ (d) efficiency
    ax = fig.add_subplot(outer[1, 1])
    for axis, sp, kap in (("PD", "miR-7", kap_a), ("PD", "miR-7 upper bound", kap_a), ("PD", "median miRNA", kap_a),
                          ("AD", "miR-29b", kap_ad), ("AD", "miR-29c", kap_ad), ("AD", "median miRNA", kap_ad)):
        rr = [scen(rows, axis, "scan midpoint", sp, T, kap) for T in intervals]
        c = SP_COLOUR[sp][0] if sp != "median miRNA" else (vs.INK2 if axis == "PD" else vs.TIER_RAMP[0])
        ls = SP_COLOUR[sp][1]
        ok = [(T, r) for T, r in zip(intervals, rr) if r["reachable"] == "True"]
        ax.plot([T for T, _ in ok], [float(r["efficiency"]) for _, r in ok], color=c, ls=ls, lw=1.4, marker="o", ms=3.4,
                mec="white", mew=0.4)
        for T, r in zip(intervals, rr):
            if r["reachable"] != "True":
                ax.plot([T], [0.0], marker="o", ms=4.2, mfc="white", mec=c, mew=0.9, ls="none")
    ax.set_xscale("log")
    ax.set_xticks(intervals)
    ax.set_xticklabels(["8", "12", "24", "48", "168"])
    ax.set_ylim(-0.04, 1.04)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel(tr(lang, "Dosing interval (h)", "Intervalo entre doses (h)"))
    ax.set_ylabel(tr(lang, "Efficiency", "Eficiência"))
    ax.set_title(tr(lang, "Closed form against the cycle (midpoint of ρ)", "Forma fechada contra o ciclo (ponto médio de ρ)"), fontsize=8, loc="left")
    style_axes(ax)
    label_panel(ax, "d")
    ax.legend([plt.Line2D([], [], color=SP_COLOUR[s][0] if s != "median miRNA" else vs.INK2, ls=SP_COLOUR[s][1], lw=1.4)
               for s in ("miR-7", "miR-7 upper bound", "median miRNA")] +
              [plt.Line2D([], [], color=vs.AD, ls="-", lw=1.4), plt.Line2D([], [], color=vs.AD, ls="--", lw=1.4),
               plt.Line2D([], [], color=vs.TIER_RAMP[0], ls="-", lw=1.4)],
              [tr(lang, "PD miR-7", "DP miR-7"), tr(lang, "PD miR-7 (5 h)", "DP miR-7 (5 h)"), tr(lang, "PD median", "DP mediano"),
               tr(lang, "AD miR-29b", "DA miR-29b"), tr(lang, "AD miR-29c", "DA miR-29c"), tr(lang, "AD median", "DA mediano")],
              frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=3, handlelength=1.8, columnspacing=1.2,
              labelspacing=0.35)

    fig.savefig(fig_path(FIG, "mech_target_response", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 7
def fine_mesh(Z, t_grid, rho_grid, kt=4, kr=3):
    """
    EN | Display grid for one map: the solved grid, refined by linear interpolation of log10(peak)
         over (log10 interval, log10 rho). The cells that no dose reaches stay marked as on the
         solved grid (nearest cell), so that no value is invented where the solver found none.
    PT | Grade de exibicao de um mapa: a grade resolvida, refinada por interpolacao linear de
         log10(pico) sobre (log10 intervalo, log10 rho). As celulas que nenhuma dose alcanca
         continuam marcadas como na grade resolvida (celula mais proxima), de modo que nenhum
         valor e inventado onde o solver nao achou nenhum.
    """
    from scipy.interpolate import RegularGridInterpolator
    lt, lr = np.log10(t_grid), np.log10(rho_grid)
    unreached = np.isnan(Z)
    zl = np.where(unreached, 3.0, np.log10(np.clip(Z, 1.0, 1000.0)))
    lt_f = np.linspace(lt[0], lt[-1], (len(lt) - 1) * kt + 1)
    lr_f = np.linspace(lr[0], lr[-1], (len(lr) - 1) * kr + 1)
    TT, RR = np.meshgrid(lt_f, lr_f, indexing="ij")
    pts = np.column_stack([TT.ravel(), RR.ravel()])
    lin = RegularGridInterpolator((lt, lr), zl, method="linear")(pts).reshape(TT.shape)
    near = RegularGridInterpolator((lt, lr), unreached.astype(float), method="nearest")(pts).reshape(TT.shape) > 0.5
    return 10 ** lr_f, 10 ** lt_f, np.ma.masked_array(10 ** lin, mask=near), near


def fmt_range(lo, hi):
    """EN/PT: 'a-b x' with one decimal below 10 and none above, as the text quotes it | faixa com uma casa abaixo de 10."""
    if hi < 10.0:
        return f"{lo:.1f}–{hi:.1f}×"
    return f"{lo:.0f}–{hi:.0f}×"


def key_line(fig, x, y, label, colour, ls="-", lw=1.0, casing=False, length=0.030, fontsize=6.6):
    """EN/PT: a short line sample with its label to the right, in figure fractions | amostra de linha com o rotulo."""
    from matplotlib.lines import Line2D
    if casing:
        fig.add_artist(Line2D([x, x + length], [y, y], color="#b9b6ae", lw=lw + 1.6, transform=fig.transFigure, solid_capstyle="butt"))
        fig.add_artist(Line2D([x, x + length], [y, y], color="white", lw=lw + 0.6, transform=fig.transFigure, solid_capstyle="butt"))
    fig.add_artist(Line2D([x, x + length], [y, y], color=colour, lw=lw, ls=ls, transform=fig.transFigure, solid_capstyle="butt"))
    return fig.text(x + length + 0.008, y, label, fontsize=fontsize, color=vs.INK2, va="center", ha="left")


def fig_design(d, lang):
    ds = d["ds"]
    rho_grid = np.array(ds["rho_grid"])
    t_grid = np.array(ds["interval_grid_h"])
    rows = d["ds_rows"]
    cap = d["tgt"]["dose_cap_in_baseline_pools"]
    rs = d["tgt"]["repression_strength"]
    arr = {}
    for r in rows:
        key = (r["axis"], r["species"], round(float(r["target_knockdown"]), 4))
        a = arr.setdefault(key, {"peak": np.full((len(t_grid), len(rho_grid)), np.nan)})
        ti = int(np.argmin(np.abs(t_grid - float(r["interval_h"]))))
        ri = int(np.argmin(np.abs(np.log10(rho_grid) - math.log10(float(r["rho"])))))
        if r["reachable"] == "True":
            a["peak"][ti, ri] = float(r["peak_pool_over_baseline"])
    kap_a, kap_b = kv(d, "K058"), kv(d, "K063")
    kap_ad = d["ode"]["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
    layout = [  # EN/PT: row-wise panels: (axis, species, knockdown, title EN, title PT)
        ("PD", "miR-7", kap_a, f"PD · miR-7 (1.7 h) · {100 * kap_a:.0f}%", f"DP · miR-7 (1.7 h) · {100 * kap_a:.0f}%"),
        ("PD", "miR-7 upper bound", kap_a, f"PD · miR-7 (5 h) · {100 * kap_a:.0f}%", f"DP · miR-7 (5 h) · {100 * kap_a:.0f}%"),
        ("PD", "median miRNA", kap_a, f"PD · median (34 h) · {100 * kap_a:.0f}%", f"DP · mediano (34 h) · {100 * kap_a:.0f}%"),
        ("PD", "miR-7", kap_b, f"PD · miR-7 (1.7 h) · {100 * kap_b:.0f}%", f"DP · miR-7 (1.7 h) · {100 * kap_b:.0f}%"),
        ("PD", "miR-7 upper bound", kap_b, f"PD · miR-7 (5 h) · {100 * kap_b:.0f}%", f"DP · miR-7 (5 h) · {100 * kap_b:.0f}%"),
        ("PD", "median miRNA", kap_b, f"PD · median (34 h) · {100 * kap_b:.0f}%", f"DP · mediano (34 h) · {100 * kap_b:.0f}%"),
        ("AD", "miR-29b", kap_ad, f"AD · miR-29b (7 h) · {100 * kap_ad:.1f}%", f"DA · miR-29b (7 h) · {100 * kap_ad:.1f}%"),
        ("AD", "miR-29c", kap_ad, f"AD · miR-29c (10.6 h) · {100 * kap_ad:.1f}%", f"DA · miR-29c (10.6 h) · {100 * kap_ad:.1f}%"),
        ("AD", "median miRNA", kap_ad, f"AD · median (34 h) · {100 * kap_ad:.1f}%", f"DA · mediano (34 h) · {100 * kap_ad:.1f}%"),
    ]
    FW, FH = 174.0, 178.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    gs = fig.add_gridspec(3, 3, left=0.082, right=0.985, top=0.950, bottom=0.172, wspace=0.34, hspace=0.50)
    levels = [3, 10, 30, 100, 300]
    norm = LogNorm(vmin=1.0, vmax=1000.0)
    letters = "abcdefghi"
    rows_t = d["rows"]
    plt.rcParams["hatch.linewidth"] = 0.4
    plt.rcParams["hatch.color"] = "#8a8883"
    for k, (axis, sp, kap, ten, ptt) in enumerate(layout):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        key = (axis, sp, round(kap, 4))
        Z = arr[key]["peak"]
        cmap = vm.PD_MAP if axis == "PD" else vm.AD_MAP
        Rf, Tf, Zf, near = fine_mesh(Z, t_grid, rho_grid)
        ax.set_facecolor("#e6e4df")
        ax.pcolormesh(Rf, Tf, Zf, cmap=cmap, norm=norm, shading="nearest", rasterized=True, zorder=1)
        if near.any():
            ax.contourf(Rf, Tf, near.astype(float), levels=[0.5, 1.5], colors=["#e6e4df"], zorder=1.5)
            ax.contourf(Rf, Tf, near.astype(float), levels=[0.5, 1.5], colors="none", hatches=["//////"], zorder=2)
        cs = ax.contour(Rf, Tf, Zf, levels=levels, colors=[vs.INK], linewidths=0.5, zorder=4)
        if k % 3 == 0:
            spots = []      # EN/PT: one spot per level, on the stretch nearest rho = 0.03 and clear of the 24 h and 168 h lines | um ponto por nivel
            for li_, lv_ in enumerate(cs.levels):
                if lv_ not in (10, 100):
                    continue
                best = None
                for seg_ in cs.allsegs[li_]:
                    for x_, y_ in seg_:
                        if min(abs(math.log(y_ / 24.0)), abs(math.log(y_ / 168.0))) < 0.25 or not (3.0 < y_ < 150.0):
                            continue
                        dx_ = abs(math.log10(x_) - math.log10(0.03))
                        if best is None or dx_ < best[0]:
                            best = (dx_, x_, y_)
                if best:
                    spots.append((best[1], best[2]))
            if spots:
                lab = ax.clabel(cs, levels=[10, 100], fmt=lambda v: f"{v:g}×", fontsize=6.5, inline=True, inline_spacing=2, manual=spots)
                for t_ in lab:
                    vm.halo(t_, 1.6)
        for T_ in (24.0, 168.0):
            ax.axhline(T_, color="white", lw=1.5, zorder=5)
            ax.axhline(T_, color=vs.INK, lw=0.55, ls=(0, (3, 2)), zorder=5)
        e = rs[axis]
        for x_ in (e["rho_scan_low"], e["rho_scan_high"]):
            ax.axvline(x_, color=vs.INK, lw=0.7, ls=(0, (1.5, 2.0)), zorder=5)
        if axis == "PD":
            ax.axvline(e["rho_floor_K065"], color=vs.INK, lw=1.1, zorder=5)
        # the daily read-out inside the declared scan of rho: the figures that the text quotes
        reads = []
        for case in ("scan low", "scan midpoint", "scan high"):
            r = scen(rows_t, axis, case, sp, 24.0, kap)
            if r["reachable"] == "True":
                reads.append(float(r["peak_pool_over_baseline"]))
        if len(reads) == 3:
            xm = math.sqrt(e["rho_scan_low"] * e["rho_scan_high"])
            t = ax.text(xm, 24.0 * 1.13, fmt_range(min(reads), max(reads)), fontsize=6.6, fontweight="bold", color=vs.INK,
                        ha="center", va="bottom", zorder=7)
            vm.halo(t, 2.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.01, 1e4)
        ax.set_ylim(t_grid.min(), t_grid.max())
        ax.set_xticks([0.01, 1, 100, 10000])
        ax.set_xticklabels(["0.01", "1", "100", "10⁴"])
        ax.set_yticks([2, 4, 8, 24, 48, 168])
        ax.set_yticklabels(["2", "4", "8", "24", "48", "168"])
        ax.minorticks_off()
        ax.spines[["top", "right"]].set_visible(True)
        for s_ in ax.spines.values():
            s_.set_color(vs.RULE)
        ax.set_title(tr(lang, ten, ptt), fontsize=7.2, loc="left", pad=3)
        if k % 3 == 0:
            ax.set_ylabel(tr(lang, "Dosing interval (h)", "Intervalo entre doses (h)"), fontsize=7)
        if k // 3 == 2:
            ax.set_xlabel("ρ", fontsize=8, labelpad=1)
        vs.panel_letter(ax, letters[k], x=-0.30, y=1.07, size=9)

    # ------------------------------------------------------------------ colour bars and key
    g = np.logspace(0, 3, 256).reshape(1, -1)
    for x0, cm, lab in ((0.082, vm.PD_MAP, tr(lang, "PD · peak of the pool needed", "DP · pico do pool exigido")),
                        (0.392, vm.AD_MAP, tr(lang, "AD · peak of the pool needed", "DA · pico do pool exigido"))):
        cax = fig.add_axes([x0, 0.092, 0.240, 0.013])
        cax.imshow(g, aspect="auto", cmap=cm, norm=norm, extent=(0, 3, 0, 1))
        cax.set_yticks([])
        cax.set_xticks([0, math.log10(3), 1, 2, 3])
        cax.set_xticklabels(["1×", "3×", "10×", "100×", "≥ 1000×"], fontsize=6.5)
        cax.tick_params(length=2, pad=1.5)
        for s_ in cax.spines.values():
            s_.set_color(vs.RULE)
        cax.set_title(lab, fontsize=6.8, loc="left", pad=2)
    hax = fig.add_axes([0.700, 0.089, 0.030, 0.021])
    hax.set_xlim(0, 1)
    hax.set_ylim(0, 1)
    hax.add_patch(Rectangle((0, 0), 1, 1, fc="#e6e4df", ec="#8a8883", lw=0.4, hatch="//////"))
    hax.axis("off")
    fig.text(0.740, 0.0995, tr(lang, f"no dose up to {cap:,.0f} × baseline\nreaches the target", f"nenhuma dose até {cap:,.0f} × o basal\nalcança a meta"),
             fontsize=6.6, color=vs.INK2, va="center", ha="left", linespacing=1.15)
    y1, y2 = 0.050, 0.022
    key_line(fig, 0.082, y1, tr(lang, "24 h (daily) and 168 h (weekly)", "24 h (diária) e 168 h (semanal)"), vs.INK, ls=(0, (3, 2)), lw=0.55, casing=True)
    key_line(fig, 0.392, y1, tr(lang, "limits of the declared scan of ρ", "limites da varredura declarada de ρ"), vs.INK, ls=(0, (1.5, 2.0)), lw=0.7)
    key_line(fig, 0.700, y1, tr(lang, "PD floor from K065", "piso do K065 (DP)"), vs.INK, ls="-", lw=1.1)
    fig.text(0.082, y2, tr(lang, "Black lines: 3×, 10×, 30×, 100× and 300×. Bold figure on the daily line: range of the peak across the scan.",
                           "Linhas pretas: 3×, 10×, 30×, 100× e 300×. Número em negrito na linha diária: faixa do pico na varredura."),
             fontsize=6.5, color=vs.INK2, va="center", ha="left")
    fig.savefig(fig_path(FIG, "mech_design_space", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 2
def fig_provenance(d, lang):
    rows = d["kin_rows"]
    sys_of, blk_of = pc.classes(rows)
    sys_keys = [k for k, *_ in pc.SYSTEMS]
    blk_keys = [k for k, *_ in pc.BLOCKS]
    sys_idx = {k: i for i, k in enumerate(sys_keys)}
    blk_idx = {k: i for i, k in enumerate(blk_keys)}
    src_of = {r["param_id"]: pc.source_of(r) for r in rows}
    axis_of = {r["param_id"]: pc.AXIS[r["axis"]] for r in rows}
    sources = sorted(set(src_of.values()), key=lambda s: (s == ("", ""), s[1], s[0]))
    n_rows = len(rows)

    # barycentre of each source over the systems it feeds: the order that keeps the ribbons from crossing
    def bary(s):
        ids = [p for p, v in src_of.items() if v == s]
        return sum(sys_idx[sys_of[p]] for p in ids) / len(ids), -len(ids)
    sources.sort(key=lambda s: (s == ("", ""), bary(s)))
    # one flow per (source, system, axis) and per (system, block, axis)
    f1, f2 = {}, {}
    for r in rows:
        p = r["param_id"]
        f1[(src_of[p], sys_of[p], axis_of[p])] = f1.get((src_of[p], sys_of[p], axis_of[p]), 0) + 1
        f2[(sys_of[p], blk_of[p], axis_of[p])] = f2.get((sys_of[p], blk_of[p], axis_of[p]), 0) + 1
    AX_ORDER = {"AD": 0, "both": 1, "PD": 2}

    FW, FH = 174.0, 132.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    top = 100.0 * FH / FW
    ax.set_xlim(0, 100)
    ax.set_ylim(0, top)
    ax.axis("off")
    U = 0.80                                  # height of one table row, in axis units
    y_hi, g_src = top - 8.2, 0.70
    count_src = {s: sum(1 for v in src_of.values() if v == s) for s in sources}
    total = sum(count_src.values()) * U + g_src * (len(sources) - 1)
    y_lo = y_hi - total

    def stack(keys, counts, y_top, y_bot):
        """EN/PT: vertical positions of a column of nodes, equal gaps, from y_top down | posicoes dos nos."""
        tot = sum(counts[k] for k in keys) * U
        gap = (y_top - y_bot - tot) / max(len(keys) - 1, 1)
        pos, y = {}, y_top
        for k in keys:
            h = counts[k] * U
            pos[k] = (y - h, y)
            y -= h + gap
        return pos

    count_sys = {k: sum(1 for v in sys_of.values() if v == k) for k in sys_keys}
    count_blk = {k: sum(1 for v in blk_of.values() if v == k) for k in blk_keys}
    pos_src = stack(sources, count_src, y_hi, y_lo)
    pos_sys = stack(sys_keys, count_sys, y_hi, y_lo)
    pos_blk = stack(blk_keys, count_blk, y_hi, y_lo)
    xs0, xs1 = 19.0, 20.8                     # source bars
    xm0, xm1 = 38.0, 61.0                     # system pills
    xb0, xb1 = 76.0, 77.8                     # block bars
    colour = {"AD": vs.AD, "PD": vs.PD, "both": "#9a9892"}
    alpha = {"AD": 0.55, "PD": 0.55, "both": 0.45}

    def lay(flows, order_a, order_b, pos_a, pos_b):
        """EN/PT: stack the ribbons inside every node, ordered by the partner node, then by axis | empilha as fitas em cada no."""
        cur_a = {k: pos_a[k][1] for k in pos_a}      # next free top edge, leaving a node
        cur_b = {k: pos_b[k][1] for k in pos_b}
        y_a, y_b = {}, {}
        for (a, b, axk), n in sorted(flows.items(), key=lambda kv: (order_a[kv[0][0]], order_b[kv[0][1]], AX_ORDER[kv[0][2]])):
            y_a[(a, b, axk)] = (cur_a[a] - n * U, cur_a[a])
            cur_a[a] -= n * U
        for (a, b, axk), n in sorted(flows.items(), key=lambda kv: (order_b[kv[0][1]], order_a[kv[0][0]], AX_ORDER[kv[0][2]])):
            y_b[(a, b, axk)] = (cur_b[b] - n * U, cur_b[b])
            cur_b[b] -= n * U
        return [(k, y_a[k], y_b[k]) for k in flows]

    src_order = {s: i for i, s in enumerate(sources)}
    GAP_FILL = "#e6e4df"
    plt.rcParams["hatch.linewidth"] = 0.4
    plt.rcParams["hatch.color"] = "#8a8883"
    def by_reach(items):
        """EN/PT: shortest bands first, so that a long band crossing many is drawn over them and can be followed | bandas curtas primeiro."""
        return sorted(items, key=lambda it: abs(sum(it[1]) / 2 - sum(it[2]) / 2))

    for (a, b, axk), ya, yb in by_reach(lay(f1, src_order, sys_idx, pos_src, pos_sys)):
        gap_row = (a == ("", ""))
        p_ = vm.ribbon(ax, xs1, ya[0], ya[1], xm0, yb[0], yb[1], GAP_FILL if gap_row else colour[axk], alpha=0.9 if gap_row else alpha[axk], z=1, lw=0.35)
        if gap_row:
            p_.set_hatch("////")
            p_.set_edgecolor("#8a8883")
    for (a, b, axk), ya, yb in by_reach(lay(f2, sys_idx, blk_idx, pos_sys, pos_blk)):
        gap_row = (a == "none")
        p_ = vm.ribbon(ax, xm1, ya[0], ya[1], xb0, yb[0], yb[1], GAP_FILL if gap_row else colour[axk], alpha=0.9 if gap_row else alpha[axk], z=1, lw=0.35)
        if gap_row:
            p_.set_hatch("////")
            p_.set_edgecolor("#8a8883")

    # ---- nodes
    for s in sources:
        y0, y1 = pos_src[s]
        ax.add_patch(Rectangle((xs0, y0), xs1 - xs0, y1 - y0, fc=vs.INK2, ec="none", zorder=3))
        name = f"{s[0]} {s[1]} · {count_src[s]}" if s[0] else tr(lang, "no source (declared gap) · 1", "sem fonte (lacuna declarada) · 1")
        ax.text(xs0 - 0.8, (y0 + y1) / 2, name, ha="right", va="center", fontsize=6.8, color=vs.INK)
    for k, en, pt in pc.SYSTEMS:
        y0, y1 = pos_sys[k]
        i = sys_idx[k]
        face = vm.SYSTEM_RAMP[i] if k != "none" else GAP_FILL
        ax.add_patch(vs.rounded(ax, xm0, y0, xm1 - xm0, y1 - y0, fc=face, ec="white", r=0.5, lw=0.6, z=3))
        n = count_sys[k]
        if k == "none":
            ax.text((xm0 + xm1) / 2, y1 + 0.35, f"{tr(lang, en, pt)} · {n}", ha="center", va="bottom", fontsize=6.8, color=vs.INK2)
        else:
            ax.text((xm0 + xm1) / 2, (y0 + y1) / 2 + 0.55, tr(lang, en, pt), ha="center", va="center", fontsize=6.8, fontweight="bold",
                    color=vm.SYSTEM_TEXT[i], zorder=4)
            ax.text((xm0 + xm1) / 2, (y0 + y1) / 2 - 0.85, tr(lang, f"{n} rows", f"{n} linhas") if n != 1 else tr(lang, "1 row", "1 linha"),
                    ha="center", va="center", fontsize=6.6, color=vm.SYSTEM_TEXT[i], zorder=4)
    for k, en, pt in pc.BLOCKS:
        y0, y1 = pos_blk[k]
        ax.add_patch(Rectangle((xb0, y0), xb1 - xb0, y1 - y0, fc=vs.INK2, ec="none", zorder=3))
        ax.text(xb1 + 0.9, (y0 + y1) / 2, textwrap.fill(tr(lang, en, pt), 26) + f" · {count_blk[k]}", ha="left", va="center", fontsize=6.8,
                color=vs.INK, linespacing=1.15)
    # ---- headings and key
    ax.text(xs0 - 0.8, top - 1.0, tr(lang, "Primary source", "Fonte primária"), ha="right", va="top", fontsize=7.4, fontweight="bold", color=vs.INK)
    ax.text((xm0 + xm1) / 2, top - 1.0, tr(lang, "System where it was measured", "Sistema em que foi medido"), ha="center", va="top", fontsize=7.4,
            fontweight="bold", color=vs.INK)
    ax.text((xm0 + xm1) / 2, top - 3.3, tr(lang, "darker: closer to the living human brain", "mais escuro: mais perto do cérebro humano vivo"),
            ha="center", va="top", fontsize=6.6, color=vs.INK2)
    ax.text(99.0, top - 1.0, tr(lang, "Block of the model it feeds", "Bloco do modelo que alimenta"), ha="right", va="top", fontsize=7.4,
            fontweight="bold", color=vs.INK)
    ax.text(99.0, top - 3.3, tr(lang, f"band width: number of the {n_rows} table rows", f"largura da faixa: número das {n_rows} linhas da tabela"),
            ha="right", va="top", fontsize=6.6, color=vs.INK2)
    yk = 1.1
    x = 3.0
    for lab_en, lab_pt, c, a in (("Alzheimer's disease", "Doença de Alzheimer", vs.AD, 0.7), ("Parkinson's disease", "Doença de Parkinson", vs.PD, 0.7),
                                 ("both diseases", "as duas doenças", "#9a9892", 0.6), ("declared gap", "lacuna declarada", GAP_FILL, 1.0)):
        ax.add_patch(Rectangle((x, yk - 0.7), 3.2, 1.4, fc=vm.blend(c, a), ec="#8a8883" if c == GAP_FILL else "none", lw=0.4,
                               hatch="////" if c == GAP_FILL else None))
        t = ax.text(x + 4.2, yk, tr(lang, lab_en, lab_pt), ha="left", va="center", fontsize=6.8, color=vs.INK2)
        x += 4.2 + 1.15 * len(tr(lang, lab_en, lab_pt)) * 0.74 + 3.0
    fig.savefig(fig_path(FIG, "mech_provenance", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 1
def arrow(ax, p0, p1, colour, style="-|>", lw=1.1, ls="-", rad=0.0, mutation=9):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=mutation, lw=lw, color=colour, ls=ls,
                        connectionstyle=f"arc3,rad={rad}", zorder=4, shrinkA=0, shrinkB=0)
    ax.add_patch(a)


def tbar(ax, p0, p1, colour, lw=1.1, half=1.5):
    """EN/PT: dashed line ending in a bar, the sign for repression | linha tracejada que termina em barra."""
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=colour, lw=lw, ls=(0, (3, 2)), zorder=4, solid_capstyle="butt")
    ax.plot([p1[0], p1[0]], [p1[1] - half, p1[1] + half], color=colour, lw=lw + 0.4, zorder=4, solid_capstyle="butt")


def fig_schematic(d, lang):
    FW, FH = 174.0, 143.0
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    top = 100.0 * FH / FW
    ax.set_xlim(0, 100)
    ax.set_ylim(0, top)
    ax.axis("off")
    dd = d["dos"]["decay_constants"]
    hl = lambda k: LN2 / dd[k]["decay_per_hour"]
    kin = d["kin"]
    ode = d["ode"]
    prov = d["prov"]

    # ------------------------------------------------------------------ the four deposited structures
    ax.text(0.8, top - 0.6, tr(lang, "DEPOSITED STRUCTURES", "ESTRUTURAS DEPOSITADAS"), fontsize=8.2, fontweight="bold",
            color=vs.INK2, va="top")
    ax.text(28.0, top - 0.9, tr(lang, "renders of the deposited coordinates, as illustration: no result is read from them",
                                "renders das coordenadas depositadas, como ilustração: nenhum resultado é lido deles"),
            fontsize=6.6, color=vs.NOTE, va="top")
    cells = [
        ("6N4O", "a2_ago2_guide_target.png", tr(lang, "Argonaute 2 with a miRNA guide", "Argonauta 2 com um miRNA guia"),
         tr(lang, "guide, seed (red) and target", "guia, semente (vermelho) e alvo"), ("AD", "PD")),
        ("4D8C", "a2_bace1_inhibitor.png", tr(lang, "BACE1 with an inhibitor", "BACE1 com um inibidor"),
         tr(lang, "cut open; catalytic Asp pair in red", "aberta; Asp catalíticos em vermelho"), ("AD",)),
        ("5OQV", "a2_abeta42_fibril_axis.png", tr(lang, "Aβ42 fibril", "Fibrila de Aβ42"),
         tr(lang, "seen down the fibril axis", "vista pelo eixo da fibrila"), ("AD",)),
        ("6CU7", "a2_asyn_fibril_axis.png", tr(lang, "α-synuclein fibril", "Fibrila de α-sinucleína"),
         tr(lang, "seen down the fibril axis", "vista pelo eixo da fibrila"), ("PD",)),
    ]
    img_h = 17.5                                   # axis units: the image box is 22 x 17.5 units
    y_img_top = top - 4.0
    for (pdb, fn, title, desc, tags), xc in zip(cells, (12.5, 37.5, 62.5, 87.5)):
        vm.place_image_fit(fig, f"{STRUCT_FIG}/{fn}", xc / 100.0, (y_img_top - img_h / 2) / top, 22.0 * FW / 100.0, img_h * FW / 100.0, FW, FH)
        ax.text(xc, y_img_top - img_h - 0.9, title, ha="center", va="top", fontsize=6.9, fontweight="bold", color=vs.INK)
        ax.text(xc, y_img_top - img_h - 3.0, structure_line(prov[pdb], lang), ha="center", va="top", fontsize=6.6, color=vs.INK2)
        ax.text(xc, y_img_top - img_h - 4.9, desc, ha="center", va="top", fontsize=6.6, color=vs.INK2)
        ys = y_img_top - img_h - 7.9
        if len(tags) == 2:
            vm.chip(ax, xc - 2.8, ys, "AD", vs.AD, vs.AD_SOFT)
            vm.chip(ax, xc + 2.8, ys, "PD", vs.PD, vs.PD_SOFT)
        else:
            vm.chip(ax, xc, ys, tags[0], vs.AD if tags[0] == "AD" else vs.PD, vs.AD_SOFT if tags[0] == "AD" else vs.PD_SOFT)
    ax.plot([0.8, 99.2], [y_img_top - img_h - 10.6] * 2, color=vs.RULE, lw=0.6)

    # ------------------------------------------------------------------ the two cascades
    W, H, GAP, X0 = 14.0, 18.0, 2.4, 17.0

    def xs(i):
        return X0 + i * (W + GAP)

    def f1(x):
        return f"{x:.1f}"

    box_top = {"A": y_img_top - img_h - 15.0}
    box_top["P"] = box_top["A"] - H - 5.6
    yA, yP = box_top["A"] - H / 2, box_top["P"] - H / 2
    n_ord = d["tgt"]["nucleation_scaling"]["order_primary_nucleation"]
    m_star = ode["human_abeta_aggregate_load"]["M_star_nM [K014]"]
    apq = ode["enzyme_substrate_stoichiometry"]["APP_per_BACE1"]
    t_b, t_c = f1(hl("miR-29b")), f1(hl("miR-29c"))
    t_b = t_b[:-2] if t_b.endswith((".0", ",0")) else t_b
    drawn = []

    def cbox(i, yc, title, sub, colour, soft, glyph, dashed=False):
        """EN/PT: a step of the cascade: glyph, title and the measured numbers | um passo da cascata: desenho, titulo e numeros medidos."""
        x0, y0 = xs(i), yc - H / 2
        p = vs.rounded(ax, x0, y0, W, H, fc=soft, ec=colour, r=1.2, lw=1.1 if not dashed else 0.9, z=2)
        if dashed:
            p.set_linestyle((0, (3, 2)))
        gy = y0 + H - 3.1
        glyph(x0 + W / 2, gy, colour, soft)
        t1 = ax.text(x0 + W / 2, y0 + H - 6.3, title, ha="center", va="center", fontsize=7.2, color=vs.INK, fontweight="bold", zorder=3)
        t2 = ax.text(x0 + W / 2, y0 + H - 8.0, sub, ha="center", va="top", fontsize=6.3, color=vs.INK2, zorder=3, linespacing=1.2)
        drawn.append((p, t1, t2))

    # glyphs: all sized in axis units, centred on (x, y)
    g_mirna = lambda x, y, c, bg: vm.mirna_glyph(ax, x - 4.0, y, 8.0, c)
    g_mrna = lambda x, y, c, bg: vm.mrna_glyph(ax, x - 4.6, y, 9.2, c)
    g_prot = lambda x, y, c, bg: vm.protein_glyph(ax, x, y, 2.5, c, over=bg)
    g_coil = lambda x, y, c, bg: vm.coil_glyph(ax, x, y, 5.6, c)
    g_nuc = lambda x, y, c, bg: vm.nucleus_glyph(ax, x, y, 1.7, c, over=bg)
    g_layers = lambda x, y, c, bg: vm.layers_glyph(ax, x - 3.2, y + 0.6, 6.4, c, n=3, gap=0.19, over=bg)

    # ---- lane titles and a rule between the lanes
    for key, c, soft, name in (("A", vs.AD, vs.AD_SOFT, tr(lang, "ALZHEIMER", "ALZHEIMER")), ("P", vs.PD, vs.PD_SOFT, tr(lang, "PARKINSON", "PARKINSON"))):
        vs.rounded(ax, 0.6, box_top[key] - H, 3.0, H, fc=soft, ec=c, r=1.0, lw=1.0, z=2)
        ax.text(2.1, box_top[key] - H / 2, name, rotation=90, ha="center", va="center", fontsize=7.4, fontweight="bold", color=vs.INK, zorder=3)
    # ---- Alzheimer lane
    items = [
        ("miR-29a, b, c", tr(lang, f"b {t_b} h, c {t_c} h\na: no decay in 12 h\nK017–K019", f"b {t_b} h, c {t_c} h\na: sem queda em 12 h\nK017–K019"), g_mirna),
        ("BACE1 mRNA", tr(lang, f"half-life {f1(LN2 / float(kin['K068']['value_si']))} h\n5 isoforms pooled\nK067, K068", f"meia-vida {f1(LN2 / float(kin['K068']['value_si']))} h\n5 isoformas juntas\nK067, K068"), g_mrna),
        (tr(lang, "BACE1 protein", "Proteína BACE1"), tr(lang, f"half-life ≈{LN2 / float(kin['K025']['value_si']):.0f} h\nBACE1:APP = 1:{apq:.0f}\nK025, K055, K056", f"meia-vida ≈{LN2 / float(kin['K025']['value_si']):.0f} h\nBACE1:APP = 1:{apq:.0f}\nK025, K055, K056"), g_prot),
        (tr(lang, "Aβ42 monomers", "Aβ42 monômero"), tr(lang, f"production {kin['K002']['value_as_written']} %/h\nclearance {kin['K004']['value_as_written']} %/h\nin AD; K001–K004", f"produção {kin['K002']['value_as_written']} %/h\ndepuração {kin['K004']['value_as_written']} %/h\nna DA; K001–K004"), g_coil),
        (tr(lang, "Aβ42 nucleation", "Nucleação Aβ42"), tr(lang, f"orders {n_ord:.0f} and {n_ord:.0f}\nM* {m_star:.0f} nM\nrate constants: gap\nK012–K014, K037", f"ordens {n_ord:.0f} e {n_ord:.0f}\nM* {m_star:.0f} nM\nconstantes: lacuna\nK012–K014, K037"), g_nuc),
    ]
    for i, (t_, s_, g_) in enumerate(items):
        cbox(i, yA, t_, s_, vs.AD, vs.AD_SOFT, g_, dashed=(i == 4))
    # ---- Parkinson lane
    items = [
        ("miR-7", tr(lang, f"half-life {f1(hl('miR-7'))} h\nupper bound {hl('miR-7 upper bound'):.0f} h\nCyrano degrades it\nK020–K023", f"meia-vida {f1(hl('miR-7'))} h\nlimite superior {hl('miR-7 upper bound'):.0f} h\nCyrano a degrada\nK020–K023"), g_mirna),
        ("SNCA mRNA", tr(lang, f"half-life {f1(LN2 / float(kin['K066']['value_si']))} h\nrat neurons\nK066", f"meia-vida {f1(LN2 / float(kin['K066']['value_si']))} h\nneurônios de rato\nK066"), g_mrna),
        (tr(lang, "α-synuclein", "α-sinucleína"), tr(lang, f"{f1(float(kin['K054']['value_si']))} µM in a bouton\nhalf-life ≥ {LN2 / float(kin['K026']['value_si']):.0f} h\nK026, K054", f"{f1(float(kin['K054']['value_si']))} µM num botão\nmeia-vida ≥ {LN2 / float(kin['K026']['value_si']):.0f} h\nK026, K054"), g_coil),
        (tr(lang, "Elongation", "Alongamento"), tr(lang, f"saturating\nm½ {float(kin['K010']['value_si']) * 1e6:.0f}–{float(kin['K009']['value_si']) * 1e6:.0f} µM\nK009, K010", f"saturante\nm½ {float(kin['K010']['value_si']) * 1e6:.0f}–{float(kin['K009']['value_si']) * 1e6:.0f} µM\nK009, K010"), g_layers),
        (tr(lang, "Nucleation", "Nucleação"), tr(lang, "gated by pH:\nnone at neutral pH\nK011, K042", "controlada por pH:\nnenhuma em pH 7\nK011, K042"), g_nuc),
    ]
    for i, (t_, s_, g_) in enumerate(items):
        cbox(i, yP, t_, s_, vs.PD, vs.PD_SOFT, g_, dashed=(i == 4))
    # ---- arrows inside the lanes: repression (dashed bar), then production steps
    for yy, c in ((yA, vs.AD), (yP, vs.PD)):
        tbar(ax, (xs(0) + W + 0.2, yy), (xs(1) - 0.3, yy), vs.INK3)
        for i in range(1, 4):
            arrow(ax, (xs(i) + W + 0.2, yy), (xs(i + 1) - 0.3, yy), c, lw=1.2)
    # ---- what has been measured at the repression step
    ax.text(xs(1) - GAP / 2 + 1.0, box_top["A"] + 1.0, tr(lang, "repression strength not measured\nBACE1 about 50% in cells (K047)", "força da repressão não medida\nBACE1 cerca de 50% em células (K047)"),
            ha="center", va="bottom", fontsize=6.3, color=vs.INK2, linespacing=1.15)
    ax.text(xs(1) - GAP / 2 + 1.0, box_top["P"] + 1.0, tr(lang, f"repression strength not measured\nSNCA {100 * float(kin['K058']['value_si']):.0f}% and {100 * float(kin['K063']['value_si']):.0f}% in cells (K058, K063)", f"força da repressão não medida\nSNCA {100 * float(kin['K058']['value_si']):.0f}% e {100 * float(kin['K063']['value_si']):.0f}% em células (K058, K063)"),
            ha="center", va="bottom", fontsize=6.3, color=vs.INK2, linespacing=1.15)
    # ---- the mimic bolus feeds the repressor pool
    for yy, c in ((yA, vs.AD), (yP, vs.PD)):
        bx, bw, bh = 4.6, 9.0, 11.4
        vs.rounded(ax, bx, yy - bh / 2, bw, bh, fc="white", ec=c, r=1.2, lw=1.1, z=2)
        vm.pulse_glyph(ax, bx + 1.5, yy + bh / 2 - 4.2, bw - 3.0, 2.6, c, over="#ffffff")
        ax.text(bx + bw / 2, yy + 0.5, tr(lang, "Mimic", "Mimético"), ha="center", va="center", fontsize=7.0, fontweight="bold", color=vs.INK, zorder=3)
        ax.text(bx + bw / 2, yy - 1.1, tr(lang, "bolus D\nevery T h", "bolus D\na cada T h"), ha="center", va="top", fontsize=6.3, color=vs.INK2, linespacing=1.2, zorder=3)
        arrow(ax, (bx + bw + 0.2, yy), (xs(0) - 0.3, yy), c, lw=1.3)
    ax.text(99.2, 0.5, tr(lang, "dashed outline: no measured rate constants      dashed bar: repression", "contorno tracejado: sem constantes medidas      barra tracejada: repressão"),
            ha="right", va="bottom", fontsize=6.3, color=vs.NOTE)
    # ---- guard: every text must lie inside its box
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    bad = []
    for p, t1, t2 in drawn:
        bb = p.get_window_extent(rend)
        for t in (t1, t2):
            tb = t.get_window_extent(rend)
            if tb.x0 < bb.x0 + 2 or tb.x1 > bb.x1 - 2 or tb.y0 < bb.y0 + 2 or tb.y1 > bb.y1 - 2:
                bad.append(f"{t.get_text()!r}: text {tb.width:.0f} x {tb.height:.0f} px, box {bb.width:.0f} x {bb.height:.0f} px")
    if bad:
        sys.exit(f"EN/PT: text overflows its box on the schematic ({lang}):\n  " + "\n  ".join(bad))
    fig.savefig(fig_path(FIG, "mech_model_schematic", lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


FIGURES = OrderedDict([
    ("mech_model_schematic", lambda d, lang, t36: fig_schematic(d, lang)),
    ("mech_provenance", lambda d, lang, t36: fig_provenance(d, lang)),
    ("mech_alzheimer", lambda d, lang, t36: fig_alzheimer(d, lang)),
    ("mech_parkinson", lambda d, lang, t36: fig_parkinson(d, lang)),
    ("mech_turnover", lambda d, lang, t36: fig_turnover(d, lang)),
    ("mech_target_response", lambda d, lang, t36: fig_target(d, lang, t36)),
    ("mech_design_space", lambda d, lang, t36: fig_design(d, lang)),
])


def main():
    os.makedirs(FIG, exist_ok=True)
    d = read_inputs()
    pc.write_table(d["kin_rows"])
    t36 = load_module("t36", "36_mimic_target_response.py")
    for stem, fn in FIGURES.items():
        if ONLY and stem not in ONLY:
            continue
        for lang in LANGS:
            fn(d, lang, t36)
    print("=" * 78)
    print("EN | Mechanistic figures | PT | Figuras mecanisticas")
    print("=" * 78)
    for stem in FIGURES:
        if ONLY and stem not in ONLY:
            continue
        for lang in LANGS:
            print(f"EN/PT -> {fig_path(FIG, stem, lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
