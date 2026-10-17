#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | A scheme in two parts. (a) The molecular mechanism of the two diseases as the models
     take it: for each disease, where a microRNA represses the mRNA of its target, how the
     target protein is made, and what the protein does downstream (Abeta42 pool and its
     clearance in Alzheimer's disease; the saturating elongation of alpha-synuclein fibrils
     in Parkinson's disease), with the renders of the deposited structures at the steps they
     belong to. (b) How the mechanistic algorithm runs, in six steps: the parameter table is
     read under a rule that stops the run on a gap; the closed form fixes the knockdown that
     is needed and the mimic factor at a repression strength; the repeated dose is followed
     through the ODE; the dose that brings the cycle-averaged protein down by the target is
     solved at the periodic steady state; the solution is repeated over a grid of repression
     strength and dosing interval and the series are ranked; and the checks that every step
     has to pass. Each step carries a tag for what it rests on (measured, scanned, assumed,
     computed). Every number is read from the parameter table and the result tables, and the
     small plots are drawn from the same tables as Figs. 3 to 8; nothing is typed in.

PT | Um esquema em duas partes. (a) O mecanismo molecular das duas doencas como os modelos o
     tomam: em cada doenca, onde um microRNA reprime o mRNA do alvo, como a proteina-alvo e
     feita e o que ela faz adiante (reservatorio de Abeta42 e sua depuracao na doenca de
     Alzheimer; alongamento saturante das fibrilas de alfa-sinucleina na doenca de Parkinson),
     com os renders das estruturas depositadas nos passos a que pertencem. (b) Como o
     algoritmo mecanicista roda, em seis passos: a tabela de parametros e lida sob uma regra
     que para a execucao diante de uma lacuna; a forma fechada fixa a queda necessaria e o
     fator de mimetico numa forca de repressao; a dose repetida e seguida pelo sistema de EDO;
     a dose que leva a proteina media do ciclo ao alvo e resolvida no estado periodico; a
     solucao e repetida numa grade de forca de repressao e intervalo e as series sao
     ordenadas; e as verificacoes que cada passo precisa cumprir. Cada passo leva uma marca do
     que o sustenta (medido, varrido, presumido, calculado). Todo numero e lido da tabela de
     parametros e das tabelas de resultado, e os graficos pequenos sao desenhados das mesmas
     tabelas das Figs. 3 a 8; nada e digitado.

    python scripts/48_mechanism_algorithm_scheme.py [--preview]     # needs scripts/12, 13, 17, 36, 39, 40 and 46
"""

import importlib.util
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m38 = load("m38", "38_mechanistic_figures.py")
import _viz_style as vs
import _viz_mech as vm
from _bilingual import LANGS, t as tr, fig_path

FIG, TAB, MM, DPI = m38.FIG, m38.TAB, m38.MM, m38.DPI
STRUCT = f"{FIG}/structures"
STEM = "mech_mechanism_algorithm"
FW, FH = 174.0, 226.0
TOP = 100.0 * FH / FW
LN2 = math.log(2.0)
BAND = "#f4f3f0"

# ----------------------------------------------------------------------------- layout (axis units: 1 unit = 1.74 mm)
A_TOP = TOP - 5.2                         # top of the two disease columns
COL_W, COL_GAP = 48.8, 1.6
COL_X = {"AD": 0.8, "PD": 0.8 + COL_W + COL_GAP}
BAND_H, BAND_GAP, HEAD_H = 13.4, 1.7, 9.6
A_BOT = A_TOP - HEAD_H - 3 * BAND_H - 2 * BAND_GAP - 0.8
B_HEAD = A_BOT - 3.2
CARD_W, CARD_H, CARD_GAP = 31.6, 25.6, 2.0
CARD_X = [0.8, 0.8 + CARD_W + CARD_GAP, 0.8 + 2 * (CARD_W + CARD_GAP)]
ROW1_TOP = B_HEAD - 4.4
ROW_GAP = 5.4
ROW2_TOP = ROW1_TOP - CARD_H - ROW_GAP
LEG_Y = ROW2_TOP - CARD_H - 3.0


def f1(x):
    s = f"{x:.1f}"
    return s[:-2] if s.endswith(".0") else s


_SUPS = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}


def sup(s):
    """EN/PT: Unicode superscripts to mathtext, because the figure font has no superscript minus | sobrescritos Unicode para mathtext."""
    import re
    return re.sub(r"([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)", lambda m: "$^{" + "".join(_SUPS[c] for c in m.group(1)) + "}$", s)


def draw(d, lang):
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(LEG_Y - 2.4, TOP)
    ax.axis("off")
    y_lo = LEG_Y - 2.4
    span = TOP - y_lo
    kin = d["kin"]
    ode, tgt, ds = d["ode"], d["tgt"], d["ds"]
    kap_ad = ode["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
    kap_a, kap_b = m38.kv(d, "K058"), m38.kv(d, "K063")
    dd = d["dos"]["decay_constants"]
    hl = lambda k: LN2 / dd[k]["decay_per_hour"]
    boxes = []      # (patch, texts) guard: a text must lie inside its box

    def mini(x0, y0, w, h):
        a = fig.add_axes([x0 / 100.0, (y0 - y_lo) / span, w / 100.0, h / span])
        a.patch.set_alpha(0.0)
        return a

    def mini_style(a, grid=True):
        m38.style_axes(a, grid=grid)
        a.tick_params(labelsize=6.5, length=2.0, pad=1.2)
        for s_ in a.spines.values():
            s_.set_linewidth(0.6)

    def txt(x, y, s, **kw):
        kw.setdefault("fontsize", 6.6)
        kw.setdefault("color", vs.INK2)
        kw.setdefault("zorder", 6)
        return ax.text(x, y, s, **kw)

    def arrow(p0, p1, colour, lw=1.2, ls="-"):
        m38.arrow(ax, p0, p1, colour, lw=lw, ls=ls)

    def img(stem, cx, cy, bw_mm, bh_mm):
        vm.place_image_fit(fig, f"{STRUCT}/{stem}.png", cx / 100.0, (cy - y_lo) / span, bw_mm, bh_mm, FW, FH)

    # ======================================================================== (a) the molecular mechanism
    ax.text(0.4, TOP - 0.3, "a", fontsize=10, fontweight="bold", color=vs.INK, va="top")
    ax.text(3.4, TOP - 0.9, tr(lang, "The molecular mechanism, as the models take it", "O mecanismo molecular, como os modelos o tomam"),
            fontsize=7.6, fontweight="bold", color=vs.INK, va="top")
    hdr = {
        "AD": (tr(lang, "Alzheimer's disease", "Doença de Alzheimer"),
               tr(lang, "The Aβ42 pool is set by clearance, so a miR-29 mimic\nacts on the production side, through BACE1",
                  "O reservatório de Aβ42 é fixado pela depuração; um mimético\nde miR-29 age do lado da produção, via BACE1")),
        "PD": (tr(lang, "Parkinson's disease", "Doença de Parkinson"),
               tr(lang, "α-synuclein sits on the rising part of its elongation\ncurve, so a knockdown reaches the fibril less than 1:1",
                  "A α-sinucleína está na parte ascendente da curva de\nalongamento; a queda chega à fibrila menos que 1:1")),
    }
    for key in ("AD", "PD"):
        c, soft = (vs.AD, vs.AD_SOFT) if key == "AD" else (vs.PD, vs.PD_SOFT)
        x0 = COL_X[key]
        vs.rounded(ax, x0, A_BOT, COL_W, A_TOP - A_BOT, fc="white", ec=c, r=1.4, lw=1.1, z=1)
        w_ = vm.chip(ax, x0 + 4.0, A_TOP - 3.0, key, c, soft)
        txt(x0 + 7.6, A_TOP - 3.0, hdr[key][0], fontsize=8.0, fontweight="bold", color=vs.DISEASE_TEXT[key], va="center", ha="left")
        txt(x0 + 2.0, A_TOP - 6.9, hdr[key][1], fontsize=6.6, va="center", ha="left", linespacing=1.2)

    def band(key, i):
        x0 = COL_X[key] + 1.0
        y1 = A_TOP - HEAD_H - i * (BAND_H + BAND_GAP)
        y0 = y1 - BAND_H
        c = vs.AD if key == "AD" else vs.PD
        ax.add_patch(FancyBboxPatch((x0, y0), COL_W - 2.0, BAND_H, boxstyle="round,pad=0,rounding_size=0.9", fc=BAND, ec="none", zorder=2))
        return x0, y0, y1, c

    for key in ("AD", "PD"):
        soft = vs.AD_SOFT if key == "AD" else vs.PD_SOFT
        # -------------------------------------------------------------- scene 1: repression
        x0, y0, y1, c = band(key, 0)
        yc = (y0 + y1) / 2
        txt(x0 + 0.9, y1 - 1.4, tr(lang, "1  Repression", "1  Repressão"), fontweight="bold", va="center", ha="left")
        ax.add_patch(Circle((x0 + 3.6, yc - 0.3), 2.2, fc="white", ec=c, lw=0.9, zorder=4))
        ax.add_patch(Circle((x0 + 3.6, yc - 0.3), 0.75, fc=soft, ec="none", zorder=5))
        vm.mrna_glyph(ax, x0 + 6.6, yc + 0.9, 13.0, c)
        mr = "BACE1 mRNA" if key == "AD" else "SNCA mRNA"
        txt(x0 + 17.0, yc + 4.9, mr, ha="center", va="center", fontsize=6.6, color=vs.INK)
        txt(x0 + 20.4, yc + 0.9, "3′UTR", ha="left", va="center", fontsize=6.5, color=vs.NOTE)
        vm.mirna_glyph(ax, x0 + 13.0, yc - 3.4, 6.4, c)
        m38.tbar(ax, (x0 + 17.4, yc - 2.7), (x0 + 17.4, yc - 0.1), vs.INK3, half=1.1)
        mi = "miR-29" if key == "AD" else "miR-7"
        txt(x0 + 12.0, yc - 3.0, mi, ha="right", va="center", fontsize=6.6, color=vs.INK)
        if key == "AD":
            t_b, t_c = f1(hl("miR-29b")), f1(hl("miR-29c"))
            txt(x0 + 12.0, yc - 5.3, tr(lang, f"t½ {t_b} and {t_c} h", f"t½ {t_b} e {t_c} h"), ha="right", va="center", fontsize=6.5, color=vs.NOTE)
        else:
            txt(x0 + 12.0, yc - 5.3, tr(lang, f"t½ {f1(hl('miR-7'))} h", f"t½ {f1(hl('miR-7'))} h"), ha="right", va="center", fontsize=6.5, color=vs.NOTE)
        img("a2_ago2_guide_target", x0 + 30.4, yc + 0.8, 17.0, 16.0)
        txt(x0 + 30.4, y0 + 0.9, "AGO2 · 6N4O", ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
        bx = x0 + 39.2
        vs.rounded(ax, bx, yc - 3.7, 6.9, 7.4, fc="white", ec=c, r=0.9, lw=1.0, z=4)
        vm.pulse_glyph(ax, bx + 0.9, yc + 0.1, 5.1, 2.0, c, over="#ffffff")
        txt(bx + 3.45, yc - 1.7, tr(lang, "mimic", "mimético"), ha="center", va="center", fontsize=6.5, color=vs.INK, fontweight="bold")
        arrow((bx - 0.2, yc), (x0 + 36.6, yc), c, lw=1.1)
        # -------------------------------------------------------------- scene 2: production of the target
        x0, y0, y1, c = band(key, 1)
        yc = (y0 + y1) / 2
        txt(x0 + 0.9, y1 - 1.4, tr(lang, "2  Production", "2  Produção"), fontweight="bold", va="center", ha="left")
        vm.mrna_glyph(ax, x0 + 1.6, yc - 1.4, 9.0, c)
        ax.add_patch(Ellipse((x0 + 6.2, yc - 0.9), 2.6, 1.9, fc=vs.INK3, ec="white", lw=0.4, zorder=5))
        ax.add_patch(Ellipse((x0 + 6.2, yc + 0.3), 2.0, 1.2, fc=vs.INK3, ec="white", lw=0.4, zorder=5))
        arrow((x0 + 11.2, yc - 0.6), (x0 + 14.2, yc - 0.6), c)
        if key == "AD":
            img("a2_bace1_inhibitor", x0 + 19.6, yc + 0.8, 13.2, 16.0)
            txt(x0 + 18.4, y0 + 0.9, "BACE1 · 4D8C", ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
            for yy in (yc + 1.6, yc - 1.0):
                ax.plot([x0 + 26.0, x0 + 36.0], [yy, yy], color=vs.INK3, lw=1.6, solid_capstyle="round", zorder=4)
            ax.plot([x0 + 31.5, x0 + 31.5], [yc - 3.2, yc + 4.0], color=c, lw=1.5, zorder=5, solid_capstyle="round")
            ax.plot([x0 + 30.6, x0 + 32.4], [yc + 1.9, yc + 3.1], color="#d62839", lw=1.1, zorder=6)
            ax.plot([x0 + 30.6, x0 + 32.4], [yc + 3.1, yc + 1.9], color="#d62839", lw=1.1, zorder=6)
            arrow((x0 + 25.6, yc + 0.4), (x0 + 30.2, yc + 2.3), vs.INK2, lw=0.9)
            txt(x0 + 32.4, y0 + 0.9, tr(lang, "APP cut at the β site", "APP cortada no sítio β"), ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
            vm.coil_glyph(ax, x0 + 41.6, yc + 1.2, 5.4, c, seed=5)
            arrow((x0 + 33.8, yc + 2.2), (x0 + 38.2, yc + 1.4), c, lw=1.0)
            txt(x0 + 41.6, yc - 2.0, "Aβ42", ha="center", va="top", fontsize=6.6, color=vs.INK)
            txt(x0 + 46.2, y1 - 1.4, f"BACE1 t½ ≈{LN2 / float(kin['K025']['value_si']):.0f} h", ha="right", va="center", fontsize=6.5, color=vs.NOTE)
        else:
            vm.coil_glyph(ax, x0 + 19.4, yc - 0.2, 7.0, c, seed=7)
            txt(x0 + 19.4, yc - 4.2, tr(lang, "α-synuclein monomer", "monômero de α-sinucleína"), ha="center", va="top", fontsize=6.6, color=vs.INK)
            arrow((x0 + 24.0, yc - 0.4), (x0 + 28.0, yc - 0.4), c)
            txt(x0 + 29.0, yc + 1.8, tr(lang, f"{f1(float(kin['K054']['value_si']))} µM in a bouton", f"{f1(float(kin['K054']['value_si']))} µM num botão"), ha="left", va="center", fontsize=6.6, color=vs.INK)
            txt(x0 + 29.0, yc - 1.6, tr(lang, f"half-life ≥ {LN2 / float(kin['K026']['value_si']):.0f} h", f"meia-vida ≥ {LN2 / float(kin['K026']['value_si']):.0f} h"), ha="left", va="center", fontsize=6.6, color=vs.INK2)
        # -------------------------------------------------------------- scene 3: what the protein does downstream
        x0, y0, y1, c = band(key, 2)
        yc = (y0 + y1) / 2
        txt(x0 + 0.9, y1 - 1.4, tr(lang, "3  Downstream", "3  A jusante"), fontweight="bold", va="center", ha="left")
        if key == "AD":
            pa, pb = float(kin["K001"]["value_si"]), float(kin["K002"]["value_si"])
            ca, cb = float(kin["K003"]["value_si"]), float(kin["K004"]["value_si"])
            R = (pb / cb) / (pa / ca)
            tx, tw, th = x0 + 4.0, 7.2, 6.6
            fl = y0 + 1.8
            lvl = th * min(R / 1.6, 1.0)
            ax.add_patch(Rectangle((tx, fl), tw, lvl, fc=vs.AD_TEXT, ec="none", zorder=3))
            ax.plot([tx, tx, tx + tw, tx + tw], [fl + th + 0.8, fl, fl, fl + th + 0.8], color=vs.INK2, lw=1.1, zorder=4)
            ax.plot([tx - 0.5, tx + tw + 0.5], [fl + th / 1.6] * 2, color=vs.INK, lw=0.7, ls=(0, (2, 1.5)), zorder=5)
            vm.block_arrow(ax, x0 + 0.9, fl + th * 0.62, tx - 0.2, fl + th * 0.62, 14 * pb, c, head=1.0, z=5)
            vm.block_arrow(ax, tx + tw + 0.2, fl + 1.6, tx + tw + 3.6, fl + 1.6, 14 * cb, vs.INK2, head=1.0, z=5)
            txt(tx + tw / 2, fl + lvl / 2, f"{R:.2f}", ha="center", va="center", fontsize=7.0, fontweight="bold", color="white")
            txt(x0 + 15.6, y0 + 9.6, tr(lang, f"production {kin['K002']['value_as_written']} %/h", f"produção {kin['K002']['value_as_written']} %/h"), ha="left", va="center", fontsize=6.5, color=vs.INK)
            txt(x0 + 15.6, y0 + 6.6, tr(lang, f"clearance {kin['K004']['value_as_written']} %/h", f"depuração {kin['K004']['value_as_written']} %/h"), ha="left", va="center", fontsize=6.5, color=vs.INK)
            txt(x0 + 15.6, y0 + 3.6, tr(lang, f"controls {kin['K001']['value_as_written']} and {kin['K003']['value_as_written']}", f"controles {kin['K001']['value_as_written']} e {kin['K003']['value_as_written']}"), ha="left", va="center", fontsize=6.5, color=vs.NOTE)
            arrow((x0 + 32.4, yc), (x0 + 34.4, yc), c, lw=1.1)
            img("a2_abeta42_fibril_axis", x0 + 40.6, yc + 0.8, 18.0, 16.0)
            txt(x0 + 40.6, y0 + 0.9, "Aβ42 · 5OQV", ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
        else:
            a = mini(x0 + 5.6, y0 + 3.7, 14.6, 7.0)
            S = np.linspace(0, 120, 200)
            m1, m2 = m38.kv(d, "K010") * 1e6, m38.kv(d, "K009") * 1e6
            s0 = float(kin["K054"]["value_si"])
            for m_ in (m1, m2):
                a.plot(S, S / (m_ + S), color=c, lw=1.0)
            a.axvline(s0, color=vs.INK3, lw=0.7, ls=(0, (2, 1.5)))
            a.scatter([s0], [s0 / (0.5 * (m1 + m2) + s0)], s=9, color=c, zorder=5, edgecolor="white", linewidth=0.4)
            a.set_xlim(0, 120); a.set_ylim(0, 0.75)
            a.set_xticks([0, 60, 120]); a.set_xticklabels(["0", "60", "120 µM"]); a.set_yticks([0, 0.5])
            mini_style(a)
            a.set_ylabel(tr(lang, "rate", "taxa"), fontsize=6.5, labelpad=0.5, color=vs.INK2)
            fr_lo, fr_hi = s0 / (m2 + s0), s0 / (m1 + s0)
            txt(x0 + 22.4, yc + 2.4, tr(lang, f"{100 * fr_lo:.0f}–{100 * fr_hi:.0f}% of the\nmaximum rate", f"{100 * fr_lo:.0f}–{100 * fr_hi:.0f}% da\ntaxa máxima"), ha="left", va="center", fontsize=6.5, color=vs.INK, linespacing=1.1)
            sl = lambda kd, m_: 100 * (1 - ((1 - kd) * s0 / (m_ + (1 - kd) * s0)) / (s0 / (m_ + s0)))
            lo_a, hi_a = sorted((sl(kap_a, m1), sl(kap_a, m2)))
            txt(x0 + 22.4, yc - 2.6, tr(lang, f"−{100 * kap_a:.0f}% protein:\n−{lo_a:.1f}–{hi_a:.1f}% rate", f"−{100 * kap_a:.0f}% de proteína:\n−{lo_a:.1f}–{hi_a:.1f}% na taxa"), ha="left", va="center", fontsize=6.5, color=vs.INK, linespacing=1.1)
            arrow((x0 + 34.2, yc), (x0 + 35.6, yc), c, lw=1.1)
            img("a2_asyn_fibril_axis", x0 + 41.0, yc + 0.8, 18.0, 16.0)
            txt(x0 + 41.0, y0 + 0.9, "α-syn · 6CU7", ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
        for i in (0, 1):
            xb = COL_X[key] + COL_W / 2
            yb = A_TOP - HEAD_H - i * (BAND_H + BAND_GAP) - BAND_H
            arrow((xb, yb + 0.1), (xb, yb - BAND_GAP + 0.15), vs.INK3, lw=1.0, ls="-")

    # ======================================================================== (b) the algorithm
    ax.text(0.4, B_HEAD + 0.6, "b", fontsize=10, fontweight="bold", color=vs.INK, va="center")
    ax.text(3.4, B_HEAD + 0.6, tr(lang, "How the mechanistic algorithm runs", "Como o algoritmo mecanicista roda"), fontsize=7.6, fontweight="bold", color=vs.INK, va="center")

    STATUS = {"measured": ("-", tr(lang, "measured", "medido")), "scanned": ((0, (3, 2)), tr(lang, "scanned", "varrido")),
              "assumed": ((0, (1, 1.6)), tr(lang, "assumed", "presumido")), "computed": ("-", tr(lang, "computed", "calculado"))}

    def card(row, col, n, title, status, lines, graphic=None):
        x0 = CARD_X[col]
        y1 = ROW1_TOP if row == 0 else ROW2_TOP
        y0 = y1 - CARD_H
        p = vs.rounded(ax, x0, y0, CARD_W, CARD_H, fc="white", ec=vs.INK2, r=1.2, lw=0.9, z=2)
        ax.add_patch(Circle((x0 + 2.6, y1 - 2.6), 1.55, fc=vs.INK, ec="none", zorder=5))
        ax.text(x0 + 2.6, y1 - 2.65, str(n), ha="center", va="center", fontsize=7.0, fontweight="bold", color="white", zorder=6)
        t1 = txt(x0 + 5.0, y1 - 2.6, title, ha="left", va="center", fontsize=7.0, fontweight="bold", color=vs.INK)
        ls_, lab_ = STATUS[status]
        ax.plot([x0 + 1.0, x0 + CARD_W - 1.0], [y1 - 4.7, y1 - 4.7], color=vs.RULE, lw=0.6, zorder=3)
        if status != "computed":
            ax.plot([x0 + 1.6, x0 + 4.6], [y1 - 6.5, y1 - 6.5], color=vs.INK2, lw=1.4, ls=ls_, zorder=4)
        t2 = txt(x0 + (5.4 if status != "computed" else 1.6), y1 - 6.5, lab_, ha="left", va="center", fontsize=6.5, color=vs.INK2, fontweight="bold")
        t3 = txt(x0 + 1.4, y0 + 0.9, lines, ha="left", va="bottom", fontsize=6.5, color=vs.INK2, linespacing=1.2)
        boxes.append((p, [t1, t2, t3]))
        return x0, y0, y1

    gx = lambda x0, off, w: (x0 + off, w)

    # ---- 1. the parameter table
    x0, y0, y1 = card(0, 0, 1, tr(lang, "Read the parameter table", "Ler a tabela de parâmetros"), "measured",
                      tr(lang, "A row loads only if measured, or derived\nand requested by name.",
                         "Uma linha só carrega se medida, ou derivada\ne pedida pelo nome."))
    kinds = [("numeric", tr(lang, "measured", "medidas"), vs.INK), ("derived", tr(lang, "derived", "derivadas"), vs.INK2),
             ("qualitative_constraint", tr(lang, "qualitative", "qualitativas"), vs.INK3), ("declared_gap", tr(lang, "gap", "lacuna"), "white")]
    counts = {k: sum(1 for r in d["kin_rows"] if r["kind"] == k) for k, _, _ in kinds}
    tot = sum(counts.values())
    bx, bw, by = x0 + 1.6, CARD_W - 3.2, y1 - 13.2
    xx = bx
    for k, lab, col in kinds:
        w_ = bw * counts[k] / tot
        ax.add_patch(Rectangle((xx, by), max(w_, 0.5), 3.0, fc=col, ec=vs.INK2 if k == "declared_gap" else "white", lw=0.7 if k == "declared_gap" else 0.5,
                               zorder=4, hatch="////" if k == "declared_gap" else None))
        xx += w_
    txt(bx, by + 4.4, tr(lang, f"{tot} rows from {int(float(d['nums']['kin.sources']))} primary sources", f"{tot} linhas de {int(float(d['nums']['kin.sources']))} fontes primárias"), ha="left", va="center", fontsize=6.5, color=vs.NOTE)
    txt(bx, by - 1.5, f"{counts['numeric']} {kinds[0][1]}, {counts['derived']} {kinds[1][1]}, {counts['qualitative_constraint']} {kinds[2][1]}", ha="left", va="center", fontsize=6.5, color=vs.INK2)
    ax.add_patch(Polygon([(bx + 0.9 + 0.95 * math.cos(math.radians(22.5 + 45 * k)), by - 4.1 + 0.95 * math.sin(math.radians(22.5 + 45 * k))) for k in range(8)], closed=True, fc="#d03b3b", ec="none", zorder=5))
    txt(bx + 2.4, by - 4.1, tr(lang, f"{counts['declared_gap']} declared gap: stops the run", f"{counts['declared_gap']} lacuna declarada: para a execução"), ha="left", va="center", fontsize=6.5, color=vs.INK2)

    # ---- 2. closed form
    x0, y0, y1 = card(0, 1, 2, tr(lang, "Closed form at steady state", "Forma fechada, estado estacionário"), "scanned",
                      tr(lang, "κ = 1 − 1/R, from the table alone\nf = 1 + κ/(1−κ) · (1+ρ)/ρ,  ρ scanned",
                         "κ = 1 − 1/R, só com a tabela\nf = 1 + κ/(1−κ) · (1+ρ)/ρ,  ρ varrido"))
    a = mini(x0 + 4.4, y0 + 7.4, CARD_W - 7.2, 8.4)
    rho = np.logspace(-2, 3, 200)
    for kap, col, lab in ((kap_ad, vs.AD, "AD"), (kap_a, vs.PD, "PD")):
        a.plot(rho, m38.fold_for(kap, rho), color=col, lw=1.1)
        a.axhline(1.0 / (1.0 - kap), color=col, lw=0.6, ls=(0, (1, 1.6)))
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(0.01, 1000); a.set_ylim(1, 100)
    a.set_xticks([0.01, 1, 100]); a.set_xticklabels(["0.01", "1", "100"]); a.set_yticks([1, 10, 100]); a.set_yticklabels(["1", "10", "100"])
    a.minorticks_off()
    mini_style(a)
    a.set_xlabel("ρ", fontsize=6.5, labelpad=0.3, color=vs.INK2)
    a.set_ylabel(tr(lang, "mimic factor f", "fator f"), fontsize=6.5, labelpad=0.5, color=vs.INK2)

    # ---- 3. the repeated dose as an ODE
    x0, y0, y1 = card(0, 2, 3, tr(lang, "Follow the repeated dose", "Seguir a dose repetida"), "assumed",
                      tr(lang, "bolus D every T h on the repressor E, then\nmRNA μ, protein π (and Aβ42 a). The mimic\nis cleared like the endogenous species.",
                         "bolus D a cada T h no repressor E, depois\nmRNA μ, proteína π (e Aβ42 a). O mimético\né depurado como a espécie endógena."))
    meas, _ = d["t36"].load_module("ode12", "12_ode_models_calibrated.py").load_measured()
    T24 = 24.0
    midr = tgt["repression_strength"]["PD"]["rho_scan_mid"]
    r = m38.scen(d["rows"], "PD", "scan midpoint", "miR-7", T24, kap_a)
    dx = LN2 / float(r["half_life_h"])
    cas = d["t36"].Cascade(meas["d_mRNA_SNCA"], midr, meas["d_aSyn"], dx)
    tt, yy = m38.trajectory(cas, float(r["dose"]), T24, dx)
    gh = 3.1
    for k_, (lab, col, ylim) in enumerate((("E", vs.PD, (0.8, 60)), ("μ", vs.INK2, (-0.05, 1.15)), ("π", vs.INK2, (0.55, 1.05)))):
        a = mini(x0 + 5.6, y1 - 8.0 - gh * (k_ + 1) - 0.1 * k_, CARD_W - 8.4, gh)
        yv = 1 + yy[0] if k_ == 0 else yy[k_]
        a.plot(tt, yv, color=col, lw=1.0)
        if k_ == 0:
            a.set_yscale("log")
        a.set_ylim(*ylim); a.set_xlim(0, 48)
        a.set_xticks([]); a.set_yticks([]); a.minorticks_off()
        for s_ in ("top", "right", "left", "bottom"):
            a.spines[s_].set_visible(False)
        if k_ == 2:
            a.axhline(1.0 - kap_a, color=vs.INK3, lw=0.7, ls=(0, (2, 1.5)))
        ax.text(x0 + 4.6, y1 - 8.0 - gh * (k_ + 1) - 0.1 * k_ + gh / 2, lab, ha="right", va="center", fontsize=6.6, color=vs.INK, fontstyle="italic")

    # ---- 4. solve for the dose
    x0, y0, y1 = card(1, 0, 4, tr(lang, "Solve for the dose", "Resolver a dose"), "computed",
                      tr(lang, "One cycle is an affine map; solve for D so\nthat the cycle mean of the protein falls by κ.",
                         "Um ciclo é um mapa afim; resolver D para\nque a média do ciclo da proteína caia κ."))
    a = mini(x0 + 4.6, y0 + 8.6, CARD_W - 7.6, 7.8)
    intervals = [8.0, 12.0, 24.0, 48.0, 168.0]
    for sp, col, ls_ in (("miR-7", vs.PD, "-"), ("median miRNA", vs.INK2, "-")):
        pk = []
        for Ti in intervals:
            rr = m38.scen(d["rows"], "PD", "scan midpoint", sp, Ti, kap_a)
            pk.append(float(rr["peak_pool_over_baseline"]) if rr["reachable"] == "True" else np.nan)
        a.plot(intervals, pk, color=col, lw=1.0, marker="o", ms=2.2)
        if any(np.isnan(pk)):
            i0 = int(np.argmax(np.isnan(pk)))
            a.plot([intervals[i0]], [1.2e4], marker="o", ms=3.0, mfc="white", mec=col, mew=0.7, ls="none", clip_on=False)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_ylim(1, 3e4); a.set_xlim(7, 190)
    a.set_xticks([8, 24, 168]); a.set_xticklabels(["8", "24", "168"]); a.set_yticks([1, 100, 10000]); a.set_yticklabels(["1", "100", "10⁴"])
    a.minorticks_off()
    mini_style(a)
    a.set_xlabel(tr(lang, "interval T (h)", "intervalo T (h)"), fontsize=6.5, labelpad=0.3, color=vs.INK2)
    a.set_ylabel(tr(lang, "peak of E", "pico de E"), fontsize=6.5, labelpad=0.5, color=vs.INK2)

    # ---- 5. sweep and rank
    gridp = len(ds["rho_grid"]) * len(ds["interval_grid_h"])
    x0, y0, y1 = card(1, 1, 5, tr(lang, "Sweep ρ and T, then rank", "Varrer ρ e T, depois ordenar"), "computed",
                      tr(lang, f"{len(ds['rho_grid'])} values of ρ × {len(ds['interval_grid_h'])} intervals × 9 series;\nrank by the longest T within a 3×, 10×\nor 30× peak.",
                         f"{len(ds['rho_grid'])} valores de ρ × {len(ds['interval_grid_h'])} intervalos × 9 séries;\nordenar pelo maior T com pico de até\n3×, 10× ou 30×."))
    t_grid = np.array(ds["interval_grid_h"]); rho_grid = np.array(ds["rho_grid"])
    Z = np.full((len(t_grid), len(rho_grid)), np.nan)
    for rw in d["ds_rows"]:
        if rw["axis"] == "PD" and rw["species"] == "miR-7 upper bound" and abs(float(rw["target_knockdown"]) - kap_a) < 1e-3 and rw["reachable"] == "True":
            ti = int(np.argmin(np.abs(t_grid - float(rw["interval_h"]))))
            ri = int(np.argmin(np.abs(np.log10(rho_grid) - math.log10(float(rw["rho"])))))
            Z[ti, ri] = float(rw["peak_pool_over_baseline"])
    Rf, Tf, Zf, near = m38.fine_mesh(Z, t_grid, rho_grid)
    a = mini(x0 + 4.2, y0 + 8.6, 12.4, 7.8)
    a.pcolormesh(Rf, Tf, Zf, cmap=vm.PD_MAP, norm=LogNorm(vmin=1.0, vmax=1000.0), shading="nearest", rasterized=True)
    a.contour(Rf, Tf, Zf, levels=[10], colors=[vs.INK], linewidths=0.8)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.01, 1e4); a.set_ylim(t_grid.min(), t_grid.max())
    a.set_xticks([0.01, 1, 100]); a.set_xticklabels(["0.01", "1", "100"]); a.set_yticks([2, 24, 168]); a.set_yticklabels(["2", "24", "168"])
    a.minorticks_off()
    mini_style(a, grid=False)
    a.set_xlabel("ρ", fontsize=6.5, labelpad=0.3, color=vs.INK2)
    a.set_ylabel("T (h)", fontsize=6.5, labelpad=0.5, color=vs.INK2)
    rk = d["rank"]
    base_x = x0 + 23.2
    names = [("miR-7 upper bound", "5 h"), ("miR-7", "1.7 h")]
    ax.text(base_x + 2.4, y1 - 9.6, tr(lang, "T within 10×", "T dentro de 10×"), ha="center", va="center", fontsize=6.5, color=vs.NOTE)
    for j, (sp, hlab) in enumerate(names):
        rows_sp = [rr for rr in d["rank_rows"] if rr["axis"] == "PD" and rr["species"] == sp and rr["threshold_x"] == "10" and rr["anchor"] == "scan midpoint" and abs(float(rr["target_knockdown"]) - kap_a) < 1e-3]
        val = float(rows_sp[0]["interval_h"])
        yb = y1 - 12.4 - j * 2.9
        ax.add_patch(Rectangle((base_x, yb - 0.9), 4.0 * min(val, 60.0) / 60.0, 1.8, fc=vs.PD if j == 0 else vm.blend(vs.PD, 0.55), ec="none", zorder=4))
        txt(base_x + 4.0 * min(val, 60.0) / 60.0 + 0.5, yb, f"{val:.0f} h", ha="left", va="center", fontsize=6.5, color=vs.INK)
        txt(base_x - 0.4, yb, hlab, ha="right", va="center", fontsize=6.5, color=vs.INK2)

    # ---- 6. checks
    x0, y0, y1 = card(1, 2, 6, tr(lang, "Check every step", "Verificar cada passo"), "computed", "")
    n = d["nums"]
    checks = [
        tr(lang, f"closed form against the ODE: {sup(n['ad.ode_check'])}", f"forma fechada contra a EDO: {sup(n['ad.ode_check'])}"),
        tr(lang, f"six independent calculations of\nthe dosing code, worst {n['dose.verif_worst']} of tolerance", f"seis cálculos independentes do\ncódigo de dose, pior {n['dose.verif_worst']} da tolerância"),
        tr(lang, f"grid against the scenarios at 24 h,\n{n['ds.verif_points']} points: {sup(n['ds.verif_err'])}", f"grade contra os cenários a 24 h,\n{n['ds.verif_points']} pontos: {sup(n['ds.verif_err'])}"),
        tr(lang, "consistency checks between tables,\nmodels and figures (scripts/08)", "verificações de consistência entre tabelas,\nmodelos e figuras (scripts/08)"),
    ]
    yy_ = y1 - 9.4
    for s in checks:
        ax.add_patch(Circle((x0 + 3.0, yy_ - 0.2 - (s.count("\n") * 0.9)), 1.15, fc="#0b7a0b", ec="none", zorder=5))
        ax.plot([x0 + 2.35, x0 + 2.85, x0 + 3.75], [yy_ - 0.1 - s.count("\n") * 0.9, yy_ - 0.65 - s.count("\n") * 0.9, yy_ + 0.45 - s.count("\n") * 0.9], color="white", lw=0.9, zorder=6, solid_capstyle="round")
        tx = txt(x0 + 5.4, yy_, s, ha="left", va="center", fontsize=6.5, color=vs.INK2, linespacing=1.15)
        boxes[-1][1].append(tx)
        yy_ -= 3.0 + 1.8 * s.count("\n")

    # ---- arrows between the cards (a snake: 1 → 2 → 3, then down and back to 4 → 5 → 6)
    for row, y1c in ((0, ROW1_TOP), (1, ROW2_TOP)):
        yc = y1c - CARD_H / 2
        for i in (0, 1):
            arrow((CARD_X[i] + CARD_W + 0.1, yc), (CARD_X[i + 1] - 0.1, yc), vs.INK, lw=1.3)
    xm = CARD_X[2] + CARD_W / 2
    ym = ROW1_TOP - CARD_H - ROW_GAP / 2
    ax.plot([xm, xm], [ROW1_TOP - CARD_H, ym], color=vs.INK, lw=1.3, zorder=4, solid_capstyle="butt")
    ax.plot([xm, CARD_X[0] + CARD_W / 2], [ym, ym], color=vs.INK, lw=1.3, zorder=4, solid_capstyle="butt")
    arrow((CARD_X[0] + CARD_W / 2, ym), (CARD_X[0] + CARD_W / 2, ROW2_TOP + 0.1), vs.INK, lw=1.3)

    # ---- the tags
    ky = LEG_Y
    txt(0.8, ky + 0.2, tr(lang, "A step rests on:", "O passo se apoia em:"), ha="left", va="center", fontweight="bold", color=vs.INK2)
    items = [("measured", tr(lang, "from the source sentences", "das frases de origem"), "-"),
             ("scanned", tr(lang, "unmeasured, swept over a range", "não medido, varrido numa faixa"), (0, (3, 2))),
             ("assumed", tr(lang, "stated", "declarada"), (0, (1, 1.6)))]
    xk = 17.5
    for k, lab, ls_ in items:
        ax.plot([xk, xk + 3.0], [ky + 0.2, ky + 0.2], color=vs.INK2, lw=1.4, ls=ls_)
        t_ = txt(xk + 3.8, ky + 0.2, f"{STATUS[k][1]}: {lab}", ha="left", va="center", fontsize=6.5, color=vs.INK2)
        fig.canvas.draw()
        bb = t_.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
        xk = bb.x1 + 2.2
    txt(0.8, ky - 2.6, tr(lang, "Renders: deposited structures, as illustration; no number is read from them. Small plots: drawn from the result tables of Figs. 3 to 8.",
                          "Renders: estruturas depositadas, como ilustração; nenhum número é lido deles. Gráficos pequenos: desenhados das tabelas de resultado das Figs. 3 a 8."),
        ha="left", va="center", fontsize=6.5, color=vs.NOTE)

    # ---- guards
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    bad = []
    for p, ts in boxes:
        bb = p.get_window_extent(rend)
        for t in ts:
            tb = t.get_window_extent(rend)
            if tb.x0 < bb.x0 + 1 or tb.x1 > bb.x1 - 1 or tb.y0 < bb.y0 + 1 or tb.y1 > bb.y1 - 1:
                bad.append(f"{t.get_text()[:34]!r}: text {tb.width:.0f} x {tb.height:.0f} px, box {bb.width:.0f} x {bb.height:.0f} px")
    W_px, H_px = fig.get_size_inches() * fig.dpi
    for tx in fig.findobj(matplotlib.text.Text):
        if tx.get_text().strip() and tx.get_visible():
            tb = tx.get_window_extent(rend)
            if tb.x1 > W_px + 1 or tb.x0 < -1:
                bad.append(f"outside the page: {tx.get_text()[:34]!r}")
    if bad:
        sys.exit(f"EN/PT: text problems on the scheme ({lang}):\n  " + "\n  ".join(bad))
    fig.savefig(fig_path(FIG, STEM, lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    d = m38.read_inputs()
    raw = json.load(open(f"{TAB}/manuscript2_numbers.json", encoding="utf-8"))
    d["nums"] = raw.get("values", raw)
    d["t36"] = m38.load_module("t36", "36_mimic_target_response.py")
    d["rank"] = json.load(open(f"{TAB}/mimic_series_ranking.json", encoding="utf-8"))
    import csv
    d["rank_rows"] = list(csv.DictReader(open(f"{TAB}/mimic_series_ranking.csv", encoding="utf-8")))
    os.makedirs(FIG, exist_ok=True)
    for lang in LANGS:
        draw(d, lang)
    print("=" * 78)
    print("EN | Scheme of the mechanism and the algorithm | PT | Esquema do mecanismo e do algoritmo")
    print("=" * 78)
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG, STEM, lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
