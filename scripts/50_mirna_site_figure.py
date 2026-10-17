#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | The pairing of the microRNAs with the 3'UTR at the level of sequence (scripts/49). For each
     disease, the 3'UTR is drawn to scale (nucleotides) with the seed-match sites that the
     scan found, and each site is opened below as a duplex: the 3'UTR read 5' to 3' above, the
     microRNA read 3' to 5' below, every Watson-Crick pair as a bar and every G-U pair as a dot,
     the seed (nucleotides 2 to 8) shaded and nucleotide 1 marked. Every letter is read from
     the archived sequences; nothing is typed in. The figure shows complementarity, not binding
     energy, and it does not say that a site is functional.

PT | O pareamento dos microRNAs com o 3'UTR em nivel de sequencia (scripts/49). Em cada doenca, o
     3'UTR e desenhado em escala (nucleotideos) com os sitios de pareamento da semente que a
     varredura achou, e cada sitio e aberto abaixo como duplex: o 3'UTR lido de 5' a 3' em
     cima, o microRNA lido de 3' a 5' embaixo, cada par Watson-Crick como barra e cada par G-U
     como ponto, a semente (nucleotideos 2 a 8) sombreada e o nucleotideo 1 marcado. Toda
     letra e lida das sequencias arquivadas; nada e digitado. A figura mostra
     complementaridade, nao energia de ligacao, e nao diz que um sitio e funcional.

    python scripts/50_mirna_site_figure.py [--preview]     # needs scripts/49
"""

import csv
import importlib.util
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

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
STEM = "mech_mirna_site_pairing"
FW, FH = 174.0, 152.0
TOP = 100.0 * FH / FW
CW = 1.50                         # width of one nucleotide cell, axis units
BAND = "#f4f3f0"
SHORT = {"hsa-miR-29a-3p": "miR-29a", "hsa-miR-29b-3p": "miR-29b", "hsa-miR-29c-3p": "miR-29c", "hsa-miR-7-5p": "miR-7", "hsa-miR-153-3p": "miR-153"}
TYPE_H = {"8mer": 3.4, "7mer-m8": 2.5, "7mer-A1": 1.7}      # lollipop heights: stronger seed match, taller
PANELS = [("BACE1", "AD", ["hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p"]),
          ("SNCA", "PD", ["hsa-miR-7-5p", "hsa-miR-153-3p"])]


def fmt(n):
    return f"{n:,}"


def draw(rows, summ, lang):
    fig = plt.figure(figsize=(FW * MM, FH * MM))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, TOP)
    ax.axis("off")
    nt_max = max(summ["utr3"][g]["length_nt"] for g, _, _ in PANELS)
    x_bar0, x_bar1 = 6.0, 98.0
    sc = (x_bar1 - x_bar0) / nt_max
    guard = []
    y = TOP - 0.8
    for gene, dz, mirs in PANELS:
        c, soft = (vs.AD, vs.AD_SOFT) if dz == "AD" else (vs.PD, vs.PD_SOFT)
        ctext = vs.DISEASE_TEXT[dz]
        info = summ["utr3"][gene]
        # ---- header
        ph = 40.5 if dz == "AD" else 35.5
        vs.rounded(ax, 0.8, y - ph, 98.4, ph, fc="white", ec=c, r=1.4, lw=1.0, z=1)
        ax.text(2.6, y - 3.0, dz, fontsize=7.4, fontweight="bold", color=ctext, va="center", ha="left")
        ax.text(7.6, y - 3.0, tr(lang, f"{gene} 3′UTR · RefSeq {info['transcript']} · {fmt(info['length_nt'])} nt",
                                  f"3′UTR de {gene} · RefSeq {info['transcript']} · {fmt(info['length_nt'])} nt"), fontsize=7.6, fontweight="bold", color=ctext, va="center", ha="left")
        # ---- the 3'UTR to scale
        yb = y - 12.4
        L = info["length_nt"] * sc
        ax.add_patch(FancyBboxPatch((x_bar0, yb - 0.9), L, 1.8, boxstyle="round,pad=0,rounding_size=0.6", fc=soft, ec=c, lw=0.9, zorder=3))
        for t in range(0, info["length_nt"] + 1, 500):
            ax.plot([x_bar0 + t * sc] * 2, [yb - 0.9, yb - 1.7], color=vs.INK3, lw=0.6, zorder=3)
            if t % 1000 == 0:
                ax.text(x_bar0 + t * sc, yb - 2.4, fmt(t), ha="center", va="top", fontsize=6.5, color=vs.NOTE)
        ax.text(x_bar0 - 0.5, yb, tr(lang, "stop", "parada"), ha="right", va="center", fontsize=6.5, color=vs.NOTE)
        ax.text(x_bar0 + L + 0.5, yb, "poly(A)", ha="left", va="center", fontsize=6.5, color=vs.NOTE) if x_bar0 + L < 91 else ax.text(x_bar0 + L, yb + 1.6, "poly(A)", ha="right", va="bottom", fontsize=6.5, color=vs.NOTE)
        sites = [r for r in rows if r["gene"] == gene and r["mirna"] in mirs]
        # one lollipop per distinct site (the three miR-29 paralogues share their sites)
        seen = {}
        for r in sites:
            k = (int(r["utr_start"]), r["site_type"])
            seen.setdefault(k, []).append(r["mirna"])
        order = sorted(seen)
        for (st, kind), ms in seen.items():
            xs = x_bar0 + st * sc
            ax.plot([xs, xs], [yb + 0.9, yb + 0.9 + TYPE_H[kind]], color=c, lw=1.1, zorder=4)
            ax.add_patch(Circle((xs, yb + 0.9 + TYPE_H[kind]), 0.62, fc=c if kind != "7mer-A1" else "white", ec=c, lw=0.9, zorder=5))
            ms_lab = "/".join(SHORT[m].replace("miR-", "") for m in ms) if gene == "BACE1" else SHORT[ms[0]]
            lab = f"{('miR-' + ms_lab) if gene == 'BACE1' else ms_lab}"
            ax.text(xs, yb + 0.9 + TYPE_H[kind] + 1.3, f"{lab} · {fmt(st)}", ha="center", va="bottom", fontsize=6.5, color=ctext, fontweight="bold")
        # ---- the duplexes
        n_rows = len(mirs)
        blocks = order
        bw = 22 * 0 + (23 * CW + 12.0)
        xb0 = [2.8, 2.8 + bw + 2.0]
        ytop = yb - 5.6
        for bi, (st, kind) in enumerate(blocks):
            x0 = xb0[bi] if bi < len(xb0) else xb0[-1]
            xl = x0 + 9.6                       # first nucleotide cell
            rr = [r for r in sites if (int(r["utr_start"]), r["site_type"]) == (st, kind)]
            rr.sort(key=lambda r: mirs.index(r["mirna"]))
            Lm = max(int(r["mirna_length_nt"]) for r in rr)
            ncell = Lm
            # a link from the site on the bar to its block
            xs = x_bar0 + st * sc
            ax.text(x0 + 0.4, ytop - 0.6, tr(lang, f"{kind} · 3′UTR nt {fmt(st)}–{fmt(int(rr[0]['utr_end']))}", f"{kind} · 3′UTR nt {fmt(st)}–{fmt(int(rr[0]['utr_end']))}"),
                    fontsize=6.8, fontweight="bold", color=ctext, ha="left", va="top")
            yt = ytop - 6.2                      # target row
            yp = yt - 1.9                        # pairing marks
            # seed shading: nucleotides 2 to 8 are the cells (ncell-2) .. (ncell-8) from the left, nt 1 is the last cell
            def cx(k):                           # x centre of the cell facing nucleotide k (k = 1 is the rightmost cell)
                return xl + (ncell - k) * CW
            yb0 = yt + 1.2
            yb1 = yt - 2.4 - 2.4 * len(rr) - 0.2
            ax.add_patch(Rectangle((cx(8) - CW / 2, yb1), 7 * CW, yb0 - yb1, fc=vm.blend(c, 0.13), ec="none", zorder=2))
            ax.add_patch(Rectangle((cx(1) - CW / 2, yb1), CW, yb0 - yb1, fc="#ece9e3", ec="none", zorder=2))
            ax.text(cx(1) - CW * 0.95, yb0 + 0.25, tr(lang, "seed (nt 2–8)", "semente (nt 2–8)"), ha="right", va="bottom", fontsize=6.5, color=ctext)
            ax.text(cx(1), yb0 + 0.25, "A1" if kind in ("8mer", "7mer-A1") else "1", ha="center", va="bottom", fontsize=6.5, color=vs.NOTE)
            # target row
            tgt = max(rr, key=lambda r: int(r["mirna_length_nt"]))["target_5to3"]
            ax.text(xl - 1.2 * CW, yt, "5′", ha="right", va="center", fontsize=6.5, color=vs.NOTE)
            ax.text(xl + ncell * CW - 0.3 * CW, yt, "3′", ha="left", va="center", fontsize=6.5, color=vs.NOTE)
            ax.text(x0 + 0.4, yt, tr(lang, "3′UTR", "3′UTR"), ha="left", va="center", fontsize=6.5, color=vs.INK2, fontweight="bold")
            for j, ch in enumerate(tgt):
                ax.text(xl + j * CW, yt, ch, ha="center", va="center", fontsize=6.6, color=vs.INK, zorder=5)
            # microRNA rows
            for ri, r in enumerate(rr):
                yr = yt - 3.7 - 2.4 * ri
                m = summ["mirna"][r["mirna"]]["sequence"]
                marks = r["pairing_marks"]
                Lr = len(m)
                off = ncell - Lr                    # a shorter microRNA leaves its 3' end cells empty
                ax.text(x0 + 0.4, yr, SHORT[r["mirna"]], ha="left", va="center", fontsize=6.5, color=ctext, fontweight="bold")
                if ri == 0:
                    ax.text(xl - 1.2 * CW, yr, "3′", ha="right", va="center", fontsize=6.5, color=vs.NOTE)
                for j, ch in enumerate(m[::-1]):
                    col_idx = off + j
                    mk = marks[j] if j < len(marks) else " "
                    ax.text(xl + col_idx * CW, yr, ch, ha="center", va="center", fontsize=6.6, color=ctext if mk != " " else vs.NOTE,
                            fontweight="bold" if mk != " " else "normal", zorder=5)
                    if mk == "|":
                        ax.plot([xl + col_idx * CW] * 2, [yt - 0.95, yr + 0.95], color=c, lw=0.9, zorder=4) if ri == 0 else ax.plot([xl + col_idx * CW] * 2, [yr + 1.45, yr + 0.95], color=c, lw=0.9, zorder=4)
                    elif mk == ":":
                        ax.add_patch(Circle((xl + col_idx * CW, (yt + yr) / 2 if ri == 0 else yr + 1.2), 0.22, fc=vs.INK2, ec="none", zorder=4))
                ax.text(xl + ncell * CW - 0.3 * CW, yr, "5′", ha="left", va="center", fontsize=6.5, color=vs.NOTE)
            nt_lab = [(1, "1"), (8, "8"), (13, "13"), (16, "16")]
            for k, lab in nt_lab:
                ax.text(cx(k), yb1 - 0.3, lab, ha="center", va="top", fontsize=6.5, color=vs.NOTE)
            ax.text(x0 + 0.4, yb1 - 0.3, tr(lang, "miRNA nt", "nt do miRNA"), ha="left", va="top", fontsize=6.5, color=vs.NOTE)
            pairs = " · ".join(f"{SHORT[r['mirna']]} {r['pairs_total']}/{r['mirna_length_nt']}" for r in rr)
            guard.append(ax.text(x0 + 0.4, yb1 - 2.9, tr(lang, f"pairs: {pairs}", f"pares: {pairs}"), ha="left", va="top", fontsize=6.5, color=vs.INK2))
        y = y - ph - 1.6
    # ---- the key
    ky = 6.3
    items = [("8mer", tr(lang, "8mer: seed match and A1", "8mer: semente e A1")), ("7mer-m8", tr(lang, "7mer-m8: seed match", "7mer-m8: semente")), ("7mer-A1", tr(lang, "7mer-A1: nt 2–7 and A1", "7mer-A1: nt 2–7 e A1"))]
    xk = 2.4
    for kind, lab in items:
        ax.plot([xk, xk], [ky - 0.9, ky - 0.9 + TYPE_H[kind] * 0.55], color=vs.INK2, lw=1.1)
        ax.add_patch(Circle((xk, ky - 0.9 + TYPE_H[kind] * 0.55), 0.55, fc=vs.INK2 if kind != "7mer-A1" else "white", ec=vs.INK2, lw=0.9, zorder=5))
        t = ax.text(xk + 1.4, ky, lab, ha="left", va="center", fontsize=6.5, color=vs.INK2)
        fig.canvas.draw()
        xk = t.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted()).x1 + 3.4
    ax.plot([xk, xk], [ky - 0.9, ky + 0.9], color=vs.INK2, lw=0.9)
    t = ax.text(xk + 1.0, ky, tr(lang, "Watson–Crick pair", "par Watson–Crick"), ha="left", va="center", fontsize=6.5, color=vs.INK2)
    fig.canvas.draw()
    xk = t.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted()).x1 + 2.8
    ax.add_patch(Circle((xk, ky), 0.25, fc=vs.INK2, ec="none"))
    ax.text(xk + 0.9, ky, tr(lang, "G–U pair", "par G–U"), ha="left", va="center", fontsize=6.5, color=vs.INK2)
    ax.text(2.4, ky - 3.4, tr(lang, "Complementarity only: no binding energy, and a site is not shown to be functional.\nPositions count from the first nucleotide after the stop codon (hg38).",
                              "Só complementaridade: sem energia de ligação, e um sítio não é mostrado como funcional.\nPosições contadas do primeiro nucleotídeo depois do códon de parada (hg38)."),
            ha="left", va="center", fontsize=6.5, color=vs.NOTE, linespacing=1.25)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W_px, H_px = fig.get_size_inches() * fig.dpi
    bad = [t.get_text()[:40] for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()
           and (t.get_window_extent(rend).x1 > W_px - 2 or t.get_window_extent(rend).x0 < 2)]
    if bad:
        sys.exit(f"EN/PT: text outside the page in {lang}: {bad}")
    fig.savefig(fig_path(FIG, STEM, lang), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main():
    p = f"{TAB}/mirna_site_pairing.csv"
    if not os.path.exists(p):
        sys.exit("EN/PT: missing results/tables/mirna_site_pairing.csv; run scripts/49 first")
    rows = list(csv.DictReader(open(p, encoding="utf-8")))
    summ = json.load(open(f"{TAB}/mirna_site_pairing.json", encoding="utf-8"))
    os.makedirs(FIG, exist_ok=True)
    for lang in LANGS:
        draw(rows, summ, lang)
    print("=" * 78)
    print("EN | Figure of the miRNA-3'UTR pairing | PT | Figura do pareamento miRNA-3'UTR")
    print("=" * 78)
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG, STEM, lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
