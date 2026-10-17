#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Shared visual style for the manuscript figures: one palette, one type scale and a
     few drawing helpers, so that every figure reads as part of the same paper.

     Colour is assigned by the job it does. Disease identity is the only categorical
     encoding (AD blue, PD orange, two slots of a palette checked for colour-vision
     deficiency with the data-viz validator, all-pairs mode). Strength of evidence is a
     one-hue ordinal ramp, dark = strongest (validated as ordinal: monotone lightness,
     visible steps, light end at least 2:1 against the page). Risk of bias uses the
     fixed status colours and always carries a letter (L, U, H) as well, so colour is
     never the only signal. Marker type (single microRNA, panel) is carried by shape and
     fill, not by a third hue.

PT | Estilo visual compartilhado das figuras do manuscrito: uma paleta, uma escala
     tipografica e alguns auxiliares de desenho, para que toda figura pareca parte do
     mesmo artigo.

     A cor e atribuida pela funcao que cumpre. A identidade da doenca e a unica
     codificacao categorica (DA azul, DP laranja, dois lugares de uma paleta checada para
     daltonismo com o validador de visualizacao de dados, modo todos os pares). A forca da
     evidencia e uma rampa ordinal de um so matiz, escuro = mais forte. O risco de vies usa
     as cores de status fixas e sempre leva uma letra (L, U, H), de modo que a cor nunca e
     o unico sinal. O tipo de marcador (miRNA isolado, painel) e dado por forma e
     preenchimento, nao por um terceiro matiz.
"""

import matplotlib
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch

# --- ink and surface (light, print) ---------------------------------------------------
PAGE = "#ffffff"
INK = "#1b1b1a"
INK2 = "#52514e"
INK3 = "#8a8984"          # EN/PT: rules, arrows, markers | filetes, setas, marcadores
NOTE = "#62615d"          # EN/PT: grey lettering, 5.4 to 6.2:1 on the page and on BAND | texto cinza
GRID = "#e7e6e2"
BAND = "#f4f3f0"
RULE = "#cfcdc7"

# --- categorical: disease identity (validated, slots 1 and 2) -------------------------
AD = "#2a78d6"
PD = "#eb6834"
DISEASE = {"AD": AD, "PD": PD}
# EN | Lettering in the disease hues, darker so that it reaches 4.5:1 on the page and on BAND (the journal rule);
#      marks, lines and fills keep AD and PD above. Measured on the vector PDFs by scripts/45.
# PT | Texto nos matizes de doenca, mais escuro para chegar a 4,5:1 na pagina e na faixa; marcas, linhas e
#      preenchimentos mantem AD e PD acima.
AD_TEXT = "#1f63b8"
PD_TEXT = "#b8400f"
DISEASE_TEXT = {"AD": AD_TEXT, "PD": PD_TEXT}
AD_SOFT = "#cde2fb"
PD_SOFT = "#f9d3c1"

# --- ordinal ramp for strength of evidence, strongest first (validated, ordinal) ------
TIER_RAMP = ["#0d366b", "#124785", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef"]
TIER_NONE = "#cfcdc7"

# --- status (fixed palette): risk of bias ---------------------------------------------
ROB_FILL = {"low": "#0b7a0b", "unclear": "#fab219", "high": "#d03b3b"}
ROB_LETTER = {"low": "L", "unclear": "U", "high": "H"}
ROB_TEXT = {"low": "#ffffff", "unclear": "#1b1b1a", "high": "#ffffff"}


def apply():
    """EN/PT: set the global matplotlib style | define o estilo global do matplotlib."""
    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.edgecolor": RULE,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": PAGE,
        "axes.facecolor": PAGE,
        "savefig.facecolor": PAGE,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })


def rounded(ax, x, y, w, h, fc, ec="none", r=0.12, lw=0.0, z=2):
    """EN/PT: rounded rectangle in data coordinates | retangulo arredondado em coordenadas de dados."""
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z)
    ax.add_patch(p)
    return p


def panel_letter(ax, letter, x=-0.02, y=1.02, size=11):
    """EN/PT: bold panel letter above the top-left corner | letra de painel acima do canto superior esquerdo."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size, fontweight="bold",
            color=INK, ha="left", va="bottom")


apply()
