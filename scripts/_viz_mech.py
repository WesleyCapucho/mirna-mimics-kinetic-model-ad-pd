#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Drawing helpers for the figures of the mechanistic study (Article 2), on top of the shared
     style of scripts/_viz_style.py: lettering that stays legible over artwork, a vertical
     gradient fill, curved ribbons for the flow diagram, a tank drawn from measured rates,
     molecular renders placed as images at a stated resolution, and a one-hue ordinal ramp.

     Nothing here knows a number from the study. Every size, width and height is passed in by
     the caller from a table, so that a drawing whose geometry carries data (the tanks, the
     ribbons) cannot drift from the table it is drawn from.

PT | Auxiliares de desenho para as figuras do estudo mecanistico (Artigo 2), sobre o estilo
     compartilhado do scripts/_viz_style.py: texto que continua legivel sobre a arte,
     preenchimento em gradiente vertical, fitas curvas para o diagrama de fluxo, um tanque
     desenhado a partir de taxas medidas, renders moleculares colocados como imagens numa
     resolucao declarada e uma rampa ordinal de um so matiz.

     Nada aqui conhece um numero do estudo. Todo tamanho, largura e altura vem do chamador, a
     partir de uma tabela, para que um desenho cuja geometria carrega dados (os tanques, as
     fitas) nao possa se afastar da tabela de onde sai.
"""

import os

import numpy as np
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap, to_rgb

import _viz_style as vs

DPI_IMAGE = 600                      # EN/PT: resolution of an embedded render at its printed size

# EN | Without this, the PDF and EPS backends resample every placed render to the page resolution
#      (100 dpi were measured on the vector files of the figures, against 600 dpi in the arrays), so
#      each render is kept as its own image, at its own resolution, as the journal asks for.
# PT | Sem isto, os backends de PDF e EPS reamostram cada render colocado para a resolucao da pagina
#      (100 dpi medidos nos arquivos vetoriais das figuras, contra 600 dpi nos arrays), entao cada
#      render e mantido como imagem propria, na resolucao propria, como a revista pede.
import matplotlib as _mpl
_mpl.rcParams["image.composite_image"] = False
_mpl.rcParams["figure.dpi"] = DPI_IMAGE      # EN/PT: the vector backends resample a render to the figure dpi | os backends vetoriais reamostram ao dpi da figura

# EN/PT: one-hue ordinal ramp for "how close the system is to the human brain in vivo", deepest first
SYSTEM_RAMP = ["#0b3b3c", "#145558", "#1d7573", "#34a096", "#4fb0a8", "#74c2ba"]   # checked with the dataviz validator, --ordinal: all pass
SYSTEM_TEXT = ["#ffffff", "#ffffff", "#ffffff", "#1b1b1a", "#1b1b1a", "#1b1b1a"]

# EN/PT: sequential ramps for the design-space maps, one per disease, same lightness steps
def ramp(colour_deep, colour_pale="#f6f5f2", n=256, name="ramp"):
    return LinearSegmentedColormap.from_list(name, [colour_pale, colour_deep], N=n)


AD_MAP = LinearSegmentedColormap.from_list("ad_map", ["#f3f7fd", "#cde2fb", "#86b6ef", "#2a78d6", "#124785", "#0d366b"])
PD_MAP = LinearSegmentedColormap.from_list("pd_map", ["#fdf5f0", "#f9d3c1", "#f5ad8d", "#eb6834", "#a63d12", "#6e2509"])


def blend(colour, alpha, over="#ffffff"):
    """
    EN | The opaque colour that `colour` at opacity `alpha` makes over `over`. Figures that go to a
         journal as EPS use these instead of transparency, which PostScript cannot hold: cairo
         would turn every translucent patch into a bitmap, and the line art would stop being line art.
    PT | A cor opaca que `colour` com opacidade `alpha` forma sobre `over`. Figuras que vao ao
         periodico em EPS usam isto no lugar de transparencia, que o PostScript nao guarda: o cairo
         transformaria cada trecho translucido em bitmap, e o desenho de linhas deixaria de ser vetorial.
    """
    c, o = np.array(to_rgb(colour)), np.array(to_rgb(over))
    return tuple(float(x) for x in alpha * c + (1.0 - alpha) * o)


def halo(text, width=2.2, colour="white"):
    """EN/PT: white outline round lettering so that it reads over lines and fills | contorno branco no texto."""
    text.set_path_effects([pe.withStroke(linewidth=width, foreground=colour)])
    return text


def gradient_rect(ax, x0, x1, y0, y1, c_bottom, c_top, z=2, alpha=1.0, clip=None):
    """EN/PT: rectangle filled with a vertical gradient | retangulo com gradiente vertical."""
    n = 64
    g = np.linspace(0.0, 1.0, n).reshape(n, 1)
    cb, ct = np.array(to_rgb(c_bottom)), np.array(to_rgb(c_top))
    img = np.zeros((n, 1, 3))
    for k in range(3):
        img[:, 0, k] = cb[k] + (ct[k] - cb[k]) * g[:, 0]
    if alpha != 1.0:
        img = np.concatenate([img, np.full((n, 1, 1), alpha)], axis=2)
    im = ax.imshow(img, extent=(x0, x1, y0, y1), origin="lower", aspect="auto", zorder=z, interpolation="bilinear")
    if clip is not None:
        im.set_clip_path(clip)
    return im


def ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, colour, alpha=0.55, z=1, lw=0.0, bend=0.5):
    """
    EN | Curved flow band from the vertical span [y0a, y0b] at x0 to [y1a, y1b] at x1, drawn
         with two cubic Bezier edges whose control points sit at a fraction `bend` of the
         horizontal distance.
    PT | Faixa curva de fluxo do trecho vertical [y0a, y0b] em x0 ao trecho [y1a, y1b] em x1,
         desenhada com duas bordas de Bezier cubica cujos pontos de controle ficam numa
         fracao `bend` da distancia horizontal.
    """
    dx = (x1 - x0) * bend
    verts = [(x0, y0a), (x0 + dx, y0a), (x1 - dx, y1a), (x1, y1a),
             (x1, y1b), (x1 - dx, y1b), (x0 + dx, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    p = PathPatch(Path(verts, codes), fc=blend(colour, alpha), ec="white" if lw else "none", lw=lw, zorder=z)
    ax.add_patch(p)
    return p


def block_arrow(ax, x0, y0, x1, y1, width, colour, head=None, z=4, alpha=1.0):
    """
    EN | Straight arrow whose shaft is `width` wide (data units): the width carries a rate.
         The head is `head` long and 1.9 times the shaft wide.
    PT | Seta reta cujo corpo tem `width` de largura (unidades de dados): a largura carrega
         uma taxa. A ponta tem comprimento `head` e 1,9 vez a largura do corpo.
    """
    v = np.array([x1 - x0, y1 - y0], dtype=float)
    L = float(np.hypot(*v))
    u = v / L
    n = np.array([-u[1], u[0]])
    head = head if head is not None else max(width * 1.1, 0.6)
    head = min(head, 0.6 * L)
    hw = width * 0.95
    body_end = np.array([x1, y1]) - u * head
    pts = [np.array([x0, y0]) + n * width / 2, body_end + n * width / 2, body_end + n * hw,
           np.array([x1, y1]), body_end - n * hw, body_end - n * width / 2, np.array([x0, y0]) - n * width / 2]
    p = Polygon(pts, closed=True, fc=colour, ec="none", zorder=z, alpha=alpha)
    ax.add_patch(p)
    return p


def downsample(path, printed_mm, dpi=DPI_IMAGE):
    """
    EN | Load a render and shrink it so that, printed `printed_mm` wide, it has `dpi` pixels per
         inch: more pixels only inflate the EPS and PDF without adding anything a reader can see.
    PT | Carrega um render e o reduz para que, impresso com `printed_mm` de largura, tenha `dpi`
         pixels por polegada: mais pixels so inflam o EPS e o PDF sem acrescentar nada visivel.
    """
    from PIL import Image
    im = Image.open(path).convert("RGBA")
    target = int(round(printed_mm / 25.4 * dpi))
    if im.width > target:
        h = int(round(im.height * target / im.width))
        im = im.resize((target, h), Image.LANCZOS)
    # EN/PT: flatten onto white, so that the PDF and EPS carry no soft mask | achata sobre branco, sem mascara no PDF e no EPS
    flat = Image.new("RGB", im.size, (255, 255, 255))
    flat.paste(im, mask=im.getchannel("A"))
    return np.asarray(flat)


def place_image(fig, path, left, bottom, width_mm, fig_w_mm, fig_h_mm, anchor="lb", zorder=5):
    """
    EN | Put a render on the figure with its lower-left (or centre, anchor="c") corner at
         (left, bottom) in figure fractions and a printed width of `width_mm`; returns the
         axes and the printed height in figure fractions.
    PT | Coloca um render na figura com o canto inferior esquerdo (ou o centro, anchor="c") em
         (left, bottom), em fracoes da figura, e largura impressa de `width_mm`; devolve os
         eixos e a altura impressa em fracoes da figura.
    """
    arr = downsample(path, width_mm)
    h_mm = width_mm * arr.shape[0] / arr.shape[1]
    w, h = width_mm / fig_w_mm, h_mm / fig_h_mm
    if anchor == "c":
        left, bottom = left - w / 2, bottom - h / 2
    ax = fig.add_axes([left, bottom, w, h], zorder=zorder)
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")
    return ax, h


def place_image_fit(fig, path, cx, cy, box_w_mm, box_h_mm, fig_w_mm, fig_h_mm, zorder=5):
    """
    EN | Put a render, as large as it can be while staying inside a box of box_w_mm x box_h_mm
         centred at (cx, cy) in figure fractions; returns the axes and the printed size in mm.
    PT | Coloca um render, o maior possivel dentro de uma caixa de box_w_mm x box_h_mm centrada em
         (cx, cy), em fracoes da figura; devolve os eixos e o tamanho impresso em mm.
    """
    from PIL import Image
    with Image.open(path) as im:
        aspect = im.height / im.width
    w_mm = min(box_w_mm, box_h_mm / aspect)
    arr = downsample(path, w_mm)
    h_mm = w_mm * arr.shape[0] / arr.shape[1]
    ax = fig.add_axes([cx - w_mm / fig_w_mm / 2, cy - h_mm / fig_h_mm / 2, w_mm / fig_w_mm, h_mm / fig_h_mm], zorder=zorder)
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")
    return ax, w_mm, h_mm


def mrna_glyph(ax, x, y, w, colour, lw=1.3):
    """EN/PT: a small messenger RNA: cap, wavy body, poly(A) tail | um mRNA pequeno: capuz, corpo ondulado e cauda poli(A)."""
    t = np.linspace(0, 1, 120)
    xs = x + w * t
    ys = y + 0.16 * w * np.sin(2 * np.pi * 2.5 * t) * (0.4 + 0.6 * np.sin(np.pi * t))
    ax.plot(xs, ys, color=colour, lw=lw, solid_capstyle="round", zorder=4)
    ax.plot([x], [y], "o", ms=3.6, color=colour, zorder=5)
    for k in range(5):
        ax.plot([x + w * (0.84 + 0.032 * k)] * 2, [y - 0.045 * w, y + 0.045 * w], color=colour, lw=0.9, zorder=4)


def coil_glyph(ax, x, y, w, colour, lw=1.3, seed=3):
    """EN/PT: a disordered chain, drawn as a smoothed random walk; stands for a monomer with no single structure | cadeia desordenada."""
    rng = np.random.default_rng(seed)
    pts = np.cumsum(rng.normal(0, 1, (26, 2)), axis=0)
    pts -= pts.mean(axis=0)
    pts /= np.abs(pts).max()
    k = np.ones(3) / 3
    sm = np.column_stack([np.convolve(np.r_[pts[0, i], pts[:, i], pts[-1, i]], k, mode="valid") for i in range(2)])
    t = np.linspace(0, 1, len(sm))
    tt = np.linspace(0, 1, 160)
    xs, ys = np.interp(tt, t, sm[:, 0]), np.interp(tt, t, sm[:, 1])
    ax.plot(x + xs * w * 0.5, y + ys * w * 0.5, color=colour, lw=lw, solid_capstyle="round", zorder=4)


def nucleus_glyph(ax, x, y, r, colour, over="#ffffff"):
    """EN/PT: a handful of molecules that gather: the sign for nucleation | poucas moleculas que se juntam: sinal de nucleacao."""
    pos = [(0.0, 0.0), (1.0, 0.25), (0.45, 0.95), (-0.5, 0.8), (-0.9, -0.1), (0.2, -0.95), (1.0, -0.6)]
    for i, (a, b) in enumerate(pos):
        ax.add_patch(plt_circle((x + a * r, y + b * r), r * 0.42, blend(colour, 0.25 + 0.75 * (i < 4), over)))


def plt_circle(xy, r, colour):
    from matplotlib.patches import Circle
    return Circle(xy, r, fc=colour, ec="white", lw=0.5, zorder=4)



def mirna_glyph(ax, x, y, w, colour, seed_colour="#d62839", lw=1.5):
    """EN/PT: a short miRNA: a gentle wave whose first stretch (the seed) is red | miRNA curto, com a semente em vermelho."""
    t = np.linspace(0, 1, 90)
    xs = x + w * t
    ys = y + 0.10 * w * np.sin(2 * np.pi * 1.6 * t)
    cut = int(0.30 * len(t))
    ax.plot(xs[:cut + 1], ys[:cut + 1], color=seed_colour, lw=lw + 0.3, solid_capstyle="round", zorder=4)
    ax.plot(xs[cut:], ys[cut:], color=colour, lw=lw, solid_capstyle="round", zorder=4)


def protein_glyph(ax, x, y, r, colour, over="#ffffff"):
    """EN/PT: a folded protein drawn as a few overlapping discs | proteina enovelada, alguns discos sobrepostos."""
    from matplotlib.patches import Circle
    for dx, dy, k in ((-0.55, -0.15, 0.78), (0.35, 0.25, 0.9), (0.55, -0.35, 0.66), (-0.1, 0.45, 0.62)):
        ax.add_patch(Circle((x + dx * r, y + dy * r), k * r * 0.62, fc=blend(colour, 0.88, over), ec="white", lw=0.5, zorder=4))


def layers_glyph(ax, x, y, w, colour, n=4, gap=0.34, ghost=True, over="#ffffff"):
    """EN/PT: layers of a fibril, the top one ghosted: where the next monomer lands | camadas de fibrila, a de cima translucida."""
    h = w * 0.10
    for i in range(n):
        ax.add_patch(Rectangle((x, y - i * w * gap), w, h, fc=blend(colour, 0.95 - 0.1 * i, over), ec="white", lw=0.4, zorder=4))
    if ghost:
        ax.add_patch(Rectangle((x, y + w * gap), w, h, fc=blend(colour, 0.22, over), ec=colour, lw=0.6, ls=(0, (2, 1.4)), zorder=4))


def pulse_glyph(ax, x, y, w, h, colour, lw=1.3, over="#ffffff"):
    """EN/PT: repeated boluses: sharp rise, slow decay | bolus repetidos: subida brusca e queda lenta."""
    xs, ys = [], []
    for k in range(3):
        t = np.linspace(0, 1, 16)
        xs += list(x + (k + t) * w / 3.0)
        ys += list(y + h * np.exp(-3.0 * t))
        xs.append(x + (k + 1) * w / 3.0)
        ys.append(y)
    ax.plot(xs, ys, color=colour, lw=lw, zorder=4, solid_joinstyle="miter")
    ax.plot([x, x + w], [y, y], color=blend(colour, 0.5, over), lw=0.6, zorder=3)


def chip(ax, x, y, text, colour, soft, fontsize=6.5, pad=0.55):
    """EN/PT: small rounded label, soft fill and coloured edge | rotulo pequeno, fundo suave e borda colorida."""
    t = ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color=vs.INK, fontweight="bold", zorder=6)
    ax.figure.canvas.draw()
    rend = ax.figure.canvas.get_renderer()
    bb = t.get_window_extent(rend).transformed(ax.transData.inverted())
    w, h = bb.width + 2 * pad * 1.4, bb.height + 2 * pad * 0.7
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={h / 2}", fc=soft, ec=colour,
                                lw=0.8, zorder=5))
    return w
