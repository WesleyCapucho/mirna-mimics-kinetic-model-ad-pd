#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Fit the main-text figures to the artwork limits of the target journal and record
     what was delivered. Installed through scripts/_bilingual.py, so no figure script
     needs to change.

     Limits applied (Molecular Neurobiology author instructions, read on 2026-10-09
     through a search tool because the Springer site is blocked from this environment;
     see docs/en/JOURNAL_COMPLIANCE.md):
       - width at most 174 mm and height at most 234 mm;
       - lettering in Helvetica or Arial, 8 to 12 pt at final size, consistent;
       - every line at least 0.1 mm (0.3 pt) wide;
       - vector art as EPS, fonts embedded.
     A figure that is larger than the page area is scaled down as a whole (text, lines,
     markers and arrow heads by the same factor), so its layout is unchanged, and any
     lettering that ends up below 6.5 pt is raised to 6.5 pt. The size
     actually delivered, the smallest lettering and the thinnest line are written to
     results/tables/journal_figure_spec.csv, and scripts/08_verify_consistency.py checks
     that file.

PT | Ajusta as figuras do texto principal aos limites de arte do periodico-alvo e
     registra o que foi entregue. Instalado por scripts/_bilingual.py, de modo que
     nenhum script de figura precisa mudar.

     Limites aplicados (instrucoes da Molecular Neurobiology, lidas em 2026-10-09 por
     uma ferramenta de busca porque o site da Springer esta bloqueado neste ambiente;
     ver docs/pt-BR/CONFORMIDADE_REVISTA.md):
       - largura de no maximo 174 mm e altura de no maximo 234 mm;
       - texto em Helvetica ou Arial, 8 a 12 pt no tamanho final, consistente;
       - toda linha com pelo menos 0,1 mm (0,3 pt);
       - arte vetorial em EPS, com fontes embutidas.
     Figura maior que a area da pagina e reduzida por inteiro (texto, linhas, marcadores
     e pontas de seta pelo mesmo fator), de modo que o desenho nao muda, e todo texto que
     ficar abaixo de 6,5 pt e elevado a 6,5 pt. O tamanho
     entregue, o menor texto e a linha mais fina vao para
     results/tables/journal_figure_spec.csv, e scripts/08_verify_consistency.py confere
     esse arquivo.
"""

import csv
import os
import shutil
import subprocess

import numpy as np
from matplotlib.collections import Collection
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.text import Annotation, Text

MAX_W_MM = 174.0
MAX_H_MM = 234.0
MIN_LINE_PT = 0.3
FLOOR_PT = 6.5                    # smallest lettering allowed in a main-text figure
PAD_IN = 0.1                      # bbox_inches="tight" default padding, each side

MAIN_STEMS = (
    "prisma_flow_diagram", "design_profile", "quadas2_summary",
    "forest_plot_primary", "robustness_specification", "molecular_evidence_map",
    # Article 2 (mechanistic study)
    "mech_model_schematic", "mech_provenance", "mech_alzheimer", "mech_parkinson", "mech_turnover",
    "mech_target_response", "mech_design_space", "mech_series_ranking", "mech_mechanism_algorithm",
)
SPEC_CSV = "results/tables/journal_figure_spec.csv"
FIELDS = ["file", "width_mm", "height_mm", "fits_174x234_mm", "fonts", "min_text_pt",
          "texts_below_8pt", "texts_below_6pt", "min_line_pt"]


def stem_of(path):
    """EN/PT: 'results/figures/<stem>.<lang>.png' -> '<stem>' | extrai o nome base."""
    return os.path.basename(path).split(".")[0]


def is_main(path):
    return stem_of(path) in MAIN_STEMS


def tight_inches(fig):
    """EN/PT: size in inches that savefig(bbox_inches='tight') will write | tamanho que o savefig vai gravar."""
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    return bb.width + 2 * PAD_IN, bb.height + 2 * PAD_IN


def _texts(fig):
    return [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]


def _tick_label_ids(fig):
    ids = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            ids.update(id(t) for t in axis.get_ticklabels(which="both"))
    return ids


def _map_tick_label_sizes(fig, f):
    """EN/PT: tick labels are rebuilt by matplotlib, so size them through tick_params | rotulos de eixo sao recriados; dimensiona por tick_params."""
    for ax in fig.axes:
        for axis, name in ((ax.xaxis, "x"), (ax.yaxis, "y")):
            labels = axis.get_ticklabels(which="major")
            if labels:
                ax.tick_params(axis=name, which="both", labelsize=f(labels[0].get_fontsize()))


def scale_figure(fig, s):
    """EN/PT: scale every size-bearing property by s, keep the layout | escala tudo que tem tamanho, mantendo o desenho."""
    ticks = _tick_label_ids(fig)
    _map_tick_label_sizes(fig, lambda x: x * s)
    for t in fig.findobj(Text):
        if id(t) not in ticks:
            t.set_fontsize(t.get_fontsize() * s)
    for ln in fig.findobj(Line2D):
        ln.set_linewidth(ln.get_linewidth() * s)
        ln.set_markersize(ln.get_markersize() * s)
        ln.set_markeredgewidth(ln.get_markeredgewidth() * s)
    for p in fig.findobj(Patch):
        p.set_linewidth(p.get_linewidth() * s)
        if isinstance(p, FancyArrowPatch):
            p.set_mutation_scale(p.get_mutation_scale() * s)
            p.shrinkA *= s
            p.shrinkB *= s
    for c in fig.findobj(Collection):
        lw = np.asarray(c.get_linewidth(), dtype=float)
        if lw.size:
            c.set_linewidth(lw * s)
        if hasattr(c, "get_sizes"):
            sz = np.asarray(c.get_sizes(), dtype=float)
            if sz.size:
                c.set_sizes(sz * s * s)
    for a in fig.findobj(Annotation):
        if a.arrow_patch is not None:
            a.arrow_patch.set_mutation_scale(a.arrow_patch.get_mutation_scale() * s)
            a.arrow_patch.set_linewidth(a.arrow_patch.get_linewidth() * s)
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * s, h * s)


def raise_floor(fig, floor=FLOOR_PT):
    """EN/PT: no lettering below the floor | nenhum texto abaixo do piso."""
    ticks = _tick_label_ids(fig)
    _map_tick_label_sizes(fig, lambda x: max(x, floor))
    for t in fig.findobj(Text):
        if id(t) not in ticks and t.get_text().strip() and t.get_fontsize() < floor:
            t.set_fontsize(floor)


def fit(fig):
    """EN/PT: scale the figure down if it exceeds 174 x 234 mm, then keep the lettering floor | reduz a figura se passar de 174 x 234 mm e mantem o piso de texto."""
    s_total = 1.0
    for _ in range(4):
        w, h = tight_inches(fig)
        s = min(1.0, (MAX_W_MM / 25.4 - 2 * PAD_IN) / (w - 2 * PAD_IN),
                (MAX_H_MM / 25.4 - 2 * PAD_IN) / (h - 2 * PAD_IN))
        if s >= 0.9999:
            break
        scale_figure(fig, s * 0.998)
        raise_floor(fig)
        s_total *= s
    raise_floor(fig)
    return s_total


def measure(fig):
    """EN/PT: width, height, fonts, smallest text and thinnest line of the figure | medidas da figura."""
    w, h = tight_inches(fig)
    texts = _texts(fig)
    sizes = [t.get_fontsize() for t in texts]
    lws = []
    for ln in fig.findobj(Line2D):
        if ln.get_visible() and ln.get_linestyle() not in ("None", "none", "", " ") and len(ln.get_xdata()) > 1:
            lws.append(ln.get_linewidth())
    for p in fig.findobj(Patch):
        if not p.get_visible() or p.get_linewidth() <= 0:
            continue
        ec = p.get_edgecolor()
        if len(ec) == 4 and ec[3] > 0:
            lws.append(p.get_linewidth())
    return {
        "width_mm": round(w * 25.4, 1),
        "height_mm": round(h * 25.4, 1),
        "fonts": "; ".join(sorted({t.get_fontname() for t in texts})),
        "min_text_pt": round(min(sizes), 2) if sizes else "",
        "texts_below_8pt": sum(1 for x in sizes if x < 8.0),
        "texts_below_6pt": sum(1 for x in sizes if x < 6.0),
        "min_line_pt": round(min(lws), 2) if lws else "",
    }


def record(path, m):
    """EN/PT: upsert one row of the spec table | grava/atualiza uma linha da tabela de especificacoes."""
    os.makedirs(os.path.dirname(SPEC_CSV), exist_ok=True)
    rows = {}
    if os.path.exists(SPEC_CSV):
        for r in csv.DictReader(open(SPEC_CSV, encoding="utf-8")):
            rows[r["file"]] = r
    fits = m["width_mm"] <= MAX_W_MM + 0.05 and m["height_mm"] <= MAX_H_MM + 0.05
    rows[path] = {"file": path, **m, "fits_174x234_mm": "yes" if fits else "no"}
    with open(SPEC_CSV, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS)
        wr.writeheader()
        for k in sorted(rows):
            wr.writerow({f: rows[k].get(f, "") for f in FIELDS})


def eps_from_pdf(pdf_path):
    """EN/PT: EPS (vector, fonts embedded) from the figure's PDF | EPS vetorial, fontes embutidas, a partir do PDF."""
    exe = shutil.which("pdftocairo")
    if not exe:
        return None
    out = pdf_path[:-4] + ".eps"
    subprocess.run([exe, "-eps", pdf_path, out], check=True)
    return out
