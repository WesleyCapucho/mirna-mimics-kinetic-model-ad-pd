#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | One figure out of the four structures scripts/13 already rendered, so a reader sees
     the mechanistic story in one look instead of four unconnected PNGs.
PT | Uma figura a partir das quatro estruturas que o scripts/13 ja renderizou, para que
     o leitor veja a historia mecanistica de uma vez, em vez de quatro PNGs soltos.

EN | This composites; it does not render. Every image here is the exact PNG scripts/13
     produced from the deposited coordinates - nothing is redrawn, and no PyMOL session
     is opened. Every caption fact (PDB id, method, resolution, the catalytic dyad
     residues, the protofilament count and spacing) is read from
     results/tables/structure_figure_provenance.json, which scripts/13 itself wrote from
     the coordinate files, so a caption cannot say something the render's own provenance
     does not. The only new content is the layout: which image sits where, and the one
     line of English or Portuguese grouping it by axis.
PT | Isto compoe; nao renderiza. Toda imagem aqui e o PNG exato que o scripts/13 produziu
     das coordenadas depositadas - nada e redesenhado, e nenhuma sessao do PyMOL e
     aberta. Todo fato de legenda (PDB id, metodo, resolucao, os residuos da diade
     catalitica, o numero e espacamento dos protofilamentos) e lido de
     results/tables/structure_figure_provenance.json, que o proprio scripts/13 escreveu
     a partir dos arquivos de coordenadas, entao uma legenda nao pode dizer algo que a
     propria proveniencia do render nao diz. O unico conteudo novo e o layout: qual
     imagem fica onde, e a linha em ingles ou portugues que a agrupa por eixo.

EN | Why no arrows between the panels. An earlier draft connected the four structures
     with arrows reading like a pathway (AGO2 -> BACE1 -> fibril), which overstates what
     a static crystal structure and a cryo-EM fibril actually show: BACE1 does not
     structurally become the Abeta fibril, it is the enzyme whose product aggregates
     into it, and no structure in this set captures that transition. The panels are
     therefore grouped by which axis they belong to (a coloured accent bar, not a
     directional arrow) and by what each one is common to, not chained as if one state
     turned into the next.
PT | Por que nao ha setas entre os paineis. Um rascunho anterior ligava as quatro
     estruturas com setas que se liam como via (AGO2 -> BACE1 -> fibrila), o que afirma
     mais do que uma estrutura cristalografica estatica e uma fibrila de crio-EM de fato
     mostram: a BACE1 nao se torna estruturalmente a fibrila de Abeta, e a enzima cujo
     produto se agrega nela, e nenhuma estrutura deste conjunto capta essa transicao. Os
     paineis ficam por isso agrupados pelo eixo a que pertencem (uma barra de cor, nao
     uma seta direcional) e pelo que cada um tem em comum, e nao encadeados como se um
     estado virasse o proximo.

    python scripts/20_structure_story_panel.py
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bilingual import LANGS, t as tr, fig_path

PROV = "results/tables/structure_figure_provenance.json"
STRUCT_DIR = "results/figures/structures"
FIG_DIR = "results/figures"
TAB_DIR = "results/tables"

AD_HUE, PD_HUE, SHARED_HUE = "#2a78d6", "#eb6834", "#8a8a86"


# EN/PT: the only two experimental methods this project's four depositions use. A
# translation map rather than a generic translator, so a fifth structure with a method
# not listed here fails loudly instead of rendering in the wrong language.
METHOD_PT = {
    "X-RAY DIFFRACTION": "Difração de raios X",
    "ELECTRON MICROSCOPY": "Microscopia eletrônica",
}


def caption(lang, prov, pid, extra=""):
    p = prov[pid]
    res = f"{p['resolution_angstrom']:.2g} Å" if p.get("resolution_angstrom") else "n/a"
    crit = p.get("resolution_method", "")
    if crit:
        res += f" ({crit.replace(' CUT-OFF', '').title().replace('Fsc', 'FSC')})"
    if lang == "en":
        method = p["method"].title()
    else:
        if p["method"] not in METHOD_PT:
            sys.exit(f"EN/PT: no Portuguese translation for method {p['method']!r} "
                     f"({pid}) - add it to METHOD_PT")
        method = METHOD_PT[p["method"]]
    head = f"PDB {pid} · {method} · {res}"
    role = tr(lang, p["role_en"], p["role_pt"])
    if extra:
        return f"{head}\n{role}\n{extra}"
    return f"{head}\n{role}"


def build_captions(lang, prov):
    dyad = prov["4D8C"]["catalytic_dyad_resi_in_this_deposition"]
    dyad_txt = tr(lang, f"Catalytic dyad: Asp{dyad['DTGS']} (DTGS) / Asp{dyad['DSGT']} (DSGT)",
                  f"Díade catalítica: Asp{dyad['DTGS']} (DTGS) / Asp{dyad['DSGT']} (DSGT)")

    def fibril_extra(pid):
        # EN/PT: the fact worth surfacing here is the count and spacing scripts/13 itself
        # measured from the coordinates; the deposition-vs-abstract resolution discrepancy
        # (6CU7 only) is condensed to one short clause rather than the internal dev-note
        # verbatim, which ran past the panel's width and was cut off mid-sentence.
        cs = prov[pid]["chain_split"]
        lo, hi = cs["layer_spacing_angstrom"]
        base = tr(lang,
                  f"{cs['n_chains']} chains, 2 protofilaments · layer spacing {lo}–{hi} Å",
                  f"{cs['n_chains']} cadeias, 2 protofilamentos · espaçamento entre "
                  f"camadas {lo}–{hi} Å".replace(".", ","))
        if pid == "6CU7":
            base += tr(lang, "\n(abstract reports 3.7 Å for both polymorphs)",
                       "\n(o resumo do artigo informa 3,7 Å para os dois polimorfos)")
        return base

    return {
        "6N4O": caption(lang, prov, "6N4O"),
        "4D8C": caption(lang, prov, "4D8C", dyad_txt),
        "6CU7": caption(lang, prov, "6CU7", fibril_extra("6CU7")),
        "5OQV": caption(lang, prov, "5OQV", fibril_extra("5OQV")),
    }


def panel(fig, rect, img_path, title, hue, cap_text):
    """EN/PT: one image cell - accent bar, image, caption - inside a figure-fraction rect."""
    x0, y0, w, h = rect
    bar_h = 0.028
    ax_bar = fig.add_axes([x0, y0 + h - bar_h, w, bar_h])
    ax_bar.set_facecolor(hue)
    ax_bar.set_xticks([]); ax_bar.set_yticks([])
    for s in ax_bar.spines.values():
        s.set_visible(False)
    ax_bar.text(0.015, 0.5, title, transform=ax_bar.transAxes, ha="left", va="center",
                fontsize=9.5, color="white", weight="bold")

    img_h = h - bar_h - 0.115
    ax_img = fig.add_axes([x0 + 0.01, y0 + 0.115, w - 0.02, img_h])
    ax_img.imshow(mpimg.imread(img_path))
    ax_img.axis("off")  # EN/PT: no frame on the image; the accent bar carries the border colour

    fig.text(x0 + w / 2, y0 + 0.008, cap_text, ha="center", va="bottom",
              fontsize=7.6, color="#2b2a28", linespacing=1.35)


def plot(prov, path, lang):
    cap = build_captions(lang, prov)
    fig = plt.figure(figsize=(13.5, 10.2))
    fig.text(0.5, 0.975,
              tr(lang, "The structural mechanism, in one figure",
                 "O mecanismo estrutural, numa figura só"),
              ha="center", fontsize=15, weight="bold")
    fig.text(0.5, 0.952,
              tr(lang, "Four deposited structures, composited from the renders in "
                       "results/figures/structures/ - no new render, no new claim",
                 "Quatro estruturas depositadas, compostas a partir dos renders em "
                       "results/figures/structures/ - nenhum render novo, nenhuma "
                       "afirmação nova"),
              ha="center", fontsize=8.8, color="#52514e", style="italic")

    gx, gy, gw, gh, m = 0.03, 0.045, 0.465, 0.435, 0.02
    panel(fig, (gx, gy + gh + m, gw, gh), f"{STRUCT_DIR}/ago2_guide_target.png",
          tr(lang, "SEED-PAIRING MECHANISM · COMMON TO BOTH AXES",
             "MECANISMO DE PAREAMENTO DA SEMENTE · COMUM AOS DOIS EIXOS"),
          SHARED_HUE, cap["6N4O"])
    panel(fig, (gx + gw + m, gy + gh + m, gw, gh),
          f"{STRUCT_DIR}/bace1_inhibitor_active_site.png",
          tr(lang, "ALZHEIMER'S AXIS · THE ENZYME miR-29 REPRESSES",
             "EIXO DE ALZHEIMER · A ENZIMA QUE O miR-29 REPRIME"),
          AD_HUE, cap["4D8C"])
    panel(fig, (gx, gy, gw, gh), f"{STRUCT_DIR}/abeta42_fibril_axis.png",
          tr(lang, "ALZHEIMER'S AXIS · AGGREGATED END STATE",
             "EIXO DE ALZHEIMER · ESTADO AGREGADO FINAL"),
          AD_HUE, cap["5OQV"])
    panel(fig, (gx + gw + m, gy, gw, gh), f"{STRUCT_DIR}/alpha_synuclein_fibril_axis.png",
          tr(lang, "PARKINSON'S AXIS · AGGREGATED END STATE",
             "EIXO DE PARKINSON · ESTADO AGREGADO FINAL"),
          PD_HUE, cap["6CU7"])

    fig.text(0.5, 0.008,
              tr(lang, "No structure here shows BACE1 becoming a fibril, or SNCA bound to "
                       "AGO2: this is illustration of separately deposited structures, not "
                       "a captured pathway. Full provenance: results/tables/structure_figure_provenance.json.",
                 "Nenhuma estrutura aqui mostra a BACE1 virando fibrila, ou a SNCA ligada "
                       "ao AGO2: é ilustração de estruturas depositadas separadamente, não "
                       "uma via capturada. Proveniência completa: "
                       "results/tables/structure_figure_provenance.json."),
              ha="center", fontsize=6.6, color="#6b6a66")

    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    if not os.path.exists(PROV):
        sys.exit(f"EN/PT: missing {PROV} - run scripts/13_structure_figures.py first")
    prov = json.load(open(PROV, encoding="utf-8"))
    needed = {
        "6N4O": f"{STRUCT_DIR}/ago2_guide_target.png",
        "4D8C": f"{STRUCT_DIR}/bace1_inhibitor_active_site.png",
        "6CU7": f"{STRUCT_DIR}/alpha_synuclein_fibril_axis.png",
        "5OQV": f"{STRUCT_DIR}/abeta42_fibril_axis.png",
    }
    for pid, fp in needed.items():
        if pid not in prov:
            sys.exit(f"EN/PT: {pid} missing from {PROV}")
        if not os.path.exists(fp):
            sys.exit(f"EN/PT: missing {fp} - run scripts/13_structure_figures.py first")

    for lang in LANGS:
        plot(prov, fig_path(FIG_DIR, "structure_story_panel", lang), lang)

    print("=" * 78)
    print("EN | Structure story panel | PT | Painel estrutural composto")
    print("=" * 78)
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG_DIR, 'structure_story_panel', lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
