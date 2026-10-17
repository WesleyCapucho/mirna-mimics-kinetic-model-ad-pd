#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Check figure PDFs against the artwork rules of the target journal (Molecular Neurobiology,
     Springer Nature), measured on the vector PDF of each figure rather than on the script that
     drew it:

       - page size, at most 174 x 234 mm;
       - lettering contrast, at least 4.5:1 for every text span against the colour that dominates
         the pixels under it (the span's own colour is read from the PDF; the page is rendered at
         150 dpi, its pixels are binned to steps of 8 levels, and the most frequent bin that is
         not the text colour is taken as the background);
       - thinnest stroked line, at least 0.3 pt;
       - effective resolution of every embedded image, at least 600 dpi;
     and reports, without failing on them: lettering size (the journal asks 8 to 12 pt, and
     Figs. 1-6 carry text from 6.5 pt, see README), fonts (Arial or Helvetica) and transparent
     objects. The exit status is 1 when a rule that is checked fails.

     With no argument it checks the PDFs of Figs. 1-6 in both languages and of Fig. 7. With
     file names it checks those. `--report FILE` also writes the text to FILE.

     The contrast method is the one that found the failures v1.2.1 fixed. In v1.2.0 the lowest
     value was 1.99:1 (the light-blue numerals of Fig. 1); the orange of PD, the grey notes and
     the white on the green of low risk lay between 2.8:1 and 3.5:1, and the blue of AD between
     3.9:1 and 4.2:1.

     Needs PyMuPDF and numpy (pip install pymupdf).

PT | Confere PDFs de figuras contra as regras de arte da revista-alvo (Molecular Neurobiology,
     Springer Nature), medidas no PDF vetorial de cada figura e nao no script que a desenhou:

       - tamanho da pagina, no maximo 174 x 234 mm;
       - contraste do texto, ao menos 4,5:1 para cada trecho de texto contra a cor que domina os
         pixels sob ele (a cor do trecho e lida do PDF; a pagina e desenhada a 150 dpi, seus
         pixels sao agrupados em passos de 8 niveis, e o grupo mais frequente que nao seja a cor
         do texto e tomado como fundo);
       - linha tracada mais fina, ao menos 0,3 pt;
       - resolucao efetiva de toda imagem embutida, ao menos 600 dpi;
     e informa, sem reprovar por isso: tamanho do texto (a revista pede 8 a 12 pt, e as
     Figs. 1-6 tem texto a partir de 6,5 pt, ver README), fontes (Arial ou Helvetica) e objetos
     transparentes. O codigo de saida e 1 quando uma regra checada falha.

     Sem argumento confere os PDFs das Figs. 1-6 nos dois idiomas e o da Fig. 7. Com nomes de
     arquivo confere esses. `--report ARQUIVO` tambem grava o texto em ARQUIVO.

     O metodo de contraste e o que achou as falhas que a v1.2.1 corrigiu. Na v1.2.0 o menor valor
     era 1,99:1 (os numerais azul-claros da Fig. 1); o laranja da DP, as notas cinza e o branco
     sobre o verde de baixo risco ficavam entre 2,8:1 e 3,5:1, e o azul da DA entre 3,9:1 e 4,2:1.

     Precisa de PyMuPDF e numpy (pip install pymupdf).

    python scripts/45_check_figure_artwork.py [--report FILE] [figure.pdf ...]
"""

import os
import sys
from collections import Counter, defaultdict

import numpy as np

try:
    import pymupdf as fitz
except ImportError:                     # EN/PT: older PyMuPDF names the module fitz | PyMuPDF antigo chama o modulo fitz
    import fitz

FIG_DIR = "results/figures"
# EN/PT: manuscript numbering of the figures | numeracao das figuras do manuscrito
FIGURES = (("prisma_flow_diagram", 1), ("design_profile", 2), ("quadas2_summary", 3), ("forest_plot_primary", 4),
           ("robustness_specification", 5), ("molecular_evidence_map", 6))
FIG7 = "mirna_mechanism_schematic"

MAX_W_MM, MAX_H_MM = 174.0, 234.0
MIN_CONTRAST = 4.5
MIN_LINE_PT = 0.3
MIN_DPI = 600
DPI_TOLERANCE = 0.5            # EN/PT: matplotlib resamples a render to exactly 600 dpi, and the division lands at 599.99995 | o matplotlib reamostra a 600 dpi exatos
ZOOM = 150 / 72.0


def lum(c):
    """EN/PT: WCAG relative luminance of an (r, g, b) triple in 0..255 | luminancia relativa WCAG."""
    out = []
    for v in c:
        v = v / 255.0
        out.append(v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4)
    return 0.2126 * out[0] + 0.7152 * out[1] + 0.0722 * out[2]


def ratio(a, b):
    la, lb = lum(a), lum(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


def measure(path):
    """EN/PT: every measurement of one figure PDF, as a dict | todas as medidas de um PDF de figura."""
    doc = fitz.open(path)
    page = doc[0]
    m = {"file": path, "w_mm": page.rect.width * 25.4 / 72, "h_mm": page.rect.height * 25.4 / 72}
    pix = page.get_pixmap(matrix=fitz.Matrix(ZOOM, ZOOM), alpha=False)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    sizes, fonts, bad = [], Counter(), defaultdict(list)
    worst = 99.0
    for blk in page.get_text("dict")["blocks"]:
        for ln in blk.get("lines", []):
            for sp in ln.get("spans", []):
                t = sp["text"].strip()
                if not t:
                    continue
                sizes.append(sp["size"])
                fonts[sp["font"]] += 1
                col = sp["color"]
                fg = ((col >> 16) & 255, (col >> 8) & 255, col & 255)
                x0, y0, x1, y1 = [int(round(v * ZOOM)) for v in sp["bbox"]]
                reg = img[max(y0, 0):max(y1, y0 + 1), max(x0, 0):max(x1, x0 + 1)].reshape(-1, 3)
                if reg.size == 0:
                    continue
                cnt = Counter(map(tuple, ((reg // 8) * 8).tolist()))
                bg = next((c for c, _ in cnt.most_common(6) if ratio(c, fg) > 1.15), cnt.most_common(1)[0][0])
                r = ratio(fg, bg)
                worst = min(worst, r)
                if r < MIN_CONTRAST:
                    bad[(fg, bg)].append((round(r, 2), t[:30]))
    m.update(n_spans=len(sizes), min_pt=min(sizes) if sizes else float("nan"), max_pt=max(sizes) if sizes else float("nan"),
             below_8=sum(1 for s in sizes if s < 7.95), above_12=sum(1 for s in sizes if s > 12.05), fonts=dict(fonts),
             worst=worst, bad=bad, n_bad=sum(len(v) for v in bad.values()))
    lws = [d["width"] for d in page.get_drawings()
           if d.get("width") not in (None, 0) and d.get("stroke_opacity", 1) != 0 and d.get("color") is not None]
    m.update(thinnest=min(lws) if lws else float("nan"), n_stroked=len(lws))
    images = []
    for info in page.get_images(full=True):
        for r in page.get_image_rects(info[0]):
            images.append((info[2], info[3], r.width * 25.4 / 72, r.height * 25.4 / 72, info[2] / (r.width / 72)))
    m["images"] = images
    m["transparent"] = any((d.get("fill_opacity", 1) not in (1, None)) or (d.get("stroke_opacity", 1) not in (1, None))
                           for d in page.get_drawings())
    return m


def verdict(m):
    """EN/PT: the rules that fail, as text | as regras que falham, em texto."""
    fails = []
    if m["w_mm"] > MAX_W_MM + 0.05 or m["h_mm"] > MAX_H_MM + 0.05:
        fails.append("page over 174 x 234 mm")
    if m["n_bad"]:
        fails.append(f"{m['n_bad']} text spans below {MIN_CONTRAST}:1")
    if m["thinnest"] == m["thinnest"] and m["thinnest"] < MIN_LINE_PT:
        fails.append(f"a line thinner than {MIN_LINE_PT} pt")
    if any(im[4] < MIN_DPI - DPI_TOLERANCE for im in m["images"]):
        fails.append(f"an image below {MIN_DPI} dpi")
    return fails


def describe(m):
    """EN/PT: the lines of the report for one figure | as linhas do relatorio de uma figura."""
    out = [f"page: {m['w_mm']:.1f} x {m['h_mm']:.1f} mm | limit {MAX_W_MM:.0f} x {MAX_H_MM:.0f} mm -> "
           + ("OK" if m["w_mm"] <= MAX_W_MM + 0.05 and m["h_mm"] <= MAX_H_MM + 0.05 else "TOO LARGE"),
           f"text spans: {m['n_spans']} | sizes {m['min_pt']:.1f} to {m['max_pt']:.1f} pt | spans below 8 pt: {m['below_8']}"
           f" | above 12 pt: {m['above_12']} | fonts: {m['fonts']}",
           f"contrast: worst {m['worst']:.2f}:1 | spans below {MIN_CONTRAST}:1: {m['n_bad']}"]
    for (fg, bg), items in m["bad"].items():
        out.append("   text rgb%s on rgb%s -> %.2f:1 e.g. %s" % (fg, bg, items[0][0], [i[1] for i in items[:4]]))
    out.append(f"thinnest stroked line: {m['thinnest']:.2f} pt (limit {MIN_LINE_PT} pt) | stroked paths: {m['n_stroked']}")
    out.append("embedded images:" + ("" if m["images"] else " none"))
    for pw, ph, wmm, hmm, dpi in m["images"]:
        out.append("   %dx%d px placed at %.1f x %.1f mm -> %.0f dpi" % (pw, ph, wmm, hmm, dpi))
    out.append(f"transparent drawing objects: {m['transparent']}")
    names = " ".join(m["fonts"])
    if names and not all(("Arial" in f or "Helvetica" in f) for f in m["fonts"]):
        out.append(f"note: fonts other than Arial or Helvetica: {names}")
    fails = verdict(m)
    out.append("RULES CHECKED: " + ("PASS" if not fails else "FAIL (" + "; ".join(fails) + ")"))
    return out


def label(path):
    stem = os.path.basename(path).split(".")[0]
    for s, n in FIGURES:
        if s == stem:
            return f"Fig. {n}"
    return "Fig. 7" if stem == FIG7 else stem


def main():
    args = sys.argv[1:]
    report = None
    if "--report" in args:
        i = args.index("--report")
        if i + 1 >= len(args):
            sys.exit("EN/PT: usage: python scripts/45_check_figure_artwork.py [--report FILE] [figure.pdf ...]")
        report = args[i + 1]
        del args[i:i + 2]
    files = args or ([f"{FIG_DIR}/{s}.{lang}.pdf" for s, _ in FIGURES for lang in ("en", "pt-BR")]
                     + [f"{FIG_DIR}/{FIG7}.en.pdf"])
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        sys.exit(f"EN/PT: missing {missing}")
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    lines = ["=" * 78, "EN | Artwork check of the figure PDFs | PT | Conferencia de arte dos PDFs de figura", "=" * 78,
             f"PyMuPDF {getattr(fitz, 'VersionBind', '?')} | contrast >= {MIN_CONTRAST}:1, lines >= {MIN_LINE_PT} pt, "
             f"images >= {MIN_DPI} dpi, page <= {MAX_W_MM:.0f} x {MAX_H_MM:.0f} mm"]
    results = []
    for f in files:
        m = measure(f)
        results.append(m)
        lines += ["", f"######## {label(f)}: {f.replace(os.sep, '/')}"] + describe(m)
    failing = [m["file"] for m in results if verdict(m)]
    lines += ["", "=" * 78,
              f"{len(results)} PDFs | worst contrast over all text: {min(m['worst'] for m in results):.2f}:1 | "
              f"text spans below {MIN_CONTRAST}:1: {sum(m['n_bad'] for m in results)} | "
              f"spans below 8 pt: {sum(m['below_8'] for m in results)} of {sum(m['n_spans'] for m in results)} | "
              + ("all checked rules pass" if not failing else f"FAILING: {failing}")]
    text = "\n".join(lines)
    print(text)
    if report:
        with open(report, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text + "\n")
        print(f"EN/PT -> {report}")
    return 1 if failing else 0


if __name__ == "__main__":
    sys.exit(main())
