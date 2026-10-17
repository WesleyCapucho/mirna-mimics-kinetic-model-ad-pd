#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Shared helper for producing every figure twice, once in English and once in
     Brazilian Portuguese.
PT | Auxiliar compartilhado para produzir cada figura duas vezes, uma em ingles e uma em
     portugues do Brasil.

EN | Why separate versions rather than bilingual labels. Until now the figures carried
     both languages in one string ("hours | horas"), which keeps a single file but reads
     badly in either language and cannot be dropped into a manuscript without editing.
     A journal figure has one language. This module makes each script emit
     <name>.en.png and <name>.pt-BR.png from the same code path, so the two versions
     cannot drift apart: there is one plotting function and one set of numbers, and only
     the strings change.
PT | Por que versoes separadas em vez de rotulos bilingues. Ate agora as figuras traziam
     os dois idiomas numa string so ("hours | horas"), o que mantem um arquivo unico mas
     le mal nos dois idiomas e nao entra num manuscrito sem edicao. Figura de periodico
     tem um idioma. Este modulo faz cada script emitir <nome>.en.png e <nome>.pt-BR.png
     pelo mesmo caminho de codigo, de modo que as duas versoes nao podem divergir: ha uma
     funcao de plotagem e um conjunto de numeros, e so as strings mudam.

EN | Numbers are never translated. Decimal separators stay as Python writes them, because
     a figure whose numbers change with the language is a figure whose numbers cannot be
     checked against the tables. Only words are switched.
PT | Numeros nunca sao traduzidos. Separadores decimais ficam como o Python escreve,
     porque uma figura cujos numeros mudam com o idioma e uma figura cujos numeros nao
     podem ser conferidos contra as tabelas. So as palavras mudam.
"""

LANGS = ("en", "pt-BR")


def t(lang, en, pt):
    """EN/PT: pick the English or the Portuguese string."""
    return en if lang == "en" else pt


def fig_path(directory, stem, lang, ext="png"):
    """
    EN | Path for one language's copy of a figure: results/figures/<stem>.<lang>.<ext>.
    PT | Caminho da copia de uma figura num idioma: results/figures/<stem>.<idioma>.<ext>.
    """
    return f"{directory}/{stem}.{lang}.{ext}"


def both(directory, stem, ext="png"):
    """EN/PT: iterate over (lang, path) for every language."""
    for lang in LANGS:
        yield lang, fig_path(directory, stem, lang, ext)


def _install_vector_twin():
    """
    EN | Journals ask for vector artwork where possible. Every figure script saves a PNG
         through Figure.savefig; this patch writes a PDF with embedded fonts next to each
         PNG from the same figure object, so the two files cannot differ. It is installed
         when this module is imported, which every figure script does.
    PT | Periodicos pedem arte vetorial quando possivel. Todo script de figura grava um PNG
         por Figure.savefig; este remendo grava um PDF com fontes embutidas ao lado de cada
         PNG, a partir do mesmo objeto de figura, de modo que os dois arquivos nao diferem.
         Instalado na importacao deste modulo, que todo script de figura faz.
    """
    try:
        import matplotlib
        from matplotlib.figure import Figure
    except ImportError:
        return
    if getattr(Figure.savefig, "_vector_twin", False):
        return
    matplotlib.rcParams["pdf.fonttype"] = 42
    original = Figure.savefig

    def savefig(self, fname, *args, **kwargs):
        is_png = isinstance(fname, str) and fname.lower().endswith(".png")
        try:
            import _journal_fit as jf
        except ImportError:
            jf = None
        main = bool(jf and is_png and jf.is_main(fname))
        if main:
            jf.fit(self)                      # main-text figures: at most 174 x 234 mm
        result = original(self, fname, *args, **kwargs)
        if is_png:
            kw = {k: v for k, v in kwargs.items() if k not in ("dpi", "format", "facecolor")}
            if "facecolor" in kwargs:
                kw["facecolor"] = kwargs["facecolor"]
            pdf = fname[:-4] + ".pdf"
            original(self, pdf, *args, format="pdf", **kw)
            if main:
                jf.record(fname, jf.measure(self))
                jf.eps_from_pdf(pdf)          # EPS twin for the journal
        return result

    savefig._vector_twin = True
    Figure.savefig = savefig


_install_vector_twin()
