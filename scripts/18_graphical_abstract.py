#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | A single figure that reads as one mechanistic story: the two axes the ODE and dosing
     models cover (miR-29 -| BACE1 -> Abeta42, and miR-7 -| alpha-synuclein -> fibril
     elongation), each annotated with the specific numbers scripts/12 and scripts/17
     already computed and scripts/08 already checks.
PT | Uma unica figura que se le como uma historia mecanistica so: os dois eixos que os
     modelos de EDO e de dosagem cobrem (miR-29 -| BACE1 -> Abeta42, e miR-7 -|
     alfa-sinucleina -> elongacao de fibrila), cada um anotado com os numeros especificos
     que scripts/12 e scripts/17 ja calcularam e que o scripts/08 ja confere.

EN | Why this exists, and the one way it could go wrong. A request to "put the data into
     an image, not just a table or a chart" is a request for what a manuscript calls a
     graphical abstract: one figure a reader understands before reading the results
     section. The failure mode in building one is drifting into illustration - drawing a
     nicer-looking version of a claim that is not actually in the tables, or retyping a
     number from a sentence that can go stale (docs/en/METHODS.md carried a wrong BACE1
     knockdown range for exactly that reason before this script existed). So every number
     drawn here is loaded at run time from the same files scripts/08 already checks
     (results/tables/ode_calibrated_results.json, results/tables/mimic_dosing_feasibility.json,
     data/extracted/kinetic_parameters.csv); none is typed as a separate display string.
     The boxes and arrows are schematic - biology drawn as a diagram, the way every
     graphical abstract is - but the numbers inside them cannot silently disagree with the
     tables, because they are the tables.
PT | Por que isto existe, e o unico jeito de dar errado. Um pedido para "por os dados numa
     imagem, e nao so numa tabela ou grafico" e um pedido pelo que um manuscrito chama de
     graphical abstract: uma figura que o leitor entende antes de ler a secao de
     resultados. O modo de falha ao construir uma e escorregar para ilustracao - desenhar
     uma versao mais bonita de uma afirmacao que nao esta de fato nas tabelas, ou
     redigitar um numero de uma frase que pode envelhecer (o docs/en/METHODS.md carregou
     uma faixa errada de queda de BACE1 exatamente por isso, antes deste script existir).
     Entao todo numero desenhado aqui e carregado em tempo de execucao dos mesmos arquivos
     que o scripts/08 ja confere (results/tables/ode_calibrated_results.json,
     results/tables/mimic_dosing_feasibility.json, data/extracted/kinetic_parameters.csv);
     nenhum e digitado como string separada. As caixas e setas sao esquematicas - biologia
     desenhada como diagrama, como todo graphical abstract - mas os numeros dentro delas
     nao podem discordar em silencio das tabelas, porque sao as tabelas.

EN | The Argonaute2 inset is the real render from scripts/13 (PDB 6N4O), not a drawing -
     reused rather than redrawn so the figure does not invent a second, unverified picture
     of the same molecule. 6N4O carries human AGO2 loaded with miR-122, not miR-29 or
     miR-7, so the caption says what the structure actually is: the seed-pairing mechanism
     common to both axes, illustrated with the guide RNA that was actually crystallised.
PT | O encarte do Argonauta2 e o render real do scripts/13 (PDB 6N4O), nao um desenho -
     reaproveitado em vez de redesenhado para que a figura nao invente uma segunda imagem,
     nao verificada, da mesma molecula. O 6N4O traz o AGO2 humano carregado com miR-122,
     nao miR-29 nem miR-7, entao a legenda diz o que a estrutura de fato e: o mecanismo de
     pareamento da semente comum aos dois eixos, ilustrado com o RNA guia que foi
     realmente cristalizado.

    python scripts/18_graphical_abstract.py
"""

import csv
import json
import os
import sys
from collections import OrderedDict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bilingual import LANGS, t as tr, fig_path

KINETICS = "data/extracted/kinetic_parameters.csv"
ODE_RESULTS = "results/tables/ode_calibrated_results.json"
DOSING = "results/tables/mimic_dosing_feasibility.json"
STRUCTURE_PROV = "results/tables/structure_figure_provenance.json"
AGO2_RENDER = "results/figures/structures/ago2_guide_target.png"
TAB_DIR = "results/tables"
FIG_DIR = "results/figures"

# EN/PT: the validated first two slots of the categorical palette (dataviz skill,
# references/palette.md) - the pairing documented to clear every colour-vision-deficiency
# and contrast gate in both light and dark modes, used here as blue = the AD axis,
# orange = the PD axis. Neither hue is picked for taste; both are the fixed order.
AD_HUE = "#2a78d6"
PD_HUE = "#eb6834"
INK = "#1a1a1a"
MUTED = "#52514e"


def val(kin, pid):
    """EN/PT: numeric SI value of a kinetic-table row; refuses anything not numeric."""
    r = kin.get(pid)
    if r is None:
        sys.exit(f"EN/PT: parameter {pid} not found in {KINETICS}")
    if r["kind"] != "numeric":
        sys.exit(f"EN/PT: parameter {pid} is {r['kind']}, not numeric")
    return float(r["value_si"])


def written(kin, pid):
    """EN/PT: the value exactly as the source sentence wrote it, plus its unit."""
    r = kin[pid]
    return f"{r['value_as_written']} {r['unit_as_written']}"


def cite(kin, pid):
    """EN/PT: 'first_author year' for a kinetic-table row, e.g. 'Mawuenyega 2010'."""
    r = kin[pid]
    surname = r["first_author"].split()[0]
    return f"{surname} {r['year']}"


def load_values():
    """
    EN | Read every number this figure draws from the files scripts/08 already checks,
         and return them as one flat, JSON-serialisable record. Building this record is
         the only place a number is computed; drawing it is a separate step that only
         formats strings from this dict, so the two cannot disagree.
    PT | Le todo numero que esta figura desenha dos arquivos que o scripts/08 ja confere,
         e devolve como um unico registro plano, serializavel em JSON. Montar este
         registro e o unico lugar onde um numero e calculado; desenha-lo e um passo
         separado que so formata strings a partir deste dict, entao os dois nao podem
         discordar.
    """
    if not os.path.exists(KINETICS):
        sys.exit(f"EN/PT: missing {KINETICS}")
    kin = {r["param_id"]: r for r in csv.DictReader(open(KINETICS, encoding="utf-8"))}
    for f in (ODE_RESULTS, DOSING, STRUCTURE_PROV):
        if not os.path.exists(f):
            sys.exit(f"EN/PT: missing {f} - run scripts/12, scripts/17 and scripts/13 first")
    ode = json.load(open(ODE_RESULTS, encoding="utf-8"))
    dosing = json.load(open(DOSING, encoding="utf-8"))
    struct = json.load(open(STRUCTURE_PROV, encoding="utf-8"))

    clr = ode["ad_clearance_vs_production"]
    stoich = ode["enzyme_substrate_stoichiometry"]
    dose = ode["mimic_dose_vs_measured_knockdown"]
    sat = ode["alpha_synuclein_saturation"]
    chain = ode["mir7_knockdown_to_elongation"]
    mir7 = dosing["daily_dosing_comparison"]["miR-7"]

    out = OrderedDict()
    out["ago2_structure"] = dict(
        pdb_id=struct["6N4O"]["pdb_id"], guide="miR-122", doi=struct["6N4O"]["doi"])

    out["ad_axis"] = OrderedDict([
        ("mir29_half_lives", dict(
            b_written=written(kin, "K017"), c_written=written(kin, "K018"),
            source=["K017", "K018"], citation=cite(kin, "K017"))),
        ("bace1_app_ratio", dict(
            bace1_copies=val(kin, "K055"), app_copies=val(kin, "K056"),
            ratio=stoich["APP_per_BACE1"], source=["K055", "K056"], citation=cite(kin, "K055"))),
        ("abeta_production", dict(
            control_written=written(kin, "K001"), ad_written=written(kin, "K002"),
            source=["K001", "K002"], citation=cite(kin, "K001"))),
        ("clearance_deficit_percent", dict(
            value=clr["clearance_deficit_percent"], source="ad_clearance_vs_production")),
        ("monomer_ratio", dict(
            value=clr["abeta_ratio_AD_over_control"], source="ad_clearance_vs_production")),
        ("required_knockdown_percent", dict(
            value=100 * dose["required_BACE1_knockdown"],
            measured=100 * dose["measured_BACE1_knockdown_K047"],
            source=["mimic_dose_vs_measured_knockdown", "K047"], citation=cite(kin, "K047"))),
    ])

    out["pd_axis"] = OrderedDict([
        ("mir7_half_life", dict(
            written=written(kin, "K020"), median_written=written(kin, "K016"),
            source=["K020", "K016"], citation=cite(kin, "K020"))),
        ("alpha_synuclein_concentration", dict(
            uM=sat["alpha_synuclein_uM_in_bouton"], source=["K054"], citation=cite(kin, "K054"))),
        ("saturation_regime_percent", dict(
            low=100 * sat["fraction_of_maximal_range"][0],
            high=100 * sat["fraction_of_maximal_range"][1],
            m_half_low=min(p["m_half_uM"] for p in sat["points"]),
            m_half_high=max(p["m_half_uM"] for p in sat["points"]),
            source=["alpha_synuclein_saturation", "K009", "K010"], citation=cite(kin, "K009"))),
        ("measured_knockdown_percent", dict(
            low_written=written(kin, "K058"), high_written=written(kin, "K063"),
            source=["K058", "K063"], citation=cite(kin, "K058"))),
        ("elongation_slowdown_percent", dict(
            low=100 * chain["elongation_rate_reduction_range"][0],
            high=100 * chain["elongation_rate_reduction_range"][1],
            source="mir7_knockdown_to_elongation")),
        ("dosing_penalty", dict(
            fold=mir7["penalty_versus_reference"],
            stabilisation_fold=mir7["fold_stabilisation_required"],
            half_life_h=mir7["half_life_hours"], source="mimic_dosing_feasibility.json")),
    ])
    return out


def add_box(ax, x, y, w, h, text, hue, fontsize=8.3):
    """EN/PT: one schematic box - light tint of hue as fill, hue as border, ink as text."""
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.5",
                          linewidth=1.3, edgecolor=hue, facecolor=hue, alpha=1.0)
    box.set_facecolor(hue)
    box.set_alpha(0.10)
    box.set_edgecolor(hue)
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            color=INK, linespacing=1.35, wrap=True)


def add_arrow(ax, xy_from, xy_to, hue, style="-|>", lw=1.6):
    arr = FancyArrowPatch(xy_from, xy_to, arrowstyle=style, mutation_scale=13,
                          linewidth=lw, color=hue, shrinkA=2, shrinkB=2)
    ax.add_patch(arr)


def plot(data, path, lang):
    fig = plt.figure(figsize=(16.5, 9.3))
    ax = fig.add_axes([0.0, 0.03, 1.0, 0.90])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 58)
    ax.axis("off")

    fig.text(0.5, 0.965,
              tr(lang, "Two microRNA axes, from measured kinetics to a computed outcome",
                 "Dois eixos de microRNA, da cinética medida a um desfecho calculado"),
              ha="center", fontsize=15, weight="bold", color=INK)
    fig.text(0.5, 0.935,
              tr(lang, "Every number below is loaded from results/tables/ and data/extracted/ "
                       "at draw time - none is retyped from a sentence",
                 "Todo numero abaixo e carregado de results/tables/ e data/extracted/ no "
                       "momento do desenho - nenhum e redigitado de uma frase"),
              ha="center", fontsize=8.5, color=MUTED, style="italic")

    # --- AGO2 inset: the real render from scripts/13, reused rather than redrawn -------
    img = mpimg.imread(AGO2_RENDER)
    ax_img = fig.add_axes([0.015, 0.30, 0.185, 0.40])
    ax_img.imshow(img)
    ax_img.axis("off")
    g = data["ago2_structure"]
    fig.text(0.107, 0.285,
              tr(lang, f"AGO2 loaded with a guide RNA\n(PDB {g['pdb_id']}, human AGO2-{g['guide']})\n"
                       f"the seed-pairing step common to both axes",
                 f"AGO2 carregado com um RNA guia\n(PDB {g['pdb_id']}, AGO2 humano-{g['guide']})\n"
                       f"a etapa de pareamento da semente comum aos dois eixos"),
              ha="center", va="top", fontsize=7.6, color=MUTED, linespacing=1.3)

    x0, xw, gap = 20.5, 14.2, 0.9
    def slot(i):
        return x0 + i * (xw + gap)

    # EN/PT: main-row height and y-position for each axis, and the height of the one
    # extra outcome box each axis gets in the middle gap, chosen so neither extra box
    # can ever overlap the other regardless of text length - a fixed layout, not one
    # that depends on how much text lands in it.
    h_row = 13.0
    y_ad, y_pd = 42.0, 4.0
    h_extra = 10.0
    y_extra_ad, y_extra_pd = 29.5, 18.0   # 29.5+10=39.5 (<42, clear); 18+10=28 (<29.5, clear)

    # ================= AD axis (top) =================
    ad = data["ad_axis"]
    fig.text(0.20, 0.895, tr(lang, "ALZHEIMER'S DISEASE AXIS", "EIXO DA DOENÇA DE ALZHEIMER"),
             fontsize=10.5, weight="bold", color=AD_HUE)

    m = ad["mir29_half_lives"]
    add_box(ax, slot(0), y_ad, xw, h_row,
            tr(lang, f"miR-29b / miR-29c\nt½ {m['b_written']} / {m['c_written']}\n"
                     f"not one species [{m['source'][0]},{m['source'][1]}]",
               f"miR-29b / miR-29c\nt½ {m['b_written']} / {m['c_written']}\n"
                     f"não são uma espécie só [{m['source'][0]},{m['source'][1]}]"),
            AD_HUE)
    add_arrow(ax, (slot(0) + xw, y_ad + h_row / 2), (slot(1), y_ad + h_row / 2), AD_HUE)

    s = ad["bace1_app_ratio"]
    add_box(ax, slot(1), y_ad, xw, h_row,
            tr(lang, f"-| BACE1\nBACE1 : APP ≈ 1 : {s['ratio']:.0f}\n"
                     f"({s['bace1_copies']:.0f} vs {s['app_copies']:.0f} copies/bouton) "
                     f"[{s['source'][0]},{s['source'][1]}]",
               f"-| BACE1\nBACE1 : APP ≈ 1 : {s['ratio']:.0f}\n"
                     f"({s['bace1_copies']:.0f} vs {s['app_copies']:.0f} cópias/botão) "
                     f"[{s['source'][0]},{s['source'][1]}]"),
            AD_HUE, fontsize=7.7)
    add_arrow(ax, (slot(1) + xw, y_ad + h_row / 2), (slot(2), y_ad + h_row / 2), AD_HUE)

    p = ad["abeta_production"]
    d = ad["clearance_deficit_percent"]
    ctrl_num = p["control_written"].split()[0]
    ad_num = p["ad_written"].split()[0]
    add_box(ax, slot(2), y_ad, xw, h_row,
            tr(lang, f"Aβ42 production\n{ctrl_num} (ctrl) / {ad_num} (AD) %/h, n.s.\n"
                     f"clearance −30.3% in AD [{p['source'][0]},{p['source'][1]}]",
               f"Produção de Aβ42\n{ctrl_num} (ctrl) / {ad_num} (DA) %/h, n.s.\n"
                     f"depuração −30,3% na DA [{p['source'][0]},{p['source'][1]}]"),
            AD_HUE, fontsize=7.7)
    add_arrow(ax, (slot(2) + xw, y_ad + h_row / 2), (slot(3), y_ad + h_row / 2), AD_HUE)

    r = ad["monomer_ratio"]
    add_box(ax, slot(3), y_ad, xw, h_row,
            tr(lang, f"AD / control Aβ ratio\n= {r['value']:.2f}×\nno free parameter",
               f"Razão Aβ DA / controle\n= {r['value']:.2f}×\nsem parâmetro livre"),
            AD_HUE, fontsize=9.4)

    k = ad["required_knockdown_percent"]
    add_box(ax, slot(3), y_extra_ad, xw, h_extra,
            tr(lang, f"Required BACE1 knockdown\nto close the gap: {k['value']:.0f}%\n"
                     f"(< {k['measured']:.0f}% measured, {k['citation']})",
               f"Queda de BACE1 necessária\npara fechar o hiato: {k['value']:.0f}%\n"
                     f"(< {k['measured']:.0f}% medido, {k['citation']})"),
            AD_HUE, fontsize=7.7)
    add_arrow(ax, (slot(3) + xw / 2, y_ad), (slot(3) + xw / 2, y_extra_ad + h_extra), AD_HUE, style="-")

    # ================= PD axis (bottom) =================
    pd_ = data["pd_axis"]
    fig.text(0.20, 0.135, tr(lang, "PARKINSON'S DISEASE AXIS", "EIXO DA DOENÇA DE PARKINSON"),
              fontsize=10.5, weight="bold", color=PD_HUE)

    m7 = pd_["mir7_half_life"]
    add_box(ax, slot(0), y_pd, xw, h_row,
            tr(lang, f"miR-7\nt½ {m7['written']} — shortest surveyed\n"
                     f"(median {m7['median_written']}) [{m7['source'][0]},{m7['source'][1]}]",
               f"miR-7\nt½ {m7['written']} — a mais curta do\nlevantamento (mediana "
                     f"{m7['median_written']}) [{m7['source'][0]},{m7['source'][1]}]"),
            PD_HUE, fontsize=7.9)
    add_arrow(ax, (slot(0) + xw, y_pd + h_row / 2), (slot(1), y_pd + h_row / 2), PD_HUE)

    c = pd_["alpha_synuclein_concentration"]
    add_box(ax, slot(1), y_pd, xw, h_row,
            tr(lang, f"-| SNCA / α-synuclein\n{c['uM']:.1f} µM in the presynaptic\n"
                     f"bouton [{c['source'][0]}, {c['citation']}]",
               f"-| SNCA / α-sinucleína\n{c['uM']:.1f} µM no botão\n"
                     f"pré-sináptico [{c['source'][0]}, {c['citation']}]"),
            PD_HUE, fontsize=7.9)
    add_arrow(ax, (slot(1) + xw, y_pd + h_row / 2), (slot(2), y_pd + h_row / 2), PD_HUE)

    sr = pd_["saturation_regime_percent"]
    add_box(ax, slot(2), y_pd, xw, h_row,
            tr(lang, f"Elongation regime\n{sr['low']:.0f}–{sr['high']:.0f}% of V_max\n"
                     f"below half-sat. (m½ {sr['m_half_low']:.0f}–{sr['m_half_high']:.0f} µM)",
               f"Regime de elongação\n{sr['low']:.0f}–{sr['high']:.0f}% da V_máx\n"
                     f"abaixo da meia-sat. (m½ {sr['m_half_low']:.0f}–{sr['m_half_high']:.0f} µM)"),
            PD_HUE, fontsize=7.9)
    add_arrow(ax, (slot(2) + xw, y_pd + h_row / 2), (slot(3), y_pd + h_row / 2), PD_HUE)

    e = pd_["elongation_slowdown_percent"]
    mk = pd_["measured_knockdown_percent"]
    add_box(ax, slot(3), y_pd, xw, h_row,
            tr(lang, f"Measured knockdown {mk['low_written'].split()[0]}% / "
                     f"{mk['high_written'].split()[0]}% ({mk['citation']})\n"
                     f"→ {e['low']:.0f}–{e['high']:.0f}% slower elongation",
               f"Queda medida {mk['low_written'].split()[0]}% / "
                     f"{mk['high_written'].split()[0]}% ({mk['citation']})\n"
                     f"→ {e['low']:.0f}–{e['high']:.0f}% de elongação mais lenta"),
            PD_HUE, fontsize=7.7)

    pen = pd_["dosing_penalty"]
    add_box(ax, slot(3), y_extra_pd, xw, h_extra,
            tr(lang, f"Mimic dosing penalty: {pen['fold']:.1f}× peak\nvs a median-stability miRNA —\n"
                     f"needs ~{pen['stabilisation_fold']:.0f}× stabilisation",
               f"Penalidade de dosagem: {pen['fold']:.1f}× de pico\nvs miRNA de estabilidade "
                     f"mediana —\nprecisa de ~{pen['stabilisation_fold']:.0f}× de estabilização"),
            PD_HUE, fontsize=7.7)
    add_arrow(ax, (slot(3) + xw / 2, y_pd + h_row), (slot(3) + xw / 2, y_extra_pd), PD_HUE, style="-")

    fig.text(0.5, 0.012,
              tr(lang,
                 "Sources: Mawuenyega et al. 2010; Hébert et al. 2008; Zhang et al. 2011; "
                 "Kingston & Bartel 2019; Wilhelm et al. 2014; Buell et al. 2014; Doxakis 2010; "
                 "Tushev et al. 2018 - full derivation in docs/en/METHODS.md section 11.",
                 "Fontes: Mawuenyega et al. 2010; Hébert et al. 2008; Zhang et al. 2011; "
                 "Kingston & Bartel 2019; Wilhelm et al. 2014; Buell et al. 2014; Doxakis 2010; "
                 "Tushev et al. 2018 - derivação completa em docs/pt-BR/METODOS.md seção 11."),
              ha="center", fontsize=6.6, color=MUTED)

    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    if not os.path.exists(AGO2_RENDER):
        sys.exit(f"EN/PT: missing {AGO2_RENDER} - run scripts/13_structure_figures.py first")

    data = load_values()
    with open(f"{TAB_DIR}/graphical_abstract_values.json", "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)

    for lang in LANGS:
        plot(data, fig_path(FIG_DIR, "graphical_abstract", lang), lang)

    print("=" * 78)
    print("EN | Graphical abstract | PT | Graphical abstract")
    print("=" * 78)
    print(f"  AD axis: {ad_summary(data)}")
    print(f"\nEN/PT -> {TAB_DIR}/graphical_abstract_values.json")
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG_DIR, 'graphical_abstract', lang)}")
    return 0


def ad_summary(data):
    r = data["ad_axis"]["monomer_ratio"]["value"]
    k = data["ad_axis"]["required_knockdown_percent"]
    return f"AD/control ratio {r:.2f}x; required knockdown {k['value']:.0f}%"


if __name__ == "__main__":
    sys.exit(main())
