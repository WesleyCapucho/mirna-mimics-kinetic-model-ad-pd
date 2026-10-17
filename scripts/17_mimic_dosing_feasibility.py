#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | What the measured half-lives imply for dosing a microRNA mimic: the peak-to-average
     penalty, the fraction of a dosing interval spent near the intended effect, and how
     much a miR-7 mimic would have to be stabilised to be dosed like an ordinary one.
PT | O que as meias-vidas medidas implicam para a dosagem de um mimetico de microRNA: a
     penalidade de pico sobre media, a fracao do intervalo passada perto do efeito
     pretendido, e o quanto um mimetico de miR-7 teria de ser estabilizado para ser dosado
     como um comum.

EN | The question this answers, and why nobody has answered it with these numbers.
     Half the literature on miR-7 and miR-29 as therapeutic targets treats them as
     interchangeable cargo: pick the miRNA whose target matters and deliver a mimic. The
     measured half-lives say they are not interchangeable at all. miR-7 turns over with a
     1.7 h half-life (K020), the shortest in a genome-wide survey whose median is 34 h
     (K016), because the lncRNA Cyrano drives it to target-directed degradation. This
     script turns that single fact into the quantity a formulation chemist would ask for.
PT | A pergunta que isto responde, e por que ninguem a respondeu com estes numeros.
     Metade da literatura sobre miR-7 e miR-29 como alvos terapeuticos os trata como carga
     intercambiavel: escolha o miRNA cujo alvo importa e entregue um mimetico. As
     meias-vidas medidas dizem que eles nao sao nada intercambiaveis. O miR-7 se renova com
     meia-vida de 1,7 h (K020), a mais curta de um levantamento genomico cuja mediana e
     34 h (K016), porque o lncRNA Cyrano o leva a degradacao dirigida pelo alvo. Este
     script transforma esse unico fato na grandeza que um quimico de formulacao pediria.

EN | The mathematics, and why it needs no fitted parameter. Under repeated dosing at
     interval T of a species cleared with first-order constant d, the concentration above
     baseline reaches a steady cycle. Just after a dose it is D / (1 - exp(-d T)); averaged
     over the interval it is D / (d T). Their ratio

         R(dT) = d T / (1 - exp(-d T))

     is the peak-to-average penalty, and D cancels: it depends only on the product of the
     decay constant and the dosing interval. Nothing is fitted, and nothing depends on how
     much repression a given dose produces. That is what makes this comparable across
     miRNAs whose potencies are not known.
PT | A matematica, e por que ela nao precisa de parametro ajustado. Sob dosagem repetida
     em intervalo T de uma especie eliminada com constante de primeira ordem d, a
     concentracao acima do basal atinge um ciclo estacionario. Logo apos a dose ela vale
     D / (1 - exp(-d T)); na media do intervalo vale D / (d T). A razao entre elas

         R(dT) = d T / (1 - exp(-d T))

     e a penalidade de pico sobre media, e D se cancela: depende so do produto da constante
     de decaimento pelo intervalo. Nada e ajustado, e nada depende de quanta repressao uma
     dose produz. E isso que torna a comparacao possivel entre miRNAs cujas potencias nao
     se conhecem.

EN | The assumption, stated plainly. The delivered mimic is taken to turn over at the same
     first-order rate as the endogenous species. A chemically stabilised mimic would not,
     which is exactly why the last calculation here is the useful one: it reports how much
     stabilisation is needed, and so states a design target rather than a verdict.
PT | A suposicao, dita com todas as letras. O mimetico entregue e tratado como tendo a
     mesma renovacao de primeira ordem da especie endogena. Um mimetico quimicamente
     estabilizado nao teria, e e justamente por isso que o ultimo calculo aqui e o util:
     ele diz de quanta estabilizacao se precisa, e portanto enuncia uma meta de projeto e
     nao um veredito.

    python scripts/17_mimic_dosing_feasibility.py
"""

import csv
import json
import math
import os
import sys
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bilingual import LANGS, t as tr, fig_path

KINETICS = "data/extracted/kinetic_parameters.csv"
TAB_DIR = "results/tables"
FIG_DIR = "results/figures"
LN2 = math.log(2.0)

# EN/PT: the species compared, each by the param_id of its measured decay constant
SPECIES = [
    ("miR-7", "K020", "miR-7", "miR-7"),
    ("miR-7 upper bound", "K021", "miR-7 upper bound", "miR-7 limite superior"),
    ("miR-29b", "K017", "miR-29b", "miR-29b"),
    ("miR-29c", "K018", "miR-29c", "miR-29c"),
    ("median miRNA", "K016", "median miRNA", "miRNA mediano"),
]
REFERENCE = "median miRNA"          # EN/PT: the yardstick an ordinary mimic would meet
INTERVALS_H = [8.0, 12.0, 24.0, 48.0, 168.0]


def load_decays():
    table = {r["param_id"]: r for r in csv.DictReader(open(KINETICS, encoding="utf-8"))}
    out = OrderedDict()
    for key, pid, lab_en, lab_pt in SPECIES:
        r = table.get(pid)
        if r is None:
            sys.exit(f"EN/PT: parameter {pid} not found in {KINETICS}")
        if r["kind"] != "numeric":
            sys.exit(f"EN/PT: parameter {pid} is {r['kind']}, not a numeric value")
        d = float(r["value_si"])
        out[key] = dict(param_id=pid, decay_per_hour=d, half_life_hours=LN2 / d,
                        label_en=lab_en, label_pt=lab_pt)
    return out


def peak_to_average(d, T):
    """
    EN | Peak concentration divided by the interval average, at steady state. The limit as
         d T goes to zero is 1 (a species that never decays), and it grows without bound
         as the interval outruns the half-life.
    PT | Concentracao de pico dividida pela media do intervalo, no estado estacionario. O
         limite quando d T tende a zero e 1 (especie que nao decai), e cresce sem limite
         quando o intervalo ultrapassa a meia-vida.
    """
    x = d * T
    if x < 1e-9:
        return 1.0
    return x / (1.0 - math.exp(-x))


def fraction_above_half_peak(d, T):
    """
    EN | Fraction of a dosing interval during which the level is still above half its
         post-dose peak. That is one half-life, capped at the interval.
    PT | Fracao do intervalo em que o nivel ainda esta acima da metade do pico pos-dose.
         E uma meia-vida, limitada pelo intervalo.
    """
    return min(1.0, (LN2 / d) / T)


def interval_for_ratio(d, target_ratio):
    """
    EN | The dosing interval at which a species with decay constant d reaches a given
         peak-to-average ratio. Solved by bisection because R is monotone in T.
    PT | O intervalo de dosagem em que uma especie de constante d atinge uma dada razao
         pico sobre media. Resolvido por bissecao porque R e monotona em T.
    """
    if target_ratio <= 1.0:
        return 0.0
    lo, hi = 1e-6, 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if peak_to_average(d, mid) < target_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def stabilisation_needed(d_fast, d_slow, T):
    """
    EN | By how much the half-life of the fast species must be multiplied so that, at the
         same dosing interval, its peak-to-average penalty equals the slow species'.
         Because R depends only on d T, the answer is simply the ratio of decay constants:
         matching the penalty at a fixed T means matching d, so the half-life must rise by
         d_fast / d_slow. Reported anyway through the solver, so the claim is computed and
         not asserted.
    PT | Por quanto a meia-vida da especie rapida tem de ser multiplicada para que, no mesmo
         intervalo, sua penalidade de pico sobre media iguale a da especie lenta. Como R
         depende so de d T, a resposta e a razao das constantes: igualar a penalidade num T
         fixo significa igualar d, entao a meia-vida tem de subir por d_fast / d_slow.
         Ainda assim e obtida pelo solver, para que a afirmacao seja calculada e nao
         declarada.
    """
    target = peak_to_average(d_slow, T)
    t_equiv = interval_for_ratio(d_fast, target)
    # EN/PT: d_required is the decay the fast species would need to reach target at T
    d_required = d_fast * t_equiv / T
    return dict(target_ratio=target, required_decay_per_hour=d_required,
                required_half_life_hours=LN2 / d_required,
                fold_stabilisation=(LN2 / d_required) / (LN2 / d_fast))


def plot(dec, rows, path, lang):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))

    ax = axes[0]
    T = np.linspace(0.5, 72, 400)
    # EN | The two miR-7 rows are the two ends of one measurement, not two species, so
    #      the upper bound is drawn as a dashed line in the same colour: the gap between
    #      the pair is the uncertainty in the miR-7 half-life, and hiding it would make
    #      the penalty look better determined than the source allows. The right panel
    #      lists both for the same reason.
    # PT | As duas linhas de miR-7 sao os dois extremos de uma mesma medida, nao duas
    #      especies, entao o limite superior sai tracejado na mesma cor: a distancia entre
    #      o par e a incerteza da meia-vida do miR-7, e esconde-la faria a penalidade
    #      parecer mais bem determinada do que a fonte permite. O painel da direita lista
    #      os dois pela mesma razao.
    mir7_colour = None
    for key, info in dec.items():
        y = [peak_to_average(info["decay_per_hour"], x) for x in T]
        upper = key == "miR-7 upper bound"
        line, = ax.plot(T, y, lw=1.4 if upper else 1.8,
                        ls="--" if upper else "-",
                        color=mir7_colour if upper else None,
                        label=f"{tr(lang, info['label_en'], info['label_pt'])} "
                              f"({info['half_life_hours']:.1f} h)")
        if key == "miR-7":
            mir7_colour = line.get_color()
    ax.axvline(24, color="#888888", ls="--", lw=1)
    ax.text(24.6, ax.get_ylim()[1] * 0.93, tr(lang, "daily", "diário"),
            fontsize=8, color="#666666")
    ax.set_yscale("log")
    ax.set_xlabel(tr(lang, "dosing interval (hours)", "intervalo entre doses (horas)"))
    ax.set_ylabel(tr(lang, "peak / average", "pico / média"))
    ax.set_title(tr(lang, "Cost of instability: peak needed per unit of average effect",
                    "Custo da instabilidade: pico necessário por efeito médio"), fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    labels, vals = [], []
    for r in rows:
        if r["dosing_interval_hours"] != 24.0:
            continue
        info = dec[r["species"]]
        labels.append(tr(lang, info["label_en"], info["label_pt"]))
        vals.append(100.0 * r["fraction_of_interval_above_half_peak"])
    ypos = range(len(labels))
    ax.barh(list(ypos), vals, color="#3B6EA5", height=0.6)
    for i, v in zip(ypos, vals):
        ax.text(v + 1.2, i, f"{v:.0f}%", va="center", fontsize=9)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel(tr(lang, "per cent of a 24 h interval above half the post-dose peak",
                     "por cento de um intervalo de 24 h acima da metade do pico"))
    ax.set_title(tr(lang, "How much of the day the dose is still working",
                    "Quanto do dia a dose ainda está agindo"), fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    dec = load_decays()

    rows = []
    for key, info in dec.items():
        for T in INTERVALS_H:
            rows.append(OrderedDict([
                ("species", key), ("param_id", info["param_id"]),
                ("half_life_hours", round(info["half_life_hours"], 4)),
                ("dosing_interval_hours", T),
                ("peak_over_average", round(peak_to_average(info["decay_per_hour"], T), 4)),
                ("fraction_of_interval_above_half_peak",
                 round(fraction_above_half_peak(info["decay_per_hour"], T), 4)),
            ]))

    ref = dec[REFERENCE]["decay_per_hour"]
    daily = 24.0
    penalty = OrderedDict()
    for key, info in dec.items():
        if key == REFERENCE:
            continue
        r_sp = peak_to_average(info["decay_per_hour"], daily)
        r_ref = peak_to_average(ref, daily)
        st = stabilisation_needed(info["decay_per_hour"], ref, daily)
        penalty[key] = OrderedDict([
            ("half_life_hours", round(info["half_life_hours"], 4)),
            ("peak_over_average_daily", round(r_sp, 4)),
            ("reference_peak_over_average_daily", round(r_ref, 4)),
            ("penalty_versus_reference", round(r_sp / r_ref, 4)),
            ("interval_matching_reference_penalty_hours",
             round(interval_for_ratio(info["decay_per_hour"], r_ref), 4)),
            ("half_life_required_for_daily_dosing_hours",
             round(st["required_half_life_hours"], 4)),
            ("fold_stabilisation_required", round(st["fold_stabilisation"], 4)),
        ])

    print("=" * 78)
    print("EN | Mimic dosing feasibility | PT | Viabilidade de dosagem de mimetico")
    print("=" * 78)
    print(f"  reference species | especie de referencia : {REFERENCE} "
          f"({dec[REFERENCE]['half_life_hours']:.1f} h, {dec[REFERENCE]['param_id']})")
    print("\n  Daily dosing (24 h):")
    print(f"    {'species':22s} {'t_half':>8s} {'peak/avg':>9s} {'penalty':>9s} "
          f"{'interval to match':>18s} {'stabilisation':>14s}")
    for key, p in penalty.items():
        print(f"    {key:22s} {p['half_life_hours']:7.1f}h {p['peak_over_average_daily']:9.2f} "
              f"{p['penalty_versus_reference']:8.1f}x "
              f"{p['interval_matching_reference_penalty_hours']:16.2f}h "
              f"{p['fold_stabilisation_required']:12.1f}x")
    print(f"\n  Fraction of a 24 h interval above half the post-dose peak:")
    for r in rows:
        if r["dosing_interval_hours"] == 24.0:
            print(f"    {r['species']:22s} {100 * r['fraction_of_interval_above_half_peak']:5.1f}%")

    with open(f"{TAB_DIR}/mimic_dosing_feasibility.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    payload = OrderedDict([
        ("model", "steady-state repeated dosing of a species cleared with first-order "
                  "constant d: peak/average = dT / (1 - exp(-dT)), independent of dose"),
        ("reference_species", REFERENCE),
        ("decay_constants", {k: dict(param_id=v["param_id"],
                                     decay_per_hour=v["decay_per_hour"],
                                     half_life_hours=v["half_life_hours"])
                             for k, v in dec.items()}),
        ("daily_dosing_comparison", penalty),
        ("assumption_en",
         "The delivered mimic is assumed to be cleared at the same first-order rate as the "
         "endogenous species. A chemically stabilised mimic would not be, which is why the "
         "result is stated as a required fold stabilisation: it is a design target, not a "
         "verdict on feasibility."),
        ("assumption_pt",
         "Supoe-se que o mimetico entregue seja eliminado na mesma taxa de primeira ordem "
         "da especie endogena. Um mimetico quimicamente estabilizado nao seria, e por isso "
         "o resultado e enunciado como fator de estabilizacao necessario: e meta de "
         "projeto, nao veredito de viabilidade."),
        ("reading_en",
         "Dosed once a day, a miR-7 mimic that turns over like endogenous miR-7 must reach "
         "a peak roughly eight times higher than a mimic of an ordinary miRNA to hold the "
         "same average level, and it spends about 7% of each day above half its own peak. "
         "To be dosed daily on the same terms as an ordinary miRNA its half-life would have "
         "to rise twentyfold, from 1.7 h to about 34 h. The miR-29 paralogues sit between "
         "the two, and differ from each other, which is a further reason not to treat the "
         "family as one agent."),
        ("reading_pt",
         "Dosado uma vez por dia, um mimetico de miR-7 que se renova como o miR-7 endogeno "
         "precisa atingir um pico cerca de oito vezes maior que o mimetico de um miRNA comum "
         "para sustentar o mesmo nivel medio, e passa cerca de 7% de cada dia acima da metade "
         "do proprio pico. Para ser dosado diariamente nas mesmas condicoes de um miRNA comum, "
         "sua meia-vida teria de subir vinte vezes, de 1,7 h para cerca de 34 h. Os parálogos "
         "de miR-29 ficam entre os dois, e diferem entre si, o que e mais uma razao para nao "
         "tratar a familia como um agente so."),
    ])
    with open(f"{TAB_DIR}/mimic_dosing_feasibility.json", "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)

    for lang in LANGS:
        plot(dec, rows, fig_path(FIG_DIR, "mimic_dosing_feasibility", lang), lang)

    print(f"\nEN/PT -> {TAB_DIR}/mimic_dosing_feasibility.csv")
    print(f"EN/PT -> {TAB_DIR}/mimic_dosing_feasibility.json")
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG_DIR, 'mimic_dosing_feasibility', lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
