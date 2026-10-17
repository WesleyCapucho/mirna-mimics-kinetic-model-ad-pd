#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Check that the kinetic parameter table, the raw sources it was read from, the model
     results, the design-space grid, the structure renders and the figures still agree with
     one another.
PT | Verifica se a tabela de parametros cineticos, as fontes brutas de onde foi lida, os
     resultados dos modelos, a grade do espaco de projeto, as renders de estrutura e as
     figuras continuam concordando entre si.

EN | This repository is a reproducibility package: raw data, scripts, tables and figures.
     Nothing here is narrated in prose, so the failure mode to guard against is a silent
     drift between artefacts: a parameter edited by hand, a table regenerated from a stale
     input, a sentence quoting a range the scan no longer produces. This script recomputes
     each derived number from its source and fails loudly on any disagreement. Run it after
     every pipeline change, before committing.
PT | Este repositorio e um pacote de reprodutibilidade: dados brutos, scripts, tabelas e
     figuras. Nada aqui e narrado em prosa, entao o modo de falha a evitar e um descolamento
     silencioso entre artefatos: um parametro editado a mao, uma tabela regerada a partir de
     entrada velha, uma frase que cita um intervalo que a varredura nao produz mais. Este
     script recalcula cada numero derivado a partir da fonte e falha ruidosamente em qualquer
     divergencia. Rode apos toda mudanca no pipeline, antes de versionar.

    python scripts/08_verify_consistency.py
"""

import csv
import os
import json
import math
import re
import sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP

KINETICS = "data/extracted/kinetic_parameters.csv"
SCHWANHAUSSER = "data/raw/kinetics_2026/schwanhausser_2011_supplementary_table.xls"
TUSHEV = "data/raw/kinetics_2026/tushev_2018_table_S1.xls"
ODE_RESULTS = "results/tables/ode_calibrated_results.json"
DOSING_JSON = "results/tables/mimic_dosing_feasibility.json"
DOSING_CSV = "results/tables/mimic_dosing_feasibility.csv"
ANIM_JSON = "results/tables/mechanism_animations.json"
TARGET_JSON = "results/tables/mimic_target_response.json"
TARGET_CSV = "results/tables/mimic_target_response.csv"


def main():
    fails, checks = [], 0

    def check(ok, msg):
        nonlocal checks
        checks += 1
        if not ok:
            fails.append(msg)

    # EN | A claim about the run itself ("N checks pass" in the README) cannot be a check,
    #      because counting it changes the number it asserts; claim() records a failure
    #      without touching the counter.
    # PT | Uma afirmacao sobre a propria execucao nao pode ser uma verificacao, porque conta-la
    #      muda o numero que afirma; o claim() registra a falha sem mexer no contador.
    def claim(ok, msg):
        if not ok:
            fails.append(msg)

    print("EN | Verifying internal consistency of parameters, model results and figures")
    print("PT | Verificando consistencia interna de parametros, resultados dos modelos e figuras")
    print("=" * 76)

    # --- 8. Kinetic parameters carry their own provenance ------------------
    # EN | Same discipline as the accuracy table: a parameter is only usable if
    #      the sentence that supports it is stored with it, and if the number
    #      recorded is the number that sentence actually states. The second
    #      condition is the one that catches a transcription slip, so it is
    #      checked mechanically rather than trusted.
    # PT | Mesma disciplina da tabela de acuracia: um parametro so serve se a
    #      frase que o sustenta estiver guardada junto, e se o numero registrado
    #      for o numero que essa frase de fato declara. A segunda condicao e a
    #      que pega erro de transcricao, entao e checada mecanicamente e nao
    #      confiada.
    try:
        kin = list(csv.DictReader(open(KINETICS, encoding="utf-8")))
    except FileNotFoundError:
        kin = []
    if kin:
        kids = [r["param_id"] for r in kin]
        check(len(set(kids)) == len(kids), "kinetics: duplicate param_id")
        for r in kin:
            pid = r["param_id"]
            if r["kind"] == "declared_gap":
                # EN | A parameter that was searched for and not found. It must say how it
                #      was searched and must not borrow a source, a quote or a number -
                #      attaching the DOI of a paper that was rejected would make a gap look
                #      like evidence.
                # PT | Um parametro procurado e nao encontrado. Precisa dizer como foi
                #      buscado e nao pode emprestar fonte, citacao nem numero - anexar o DOI
                #      de um artigo rejeitado faria uma lacuna parecer evidencia.
                check(not r["value_si"].strip() and not r["verbatim_quote"].strip()
                      and not r["pmid"].strip() and not r["doi"].strip(),
                      f"kinetics {pid}: a declared gap must not carry a value, quote or source")
                check(len(r["note"].strip()) > 80,
                      f"kinetics {pid}: a declared gap must record how it was searched")
                continue
            if r["kind"] == "derived":
                # EN | A value this project COMPUTED from a source table rather than read
                #      in a sentence. It must name its source and document the derivation,
                #      and it must not claim a verbatim quote it does not have: quoting a
                #      sentence that never contained the number is how a computed figure
                #      starts looking like a measured one.
                # PT | Valor que este projeto CALCULOU de uma tabela-fonte em vez de ler
                #      numa frase. Precisa nomear a fonte e documentar a derivacao, e nao
                #      pode alegar citacao verbatim que nao tem: citar uma frase que nunca
                #      continha o numero e como uma cifra calculada passa a parecer medida.
                check(not r["verbatim_quote"].strip(),
                      f"kinetics {pid}: a derived value must not carry a verbatim quote")
                check(bool(r["source"].strip()) and bool(r["value_si"].strip()),
                      f"kinetics {pid}: a derived value needs a source and an SI value")
                check("DERIVED" in r["note"] or "derived" in r["note"],
                      f"kinetics {pid}: a derived value must document its derivation")
                continue
            check(bool(r["verbatim_quote"].strip()), f"kinetics {pid}: no verbatim quote")
            check(bool(r["pmid"].strip() or r["doi"].strip()), f"kinetics {pid}: no PMID or DOI")
            check(r["kind"] in ("numeric", "qualitative_constraint", "declared_gap", "derived"),
                  f"kinetics {pid}: unknown kind {r['kind']!r}")
            check(bool(r["species"].strip()), f"kinetics {pid}: species not recorded")
            if r["kind"] != "numeric":
                # EN | A measured absence ("undetectable") is a real constraint on the
                #      model but is not a number, and must not be stored as a fitted
                #      zero. Only the wording is checked against the source.
                # PT | Uma ausencia medida ("undetectable") e uma restricao real ao
                #      modelo, mas nao e um numero, e nao pode ser guardada como zero
                #      ajustado. So a redacao e conferida contra a fonte.
                check(r["value_as_written"].strip() in r["verbatim_quote"],
                      f"kinetics {pid}: qualitative value does not appear in its quote")
                check(not r["value_si"].strip(),
                      f"kinetics {pid}: a qualitative constraint must not carry an SI value")
                continue
            check(bool(r["unit_si"].strip()), f"kinetics {pid}: no SI unit")
            check(bool(r["condition"].strip()),
                  f"kinetics {pid}: no experimental condition - a rate constant without "
                  "its pH, temperature and buffer is not a usable parameter")
            try:
                float(r["value_si"])
            except ValueError:
                check(False, f"kinetics {pid}: value_si {r['value_si']!r} is not numeric")
                continue
            # EN | The number as printed in the source must appear in the source's own
            #      sentence. value_si may differ from it (unit conversion, percent to
            #      rate constant, scientific notation), so the two are checked apart:
            #      value_as_written is the audit trail, value_si is what the model runs.
            # PT | O numero como impresso na fonte precisa aparecer na frase da propria
            #      fonte. value_si pode diferir dele (conversao de unidade, porcentagem
            #      para constante, notacao cientifica), entao os dois sao checados em
            #      separado: value_as_written e a trilha de auditoria, value_si e o que
            #      o modelo roda.
            quote = r["verbatim_quote"]
            written = r["value_as_written"].strip()
            trimmed = written.rstrip("0").rstrip(".") if "." in written else written
            check(written in quote or trimmed in quote,
                  f"kinetics {pid}: value_as_written {written!r} does not appear in its "
                  "verbatim quote")

    # EN | The two derived rows and the K030 gap rest on one archived spreadsheet
    #      (Schwanhausser et al. 2011, supplementary table). They were first worked
    #      out by hand; here they are recomputed from the file, so that anyone can see
    #      the medians and the absence of BACE1/SNCA come from the data and not from
    #      a note.
    # PT | As duas linhas derivadas e a lacuna K030 dependem de uma planilha arquivada
    #      (Schwanhausser et al. 2011, tabela suplementar). Foram calculadas a mao
    #      primeiro; aqui sao recalculadas a partir do arquivo, para que qualquer um
    #      veja que as medianas e a ausencia de BACE1/SNCA vem dos dados e nao de uma
    #      nota.
    kin_by_id = {r["param_id"]: r for r in kin}
    if {"K030", "K040", "K041"} <= set(kin_by_id):
        try:
            import pandas as pd
            sch = pd.read_excel(SCHWANHAUSSER)
        except Exception as e:  # EN/PT: missing xlrd or missing file is a failure
            sch = None
            check(False, f"kinetics: cannot read {SCHWANHAUSSER}: {e}")
        if sch is not None:
            for pid, col in (("K040", "mRNA half-life average [h]"),
                             ("K041", "Protein half-life average [h]")):
                vals = pd.to_numeric(sch[col], errors="coerce").dropna()
                r = kin_by_id[pid]
                check(int(r["n"]) == len(vals),
                      f"kinetics {pid}: n={r['n']} but the table has {len(vals)} values")
                check(abs(float(r["value_as_written"]) - float(vals.median())) < 1e-9,
                      f"kinetics {pid}: median {vals.median()} != {r['value_as_written']}")
                check(abs(float(r["value_si"]) - math.log(2) / float(vals.median())) < 1e-6,
                      f"kinetics {pid}: value_si is not ln2 / median half-life")
            genes = {g.strip().lower() for cell in sch["Gene Names"].dropna()
                     for g in str(cell).split(";")}
            for g in ("bace1", "snca", "sncb", "sncg", "app"):
                check(g not in genes,
                      f"kinetics K030: {g} IS in the table - the declared gap is wrong")

    # EN | Mimic washout: the reported time must be the time for a 10-fold bolus
    #      (excess 9x) to decay within 10% of baseline, t = ln(90) * t_half / ln 2.
    #      This once used ln(9) and halved every washout time.
    # PT | Eliminacao do mimetico: o tempo reportado deve ser o de um bolus de 10x
    #      (excesso 9x) decair a menos de 10% do basal, t = ln(90) * t_meia / ln 2.
    #      Isto ja usou ln(9) e cortou pela metade todos os tempos de eliminacao.
    try:
        ode = json.load(open(ODE_RESULTS, encoding="utf-8"))
    except FileNotFoundError:
        ode = None
    if ode is not None:
        for lab, w in ode["mimic_washout"].items():
            expect = math.log(90.0) * w["half_life_hours"] / math.log(2.0)
            check(abs(w["hours_to_return_to_baseline"] - expect) < 1e-6 * expect,
                  f"ode washout {lab}: {w['hours_to_return_to_baseline']:.3f} h, "
                  f"expected {expect:.3f} h")

    # EN | The dose comparison is analytic, so it can be checked in closed form rather
    #      than trusted: a mimic factor m gives a fractional BACE1 knockdown of
    #      1 - (d + 3k)/(d + 3k m), and the inverse must return the measured knockdown.
    #      Both ends of each reported range must also be consistent with the free
    #      parameter ranges the same file declares.
    # PT | A comparacao de doses e analitica, entao pode ser checada em forma fechada em
    #      vez de aceita: um fator de mimetico m da uma queda fracional de BACE1 de
    #      1 - (d + 3k)/(d + 3k m), e a inversa tem de devolver a queda medida. As duas
    #      pontas de cada faixa reportada tambem precisam ser coerentes com as faixas de
    #      parametro livre que o mesmo arquivo declara.
    if ode is not None and "mimic_dose_vs_measured_knockdown" in ode:
        dose = ode["mimic_dose_vs_measured_knockdown"]
        fr = ode["free_parameters"]
        target = dose["measured_BACE1_knockdown_K047"]
        kin_by_pid = {r["param_id"]: r for r in kin}
        check(abs(target - float(kin_by_pid["K047"]["value_si"])) < 1e-12,
              "ode dose: measured knockdown does not match K047")

        def kd(d, k, m):
            return 1.0 - (d + 3.0 * k) / (d + 3.0 * k * m)

        # EN/PT: d_mRNA is measured now (K068), so the corners span k_repress only
        d_meas = dose["d_mRNA_BACE1_measured"]
        check(abs(d_meas - float(kin_by_id["K068"]["value_si"])) < 1e-12,
              "ode dose: the BACE1 mRNA decay used is not K068")
        corners = [(d_meas, k)
                   for k in (fr["k_repress"]["low"], fr["k_repress"]["high"])]

        # EN | The knockdown the clearance deficit requires is rebuilt here from the four
        #      Abeta42 rows (K001-K004), not read back from the JSON: bringing the AD
        #      monomer pool (k_prod / k_clear) back to the control pool means lowering
        #      BACE1 by 1 - 1/R. It must be one number, the same at every k_repress. An
        #      earlier grid search reported it as a 30-33% range.
        # PT | A queda que o deficit de depuracao exige e refeita aqui a partir das quatro
        #      linhas de Abeta42 (K001-K004), nao relida do JSON: levar o pool de monomero
        #      da AD (k_prod / k_clear) de volta ao do controle significa baixar a BACE1
        #      em 1 - 1/R. Tem de ser um numero so, igual em todo k_repress. Uma busca em
        #      grade anterior o reportava como faixa de 30-33%.
        ab = {pid: float(kin_by_pid[pid]["value_si"]) for pid in ("K001", "K002", "K003", "K004")}
        ratio_ab = (ab["K002"] / ab["K004"]) / (ab["K001"] / ab["K003"])
        kappa = 1.0 - 1.0 / ratio_ab
        check(abs(dose["required_BACE1_knockdown"] - kappa) < 1e-12,
              f"ode dose: required knockdown {dose['required_BACE1_knockdown']:.6f} != "
              f"1 - 1/R = {kappa:.6f} from K001-K004")
        check(abs(ode["ad_clearance_vs_production"]["required_bace1_knockdown"] - kappa) < 1e-12,
              "ode dose: the clearance block and the dose block disagree on the required knockdown")
        # EN/PT: the mimic factor that delivers kappa at each end of the declared k_repress
        #        range, from the inverse written out independently of the script
        def mimic_for(d, k, t):
            return 1.0 + t / (1.0 - t) * (d + 3.0 * k) / (3.0 * k)
        m_lo, m_hi = dose["mimic_fold_to_offset_clearance"]
        (_, k_lo), (_, k_hi) = corners
        check(abs(m_hi - mimic_for(d_meas, k_lo, kappa)) < 1e-9
              and abs(m_lo - mimic_for(d_meas, k_hi, kappa)) < 1e-9,
              f"ode dose: mimic fold range {m_lo:.4f}-{m_hi:.4f} does not match the closed "
              f"form at the ends of the k_repress range")
        check(abs(kd(d_meas, k_lo, m_hi) - kappa) < 1e-9 and abs(kd(d_meas, k_hi, m_lo) - kappa) < 1e-9,
              "ode dose: the mimic factor at an end of the range does not give the required knockdown")
        clr = ode["ad_clearance_vs_production"]
        check(clr["mimic_fold_needed_range"] == [m_lo, m_hi]
              or all(abs(a - b) < 1e-12 for a, b in zip(clr["mimic_fold_needed_range"], [m_lo, m_hi])),
              "ode dose: the clearance block and the dose block disagree on the mimic fold range")
        check(m_lo <= clr["mimic_fold_needed"] <= m_hi,
              "ode dose: the midpoint mimic fold lies outside its own range")
        check(clr["mimic_fold_ode_check_max_relative_error"] < 1e-6,
              "ode dose: the integrated AD pool does not reach the control pool at the closed-form fold")
        for m_end in dose["mimic_fold_to_reproduce_measured_knockdown"]:
            got = [kd(d, k, m_end) for d, k in corners]
            check(min(abs(g - target) for g in got) < 1e-9,
                  f"ode dose: mimic fold {m_end:.4f} does not give the measured knockdown")
        check(abs(dose["required_over_measured_knockdown"] - kappa / target) < 1e-12,
              "ode dose: the required-over-measured ratio is not kappa / K047")
        check(dose["required_knockdown_below_measured"] == (kappa < target),
              "ode dose: the below-measured verdict disagrees with the numbers")

    # EN | Wilhelm et al. print copy numbers and concentrations in the same table and
    #      state in the legend that the concentrations were computed over the synaptic
    #      volume MINUS the mitochondrial volume. Both halves of that statement are
    #      recomputed here: the alpha-synuclein split that K054 derives, and the volume
    #      implied by each copy number / concentration pair, which must agree with the
    #      two geometric rows read from a different figure of the same paper.
    # PT | Wilhelm et al. imprimem numero de copias e concentracao na mesma tabela e
    #      dizem na legenda que as concentracoes foram calculadas sobre o volume
    #      sinaptico MENOS o volume mitocondrial. As duas metades disso sao recalculadas
    #      aqui: a separacao da alfa-sinucleina que o K054 deriva, e o volume implicado
    #      por cada par copias / concentracao, que tem de bater com as duas linhas
    #      geometricas lidas de outra figura do mesmo artigo.
    if {"K050", "K052", "K053", "K054", "K057"} <= set(kin_by_id):
        AVOGADRO = 6.02214076e23
        ratio = float(kin_by_id["K053"]["value_si"])
        share = ratio / (1.0 + ratio)
        total_uM = 43.57                      # printed beside K052 in the same row
        check(f"{total_uM}" in kin_by_id["K052"]["verbatim_quote"],
              "kinetics K052: the quoted row no longer carries the 43.57 uM used by K054")
        check(abs(float(kin_by_id["K054"]["value_si"]) - total_uM * share) < 1e-3,
              f"kinetics K054: {kin_by_id['K054']['value_si']} != "
              f"{total_uM} x {share:.6f} = {total_uM * share:.4f}")
        v_cyto = float(kin_by_id["K050"]["value_si"]) - float(kin_by_id["K057"]["value_si"])
        for pid, conc_uM in (("K052", 43.57), ("K055", 0.77), ("K056", 41.96)):
            copies = float(kin_by_id[pid]["value_si"])
            implied = copies / (conc_uM * 1e-6 * AVOGADRO) * 1e15     # litres -> um^3
            check(abs(implied - v_cyto) / v_cyto < 0.02,
                  f"kinetics {pid}: copies and concentration imply {implied:.4f} um^3, "
                  f"but the bouton minus its mitochondria is {v_cyto:.4f} um^3")

    # EN | The saturation and stoichiometry blocks are pure arithmetic on measured rows,
    #      so they are recomputed here from the parameter table rather than trusted.
    # PT | Os blocos de saturacao e estequiometria sao aritmetica pura sobre linhas
    #      medidas, entao sao recalculados aqui a partir da tabela em vez de aceitos.
    if ode is not None and "alpha_synuclein_saturation" in ode:
        sat = ode["alpha_synuclein_saturation"]
        C = float(kin_by_id["K054"]["value_si"])
        check(abs(sat["alpha_synuclein_uM_in_bouton"] - C) < 1e-9,
              "ode saturation: concentration does not match K054")
        for pt in sat["points"]:
            mh = pt["m_half_uM"]
            check(abs(pt["fraction_of_maximal_elongation_rate"] - C / (mh + C)) < 1e-9,
                  f"ode saturation: fraction at m_half {mh} is not C/(m_half+C)")
            check(abs(pt["percent_rate_change_per_percent_concentration_change"]
                      - mh / (mh + C)) < 1e-9,
                  f"ode saturation: elasticity at m_half {mh} is not m_half/(m_half+C)")
            check(abs(pt["m_half_uM"] - float(kin_by_id["K009"]["value_si"]) * 1e6) < 1e-6
                  or abs(pt["m_half_uM"] - float(kin_by_id["K010"]["value_si"]) * 1e6) < 1e-6,
                  f"ode saturation: m_half {mh} is neither K009 nor K010")
        check(sat["below_half_saturation"] == (max(
            p["fraction_of_maximal_elongation_rate"] for p in sat["points"]) < 0.5),
            "ode saturation: the below-half-saturation verdict disagrees with its own numbers")
    if ode is not None and "enzyme_substrate_stoichiometry" in ode:
        st = ode["enzyme_substrate_stoichiometry"]
        nb, na = float(kin_by_id["K055"]["value_si"]), float(kin_by_id["K056"]["value_si"])
        check(abs(st["BACE1_copies_per_bouton"] - nb) < 1e-9
              and abs(st["APP_copies_per_bouton"] - na) < 1e-9,
              "ode stoichiometry: copy numbers do not match K055 and K056")
        check(abs(st["APP_per_BACE1"] - na / nb) < 1e-9,
              "ode stoichiometry: the ratio is not APP over BACE1")

    # EN | The three-source chain is arithmetic on measured rows, so every step is
    #      recomputed from the parameter table: the knockdown against its row, the
    #      concentration against K054, the saturation form against K009 and K010, and
    #      the additivity comparison against K058, K059 and K060.
    # PT | A cadeia de tres fontes e aritmetica sobre linhas medidas, entao cada passo e
    #      recalculado a partir da tabela: a queda contra a linha dela, a concentracao
    #      contra o K054, a forma de saturacao contra K009 e K010, e a comparacao de
    #      aditividade contra K058, K059 e K060.
    if ode is not None and "mir7_knockdown_to_elongation" in ode:
        ch = ode["mir7_knockdown_to_elongation"]
        C0 = float(kin_by_id["K054"]["value_si"])
        check(abs(ch["alpha_synuclein_uM_in_bouton"] - C0) < 1e-9,
              "ode chain: concentration does not match K054")
        known_kd = {float(kin_by_id[pid]["value_si"]) for pid in ("K058", "K063")}
        known_mh = {float(kin_by_id[pid]["value_si"]) * 1e6 for pid in ("K009", "K010")}
        for c in ch["cases"]:
            kd, mh = c["knockdown"], c["m_half_uM"]
            check(any(abs(kd - k) < 1e-12 for k in known_kd),
                  f"ode chain: knockdown {kd} is not a row of the parameter table")
            check(any(abs(mh - m) < 1e-6 for m in known_mh),
                  f"ode chain: m_half {mh} is neither K009 nor K010")
            C1 = C0 * (1.0 - kd)
            check(abs(c["alpha_synuclein_after_uM"] - C1) < 1e-9,
                  "ode chain: post-knockdown concentration is not C0 x (1 - knockdown)")
            expect = 1.0 - (C1 / (mh + C1)) / (C0 / (mh + C0))
            check(abs(c["elongation_rate_reduction"] - expect) < 1e-9,
                  f"ode chain: elongation reduction {c['elongation_rate_reduction']:.6f} "
                  f"!= {expect:.6f} from the saturation form")
            check(c["elongation_rate_reduction"] < kd,
                  "ode chain: a saturating curve cannot give a more than proportional "
                  "slowdown; this one did")
        check(ch["less_than_proportional_everywhere"] ==
              all(c["elongation_rate_reduction"] < c["knockdown"] for c in ch["cases"]),
              "ode chain: the less-than-proportional verdict disagrees with its own cases")
        ad = ch["additivity"]
        a_, b_ = float(kin_by_id["K058"]["value_si"]), float(kin_by_id["K059"]["value_si"])
        check(abs(ad["pair_measured"] - float(kin_by_id["K060"]["value_si"])) < 1e-12,
              "ode chain: the measured pair does not match K060")
        check(abs(ad["pair_if_independent"] - (1.0 - (1.0 - a_) * (1.0 - b_))) < 1e-12,
              "ode chain: the independent-action prediction is not 1 - (1-a)(1-b)")

    # EN | scripts/36 (what a repeated mimic dose does to the protein it targets). The
    #      closed forms are rebuilt here from the parameter table, not read back from the
    #      JSON, and the scenario table is held to the identities the model has to obey:
    #      the knockdown reached is the one asked for; the dose does not depend on the
    #      protein half-life, because over a periodic cycle the mean protein equals the
    #      mean mRNA; peak over mean excess is scripts/17's penalty; efficiency never
    #      exceeds 1; the dose rises with the dosing interval and falls with the stability
    #      of the delivered species; and a target that is out of reach stays out of reach
    #      at longer intervals. The internal verification of scripts/36 must also have
    #      passed.
    # PT | scripts/36 (o que uma dose repetida de mimetico faz com a proteina que ele
    #      mira). As formas fechadas sao refeitas aqui a partir da tabela de parametros,
    #      nao relidas do JSON, e a tabela de cenarios e submetida as identidades que o
    #      modelo precisa obedecer: a queda alcancada e a pedida; a dose nao depende da
    #      meia-vida da proteina, porque num ciclo periodico a proteina media iguala o
    #      mRNA medio; excesso de pico sobre medio e a penalidade do scripts/17; a
    #      eficiencia nunca passa de 1; a dose sobe com o intervalo e cai com a
    #      estabilidade da especie entregue; e uma meta fora de alcance continua fora em
    #      intervalos maiores. A verificacao interna do scripts/36 tambem precisa ter
    #      passado.
    try:
        tgt = json.load(open(TARGET_JSON, encoding="utf-8"))
        tgt_rows = list(csv.DictReader(open(TARGET_CSV, encoding="utf-8")))
    except FileNotFoundError:
        tgt, tgt_rows = None, []
    if tgt is not None and ode is not None and kin_by_id:
        def pv(pid):
            return float(kin_by_id[pid]["value_si"])

        def rel(a, b):
            return abs(a - b) / max(abs(b), 1e-300)

        def fold(kappa, rho):
            return 1.0 + kappa / (1.0 - kappa) * (1.0 + rho) / rho

        fr_ = ode["free_parameters"]["k_repress"]
        k_lo, k_hi = fr_["low"], fr_["high"]
        k_mid = math.sqrt(k_lo * k_hi)
        d_ad, d_pd = pv("K068"), pv("K066")
        r_ab = (pv("K002") / pv("K004")) / (pv("K001") / pv("K003"))
        kappa_ad = 1.0 - 1.0 / r_ab

        check(all(v["passed"] for v in tgt["verification"].values()),
              "mimic target: an internal verification of scripts/36 did not pass")
        check(tgt["intervals_h"] == [8.0, 12.0, 24.0, 48.0, 168.0],
              "mimic target: the dosing intervals are not those of scripts/17")
        check(len(tgt_rows) == tgt["scenarios"],
              f"mimic target: {len(tgt_rows)} rows, the JSON says {tgt['scenarios']}")
        check(sum(1 for r in tgt_rows if r["reachable"] != "True") == tgt["scenarios_unreachable"],
              "mimic target: the unreachable count disagrees with the table")

        # EN/PT: repression strength rebuilt from K068, K066, K065 and the declared k_repress range
        rs = tgt["repression_strength"]
        want = {("AD", "rho_scan_low"): 3.0 * k_lo / d_ad, ("AD", "rho_scan_mid"): 3.0 * k_mid / d_ad,
                ("AD", "rho_scan_high"): 3.0 * k_hi / d_ad, ("PD", "rho_scan_low"): k_lo / d_pd,
                ("PD", "rho_scan_mid"): k_mid / d_pd, ("PD", "rho_scan_high"): k_hi / d_pd}
        for (ax_, key_), v in want.items():
            check(rel(rs[ax_][key_], v) < 1e-12, f"mimic target: {ax_} {key_} is not n k / d from the table")
        check(abs(rs["PD"]["rho_floor_K065"] - pv("K065")) < 1e-12 and
              rs["PD"]["floor_below_scan"] == (pv("K065") < k_lo / d_pd),
              "mimic target: the rho floor from K065 does not match the row, or its position against the scan is wrong")
        for ax_, tmap in (("AD", {"required to restore control Abeta42 [K001-K004]": kappa_ad}),
                          ("PD", {"miR-7 alone, SNCA construct with 3'UTR [K058]": pv("K058"),
                                  "miR-7 + miR-153, cortical neurons [K063]": pv("K063")})):
            for lab, kap in tmap.items():
                t_ = rs[ax_]["targets"][lab]
                check(abs(t_["knockdown"] - kap) < 1e-12 and
                      rel(t_["proportional_limit_fold"], 1.0 / (1.0 - kap)) < 1e-12,
                      f"mimic target: {ax_} target '{lab}' does not match its table rows")
                for key_, rho_ in (("fold_at_scan_low", rs[ax_]["rho_scan_low"]),
                                   ("fold_at_scan_mid", rs[ax_]["rho_scan_mid"]),
                                   ("fold_at_scan_high", rs[ax_]["rho_scan_high"])):
                    check(rel(t_[key_], fold(kap, rho_)) < 1e-12,
                          f"mimic target: {ax_} '{lab}' {key_} is not the closed form")
        for lab, t_ in rs["PD"]["targets"].items():
            check(rel(t_["fold_at_floor"], fold(t_["knockdown"], pv("K065"))) < 1e-12,
                  f"mimic target: PD '{lab}' fold at the K065 floor is not the closed form")
        # EN/PT: the AD closed form must also agree with scripts/12's range for the same fold
        lo_, hi_ = ode["mimic_dose_vs_measured_knockdown"]["mimic_fold_to_offset_clearance"]
        check(rel(lo_, fold(kappa_ad, rs["AD"]["rho_scan_high"])) < 1e-9 and
              rel(hi_, fold(kappa_ad, rs["AD"]["rho_scan_low"])) < 1e-9,
              "mimic target: scripts/12 and scripts/36 disagree on the AD mimic fold at the ends of the scan")

        # EN/PT: nucleation scaling, with the two orders asked for by what each row says it is
        nuc = tgt["nucleation_scaling"]
        check("primary" in kin_by_id["K013"]["parameter"] and "secondary" in kin_by_id["K012"]["parameter"],
              "mimic target: K013 must be the primary and K012 the secondary nucleation order")
        check(nuc["order_primary_nucleation"] == pv("K013") and nuc["order_secondary_nucleation"] == pv("K012"),
              "mimic target: the nucleation orders are not read from K013 (primary) and K012 (secondary)")
        check(rel(nuc["monomer_ratio_AD_over_control"], r_ab) < 1e-12 and
              abs(nuc["required_monomer_reduction"] - kappa_ad) < 1e-12,
              "mimic target: the monomer ratio or the required reduction does not match K001-K004")
        for kind_, n_ in (("primary", pv("K013")), ("secondary", pv("K012"))):
            check(rel(nuc[f"{kind_}_nucleation_rate_ratio_AD_over_control"], r_ab ** n_) < 1e-12 and
                  abs(nuc[f"{kind_}_nucleation_rate_reduction_after_offset"]
                      - (1.0 - (1.0 - kappa_ad) ** n_)) < 1e-12,
                  f"mimic target: {kind_} nucleation scaling is not ratio^order")

        # EN/PT: the scenario table
        half = {pid: math.log(2.0) / pv(pid) for pid in ("K016", "K017", "K018", "K020", "K021")}
        decay = {pid: pv(pid) for pid in half}
        prot_ok = {"AD": {round(math.log(2.0) / pv("K025"), 6), round(math.log(2.0) / pv("K029"), 6)},
                   "PD": {round(math.log(2.0) / pv("K026"), 6), round(math.log(2.0) / pv("K029"), 6)}}
        target_ok = {"AD": [kappa_ad], "PD": [pv("K058"), pv("K063")]}
        bad = {k_: [] for k_ in ("rho", "target", "species", "protein", "reached", "eff", "peak",
                                 "pool", "linear", "unreach")}
        for r in tgt_rows:
            ax_ = r["axis"]
            n_sp, d_m = (3.0, d_ad) if ax_ == "AD" else (1.0, d_pd)
            dx = decay[r["param_id"]]
            T = float(r["interval_h"])
            kap = float(r["target_knockdown"])
            tag = f"{ax_}/{r['rho_case']}/{r['species']}/{T:g}h/{r['target'][:12]}"
            if rel(float(r["rho"]), n_sp * float(r["k_repress"]) / d_m) > 1e-9:
                bad["rho"].append(tag)
            if not any(abs(kap - t_) < 1e-12 for t_ in target_ok[ax_]):
                bad["target"].append(tag)
            if rel(float(r["half_life_h"]), half[r["param_id"]]) > 1e-9:
                bad["species"].append(tag)
            if round(float(r["protein_half_life_h"]), 6) not in prot_ok[ax_]:
                bad["protein"].append(tag)
            e_lin = fold(kap, float(r["rho"])) - 1.0
            if rel(float(r["linear_excess"]), e_lin) > 1e-9 or rel(float(r["linear_dose"]), e_lin * dx * T) > 1e-9:
                bad["linear"].append(tag)
            if r["reachable"] == "True":
                dose = float(r["dose"])
                if abs(float(r["knockdown_reached"]) - kap) > 1e-8:
                    bad["reached"].append(tag)
                eff = float(r["efficiency"])
                if rel(eff, float(r["linear_dose"]) / dose) > 1e-9 or eff > 1.0 + 1e-9:
                    bad["eff"].append(tag)
                x = dx * T
                if rel(float(r["peak_excess"]) / float(r["mean_excess"]), x / (1.0 - math.exp(-x))) > 1e-9 \
                        or rel(float(r["mean_excess"]), dose / x) > 1e-9:
                    bad["peak"].append(tag)
                if rel(float(r["peak_pool_over_baseline"]), 1.0 + float(r["peak_excess"])) > 1e-9 \
                        or rel(float(r["trough_pool_over_baseline"]),
                               1.0 + float(r["peak_excess"]) * math.exp(-x)) > 1e-9:
                    bad["pool"].append(tag)
            elif r["dose"] != "" or float(r["knockdown_reached"]) >= kap:
                bad["unreach"].append(tag)
        for k_, what in (("rho", "rho is not n k / d from the table"),
                         ("target", "target knockdown is not a table row (K001-K004, K058, K063)"),
                         ("species", "species half-life is not its table row"),
                         ("protein", "protein half-life is not a table row (K025, K026, K029)"),
                         ("linear", "linear benchmark is not the closed form"),
                         ("reached", "knockdown reached differs from the target"),
                         ("eff", "efficiency is not linear dose over dose, or exceeds 1"),
                         ("peak", "peak over mean excess is not scripts/17's ratio"),
                         ("pool", "peak or trough pool is not 1 + excess"),
                         ("unreach", "an unreachable row carries a dose or reaches its target")):
            check(not bad[k_], f"mimic target: {what}; first: {bad[k_][:3]}")

        # EN/PT: the K065 floor case is present for every PD scenario that the scan has
        floor_rows = [r for r in tgt_rows if r["rho_case"] == "K065 lower bound"]
        check(len(floor_rows) == 60 and all(r["axis"] == "PD" for r in floor_rows),
              f"mimic target: {len(floor_rows)} rows at the K065 floor, expected 60 (3 species x 5 intervals x 2 protein cases x 2 targets, all PD)")

        # EN | The rho sweep: the grid is geometric, every reachable row obeys the same
        #      identities as the scenario table, and the per-scenario summary is the
        #      minimum and the end points of the rows it summarises.
        # PT | A varredura de rho: a grade e geometrica, toda linha alcancavel obedece as
        #      mesmas identidades da tabela de cenarios, e o resumo por cenario e o minimo
        #      e as pontas das linhas que resume.
        try:
            sw_rows = list(csv.DictReader(open("results/tables/mimic_target_rho_sweep.csv", encoding="utf-8")))
        except FileNotFoundError:
            sw_rows = []
        sw = tgt.get("rho_sweep")
        check(bool(sw_rows) and sw is not None, "mimic target: the rho sweep table or its JSON block is missing")
        if sw_rows and sw is not None:
            grid = sorted({float(r["rho"]) for r in sw_rows})
            ratios_ = [b_ / a_ for a_, b_ in zip(grid, grid[1:])]
            check(len(grid) == sw["points"] and max(ratios_) / min(ratios_) < 1.0 + 1e-9,
                  "mimic target: the rho sweep grid is not geometric or its size disagrees with the JSON")
            check(abs(grid[0] - sw["rho_min"]) < 1e-12 and abs(grid[-1] / sw["rho_max"] - 1.0) < 1e-9,
                  "mimic target: the rho sweep ends do not match the JSON")
            bad_sw = []
            for r in sw_rows:
                if r["reachable"] != "True":
                    continue
                ax_ = r["axis"]
                dx = decay[r["param_id"]]
                T = float(r["interval_h"])
                kap = float(r["target_knockdown"])
                rho_ = float(r["rho"])
                dose_ = float(r["dose"])
                x = dx * T
                e_peak = dose_ / (1.0 - math.exp(-x))
                tag = f"{ax_}/{r['species']}/rho {rho_:.3g}"
                if rel(float(r["peak_pool_over_baseline"]), 1.0 + e_peak) > 1e-9:
                    bad_sw.append(tag + " peak")
                if rel(float(r["linear_dose"]), (fold(kap, rho_) - 1.0) * x) > 1e-9:
                    bad_sw.append(tag + " linear")
                eff_ = float(r["efficiency"])
                if rel(eff_, float(r["linear_dose"]) / dose_) > 1e-9 or eff_ > 1.0 + 1e-9:
                    bad_sw.append(tag + " efficiency")
                if abs(float(r["knockdown_reached"]) - kap) > 1e-8:
                    bad_sw.append(tag + " knockdown")
            check(not bad_sw, f"mimic target: rho sweep rows break an identity; first: {bad_sw[:3]}")
            for e_ in sw["summary"]:
                pts_ = [r for r in sw_rows if r["axis"] == e_["axis"] and r["species"] == e_["species"]
                        and r["target"] == e_["target"]]
                pts_.sort(key=lambda r: float(r["rho"]))
                ok_ = [float(r["peak_pool_over_baseline"]) for r in pts_ if r["reachable"] == "True"]
                check(ok_ and abs(min(ok_) - e_["min_peak_pool_over_baseline"]) < 1e-9 * max(ok_)
                      and abs(float(pts_[0]["peak_pool_over_baseline"] or "nan") - e_["peak_at_lowest_rho"]) < 1e-6 * e_["peak_at_lowest_rho"]
                      and abs(float(pts_[-1]["peak_pool_over_baseline"] or "nan") - e_["peak_at_highest_rho"]) < 1e-6 * e_["peak_at_highest_rho"],
                      f"mimic target: the sweep summary for {e_['axis']} {e_['species']} {e_['target'][:14]} does not match its rows")

        # EN/PT: dose independent of the protein half-life (mean protein = mean mRNA over a cycle)
        groups = {}
        for r in tgt_rows:
            key_ = (r["axis"], r["rho_case"], r["species"], r["interval_h"], r["target"])
            groups.setdefault(key_, []).append(r)
        miss = [k_ for k_, g_ in groups.items()
                if len({x["dose"] == "" for x in g_}) > 1
                or any(x["dose"] != "" and rel(float(x["dose"]), float(g_[0]["dose"])) > 1e-7 for x in g_)]
        check(not miss, f"mimic target: the dose depends on the protein half-life in {miss[:2]}")

        # EN/PT: dose rises with the interval, and unreachable stays unreachable
        series = {}
        for r in tgt_rows:
            if not r["protein_case"].startswith(("BACE1 in", "alpha-syn")):
                continue
            series.setdefault((r["axis"], r["rho_case"], r["species"], r["target"]), []).append(r)
        nonmono, flip = [], []
        for k_, g_ in series.items():
            g_.sort(key=lambda x: float(x["interval_h"]))
            doses = [float(x["dose"]) for x in g_ if x["dose"] != ""]
            if any(b_ <= a_ for a_, b_ in zip(doses, doses[1:])):
                nonmono.append(k_)
            flags = [x["reachable"] == "True" for x in g_]
            if any((not a_) and b_ for a_, b_ in zip(flags, flags[1:])):
                flip.append(k_)
        check(not nonmono, f"mimic target: dose does not rise with the interval in {nonmono[:2]}")
        check(not flip, f"mimic target: a target that is out of reach becomes reachable at a longer interval in {flip[:2]}")

        # EN/PT: the dose falls as the delivered species gets more stable
        stab = {}
        for r in tgt_rows:
            if r["dose"] == "" or not r["protein_case"].startswith(("BACE1 in", "alpha-syn")):
                continue
            stab.setdefault((r["axis"], r["rho_case"], r["interval_h"], r["target"]), []).append(r)
        notfall = []
        for k_, g_ in stab.items():
            g_.sort(key=lambda x: float(x["half_life_h"]))
            doses = [float(x["dose"]) for x in g_]
            if any(b_ > a_ * (1.0 + 1e-9) for a_, b_ in zip(doses, doses[1:])):
                notfall.append(k_)
        check(not notfall, f"mimic target: the dose does not fall with the species' stability in {notfall[:2]}")

    # EN | The Tushev rows are read straight out of a 24,435-row table, so each one is
    #      looked up again here by gene symbol and compared field by field, and the
    #      derived pooled decay constant is recomputed from all five BACE1 isoforms.
    #      The pooling is the part worth guarding: the abundance-weighted mean of the
    #      RATE constants is the right quantity, and the mean of the half-lives is not.
    # PT | As linhas do Tushev sao lidas de uma tabela de 24.435 linhas, entao cada uma e
    #      procurada de novo aqui pelo simbolo do gene e comparada campo a campo, e a
    #      constante de decaimento derivada e recalculada a partir das cinco isoformas de
    #      BACE1. O pooling e a parte que vale guardar: a media ponderada das CONSTANTES
    #      DE VELOCIDADE e a grandeza certa, e a media das meias-vidas nao e.
    if {"K066", "K067", "K068", "K069"} <= set(kin_by_id):
        try:
            import pandas as pd
            tu = pd.read_excel(TUSHEV, sheet_name="PASSData", header=0)
            tu.columns = [str(c) for c in tu.columns]
        except Exception as e:
            tu = None
            check(False, f"kinetics: cannot read {TUSHEV}: {e}")
        if tu is not None:
            sym = tu["gene.symbol"].astype(str).str.lower()
            hl = pd.to_numeric(tu["half.life[hours]"], errors="coerce")
            rpm = pd.to_numeric(tu["rpm.neuron.culture"], errors="coerce")
            for pid, gene, single in (("K066", "snca", True), ("K069", "app", True)):
                rowsg = tu[sym == gene]
                check(len(rowsg) == 1,
                      f"kinetics {pid}: {gene} has {len(rowsg)} isoforms, the row assumes one")
                if len(rowsg) == 1:
                    t_printed = float(hl[rowsg.index[0]])
                    check(abs(float(kin_by_id[pid]["value_as_written"]) - t_printed) < 1e-4,
                          f"kinetics {pid}: half-life {kin_by_id[pid]['value_as_written']} "
                          f"!= {t_printed} in the table")
                    check(abs(float(kin_by_id[pid]["value_si"]) - math.log(2) / t_printed) < 1e-7,
                          f"kinetics {pid}: value_si is not ln2 over the printed half-life")
            b = tu[sym == "bace1"]
            check(len(b) == 5, f"kinetics K067: bace1 has {len(b)} isoforms, the note says five")
            if len(b):
                dom = rpm[b.index].idxmax()
                check(abs(float(kin_by_id["K067"]["value_as_written"]) - float(hl[dom])) < 1e-4,
                      "kinetics K067: this is not the most abundant BACE1 isoform's half-life")
                k_eff = float((rpm[b.index] * (math.log(2) / hl[b.index])).sum()
                              / rpm[b.index].sum())
                check(abs(float(kin_by_id["K068"]["value_si"]) - k_eff) < 1e-7,
                      f"kinetics K068: pooled constant {kin_by_id['K068']['value_si']} != {k_eff:.8f}")
                # EN/PT: and it must NOT be ln2 over the mean of the half-lives
                wrong = math.log(2) / float((rpm[b.index] * hl[b.index]).sum() / rpm[b.index].sum())
                check(abs(float(kin_by_id["K068"]["value_si"]) - wrong) > 1e-4,
                      "kinetics K068: the pooled constant equals the mean-of-half-lives form, "
                      "which is the wrong average for a decaying pool")

    # --- 8e. Mimic dosing feasibility ---------------------------------------
    # EN | Closed form throughout, so every published number is recomputed from the decay
    #      constants in the parameter table. Two properties are asserted as properties of
    #      the function rather than of this run: the peak-to-average ratio never drops
    #      below one, and a species that decays faster always pays a larger penalty at the
    #      same interval. A version that violated either would be wrong regardless of what
    #      the data said.
    # PT | Tudo em forma fechada, entao cada numero publicado e recalculado a partir das
    #      constantes de decaimento da tabela de parametros. Duas propriedades sao exigidas
    #      como propriedades da funcao e nao desta execucao: a razao pico sobre media nunca
    #      cai abaixo de um, e uma especie que decai mais rapido sempre paga penalidade
    #      maior no mesmo intervalo. Uma versao que violasse qualquer uma estaria errada
    #      independentemente do que os dados dissessem.
    try:
        dos = json.load(open(DOSING_JSON, encoding="utf-8"))
        dos_rows = list(csv.DictReader(open(DOSING_CSV, encoding="utf-8")))
    except FileNotFoundError:
        dos, dos_rows = None, None
    if dos is not None:
        ln2 = math.log(2.0)

        def p2a(d, T):
            x = d * T
            return 1.0 if x < 1e-9 else x / (1.0 - math.exp(-x))

        for name, rec in dos["decay_constants"].items():
            row = kin_by_id.get(rec["param_id"])
            check(row is not None and abs(float(row["value_si"]) - rec["decay_per_hour"]) < 1e-12,
                  f"dosing: {name} decay does not match {rec['param_id']}")
            check(abs(rec["half_life_hours"] - ln2 / rec["decay_per_hour"]) < 1e-9,
                  f"dosing: {name} half-life is not ln2 over its decay constant")
        for r in dos_rows:
            d = ln2 / float(r["half_life_hours"])
            T = float(r["dosing_interval_hours"])
            check(abs(float(r["peak_over_average"]) - p2a(d, T)) < 1e-3,
                  f"dosing {r['species']} at {T} h: peak/average is not dT/(1-exp(-dT))")
            check(abs(float(r["fraction_of_interval_above_half_peak"])
                      - min(1.0, (ln2 / d) / T)) < 1e-3,
                  f"dosing {r['species']} at {T} h: time above half peak is not one half-life")
            check(float(r["peak_over_average"]) >= 1.0,
                  f"dosing {r['species']}: a peak below the average is impossible")
        ref = dos["decay_constants"][dos["reference_species"]]["decay_per_hour"]
        for name, p in dos["daily_dosing_comparison"].items():
            d = ln2 / p["half_life_hours"]
            check(abs(p["peak_over_average_daily"] - p2a(d, 24.0)) < 1e-3,
                  f"dosing {name}: daily peak/average disagrees")
            check(abs(p["penalty_versus_reference"]
                      - p2a(d, 24.0) / p2a(ref, 24.0)) < 1e-3,
                  f"dosing {name}: the penalty is not its ratio to the reference")
            check(d > ref and p["penalty_versus_reference"] > 1.0,
                  f"dosing {name}: a species decaying faster than the reference must pay "
                  "a penalty above one")
            check(abs(p["fold_stabilisation_required"] - d / ref) < 1e-2,
                  f"dosing {name}: matching the reference at a fixed interval means "
                  "matching its decay constant, so the stabilisation is d/d_ref")

    # --- 8g. Mechanism animations recompute the same simulation ------------
    # EN | An animation is a movie of a simulation, and the only thing worth checking
    #      about it automatically is whether it is a movie of THIS simulation: its final
    #      frame has to land on the same numbers ode_calibrated_results.json and
    #      mimic_dosing_feasibility.json already carry. This does not open the GIF files
    #      - matplotlib figures are not meaningfully diffable - it recomputes the
    #      manifest scripts/19 wrote its numbers from and compares that manifest to the
    #      already-verified tables, the same one step removed from the pixels that every
    #      other figure check in this file uses.
    # PT | Uma animacao e o filme de uma simulacao, e a unica coisa que vale a pena
    #      conferir automaticamente nela e se e o filme DESTA simulacao: seu ultimo
    #      quadro tem de cair nos mesmos numeros que o ode_calibrated_results.json e o
    #      mimic_dosing_feasibility.json ja guardam. Isto nao abre os arquivos GIF -
    #      figuras do matplotlib nao sao comparaveis de forma significativa - recalcula
    #      o manifesto do qual o scripts/19 tirou seus numeros e compara esse manifesto
    #      as tabelas ja verificadas, o mesmo passo removido dos pixels que toda outra
    #      checagem de figura deste arquivo usa.
    try:
        anim = json.load(open(ANIM_JSON, encoding="utf-8"))
    except FileNotFoundError:
        anim = None
    if anim is not None and ode is not None:
        clr = ode["ad_clearance_vs_production"]
        check(abs(anim["ad_monomer_accumulation"]["ratio_AD_over_control"]
                  - clr["abeta_ratio_AD_over_control"]) < 1e-3,
              "animation manifest: AD/control ratio disagrees with ode_calibrated_results.json")

        wo = ode["mimic_washout"]
        m = anim["mimic_washout"]
        check(m["fast_species"] in wo and m["slow_species"] in wo,
              "animation manifest: washout species are not in ode_calibrated_results.json")
        if m["fast_species"] in wo and m["slow_species"] in wo:
            check(abs(m["fast_hours_to_baseline"]
                      - wo[m["fast_species"]]["hours_to_return_to_baseline"]) < 1e-3,
                  "animation manifest: fast-species washout time disagrees with the ODE table")
            check(abs(m["slow_hours_to_baseline"]
                      - wo[m["slow_species"]]["hours_to_return_to_baseline"]) < 1e-3,
                  "animation manifest: slow-species washout time disagrees with the ODE table")

    if anim is not None and dos is not None:
        d = anim["dosing_sawtooth"]
        mir7 = dos["daily_dosing_comparison"]["miR-7"]
        check(abs(d["peak_over_average_fast"] - mir7["peak_over_average_daily"]) < 1e-3,
              "animation manifest: fast-species peak/average disagrees with mimic_dosing_feasibility.json")
        check(abs(d["peak_over_average_slow"] - mir7["reference_peak_over_average_daily"]) < 1e-3,
              "animation manifest: reference peak/average disagrees with mimic_dosing_feasibility.json")
        check(abs(d["penalty_ratio"] - mir7["penalty_versus_reference"]) < 1e-3,
              "animation manifest: penalty ratio disagrees with mimic_dosing_feasibility.json")

    anim_stems = ("anim_ad_monomer", "anim_mimic_washout", "anim_dosing_sawtooth")
    for stem in anim_stems:
        for lang in ("en", "pt-BR"):
            fp = f"results/figures/{stem}.{lang}.gif"
            claim(os.path.isfile(fp), f"animation: {fp} is missing")

    # EN | The structure story panel composites scripts/13's own renders and fails loudly
    #      at generation time (missing provenance, missing PNG) rather than at check time,
    #      matching scripts/13's own stance that these figures are illustration, not a
    #      result, and so carry a lighter checking bar than a quantitative table. All
    #      scripts/08 confirms is that the composite was actually produced.
    # PT | O painel estrutural composto reaproveita os renders do proprio scripts/13 e
    #      falha alto na geracao (proveniencia ou PNG faltando), nao na checagem, no mesmo
    #      espirito do proprio scripts/13 de que essas figuras sao ilustracao, nao
    #      resultado, e por isso carregam uma barra de checagem mais leve que uma tabela
    #      quantitativa. Tudo que o scripts/08 confirma e que o composto foi de fato
    #      produzido.
    for lang in ("en", "pt-BR"):
        fp = f"results/figures/structure_story_panel.{lang}.png"
        claim(os.path.isfile(fp), f"structure story panel: {fp} is missing")

    # EN | The graphical abstract's every number is loaded at draw time from files
    #      scripts/08 already checks (see the animation and dosing sections above), so
    #      what remains to confirm here is only that the manifest and the two images
    #      were actually produced together.
    # PT | Todo numero do graphical abstract e carregado no momento do desenho de
    #      arquivos que o scripts/08 ja confere (ver as secoes de animacao e dosagem
    #      acima), entao o que resta confirmar aqui e so que o manifesto e as duas
    #      imagens foram de fato produzidos juntos.
    claim(os.path.isfile("results/tables/graphical_abstract_values.json"),
          "graphical abstract: results/tables/graphical_abstract_values.json is missing")
    for lang in ("en", "pt-BR"):
        fp = f"results/figures/graphical_abstract.{lang}.png"
        claim(os.path.isfile(fp), f"graphical abstract: {fp} is missing")

    # --- 9b. A prose claim about a scanned range must match the scan -------
    # EN | The methods prose once said the required BACE1 knockdown was "19-33%",
    #      copied from before d_mRNA was measured (K068) and the free-parameter grid
    #      narrowed. Nothing recomputed it when the grid changed, and it sat wrong in
    #      both language files. It is corrected here and pinned to the JSON, exactly as
    #      the README's own numbers were pinned in section 10 below: the whole point of
    #      writing the range into ode_calibrated_results.json was so a sentence quoting
    #      it could be checked instead of trusted.
    # PT | A prosa dos metodos dizia "19-33%" para a queda de BACE1 necessaria, copiada
    #      de antes de o d_mRNA ser medido (K068) e a grade de parametros livres se
    #      estreitar. Nada recalculou quando a grade mudou, e o texto ficou errado nos
    #      dois idiomas. E corrigido aqui e preso ao JSON, exatamente como os numeros do
    #      README sao presos na secao 10 abaixo: o motivo de escrever a faixa em
    #      ode_calibrated_results.json era justamente poder conferir uma frase que a cita,
    #      em vez de confiar nela.
    if ode is not None and "mimic_dose_vs_measured_knockdown" in ode:
        kappa_doc = ode["mimic_dose_vs_measured_knockdown"]["required_BACE1_knockdown"]
        # EN/PT: the Portuguese text writes the decimal with a comma, the English with a point
        anchor = {"docs/en/METHODS.md": ("knockdown needed to offset", f"{100 * kappa_doc:.1f}%"),
                  "docs/pt-BR/METODOS.md": ("queda de BACE1 necess\u00e1ria",
                                            f"{100 * kappa_doc:.1f}%".replace(".", ","))}
        for name, (needle, expect) in anchor.items():
            try:
                text = open(name, encoding="utf-8").read()
            except FileNotFoundError:
                continue
            claim(needle in text, f"{name}: the sentence this check anchors to is gone; "
                                  "update the anchor before trusting its pass")
            claim(expect in text,
                  f"{name}: does not quote the current required BACE1 knockdown {expect} "
                  f"computed from {ODE_RESULTS}")

    # --- 9b2. The prose about scripts/36 must quote what scripts/36 computed ---
    # EN | Same rule as 9b: a number written into the methods about the repeated-dosing
    #      analysis is looked up in its JSON, in the decimal form each language uses, so
    #      a rerun that moves a value cannot leave the sentence behind.
    # PT | Mesma regra do 9b: um numero escrito nos metodos sobre a analise de dosagem
    #      repetida e conferido no JSON dela, na forma decimal de cada idioma, para que
    #      uma nova rodada que mova um valor nao deixe a frase para tras.
    if tgt is not None:
        rs_, nuc_ = tgt["repression_strength"], tgt["nucleation_scaling"]
        quotes = {
            "docs/en/METHODS.md": (
                f"\u03c1 = {rs_['PD']['rho_floor_K065']:.2f}",
                f"between {rs_['PD']['rho_scan_low']:.2f} and {round(rs_['PD']['rho_scan_high'])} in the PD arm",
                f"between {rs_['AD']['rho_scan_low']:.1f} and {round(rs_['AD']['rho_scan_high'])} in the AD arm",
                f"nucleation rate ratio of {nuc_['primary_nucleation_rate_ratio_AD_over_control']:.3f}",
                f"{100 * nuc_['primary_nucleation_rate_reduction_after_offset']:.1f}% fall in the nucleation rate"),
            "docs/pt-BR/METODOS.md": (
                f"\u03c1 = {rs_['PD']['rho_floor_K065']:.2f}".replace(".", ","),
                f"entre {rs_['PD']['rho_scan_low']:.2f} e {round(rs_['PD']['rho_scan_high'])} no bra\u00e7o da DP".replace(".", ","),
                f"entre {rs_['AD']['rho_scan_low']:.1f} e {round(rs_['AD']['rho_scan_high'])} no bra\u00e7o da DA".replace(".", ","),
                f"raz\u00e3o de taxa de nuclea\u00e7\u00e3o de {nuc_['primary_nucleation_rate_ratio_AD_over_control']:.3f}".replace(".", ","),
                f"{100 * nuc_['primary_nucleation_rate_reduction_after_offset']:.1f}% na taxa de nuclea\u00e7\u00e3o".replace(".", ",")),
        }
        for name, needles in quotes.items():
            try:
                text = open(name, encoding="utf-8").read()
            except FileNotFoundError:
                continue
            for needle in needles:
                claim(needle in text,
                      f"{name}: does not quote '{needle}' from {TARGET_JSON}")

    # --- 9b3. The design-space grid and the molecular renders keep to their sources ---
    # EN | The design-space maps of Article 2 (scripts/40) are drawn from a grid of 8379 solved
    #      points and the structure panels (scripts/39) from four deposited entries. The grid
    #      has to be complete and free of duplicates, to reproduce scripts/36 where the two
    #      overlap (daily dosing), to be monotone in the dosing interval, and to count the
    #      points that no dose reaches as many times as the table lists them; every render has
    #      to exist, point at an entry that has provenance, and, where a layer was constructed
    #      for the picture, to have been built from an exact superposition of two deposited
    #      chains.
    # PT | Os mapas de espaco de projeto do Artigo 2 (scripts/40) saem de uma grade de 8379
    #      pontos resolvidos e os paineis de estrutura (scripts/39) de quatro entradas
    #      depositadas. A grade tem de estar completa e sem duplicatas, reproduzir o
    #      scripts/36 onde os dois se sobrepoem (dose diaria), ser monotona no intervalo e
    #      contar os pontos que nenhuma dose alcanca tantas vezes quanto a tabela os lista;
    #      todo render tem de existir, apontar para uma entrada com proveniencia e, onde uma
    #      camada foi construida para a figura, ter sido feita de uma superposicao exata de
    #      duas cadeias depositadas.
    try:
        ds_ = json.load(open("results/tables/mimic_design_space.json", encoding="utf-8"))
        ds_rows_ = list(csv.DictReader(open("results/tables/mimic_design_space.csv", encoding="utf-8")))
        n_expect_ = len(ds_["series"]) * len(ds_["rho_grid"]) * len(ds_["interval_grid_h"])
        check(len(ds_rows_) == n_expect_, f"design space: {len(ds_rows_)} rows in the table, grid says {n_expect_}")
        keys_ = {(r["axis"], r["species"], r["target_knockdown"], r["rho"], r["interval_h"]) for r in ds_rows_}
        check(len(keys_) == len(ds_rows_), "design space: a grid point is listed twice")
        v36_ = ds_["verification"]["against_sweep_36_at_24h"]
        check(v36_["points"] > 0 and v36_["max_relative_error"] <= v36_["tolerance"],
              f"design space: differs from scripts/36 at daily dosing by {v36_['max_relative_error']:.2e} (tolerance {v36_['tolerance']:.0e})")
        mono_ = ds_["verification"]["monotone_in_interval"]
        check(mono_["max_relative_fall"] <= 1e-9 and mono_["reached_after_unreached"] == 0,
              "design space: the peak needed falls with the dosing interval, or a target is reached after one that was not")
        n_unreached_ = sum(1 for r in ds_rows_ if r["reachable"] != "True")
        check(n_unreached_ == ds_["unreached_points"],
              f"design space: table lists {n_unreached_} unreached points, JSON says {ds_['unreached_points']}")
        if tgt is not None:
            check(ds_["dose_cap_in_baseline_pools"] == tgt["dose_cap_in_baseline_pools"],
                  "design space: dose cap differs from the cap of scripts/36")
        # the classes of Fig. 2: the saved table is what scripts/_provenance.py computes from the kinetic table now
        sys.path.insert(0, "scripts")
        import _provenance as pc_
        kin_rows_ = list(csv.DictReader(open(KINETICS, encoding="utf-8")))
        saved_ = list(csv.DictReader(open(pc_.OUT, encoding="utf-8")))
        check(saved_ == [{k_: str(v_) for k_, v_ in e_.items()} for e_ in pc_.table(kin_rows_)],
              f"{pc_.OUT} differs from the classification that scripts/_provenance.py gives for the kinetic table; run scripts/38")
        check(len(saved_) == len(kin_rows_) and all(r_["system"] != "none" for r_ in saved_ if r_["kind"] != "declared_gap"),
              "provenance classes: a row that is not the declared gap has no system, or a row is missing")
        rn_ = json.load(open("results/tables/structure_renders_article2.json", encoding="utf-8"))
        pv_ = json.load(open("results/tables/structure_figure_provenance.json", encoding="utf-8"))
        for name_, info_ in rn_["renders"].items():
            check(os.path.exists(f"results/figures/structures/{name_}.png"), f"structure render {name_}.png is missing")
            check(info_["pdb_id"] in pv_ and pv_[info_["pdb_id"]].get("method") and pv_[info_["pdb_id"]].get("resolution_angstrom"),
                  f"structure render {name_}: entry {info_['pdb_id']} has no method and resolution in structure_figure_provenance.json")
            for lt_ in info_.get("layer_transformation", []):
                check(lt_["rms_ca_angstrom"] < 0.05 and lt_["ca_pairs"] >= 40,
                      f"structure render {name_}: the constructed layer rests on a fit of {lt_['rms_ca_angstrom']} A over {lt_['ca_pairs']} pairs")
    except FileNotFoundError as exc_:
        claim(False, f"{exc_.filename} missing; run scripts/39_structure_renders_article2.py and scripts/40_mimic_design_space.py")

    # --- 9d. The ranking of the mimic series rests on the design-space summary ---
    # EN | scripts/46 ranks the series by the longest dosing interval within 3, 10 or 30 times
    #      the baseline pool. The ranking is recomputed here from the same JSON by a separate
    #      loop, the reference series is kept out of the ranks, ties share a rank, and the
    #      featured series of each axis is checked against the file's own statement that first
    #      place does not change with the threshold, the anchor or the target.
    # PT | O scripts/46 ordena as series pelo maior intervalo de dose dentro de 3, 10 ou 30
    #      vezes o reservatorio basal. A ordenacao e recalculada aqui a partir do mesmo JSON por
    #      um laco separado, a serie de referencia fica fora das posicoes, empates dividem a
    #      posicao, e a serie em destaque de cada eixo e conferida com a afirmacao do proprio
    #      arquivo de que o primeiro lugar nao muda com o limiar, a ancora ou o alvo.
    try:
        rank_json = json.load(open("results/tables/mimic_series_ranking.json", encoding="utf-8"))
        rank_rows = list(csv.DictReader(open("results/tables/mimic_series_ranking.csv", encoding="utf-8")))
        space = json.load(open("results/tables/mimic_design_space.json", encoding="utf-8"))
        n_keys = len(space["interval_within_threshold"])
        check(len(rank_rows) == n_keys * 3, f"series ranking: {len(rank_rows)} rows for {n_keys} series-anchor pairs x 3 thresholds")
        by_group = {}
        for r_ in rank_rows:
            by_group.setdefault((r_["axis"], r_["target_knockdown"], r_["anchor"], r_["threshold_x"]), []).append(r_)
        tops_by_axis = {}
        for k_, g_ in by_group.items():
            studied_ = [r_ for r_ in g_ if r_["reference"] == "False"]
            ref_ = [r_ for r_ in g_ if r_["reference"] == "True"]
            check(len(ref_) == 1 and ref_[0]["rank_among_studied"] == "", f"series ranking: the reference series is ranked in {k_}")
            best_ = max(float(r_["interval_h"]) for r_ in studied_)
            for r_ in studied_:
                want_ = 1 + sum(1 for o_ in studied_ if float(o_["interval_h"]) > float(r_["interval_h"]))
                check(int(r_["rank_among_studied"]) == want_, f"series ranking: {k_} {r_['species']} has rank {r_['rank_among_studied']}, expected {want_}")
            winners_ = sorted(r_["species"] for r_ in studied_ if float(r_["interval_h"]) == best_)
            if len(winners_) == 1:
                tops_by_axis.setdefault(k_[0], set()).add(winners_[0])
        for axis_, f_ in rank_json["featured_by_axis"].items():
            check(sorted(tops_by_axis.get(axis_, set())) == f_["first_place_in_every_combination"],
                  f"series ranking: first places of {axis_} are {sorted(tops_by_axis.get(axis_, set()))}, the file says {f_['first_place_in_every_combination']}")
            check(f_["first_place_changes"] == (len(tops_by_axis.get(axis_, set())) > 1), f"series ranking: {axis_} first_place_changes is wrong")
            # EN/PT: the featured series is a measured species and, in the primary reading, has the longest half-life among them
            half_ = {s_["species"]: s_["half_life_h"] for s_ in space["series"] if s_["axis"] == axis_ and s_["species"] != "median miRNA"}
            check(f_["species"] == [max(half_, key=half_.get)], f"series ranking: the featured series of {axis_} is {f_['species']}, the most stable measured species is {max(half_, key=half_.get)}")
        check(rank_json["reference_at_least_as_long_as_first_place"]["of"] == len(by_group) - sum(1 for k_ in by_group if not any(r_["reference"] == "True" for r_ in by_group[k_])),
              "series ranking: the count of reference comparisons does not match the groups")
        # EN/PT: the combinations with a single first place and the rounding of the intervals quoted in the text | combinacoes com um so primeiro lugar e o arredondamento dos intervalos citados no texto
        single_ = sum(1 for k_, g_ in by_group.items()
                      if sum(1 for r_ in g_ if r_["reference"] == "False" and r_["rank_among_studied"] == "1") == 1)
        nums_ = json.load(open("results/tables/manuscript2_numbers.json", encoding="utf-8"))
        nums_ = nums_.get("values", nums_)
        check(int(nums_["rank.combos_single"]) == single_ and int(nums_["rank.combos"]) - int(nums_["rank.ties"]) == single_,
              f"series ranking: {single_} combinations have a single first place; the numbers file says {nums_['rank.combos_single']}")
        def near_(x_):
            return str(int(Decimal(repr(round(x_, 8))).quantize(Decimal(1), rounding=ROUND_HALF_UP))) if x_ >= 10 else str(Decimal(repr(round(x_, 8))).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
        for tok_, ax_, sp_, kd_ in (("ds.AD.c.ad.i10_rng", "AD", "miR-29c", "0.2921"), ("ds.AD.b.ad.i10_rng", "AD", "miR-29b", "0.2921"),
                                    ("ds.PD.m7ub.a.i10_rng", "PD", "miR-7 upper bound", "0.3"), ("ds.PD.m7.a.i10_rng", "PD", "miR-7", "0.3")):
            v_ = [float(r_["interval_h"]) for r_ in rank_rows if r_["axis"] == ax_ and r_["species"] == sp_ and r_["threshold_x"] == "10"
                  and abs(float(r_["target_knockdown"]) - float(kd_)) < 1e-3 and r_["anchor"] != "K065 floor"]
            want_ = near_(min(v_)) if near_(min(v_)) == near_(max(v_)) else f"{near_(min(v_))} to {near_(max(v_))}"
            check(nums_[tok_] == want_, f"series ranking: the text quotes {nums_[tok_]} h for {sp_} ({ax_}) at 10 times; the table gives {want_} h, rounded to the nearest")
    except FileNotFoundError:
        claim(False, "results/tables/mimic_series_ranking.json missing; run scripts/46 after scripts/40")

    # --- 9e. The miRNA-3'UTR pairing rests on archived sequences and is recomputed here ---
    # EN | scripts/49 scans the 3'UTR of BACE1 (NM_012104) and SNCA (NM_000345) for seed-match
    #      sites of five microRNAs. Here the archived sequences are checked against the SHA-256
    #      in provenance.json, the scan is repeated by a separate loop with regular expressions,
    #      and the sites, their types and the pair counts must equal the table.
    # PT | O scripts/49 varre o 3'UTR de BACE1 e SNCA atras de sitios de pareamento da semente
    #      de cinco microRNAs. Aqui as sequencias arquivadas sao conferidas contra o SHA-256 de
    #      provenance.json, a varredura e refeita por outro laco com expressoes regulares, e
    #      os sitios, seus tipos e as contagens de pares precisam ser iguais a tabela.
    try:
        import hashlib as hl_
        prov_ = json.load(open("data/raw/sequences_2026/provenance.json", encoding="utf-8"))
        fa_ = {}
        nm_ = None
        for ln_ in open("data/extracted/utr3_sequences.fa", encoding="utf-8"):
            ln_ = ln_.strip()
            if ln_.startswith(">"):
                nm_ = ln_[1:].split("|")[0]
                fa_[nm_] = ""
            elif nm_:
                fa_[nm_] += ln_
        for g_ in ("BACE1", "SNCA"):
            check(hl_.sha256(fa_[g_].encode()).hexdigest() == prov_["sha256"][g_], f"3'UTR of {g_}: the archived sequence no longer matches provenance.json")
            check(set(fa_[g_]) <= set("ACGU"), f"3'UTR of {g_}: a character other than A, C, G, U")
        check(len(fa_["BACE1"]) == 3880 and len(fa_["SNCA"]) == 2529, "3'UTR lengths changed (BACE1 3,880 nt, SNCA 2,529 nt)")
        mat_ = {}
        for ln_ in open("data/raw/sequences_2026/mirbase_mature_selected.fa", encoding="utf-8"):
            ln_ = ln_.strip()
            if ln_.startswith(">"):
                nm_ = ln_[1:].split()[0]
                mat_[nm_] = ""
            else:
                mat_[nm_] += ln_
        comp_ = str.maketrans("ACGU", "UGCA")
        rc_ = lambda x_: x_.translate(comp_)[::-1]
        check(len({mat_[m_][1:8] for m_ in ("hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p")}) == 1, "the three miR-29 paralogues no longer share the seed (nt 2-8)")
        tab_ = list(csv.DictReader(open("results/tables/mirna_site_pairing.csv", encoding="utf-8")))
        import re as re_
        want_ = 0
        for g_ in ("BACE1", "SNCA"):
            for m_, seq_ in mat_.items():
                s78_, s27_ = rc_(seq_[1:8]), rc_(seq_[1:7])
                u_ = fa_[g_]
                pos78_ = [m.start() for m in re_.finditer(f"(?={s78_})", u_)]
                pos27_ = [m.start() for m in re_.finditer(f"(?={s27_}A)", u_) if m.start() >= 1 and u_[m.start() - 1:m.start() + 6] != s78_]
                got_ = sorted((int(r_["utr_start"]), r_["site_type"]) for r_ in tab_ if r_["gene"] == g_ and r_["mirna"] == m_)
                exp_ = sorted([(p_ + 1, "8mer" if u_[p_ + 7:p_ + 8] == "A" else "7mer-m8") for p_ in pos78_] + [(p_ + 1, "7mer-A1") for p_ in pos27_])
                check(got_ == exp_, f"miRNA sites of {m_} in the 3'UTR of {g_}: the table has {got_}, a separate scan finds {exp_}")
                want_ += len(exp_)
        check(len(tab_) == want_, f"mirna_site_pairing.csv has {len(tab_)} sites, a separate scan finds {want_}")
        for r_ in tab_:      # the pair marks recounted from the two sequences
            m_ = mat_[r_["mirna"]]
            a1_ = int(r_["a1_position"]) - 1
            pairs_ = 0
            for k_ in range(1, len(m_) + 1):
                p_ = a1_ - (k_ - 1)
                if 0 <= p_ < len(fa_[r_["gene"]]) and (m_[k_ - 1], fa_[r_["gene"]][p_]) in {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}:
                    pairs_ += 1
            check(pairs_ == int(r_["pairs_total"]), f"{r_['mirna']} at {r_['gene']} {r_['utr_start']}: {pairs_} pairs recounted, the table says {r_['pairs_total']}")
    except FileNotFoundError as exc_:
        claim(False, f"{exc_.filename} missing; run scripts/49_mirna_site_pairing.py")

    # --- 9b. Main-text figures stay inside the journal's artwork limits ---
    # EN | scripts/_journal_fit.py records the size, lettering and line width of every
    #      main-text figure it writes. A figure regenerated outside the limits (174 x 234
    #      mm, no lettering under 6.5 pt, no line under 0.3 pt, Arial-compatible font, an
    #      EPS twin) would otherwise reach the submission unnoticed.
    # PT | scripts/_journal_fit.py registra tamanho, texto e espessura de linha de toda
    #      figura do texto principal. Figura regenerada fora dos limites (174 x 234 mm,
    #      sem texto abaixo de 6,5 pt, sem linha abaixo de 0,3 pt, fonte compativel com
    #      Arial, par EPS) chegaria ao envio sem ser notada.
    try:
        spec_rows = {r["file"]: r for r in csv.DictReader(open("results/tables/journal_figure_spec.csv", encoding="utf-8"))}
        for stem in (# scripts/38_mechanistic_figures.py
                     "mech_model_schematic", "mech_provenance", "mech_alzheimer", "mech_parkinson", "mech_turnover",
                     "mech_target_response", "mech_design_space", "mech_series_ranking", "mech_mechanism_algorithm", "mech_mirna_site_pairing"):
            for lang_ in ("en", "pt-BR"):
                f_ = f"results/figures/{stem}.{lang_}.png"
                r_ = spec_rows.get(f_)
                check(r_ is not None, f"journal artwork: {f_} has no row in journal_figure_spec.csv")
                if r_ is None:
                    continue
                check(float(r_["width_mm"]) <= 174.05 and float(r_["height_mm"]) <= 234.05,
                      f"journal artwork: {f_} is {r_['width_mm']} x {r_['height_mm']} mm, over 174 x 234 mm")
                check(float(r_["min_text_pt"]) >= 6.5 - 1e-6, f"journal artwork: {f_} has lettering of {r_['min_text_pt']} pt")
                check(float(r_["min_line_pt"]) >= 0.3, f"journal artwork: {f_} has a line of {r_['min_line_pt']} pt")
                check(all(any(k_ in fnt for k_ in ("Arial", "Helvetica", "Liberation Sans")) for fnt in r_["fonts"].split("; ")),
                      f"journal artwork: {f_} uses {r_['fonts']}")
                check(os.path.exists(f_[:-4] + ".eps"), f"journal artwork: {f_[:-4]}.eps is missing")
    except FileNotFoundError:
        claim(False, "results/tables/journal_figure_spec.csv missing; run the main-text figure scripts (14, 21, 27, 28, 32, 33, and 38 after 12, 13, 17, 36, 39 and 40)")

    # --- Lettering of the main-text figures reaches 4.5:1 and the renders keep 600 dpi in the vector files ---
    # EN | scripts/45 measures every text span of a figure PDF against the colour under it, the
    #      thinnest line, the page, and the effective resolution of every placed render. The figures
    #      were found with renders at 100 dpi in the PDF and EPS (600 dpi in the PNG), because the
    #      vector backends resample an image to the figure dpi; the check keeps that from returning.
    #      It needs PyMuPDF, an optional dependency, and says so when it is missing.
    # PT | O scripts/45 mede cada trecho de texto de um PDF de figura contra a cor sob ele, a linha
    #      mais fina, a pagina e a resolucao efetiva de cada render colocado. As figuras tinham
    #      renders a 100 dpi no PDF e no EPS (600 dpi no PNG), porque os backends vetoriais
    #      reamostram uma imagem para o dpi da figura; a conferencia impede que isso volte.
    #      Precisa do PyMuPDF, dependencia opcional, e avisa quando falta.
    try:
        import importlib.util
        spec_ = importlib.util.spec_from_file_location("art45", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                              "45_check_figure_artwork.py"))
        art45 = importlib.util.module_from_spec(spec_)
        spec_.loader.exec_module(art45)
    except ImportError:
        art45 = None
        print("EN/PT: PyMuPDF is not installed; the contrast and image resolution of the figures were not measured "
              "(pip install pymupdf) | o contraste e a resolucao de imagem das figuras nao foram medidos")
    if art45 is not None:
        for stem_ in ("mech_model_schematic", "mech_provenance", "mech_alzheimer", "mech_parkinson", "mech_turnover",
                      "mech_target_response", "mech_design_space", "mech_series_ranking", "mech_mechanism_algorithm", "mech_mirna_site_pairing"):
            for lang_ in ("en", "pt-BR"):
                f_ = f"results/figures/{stem_}.{lang_}.pdf"
                if not os.path.exists(f_):
                    claim(False, f"{f_} is missing; run the figure scripts")
                    continue
                m_ = art45.measure(f_)
                check(m_["n_bad"] == 0, f"{f_}: {m_['n_bad']} text spans below 4.5:1 (worst {m_['worst']:.2f}:1)")
                check(m_["w_mm"] <= 174.05 and m_["h_mm"] <= 234.05, f"{f_}: page {m_['w_mm']:.1f} x {m_['h_mm']:.1f} mm")
                check(m_["thinnest"] != m_["thinnest"] or m_["thinnest"] >= 0.3, f"{f_}: a line of {m_['thinnest']:.2f} pt")
                check(all(im_[4] >= 600 - 0.5 for im_ in m_["images"]), f"{f_}: an image below 600 dpi in the vector file")

    # --- README states the number of checks that ran ---------------------------------
    for name, count_claim in (("README.md", "{n} checks currently pass."),
                              ("README.pt-BR.md", "Atualmente {n} verificacoes passam.")):
        try:
            text = open(name, encoding="utf-8").read()
        except FileNotFoundError:
            continue
        stated = count_claim.format(n=checks)
        claim(stated in text or stated.replace("verificacoes", "verificações") in text,
              f"{name}: the README claims a different number of checks than the {checks} that just ran")

    print(f"Checks run / verificacoes: {checks}")
    if fails:
        print(f"\nFAILURES / FALHAS: {len(fails)}")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("\nEN | Parameters, model results and figures are mutually consistent.")
    print("PT | Parametros, resultados dos modelos e figuras sao consistentes entre si.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
