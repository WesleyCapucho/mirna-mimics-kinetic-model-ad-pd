#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Recalibrated ODE models of the miR-29/BACE1/Abeta and miR-7/SNCA/alpha-synuclein
     axes, driven by the measured constants in data/extracted/kinetic_parameters.csv.
PT | Modelos EDO recalibrados dos eixos miR-29/BACE1/Abeta e miR-7/SNCA/alfa-sinucleina,
     guiados pelas constantes medidas em data/extracted/kinetic_parameters.csv.

EN | What changed, and why. The first models of this work (an earlier script of the companion
     repository, 10.5281/zenodo.23438256) came from the source monograph and
     used illustrative parameters: one lumped "miRNA" variable, first-order aggregation,
     and a disease effect entering through production. Every one of those choices is
     contradicted by a measured value now in the parameter table:

       1. There is no single miR-29. In one assay miR-29b decays with a 7 h half-life,
          miR-29c 10.6 h, and miR-29a does not decay measurably in 12 h (K017-K019).
          The three paralogues are separate state variables here.
       2. miR-7 is not a stable miRNA. Its half-life is 1.7 h, the shortest in a
          genome-wide survey whose median is 34 h (K016, K020), because the lncRNA
          Cyrano drives target-directed degradation. A single mimic dose is therefore
          back within 10% of baseline in 11-32 h (K020, K021), where the median
          miRNA would take about 9 days - a gap a generic decay rate hides.
       3. alpha-synuclein fibrils do not grow first-order in monomer. Elongation
          saturates with a half-maximal concentration near 46-50 uM (K009, K010).
       4. At neutral pH, primary and secondary nucleation of alpha-synuclein are
          undetectable (K011). The autocatalytic term must be off unless the
          compartment is acidic. The old model had no such gate.
       5. The Abeta42 elongation constant is NOT in the parameter table: it could not
          be read from any retrievable sentence of its primary source. It is therefore
          a free parameter, and the Abeta aggregation block is illustrative. Borrowing
          the alpha-synuclein value instead would be wrong by roughly two orders of
          magnitude against the measured clearance rate.
       6. BACE1 mRNA is not one species either. It has five 3'UTR isoforms whose
          half-lives run from 2.8 to 23.9 h (K067), so the pool's decay constant is the
          abundance-weighted mean of the RATE constants, not of the half-lives (K068).
          SNCA has a single isoform at 6.53 h (K066). Both were a declared gap until the
          Tushev supplementary table was read, and both are now measured, which is why
          mRNA decay is no longer a free parameter.
       7. In humans, Abeta42 production does not differ between AD and control
          (6.7 vs 6.6 %/h, p = 0.96); clearance does (7.6 vs 5.3 %/h, p = 0.03)
          (K001-K004). The disease lever is clearance. The old model moved production.

PT | O que mudou, e por que. Os primeiros modelos deste trabalho (um script anterior do repositorio
     companheiro, 10.5281/zenodo.23438256) vieram da monografia e usavam
     parametros ilustrativos: uma variavel "miRNA" agregada, agregacao de primeira
     ordem e efeito de doenca entrando pela producao. Cada uma dessas escolhas e
     contrariada por um valor medido que agora esta na tabela de parametros:

       1. Nao existe um unico miR-29. Num mesmo ensaio o miR-29b decai com meia-vida
          de 7 h, o miR-29c 10,6 h e o miR-29a nao decai de forma mensuravel em 12 h
          (K017-K019). Os tres paralogos sao variaveis de estado separadas aqui.
       2. O miR-7 nao e um miRNA estavel. Sua meia-vida e 1,7 h, a mais curta de um
          levantamento genomico cuja mediana e 34 h (K016, K020), porque o lncRNA
          Cyrano promove degradacao dirigida pelo alvo. Uma dose unica de mimetico volta,
          portanto, a menos de 10% do basal em 11-32 h (K020, K021), quando o miRNA
          mediano levaria cerca de 9 dias - diferenca que uma taxa generica esconde.
       3. Fibrilas de alfa-sinucleina nao crescem em primeira ordem no monomero. A
          elongacao satura, com concentracao de meia-saturacao perto de 46-50 uM
          (K009, K010).
       4. Em pH neutro, nucleacao primaria e secundaria da alfa-sinucleina sao
          indetectaveis (K011). O termo autocatalitico precisa ficar desligado a menos
          que o compartimento seja acido. O modelo antigo nao tinha essa trava.
       5. O mRNA de BACE1 tambem nao e uma especie so. Ele tem cinco isoformas de 3'UTR
          com meias-vidas de 2,8 a 23,9 h (K067), entao a constante de decaimento do pool
          e a media ponderada das CONSTANTES DE VELOCIDADE, nao das meias-vidas (K068). O
          SNCA tem uma isoforma so, de 6,53 h (K066). Os dois eram lacuna declarada ate a
          tabela suplementar do Tushev ser lida, e os dois agora sao medidos, e por isso o
          decaimento de mRNA deixou de ser parametro livre.
       6. Em humanos, a producao de Abeta42 nao difere entre AD e controle
          (6,7 vs 6,6 %/h, p = 0,96); a depuracao difere (7,6 vs 5,3 %/h, p = 0,03)
          (K001-K004). A alavanca da doenca e a depuracao. O modelo antigo mexia na
          producao.

EN | What this model is and is not. Where a constant was measured it is used and cited
     by param_id. Where it was not - the translation rate and the strength of miRNA
     repression - it is a FREE parameter,
     named as such, and the conclusions are reported across a scan of its plausible
     range rather than at one convenient value. No free parameter was tuned to make a
     result come out. The output is therefore semi-quantitative: it supports statements
     about direction, ordering and order of magnitude, not about absolute concentrations.
PT | O que este modelo e e o que nao e. Onde a constante foi medida, ela e usada e
     citada pelo param_id. Onde nao foi - taxa de traducao e forca da repressao por
     miRNA - ela e parametro LIVRE, nomeado como
     tal, e as conclusoes sao reportadas ao longo de uma varredura da faixa plausivel,
     nao num valor conveniente. Nenhum parametro livre foi ajustado para produzir um
     resultado. A saida e, portanto, semiquantitativa: sustenta afirmacoes sobre
     direcao, ordenamento e ordem de grandeza, nao sobre concentracoes absolutas.

    python scripts/12_ode_models_calibrated.py
"""

import csv
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# EN/PT: imported as tr because t is the time vector throughout this script
from _bilingual import LANGS, t as tr, fig_path

KINETICS = "data/extracted/kinetic_parameters.csv"
TAB_DIR = "results/tables"
FIG_DIR = "results/figures"
LN2 = math.log(2.0)


# --------------------------------------------------------------------------
# Measured parameters | Parametros medidos
# --------------------------------------------------------------------------
def load_measured():
    """
    EN | Read the measured constants by param_id. Nothing in this script invents a
         number: if an id is missing the script stops rather than falling back to a
         default, because a silent default is how an illustrative parameter becomes a
         quoted result.
    PT | Le as constantes medidas por param_id. Nada neste script inventa numero: se um
         id faltar, o script para em vez de recorrer a um padrao, porque um padrao
         silencioso e justamente como um parametro ilustrativo vira resultado citado.
    """
    if not os.path.exists(KINETICS):
        sys.exit(f"EN/PT: missing {KINETICS}")
    table = {r["param_id"]: r for r in csv.DictReader(open(KINETICS, encoding="utf-8"))}

    def val(pid):
        r = table.get(pid)
        if r is None:
            sys.exit(f"EN/PT: parameter {pid} not found in {KINETICS}")
        if r["kind"] != "numeric":
            sys.exit(f"EN/PT: parameter {pid} is {r['kind']}, not a numeric value")
        return float(r["value_si"])

    def derived(pid):
        """
        EN | Load a value we computed ourselves from a primary table. Kept separate from
             val() on purpose: a derived number must be asked for by name, never picked
             up by accident, and its row has to say how it was derived. scripts/08
             recomputes every derived row from its sources.
        PT | Carrega um valor que nos mesmos calculamos a partir de uma tabela primaria.
             Separado de val() de proposito: um numero derivado tem de ser pedido pelo
             nome, nunca apanhado por acaso, e a linha dele precisa dizer como foi
             derivado. O scripts/08 recalcula cada linha derivada a partir das fontes.
        """
        r = table.get(pid)
        if r is None:
            sys.exit(f"EN/PT: parameter {pid} not found in {KINETICS}")
        if r["kind"] != "derived":
            sys.exit(f"EN/PT: parameter {pid} is {r['kind']}, not a derived value")
        if "DERIVED" not in r["note"] and "derived" not in r["note"]:
            sys.exit(f"EN/PT: parameter {pid} does not document its derivation")
        return float(r["value_si"])

    m = {
        # Abeta42 handling, human in vivo (Mawuenyega 2010)
        "k_prod_Ab_ctrl": val("K001"),      # 1/h
        "k_prod_Ab_ad":   val("K002"),      # 1/h
        "k_clear_Ab_ctrl": val("K003"),     # 1/h
        "k_clear_Ab_ad":   val("K004"),     # 1/h
        # alpha-synuclein fibril elongation (Buell 2014); M^-1 s^-1 -> uM^-1 h^-1
        "kplus_high": val("K007") * 3600.0 / 1e6,
        "kplus_low":  val("K008") * 3600.0 / 1e6,
        "m_half_uM":  val("K009") * 1e6,    # M -> uM
        "m_half_uM_b": val("K010") * 1e6,
        # Abeta42 aggregation reaction orders (Cohen 2013).
        # EN | K013 is the PRIMARY order and K012 the SECONDARY one. The two were read the
        #      other way round until scripts/36 asked for them by meaning; both rows are 2,
        #      so no earlier output changes, but the names now match the rows.
        # PT | K013 e a ordem PRIMARIA e K012 a SECUNDARIA. Os dois eram lidos ao
        #      contrario ate o scripts/36 pedi-los pelo significado; as duas linhas valem
        #      2, entao nenhuma saida anterior muda, mas os nomes agora batem com as linhas.
        "n_c": val("K013"),
        "n_2": val("K012"),
        "M_star_uM": val("K014") * 1e6,     # M -> uM
        # miRNA decay constants, 1/h
        "d_miR_median": val("K016"),
        "d_miR29b": val("K017"),
        "d_miR29c": val("K018"),
        "d_miR7":   val("K020"),
        "d_miR7_slow": val("K021"),
        # mRNA decay of the two target genes, measured in neurons (Tushev 2018)
        "d_mRNA_SNCA": val("K066"),
        # absolute abundances in the presynaptic bouton (Wilhelm 2014)
        "C_aSyn_uM": derived("K054"),
        "d_mRNA_BACE1": derived("K068"),           # alpha only, derived from K052 and K053
        "N_BACE1": val("K055"),
        "N_APP": val("K056"),
        # measured size of miR-7 repression of alpha-synuclein (Doxakis 2010)
        "kd_SNCA_miR7_construct": val("K058"),
        "kd_SNCA_miR153_construct": val("K059"),
        "kd_SNCA_pair_construct": val("K060"),
        "kd_SNCA_neurons": val("K063"),
        # measured size of miR-29 repression of BACE1 (Hebert 2008)
        "BACE1_knockdown_miR29": val("K047"),
        # miR-7 neuronal steady state (Kleaveland 2018)
        "miR7_copies": val("K022"),
        "miR7_fold_noCyrano": val("K023"),
        # protein decay, 1/h
        "d_BACE1_hek": val("K024"),
        "d_BACE1_neuron": val("K025"),
        "d_aSyn": val("K026"),
        "d_proteome_brain": val("K029"),
        # human ageing observations
        "aSyn_protein_aging": val("K027"),  # +1.0 = doubling
        "SNCA_mRNA_aging": val("K028"),     # -0.6
        # Abeta42 aggregate load in human brain, Cohen 2013 Table S2; M
        "Ab_load_ctrl_median": val("K032"),
        "Ab_load_AD_median":   val("K035"),
        "Ab_load_ctrl_lq":     val("K031"),
        "Ab_load_AD_uq":       val("K036"),
        # EN/PT: K040/K041 (genome-wide medians) are deliberately NOT loaded here. They
        #        are kind=derived - computed by this project from a supplementary table,
        #        not read in a sentence - and they were only ever context for the mRNA
        #        decay rate, which is now measured per gene (K066, K068). The val() guard
        #        refuses non-numeric kinds, which is what caught an earlier attempt to
        #        feed one straight into the model.
    }
    # EN/PT: miR-29a does not decay measurably within the assay (K019, qualitative).
    #        Represented by the slowest rate the data can support, the whole-population
    #        median, and flagged in the output rather than invented as a number.
    m["d_miR29a_proxy"] = m["d_miR_median"]
    return m, table


# --------------------------------------------------------------------------
# Free parameters | Parametros livres
# --------------------------------------------------------------------------
# EN | Every entry here is UNMEASURED in the parameter table. The scan range is what
#      the scripts explore; the midpoint is used only to draw the time courses.
#      mRNA decay used to live here. It does not any more: Tushev et al. 2018 measured
#      it in neurons for both genes (SNCA 6.53 h, K066; BACE1 pooled over its five
#      3'UTR isoforms, K068), so each arm now carries its own measured rate and the
#      scan is one parameter shorter.
# PT | Toda entrada aqui e NAO MEDIDA na tabela de parametros. A faixa e o que os
#      scripts varrem; o ponto medio serve apenas para desenhar as curvas.
#      O decaimento de mRNA morava aqui. Nao mora mais: Tushev et al. 2018 o mediram em
#      neuronios para os dois genes (SNCA 6,53 h, K066; BACE1 agrupado sobre suas cinco
#      isoformas de 3'UTR, K068), entao cada braco carrega agora sua taxa medida e a
#      varredura tem um parametro a menos.
FREE = {
    "k_repress": dict(low=0.1, high=3.0,
                      why_en="Strength of miRNA repression per unit miRNA. Not measured as a "
                              "rate anywhere in the table.",
                      why_pt="Forca da repressao por unidade de miRNA. Nao medida como taxa em "
                              "nenhum ponto da tabela."),
    "k_translate": dict(low=0.2, high=5.0,
                        why_en="Protein output per unit mRNA per hour.",
                        why_pt="Producao de proteina por unidade de mRNA por hora."),
}


# EN | Aggregation rate constants that no source in the table gives as a number.
#      Unlike FREE, they are not swept in the main scan, because no result that is
#      reported as a finding depends on them: the Abeta block has no feedback on the
#      monomer pool, and for alpha-synuclein only the DIRECTION of the pH switch is
#      sourced (K011, K042). experiment_ph_gate sweeps k2_aSyn to show that the size
#      of the switch is set by these numbers and not by data.
# PT | Constantes de agregacao que nenhuma fonte da tabela da como numero. Diferente
#      de FREE, nao entram na varredura principal, porque nenhum resultado reportado
#      como achado depende delas: o bloco do Abeta nao realimenta o pool de monomero,
#      e para a alfa-sinucleina so a DIRECAO da chave de pH tem fonte (K011, K042).
#      experiment_ph_gate varre k2_aSyn para mostrar que o tamanho da chave e fixado
#      por esses numeros e nao por dados.
ILLUSTRATIVE = {
    "kn_Ab": dict(value=1e-6, why_en="Abeta42 primary nucleation. Not separable from k+ with "
                  "the published data (K037).", why_pt="Nucleacao primaria do Abeta42. Nao "
                  "separavel de k+ com os dados publicados (K037)."),
    "k2_Ab": dict(value=1e-4, why_en="Abeta42 secondary nucleation. Not separable (K037).",
                  why_pt="Nucleacao secundaria do Abeta42. Nao separavel (K037)."),
    "kplus_Ab": dict(value=1e-3, why_en="Abeta42 elongation. Not separable (K037); the "
                     "alpha-synuclein value must not be borrowed.", why_pt="Alongamento do "
                     "Abeta42. Nao separavel (K037); o valor da alfa-sinucleina nao pode ser "
                     "emprestado."),
    "kn_aSyn": dict(value=1e-8, why_en="alpha-synuclein primary nucleation at acidic pH. Only "
                    "qualitative in the source (K042).", why_pt="Nucleacao primaria da "
                    "alfa-sinucleina em pH acido. So qualitativa na fonte (K042)."),
    "k2_aSyn": dict(value=1e-6, why_en="alpha-synuclein secondary nucleation at acidic pH. "
                    "Only qualitative in the source (K042).", why_pt="Nucleacao secundaria da "
                    "alfa-sinucleina em pH acido. So qualitativa na fonte (K042)."),
}


def midpoint(name):
    """
    EN | The central value of a free parameter: a measured value where one exists for a
         comparable system, otherwise the geometric mean of the declared range. This is
         still not a fit - the default is read from the parameter table, not chosen to
         make a result come out - but it is better than the middle of an interval.
    PT | O valor central de um parametro livre: um valor medido quando existe para um
         sistema comparavel, senao a media geometrica da faixa declarada. Isto continua
         nao sendo ajuste - o padrao vem da tabela de parametros, nao e escolhido para
         produzir um resultado - mas e melhor que o meio de um intervalo.
    """
    f = FREE[name]
    if "default" in f:
        assert f["low"] <= f["default"] <= f["high"], f"{name}: default outside its range"
        return f["default"]
    return math.sqrt(f["low"] * f["high"])   # EN/PT: geometric mean of the scan range


# --------------------------------------------------------------------------
# AD arm | Braco AD : miR-29 paralogues -> BACE1 -> Abeta42
# --------------------------------------------------------------------------
def ad_rhs(t, y, p):
    """
    EN | State: [miR-29a, miR-29b, miR-29c, BACE1 mRNA, BACE1 protein, Abeta monomer,
         fibril number, fibril mass]. Abeta monomer is in relative units; the
         aggregation block uses the measured reaction orders and the critical fibril
         concentration, both from Cohen 2013.
    PT | Estado: [miR-29a, miR-29b, miR-29c, mRNA de BACE1, proteina BACE1, monomero de
         Abeta, numero de fibrilas, massa de fibrilas]. O monomero esta em unidades
         relativas; o bloco de agregacao usa as ordens de reacao medidas e a
         concentracao critica de fibrila, ambas de Cohen 2013.
    """
    a29a, a29b, a29c, mR, pB, A, P, M = y

    # --- miR-29 paralogues, each with its own measured decay -----------------
    d_a = a29a * p["d_miR29a"]
    d_b = a29b * p["d_miR29b"]
    d_c = a29c * p["d_miR29c"]
    dm29a = p["s29a"] - d_a
    dm29b = p["s29b"] - d_b
    dm29c = p["s29c"] - d_c

    # --- BACE1 mRNA: repressed by the summed miR-29 pool ---------------------
    repression = p["k_repress"] * (a29a + a29b + a29c)
    dmR = p["s_mRNA"] - (p["d_mRNA"] + repression) * mR

    # --- BACE1 protein -------------------------------------------------------
    dpB = p["k_translate"] * mR - p["d_BACE1"] * pB

    # --- Abeta42 monomer pool --------------------------------------------
    # EN | A is the MONOMER pool, which is exactly what Mawuenyega et al. measured:
    #      fractional production and clearance of labelled Abeta42 in human CSF. It is
    #      not a total that also contains fibril mass. An earlier version subtracted the
    #      fibril mass from it and then subtracted the elongation flux again, counting
    #      the same molecules twice and driving the AD trajectory BELOW control.
    # PT | A e o pool de MONOMERO, exatamente o que Mawuenyega et al. mediram: producao e
    #      depuracao fracionais de Abeta42 marcado no liquor humano. Nao e um total que
    #      contenha tambem massa de fibrila. Uma versao anterior subtraia a massa de
    #      fibrila dele e depois subtraia o fluxo de elongacao de novo, contando as
    #      mesmas moleculas duas vezes e levando a trajetoria de AD ABAIXO do controle.
    dA = p["k_prod_Ab"] * pB - p["k_clear_Ab"] * A

    # --- Abeta42 aggregation: reaction orders measured, rate constants NOT ----
    # EN | Cohen 2013 gives n_c, n_2 and the critical fibril concentration M_star, but
    #      its elongation constant could not be read from any retrievable sentence of
    #      the primary source, so k_plus, k_n and k_2 are FREE here. An earlier version
    #      borrowed the alpha-synuclein k_plus (K007), which is about a hundred times
    #      the Abeta clearance rate and swamped the monomer pool. Constants are not
    #      transferable between proteins, and this block is reported as illustrative.
    # PT | Cohen 2013 da n_c, n_2 e a concentracao critica de fibrila M_star, mas sua
    #      constante de elongacao nao pode ser lida em nenhuma frase recuperavel da fonte
    #      primaria, entao k_plus, k_n e k_2 sao LIVRES aqui. Uma versao anterior tomou
    #      emprestado o k_plus da alfa-sinucleina (K007), cerca de cem vezes a taxa de
    #      depuracao do Abeta, e engoliu o pool de monomero. Constantes nao sao
    #      transferiveis entre proteinas, e este bloco e reportado como ilustrativo.
    nucleation = p["k_n"] * A ** p["n_c"]
    secondary = p["k_2"] * A ** p["n_2"] * M if M > p["M_star"] else 0.0
    dP = nucleation + secondary
    dM = 2.0 * p["k_plus_Ab_FREE"] * A * P

    return [dm29a, dm29b, dm29c, dmR, dpB, dA, dP, dM]


def ad_monomer_rhs(t, y, p):
    """
    EN | The first six equations of ad_rhs: the miR-29 paralogues, BACE1 and the Abeta42
         monomer. The aggregation block has no feedback on any of them, so the monomer
         trajectory is the same with or without it. This form exists so that a question
         about the monomer is not at the mercy of the illustrative block, which does not
         conserve mass and diverges at two corners of the free-parameter scan. It calls
         ad_rhs rather than restating the equations, so there is one copy of them.
    PT | As seis primeiras equacoes de ad_rhs: os paralogos de miR-29, a BACE1 e o
         monomero de Abeta42. O bloco de agregacao nao realimenta nenhuma delas, entao a
         trajetoria do monomero e a mesma com ou sem ele. Esta forma existe para que uma
         pergunta sobre o monomero nao dependa do bloco ilustrativo, que nao conserva
         massa e diverge em dois cantos da varredura de parametros livres. Chama ad_rhs em
         vez de reescrever as equacoes, entao ha uma copia so delas.
    """
    return ad_rhs(t, list(y) + [0.0, 0.0], p)[:6]


def run_ad(measured, condition, free=None, t_end=2000.0, mimic=1.0, aggregation=True):
    """
    EN | One AD trajectory. `mimic` multiplies miR-29 synthesis. With aggregation=False
         only the monomer cascade is integrated (six states instead of eight); sol.y[5]
         is the monomer pool either way.
    PT | Uma trajetoria de AD. `mimic` multiplica a sintese de miR-29. Com
         aggregation=False so a cascata do monomero e integrada (seis estados em vez de
         oito); sol.y[5] e o pool de monomero nos dois casos.
    """
    free = free or {k: midpoint(k) for k in FREE}
    p = {
        "d_miR29a": measured["d_miR29a_proxy"],
        "d_miR29b": measured["d_miR29b"],
        "d_miR29c": measured["d_miR29c"],
        "d_mRNA": measured["d_mRNA_BACE1"],
        "k_repress": free["k_repress"],
        "k_translate": free["k_translate"],
        "d_BACE1": measured["d_BACE1_neuron"],
        "n_c": measured["n_c"],
        "n_2": measured["n_2"],
        "M_star": measured["M_star_uM"],
        # EN/PT: ILLUSTRATIVE - no measured Abeta42 constant exists (K037). They set
        #        the aggregation block only and do not feed back on the monomer pool.
        "k_n": ILLUSTRATIVE["kn_Ab"]["value"],
        "k_2": ILLUSTRATIVE["k2_Ab"]["value"],
        "k_plus_Ab_FREE": ILLUSTRATIVE["kplus_Ab"]["value"],
    }
    # EN/PT: each paralogue held at a steady state of 1.0 before the mimic
    p["s29a"] = 1.0 * p["d_miR29a"] * mimic
    p["s29b"] = 1.0 * p["d_miR29b"] * mimic
    p["s29c"] = 1.0 * p["d_miR29c"] * mimic
    p["s_mRNA"] = 1.0

    if condition == "control":
        p["k_prod_Ab"] = measured["k_prod_Ab_ctrl"]
        p["k_clear_Ab"] = measured["k_clear_Ab_ctrl"]
    elif condition == "AD":
        p["k_prod_Ab"] = measured["k_prod_Ab_ad"]
        p["k_clear_Ab"] = measured["k_clear_Ab_ad"]
    else:
        raise ValueError(condition)

    y0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-3, 0.0]
    rhs = ad_rhs
    if not aggregation:
        y0, rhs = y0[:6], ad_monomer_rhs
    sol = solve_ivp(rhs, (0.0, t_end), y0, args=(p,),
                    t_eval=np.linspace(0.0, t_end, 1200), method="LSODA", rtol=1e-8, atol=1e-10)
    return sol, p


# --------------------------------------------------------------------------
# PD arm | Braco PD : miR-7 -> SNCA -> alpha-synuclein
# --------------------------------------------------------------------------
def pd_rhs(t, y, p):
    """
    EN | State: [miR-7, SNCA mRNA, alpha-syn monomer, fibril number, fibril mass].
         Elongation saturates with monomer (Buell 2014), and nucleation is gated on
         compartment pH: at neutral pH both nucleation terms are zero, which is a
         measured absence, not an assumption.
    PT | Estado: [miR-7, mRNA de SNCA, monomero de alfa-sin, numero de fibrilas, massa].
         A elongacao satura com o monomero (Buell 2014) e a nucleacao e condicionada ao
         pH do compartimento: em pH neutro os dois termos sao zero, o que e ausencia
         medida, e nao suposicao.
    """
    miR7, mR, S, P, M = y

    dmiR7 = p["s_miR7"] - p["d_miR7"] * miR7
    repression = p["k_repress"] * miR7
    dmR = p["s_mRNA"] - (p["d_mRNA"] + repression) * mR

    free_S = max(S - M, 0.0)
    # EN/PT: Michaelis-Menten-like saturation, m_half from Buell 2014 (K009/K010)
    elong = 2.0 * p["k_plus"] * P * free_S / (p["m_half"] + free_S)
    if p["acidic"]:
        dP = p["k_n"] * free_S ** 2 + p["k_2"] * free_S ** 2 * M
    else:
        dP = 0.0     # EN/PT: K011, undetectable at neutral pH under quiescent conditions

    dS = p["k_translate"] * mR - p["d_aSyn"] * free_S - elong
    dM = elong
    return [dmiR7, dmR, dS, dP, dM]


def run_pd(measured, acidic=False, free=None, t_end=2000.0, mimic=1.0,
           d_miR7=None, seed=1e-3, k2=None):
    free = free or {k: midpoint(k) for k in FREE}
    p = {
        "d_miR7": d_miR7 if d_miR7 is not None else measured["d_miR7"],
        "d_mRNA": measured["d_mRNA_SNCA"],
        "k_repress": free["k_repress"],
        "k_translate": free["k_translate"],
        "d_aSyn": measured["d_aSyn"],
        "k_plus": measured["kplus_high"],
        "m_half": measured["m_half_uM"],
        "k_n": ILLUSTRATIVE["kn_aSyn"]["value"],
        "k_2": ILLUSTRATIVE["k2_aSyn"]["value"] if k2 is None else k2,
        "acidic": acidic,
    }
    p["s_miR7"] = 1.0 * p["d_miR7"] * mimic
    p["s_mRNA"] = 1.0
    y0 = [1.0, 1.0, 1.0, seed, 0.0]
    sol = solve_ivp(pd_rhs, (0.0, t_end), y0, args=(p,),
                    t_eval=np.linspace(0.0, t_end, 1200), method="LSODA", rtol=1e-8, atol=1e-10)
    return sol, p


# --------------------------------------------------------------------------
# Experiments | Experimentos
# --------------------------------------------------------------------------
N_PARALOGUES = 3     # EN/PT: miR-29a, -29b and -29c, each at a steady state of 1.0


def bace1_knockdown(d, k, m):
    """
    EN | Fractional fall in steady-state BACE1 (mRNA and, because protein follows mRNA
         linearly, protein) when every miR-29 paralogue is held at m times baseline. The
         mRNA level is s / (d + 3 k m), so against the baseline s / (d + 3 k) the fall is
         1 - (d + 3 k) / (d + 3 k m). Writing rho = 3 k / d, the decay that the miRNAs add to
         the mRNA at baseline in units of its intrinsic decay, this is 1 - (1 + rho) / (1 + rho m).
    PT | Queda fracional da BACE1 no estado estacionario (mRNA e, como a proteina segue o
         mRNA linearmente, proteina) quando cada paralogo de miR-29 e mantido em m vezes o
         basal. O nivel de mRNA e s / (d + 3 k m), entao frente ao basal s / (d + 3 k) a
         queda e 1 - (d + 3 k) / (d + 3 k m). Com rho = 3 k / d, o decaimento que os miRNAs
         somam ao do mRNA no basal, em unidades do decaimento intrinseco, vale 1 - (1 + rho) / (1 + rho m).
    """
    n = N_PARALOGUES
    return 1.0 - (d + n * k) / (d + n * k * m)


def mimic_for_knockdown(d, k, target):
    """
    EN | The inverse of bace1_knockdown: the mimic factor that gives a fractional fall
         `target`. As rho grows the answer tends to 1 / (1 - target), the proportional
         limit; as rho shrinks it grows without bound, because a miRNA that adds only a
         small decay to the mRNA's own has little leverage.
    PT | A inversa de bace1_knockdown: o fator de mimetico que da uma queda fracional
         `target`. Quando rho cresce a resposta tende a 1 / (1 - target), o limite
         proporcional; quando rho encolhe ela cresce sem limite, porque um miRNA que
         soma so um pequeno decaimento ao proprio do mRNA tem pouca alavanca.
    """
    n = N_PARALOGUES
    return 1.0 + target / (1.0 - target) * (d + n * k) / (n * k)


def required_bace1_knockdown(measured):
    """
    EN | How far BACE1 has to fall to bring the AD monomer pool back to the control
         pool, assuming production is proportional to BACE1 protein. At steady state the
         monomer pool is (k_prod / k_clear) times BACE1 protein, so the AD pool matches
         the control pool when BACE1 protein is lowered by 1 - 1 / R, with R the AD over
         control ratio of k_prod / k_clear. It depends on the two measured Abeta42
         constants and on nothing else.
    PT | O quanto a BACE1 precisa cair para trazer o pool de monomero da AD de volta ao
         pool do controle, supondo producao proporcional a proteina BACE1. No estado
         estacionario o pool de monomero e (k_prod / k_clear) vezes a proteina BACE1,
         entao o pool da AD iguala o do controle quando a proteina BACE1 cai 1 - 1 / R,
         com R a razao AD sobre controle de k_prod / k_clear. Depende das duas constantes
         medidas de Abeta42 e de mais nada.
    """
    r = ((measured["k_prod_Ab_ad"] / measured["k_clear_Ab_ad"])
         / (measured["k_prod_Ab_ctrl"] / measured["k_clear_Ab_ctrl"]))
    return 1.0 - 1.0 / r


def experiment_clearance_vs_production(measured):
    """
    EN | The central AD question. Human data say the disease lever is clearance, not
         production. A miR-29 mimic acts on production. How far must BACE1 fall, and how
         large a mimic does that take, to bring AD Abeta back to the control steady
         state?

         The required fall is closed-form and has no free parameter (see
         required_bace1_knockdown). The mimic that delivers it is closed-form too and
         depends on one free parameter, k_repress (see mimic_for_knockdown); it is
         reported at the geometric midpoint of the declared range and across the range.
         Each of the three is then checked by integrating the ODE. An earlier version
         searched a grid of mimic factors with a step of 0.247 and reported the first
         one that reached the control level, 1.494, which overshoots the exact
         requirement and made the knockdown it implied look like a range of 30-33%
         when it is 29.2% at every value of the free parameter.
    PT | A questao central de AD. Os dados humanos dizem que a alavanca e a depuracao,
         nao a producao. Um mimetico de miR-29 age na producao. Quanto a BACE1 precisa
         cair, e que tamanho de mimetico isso exige, para trazer o Abeta de AD de volta
         ao estado estacionario do controle?

         A queda necessaria e de forma fechada e nao tem parametro livre (ver
         required_bace1_knockdown). O mimetico que a entrega tambem e de forma fechada e
         depende de um parametro livre, k_repress (ver mimic_for_knockdown); e reportado
         na media geometrica da faixa declarada e ao longo da faixa. Cada um dos tres e
         entao conferido integrando a EDO. Uma versao anterior varria uma grade de
         fatores de mimetico com passo de 0,247 e reportava o primeiro que alcancava o
         nivel do controle, 1,494, que ultrapassa a exigencia exata e fazia a queda
         implicada parecer uma faixa de 30-33% quando e 29,2% em todo valor do parametro
         livre.
    """
    sol_c, _ = run_ad(measured, "control", aggregation=False)
    sol_d, _ = run_ad(measured, "AD", aggregation=False)
    ab_c = sol_c.y[5][-1]
    ab_d = sol_d.y[5][-1]

    kappa = required_bace1_knockdown(measured)
    d = measured["d_mRNA_BACE1"]
    k_mid = midpoint("k_repress")
    k_lo, k_hi = FREE["k_repress"]["low"], FREE["k_repress"]["high"]
    folds = {name: mimic_for_knockdown(d, k, kappa)
             for name, k in (("low", k_lo), ("mid", k_mid), ("high", k_hi))}

    # EN/PT: integrate the ODE at each closed-form mimic factor; the AD monomer pool must
    #        land on the control pool of the same free-parameter setting
    errors = []
    for name, k in (("low", k_lo), ("mid", k_mid), ("high", k_hi)):
        free = {kk: midpoint(kk) for kk in FREE}
        free["k_repress"] = k
        s_c, _ = run_ad(measured, "control", free=free, aggregation=False)
        s_d, _ = run_ad(measured, "AD", free=free, mimic=folds[name], aggregation=False)
        errors.append(abs(s_d.y[5][-1] / s_c.y[5][-1] - 1.0))
    if max(errors) > 1e-6:
        sys.exit("EN/PT: the closed-form mimic factor does not bring the integrated AD "
                 f"monomer pool to the control pool (relative error {max(errors):.2e})")

    return dict(
        abeta_control=float(ab_c),
        abeta_AD=float(ab_d),
        abeta_ratio_AD_over_control=float(ab_d / ab_c) if ab_c else float("nan"),
        clearance_deficit_percent=float(
            100.0 * (1.0 - measured["k_clear_Ab_ad"] / measured["k_clear_Ab_ctrl"])),
        production_difference_percent=float(
            100.0 * (1.0 - measured["k_prod_Ab_ad"] / measured["k_prod_Ab_ctrl"])),
        required_bace1_knockdown=float(kappa),
        mimic_fold_needed=float(folds["mid"]),
        mimic_fold_needed_range=[float(folds["high"]), float(folds["low"])],
        mimic_fold_ode_check_max_relative_error=float(max(errors)),
    )


def experiment_mimic_washout(measured):
    """
    EN | How long does a single-dose miRNA mimic last, given the measured decay rates?
    PT | Quanto dura um mimetico de dose unica, dadas as taxas de decaimento medidas?
    """
    out = {}
    for label, d in [("miR-7 (1.7 h, K020)", measured["d_miR7"]),
                     ("miR-7 upper bound (5 h, K021)", measured["d_miR7_slow"]),
                     ("miR-29b (7 h, K017)", measured["d_miR29b"]),
                     ("miR-29c (10.6 h, K018)", measured["d_miR29c"]),
                     ("median miRNA (34 h, K016)", measured["d_miR_median"])]:
        # EN | Time for a 10-fold bolus to fall back within 10% of baseline: the excess
        #      starts at 9x baseline and must decay to 0.1x, so t = ln(9 / 0.1) / d.
        #      An earlier version used ln(9) / d, which is the time to reach TWICE
        #      baseline, and so reported washout times about half as long as the
        #      threshold drawn in the figure (dashed line at 1.1x).
        # PT | Tempo para um bolus de 10x voltar a menos de 10% do basal: o excesso
        #      comeca em 9x o basal e precisa cair a 0,1x, entao t = ln(9 / 0,1) / d.
        #      Uma versao anterior usava ln(9) / d, que e o tempo para chegar ao DOBRO
        #      do basal, e por isso reportava eliminacoes com cerca de metade do tempo
        #      do limiar desenhado na figura (linha tracejada em 1,1x).
        excess0, tolerance = 9.0, 0.1
        t = math.log(excess0 / tolerance) / d if d > 0 else float("nan")
        out[label] = dict(decay_constant_per_hour=float(d),
                          half_life_hours=float(LN2 / d),
                          hours_to_return_to_baseline=float(t),
                          days_to_return_to_baseline=float(t / 24.0))
    return out


def experiment_human_aggregate_load(measured):
    """
    EN | Where the measured human brain aggregate loads sit relative to the critical
         fibril concentration above which secondary nucleation takes over. This is a
         comparison of measured quantities, not a simulation output, so it does not
         depend on any free parameter.
    PT | Onde as cargas de agregado medidas em cerebro humano ficam em relacao a
         concentracao critica de fibrila acima da qual a nucleacao secundaria assume.
         E uma comparacao entre quantidades medidas, nao saida de simulacao, entao nao
         depende de nenhum parametro livre.
    """
    mstar = measured["M_star_uM"] * 1e-6          # uM -> M
    out = {}
    for lab, key in [("control median [K032]", "Ab_load_ctrl_median"),
                     ("control lower quartile [K031]", "Ab_load_ctrl_lq"),
                     ("AD median [K035]", "Ab_load_AD_median"),
                     ("AD upper quartile [K036]", "Ab_load_AD_uq")]:
        out[lab] = dict(load_nM=measured[key] * 1e9,
                        multiples_of_M_star=measured[key] / mstar,
                        above_M_star=bool(measured[key] > mstar))
    out["M_star_nM [K014]"] = mstar * 1e9
    out["AD_over_control_median_fold"] = (measured["Ab_load_AD_median"]
                                          / measured["Ab_load_ctrl_median"])
    return out


def experiment_ph_gate(measured):
    """
    EN | What the pH gate does to aggregate number, holding everything else fixed.
    PT | O que a trava de pH faz com o numero de agregados, com o resto fixo.
    """
    neutral, _ = run_pd(measured, acidic=False)
    acidic, _ = run_pd(measured, acidic=True)
    n0 = float(neutral.y[3][-1])
    # EN | The size of the switch, swept over four orders of magnitude of k2_aSyn.
    #      If the fold change moves with k2, the magnitude is not a result.
    # PT | O tamanho da chave, varrido em quatro ordens de grandeza de k2_aSyn. Se a
    #      razao muda com k2, a magnitude nao e resultado.
    sweep = []
    for k2 in np.geomspace(1e-8, 1e-4, 9):
        a, _ = run_pd(measured, acidic=True, k2=float(k2))
        sweep.append(dict(k2_aSyn=float(k2),
                          fold=float(a.y[3][-1] / n0) if n0 else float("inf")))
    folds = [x["fold"] for x in sweep]
    return dict(
        fibril_number_neutral_pH=n0,
        fibril_number_acidic_pH=float(acidic.y[3][-1]),
        fold_difference_at_illustrative_values=float(acidic.y[3][-1] / n0) if n0 else float("inf"),
        fibril_mass_neutral_pH=float(neutral.y[4][-1]),
        fibril_mass_acidic_pH=float(acidic.y[4][-1]),
        neutral_stays_at_seed=bool(abs(n0 - 1e-3) < 1e-12),
        acidic_above_neutral_in_every_k2=bool(all(f > 1.0 for f in folds)),
        fold_range_over_k2_sweep=[float(min(folds)), float(max(folds))],
        k2_sweep=sweep,
        reading_en="Direction sourced (K011, K042); magnitude set by illustrative constants.",
        reading_pt="Direcao com fonte (K011, K042); magnitude fixada por constantes ilustrativas.",
    )


def experiment_alpha_synuclein_saturation(measured):
    """
    EN | The one place in this project where two independent measurements meet with no
         free parameter between them. Buell et al. measured that alpha-synuclein fibril
         elongation saturates, with half-maximal rate near 46-50 uM of monomer (K009,
         K010). Wilhelm et al. measured how much alpha-synuclein is actually in a
         presynaptic bouton: 21.6 uM (K054, derived from the combined synuclein number
         and the alpha-to-beta ratio). Putting one on the other says which regime the
         protein sits in, and therefore how much a change in its level matters.
    PT | O unico ponto deste projeto em que duas medidas independentes se encontram sem
         nenhum parametro livre entre elas. Buell et al. mediram que a elongacao da
         fibrila de alfa-sinucleina satura, com meia-taxa perto de 46-50 uM de monomero
         (K009, K010). Wilhelm et al. mediram quanta alfa-sinucleina existe de fato num
         botao presinaptico: 21,6 uM (K054, derivado do numero combinado de sinucleina e
         da razao alfa para beta). Por um sobre o outro se descobre em que regime a
         proteina esta, e portanto o quanto uma mudanca no seu nivel importa.

    EN | Two systems, and the reader must be told: the saturation constant is recombinant
         protein in vitro, the concentration is a rat synaptosome. The comparison is an
         order-of-magnitude statement about regime, not a rate prediction.
    PT | Dois sistemas, e o leitor precisa saber: a constante de saturacao e proteina
         recombinante in vitro, a concentracao e um sinaptossomo de rato. A comparacao e
         uma afirmacao de ordem de grandeza sobre regime, nao predicao de taxa.
    """
    C = measured["C_aSyn_uM"]
    out = dict(alpha_synuclein_uM_in_bouton=float(C), points=[])
    for label, m_half in (("K009 (pH 7.4, 45 C): 49.8 uM", measured["m_half_uM"]),
                          ("K010 (pH 7.4, 37 C): 46 uM", measured["m_half_uM_b"])):
        frac = C / (m_half + C)          # EN/PT: rate as a fraction of maximal
        elasticity = m_half / (m_half + C)   # EN/PT: d(ln rate) / d(ln concentration)
        out["points"].append(dict(
            m_half_source=label, m_half_uM=float(m_half),
            concentration_over_m_half=float(C / m_half),
            fraction_of_maximal_elongation_rate=float(frac),
            percent_rate_change_per_percent_concentration_change=float(elasticity)))
    fr = [p["fraction_of_maximal_elongation_rate"] for p in out["points"]]
    el = [p["percent_rate_change_per_percent_concentration_change"] for p in out["points"]]
    out["below_half_saturation"] = bool(max(fr) < 0.5)
    out["fraction_of_maximal_range"] = [float(min(fr)), float(max(fr))]
    out["elasticity_range"] = [float(min(el)), float(max(el))]
    out["reading_en"] = (
        "Physiological alpha-synuclein sits below half-saturation, on the rising part of "
        "the elongation curve, so lowering it lowers elongation almost in proportion. "
        "That is what makes a miR-7 mimic a coherent idea at all. The saturation constant "
        "is from recombinant protein in vitro and the concentration from rat synaptosomes, "
        "so this is a statement about regime, not a predicted rate.")
    out["reading_pt"] = (
        "A alfa-sinucleina fisiologica fica abaixo da meia-saturacao, na parte ascendente "
        "da curva de elongacao, entao baixa-la reduz a elongacao quase em proporcao. E o "
        "que torna um mimetico de miR-7 uma ideia coerente. A constante de saturacao vem "
        "de proteina recombinante in vitro e a concentracao de sinaptossomo de rato, "
        "entao isto fala de regime, nao de taxa predita.")
    return out


def experiment_mir7_to_elongation(measured):
    """
    EN | The longest chain in this project that never touches a free parameter. Three
         independent measurements, in this order:
           Doxakis 2010  - how much miR-7 lowers alpha-synuclein (K058, K063)
           Wilhelm 2014  - how much alpha-synuclein is in a presynaptic bouton (K054)
           Buell 2014    - how fibril elongation depends on that concentration (K009, K010)
         Together they answer a question neither answers alone: what does a measured
         miR-7 knockdown do to the rate at which fibrils grow?
    PT | A cadeia mais longa deste projeto que nao toca em nenhum parametro livre. Tres
         medidas independentes, nesta ordem:
           Doxakis 2010  - quanto o miR-7 baixa a alfa-sinucleina (K058, K063)
           Wilhelm 2014  - quanta alfa-sinucleina existe num botao presinaptico (K054)
           Buell 2014    - como a elongacao da fibrila depende dessa concentracao (K009, K010)
         Juntas respondem a uma pergunta que nenhuma responde sozinha: o que uma queda
         medida de miR-7 faz com a velocidade de crescimento das fibrilas?

    EN | The answer is less than proportional, and that is the point. Because elongation
         saturates, a knockdown buys a smaller fractional slowdown than the knockdown
         itself. It is worth more than a hand-wave in either direction: the effect is
         neither negligible nor one-for-one.
    PT | A resposta e menos que proporcional, e e esse o ponto. Como a elongacao satura,
         uma queda compra um retardo fracional menor que a propria queda. Vale mais que
         um aceno em qualquer direcao: o efeito nao e desprezivel nem de um para um.
    """
    C0 = measured["C_aSyn_uM"]

    def rate(S, m_half):
        """EN/PT: elongation rate as a fraction of its maximum, Buell saturation form."""
        return S / (m_half + S)

    out = dict(alpha_synuclein_uM_in_bouton=float(C0), cases=[])
    for kd_label, kd in (("miR-7 alone, SNCA construct with 3'UTR [K058]", measured["kd_SNCA_miR7_construct"]),
                         ("miR-7 + miR-153, cortical neurons [K063]", measured["kd_SNCA_neurons"])):
        for mh_label, m_half in (("K009: 49.8 uM", measured["m_half_uM"]),
                                 ("K010: 46 uM", measured["m_half_uM_b"])):
            C1 = C0 * (1.0 - kd)
            slow = 1.0 - rate(C1, m_half) / rate(C0, m_half)
            out["cases"].append(dict(
                knockdown_source=kd_label, knockdown=float(kd),
                m_half_source=mh_label, m_half_uM=float(m_half),
                alpha_synuclein_after_uM=float(C1),
                elongation_rate_reduction=float(slow),
                reduction_over_knockdown=float(slow / kd)))
    slows = [c["elongation_rate_reduction"] for c in out["cases"]]
    out["elongation_rate_reduction_range"] = [float(min(slows)), float(max(slows))]
    out["less_than_proportional_everywhere"] = bool(
        all(c["elongation_rate_reduction"] < c["knockdown"] for c in out["cases"]))

    # EN | The paper calls the two miRNAs additive. Independent action on what the other
    #      leaves would predict 1 - (1-a)(1-b); compare that with the measured pair.
    # PT | O artigo chama os dois miRNAs de aditivos. Acao independente sobre o que o
    #      outro deixa preveria 1 - (1-a)(1-b); compare com o par medido.
    a, b = measured["kd_SNCA_miR7_construct"], measured["kd_SNCA_miR153_construct"]
    pair = measured["kd_SNCA_pair_construct"]
    out["additivity"] = dict(
        miR7_alone=float(a), miR153_alone=float(b), pair_measured=float(pair),
        pair_if_independent=float(1.0 - (1.0 - a) * (1.0 - b)),
        excess_over_independent=float(pair - (1.0 - (1.0 - a) * (1.0 - b))),
        reading_en=("Independent action predicts the pair almost exactly, so the data do "
                    "not require synergy between the two miRNAs at the protein level."),
        reading_pt=("Acao independente preve o par quase exatamente, entao os dados nao "
                    "exigem sinergia entre os dois miRNAs no nivel da proteina."))

    out["reading_en"] = (
        "Three systems in one chain: a knockdown measured in a cell line and in rat "
        "cortical neurons, a concentration measured in rat synaptosomes, and a saturation "
        "curve measured on recombinant protein in vitro. The chain is arithmetic on "
        "measured numbers with no free parameter anywhere, but it is not a prediction for "
        "a human brain, and each step carries the system it came from.")
    out["reading_pt"] = (
        "Tres sistemas numa cadeia so: uma queda medida em linhagem celular e em neuronios "
        "corticais de rato, uma concentracao medida em sinaptossomos de rato e uma curva de "
        "saturacao medida em proteina recombinante in vitro. A cadeia e aritmetica sobre "
        "numeros medidos, sem nenhum parametro livre, mas nao e predicao para cerebro "
        "humano, e cada passo carrega o sistema de onde veio.")
    return out


def experiment_enzyme_substrate_stoichiometry(measured):
    """
    EN | BACE1 and its substrate counted in the same preparation, in the same units.
    PT | BACE1 e seu substrato contados na mesma preparacao, nas mesmas unidades.
    """
    n_b, n_a = measured["N_BACE1"], measured["N_APP"]
    return dict(
        BACE1_copies_per_bouton=float(n_b),
        APP_copies_per_bouton=float(n_a),
        APP_per_BACE1=float(n_a / n_b),
        reading_en=("The enzyme is outnumbered by its substrate about fifty to one in the "
                    "presynaptic bouton, which is the arrangement in which the enzyme's own "
                    "abundance sets the flux. That is the arrangement a miR-29 mimic would "
                    "act on. Rat synaptosomes; no human equivalent is in this table."),
        reading_pt=("A enzima e superada pelo substrato em cerca de cinquenta para um no "
                    "botao presinaptico, que e o arranjo em que a abundancia da propria "
                    "enzima fixa o fluxo. E o arranjo sobre o qual um mimetico de miR-29 "
                    "agiria. Sinaptossomos de rato; nao ha equivalente humano nesta tabela."))


def experiment_mimic_versus_measured_knockdown(measured, clearance, n=9):
    """
    EN | Put the dose the model asks for next to a dose that has actually been achieved.
         The model says how much miR-29 has to rise to bring AD Abeta back to the control
         steady state (clearance experiment). Hebert et al. 2008 measured what miR-29a/b-1
         transfection does to BACE1 protein in cells: about 50% knockdown (K047). The two
         meet on BACE1, so the comparison is: how large is the BACE1 knockdown the model
         requires, and how large a mimic would be needed to reproduce the measured 50%?
         The required knockdown is exact and free of free parameters (it is 1 - 1/R, see
         required_bace1_knockdown); only the mimic factor that delivers it depends on
         k_repress.
    PT | Poe a dose que o modelo pede ao lado de uma dose de fato alcancada. O modelo diz
         quanto o miR-29 precisa subir para levar o Abeta de AD ao estado do controle
         (experimento de depuracao). Hebert et al. 2008 mediram o que a transfeccao de
         miR-29a/b-1 faz com a proteina BACE1 em celulas: cerca de 50% de queda (K047).
         Os dois se encontram na BACE1, entao a comparacao e: quao grande e a queda de
         BACE1 que o modelo exige, e que dose reproduziria os 50% medidos? A queda
         exigida e exata e livre de parametros livres (e 1 - 1/R, ver
         required_bace1_knockdown); so o fator de mimetico que a entrega depende de
         k_repress.

    EN | At steady state the BACE1 mRNA level is s_mRNA / (d_mRNA + 3 k_repress m) for a
         mimic factor m, because each of the three paralogues settles at m times baseline.
         Protein follows mRNA linearly, so the knockdown does not depend on k_translate
         and the whole comparison is analytic. Since Tushev et al. 2018 the mRNA decay is
         measured (K068), so the only free quantity left here is k_repress, reported
         across its declared range.
    PT | No estado estacionario o mRNA de BACE1 vale s_mRNA / (d_mRNA + 3 k_repress m)
         para um fator de mimetico m, porque cada um dos tres paralogos se estabiliza em m
         vezes o basal. A proteina acompanha o mRNA linearmente, entao a queda nao depende
         de k_translate e a comparacao inteira e analitica. Desde Tushev et al. 2018 o
         decaimento do mRNA e medido (K068), entao a unica grandeza livre aqui e o
         k_repress, reportado em toda a faixa declarada.
    """
    kappa = clearance["required_bace1_knockdown"]
    measured_kd = measured["BACE1_knockdown_miR29"]      # K047, fraction

    # EN/PT: d is measured now (K068), so only k_repress is swept: n*n values of it
    d_fixed = measured["d_mRNA_BACE1"]
    ks = [float(k) for k in
          np.geomspace(FREE["k_repress"]["low"], FREE["k_repress"]["high"], n * n)]
    offset = [mimic_for_knockdown(d_fixed, k, kappa) for k in ks]
    reproduce = [mimic_for_knockdown(d_fixed, k, measured_kd) for k in ks]
    # EN/PT: the knockdown that each required mimic produces must be kappa, by construction
    back = [bace1_knockdown(d_fixed, k, m) for k, m in zip(ks, offset)]
    if max(abs(b - kappa) for b in back) > 1e-12:
        sys.exit("EN/PT: the closed-form mimic factor does not return the required knockdown")

    out = dict(
        required_BACE1_knockdown=float(kappa),
        d_mRNA_BACE1_measured=float(d_fixed),
        measured_BACE1_knockdown_K047=measured_kd,
        required_over_measured_knockdown=float(kappa / measured_kd),
        required_knockdown_below_measured=bool(kappa < measured_kd),
        k_repress_values=len(ks),
        mimic_fold_to_offset_clearance=[float(min(offset)), float(max(offset))],
        mimic_fold_to_reproduce_measured_knockdown=[float(min(reproduce)), float(max(reproduce))],
    )
    out["reading_en"] = (
        "The comparison is between a modelled dose and a measured effect in a cell line, "
        "not a prediction of what a mimic would do in a human brain. It says only whether "
        "the intervention the model asks for is larger or smaller than one already shown "
        "to be achievable. The transfected miRNA level in the cell experiment is not "
        "reported as a fold over baseline, so the 50% constrains the product of repression "
        "strength and dose, not either alone.")
    out["reading_pt"] = (
        "A comparacao e entre uma dose modelada e um efeito medido em linhagem celular, "
        "nao uma predicao do que um mimetico faria num cerebro humano. Ela diz apenas se a "
        "intervencao que o modelo pede e maior ou menor que uma ja demonstrada como "
        "alcancavel. O nivel de miRNA transfectado no experimento celular nao e "
        "reportado como razao sobre o basal, entao os 50% restringem o produto da forca "
        "de repressao pela dose, e nenhum dos dois sozinho.")
    return out


def experiment_free_parameter_sensitivity(measured, n=9):
    """
    EN | Do the qualitative conclusions survive the free parameters? Each is scanned
         across its full range while the others sit at their midpoint.
    PT | As conclusoes qualitativas sobrevivem aos parametros livres? Cada um e varrido
         em toda a faixa enquanto os outros ficam no ponto medio.
    """
    rows = []
    base = {k: midpoint(k) for k in FREE}
    for name, spec in FREE.items():
        for v in np.geomspace(spec["low"], spec["high"], n):
            free = dict(base)
            free[name] = float(v)
            # EN | Monomer cascade only (aggregation=False). The question here is the
            #      monomer ratio, which the illustrative aggregation block cannot change;
            #      with that block included, the corners k_repress = 0.1 and
            #      k_translate = 5 diverged and were reported as not evaluable.
            # PT | So a cascata do monomero (aggregation=False). A pergunta aqui e a razao
            #      de monomero, que o bloco ilustrativo de agregacao nao pode mudar; com
            #      esse bloco incluido, os cantos k_repress = 0,1 e k_translate = 5
            #      divergiam e eram reportados como nao avaliaveis.
            sc, _ = run_ad(measured, "control", free=free, aggregation=False)
            sd, _ = run_ad(measured, "AD", free=free, aggregation=False)
            if not (sc.success and sd.success):
                ratio = float("nan")
            else:
                ratio = sd.y[5][-1] / sc.y[5][-1] if sc.y[5][-1] else float("nan")
            # EN | Three outcomes, not two. A non-finite ratio means the illustrative
            #      aggregation block diverged at this corner of the scan; it is NOT
            #      evidence against the conclusion. An earlier version counted such a
            #      NaN as a flip and printed "CONCLUSION NOT ROBUST" off a single
            #      numerical blow-up.
            # PT | Tres desfechos, nao dois. Um ratio nao finito significa que o bloco
            #      ilustrativo de agregacao divergiu naquele canto da varredura; NAO e
            #      evidencia contra a conclusao. Uma versao anterior contava esse NaN
            #      como inversao e imprimia "CONCLUSAO NAO ROBUSTA" por causa de uma
            #      unica explosao numerica.
            if not math.isfinite(ratio):
                verdict = "not_evaluable"
            elif ratio > 1.0:
                verdict = "AD_above_control"
            else:
                verdict = "AD_at_or_below_control"
            rows.append(dict(free_parameter=name, value=float(v),
                             abeta_AD_over_control=float(ratio),
                             verdict=verdict))
    return rows


def figures(measured, lang):
    os.makedirs(FIG_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # --- miR-29 paralogues decay apart ---------------------------------------
    ax = axes[0, 0]
    t = np.linspace(0, 48, 400)
    for lab, d, pid in [(tr(lang, "miR-29a (no measurable decay)",
                            "miR-29a (sem decaimento mensurável)"),
                         measured["d_miR29a_proxy"], "K019"),
                        ("miR-29b (7 h)", measured["d_miR29b"], "K017"),
                        ("miR-29c (10,6 h)" if lang != "en" else "miR-29c (10.6 h)",
                         measured["d_miR29c"], "K018")]:
        ax.plot(t, np.exp(-d * t), label=f"{lab}  [{pid}]")
    ax.set_title(tr(lang, "miR-29 paralogues are not one species",
                   "Parálogos de miR-29 não são uma espécie só"), fontsize=9)
    ax.set_xlabel(tr(lang, "hours", "horas"))
    ax.set_ylabel(tr(lang, "fraction remaining", "fração restante"))
    ax.legend(fontsize=7); ax.grid(alpha=.3)

    # --- miR-7 washout -------------------------------------------------------
    ax = axes[0, 1]
    t = np.linspace(0, 120, 400)
    for lab, d in [(tr(lang, "miR-7 (1.7 h) [K020]", "miR-7 (1,7 h) [K020]"),
                    measured["d_miR7"]),
                   (tr(lang, "miR-7 upper bound (5 h) [K021]",
                       "miR-7 limite superior (5 h) [K021]"), measured["d_miR7_slow"]),
                   (tr(lang, "median miRNA (34 h) [K016]",
                       "miRNA mediano (34 h) [K016]"), measured["d_miR_median"])]:
        ax.plot(t, 1.0 + 9.0 * np.exp(-d * t), label=lab)
    ax.axhline(1.1, ls="--", c="k", lw=.8)
    ax.set_yscale("log")
    ax.set_title(tr(lang, "A 10-fold mimic bolus, by measured decay",
                   "Bolus de mimético de 10x, pelo decaimento medido"), fontsize=9)
    ax.set_xlabel(tr(lang, "hours", "horas"))
    ax.set_ylabel(tr(lang, "fold over baseline", "vezes o basal"))
    ax.legend(fontsize=7); ax.grid(alpha=.3)

    # --- Abeta: clearance vs production --------------------------------------
    ax = axes[1, 0]
    sc, _ = run_ad(measured, "control")
    sd, _ = run_ad(measured, "AD")
    ax.plot(sc.t / 24.0, sc.y[5], label=tr(lang, "control", "controle"))
    ax.plot(sd.t / 24.0, sd.y[5],
            label=tr(lang, "AD (clearance -30%)", "DA (depuração -30%)"))
    ax.set_title(tr(lang, "Abeta42: the human lever is clearance [K001-K004]",
                   "Aβ42: a alavanca humana é a depuração [K001-K004]"), fontsize=9)
    ax.set_xlabel(tr(lang, "days", "dias"))
    ax.set_ylabel(tr(lang, "Abeta (relative)", "Aβ (relativo)"))
    ax.legend(fontsize=7); ax.grid(alpha=.3)

    # --- alpha-syn aggregation, pH gate --------------------------------------
    ax = axes[1, 1]
    sn, _ = run_pd(measured, acidic=False)
    sa, _ = run_pd(measured, acidic=True)
    ax.plot(sn.t / 24.0, sn.y[4],
            label=tr(lang, "neutral pH: no self-amplification [K011]",
                    "pH neutro: sem autoamplificação [K011]"))
    ax.plot(sa.t / 24.0, sa.y[4],
            label=tr(lang, "acidic pH: secondary nucleation on [K042]",
                    "pH ácido: nucleação secundária ligada [K042]"))
    ax.set_title(tr(lang, "alpha-synuclein fibril mass",
                   "Massa de fibrila de α-sinucleína"), fontsize=9)
    ax.set_xlabel(tr(lang, "days", "dias"))
    ax.set_ylabel(tr(lang, "fibril mass (relative)", "massa de fibrila (relativa)"))
    ax.legend(fontsize=7); ax.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(fig_path(FIG_DIR, "ode_calibrated_overview", lang), dpi=200)
    plt.close(fig)


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    measured, table = load_measured()

    print("=" * 78)
    print("EN | Recalibrated ODE models | PT | Modelos EDO recalibrados")
    print("=" * 78)
    n_used = sum(1 for r in table.values() if r["kind"] == "numeric")
    print(f"Measured constants available / constantes medidas : {n_used}")
    print(f"Free parameters / parametros livres               : {len(FREE)}")
    for k, v in FREE.items():
        print(f"    {k:<14} scan {v['low']:.4g} .. {v['high']:.4g}")

    clearance = experiment_clearance_vs_production(measured)
    washout = experiment_mimic_washout(measured)
    ph = experiment_ph_gate(measured)
    loads = experiment_human_aggregate_load(measured)
    sens = experiment_free_parameter_sensitivity(measured)
    dose = experiment_mimic_versus_measured_knockdown(measured, clearance)
    sat = experiment_alpha_synuclein_saturation(measured)
    stoich = experiment_enzyme_substrate_stoichiometry(measured)
    chain = experiment_mir7_to_elongation(measured)

    print("\n--- AD: clearance versus production ---")
    print(f"  Abeta AD / control            : {clearance['abeta_ratio_AD_over_control']:.3f}")
    print(f"  measured clearance deficit    : {clearance['clearance_deficit_percent']:.1f}%")
    print(f"  measured production difference: {clearance['production_difference_percent']:.1f}%")
    print(f"  required BACE1 knockdown      : {100 * clearance['required_bace1_knockdown']:.1f}%"
          f"  (1 - 1/ratio, no free parameter)")
    lo_m, hi_m = clearance["mimic_fold_needed_range"]
    print(f"  miR-29 mimic fold needed      : {clearance['mimic_fold_needed']:.3f} at the midpoint "
          f"of k_repress, {lo_m:.3f} .. {hi_m:.3f} across its range")
    print(f"  ODE check, max relative error : {clearance['mimic_fold_ode_check_max_relative_error']:.1e}")

    print("\n--- Single-dose mimic washout ---")
    for lab, d in washout.items():
        print(f"  {lab:<32} t_half {d['half_life_hours']:6.1f} h   "
              f"back to baseline in {d['days_to_return_to_baseline']:6.2f} d")

    print("\n--- Human brain Abeta42 aggregate load vs the critical concentration ---")
    print(f"  M* (secondary nucleation takes over) : {loads['M_star_nM [K014]']:.4g} nM")
    for lab in ["control lower quartile [K031]", "control median [K032]",
                "AD median [K035]", "AD upper quartile [K036]"]:
        d = loads[lab]
        print(f"  {lab:<32} {d['load_nM']:8.4g} nM   "
              f"{d['multiples_of_M_star']:8.4g} x M*   "
              f"{'above' if d['above_M_star'] else 'below'}")
    print(f"  AD / control, median                 : "
          f"{loads['AD_over_control_median_fold']:.4g} fold")

    print("\n--- alpha-synuclein pH gate ---")
    print(f"  fibril number, neutral pH     : {ph['fibril_number_neutral_pH']:.4g}")
    print(f"  fibril number, acidic pH      : {ph['fibril_number_acidic_pH']:.4g}")
    print(f"  neutral stays at seed         : {ph['neutral_stays_at_seed']}")
    print(f"  acidic > neutral, every k2    : {ph['acidic_above_neutral_in_every_k2']}")
    lo, hi = ph["fold_range_over_k2_sweep"]
    print(f"  fold range over k2 sweep      : {lo:.3g} .. {hi:.3g}")
    print("  -> direction is sourced; magnitude is not a result | "
          "direcao tem fonte; magnitude nao e resultado")

    flips = [r for r in sens if r["verdict"] == "AD_at_or_below_control"]
    nanq = [r for r in sens if r["verdict"] == "not_evaluable"]
    ok = [r for r in sens if r["verdict"] == "AD_above_control"]
    ratios = [r["abeta_AD_over_control"] for r in ok]
    print("\n--- alpha-synuclein: measured level on the measured saturation curve ---")
    print(f"  in the presynaptic bouton [K054]      : {sat['alpha_synuclein_uM_in_bouton']:.1f} uM")
    for pt in sat["points"]:
        print(f"  vs m_half {pt['m_half_uM']:.1f} uM  ->  {100*pt['fraction_of_maximal_elongation_rate']:.0f}% "
              f"of maximal elongation, {100*pt['percent_rate_change_per_percent_concentration_change']:.0f}% "
              f"rate change per 100% level change")
    print(f"  below half-saturation                 : {sat['below_half_saturation']}")

    print("\n--- miR-7 knockdown carried through to fibril elongation ---")
    for c in chain["cases"]:
        print(f"  {c['knockdown_source']:<42} {100*c['knockdown']:.0f}% down, "
              f"m_half {c['m_half_uM']:.1f} uM  ->  elongation {100*c['elongation_rate_reduction']:.0f}% slower")
    print(f"  less than proportional everywhere     : {chain['less_than_proportional_everywhere']}")
    ad = chain["additivity"]
    print(f"  additivity: pair measured {100*ad['pair_measured']:.0f}%, "
          f"independent action predicts {100*ad['pair_if_independent']:.0f}%")

    print("\n--- BACE1 and APP in the same bouton [K055, K056] ---")
    print(f"  BACE1 copies                          : {stoich['BACE1_copies_per_bouton']:.0f}")
    print(f"  APP copies                            : {stoich['APP_copies_per_bouton']:.0f}")
    print(f"  substrate per enzyme                  : {stoich['APP_per_BACE1']:.0f} to 1")

    print("\n--- Required mimic dose vs a measured knockdown [K047] ---")
    print(f"  BACE1 knockdown the model requires    : {100*dose['required_BACE1_knockdown']:.1f}%"
          f"  (same at every k_repress)")
    print(f"  measured knockdown in cells [K047]    : {100*dose['measured_BACE1_knockdown_K047']:.0f}%"
          f"  (required is {100*dose['required_over_measured_knockdown']:.0f}% of it)")
    lo, hi = dose["mimic_fold_to_offset_clearance"]
    print(f"  mimic fold that delivers the required : {lo:.3f} .. {hi:.3f}"
          f"  across {dose['k_repress_values']} values of k_repress")
    lo, hi = dose["mimic_fold_to_reproduce_measured_knockdown"]
    print(f"  mimic fold to reproduce the measured  : {lo:.2f} .. {hi:.2f}")

    print("\n--- Free-parameter sensitivity ---")
    print(f"  scans run / varreduras        : {len(sens)}")
    print(f"  AD above control              : {len(ok)}")
    print(f"  AD at or below control        : {len(flips)}")
    print(f"  not evaluable                 : {len(nanq)}")
    if ratios:
        print(f"  ratio range across all scans  : {min(ratios):.4f} .. {max(ratios):.4f}")
    print("  " + ("conclusion holds throughout | conclusao se mantem"
                  if not flips else "CONCLUSION NOT ROBUST | CONCLUSAO NAO ROBUSTA"))

    out = dict(
        generated_on="2026-10-14",
        measured_parameter_source=KINETICS,
        free_parameters={k: {kk: (vv if not isinstance(vv, float) else float(vv))
                             for kk, vv in v.items()} for k, v in FREE.items()},
        illustrative_constants=ILLUSTRATIVE,
        ad_clearance_vs_production=clearance,
        mimic_washout=washout,
        alpha_synuclein_ph_gate=ph,
        human_abeta_aggregate_load=loads,
        mimic_dose_vs_measured_knockdown=dose,
        alpha_synuclein_saturation=sat,
        enzyme_substrate_stoichiometry=stoich,
        mir7_knockdown_to_elongation=chain,
        free_parameter_sensitivity_aggregation_block_excluded=True,
        free_parameter_sensitivity_scans=len(sens),
        free_parameter_sensitivity_ad_above_control=len(ok),
        free_parameter_sensitivity_flips=len(flips),
        free_parameter_sensitivity_not_evaluable=len(nanq),
        abeta_ratio_min=float(min(ratios)) if ratios else None,
        abeta_ratio_max=float(max(ratios)) if ratios else None,
        abeta_ratio_note_en=("The AD/control monomer ratio is analytically "
                             "(k_prod_AD/k_clear_AD)/(k_prod_ctrl/k_clear_ctrl) and is therefore "
                             "independent of every free parameter. The scan confirms this "
                             "numerically rather than discovering it."),
        abeta_ratio_note_pt=("O ratio de monomero AD/controle e analiticamente "
                             "(k_prod_AD/k_clear_AD)/(k_prod_ctrl/k_clear_ctrl) e, portanto, "
                             "independente de todo parametro livre. A varredura confirma isso "
                             "numericamente, em vez de descobrir."),
    )
    json.dump(out, open(f"{TAB_DIR}/ode_calibrated_results.json", "w"),
              ensure_ascii=False, indent=1)
    with open(f"{TAB_DIR}/ode_free_parameter_sensitivity.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(sens[0].keys()))
        w.writeheader(); w.writerows(sens)

    for lang in LANGS:
        figures(measured, lang)
    print(f"\nEN/PT -> {TAB_DIR}/ode_calibrated_results.json")
    print(f"EN/PT -> {TAB_DIR}/ode_free_parameter_sensitivity.csv")
    for lang in LANGS:
        print(f"EN/PT -> {fig_path(FIG_DIR, 'ode_calibrated_overview', lang)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
