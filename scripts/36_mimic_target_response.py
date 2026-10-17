#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | What a microRNA mimic given repeatedly does to the protein it is meant to lower,
     with the measured turnover of every step in the chain.
PT | O que um mimetico de microRNA dado repetidamente faz com a proteina que ele deve
     baixar, com a renovacao medida de cada etapa da cadeia.

EN | The question. scripts/12 says how far a target has to fall (BACE1 by 29.2% to offset
     the Abeta42 clearance deficit; alpha-synuclein by 30% or 43%, the knockdowns that
     miR-7 was measured to produce) and scripts/17 says what the half-life of each miRNA
     costs in peak concentration when a mimic is dosed repeatedly. Neither says whether
     a dosing regimen that holds the AVERAGE mimic level also holds the average protein
     level. It need not: the mRNA and the protein are filters with their own half-lives,
     and repression is not linear in the miRNA. This script closes that gap, and keeps
     every constant the models already use.

PT | A pergunta. O scripts/12 diz quanto um alvo precisa cair (BACE1 em 29,2% para
     compensar o deficit de depuracao do Abeta42; alfa-sinucleina em 30% ou 43%, as
     quedas que o miR-7 produziu nas medidas) e o scripts/17 diz o que a meia-vida de cada
     miRNA custa em concentracao de pico quando um mimetico e dosado repetidamente. Nenhum
     dos dois diz se um regime que sustenta o nivel MEDIO do mimetico sustenta tambem o
     nivel medio da proteina. Nao precisa: o mRNA e a proteina sao filtros com meias-vidas
     proprias, e a repressao nao e linear no miRNA. Este script fecha essa lacuna, e
     mantem toda constante que os modelos ja usam.

EN | The model. Every production term is divided by its baseline value, so the state is
     the level relative to the untreated steady state (1.0). With r0 the baseline decay of
     the mRNA that the miRNA is responsible for, rho = r0 / d the ratio of that decay to
     the mRNA's intrinsic decay d, E(t) the mimic excess in units of the baseline pool of
     the repressor, mu the mRNA and pi the protein:

         dE/dt  = -d_x E                     (and E -> E + D at every dose, T hours apart)
         dmu/dt = (d + r0) - (d + r0 (1 + E)) mu
         dpi/dt = d_p (mu - pi)
         dal/dt = k_clear (pi - al)          (Abeta42 monomer, Alzheimer arm only)

     k_translate divides out of every one of these, so the only free parameter left is
     k_repress, and it enters through rho alone: rho = 3 k_repress / d for the three
     miR-29 paralogues of the AD arm and rho = k_repress / d for the single miR-7 of the
     PD arm. For a constant excess E the steady-state mRNA is 1 / (1 + rho E / (1 + rho)),
     so a knockdown kappa needs E = kappa / (1 - kappa) * (1 + rho) / rho, which is
     scripts/12's closed form. Every other constant is the one scripts/12 loads: the
     mRNA decay (K068, K066), the protein decay (K025, K026, and the in vivo brain
     proteome median K029 as a sensitivity case), the Abeta42 clearance in AD (K004) and
     the miRNA decay of the delivered species (K016-K021, the same assumption as
     scripts/17: the mimic is cleared like the endogenous species).

PT | O modelo. Todo termo de producao e dividido pelo seu valor basal, entao o estado e o
     nivel relativo ao estado estacionario sem tratamento (1,0). Com r0 o decaimento basal
     do mRNA de que o miRNA e responsavel, rho = r0 / d a razao entre esse decaimento e o
     decaimento intrinseco d do mRNA, E(t) o excesso de mimetico em unidades do pool basal
     do repressor, mu o mRNA e pi a proteina:

         dE/dt  = -d_x E                     (e E -> E + D a cada dose, T horas depois)
         dmu/dt = (d + r0) - (d + r0 (1 + E)) mu
         dpi/dt = d_p (mu - pi)
         dal/dt = k_clear (pi - al)          (monomero de Abeta42, so no braco de AD)

     O k_translate se cancela em todas elas, entao o unico parametro livre que sobra e o
     k_repress, e ele entra so por rho: rho = 3 k_repress / d para os tres paralogos de
     miR-29 do braco de AD e rho = k_repress / d para o miR-7 unico do braco de DP. Para
     um excesso E constante o mRNA no estado estacionario e 1 / (1 + rho E / (1 + rho)),
     entao uma queda kappa pede E = kappa / (1 - kappa) * (1 + rho) / rho, que e a forma
     fechada do scripts/12. Toda outra constante e a que o scripts/12 carrega: o
     decaimento do mRNA (K068, K066), o da proteina (K025, K026, e a mediana do proteoma
     cerebral in vivo K029 como caso de sensibilidade), a depuracao do Abeta42 na AD
     (K004) e o decaimento do miRNA da especie entregue (K016-K021, a mesma suposicao do
     scripts/17: o mimetico e eliminado como a especie endogena).

EN | What is computed. For each scenario the dose D that makes the CYCLE-AVERAGED protein
     level fall by the required kappa at the periodic steady state of repeated dosing. The
     periodic state is exact, not simulated to convergence: for a given E(t) the mRNA,
     protein and lagged variables obey a linear system, so one cycle maps the state at
     its start by an affine map, and the fixed point solves a small linear system (see
     periodic_state). Cycle averages are integrated as extra states, so they carry the
     integrator's accuracy and not a quadrature's. Reported per scenario: the dose, the
     peak and mean mimic excess, the dose the steady-state closed form would have asked
     for and the ratio of the two (efficiency), the peak-to-trough swing of the protein
     and, in the AD arm, of the Abeta42 monomer, and in the PD arm the average slowdown of
     fibril elongation on the measured saturation curve (K009, K010, K054). A target the
     regimen cannot reach at any dose up to 10^4 times the baseline pool is reported as
     unreachable, with the knockdown that dose does reach.

PT | O que e calculado. Para cada cenario, a dose D que faz o nivel de proteina MEDIO NO
     CICLO cair o kappa exigido no estado periodico da dosagem repetida. O estado
     periodico e exato, nao simulado ate convergir: para um E(t) dado, o mRNA, a proteina
     e as variaveis defasadas obedecem a um sistema linear, entao um ciclo leva o estado
     do inicio por uma aplicacao afim, e o ponto fixo resolve um sistema linear pequeno
     (ver periodic_state). As medias no ciclo sao integradas como estados extras, entao
     carregam a precisao do integrador e nao a de uma quadratura. Reportado por cenario: a
     dose, o excesso de pico e medio do mimetico, a dose que a forma fechada do estado
     estacionario pediria e a razao entre as duas (eficiencia), a oscilacao pico a vale da
     proteina e, no braco de AD, do monomero de Abeta42, e no braco de DP o retardo medio da
     elongacao da fibrila na curva de saturacao medida (K009, K010, K054). Uma meta que o
     regime nao alcanca com nenhuma dose ate 10^4 vezes o pool basal e reportada como
     inalcancavel, com a queda que essa dose alcanca.

EN | What it is not. A bolus enters the repressor pool instantly, which is the most
     favourable onset a dose can have and the least favourable peak. Repression is the
     mass-action form of the models (linear in the miRNA, saturating only because the
     mRNA is used up), so it does not saturate the way Argonaute loading and target
     occupancy can. The delivered species is assumed to act like the endogenous pool
     (equal potency per molecule, for the three miR-29 paralogues as for miR-7 and
     miR-153 lumped under miR-7). None of this is measured, and the results are
     conditional on it. They are a design bound: the dose a regimen has to reach if the
     assumptions hold.

PT | O que nao e. Um bolus entra no pool do repressor instantaneamente, o que e o inicio
     mais favoravel que uma dose pode ter e o pico menos favoravel. A repressao e a forma
     de acao de massas dos modelos (linear no miRNA, saturando so porque o mRNA se
     esgota), entao nao satura como o carregamento de Argonauta e a ocupacao de alvos
     podem saturar. A especie entregue e suposta agir como o pool endogeno (mesma
     potencia por molecula, nos tres paralogos de miR-29 como no miR-7 e no miR-153
     agrupados sob o miR-7). Nada disso e medido, e os resultados dependem disso. Sao um
     limite de projeto: a dose que um regime precisa alcancar se as suposicoes valem.

EN | Where the unmeasured repression strength sits. rho is the one free quantity, and
     the table holds one observation that bears on it, from below only. Blocking
     endogenous miR-7 and miR-153 with antisense oligonucleotides raised an SNCA 3'UTR
     reporter by 44% in rat cortical neurons (K065). For a repressor whose baseline decay
     of the mRNA is rho times the intrinsic decay, removing all of it raises the level by a factor 1 + rho, and
     a blockade cannot be more than complete, so rho is at least 0.44 for that pair. It
     is a floor and not an estimate: the blockade was probably partial, the readout is a
     reporter, and two microRNAs are lumped, so it says nothing about miR-7 alone. The
     declared scan of k_repress (0.1 to 3) puts rho between 0.94 and 28 in the PD arm,
     above that floor but not down to it, so the PD scenarios are also run at rho = 0.44,
     and a separate sweep runs every scenario at 24 h over rho from 0.01 to 10^4. Weaker
     repression raises every dose roughly in proportion to 1 / rho and stronger
     repression approaches a limit, so the sweep reports, for each scenario, the smallest
     peak over all rho: a bound that holds whatever rho is. For the AD arm no
     loss-of-function number is in the table (K046 is qualitative), and rho runs from 6.6
     to 198 over the scan.

PT | Onde fica a forca de repressao nao medida. rho e a unica grandeza livre, e a tabela
     tem uma observacao que a toca, so por baixo. Bloquear o miR-7 e o miR-153 endogenos
     com oligonucleotideos antisense elevou em 44% um repórter da 3'UTR de SNCA em
     neuronios corticais de rato (K065). Para um repressor cujo decaimento basal do mRNA
     vale rho vezes o decaimento intrinseco, retirar tudo dele eleva o nivel por um fator 1 + rho, e um
     bloqueio nao pode passar de completo, entao rho e no minimo 0,44 para esse par. E um
     piso e nao uma estimativa: o bloqueio provavelmente foi parcial, a leitura e um
     repórter, e dois microRNAs estao agrupados, entao nada diz sobre o miR-7 sozinho. A
     varredura declarada de k_repress (0,1 a 3) poe rho entre 0,94 e 28 no braco de DP,
     acima desse piso mas sem chegar a ele, entao os cenarios de DP tambem rodam em rho =
     0,44, e uma varredura a parte roda todo cenario em 24 h com rho de 0,01 a 10^4. Uma
     repressao mais fraca eleva toda dose aproximadamente em proporcao a 1 / rho e uma
     mais forte se aproxima de um limite, entao a varredura reporta, para cada cenario, o
     menor pico sobre todo rho: um limite que vale seja qual for rho. Para o braco de AD
     nao ha numero de perda de funcao na tabela (o K046 e qualitativo), e rho vai de 6,6 a
     198 na varredura.

EN | A second, smaller result. Abeta42 nucleation rates scale with the monomer level to
     the measured reaction orders (K013 primary, K012 secondary, both 2), and the ratio of
     two rates at two monomer levels needs the orders and nothing else, so it is free of
     the Abeta42 rate constants that no source gives (K037). The AD over control monomer
     ratio of 1.4126 therefore corresponds to a nucleation rate ratio of 1.4126^2, and
     the 29.2% fall in the monomer that offsets the clearance deficit to a 49.9% fall in
     the nucleation rate. This crosses systems (in vitro orders applied to the monomer
     pool of the brain) and is reported as a scaling, not a rate.

PT | Um segundo resultado, menor. As taxas de nucleacao do Abeta42 escalam com o nivel de
     monomero pelas ordens de reacao medidas (K013 primaria, K012 secundaria, ambas 2), e
     a razao entre duas taxas em dois niveis de monomero precisa das ordens e de mais
     nada, entao e livre das constantes de velocidade do Abeta42 que nenhuma fonte da
     (K037). A razao de monomero DA sobre controle de 1,4126 corresponde assim a uma
     razao de taxa de nucleacao de 1,4126^2, e a queda de 29,2% no monomero que compensa
     o deficit de depuracao a uma queda de 49,9% na taxa de nucleacao. Isto atravessa
     sistemas (ordens in vitro aplicadas ao pool de monomero do cerebro) e e reportado
     como escalonamento, nao como taxa.

    python scripts/36_mimic_target_response.py
"""

import csv
import importlib.util
import json
import math
import os
import platform
import sys
from collections import OrderedDict

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = "results/tables"
LN2 = math.log(2.0)

DOSE_CAP = 1.0e4        # EN/PT: largest dose tried, in units of the baseline repressor pool
RIPPLE_POINTS = 3000    # EN/PT: grid on which the swing of one cycle is read
RHO_SWEEP = np.geomspace(0.01, 1.0e4, 31)       # EN/PT: five points per decade


def load_module(name, filename):
    """
    EN | Import a sibling script whose filename starts with a digit. Reuse, not copy: the
         constants and the closed forms stay in one place.
    PT | Importa um script irmao cujo nome comeca com digito. Reaproveita, nao copia: as
         constantes e as formas fechadas ficam num lugar so.
    """
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        sys.exit(f"EN/PT: missing {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------
# The cascade | A cascata
# --------------------------------------------------------------------------
class Cascade:
    """
    EN | The normalised cascade of the module docstring, with the integrals needed for
         cycle averages carried as states. State layout:
             [E, mu, pi, lag_1 .. lag_L, I_pi, I_lag_1 .. I_lag_L, I_v_1 .. I_v_V]
         where lag_i are first-order lags that follow pi (the Abeta42 monomer in the AD
         arm), v_j is the relative fibril-elongation rate on a saturation curve (PD arm)
         and I_x is the running integral of x, so the cycle average is I_x(T) / T.
    PT | A cascata normalizada do docstring do modulo, com as integrais necessarias para
         as medias do ciclo carregadas como estados. Layout do estado:
             [E, mu, pi, lag_1 .. lag_L, I_pi, I_lag_1 .. I_lag_L, I_v_1 .. I_v_V]
         onde lag_i sao defasagens de primeira ordem que seguem pi (o monomero de Abeta42
         no braco de AD), v_j e a taxa relativa de elongacao da fibrila numa curva de
         saturacao (braco de DP) e I_x e a integral corrente de x, entao a media do ciclo
         e I_x(T) / T.
    """

    def __init__(self, d, rho, dp, dx, lags=(), v_maps=()):
        self.d, self.rho, self.dp, self.dx = d, rho, dp, dx
        self.r0 = rho * d
        self.lags = tuple(lags)
        self.v_maps = tuple(v_maps)             # (S0, m_half) pairs, same units
        self.L = len(self.lags)
        self.n_lin = 2 + self.L                 # mu, pi, lags: the affine part
        self.n_int = 1 + self.L + len(self.v_maps)

    def rhs(self, t, y):
        d, r0, dp, dx = self.d, self.r0, self.dp, self.dx
        E, mu, pi = y[0], y[1], y[2]
        dy = [0.0] * len(y)
        dy[0] = -dx * E
        dy[1] = (d + r0) - (d + r0 * (1.0 + E)) * mu
        dy[2] = dp * (mu - pi)
        for i, kc in enumerate(self.lags):
            dy[3 + i] = kc * (pi - y[3 + i])
        b = 3 + self.L
        dy[b] = pi
        for i in range(self.L):
            dy[b + 1 + i] = y[3 + i]
        for j, (s0, mh) in enumerate(self.v_maps):
            dy[b + 1 + self.L + j] = pi * (mh + s0) / (mh + s0 * pi)
        return dy

    def one_cycle(self, x0, e_peak, T, dense=False):
        """EN/PT: integrate one dosing interval from state x0 (mu, pi, lags) just after a dose."""
        y0 = [e_peak] + list(x0) + [0.0] * self.n_int
        sol = solve_ivp(self.rhs, (0.0, T), y0, method="LSODA", rtol=1e-10, atol=1e-13,
                        dense_output=dense)
        if not sol.success:
            raise RuntimeError(f"integration failed (E_peak {e_peak:.3g}, T {T}): {sol.message}")
        return sol

    def periodic_state(self, e_peak, T):
        """
        EN | The state at the start of a cycle that the cycle returns to. For a given
             E(t) the (mu, pi, lags) block is linear, so one cycle maps x to Phi x + c,
             and the periodic state solves (I - Phi) x = c. Phi and c come from one run
             from zero and one run per unit vector.
        PT | O estado no inicio de um ciclo ao qual o ciclo retorna. Para um E(t) dado o
             bloco (mu, pi, defasagens) e linear, entao um ciclo leva x a Phi x + c, e o
             estado periodico resolve (I - Phi) x = c. Phi e c saem de uma rodada a
             partir de zero e uma por vetor unitario.
        """
        m = self.n_lin
        zero = [0.0] * m
        c = self.one_cycle(zero, e_peak, T).y[1:1 + m, -1]
        cols = []
        for i in range(m):
            e = list(zero)
            e[i] = 1.0
            cols.append(self.one_cycle(e, e_peak, T).y[1:1 + m, -1] - c)
        phi = np.array(cols).T
        return np.linalg.solve(np.eye(m) - phi, c)

    def cycle(self, dose, T):
        """
        EN | Everything one needs about the periodic cycle for a dose D: the averages
             over the cycle and the swing of the protein and of each lagged variable.
        PT | Tudo o que se precisa do ciclo periodico para uma dose D: as medias no ciclo
             e a oscilacao da proteina e de cada variavel defasada.
        """
        if dose <= 0.0:
            return dict(mean_pi=1.0, mean_lags=[1.0] * self.L, mean_v=[1.0] * len(self.v_maps),
                        ripple_pi=0.0, ripple_lags=[0.0] * self.L, e_peak=0.0, residual=0.0)
        e_peak = dose / (1.0 - math.exp(-self.dx * T))
        x0 = self.periodic_state(e_peak, T)
        sol = self.one_cycle(x0, e_peak, T, dense=True)
        yT = sol.y[:, -1]
        residual = float(np.max(np.abs(yT[1:1 + self.n_lin] - x0)))
        b = 3 + self.L
        mean_pi = yT[b] / T
        mean_lags = [yT[b + 1 + i] / T for i in range(self.L)]
        mean_v = [yT[b + 1 + self.L + j] / T for j in range(len(self.v_maps))]
        grid = sol.sol(np.linspace(0.0, T, RIPPLE_POINTS))
        pi_t = grid[2]
        ripple_pi = float((pi_t.max() - pi_t.min()) / mean_pi)
        ripple_lags = [float((grid[3 + i].max() - grid[3 + i].min()) / mean_lags[i])
                       for i in range(self.L)]
        return dict(mean_pi=float(mean_pi), mean_lags=[float(v) for v in mean_lags],
                    mean_v=[float(v) for v in mean_v], ripple_pi=ripple_pi,
                    ripple_lags=ripple_lags, e_peak=float(e_peak), residual=residual)

    def required_dose(self, kappa, T, e_lin):
        """
        EN | The dose that brings the cycle-averaged protein to 1 - kappa, or None when
             no dose up to DOSE_CAP does. The mean protein falls monotonically with the
             dose, so a bracket and Brent's method suffice.
        PT | A dose que leva a proteina media no ciclo a 1 - kappa, ou None quando nenhuma
             dose ate DOSE_CAP o faz. A proteina media cai monotonicamente com a dose,
             entao um intervalo e o metodo de Brent bastam.
        """
        def gap(dose):
            return 1.0 - self.cycle(dose, T)["mean_pi"] - kappa

        hi = max(2.0 * e_lin * self.dx * T, 1e-6)
        while gap(hi) < 0.0:
            hi *= 2.0
            if hi > DOSE_CAP:
                return None, 1.0 - self.cycle(DOSE_CAP, T)["mean_pi"]
        dose = brentq(gap, 0.0, hi, xtol=1e-13, rtol=1e-12, maxiter=200)
        return dose, kappa


# --------------------------------------------------------------------------
# Closed forms | Formas fechadas
# --------------------------------------------------------------------------
def required_excess(kappa, rho):
    """
    EN | Mimic excess (units of the baseline pool) that holds a steady knockdown kappa:
         kappa / (1 - kappa) * (1 + rho) / rho. Its limit for large rho is the
         proportional one, kappa / (1 - kappa).
    PT | Excesso de mimetico (unidades do pool basal) que sustenta uma queda estavel
         kappa: kappa / (1 - kappa) * (1 + rho) / rho. O limite para rho grande e o
         proporcional, kappa / (1 - kappa).
    """
    return kappa / (1.0 - kappa) * (1.0 + rho) / rho


def steady_knockdown(excess, rho):
    """EN/PT: knockdown at a constant excess E | queda num excesso E constante."""
    return 1.0 - 1.0 / (1.0 + rho * excess / (1.0 + rho))


# --------------------------------------------------------------------------
# Scenarios | Cenarios
# --------------------------------------------------------------------------
def build_axes(ode, dos, measured, table):
    """
    EN | The two axes with their measured constants. Nothing here is a number typed by
         hand: every decay, concentration and knockdown is loaded by identifier from
         scripts/12 and scripts/17, which load them from the parameter table.
    PT | Os dois eixos com suas constantes medidas. Nada aqui e numero digitado a mao:
         todo decaimento, concentracao e queda e carregado por identificador dos
         scripts/12 e scripts/17, que os carregam da tabela de parametros.
    """
    dec = dos.load_decays()
    k_lo, k_hi = ode.FREE["k_repress"]["low"], ode.FREE["k_repress"]["high"]
    k_mid = ode.midpoint("k_repress")
    kappa_ad = ode.required_bace1_knockdown(measured)
    s0 = measured["C_aSyn_uM"]
    v_maps = ((s0, measured["m_half_uM"]), (s0, measured["m_half_uM_b"]))

    ad = OrderedDict(
        name="AD", n_sp=3, d=measured["d_mRNA_BACE1"], lags=(measured["k_clear_Ab_ad"],),
        v_maps=(),
        rho_cases=[(lab, k, 3.0 * k / measured["d_mRNA_BACE1"])
                   for lab, k in (("scan low", k_lo), ("scan midpoint", k_mid), ("scan high", k_hi))],
        species=[("miR-29b", "K017"), ("miR-29c", "K018"), ("median miRNA", "K016")],
        protein=[("BACE1 in neurons, 19 h [K025]", measured["d_BACE1_neuron"]),
                 ("brain proteome median in vivo, 8.1 d [K029]", measured["d_proteome_brain"])],
        targets=[("required to restore control Abeta42 [K001-K004]", kappa_ad)])

    d_pd = measured["d_mRNA_SNCA"]
    k065 = float(table["K065"]["value_si"])
    if table["K065"]["kind"] != "numeric":
        sys.exit("EN/PT: K065 is not a numeric row")
    pd = OrderedDict(
        name="PD", n_sp=1, d=d_pd, lags=(), v_maps=v_maps,
        rho_cases=[("K065 lower bound", k065 * d_pd, k065)]
                  + [(lab, k, k / d_pd)
                     for lab, k in (("scan low", k_lo), ("scan midpoint", k_mid), ("scan high", k_hi))],
        species=[("miR-7", "K020"), ("miR-7 upper bound", "K021"), ("median miRNA", "K016")],
        protein=[("alpha-synuclein, SH-SY5Y, 24 h [K026]", measured["d_aSyn"]),
                 ("brain proteome median in vivo, 8.1 d [K029]", measured["d_proteome_brain"])],
        targets=[("miR-7 alone, SNCA construct with 3'UTR [K058]", measured["kd_SNCA_miR7_construct"]),
                 ("miR-7 + miR-153, cortical neurons [K063]", measured["kd_SNCA_neurons"])])
    for ax in (ad, pd):
        ax["decay"] = {k: dec[k]["decay_per_hour"] for k, _ in ax["species"]}
        ax["half_life"] = {k: dec[k]["half_life_hours"] for k, _ in ax["species"]}
    return OrderedDict(AD=ad, PD=pd), kappa_ad, k065


def run_scenarios(axes, intervals):
    """
    EN | Solve every scenario. One row per axis, repression strength, delivered species,
         dosing interval, protein half-life and target.
    PT | Resolve todo cenario. Uma linha por eixo, forca de repressao, especie entregue,
         intervalo de dosagem, meia-vida da proteina e meta.
    """
    rows = []
    for key, ax in axes.items():
        for rho_label, k_rep, rho in ax["rho_cases"]:
            for sp, pid in ax["species"]:
                dx = ax["decay"][sp]
                for T in intervals:
                    for prot_label, dp in ax["protein"]:
                        cas = Cascade(ax["d"], rho, dp, dx, lags=ax["lags"], v_maps=ax["v_maps"])
                        for tgt_label, kappa in ax["targets"]:
                            e_lin = required_excess(kappa, rho)
                            dose, reached = cas.required_dose(kappa, T, e_lin)
                            r = OrderedDict(
                                axis=key, rho_case=rho_label, k_repress=k_rep, rho=rho,
                                species=sp, param_id=pid, half_life_h=ax["half_life"][sp],
                                interval_h=T, protein_case=prot_label,
                                protein_half_life_h=LN2 / dp, target=tgt_label,
                                target_knockdown=kappa, reachable=dose is not None,
                                knockdown_reached=reached,
                                linear_excess=e_lin, linear_dose=e_lin * dx * T)
                            if dose is None:
                                for f in ("dose", "peak_excess", "peak_pool_over_baseline",
                                          "trough_pool_over_baseline", "mean_excess",
                                          "efficiency", "protein_ripple", "abeta_ripple",
                                          "elongation_reduction_m49_8",
                                          "elongation_reduction_m46_0", "residual"):
                                    r[f] = None
                            else:
                                c = cas.cycle(dose, T)
                                if c["residual"] > 1e-8:
                                    sys.exit("EN/PT: the periodic state did not close "
                                             f"(residual {c['residual']:.1e}) for {key} {rho_label} "
                                             f"{sp} T={T} {prot_label}")
                                r["dose"] = dose
                                r["peak_excess"] = c["e_peak"]
                                r["peak_pool_over_baseline"] = 1.0 + c["e_peak"]
                                r["trough_pool_over_baseline"] = 1.0 + c["e_peak"] * math.exp(-dx * T)
                                r["mean_excess"] = dose / (dx * T)
                                r["efficiency"] = r["linear_dose"] / dose
                                r["protein_ripple"] = c["ripple_pi"]
                                r["abeta_ripple"] = c["ripple_lags"][0] if c["ripple_lags"] else None
                                r["elongation_reduction_m49_8"] = (1.0 - c["mean_v"][0]) if c["mean_v"] else None
                                r["elongation_reduction_m46_0"] = (1.0 - c["mean_v"][1]) if c["mean_v"] else None
                                r["residual"] = c["residual"]
                                r["knockdown_reached"] = 1.0 - c["mean_pi"]
                            rows.append(r)
        print(f"  {key}: {sum(1 for r in rows if r['axis'] == key)} scenarios solved")
    return rows


# --------------------------------------------------------------------------
# Verification | Verificacao
# --------------------------------------------------------------------------
def verify(axes, ode, dos, measured, intervals):
    """
    EN | Independent checks of the machinery, each against something that does not use
         it. A failure stops the script rather than writing a table.
         1. Linear limit: a tiny dose must give the knockdown of the linear response,
            rho / (1 + rho) times the mean excess.
         2. Peak over mean excess must equal scripts/17's peak-to-average ratio.
         3. Quasi-continuous dosing (T = 0.05 h) must reproduce the steady-state closed
            form.
         4. The exact periodic state must agree with brute-force simulation over many
            cycles from the untreated state.
         5. Constant exposure through scripts/12's own ODE (run_ad, mimic factor f) must
            lower BACE1 by exactly the knockdown the closed form says.
         6. scripts/12's mimic_for_knockdown and this script's required_excess must be the
            same formula.
    PT | Checagens independentes da maquinaria, cada uma contra algo que nao a usa. Uma
         falha para o script em vez de escrever uma tabela.
         1. Limite linear: uma dose minuscula tem de dar a queda da resposta linear,
            rho / (1 + rho) vezes o excesso medio.
         2. O excesso de pico sobre o medio tem de igualar a razao pico sobre media do
            scripts/17.
         3. Dosagem quase continua (T = 0,05 h) tem de reproduzir a forma fechada do
            estado estacionario.
         4. O estado periodico exato tem de concordar com a simulacao bruta em muitos
            ciclos a partir do estado sem tratamento.
         5. Exposicao constante pela propria EDO do scripts/12 (run_ad, fator f) tem de
            baixar a BACE1 exatamente o que a forma fechada diz.
         6. O mimic_for_knockdown do scripts/12 e o required_excess deste script tem de
            ser a mesma formula.
    """
    out = OrderedDict()

    def record(name, error, tol, note):
        out[name] = OrderedDict(max_abs_error=float(error), tolerance=tol,
                                passed=bool(error <= tol), note=note)
        if error > tol:
            sys.exit(f"EN/PT: verification '{name}' failed: error {error:.3e} > {tol:.0e}")

    def rho_of(ax, label):
        return next(r for lab, _, r in ax["rho_cases"] if lab == label)

    # 1. linear limit, on the fast species of each axis, at the scan midpoint
    errs = []
    for key, ax in axes.items():
        sp = ax["species"][0][0]
        dx = ax["decay"][sp]
        rho = rho_of(ax, "scan midpoint")
        cas = Cascade(ax["d"], rho, ax["protein"][0][1], dx, lags=ax["lags"])
        T, dose = 24.0, 1e-4
        e_mean = dose / (dx * T)
        got = 1.0 - cas.cycle(dose, T)["mean_pi"]
        want = rho / (1.0 + rho) * e_mean
        errs.append(abs(got - want) / want)
    record("linear_limit", max(errs), 1e-3, "relative error of the knockdown at D = 1e-4")

    # 2. peak over mean excess against scripts/17
    errs = []
    for key, ax in axes.items():
        for sp, _ in ax["species"]:
            dx = ax["decay"][sp]
            for T in intervals:
                dose = 0.37
                e_peak = dose / (1.0 - math.exp(-dx * T))
                errs.append(abs(e_peak / (dose / (dx * T)) - dos.peak_to_average(dx, T)))
    record("peak_over_mean_vs_scripts_17", max(errs), 1e-12,
           "Epeak / Emean against scripts/17 peak_to_average")

    # 3. quasi-continuous dosing against the steady-state closed form
    errs = []
    for key, ax in axes.items():
        for rho_label, k_rep, rho in ax["rho_cases"]:
            sp = ax["species"][0][0]
            dx = ax["decay"][sp]
            kappa = ax["targets"][0][1]
            T = 0.05
            dose = required_excess(kappa, rho) * dx * T
            cas = Cascade(ax["d"], rho, ax["protein"][0][1], dx, lags=ax["lags"])
            errs.append(abs(1.0 - cas.cycle(dose, T)["mean_pi"] - kappa))
    record("continuous_limit", max(errs), 1e-3,
           "knockdown at T = 0.05 h with the closed-form mean excess")

    # 4. brute force over many cycles from the untreated state
    errs = []
    for key, ax in axes.items():
        sp = ax["species"][0][0]
        dx, rho = ax["decay"][sp], rho_of(ax, "scan high")
        dp = ax["protein"][0][1]
        kappa = ax["targets"][0][1]
        T = 24.0
        cas = Cascade(ax["d"], rho, dp, dx, lags=ax["lags"])
        dose, _ = cas.required_dose(kappa, T, required_excess(kappa, rho))
        r0, d = rho * ax["d"], ax["d"]

        def rhs(t, y):
            return [-dx * y[0], (d + r0) - (d + r0 * (1.0 + y[0])) * y[1], dp * (y[1] - y[2])]
        y = [0.0, 1.0, 1.0]
        for _ in range(120):
            y[0] += dose
            sol = solve_ivp(rhs, (0.0, T), y, method="LSODA", rtol=1e-11, atol=1e-13,
                            t_eval=np.linspace(0.0, T, 4001))
            y = list(sol.y[:, -1])
        last = sol.y[2]
        mean_pi = np.trapezoid(last, sol.t) / T
        errs.append(abs(1.0 - mean_pi - kappa))
    record("brute_force_cycles", max(errs), 1e-5,
           "mean knockdown of cycle 120 from the untreated state, trapezoid on 4001 points")

    # 5. constant exposure through scripts/12's own ODE (AD arm)
    errs = []
    kappa = axes["AD"]["targets"][0][1]
    for lab, k, rho in axes["AD"]["rho_cases"]:
        free = {kk: ode.midpoint(kk) for kk in ode.FREE}
        free["k_repress"] = k
        f = 1.0 + required_excess(kappa, rho)
        s1, _ = ode.run_ad(measured, "AD", free=free, aggregation=False)
        sf, _ = ode.run_ad(measured, "AD", free=free, mimic=f, aggregation=False)
        errs.append(abs(sf.y[4][-1] / s1.y[4][-1] - (1.0 - kappa)))
    record("script12_constant_exposure", max(errs), 1e-6,
           "BACE1 protein ratio at the closed-form mimic factor against 1 - kappa")

    # 6. same formula as scripts/12
    errs = []
    for lab, k, rho in axes["AD"]["rho_cases"]:
        errs.append(abs(ode.mimic_for_knockdown(axes["AD"]["d"], k, kappa)
                        - (1.0 + required_excess(kappa, rho))))
    record("same_formula_as_script12", max(errs), 1e-12,
           "scripts/12 mimic_for_knockdown against 1 + required_excess")
    return out


def nucleation_scaling(measured):
    """
    EN | Nucleation rate ratios from the measured reaction orders (see the module
         docstring). Every input is a measured constant of scripts/12.
    PT | Razoes de taxa de nucleacao a partir das ordens de reacao medidas (ver o
         docstring do modulo). Toda entrada e uma constante medida do scripts/12.
    """
    ratio = ((measured["k_prod_Ab_ad"] / measured["k_clear_Ab_ad"])
             / (measured["k_prod_Ab_ctrl"] / measured["k_clear_Ab_ctrl"]))
    kappa = 1.0 - 1.0 / ratio
    n1, n2 = measured["n_c"], measured["n_2"]
    return OrderedDict(
        monomer_ratio_AD_over_control=float(ratio),
        required_monomer_reduction=float(kappa),
        order_primary_nucleation=float(n1),
        order_secondary_nucleation=float(n2),
        primary_nucleation_rate_ratio_AD_over_control=float(ratio ** n1),
        secondary_nucleation_rate_ratio_AD_over_control=float(ratio ** n2),
        primary_nucleation_rate_reduction_after_offset=float(1.0 - (1.0 - kappa) ** n1),
        secondary_nucleation_rate_reduction_after_offset=float(1.0 - (1.0 - kappa) ** n2),
        elongation_rate_ratio_AD_over_control=float(ratio),
        note_en=("Ratios of rates at two monomer levels depend on the reaction orders only, "
                 "so they do not need the Abeta42 rate constants that no source gives (K037). "
                 "In vitro orders applied to a brain monomer pool: a scaling, not a rate."),
        note_pt=("Razoes de taxas em dois niveis de monomero dependem so das ordens de reacao, "
                 "entao nao precisam das constantes de velocidade do Abeta42 que nenhuma fonte "
                 "da (K037). Ordens in vitro aplicadas a um pool cerebral de monomero: um "
                 "escalonamento, nao uma taxa."))


def repression_strength(axes, k065):
    """
    EN | rho over the declared scan for each axis, the proportional limit, the mimic
         factor the closed form asks for at the ends of the scan and, for PD, at the
         K065 floor.
    PT | rho na varredura declarada para cada eixo, o limite proporcional, o fator de
         mimetico que a forma fechada pede nas pontas da varredura e, para DP, no valor
         piso do K065.
    """
    out = OrderedDict()
    for key, ax in axes.items():
        entry = OrderedDict(n_species=ax["n_sp"], d_mRNA_per_hour=ax["d"])
        scan = [(lab, k, r) for lab, k, r in ax["rho_cases"] if lab != "K065 lower bound"]
        entry["rho_scan_low"] = scan[0][2]
        entry["rho_scan_mid"] = scan[1][2]
        entry["rho_scan_high"] = scan[2][2]
        targets = OrderedDict()
        for lab, kappa in ax["targets"]:
            targets[lab] = OrderedDict(
                knockdown=kappa,
                proportional_limit_fold=1.0 / (1.0 - kappa),
                fold_at_scan_low=1.0 + required_excess(kappa, entry["rho_scan_low"]),
                fold_at_scan_mid=1.0 + required_excess(kappa, entry["rho_scan_mid"]),
                fold_at_scan_high=1.0 + required_excess(kappa, entry["rho_scan_high"]))
            if key == "PD":
                targets[lab]["fold_at_floor"] = 1.0 + required_excess(kappa, k065)
        entry["targets"] = targets
        out[key] = entry
    out["PD"]["rho_floor_K065"] = k065
    out["PD"]["floor_below_scan"] = bool(k065 < out["PD"]["rho_scan_low"])
    out["PD"]["note_en"] = (
        "K065: blocking endogenous miR-7 and miR-153 with antisense oligonucleotides raised an "
        "SNCA 3'UTR reporter by 44% in rat cortical neurons. For a repressor whose baseline decay "
        "of the mRNA is rho times the intrinsic decay, removing all of it raises the level by 1 + rho, and a blockade "
        "cannot be more than complete, so rho is at least 0.44 for the pair: a floor, not an "
        "estimate. Reporter, not endogenous gene; two miRNAs lumped, so nothing is said about "
        "miR-7 alone; the only loss-of-function observation in the table.")
    out["PD"]["note_pt"] = (
        "K065: bloquear miR-7 e miR-153 endogenos com oligonucleotideos antisense elevou em 44% "
        "um repórter da 3'UTR de SNCA em neuronios corticais de rato. Para um repressor cujo "
        "decaimento basal do mRNA vale rho vezes o decaimento intrinseco, retirar tudo dele eleva o nivel por 1 + rho, "
        "e um bloqueio nao pode passar de completo, entao rho e no minimo 0,44 para o par: um "
        "piso, nao uma estimativa. Repórter, nao gene endogeno; dois miRNAs agrupados, entao nada "
        "se diz do miR-7 sozinho; a unica observacao de perda de funcao da tabela.")
    return out


def rho_sweep(axes, interval=24.0):
    """
    EN | Every scenario at one dosing interval over a wide range of the repression
         strength, so that a statement about the dose does not have to rest on the
         declared scan. For each axis, delivered species and target the sweep keeps the
         smallest peak over rho, which is a bound that holds whatever rho is, the peak at
         the weakest and at the strongest repression tried, and the full rows. Weak
         repression makes the dose grow roughly as 1 / rho; strong repression approaches a
         limit, because the mRNA then tracks 1 / (1 + E) and rho drops out.
    PT | Todo cenario em um intervalo de dosagem sobre uma faixa ampla de forca de
         repressao, para que uma afirmacao sobre a dose nao dependa da varredura
         declarada. Para cada eixo, especie entregue e meta a varredura guarda o menor pico
         sobre rho, que e um limite valido seja qual for rho, o pico na repressao mais
         fraca e na mais forte testadas, e as linhas completas. A repressao fraca faz a
         dose crescer aproximadamente como 1 / rho; a forte se aproxima de um limite,
         porque o mRNA passa a seguir 1 / (1 + E) e rho some.
    """
    rows, summary = [], []
    for key, ax in axes.items():
        dp = ax["protein"][0][1]
        for sp, pid in ax["species"]:
            dx = ax["decay"][sp]
            for tgt_label, kappa in ax["targets"]:
                pts = []
                for rho in RHO_SWEEP:
                    rho = float(rho)
                    cas = Cascade(ax["d"], rho, dp, dx, lags=ax["lags"])
                    e_lin = required_excess(kappa, rho)
                    dose, reached = cas.required_dose(kappa, interval, e_lin)
                    r = OrderedDict(axis=key, species=sp, param_id=pid, half_life_h=ax["half_life"][sp],
                                    interval_h=interval, target=tgt_label, target_knockdown=kappa,
                                    rho=rho, reachable=dose is not None,
                                    linear_dose=e_lin * dx * interval)
                    if dose is None:
                        r.update(dose=None, peak_pool_over_baseline=None, efficiency=None,
                                 protein_ripple=None, knockdown_reached=reached)
                    else:
                        c = cas.cycle(dose, interval)
                        r.update(dose=dose, peak_pool_over_baseline=1.0 + c["e_peak"],
                                 efficiency=r["linear_dose"] / dose, protein_ripple=c["ripple_pi"],
                                 knockdown_reached=1.0 - c["mean_pi"])
                    pts.append(r)
                rows.extend(pts)
                ok = [r for r in pts if r["reachable"]]
                best = min(ok, key=lambda r: r["peak_pool_over_baseline"])
                summary.append(OrderedDict(
                    axis=key, species=sp, target=tgt_label, target_knockdown=kappa,
                    interval_h=interval,
                    min_peak_pool_over_baseline=best["peak_pool_over_baseline"],
                    rho_at_min=best["rho"],
                    peak_at_lowest_rho=pts[0]["peak_pool_over_baseline"],
                    peak_at_highest_rho=pts[-1]["peak_pool_over_baseline"],
                    unreachable_points=sum(1 for r in pts if not r["reachable"])))
        print(f"  {key}: rho sweep done ({len([r for r in rows if r['axis'] == key])} points)")
    return rows, summary


SWEEP_FIELDS = ["axis", "species", "param_id", "half_life_h", "interval_h", "target", "target_knockdown",
                "rho", "reachable", "knockdown_reached", "dose", "peak_pool_over_baseline",
                "linear_dose", "efficiency", "protein_ripple"]


FIELDS = ["axis", "rho_case", "k_repress", "rho", "species", "param_id", "half_life_h",
          "interval_h", "protein_case", "protein_half_life_h", "target", "target_knockdown",
          "reachable", "knockdown_reached", "dose", "peak_excess", "peak_pool_over_baseline",
          "trough_pool_over_baseline", "mean_excess", "linear_excess", "linear_dose",
          "efficiency", "protein_ripple", "abeta_ripple", "elongation_reduction_m49_8",
          "elongation_reduction_m46_0", "residual"]


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    ode = load_module("ode12", "12_ode_models_calibrated.py")
    dos = load_module("dosing17", "17_mimic_dosing_feasibility.py")
    measured, table = ode.load_measured()
    axes, kappa_ad, k065 = build_axes(ode, dos, measured, table)
    intervals = list(dos.INTERVALS_H)

    print("=" * 78)
    print("EN | Mimic target response | PT | Resposta do alvo ao mimetico")
    print("=" * 78)

    print("\nVerification | Verificacao")
    checks = verify(axes, ode, dos, measured, intervals)
    for name, c in checks.items():
        print(f"  {name:32s} error {c['max_abs_error']:.2e}  (tolerance {c['tolerance']:.0e})  ok")

    print("\nScenarios | Cenarios")
    rows = run_scenarios(axes, intervals)

    strength = repression_strength(axes, k065)
    nuc = nucleation_scaling(measured)
    print("\nRepression-strength sweep at 24 h | Varredura da forca de repressao em 24 h")
    sweep_rows, sweep_summary = rho_sweep(axes)

    print("\n--- Repression strength rho | Forca de repressao rho ---")
    for key in ("AD", "PD"):
        e = strength[key]
        print(f"  {key}: rho {e['rho_scan_low']:.2f} .. {e['rho_scan_high']:.1f} over the scan "
              f"(midpoint {e['rho_scan_mid']:.2f})")
    print(f"  PD, K065 floor on the pair: rho >= {k065:.2f}  (floor below the scan: "
          f"{strength['PD']['floor_below_scan']})")
    for lab, t in strength["PD"]["targets"].items():
        print(f"  PD {lab}: proportional limit {t['proportional_limit_fold']:.2f}x, "
              f"scan {t['fold_at_scan_high']:.2f}..{t['fold_at_scan_low']:.2f}x, "
              f"at the floor {t['fold_at_floor']:.2f}x")

    print("\n--- Abeta42 nucleation scaling | Escalonamento da nucleacao do Abeta42 ---")
    print(f"  monomer AD/control {nuc['monomer_ratio_AD_over_control']:.4f}; nucleation rate "
          f"ratio {nuc['primary_nucleation_rate_ratio_AD_over_control']:.3f} (primary), "
          f"{nuc['secondary_nucleation_rate_ratio_AD_over_control']:.3f} (secondary)")
    print(f"  offsetting the clearance deficit lowers nucleation by "
          f"{100 * nuc['primary_nucleation_rate_reduction_after_offset']:.1f}% (primary), "
          f"{100 * nuc['secondary_nucleation_rate_reduction_after_offset']:.1f}% (secondary)")

    print("\n--- Daily dosing (24 h), primary protein half-life | Dosagem diaria ---")
    for r in rows:
        if r["interval_h"] != 24.0 or not r["protein_case"].startswith(("BACE1 in", "alpha-syn")):
            continue
        if r["axis"] == "PD" and r["rho_case"] not in ("K065 lower bound", "scan midpoint"):
            continue
        if r["axis"] == "AD" and r["rho_case"] != "scan midpoint":
            continue
        if r["dose"] is None:
            print(f"  {r['axis']} {r['rho_case']:22s} {r['species']:18s} kappa {r['target_knockdown']:.3f}"
                  f"  unreachable (reaches {100 * r['knockdown_reached']:.1f}%)")
            continue
        print(f"  {r['axis']} {r['rho_case']:22s} {r['species']:18s} kappa {r['target_knockdown']:.3f}"
              f"  dose {r['dose']:8.3f}  peak pool {r['peak_pool_over_baseline']:7.2f}x"
              f"  efficiency {r['efficiency']:.3f}  protein swing {100 * r['protein_ripple']:5.1f}%")

    print("\n--- Smallest daily peak over all rho (x baseline pool) | Menor pico diario sobre todo rho ---")
    for e in sweep_summary:
        print(f"  {e['axis']} {e['species']:18s} kappa {e['target_knockdown']:.3f}: "
              f"min {e['min_peak_pool_over_baseline']:8.2f} at rho {e['rho_at_min']:9.3g}   "
              f"weakest {e['peak_at_lowest_rho']:9.1f}   strongest {e['peak_at_highest_rho']:8.2f}")
    with open(f"{TAB_DIR}/mimic_target_rho_sweep.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SWEEP_FIELDS)
        w.writeheader()
        for r in sweep_rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in SWEEP_FIELDS})
    with open(f"{TAB_DIR}/mimic_target_response.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in FIELDS})

    payload = OrderedDict([
        ("generated_on", "2026-10-14"),
        ("software", OrderedDict(python=platform.python_version(), numpy=np.__version__,
                                  scipy=scipy.__version__)),
        ("integration", "LSODA via scipy.integrate.solve_ivp, rtol 1e-10, atol 1e-13; periodic "
                        "state from the affine one-cycle map; cycle averages integrated as states"),
        ("dose_cap_in_baseline_pools", DOSE_CAP),
        ("intervals_h", intervals),
        ("scenarios", len(rows)),
        ("scenarios_unreachable", sum(1 for r in rows if not r["reachable"])),
        ("repression_strength", strength),
        ("rho_sweep", OrderedDict(interval_h=24.0, rho_min=float(RHO_SWEEP[0]), rho_max=float(RHO_SWEEP[-1]),
                                   points=len(RHO_SWEEP), summary=sweep_summary)),
        ("nucleation_scaling", nuc),
        ("verification", checks),
        ("assumptions_en",
         "The delivered mimic is cleared at the first-order rate of the endogenous species "
         "(K016-K021) and enters the repressor pool instantly. Repression is mass action: linear "
         "in the miRNA, saturating only through depletion of the mRNA. The delivered species acts "
         "like the endogenous pool, molecule for molecule. The three miR-29 paralogues start at "
         "equal levels and are pooled; miR-153 is lumped with miR-7. k_translate cancels, and "
         "k_repress enters through rho alone."),
        ("assumptions_pt",
         "O mimetico entregue e eliminado na taxa de primeira ordem da especie endogena "
         "(K016-K021) e entra no pool do repressor instantaneamente. A repressao e acao de "
         "massas: linear no miRNA, saturando so pelo esgotamento do mRNA. A especie entregue age "
         "como o pool endogeno, molecula por molecula. Os tres paralogos de miR-29 comecam em "
         "niveis iguais e sao agrupados; o miR-153 e agrupado com o miR-7. O k_translate se "
         "cancela, e o k_repress entra so por rho."),
    ])
    with open(f"{TAB_DIR}/mimic_target_response.json", "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)

    print(f"\nEN/PT -> {TAB_DIR}/mimic_target_rho_sweep.csv  ({len(sweep_rows)} rows)")
    print(f"EN/PT -> {TAB_DIR}/mimic_target_response.csv  ({len(rows)} rows)")
    print(f"EN/PT -> {TAB_DIR}/mimic_target_response.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
