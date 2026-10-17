#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Three animations of mechanisms this project already simulated, not new ones. Each
     reuses the exact ODE right-hand sides and measured constants from scripts/12, or the
     closed-form dosing formula from scripts/17, imported directly rather than re-derived,
     so an animation cannot silently disagree with the static result it is a movie of.
PT | Tres animacoes de mecanismos que este projeto ja simulou, nao mecanismos novos. Cada
     uma reaproveita exatamente os lados-direitos de EDO e as constantes medidas do
     scripts/12, ou a formula fechada de dosagem do scripts/17, importados diretamente em
     vez de rederivados, para que uma animacao nao possa discordar em silencio do
     resultado estatico do qual e o filme.

EN | Where the line is, and why it is drawn here. "Illustrate the mechanism with a
     computational simulation" can mean two different things. One is: animate a
     simulation that was actually run - a numerical integration whose equations,
     constants and results are already in this repository and already checked by
     scripts/08. The other is: produce something that looks like a simulation - a
     smooth interpolation between two states with no physics connecting them, such as a
     protein "folding" between two crystal structures, or a fibril "growing" atom by atom
     with no force field behind it. The first is this script. The second would be
     fabrication with better production values, and this project's rule against forging
     data does not have an animation exception. Every curve drawn here is therefore
     restricted to trajectories with either a measured decay constant (the mimic washout
     and dosing panels: no free parameter enters an exponential decay) or a monomer pool
     with no free aggregation constant (the Abeta42 accumulation panel: k_prod and
     k_clear are measured, and the AD/control ratio does not depend on the illustrative
     aggregation block at all). The alpha-synuclein pH-gate trajectory is deliberately
     NOT animated here, because its rate constant at acidic pH is declared illustrative
     (K042, qualitative only) - animating it would put a specific speed on screen for a
     quantity this project explicitly refuses to report as a number, only as a direction.
PT | Onde fica a linha, e por que ela e tracada aqui. "Ilustrar o mecanismo com uma
     simulacao computacional" pode significar duas coisas diferentes. Uma e: animar uma
     simulacao que de fato rodou - uma integracao numerica cujas equacoes, constantes e
     resultados ja estao neste repositorio e ja sao conferidos pelo scripts/08. A outra e:
     produzir algo que parece uma simulacao - uma interpolacao suave entre dois estados
     sem fisica nenhuma ligando-os, como uma proteina "se dobrando" entre duas estruturas
     cristalograficas, ou uma fibrila "crescendo" atomo por atomo sem campo de forca por
     tras. A primeira e este script. A segunda seria fabricacao com producao melhor, e a
     regra deste projeto contra forjar dados nao tem excecao para animacao. Toda curva
     aqui desenhada fica por isso restrita a trajetorias com constante de decaimento
     medida (os paineis de eliminacao e dosagem: nenhum parametro livre entra numa
     exponencial) ou a um pool de monomero sem constante de agregacao livre (o painel de
     acumulo de Abeta42: k_prod e k_clear sao medidos, e a razao AD/controle nao depende
     em nada do bloco de agregacao ilustrativo). A trajetoria da comporta de pH da
     alfa-sinucleina deliberadamente NAO e animada aqui, porque sua constante de
     velocidade em pH acido e declarada ilustrativa (K042, so qualitativa) - anima-la
     poria uma velocidade especifica na tela para uma grandeza que este projeto
     explicitamente se recusa a reportar como numero, so como direcao.

EN | Every number that appears on a frame - the final AD/control ratio, the washout
     times, the peak-to-average penalty - is read at run time from scripts/12's and
     scripts/17's own functions and from results/tables/*.json, and a manifest of those
     numbers is written to results/tables/mechanism_animations.json so scripts/08 can
     check the animation's endpoint against the already-verified result, the same
     discipline used for every other figure here.
PT | Todo numero que aparece num quadro - a razao final AD/controle, os tempos de
     eliminacao, a penalidade de pico sobre media - e lido em tempo de execucao das
     proprias funcoes do scripts/12 e do scripts/17 e de results/tables/*.json, e um
     manifesto desses numeros e escrito em results/tables/mechanism_animations.json para
     que o scripts/08 confira o final da animacao contra o resultado ja verificado, a
     mesma disciplina usada em toda outra figura daqui.

    python scripts/19_mechanism_animations.py
"""

import importlib.util
import json
import math
import os
import sys
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bilingual import LANGS, t as tr, fig_path

HERE = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = "results/tables"
FIG_DIR = "results/figures"
AD_HUE, PD_HUE, REF_HUE = "#2a78d6", "#eb6834", "#52514e"
LN2 = math.log(2.0)


def load_module(name, filename):
    """
    EN | Import a sibling script whose filename starts with a digit (12_..., 17_...) and
         so cannot be named in a normal `import` statement. This is the only way to reuse
         its functions rather than copy them, which is the point: two copies of the same
         ODE would be two places for them to quietly drift apart.
    PT | Importa um script irmao cujo nome comeca com digito (12_..., 17_...) e por isso
         nao pode ser nomeado num `import` comum. E o unico jeito de reaproveitar suas
         funcoes em vez de copia-las, que e o propósito: duas copias da mesma EDO seriam
         dois lugares onde elas poderiam divergir em silencio.
    """
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        sys.exit(f"EN/PT: missing {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_manifest(ode, dosing):
    """
    EN | Run each reused simulation once, here, and record the numbers every animation
         will draw and every scripts/08 check will compare against. Drawing happens later
         and only formats these numbers as text; it does not compute anything new.
    PT | Roda cada simulacao reaproveitada uma vez, aqui, e registra os numeros que toda
         animacao vai desenhar e que todo check do scripts/08 vai comparar. O desenho
         acontece depois e so formata estes numeros como texto; nao calcula nada novo.
    """
    measured, _kin_table = ode.load_measured()

    sc, _ = ode.run_ad(measured, "control")
    sd, _ = ode.run_ad(measured, "AD")
    ratio_final = float(sd.y[5][-1] / sc.y[5][-1])

    washout = ode.experiment_mimic_washout(measured)
    fast_label = "miR-7 (1.7 h, K020)"
    slow_label = "median miRNA (34 h, K016)"

    dec = dosing.load_decays()
    d_fast = dec["miR-7"]["decay_per_hour"]
    d_slow = dec["median miRNA"]["decay_per_hour"]
    T = 24.0
    peak_fast = dosing.peak_to_average(d_fast, T)
    peak_slow = dosing.peak_to_average(d_slow, T)

    manifest = OrderedDict([
        ("ad_monomer_accumulation", dict(
            t_end_hours=float(sc.t[-1]), ratio_AD_over_control=ratio_final,
            source="ode_calibrated_results.json:ad_clearance_vs_production")),
        ("mimic_washout", dict(
            fast_species=fast_label, slow_species=slow_label,
            fast_hours_to_baseline=washout[fast_label]["hours_to_return_to_baseline"],
            slow_hours_to_baseline=washout[slow_label]["hours_to_return_to_baseline"],
            source="ode_calibrated_results.json:mimic_washout")),
        ("dosing_sawtooth", dict(
            dosing_interval_hours=T, fast_species="miR-7", slow_species="median miRNA",
            peak_over_average_fast=peak_fast, peak_over_average_slow=peak_slow,
            penalty_ratio=peak_fast / peak_slow,
            source="mimic_dosing_feasibility.json:daily_dosing_comparison")),
    ])
    return measured, manifest, dict(sc=sc, sd=sd, washout=washout,
                                    d_fast=d_fast, d_slow=d_slow, T=T,
                                    peak_fast=peak_fast, peak_slow=peak_slow)


# ---------------------------------------------------------------------------
# 1. AD axis: Abeta42 monomer accumulation, control vs AD - no free parameter
# ---------------------------------------------------------------------------
def animate_ad_monomer(sim, manifest, path, lang, n_frames=110, fps=18):
    sc, sd = sim["sc"], sim["sd"]
    t_days = sc.t / 24.0
    idx = np.linspace(0, len(t_days) - 1, n_frames).astype(int)
    y_max = 1.08 * max(sc.y[5].max(), sd.y[5].max())

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    line_c, = ax.plot([], [], color=AD_HUE, lw=2.4, alpha=0.55,
                      label=tr(lang, "control", "controle"))
    line_d, = ax.plot([], [], color=AD_HUE, lw=2.4,
                      label=tr(lang, "AD (clearance -30%)", "DA (depuração -30%)"))
    marker, = ax.plot([], [], "o", color=AD_HUE, ms=6)
    txt = ax.text(0.97, 0.06, "", transform=ax.transAxes, ha="right", fontsize=10.5,
                  color=AD_HUE, weight="bold")
    ax.set_xlim(0, t_days[-1])
    ax.set_ylim(0, y_max)
    ax.set_xlabel(tr(lang, "days", "dias"))
    ax.set_ylabel(tr(lang, "Aβ42 monomer (relative to baseline)",
                     "monômero de Aβ42 (relativo ao basal)"))
    ax.set_title(tr(lang, "Aβ42 reaches a new steady state - no free parameter",
                    "Aβ42 alcança um novo estado estacionário - sem parâmetro livre"),
                fontsize=11)
    ax.legend(loc="center right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    ratio = manifest["ad_monomer_accumulation"]["ratio_AD_over_control"]

    def update(frame):
        i = idx[frame]
        line_c.set_data(t_days[:i + 1], sc.y[5][:i + 1])
        line_d.set_data(t_days[:i + 1], sd.y[5][:i + 1])
        marker.set_data([t_days[i]], [sd.y[5][i]])
        if frame >= n_frames - 8:
            txt.set_text(tr(lang, f"AD / control = {ratio:.2f}×",
                            f"DA / controle = {ratio:.2f}×"))
        return line_c, line_d, marker, txt

    ani = FuncAnimation(fig, update, frames=n_frames, blit=False)
    ani.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. Single-dose mimic washout: measured half-lives, no simulation needed beyond
#    the same first-order decay scripts/12 already reports.
# ---------------------------------------------------------------------------
def animate_mimic_washout(sim, manifest, path, lang, n_frames=100, fps=18):
    m = manifest["mimic_washout"]
    washout = sim["washout"]
    d_fast = washout[m["fast_species"]]["decay_constant_per_hour"]
    d_slow = washout[m["slow_species"]]["decay_constant_per_hour"]
    t_end = 1.05 * washout[m["slow_species"]]["hours_to_return_to_baseline"]
    t = np.linspace(0, t_end, 600)
    y_fast = 1.0 + 9.0 * np.exp(-d_fast * t)
    y_slow = 1.0 + 9.0 * np.exp(-d_slow * t)
    idx = np.linspace(0, len(t) - 1, n_frames).astype(int)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axhline(1.1, color=REF_HUE, ls="--", lw=1.0)
    ax.text(t_end * 0.99, 1.15, tr(lang, "within 10% of baseline", "a 10% do basal"),
            ha="right", fontsize=8, color=REF_HUE)
    line_f, = ax.plot([], [], color=PD_HUE, lw=2.4,
                      label=tr(lang, "miR-7 (1.7 h)", "miR-7 (1,7 h)"))
    line_s, = ax.plot([], [], color=REF_HUE, lw=2.4,
                      label=tr(lang, "median miRNA (34 h)", "miRNA mediano (34 h)"))
    ax.set_yscale("log")
    ax.set_xlim(0, t_end)
    ax.set_ylim(0.9, 11)
    ax.set_xlabel(tr(lang, "hours after a 10-fold bolus", "horas após um bolus de 10×"))
    ax.set_ylabel(tr(lang, "fold over baseline", "vezes o basal"))
    ax.set_title(tr(lang, "A single mimic dose washes out at its own measured half-life",
                    "Uma dose única de mimético se elimina na própria meia-vida medida"),
                fontsize=11)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    def update(frame):
        i = idx[frame]
        line_f.set_data(t[:i + 1], y_fast[:i + 1])
        line_s.set_data(t[:i + 1], y_slow[:i + 1])
        return line_f, line_s

    ani = FuncAnimation(fig, update, frames=n_frames, blit=False)
    ani.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3. Repeated daily dosing: the peak-to-average penalty as a moving sawtooth.
#    Closed form, dose-free (scripts/17): C(tau)/average = d*T*exp(-d*tau)/(1-exp(-dT)).
# ---------------------------------------------------------------------------
def sawtooth(d, T, n_cycles, points_per_cycle=120):
    """EN/PT: concentration relative to ITS OWN average, so both species plot with
    average = 1 and the peak height read directly off the axis IS peak_to_average(d,T)."""
    tau = np.linspace(0, T, points_per_cycle, endpoint=False)
    one_cycle = d * T * np.exp(-d * tau) / (1.0 - math.exp(-d * T))
    t = np.concatenate([tau + k * T for k in range(n_cycles)])
    y = np.tile(one_cycle, n_cycles)
    return t, y


def animate_dosing_sawtooth(sim, manifest, path, lang, n_frames=130, fps=20, n_cycles=4):
    d = manifest["dosing_sawtooth"]
    T = d["dosing_interval_hours"]
    t_fast, y_fast = sawtooth(sim["d_fast"], T, n_cycles)
    t_slow, y_slow = sawtooth(sim["d_slow"], T, n_cycles)
    idx = np.linspace(0, len(t_fast) - 1, n_frames).astype(int)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for k in range(n_cycles + 1):
        ax.axvline(k * T, color="#cccccc", lw=0.8, zorder=0)
    line_f, = ax.plot([], [], color=PD_HUE, lw=2.0, label=tr(lang, "miR-7 mimic", "mimético de miR-7"))
    line_s, = ax.plot([], [], color=REF_HUE, lw=2.0,
                      label=tr(lang, "median-stability mimic", "mimético de estabilidade mediana"))
    txt = ax.text(0.02, 0.94, "", transform=ax.transAxes, fontsize=10, color=PD_HUE,
                  weight="bold", va="top")
    ax.set_xlim(0, n_cycles * T)
    ax.set_ylim(0, 1.08 * d["peak_over_average_fast"])
    ax.set_xlabel(tr(lang, "hours (dosed once daily)", "horas (dose diária)"))
    ax.set_ylabel(tr(lang, "concentration ÷ its own average", "concentração ÷ sua própria média"))
    ax.set_title(tr(lang, "The peak a fast-decaying mimic must reach to hold the same average",
                    "O pico que um mimético de decaimento rápido precisa atingir para "
                    "sustentar a mesma média"), fontsize=10.5)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    def update(frame):
        i = idx[frame]
        line_f.set_data(t_fast[:i + 1], y_fast[:i + 1])
        line_s.set_data(t_slow[:i + 1], y_slow[:i + 1])
        if frame >= n_frames - 10:
            txt.set_text(tr(lang, f"{d['penalty_ratio']:.1f}× penalty",
                            f"penalidade de {d['penalty_ratio']:.1f}×"))
        return line_f, line_s, txt

    ani = FuncAnimation(fig, update, frames=n_frames, blit=False)
    ani.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    ode = load_module("ode12", "12_ode_models_calibrated.py")
    dosing = load_module("dosing17", "17_mimic_dosing_feasibility.py")
    measured, manifest, sim = build_manifest(ode, dosing)

    with open(f"{TAB_DIR}/mechanism_animations.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=1)

    for lang in LANGS:
        animate_ad_monomer(sim, manifest, fig_path(FIG_DIR, "anim_ad_monomer", lang, "gif"), lang)
        animate_mimic_washout(sim, manifest, fig_path(FIG_DIR, "anim_mimic_washout", lang, "gif"), lang)
        animate_dosing_sawtooth(sim, manifest, fig_path(FIG_DIR, "anim_dosing_sawtooth", lang, "gif"), lang)

    print("=" * 78)
    print("EN | Mechanism animations | PT | Animações de mecanismo")
    print("=" * 78)
    print(f"  AD monomer ratio (final frame)   : {manifest['ad_monomer_accumulation']['ratio_AD_over_control']:.3f}")
    print(f"  Dosing penalty (final frame)     : {manifest['dosing_sawtooth']['penalty_ratio']:.2f}x")
    print(f"\nEN/PT -> {TAB_DIR}/mechanism_animations.json")
    for stem in ("anim_ad_monomer", "anim_mimic_washout", "anim_dosing_sawtooth"):
        for lang in LANGS:
            print(f"EN/PT -> {fig_path(FIG_DIR, stem, lang, 'gif')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
