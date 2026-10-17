#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | The design space of a microRNA mimic: the peak of the repressor pool, in units of its
     baseline pool, that a repeated dose has to reach to bring the cycle-averaged protein
     down by the target, over a grid of repression strengths (rho) and dosing intervals,
     for every delivered species and target of the mechanistic study.

     scripts/36 solved the scenarios the article quotes (five intervals at four values of
     rho, and one interval at 31 values of rho). This script fills the plane between them
     with the same solver, so that the figure can show where in (rho, interval) a regimen
     is within reach and where it is not, instead of two slices through it. Nothing new is
     assumed: the cascade, the constants, the species and the targets are the ones of
     scripts/36, built by its own build_axes; the protein half-life does not enter the
     dose (the mean protein equals the mean mRNA over a periodic cycle), so one protein
     half-life is used and the lagged Abeta42 state is left out, which cannot change the
     dose either.

     What the script checks before it writes anything:
       1. at 24 h, at the values of rho that the grid shares with the sweep of scripts/36
          (one per decade), the peak equals the sweep's;
       2. along the interval axis the required peak never falls as the interval grows, and
          a point that no dose up to the cap reaches stays unreached at longer intervals;
       3. the shortest intervals approach the closed form of the steady state, as they must
          when the dosing is nearly continuous.
     It then reads off, for a few thresholds of the peak, the longest interval that keeps the
     required peak within the threshold. The thresholds are descriptive (3, 10 and 30 times
     the baseline pool), not a claim about what a formulation can reach.

PT | O espaco de projeto de um mimetico de microRNA: o pico do pool do repressor, em unidades
     do pool basal, que uma dose repetida precisa atingir para levar a proteina media no
     ciclo a queda-alvo, numa grade de forcas de repressao (rho) e intervalos de dosagem,
     para toda especie entregue e meta do estudo mecanistico.

     O scripts/36 resolveu os cenarios que o artigo cita (cinco intervalos em quatro valores
     de rho, e um intervalo em 31 valores de rho). Este script preenche o plano entre eles
     com o mesmo resolvedor, para que a figura mostre onde em (rho, intervalo) um regime esta
     ao alcance e onde nao, em vez de duas fatias dele. Nada novo e suposto: a cascata, as
     constantes, as especies e as metas sao as do scripts/36, montadas pelo seu proprio
     build_axes; a meia-vida da proteina nao entra na dose (a proteina media iguala o mRNA
     medio num ciclo periodico), entao uma meia-vida de proteina e usada e o estado
     defasado do Abeta42 fica de fora, o que tambem nao pode mudar a dose.

     O que o script confere antes de gravar qualquer coisa:
       1. em 24 h, nos valores de rho que a grade divide com a varredura do scripts/36 (um por
          decada), o pico e igual ao da varredura;
       2. no eixo do intervalo o pico exigido nunca cai quando o intervalo cresce, e um ponto
          que nenhuma dose ate o teto alcanca continua sem alcance em intervalos maiores;
       3. os intervalos mais curtos se aproximam da forma fechada do estado estacionario,
          como devem quando a dosagem e quase continua.
     Depois le, para alguns limiares do pico, o maior intervalo que mantem o pico exigido
     dentro do limiar. Os limiares sao descritivos (3, 10 e 30 vezes o pool basal), nao uma
     afirmacao sobre o que uma formulacao alcanca.

    python scripts/40_mimic_design_space.py
"""

import csv
import importlib.util
import json
import math
import multiprocessing as mp
import os
import sys
import time
from collections import OrderedDict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TAB_DIR = "results/tables"
SWEEP_CSV = f"{TAB_DIR}/mimic_target_rho_sweep.csv"
OUT_CSV = f"{TAB_DIR}/mimic_design_space.csv"
OUT_JSON = f"{TAB_DIR}/mimic_design_space.json"

RHO_GRID = np.geomspace(0.01, 1.0e4, 49)                 # EN/PT: eight points per decade
T_GRID = [2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0, 30.0, 36.0, 48.0, 60.0, 72.0, 96.0,
          120.0, 144.0, 168.0]                           # EN/PT: hours; holds the five intervals of scripts/17 and 36
THRESHOLDS = (3.0, 10.0, 30.0)                           # EN/PT: peak of the repressor pool over baseline

_G = {}                                                  # EN/PT: shared with the worker processes (fork)


def load_module(name, filename):
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        sys.exit(f"EN/PT: missing {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def series_list(axes):
    """EN/PT: every (axis, species, target) of scripts/36, in its own order."""
    out = []
    for key, ax in axes.items():
        for sp, pid in ax["species"]:
            for tgt_label, kappa in ax["targets"]:
                out.append(OrderedDict(axis=key, species=sp, param_id=pid, target=tgt_label,
                                       target_knockdown=float(kappa), half_life_h=float(ax["half_life"][sp]),
                                       decay_per_hour=float(ax["decay"][sp])))
    return out


def solve_row(task):
    """EN/PT: one (series, rho) at every interval | uma (serie, rho) em todos os intervalos."""
    si, ri = task
    m36, axes, series = _G["m36"], _G["axes"], _G["series"]
    s = series[si]
    ax = axes[s["axis"]]
    rho = float(RHO_GRID[ri])
    dp = ax["protein"][0][1]
    kappa = s["target_knockdown"]
    out = []
    for T in T_GRID:
        cas = m36.Cascade(ax["d"], rho, dp, s["decay_per_hour"], lags=())
        e_lin = m36.required_excess(kappa, rho)
        dose, reached = cas.required_dose(kappa, T, e_lin)
        if dose is None:
            out.append((si, ri, T, False, None, reached))
        else:
            c = cas.cycle(dose, T)
            out.append((si, ri, T, True, 1.0 + c["e_peak"], 1.0 - c["mean_pi"]))
    return out


def interval_for_peak(Ts, peaks, level):
    """
    EN | The longest interval at which the required peak is within `level`, by linear
         interpolation of log(peak) against log(interval) between grid points. Returns
         (hours, flag) with flag "none" when even the shortest interval of the grid needs
         more than `level`, "beyond" when the longest one does not, and "interp" otherwise.
         Unreached points count as above any level.
    PT | O maior intervalo em que o pico exigido fica dentro de `level`, por interpolacao
         linear de log(pico) contra log(intervalo) entre pontos da grade. Devolve (horas,
         marca): "none" quando ate o menor intervalo da grade exige mais que `level`,
         "beyond" quando o maior nao exige, "interp" nos demais casos. Pontos sem alcance
         contam como acima de qualquer nivel.
    """
    vals = [math.inf if p is None else p for p in peaks]
    if vals[0] > level:
        return 0.0, "none"
    if vals[-1] <= level:
        return Ts[-1], "beyond"
    for i in range(len(Ts) - 1):
        if vals[i] <= level < vals[i + 1]:
            if math.isinf(vals[i + 1]):
                # EN/PT: cross into the unreached region: no finite end to interpolate to
                return Ts[i], "edge"
            f = (math.log(level) - math.log(vals[i])) / (math.log(vals[i + 1]) - math.log(vals[i]))
            return float(math.exp(math.log(Ts[i]) + f * (math.log(Ts[i + 1]) - math.log(Ts[i])))), "interp"
    return Ts[-1], "beyond"


def main():
    os.makedirs(TAB_DIR, exist_ok=True)
    m36 = load_module("m36", "36_mimic_target_response.py")
    ode = m36.load_module("ode12", "12_ode_models_calibrated.py")
    dos = m36.load_module("dosing17", "17_mimic_dosing_feasibility.py")
    measured, table = ode.load_measured()
    axes, kappa_ad, k065 = m36.build_axes(ode, dos, measured, table)
    series = series_list(axes)
    _G.update(m36=m36, axes=axes, series=series)

    print("=" * 78)
    print("EN | Mimic design space | PT | Espaco de projeto do mimetico")
    print("=" * 78)
    n_pts = len(series) * len(RHO_GRID) * len(T_GRID)
    print(f"  {len(series)} series x {len(RHO_GRID)} values of rho x {len(T_GRID)} intervals = {n_pts} points")

    t0 = time.time()
    tasks = [(si, ri) for si in range(len(series)) for ri in range(len(RHO_GRID))]
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=max(1, (os.cpu_count() or 2))) as pool:
        chunks = pool.map(solve_row, tasks, chunksize=2)
    rows = [r for ch in chunks for r in ch]
    print(f"  solved in {time.time() - t0:.0f} s")

    # EN/PT: arrays [series, rho, T]
    peak = np.full((len(series), len(RHO_GRID), len(T_GRID)), np.nan)
    reach = np.zeros_like(peak, dtype=bool)
    kd = np.zeros_like(peak)
    for si, ri, T, ok, pk, k in rows:
        ti = T_GRID.index(T)
        reach[si, ri, ti] = ok
        kd[si, ri, ti] = k
        if ok:
            peak[si, ri, ti] = pk

    # ------------------------------------------------------------------ checks
    checks = OrderedDict()
    # 1. against the sweep of scripts/36 at 24 h, at the values of rho the grids share
    sweep = list(csv.DictReader(open(SWEEP_CSV, encoding="utf-8")))
    t24 = T_GRID.index(24.0)
    worst, n_cmp = 0.0, 0
    for si, s in enumerate(series):
        for r in sweep:
            if (r["axis"] != s["axis"] or r["species"] != s["species"]
                    or abs(float(r["target_knockdown"]) - s["target_knockdown"]) > 1e-9):
                continue
            rho = float(r["rho"])
            ri = int(np.argmin(np.abs(np.log10(RHO_GRID) - math.log10(rho))))
            if abs(math.log10(RHO_GRID[ri]) - math.log10(rho)) > 1e-9:
                continue
            n_cmp += 1
            if (r["reachable"] == "True") != bool(reach[si, ri, t24]):
                sys.exit(f"EN/PT: reachability differs from the sweep for {s['species']} at rho {rho:.3g}")
            if r["reachable"] == "True":
                worst = max(worst, abs(float(r["peak_pool_over_baseline"]) - peak[si, ri, t24]) /
                            float(r["peak_pool_over_baseline"]))
    if n_cmp < len(series) * 7 or worst > 1e-6:
        sys.exit(f"EN/PT: the grid disagrees with the sweep of scripts/36 (worst relative difference "
                 f"{worst:.2e} over {n_cmp} points)")
    checks["against_sweep_36_at_24h"] = dict(points=n_cmp, max_relative_error=worst, tolerance=1e-6)
    # 2. monotone in the interval; unreached stays unreached
    worst_drop, n_viol = 0.0, 0
    for si in range(len(series)):
        for ri in range(len(RHO_GRID)):
            seen_unreached = False
            prev = None
            for ti in range(len(T_GRID)):
                if not reach[si, ri, ti]:
                    seen_unreached = True
                    continue
                if seen_unreached:
                    n_viol += 1
                if prev is not None:
                    worst_drop = max(worst_drop, (prev - peak[si, ri, ti]) / prev)
                prev = peak[si, ri, ti]
    if n_viol or worst_drop > 1e-7:
        sys.exit(f"EN/PT: the required peak falls with the interval somewhere (worst relative fall "
                 f"{worst_drop:.2e}; {n_viol} reached points after an unreached one)")
    checks["monotone_in_interval"] = dict(max_relative_fall=worst_drop, reached_after_unreached=n_viol)
    # 3. shortest interval against the closed form: peak -> 1 + excess when dosing is continuous
    worst_cf = 0.0
    for si, s in enumerate(series):
        for ri, rho in enumerate(RHO_GRID):
            if not reach[si, ri, 0]:
                continue
            e_lin = m36.required_excess(s["target_knockdown"], float(rho))
            # EN/PT: at 2 h the peak is above the closed-form pool by the peak-to-average factor of the delivered species
            x = s["decay_per_hour"] * T_GRID[0]
            r_pa = x / (1.0 - math.exp(-x))
            approx = 1.0 + e_lin * r_pa
            worst_cf = max(worst_cf, abs(peak[si, ri, 0] - approx) / approx)
    # EN/PT: not an equality (the mRNA and protein filter even a 2 h cycle a little); a loose, stated bound
    if worst_cf > 0.35:
        sys.exit(f"EN/PT: at 2 h the peak is {100 * worst_cf:.0f}% from the closed-form value; the grid "
                 "is not approaching the steady-state limit")
    checks["shortest_interval_vs_closed_form"] = dict(max_relative_difference=worst_cf, bound=0.35,
                                                      note="closed form: 1 + excess x peak-to-average factor of the delivered species")

    # ------------------------------------------------------------------ the longest interval within a threshold
    strength = ode.midpoint("k_repress")
    mids = {}
    for key, ax in axes.items():
        d = ax["d"]
        k_lo, k_hi = ode.FREE["k_repress"]["low"], ode.FREE["k_repress"]["high"]
        n = ax["n_sp"]
        mids[key] = dict(low=n * k_lo / d, mid=n * strength / d, high=n * k_hi / d)
    cases = OrderedDict()
    for si, s in enumerate(series):
        ax_key = s["axis"]
        rho_cases = OrderedDict([("scan low", mids[ax_key]["low"]), ("scan midpoint", mids[ax_key]["mid"]),
                                 ("scan high", mids[ax_key]["high"])])
        if ax_key == "PD":
            rho_cases["K065 floor"] = k065
        for lab, rho in rho_cases.items():
            # EN/PT: the metric uses the solver at exactly this rho, not the nearest grid column
            ax = axes[ax_key]
            dp = ax["protein"][0][1]
            peaks = []
            for T in T_GRID:
                cas = m36.Cascade(ax["d"], rho, dp, s["decay_per_hour"], lags=())
                dose, _ = cas.required_dose(s["target_knockdown"], T, m36.required_excess(s["target_knockdown"], rho))
                peaks.append(None if dose is None else 1.0 + cas.cycle(dose, T)["e_peak"])
            entry = OrderedDict(rho=rho, peaks=peaks)
            for lv in THRESHOLDS:
                h, flag = interval_for_peak(T_GRID, peaks, lv)
                entry[f"interval_within_{lv:g}x_h"] = h
                entry[f"interval_within_{lv:g}x_flag"] = flag
            cases[f"{ax_key}|{s['species']}|{s['target_knockdown']:.4f}|{lab}"] = entry
    print("  longest interval that keeps the required peak within a threshold (scan midpoint):")
    for k, e in cases.items():
        if k.endswith("scan midpoint"):
            print(f"    {k:46s} " + "  ".join(f"<= {lv:g}x: {e[f'interval_within_{lv:g}x_h']:6.1f} h ({e[f'interval_within_{lv:g}x_flag']})" for lv in THRESHOLDS))

    # ------------------------------------------------------------------ write
    fields = ["axis", "species", "param_id", "half_life_h", "target", "target_knockdown", "rho", "interval_h",
              "reachable", "peak_pool_over_baseline", "knockdown_reached"]
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(fields)
        for si, s in enumerate(series):
            for ri, rho in enumerate(RHO_GRID):
                for ti, T in enumerate(T_GRID):
                    wr.writerow([s["axis"], s["species"], s["param_id"], f"{s['half_life_h']:.6g}", s["target"],
                                 f"{s['target_knockdown']:.10g}", f"{float(rho):.10g}", f"{T:g}",
                                 bool(reach[si, ri, ti]),
                                 "" if not reach[si, ri, ti] else f"{peak[si, ri, ti]:.10g}",
                                 f"{kd[si, ri, ti]:.10g}"])
    out = OrderedDict(
        note_en=("Peak of the repressor pool over its baseline pool that brings the cycle-averaged protein down by the "
                 "target, on a grid of repression strength and dosing interval; the same solver and constants as "
                 "scripts/36. A point is unreached when no dose up to the cap does."),
        note_pt=("Pico do pool do repressor sobre o pool basal que leva a proteina media no ciclo a queda-alvo, numa "
                 "grade de forca de repressao e intervalo de dosagem; o mesmo resolvedor e as mesmas constantes do "
                 "scripts/36. Um ponto fica sem alcance quando nenhuma dose ate o teto o faz."),
        dose_cap_in_baseline_pools=m36.DOSE_CAP,
        rho_grid=[float(x) for x in RHO_GRID], interval_grid_h=T_GRID, thresholds=list(THRESHOLDS),
        series=[dict(s) for s in series], verification=checks, interval_within_threshold=cases,
        unreached_points=int((~reach).sum()), points=int(reach.size))
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"EN/PT -> {OUT_CSV}  ({len(series) * len(RHO_GRID) * len(T_GRID)} rows)")
    print(f"EN/PT -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
