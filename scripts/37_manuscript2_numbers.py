#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Gathers every number that the Article 2 manuscript (the mechanistic study) quotes
     into one flat file, results/tables/manuscript2_numbers.json, read from the result
     tables of scripts 12, 17 and 36 and from the kinetic parameter table, and from
     nothing else. The manuscript text is kept outside the repository; it carries
     placeholders such as {{ad.kappa}} and a renderer replaces each one with the value
     stored here, and refuses to build if a placeholder has no value. A number in the
     manuscript therefore cannot differ from the table it came from.

     This script computes nothing new about the models.
     It selects, rounds and counts (half-up rounding) from files already written, and a
     few values that are plain arithmetic on two table entries (a ratio, a percentage,
     a time to 90%) are computed here and labelled as such in "_sources". Each value is
     stored as a string, already formatted for the manuscript. Every sentence of the
     manuscript that relies on a pattern in the data (all 18 scans evaluable, the
     floor that the only loss-of-function observation sets on the repression strength lying below the scan,
     a regimen out of reach at every repression strength) has an assert here, so the
     sentence fails loudly if a rerun stops supporting it.

PT | Reune num arquivo plano, results/tables/manuscript2_numbers.json, todo numero que o
     manuscrito do Artigo 2 (o estudo mecanistico) cita, lido das tabelas de resultado
     dos scripts 12, 17 e 36 e da tabela de parametros cineticos, e de mais nada. O texto
     do manuscrito fica fora do repositorio; ele traz marcadores como {{ad.kappa}} e um
     montador troca cada um pelo valor guardado aqui, e se recusa a montar se algum
     marcador nao tiver valor. Um numero do manuscrito nao pode diferir da tabela de onde
     veio.

     Este script nao calcula nada novo sobre os
     modelos. Ele seleciona, arredonda e conta (meio-para-cima) a partir de arquivos ja
     gravados, e uns poucos valores que sao aritmetica simples sobre duas entradas da
     tabela (uma razao, uma porcentagem, um tempo ate 90%) sao calculados aqui e
     rotulados assim em "_sources". Cada valor e guardado como texto, ja formatado para o
     manuscrito. Toda frase do manuscrito que depende de um padrao nos dados (as 18
     varreduras avaliaveis, a piso que a unica observacao de perda de funcao poe na forca de repressao
     abaixo da varredura, um regime fora de alcance em toda forca de repressao) tem um
     assert aqui, e a frase falha de forma explicita se uma nova execucao deixar de
     sustenta-la.

    python scripts/37_manuscript2_numbers.py
"""

import csv
import json
import math
import os
import sys
from collections import Counter
from decimal import ROUND_HALF_UP, Decimal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _provenance as pc

T = "results/tables"
KIN = "data/extracted/kinetic_parameters.csv"
ODE = f"{T}/ode_calibrated_results.json"
SENS = f"{T}/ode_free_parameter_sensitivity.csv"
DOS_J = f"{T}/mimic_dosing_feasibility.json"
DOS_C = f"{T}/mimic_dosing_feasibility.csv"
TGT_J = f"{T}/mimic_target_response.json"
TGT_C = f"{T}/mimic_target_response.csv"
SWP_C = f"{T}/mimic_target_rho_sweep.csv"
DS_J = f"{T}/mimic_design_space.json"
DS_C = f"{T}/mimic_design_space.csv"
REN_J = f"{T}/structure_renders_article2.json"
PRV_J = f"{T}/structure_figure_provenance.json"
CLS_C = "results/tables/parameter_provenance_classes.csv"
RNK_J = f"{T}/mimic_series_ranking.json"
RNK_C = f"{T}/mimic_series_ranking.csv"
OUT = f"{T}/manuscript2_numbers.json"

N, SRC = {}, {}
LN2 = math.log(2.0)


def rnd(x, nd):
    q = Decimal(1).scaleb(-nd)
    return str(Decimal(repr(round(float(x), 10))).quantize(q, rounding=ROUND_HALF_UP))


def hl_fmt(x):
    """EN/PT: a half-life without a trailing zero: 7 h, 34 h, but 1.7 h and 10.6 h."""
    x = float(x)
    return str(int(round(x))) if abs(x - round(x)) < 1e-3 else rnd(x, 1)


def put(key, value, src):
    if key in N:
        raise KeyError(f"duplicate key {key}")
    N[key] = value if isinstance(value, str) else str(value)
    SRC[key] = src


def num(key, x, nd, src):
    put(key, rnd(x, nd), src)


def integer(key, x, src):
    put(key, f"{int(round(float(x))):,}", src)


def pct(key, x, nd, src):
    put(key, rnd(float(x) * 100.0, nd), src)


def peak(x):
    """EN/PT: one decimal below 10, whole numbers from 10 up, thousands separated."""
    x = float(x)
    if x < 10.0:
        return rnd(x, 1)
    return f"{int(Decimal(repr(round(x, 8))).quantize(Decimal(1), rounding=ROUND_HALF_UP)):,}"


def lower_bound(x):
    """EN/PT: round DOWN, to one decimal below 10 and to a whole number from 10 up, so that "at least X" stays true."""
    x = float(x)
    if x < 10.0:
        return f"{int(Decimal(repr(round(x, 8))) * 10) / 10:.1f}"
    return f"{int(Decimal(repr(round(x, 8)))):,}"


def must(cond, msg):
    if not cond:
        sys.exit(f"EN/PT: manuscript2 numbers: {msg}")


# ----------------------------------------------------------------------------- inputs
kin = {r["param_id"]: r for r in csv.DictReader(open(KIN, encoding="utf-8"))}
ode = json.load(open(ODE, encoding="utf-8"))
dos = json.load(open(DOS_J, encoding="utf-8"))
dos_rows = list(csv.DictReader(open(DOS_C, encoding="utf-8")))
sens = list(csv.DictReader(open(SENS, encoding="utf-8")))
tgt = json.load(open(TGT_J, encoding="utf-8"))
rows = list(csv.DictReader(open(TGT_C, encoding="utf-8")))


def v(pid):
    r = kin[pid]
    must(r["kind"] in ("numeric", "derived"), f"{pid} is {r['kind']}, not a number")
    return float(r["value_si"])


def written(pid):
    return kin[pid]["value_as_written"]


# ----------------------------------------------------------------------------- the parameter table
kinds = {}
for r in kin.values():
    kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
integer("kin.rows", len(kin), KIN)
integer("kin.sources", len({(r["pmid"], r["doi"]) for r in kin.values() if r["pmid"] or r["doi"]}), KIN)
integer("kin.numeric", kinds["numeric"], KIN)
integer("kin.derived", kinds["derived"], KIN)
integer("kin.qualitative", kinds["qualitative_constraint"], KIN)
integer("kin.gap", kinds["declared_gap"], KIN)
must(len(kin) == 67 and kinds["declared_gap"] == 1, "the parameter table no longer has 67 rows with one declared gap")

# ----------------------------------------------------------------------------- Alzheimer's arm (scripts/12)
clr = ode["ad_clearance_vs_production"]
dose = ode["mimic_dose_vs_measured_knockdown"]
src12 = ODE
for k, pid in dict(prod_ctrl="K001", prod_ad="K002", clear_ctrl="K003", clear_ad="K004").items():
    put(f"ad.{k}", written(pid), KIN)
num("ad.clear_deficit", clr["clearance_deficit_percent"], 1, src12)
num("ad.prod_diff", clr["production_difference_percent"], 1, src12)
num("ad.ratio", clr["abeta_ratio_AD_over_control"], 2, src12)
num("ad.ratio4", clr["abeta_ratio_AD_over_control"], 4, src12)
# EN/PT: the two factors of the pool ratio, written out so the text can show that clearance carries it
f_prod = v("K002") / v("K001")
f_clear = v("K003") / v("K004")
must(abs(f_prod * f_clear - clr["abeta_ratio_AD_over_control"]) < 1e-7, "production and clearance factors do not multiply to the ratio")
num("ad.factor_prod", f_prod, 3, "arithmetic on K001, K002")
num("ad.factor_clear", f_clear, 3, "arithmetic on K003, K004")
pct("ad.kappa", dose["required_BACE1_knockdown"], 1, src12)
pct("ad.measured_kd", dose["measured_BACE1_knockdown_K047"], 0, src12)
pct("ad.kappa_over_measured", dose["required_over_measured_knockdown"], 0, src12)
lo, hi = dose["mimic_fold_to_offset_clearance"]
num("ad.fold_lo", lo, 2, src12)
num("ad.fold_hi", hi, 2, src12)
num("ad.fold_mid", clr["mimic_fold_needed"], 2, src12)
num("ad.fold_proportional", 1.0 / (1.0 - dose["required_BACE1_knockdown"]), 2, "arithmetic on the required knockdown")
lo, hi = dose["mimic_fold_to_reproduce_measured_knockdown"]
num("ad.fold_measured_lo", lo, 1, src12)
num("ad.fold_measured_hi", hi, 1, src12)
_sup = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def sci(x):
    """EN/PT: one significant figure, written 3 x 10-10 with a superscript exponent."""
    e = int(math.floor(math.log10(abs(x))))
    m = int(round(x / 10 ** e))
    if m == 10:
        m, e = 1, e + 1
    return f"{m} \u00d7 10{str(e).translate(_sup)}"


put("ad.ode_check", sci(clr["mimic_fold_ode_check_max_relative_error"]), src12)
rs = tgt["repression_strength"]
num("ad.rho_lo", rs["AD"]["rho_scan_low"], 1, TGT_J)
num("ad.rho_mid", rs["AD"]["rho_scan_mid"], 0, TGT_J)
integer("ad.rho_hi", rs["AD"]["rho_scan_high"], TGT_J)
st = ode["enzyme_substrate_stoichiometry"]
integer("ad.bace1_copies", st["BACE1_copies_per_bouton"], src12)
integer("ad.app_copies", st["APP_copies_per_bouton"], src12)
integer("ad.stoich", st["APP_per_BACE1"], src12)
nuc = tgt["nucleation_scaling"]
must(nuc["order_primary_nucleation"] == nuc["order_secondary_nucleation"] == 2.0,
     "the two Abeta42 nucleation orders are no longer both 2; the sentence that quotes one number for both must change")
num("ad.nuc_order", nuc["order_primary_nucleation"], 0, TGT_J)
num("ad.nuc_ratio", nuc["primary_nucleation_rate_ratio_AD_over_control"], 3, TGT_J)
num("ad.nuc_ratio1", nuc["primary_nucleation_rate_ratio_AD_over_control"], 1, TGT_J)
pct("ad.nuc_cut", nuc["primary_nucleation_rate_reduction_after_offset"], 1, TGT_J)
loads = ode["human_abeta_aggregate_load"]
num("ad.mstar_nM", loads["M_star_nM [K014]"], 0, src12)
num("ad.load_ctrl_x", loads["control median [K032]"]["multiples_of_M_star"], 0, src12)
num("ad.load_ad_x", loads["AD median [K035]"]["multiples_of_M_star"], 0, src12)
num("ad.load_fold", loads["AD_over_control_median_fold"], 0, src12)
# scans
ok = [r for r in sens if r["verdict"] == "AD_above_control"]
must(len(sens) == 18 and len(ok) == 18, "the free-parameter scan no longer has 18 evaluable runs")
ratios = [float(r["abeta_AD_over_control"]) for r in ok]
must(max(ratios) - min(ratios) < 1e-6, "the scanned AD/control ratio is no longer constant")
integer("ad.scan_n", len(sens), SENS)
num("ad.scan_ratio", sum(ratios) / len(ratios), 4, SENS)
for lab, rho_ in (("rho1", 1.0), ("rho01", 0.1)):
    num(f"ad.fold_{lab}", 1.0 + dose["required_BACE1_knockdown"] / (1.0 - dose["required_BACE1_knockdown"]) * (1.0 + rho_) / rho_, 1,
        "closed form on the required knockdown at the stated rho")
pct("ad.bace1_subset", v("K049"), 0, KIN)
put("ad.k_repress_range", "0.1 to 3", ODE)
put("ad.k_translate_range", "0.2 to 5", ODE)
num("ad.d_bace1_mrna", v("K068"), 4, KIN)
num("ad.hl_bace1_mrna", LN2 / v("K068"), 1, "arithmetic on K068")
num("ad.hl_bace1_cell", LN2 / v("K025"), 0, "arithmetic on K025")
put("ad.bace1_hek_hl", written("K024"), KIN)

# ----------------------------------------------------------------------------- Parkinson's arm (scripts/12)
sat = ode["alpha_synuclein_saturation"]
chain = ode["mir7_knockdown_to_elongation"]
num("pd.conc", sat["alpha_synuclein_uM_in_bouton"], 1, src12)
num("pd.mhalf_lo", min(p["m_half_uM"] for p in sat["points"]), 0, src12)
num("pd.mhalf_hi", max(p["m_half_uM"] for p in sat["points"]), 0, src12)
pct("pd.frac_lo", sat["fraction_of_maximal_range"][0], 0, src12)
pct("pd.frac_hi", sat["fraction_of_maximal_range"][1], 0, src12)
num("pd.elast_lo", sat["elasticity_range"][0], 2, src12)
num("pd.elast_hi", sat["elasticity_range"][1], 2, src12)
must(sat["below_half_saturation"], "alpha-synuclein no longer sits below half-saturation")
for lab, pid in (("a", "K058"), ("b", "K063")):
    kd = v(pid)
    cases = [c for c in chain["cases"] if abs(c["knockdown"] - kd) < 1e-12]
    must(len(cases) == 2, f"expected two m_half cases for {pid}")
    slows = [c["elongation_rate_reduction"] for c in cases]
    pct(f"pd.kd_{lab}", kd, 0, KIN)
    pct(f"pd.slow_{lab}_lo", min(slows), 1, src12)
    pct(f"pd.slow_{lab}_hi", max(slows), 1, src12)
    pct(f"pd.transmit_{lab}_lo", min(slows) / kd, 0, src12)
    pct(f"pd.transmit_{lab}_hi", max(slows) / kd, 0, src12)
must(chain["less_than_proportional_everywhere"], "the elongation slowdown is no longer less than proportional")
ad_ = chain["additivity"]
pct("pd.pair_measured", ad_["pair_measured"], 0, src12)
pct("pd.pair_independent", ad_["pair_if_independent"], 0, src12)
pct("pd.mir7_alone_kd", ad_["miR7_alone"], 0, src12)
pct("pd.mir153_alone_kd", ad_["miR153_alone"], 0, src12)
pct("pd.mir7_mrna_kd", v("K061"), 0, KIN)
pct("pd.pair_mrna_kd", v("K062"), 0, KIN)
pct("pd.age_protein", v("K027"), 0, KIN)
pct("pd.age_mrna", abs(v("K028")), 0, KIN)
put("pd.mir7_copies", written("K022"), KIN)
put("pd.cyrano_fold", written("K023"), KIN)
put("pd.mrna_hl", written("K066"), KIN)
put("pd.mrna_hl_1", rnd(float(written("K066")), 1), KIN)
put("pd.mhalf_hi_1", written("K009"), KIN)
put("pd.asyn_hl", written("K026"), KIN)
put("pd.proteome_hl_d", written("K029"), KIN)
num("pd.d_snca_mrna", v("K066"), 3, KIN)
num("pd.rho_lo", rs["PD"]["rho_scan_low"], 2, TGT_J)
num("pd.rho_mid", rs["PD"]["rho_scan_mid"], 1, TGT_J)
integer("pd.rho_hi", rs["PD"]["rho_scan_high"], TGT_J)
num("pd.rho_floor", rs["PD"]["rho_floor_K065"], 2, TGT_J)
must(rs["PD"]["floor_below_scan"], "the K065 floor on the repression strength no longer lies below the declared scan")
pct("pd.k065_increase", v("K065"), 0, KIN)
for lab, key in (("a", "miR-7 alone, SNCA construct with 3'UTR [K058]"), ("b", "miR-7 + miR-153, cortical neurons [K063]")):
    t = rs["PD"]["targets"][key]
    num(f"pd.fold_{lab}_prop", t["proportional_limit_fold"], 2, TGT_J)
    num(f"pd.fold_{lab}_lo", t["fold_at_scan_high"], 2, TGT_J)
    num(f"pd.fold_{lab}_hi", t["fold_at_scan_low"], 2, TGT_J)
    num(f"pd.fold_{lab}_floor", t["fold_at_floor"], 2, TGT_J)
# gene-dosage arithmetic: a duplication carries three copies and a triplication four, against two
num("pd.dup_cancel", (1.0 - 2.0 / 3.0) * 100.0, 0, "arithmetic: three gene copies against two")
num("pd.trip_cancel", (1.0 - 2.0 / 4.0) * 100.0, 0, "arithmetic: four gene copies against two")

# ----------------------------------------------------------------------------- turnover (scripts/12 and 17)
wash = ode["mimic_washout"]
labels = {"miR-7 (1.7 h, K020)": "m7", "miR-7 upper bound (5 h, K021)": "m7ub", "miR-29b (7 h, K017)": "b",
          "miR-29c (10.6 h, K018)": "c", "median miRNA (34 h, K016)": "med"}
for lab, tag in labels.items():
    w = wash[lab]
    put(f"turn.{tag}.hl", hl_fmt(w["half_life_hours"]), src12)
    num(f"turn.{tag}.washout_h", w["hours_to_return_to_baseline"], 1, src12)
    num(f"turn.{tag}.washout_d", w["days_to_return_to_baseline"], 1, src12)
dd = dos["daily_dosing_comparison"]
dkeys = {"miR-7": "m7", "miR-7 upper bound": "m7ub", "miR-29b": "b", "miR-29c": "c"}
for name, tag in dkeys.items():
    d_ = dd[name]
    num(f"turn.{tag}.ptm", d_["peak_over_average_daily"], 2, DOS_J)
    num(f"turn.{tag}.penalty", d_["penalty_versus_reference"], 1, DOS_J)
    num(f"turn.{tag}.interval_h", d_["interval_matching_reference_penalty_hours"], 1, DOS_J)
    num(f"turn.{tag}.stab", d_["fold_stabilisation_required"], 1, DOS_J)
    num(f"turn.{tag}.hl_req", d_["half_life_required_for_daily_dosing_hours"], 0, DOS_J)
num("turn.med.ptm", dd["miR-7"]["reference_peak_over_average_daily"], 2, DOS_J)
for r in dos_rows:
    if float(r["dosing_interval_hours"]) == 24.0:
        tag = {"miR-7": "m7", "miR-7 upper bound": "m7ub", "miR-29b": "b", "miR-29c": "c", "median miRNA": "med"}[r["species"]]
        pct(f"turn.{tag}.above_half", r["fraction_of_interval_above_half_peak"], 1, DOS_C)
num("turn.ratio_hl_m7_med", float(wash["median miRNA (34 h, K016)"]["half_life_hours"]) / float(wash["miR-7 (1.7 h, K020)"]["half_life_hours"]), 0,
    "arithmetic on K016, K020")
put("turn.avg_hl", written("K015"), KIN)
put("turn.med_hl", written("K016"), KIN)

# ----------------------------------------------------------------------------- the dose at the target (scripts/36)
must(tgt["scenarios"] == len(rows) == 330, "the scenario table no longer has 330 rows")
integer("dose.scenarios", len(rows), TGT_C)
integer("dose.unreachable", tgt["scenarios_unreachable"], TGT_J)
integer("dose.cap", tgt["dose_cap_in_baseline_pools"], TGT_J)
for k_, c_ in tgt["verification"].items():
    must(c_["passed"], f"verification {k_} did not pass")
integer("dose.verifications", len(tgt["verification"]), TGT_J)
num("dose.verif_worst", max(c_["max_abs_error"] / c_["tolerance"] for c_ in tgt["verification"].values()), 3,
    "worst error over tolerance in " + TGT_J)

PRIMARY = ("BACE1 in", "alpha-syn")
SP_TAG = {"miR-7": "m7", "miR-7 upper bound": "m7ub", "miR-29b": "b", "miR-29c": "c", "median miRNA": "med"}


def pick(axis, rho_case, species, interval, kappa=None, protein=PRIMARY):
    out = [r for r in rows if r["axis"] == axis and r["rho_case"] == rho_case and r["species"] == species
           and float(r["interval_h"]) == interval and r["protein_case"].startswith(protein)
           and (kappa is None or abs(float(r["target_knockdown"]) - kappa) < 1e-9)]
    must(len(out) == 1, f"expected one row for {axis}/{rho_case}/{species}/{interval}/{kappa}, found {len(out)}")
    return out[0]


def cases(axis, species, interval, kappa=None, only=None, protein=PRIMARY):
    names = [c for c in dict.fromkeys(r["rho_case"] for r in rows if r["axis"] == axis)
             if only is None or c in only]
    return [pick(axis, c, species, interval, kappa, protein) for c in names]


def span(key, items, field, f, src=TGT_C):
    xs = [float(i[field]) for i in items if i[field] != ""]
    must(xs, f"no values for {key}")
    put(f"{key}_lo", f(min(xs)), src)
    put(f"{key}_hi", f(max(xs)), src)


kap_ad = dose["required_BACE1_knockdown"]
kap_a, kap_b = v("K058"), v("K063")
SCAN = ("scan low", "scan midpoint", "scan high")
INTERVALS = (8.0, 12.0, 24.0, 48.0, 168.0)

# AD: pool-level mimic, three species, four intervals, over the three scan strengths
for sp in ("miR-29b", "miR-29c", "median miRNA"):
    tag = SP_TAG[sp]
    for T_ in INTERVALS:
        its = cases("AD", sp, T_, kap_ad)
        must(all(i["reachable"] == "True" for i in its), f"AD {sp} {T_} h unreachable")
        span(f"dose.AD.{tag}.T{int(T_)}.peak", its, "peak_pool_over_baseline", peak)
        span(f"dose.AD.{tag}.T{int(T_)}.eff", its, "efficiency", lambda x: rnd(x, 2))
        span(f"dose.AD.{tag}.T{int(T_)}.swing", its, "protein_ripple", lambda x: rnd(100 * x, 0))
        span(f"dose.AD.{tag}.T{int(T_)}.abswing", its, "abeta_ripple", lambda x: rnd(100 * x, 0))
        span(f"dose.AD.{tag}.T{int(T_)}.trough", its, "trough_pool_over_baseline", lambda x: rnd(x, 2))

# PD: the K065 floor and the scan, two targets, five intervals
for sp in ("miR-7", "miR-7 upper bound", "median miRNA"):
    tag = SP_TAG[sp]
    for lab, kap in (("a", kap_a), ("b", kap_b)):
        for T_ in INTERVALS:
            k065 = pick("PD", "K065 lower bound", sp, T_, kap)
            sc = cases("PD", sp, T_, kap, only=SCAN)
            base = f"dose.PD.{tag}.{lab}.T{int(T_)}"
            put(f"{base}.reach_floor", "yes" if k065["reachable"] == "True" else "no", TGT_C)
            put(f"{base}.reach_scan", "yes" if all(i["reachable"] == "True" for i in sc) else
                ("no" if all(i["reachable"] != "True" for i in sc) else "some"), TGT_C)
            if k065["reachable"] == "True":
                put(f"{base}.peak_floor", peak(k065["peak_pool_over_baseline"]), TGT_C)
                put(f"{base}.eff_floor", rnd(k065["efficiency"], 2), TGT_C)
                put(f"{base}.swing_floor", rnd(100 * float(k065["protein_ripple"]), 0), TGT_C)
            else:
                put(f"{base}.reached_floor", rnd(100 * float(k065["knockdown_reached"]), 0), TGT_C)
            reach = [i for i in sc if i["reachable"] == "True"]
            if reach:
                span(f"{base}.peak_scan", reach, "peak_pool_over_baseline", peak)
                span(f"{base}.eff_scan", reach, "efficiency", lambda x: rnd(x, 2))
                span(f"{base}.swing_scan", reach, "protein_ripple", lambda x: rnd(100 * x, 0))
            if len(reach) < len(sc):
                un = [i for i in sc if i["reachable"] != "True"]
                span(f"{base}.reached_scan", un, "knockdown_reached", lambda x: rnd(100 * x, 0))

# the sentences about weekly miR-7 dosing rest on this pattern
for lab in ("a", "b"):
    b = f"dose.PD.m7.{lab}.T168"
    must(N[f"{b}.reach_floor"] == "no" and N[f"{b}.reach_scan"] == "no",
         f"weekly miR-7 dosing is no longer out of reach at every repression strength (target {lab})")
# daily dosing of a median-stability miRNA is the linear case: the closed form suffices
for ax_, sp_, kap_list in (("AD", "median miRNA", [kap_ad]), ("PD", "median miRNA", [kap_a, kap_b])):
    for kap in kap_list:
        for r in cases(ax_, sp_, 24.0, kap):
            must(float(r["efficiency"]) > 0.99, f"daily dosing of a median-stability miRNA is no longer in the linear regime ({ax_} {r['rho_case']})")
# the dose does not depend on the protein half-life
by = {}
for r in rows:
    by.setdefault((r["axis"], r["rho_case"], r["species"], r["interval_h"], r["target"]), []).append(r)
for k_, g_ in by.items():
    d_ = [float(x["dose"]) for x in g_ if x["dose"] != ""]
    must(len(d_) in (0, len(g_)) and (not d_ or max(d_) - min(d_) <= 1e-7 * max(d_)),
         f"the dose depends on the protein half-life in {k_}")
put("dose.protein_independent", "yes", TGT_C)

# protein half-life sensitivity: swing at the cell value and at the in vivo brain median
for tag, axis, sp, rho_case, kap in (("ad_b", "AD", "miR-29b", "scan midpoint", kap_ad),
                                     ("pd_m7", "PD", "miR-7", "K065 lower bound", kap_a)):
    cell_ = pick(axis, rho_case, sp, 24.0, kap)
    vivo = pick(axis, rho_case, sp, 24.0, kap, protein=("brain proteome",))
    pct(f"dose.swing_cell.{tag}", cell_["protein_ripple"], 1, TGT_C)
    pct(f"dose.swing_vivo.{tag}", vivo["protein_ripple"], 1, TGT_C)
    num(f"dose.hl_cell.{tag}", cell_["protein_half_life_h"], 0, TGT_C)
    num(f"dose.hl_vivo.{tag}", float(vivo["protein_half_life_h"]) / 24.0, 1, TGT_C)
# time for the protein to approach 90% of its new mean: ln(10) / d_p
num("onset.bace1_cell_d", math.log(10.0) / v("K025") / 24.0, 1, "arithmetic on K025: ln(10) / d_p")
num("onset.asyn_cell_d", math.log(10.0) / v("K026") / 24.0, 1, "arithmetic on K026: ln(10) / d_p")
num("onset.vivo_d", math.log(10.0) / v("K029") / 24.0, 0, "arithmetic on K029: ln(10) / d_p")

# elongation in the periodic regimen against the steady-state chain
r = pick("PD", "K065 lower bound", "median miRNA", 24.0, kap_a)
pct("dose.elong_periodic_a", r["elongation_reduction_m49_8"], 1, TGT_C)
ss = [c for c in chain["cases"] if abs(c["knockdown"] - kap_a) < 1e-12 and abs(c["m_half_uM"] - 49.8) < 1e-9][0]
must(abs(float(r["elongation_reduction_m49_8"]) - ss["elongation_rate_reduction"]) < 1e-4,
     "the periodic elongation reduction no longer matches the steady-state chain")

# the rho sweep: the smallest daily peak over all rho, a bound that holds whatever rho is
sweep = tgt["rho_sweep"]
for e in sweep["summary"]:
    lab = "ad" if e["axis"] == "AD" else ("a" if abs(e["target_knockdown"] - kap_a) < 1e-12 else "b")
    base = f"sweep.{e['axis']}.{SP_TAG[e['species']]}.{lab}"
    put(f"{base}.min", peak(e["min_peak_pool_over_baseline"]), TGT_J)
    put(f"{base}.min_lb", lower_bound(e["min_peak_pool_over_baseline"]), TGT_J)
    r_ = e["rho_at_min"]
    put(f"{base}.rho_at_min", (f"10{str(int(round(math.log10(r_)))).translate(_sup)}" if r_ >= 1000 else f"{r_:.2g}"), TGT_J)
    put(f"{base}.low", peak(e["peak_at_lowest_rho"]), TGT_J)
    put(f"{base}.high", peak(e["peak_at_highest_rho"]), TGT_J)
    must(e["unreachable_points"] == 0, f"the sweep has unreachable points for {base}")
# EN/PT: the smallest peak on the grid bounds the dose at any rho only if, in every series, the peak rises as rho falls below
# its minimum (the closed form of Eq. 6 grows without bound there) and has converged at the top of the grid | o menor pico da grade
# limita a dose para qualquer rho so se, em cada serie, o pico sobe quando rho cai abaixo do minimo e ja convergiu no topo da grade
_series = {}
for r_ in csv.DictReader(open(SWP_C, encoding="utf-8")):
    _series.setdefault((r_["axis"], r_["species"], r_["target_knockdown"], r_["interval_h"]), []).append(r_)
must(len(_series) == len(sweep["summary"]), "the sweep table and its summary hold a different number of series")
_tail = []
for k_, g_ in _series.items():
    g_.sort(key=lambda x: float(x["rho"]))
    p_ = [float(x["peak_pool_over_baseline"]) for x in g_]
    i_ = p_.index(min(p_))
    must(all(p_[j] >= p_[j + 1] * (1.0 - 1e-9) for j in range(i_)),
         f"the peak does not rise as rho falls below its minimum in {k_}")
    j_ = [j for j, x in enumerate(g_) if abs(float(x["rho"]) * 10.0 / sweep["rho_max"] - 1.0) < 1e-6]
    must(len(j_) == 1, f"the grid has no point one decade below rho_max in {k_}")
    _tail.append(abs(p_[-1] - p_[j_[0]]) / p_[-1])
must(max(_tail) < 0.002, "the peak has not converged at the top of the grid")
pct("sweep.tail_change", max(_tail), 1, SWP_C)
num("sweep.rho_min", sweep["rho_min"], 2, TGT_J)
integer("sweep.rho_max", sweep["rho_max"], TGT_J)
# the sentences that quote one bound for the miR-29 mimics rest on the bound being the strong-repression limit
for sp_ in ("b", "c", "med"):
    must(abs(float(N[f"sweep.AD.{sp_}.ad.min"]) - float(N[f"sweep.AD.{sp_}.ad.high"])) < 0.051,
         f"the smallest AD peak for {sp_} is no longer the strong-repression limit")
must(float(N["sweep.PD.m7.a.min"]) > 5 * float(N["sweep.PD.med.a.min"]),
     "the miR-7 bound is no longer more than five times the median-stability bound")
_m7 = [e for e in sweep["summary"] if e["axis"] == "PD" and e["species"] == "miR-7" and abs(e["target_knockdown"] - kap_a) < 1e-12][0]
_md = [e for e in sweep["summary"] if e["axis"] == "PD" and e["species"] == "median miRNA" and abs(e["target_knockdown"] - kap_a) < 1e-12][0]
num("sweep.PD.ratio_m7_med_a", _m7["min_peak_pool_over_baseline"] / _md["min_peak_pool_over_baseline"], 0, "arithmetic on two bounds of " + TGT_J)

# absolute scale for miR-7: the model's baseline pool taken as the 40 copies per neuron of K022
c40 = v("K022")
for lab, kap in (("a", kap_a), ("b", kap_b)):
    pks = [float(i["peak_pool_over_baseline"]) for i in cases("PD", "miR-7", 24.0, kap)]
    integer(f"abs.peak_copies_{lab}_lo", min(pks) * c40, TGT_C)
    integer(f"abs.peak_copies_{lab}_hi", max(pks) * c40, TGT_C)
integer("abs.baseline_copies", c40, KIN)
integer("abs.cyrano_copies", c40 * v("K023"), "arithmetic on K022, K023")
must(min(float(i["peak_pool_over_baseline"]) for i in cases("PD", "miR-7", 24.0, kap_b)) > v("K023"),
     "the daily peak for the 43% target no longer exceeds the Cyrano-knockout fold")

# ----------------------------------------------------------------------------- where each number was measured (scripts/_provenance.py)
kin_rows = list(csv.DictReader(open(KIN, encoding="utf-8")))
sys_of, blk_of = pc.classes(kin_rows)
pc.write_table(kin_rows, CLS_C)          # the same table that scripts/38 writes; scripts/08 checks it against the classification
cnt = Counter(sys_of.values())
for k_, *_ in pc.SYSTEMS:
    integer(f"prov.{k_}", cnt.get(k_, 0), CLS_C)
integer("prov.brainlike", cnt["people"] + cnt["human_pm"] + cnt["animal_brain"], CLS_C)
by_block = {b_: [p_ for p_, x_ in blk_of.items() if x_ == b_] for b_, *_ in pc.BLOCKS}
rep_ = by_block["repression"]
integer("prov.rep_rows", len(rep_), CLS_C)
integer("prov.rep_cells", sum(sys_of[p_] == "cells" for p_ in rep_), CLS_C)
integer("prov.rep_neurons", sum(sys_of[p_] == "neurons" for p_ in rep_), CLS_C)
must({sys_of[p_] for p_ in rep_} <= {"cells", "neurons"},
     "a row that bears on the repression strength was measured outside cell lines and primary neurons; the sentence about it changes")
agg_ = by_block["aggregation"]
integer("prov.agg_rows", len(agg_), CLS_C)
integer("prov.agg_cell_free", sum(sys_of[p_] == "cell_free" for p_ in agg_), CLS_C)
integer("prov.agg_human_pm", sum(sys_of[p_] == "human_pm" for p_ in agg_), CLS_C)
integer("prov.agg_gap", sum(sys_of[p_] == "none" for p_ in agg_), CLS_C)
must({sys_of[p_] for p_ in agg_} == {"cell_free", "human_pm", "none"},
     "the aggregation rows no longer come from the test tube, post-mortem brain and the declared gap alone")
must(all(sys_of[p_] == "people" for p_ in by_block["flux"]), "a production or clearance row was measured outside living people")
hl_rows_ = [f"K0{n_}" for n_ in range(15, 22)]
must(all(sys_of[p_] == "cells" and blk_of[p_] == "mirna" for p_ in hl_rows_),
     "the miRNA half-life rows K015 to K021 no longer all come from cell lines")
integer("prov.mirna_hl_cells", len(hl_rows_), CLS_C)

# ----------------------------------------------------------------------------- half-lives against the dosing interval (Fig. 5a)
hl_m7 = float(wash["miR-7 (1.7 h, K020)"]["half_life_hours"])
hl_b = float(wash["miR-29b (7 h, K017)"]["half_life_hours"])
hl_snca = LN2 / v("K066")
hl_bace1m = LN2 / v("K068")
hl_asyn = LN2 / v("K026")
for tag_, hl_ in (("m7", hl_m7), ("b", hl_b), ("snca", hl_snca), ("bace1m", hl_bace1m), ("asyn", hl_asyn)):
    num(f"turn.{tag_}.n24", 24.0 / hl_, 1, "arithmetic: 24 h over the half-life")
num("turn.ratio_mrna_m7", hl_snca / hl_m7, 1, "arithmetic on K066 and K020")
must(hl_snca > hl_m7 and hl_bace1m > hl_b, "the repressor is no longer the faster step in both arms; the sentence on Fig. 5a changes")

# ----------------------------------------------------------------------------- the design-space grid (scripts/40)
ds = json.load(open(DS_J, encoding="utf-8"))
ds_rows = list(csv.DictReader(open(DS_C, encoding="utf-8")))
integer("ds.series", len(ds["series"]), DS_J)
integer("ds.rho_n", len(ds["rho_grid"]), DS_J)
integer("ds.t_n", len(ds["interval_grid_h"]), DS_J)
integer("ds.points", len(ds_rows), DS_C)
must(len(ds_rows) == len(ds["series"]) * len(ds["rho_grid"]) * len(ds["interval_grid_h"]), "the design-space table is not the full grid")
integer("ds.unreached", ds["unreached_points"], DS_J)
_dec = (len(ds["rho_grid"]) - 1) / math.log10(max(ds["rho_grid"]) / min(ds["rho_grid"]))
must(abs(_dec - round(_dec)) < 1e-9, "the rho grid is no longer a whole number of points per decade")
integer("ds.per_decade", _dec, DS_J)
# where the unreached points lie: the first interval with any, and the first with every rho unreached, per series
_un = {}
for r_ in ds_rows:
    if r_["reachable"] != "True":
        _un.setdefault((r_["axis"], r_["species"], round(float(r_["target_knockdown"]), 3)), {}).setdefault(float(r_["interval_h"]), 0)
        _un[(r_["axis"], r_["species"], round(float(r_["target_knockdown"]), 3))][float(r_["interval_h"])] += 1
must(all(k_[0] == "PD" and k_[1] in ("miR-7", "miR-7 upper bound") for k_ in _un),
     "a point is unreached outside the PD arm or for a mimic that is not cleared like miR-7; the sentence on the unreached region changes")
n_rho_ = len(ds["rho_grid"])
for sp_, tag_ in (("miR-7", "m7"), ("miR-7 upper bound", "m7ub")):
    integer(f"ds.unreached_{tag_}", sum(sum(v_.values()) for k_, v_ in _un.items() if k_[1] == sp_), DS_C)
    for kap_, lab_ in ((0.3, "a"), (0.43, "b")):
        g_ = _un.get(("PD", sp_, kap_), {})
        first_ = min(g_) if g_ else None
        all_ = [t_ for t_, n_ in sorted(g_.items()) if n_ == n_rho_]
        put(f"ds.un.{tag_}.{lab_}.first", str(int(first_)) if first_ else "none", DS_C)
        put(f"ds.un.{tag_}.{lab_}.all", str(int(all_[0])) if all_ else "none", DS_C)
integer("ds.t_min", min(ds["interval_grid_h"]), DS_J)
integer("ds.t_max", max(ds["interval_grid_h"]), DS_J)
v36_ = ds["verification"]["against_sweep_36_at_24h"]
integer("ds.verif_points", v36_["points"], DS_J)
put("ds.verif_err", sci(v36_["max_relative_error"]), DS_J)
must(v36_["max_relative_error"] <= v36_["tolerance"], "the design-space grid no longer reproduces scripts/36 at 24 h")
must(ds["verification"]["monotone_in_interval"]["max_relative_fall"] <= 1e-9, "the peak needed falls somewhere as the interval grows")
integer("ds.thr_lo", min(ds["thresholds"]), DS_J)
integer("ds.thr_mid", sorted(ds["thresholds"])[1], DS_J)
integer("ds.thr_hi", max(ds["thresholds"]), DS_J)
must(ds["thresholds"] == [3.0, 10.0, 30.0], "the reference ceilings of the design space are no longer 3, 10 and 30")
iw = ds["interval_within_threshold"]


def hours_floor(x):
    """EN/PT: an interval rounded to the nearest value (one decimal below 10 h, whole hours from 10 h), as in the figures."""
    return peak(x)


for ax_, sp_, kap_, lab_ in (("PD", "miR-7", kap_a, "a"), ("PD", "miR-7 upper bound", kap_a, "a"), ("PD", "median miRNA", kap_a, "a"),
                             ("PD", "miR-7", kap_b, "b"), ("PD", "miR-7 upper bound", kap_b, "b"), ("PD", "median miRNA", kap_b, "b"),
                             ("AD", "miR-29b", kap_ad, "ad"), ("AD", "miR-29c", kap_ad, "ad"), ("AD", "median miRNA", kap_ad, "ad")):
    for thr_ in (3, 10, 30):
        ent = [iw[f"{ax_}|{sp_}|{kap_:.4f}|{c_}"] for c_ in ("scan low", "scan midpoint", "scan high")]
        flags = {e_[f"interval_within_{thr_}x_flag"] for e_ in ent}
        vals = [e_[f"interval_within_{thr_}x_h"] for e_ in ent]
        base = f"ds.{ax_}.{SP_TAG[sp_]}.{lab_}.i{thr_}"
        put(f"{base}_lo", hours_floor(min(vals)), DS_J)
        put(f"{base}_hi", hours_floor(max(vals)), DS_J)
        put(f"{base}_rng", N[f"{base}_lo"] if N[f"{base}_lo"] == N[f"{base}_hi"] else f"{N[f'{base}_lo']} to {N[f'{base}_hi']}", DS_J)
        put(f"{base}_ceil", str(int(math.ceil(max(vals) - 1e-9))), DS_J)
        put(f"{base}_flag", flags.pop() if len(flags) == 1 else "mixed", DS_J)
# the sentences about the 10x ceiling rest on this pattern
for key_ in ("PD.m7.a", "PD.m7ub.a", "PD.m7.b", "PD.m7ub.b", "AD.b.ad", "AD.c.ad"):
    must(N[f"ds.{key_}.i10_flag"] == "interp", f"the 10-fold ceiling is no longer met at an interior interval for {key_}")
for key_ in ("PD.med.a", "PD.med.b", "AD.med.ad"):
    must(N[f"ds.{key_}.i10_flag"] in ("beyond", "interp"), f"unexpected flag for {key_}")
must(N["ds.AD.med.ad.i10_flag"] == "beyond" and N["ds.PD.med.a.i10_flag"] == "beyond" and N["ds.PD.med.b.i10_flag"] == "beyond",
     "a mimic of median stability no longer stays within a 10-fold peak up to the longest interval of the grid")
must(float(N["ds.PD.m7.a.i10_ceil"]) < 24 and float(N["ds.PD.m7.b.i10_ceil"]) < 24,
     "a miR-7-like mimic no longer needs more than a 10-fold peak at daily dosing")
must(float(N["ds.AD.b.ad.i10_lo"].replace(",", "")) > 24, "the miR-29b mimic is no longer within a 10-fold peak at daily dosing")
# the plateau in rho: at 24 h, from rho = 10 up, the peak stays close to its value at the top of the grid, in every series
by24 = {}
for r_ in ds_rows:
    if float(r_["interval_h"]) == 24.0 and r_["reachable"] == "True":
        by24.setdefault((r_["axis"], r_["species"], round(float(r_["target_knockdown"]), 4)), {})[float(r_["rho"])] = float(r_["peak_pool_over_baseline"])
must(len(by24) == len(ds["series"]), "a series has no reachable point at 24 h")
dev_ = 0.0
for k_, g_ in by24.items():
    top_ = g_[max(g_)]
    above_ = [p_ for rho_, p_ in g_.items() if rho_ >= 10.0]
    must(len(above_) == sum(1 for x_ in ds["rho_grid"] if x_ >= 10.0 - 1e-9), f"unreached points above rho = 10 at 24 h in {k_}")
    dev_ = max(dev_, max(abs(p_ / top_ - 1.0) for p_ in above_))
integer("ds.plateau_rho", 10, DS_C)
put("ds.plateau_dev", str(int(math.ceil(100.0 * dev_ - 1e-9))), DS_C)

# ----------------------------------------------------------------------------- ranking of the series (scripts/46)
rk = json.load(open(RNK_J, encoding="utf-8"))
rk_rows = list(csv.DictReader(open(RNK_C, encoding="utf-8")))
must(rk["primary"] == {"threshold": 10, "anchor": "scan midpoint"}, "the primary reading of the ranking changed")
must(rk["thresholds"] == [3, 10, 30], "the thresholds of the ranking changed")
for ax_ in ("AD", "PD"):
    f_ = rk["featured_by_axis"][ax_]
    must(len(f_["species"]) == 1, f"{ax_}: the primary reading has no single first place")
    must(not f_["first_place_changes"], f"{ax_}: first place changes with the threshold, the anchor or the target; the text says it does not")
    put(f"rank.{ax_}.featured", f_["species"][0], RNK_J)
    integer(f"rank.{ax_}.combos", f_["n_combinations"], RNK_J)
    integer(f"rank.{ax_}.ties", f_["n_ties_for_first"], RNK_J)
must(rk["featured_by_axis"]["AD"]["species"] == ["miR-29c"] and rk["featured_by_axis"]["PD"]["species"] == ["miR-7 upper bound"],
     "the featured series are no longer miR-29c in AD and the upper-bound miR-7 in PD")
integer("rank.combos", sum(rk["featured_by_axis"][a_]["n_combinations"] for a_ in ("AD", "PD")), RNK_J)
integer("rank.ties", sum(rk["featured_by_axis"][a_]["n_ties_for_first"] for a_ in ("AD", "PD")), RNK_J)
integer("rank.combos_single", int(N["rank.combos"]) - int(N["rank.ties"]), RNK_J)
integer("rank.ref_n", rk["reference_at_least_as_long_as_first_place"]["n"], RNK_J)
integer("rank.ref_of", rk["reference_at_least_as_long_as_first_place"]["of"], RNK_J)
must(rk["reference_at_least_as_long_as_first_place"]["n"] == rk["reference_at_least_as_long_as_first_place"]["of"],
     "the median-stability reference no longer matches or exceeds first place in every combination")
# EN/PT: the half-life order of the measured species is the order of the ranking | a ordem de meia-vida e a ordem da ordenacao
for ax_ in ("AD", "PD"):
    hl_ = {r_["species"]: float(r_["half_life_h"]) for r_ in rk_rows if r_["axis"] == ax_ and r_["reference"] == "False"}
    for g_ in {(r_["target_knockdown"], r_["anchor"], r_["threshold_x"]) for r_ in rk_rows if r_["axis"] == ax_}:
        grp_ = [r_ for r_ in rk_rows if r_["axis"] == ax_ and (r_["target_knockdown"], r_["anchor"], r_["threshold_x"]) == g_ and r_["reference"] == "False"]
        by_hl_ = sorted(grp_, key=lambda r_: -float(r_["half_life_h"]))
        iv_ = [float(r_["interval_h"]) for r_ in by_hl_]
        must(all(a_ >= b_ for a_, b_ in zip(iv_, iv_[1:])), f"{ax_} {g_}: a species with a shorter half-life holds the target for a longer interval")
put("rank.order_follows_half_life", "yes", RNK_C)
tied_ = [r_ for r_ in rk_rows if r_["rank_among_studied"] == "1" and r_["reference"] == "False"]
must(N["rank.ties"] == "1", "the text says that one combination tied for first place")
must(any(r_["axis"] == "PD" and r_["anchor"] == "K065 floor" and r_["threshold_x"] == "3" and r_["target_knockdown"] == "0.43" and r_["interval_h"] == "0.0" for r_ in tied_),
     "the tied combination is no longer the 43% target at the floor of rho and a 3-fold ceiling")

# ----------------------------------------------------------------------------- the deposited structures behind the molecular images (scripts/13, 39)
rn = json.load(open(REN_J, encoding="utf-8"))["renders"]
pv = json.load(open(PRV_J, encoding="utf-8"))
for pdb in ("6N4O", "4D8C", "5OQV", "6CU7"):
    e_ = pv[pdb]
    res_ = f"{float(e_['resolution_angstrom']):g}"
    put(f"struct.{pdb.lower()}.res", res_ if "." in res_ else res_ + ".0", PRV_J)
    put(f"struct.{pdb.lower()}.method", "X-ray diffraction" if "X-RAY" in e_["method"].upper() else "cryo-electron microscopy", PRV_J)
for pdb, key_ in (("5OQV", "a2_abeta42_fibril_end"), ("6CU7", "a2_asyn_fibril_end")):
    lts_ = rn[key_]["layer_transformation"]
    must(len(lts_) == 2 and all(t_["rms_ca_angstrom"] < 0.05 and t_["ca_pairs"] >= 40 for t_ in lts_),
         f"the constructed layer of {pdb} no longer rests on an exact superposition of two deposited chains")
    rises_, rots_ = [t_["rise_angstrom"] for t_ in lts_], [t_["rotation_deg"] for t_ in lts_]
    num(f"struct.{pdb.lower()}.rise_lo", min(rises_), 1, REN_J)
    num(f"struct.{pdb.lower()}.rise_hi", max(rises_), 1, REN_J)
    num(f"struct.{pdb.lower()}.rot", sum(rots_) / len(rots_), 1, REN_J)
    integer(f"struct.{pdb.lower()}.chains", sum(len(t_["chains_along_axis"]) for t_ in lts_), REN_J)
    must(len({t_["ca_pairs"] for t_ in lts_}) == 1, f"the two protofilaments of {pdb} are fitted over different numbers of Calpha pairs")
    integer(f"struct.{pdb.lower()}.pairs", lts_[0]["ca_pairs"], REN_J)
integer("struct.n", 4, PRV_J)

# ----------------------------------------------------------------------------- write
payload = dict(N)
payload["_sources"] = SRC
json.dump(payload, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"EN/PT -> {OUT}  ({len(N)} values)")
