#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Which of the nine mimic series holds a target down for the longest dosing interval, and
     whether that answer depends on the choices that were made.

     The criterion is fixed here, before any figure uses it, and it uses only quantities that
     scripts/40 already defines: for each series, each target knockdown and each anchor of the
     repression strength rho (the low, midpoint and high ends of the scan, and the floor set by
     the one observation that bears on rho, K065), the longest dosing interval at which the
     peak of the repressor pool needed to hold the cycle-averaged protein down by the target
     stays within 3, 10 or 30 times the baseline pool (`interval_within_Tx_h` of
     results/tables/mimic_design_space.json; 0 when no interval on the grid qualifies, and 168 h,
     the longest interval on the grid, when every interval does). The primary reading is 10 times
     at the midpoint of the scan. Series are ranked within an axis and a target, among the series
     that are measured species; the series labelled "median miRNA" is a reference of median
     stability, not a mimic that was studied, and is ranked separately.

     Nothing is selected by looking at an outcome: every combination of threshold and anchor is
     ranked, the table lists all of them, and the file records whether the first place changes
     with the threshold, the anchor or the target. The featured series of a figure is the first
     place of the primary reading; the same figure shows all the others.

PT | Qual das nove series de mimetico mantem um alvo baixo pelo maior intervalo de dose, e se essa
     resposta depende das escolhas feitas.

     O criterio e fixado aqui, antes de qualquer figura usa-lo, e usa so grandezas que o
     scripts/40 ja define: para cada serie, cada reducao-alvo e cada ancora da forca de repressao
     rho (as pontas baixa, media e alta da varredura, e o piso dado pela unica observacao que
     toca rho, K065), o maior intervalo de dose em que o pico do reservatorio do repressor,
     necessario para manter a proteina media do ciclo abaixo pelo alvo, fica ate 3, 10 ou 30
     vezes o reservatorio basal (`interval_within_Tx_h` de results/tables/mimic_design_space.json;
     0 quando nenhum intervalo da grade serve, e 168 h, o maior intervalo da grade, quando todos
     servem). A leitura primaria e 10 vezes no ponto medio da varredura. As series sao
     ordenadas dentro de um eixo e de um alvo, entre as que sao especies medidas; a serie
     "median miRNA" e uma referencia de estabilidade mediana, nao um mimetico estudado, e e
     ordenada a parte.

     Nada e escolhido olhando um resultado: toda combinacao de limiar e ancora e ordenada, a tabela
     lista todas, e o arquivo registra se o primeiro lugar muda com o limiar, a ancora ou o alvo. A
     serie em destaque de uma figura e o primeiro lugar da leitura primaria; a mesma figura mostra
     todas as outras.

    python scripts/46_mimic_series_ranking.py     # needs scripts/40
"""

import csv
import json
import os
import sys

SRC = "results/tables/mimic_design_space.json"
OUT_CSV = "results/tables/mimic_series_ranking.csv"
OUT_JSON = "results/tables/mimic_series_ranking.json"
REFERENCE = "median miRNA"
THRESHOLDS = (3, 10, 30)
PRIMARY = {"threshold": 10, "anchor": "scan midpoint"}


def main():
    if not os.path.exists(SRC):
        sys.exit(f"EN/PT: missing {SRC}; run scripts/40 first")
    d = json.load(open(SRC, encoding="utf-8"))
    meta = {(s["axis"], s["species"], round(s["target_knockdown"], 4)): s for s in d["series"]}
    rows = []
    for key, v in d["interval_within_threshold"].items():
        axis, species, kd, anchor = key.split("|")
        s = meta[(axis, species, round(float(kd), 4))]
        for t in THRESHOLDS:
            h, flag = v[f"interval_within_{t}x_h"], v[f"interval_within_{t}x_flag"]
            rows.append(dict(axis=axis, target=s["target"], target_knockdown=round(float(kd), 4), species=species,
                             param_id=s["param_id"], half_life_h=round(s["half_life_h"], 2), anchor=anchor,
                             rho=round(v["rho"], 4), threshold_x=t, interval_h=round(float(h or 0.0), 2), flag=flag,
                             reference=species == REFERENCE))
    # EN/PT: rank among the measured species of each axis, target, anchor and threshold | ordena entre as especies medidas
    groups = {}
    for r in rows:
        groups.setdefault((r["axis"], r["target_knockdown"], r["anchor"], r["threshold_x"]), []).append(r)
    for g in groups.values():
        studied = [r for r in g if not r["reference"]]
        for r in studied:   # EN/PT: competition ranking, equal intervals share a rank | empates dividem a posicao
            r["rank_among_studied"] = 1 + sum(1 for o in studied if o["interval_h"] > r["interval_h"])
        for r in g:
            r.setdefault("rank_among_studied", "")
    cols = ["axis", "target", "target_knockdown", "species", "param_id", "half_life_h", "anchor", "rho", "threshold_x",
            "interval_h", "flag", "reference", "rank_among_studied"]
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["axis"], r["target_knockdown"], r["threshold_x"], r["anchor"], r["rank_among_studied"] == "", r["rank_among_studied"] or 0)):
            w.writerow({c: r[c] for c in cols})
    # EN/PT: first place of every combination, and whether it ever changes | primeiro lugar de cada combinacao
    first = {}
    for k, g in groups.items():
        top = [r for r in g if r["rank_among_studied"] == 1]
        first[k] = top[0]["species"] if len(top) == 1 else None      # None: a tie for first place
    ties = [k for k, v in first.items() if v is None]
    featured = {}
    for axis in sorted({k[0] for k in first}):
        combos = {k: v for k, v in first.items() if k[0] == axis}
        names = sorted({v for v in combos.values() if v})
        prim = {k: v for k, v in combos.items() if k[2] == PRIMARY["anchor"] and k[3] == PRIMARY["threshold"]}
        n_tie = sum(1 for v in combos.values() if v is None)
        featured[axis] = dict(species=sorted({v for v in prim.values() if v}), first_place_in_every_combination=names,
                              n_combinations=len(combos), n_ties_for_first=n_tie,
                              first_place_changes=len(names) > 1)
    ref_beats = {}
    for k, g in groups.items():
        ref = [r for r in g if r["reference"]]
        top = [r for r in g if r["rank_among_studied"] == 1]
        if ref and top:
            ref_beats[k] = ref[0]["interval_h"] >= max(r["interval_h"] for r in top)
    out = dict(
        criterion_en="Longest dosing interval at which the repressor peak needed to hold the cycle-averaged protein down by the target stays within T times the baseline pool, for T = 3, 10 and 30 and at each anchor of the repression strength; primary reading T = 10 at the midpoint of the scan. Ranked among measured species; the series 'median miRNA' is a reference and is ranked apart.",
        criterion_pt="Maior intervalo de dose em que o pico do repressor necessario para manter a proteina media do ciclo abaixo pelo alvo fica ate T vezes o reservatorio basal, para T = 3, 10 e 30 e em cada ancora da forca de repressao; leitura primaria T = 10 no ponto medio da varredura. Ordenado entre especies medidas; a serie 'median miRNA' e uma referencia e e ordenada a parte.",
        primary=PRIMARY, thresholds=list(THRESHOLDS), reference_series=REFERENCE,
        featured_by_axis=featured, combinations_with_a_tie=[list(k) for k in ties],
        reference_at_least_as_long_as_first_place=dict(n=sum(ref_beats.values()), of=len(ref_beats)),
        source=SRC)
    with open(OUT_JSON, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
        fh.write("\n")
    print("=" * 78)
    print("EN | Ranking of the mimic series | PT | Ordenacao das series de mimetico")
    print("=" * 78)
    for axis, f in featured.items():
        print(f"{axis}: featured {f['species']} | first place in every combination: {f['first_place_in_every_combination']} "
              f"({f['n_combinations']} combinations, {f['n_ties_for_first']} tied) | changes: {f['first_place_changes']}")
    print(f"ties: {len(ties)} | reference series at least as long as first place: {sum(ref_beats.values())} of {len(ref_beats)}")
    print(f"EN/PT -> {OUT_CSV}, {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
