#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Where each row of the kinetic parameter table was measured, and which block of the models
     it feeds. The two classifications are declared here, one list of parameter ids per class,
     and checked against the table: a row may appear once, none may be missing, and the system
     has to agree with what the table says about the species, the population, the condition and
     the method. scripts/38 draws them (Fig. 2 of the mechanistic study), scripts/37 counts them
     for the manuscript and scripts/08 checks the table that this module writes.

PT | Onde cada linha da tabela de parametros cineticos foi medida e a que bloco dos modelos
     ela alimenta. As duas classificacoes sao declaradas aqui, uma lista de ids por classe, e
     conferidas com a tabela: uma linha aparece uma vez, nenhuma falta, e o sistema tem de
     concordar com o que a tabela diz de especie, populacao, condicao e metodo. O scripts/38
     as desenha (Fig. 2 do estudo mecanistico), o scripts/37 as conta para o manuscrito e o
     scripts/08 confere a tabela que este modulo grava.
"""

import csv
import os
import sys

OUT = "results/tables/parameter_provenance_classes.csv"
FIELDS = ["param_id", "kind", "axis", "source", "system", "block"]

SYSTEMS = [   # EN/PT: key, label EN, label PT; the system in which a number was measured, nearest to the living human brain first
    ("people", "People, in vivo", "Pessoas, in vivo"),
    ("human_pm", "Human post-mortem brain", "Cérebro humano post-mortem"),
    ("animal_brain", "Rodent brain", "Cérebro de roedor"),
    ("neurons", "Primary neurons", "Neurônios primários"),
    ("cells", "Cell lines and fibroblasts", "Linhagens e fibroblastos"),
    ("cell_free", "Cell-free, recombinant", "Sem células, recombinante"),
    ("none", "No measurement", "Sem medida"),
]
BLOCKS = [    # EN/PT: key, label EN, label PT; the block of the model that a number feeds
    ("mirna", "miRNA pool and turnover", "Pool e renovação do miRNA"),
    ("mrna", "mRNA turnover", "Renovação do mRNA"),
    ("repression", "Repression strength", "Força da repressão"),
    ("protein", "Protein turnover and abundance", "Renovação e abundância da proteína"),
    ("flux", "Aβ42 production and clearance", "Produção e depuração de Aβ42"),
    ("aggregation", "Aggregation kinetics", "Cinética de agregação"),
]
SYSTEM_OF = {
    "people": "K001 K002 K003 K004 K005 K006",
    "human_pm": "K027 K028 K031 K032 K033 K034 K035 K036 K045 K046 K049",
    "animal_brain": "K029 K050 K052 K053 K054 K055 K056 K057",
    "neurons": "K022 K023 K025 K043 K044 K063 K064 K065 K066 K067 K068 K069",
    "cells": "K015 K016 K017 K018 K019 K020 K021 K024 K026 K040 K041 K047 K048 K058 K059 K060 K061 K062",
    "cell_free": "K007 K008 K009 K010 K011 K012 K013 K014 K038 K039 K042",
    "none": "K037",
}
BLOCK_OF = {
    "mirna": "K015 K016 K017 K018 K019 K020 K021 K022 K023 K045 K046 K049",
    "mrna": "K040 K043 K044 K066 K067 K068 K069",
    "repression": "K047 K048 K058 K059 K060 K061 K062 K063 K064 K065",
    "protein": "K024 K025 K026 K027 K028 K029 K041 K050 K052 K053 K054 K055 K056 K057",
    "flux": "K001 K002 K003 K004 K005 K006",
    "aggregation": "K007 K008 K009 K010 K011 K012 K013 K014 K031 K032 K033 K034 K035 K036 K037 K038 K039 K042",
}
SURNAME = {"Schwanhausser": "Schwanhäusser"}
AXIS = {"miR-29/BACE1/Abeta": "AD", "miR-7/SNCA/alpha-synuclein": "PD", "both": "both"}


def _invert(spec):
    out = {}
    for key, ids in spec.items():
        for pid in ids.split():
            if pid in out:
                sys.exit(f"EN/PT: {pid} is listed twice in the provenance classes")
            out[pid] = key
    return out


def _system_from_text(r):
    """EN/PT: the system, read from the species, population, condition and method of the row | o sistema, lido da propria linha."""
    sp = r["species"].lower()
    txt = " ".join(r[k].lower() for k in ("species", "population", "condition", "method"))
    if sp == "not applicable":
        return "none"
    if "recombinant" in sp:
        return "cell_free"
    if "post-mortem" in txt:
        return "human_pm"
    if "synaptosome" in txt or "brain proteome" in txt:
        return "animal_brain"
    if "cell line" in sp or "fibroblast" in sp or "mefs" in txt:
        return "cells"
    if "neuron" in txt:
        return "neurons"
    if sp == "human":
        return "people"
    return "?"


def source_of(r):
    """EN/PT: (surname, year) of the primary source of a row; empty for the declared gap | fonte primaria da linha."""
    if not r["first_author"]:
        return ("", "")
    sur = r["first_author"].split()[0]
    return (SURNAME.get(sur, sur), r["year"])


def classes(rows):
    """EN/PT: (system of each param_id, block of each param_id), after the checks described above | classes apos as conferencias."""
    sys_of, blk_of = _invert(SYSTEM_OF), _invert(BLOCK_OF)
    ids = [r["param_id"] for r in rows]
    for name, m in (("system", sys_of), ("block", blk_of)):
        missing, extra = sorted(set(ids) - set(m)), sorted(set(m) - set(ids))
        if missing or extra:
            sys.exit(f"EN/PT: provenance {name} classes do not match the table; missing {missing}, unknown {extra}")
    bad = [f"{r['param_id']}: declared {sys_of[r['param_id']]}, table text says {_system_from_text(r)}"
           for r in rows if _system_from_text(r) != sys_of[r["param_id"]]]
    if bad:
        sys.exit("EN/PT: provenance system disagrees with the table text:\n  " + "\n  ".join(bad))
    return sys_of, blk_of


def table(rows):
    """EN/PT: one dict per row of the kinetic table, with its source, system and block | um dict por linha."""
    sys_of, blk_of = classes(rows)
    out = []
    for r in rows:
        s = source_of(r)
        out.append({"param_id": r["param_id"], "kind": r["kind"], "axis": AXIS[r["axis"]],
                    "source": f"{s[0]} {s[1]}".strip(), "system": sys_of[r["param_id"]], "block": blk_of[r["param_id"]]})
    return out


def write_table(rows, path=OUT):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=FIELDS)
        wr.writeheader()
        wr.writerows(table(rows))
    return path
