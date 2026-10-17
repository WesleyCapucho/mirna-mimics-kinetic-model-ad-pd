#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Molecular renders for the figures of the mechanistic study (Article 2), drawn with
     PyMOL from the deposited coordinates archived in data/raw/structures_2026/ and
     styled to sit beside the data panels: the paper's palette, a transparent background,
     and one light set for all of them. Four molecules, six images:

         a2_ago2_guide_target     Argonaute2 holding a guide, seed paired to a target (6N4O)
         a2_bace1_inhibitor       BACE1 with an inhibitor in the active site (4D8C)
         a2_abeta42_fibril_end    Abeta42 fibril, side view, growing end up (5OQV)
         a2_abeta42_fibril_axis   the same fibril seen down its axis
         a2_asyn_fibril_end       alpha-synuclein fibril, side view, growing end up (6CU7)
         a2_asyn_fibril_axis      the same fibril seen down its axis

     What they are for, and what they are not. They illustrate the molecules at the steps
     of the two models. They test no hypothesis and they are not a result. Two things in
     them are constructed and are said so wherever the images are used:
       - the translucent layer that continues each fibril past its last chain is not in
         the deposition. It is the last chain moved by the very transformation (rise and
         twist) that relates the last two chains of the deposition, one per protofilament,
         so it shows where the next layer of the same fibril would sit. The rms deviation
         of that transformation on the C-alpha atoms and the rise per layer are written to
         results/tables/structure_renders_article2.json;
       - the Abeta42 and alpha-synuclein monomers are disordered in solution and have no
         single structure, so no monomer is drawn.
     Every title, method, resolution and citation quoted with the images is read from the
     coordinate files by scripts/13 (results/tables/structure_figure_provenance.json) and
     re-read here; the script stops if the two disagree.

PT | Renders moleculares para as figuras do estudo mecanistico (Artigo 2), feitos no PyMOL a
     partir das coordenadas depositadas e arquivadas em data/raw/structures_2026/, e
     estilizados para ficar ao lado dos paineis de dados: a paleta do artigo, fundo
     transparente e uma so luz para todos. Quatro moleculas, seis imagens (lista acima).

     Para que servem e para que nao servem. Ilustram as moleculas nos passos dos dois
     modelos. Nao testam hipotese e nao sao resultado. Duas coisas nelas sao construidas e
     isso e dito onde quer que as imagens sejam usadas:
       - a camada translucida que continua cada fibrila alem da ultima cadeia nao esta no
         deposito. E a ultima cadeia movida pela mesma transformacao (subida e torcao) que
         relaciona as duas ultimas cadeias do deposito, uma por protofilamento, e mostra
         onde ficaria a proxima camada da mesma fibrila. O desvio quadratico medio dessa
         transformacao nos atomos C-alfa e a subida por camada ficam em
         results/tables/structure_renders_article2.json;
       - os monomeros de Abeta42 e de alfa-sinucleina sao desordenados em solucao e nao tem
         estrutura unica, entao nenhum monomero e desenhado.
     Todo titulo, metodo, resolucao e citacao que acompanha as imagens e lido dos arquivos
     de coordenadas pelo scripts/13 (results/tables/structure_figure_provenance.json) e
     relido aqui; o script para se os dois discordam.

     Needs PyMOL (pip install pymol-open-source; the wheel targets the system Python 3.13):
         python3.13 scripts/39_structure_renders_article2.py [--preview]
"""

import importlib.util
import json
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STRUCT_DIR = "data/raw/structures_2026"
FIG_DIR = "results/figures/structures"
TAB_DIR = "results/tables"
PROV13 = f"{TAB_DIR}/structure_figure_provenance.json"
OUT_JSON = f"{TAB_DIR}/structure_renders_article2.json"

PREVIEW = "--preview" in sys.argv
ONLY = ([a.split("=", 1)[1] for a in sys.argv if a.startswith("--only=")] or [""])[0].split(",")
ONLY = [x for x in ONLY if x]
W, H = (900, 675) if PREVIEW else (2400, 1800)

# EN/PT: the paper's palette (scripts/_viz_style.py), written out here because the figure
# style module imports matplotlib and this script runs in the interpreter that has PyMOL
AD, PD = "#2a78d6", "#eb6834"
AD_RAMP = ["#cde2fb", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#124785", "#0d366b"]    # pale -> deep
PD_RAMP = ["#fbe0d2", "#f5b496", "#f08b5f", "#eb6834", "#c04a1c", "#963612", "#6e2509"]
GUIDE, TARGET, SEED = "#1a2456", "#4f8fe3", "#d62839"
SURF_AGO2, SURF_BACE1 = "#aab8d4", "#dce9fa"
LIGAND_C = "#2b2b2b"
CUT_AWAY = 6.0          # EN/PT: angstrom in front of the BACE1 ligand centre where the front clipping plane sits


def hex_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def load_module(name, filename):
    path = os.path.join(HERE, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def style(cmd, shadows=False):
    """EN/PT: one light and one finish for every image | uma luz e um acabamento para todas."""
    cmd.reinitialize("settings")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.40)
    cmd.set("direct", 0.50)
    cmd.set("reflect", 0.40)
    cmd.set("specular", 0.35)
    cmd.set("shininess", 40)
    cmd.set("light_count", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_flat_sheets", 0)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.55)
    cmd.set("ray_trace_mode", 0)
    cmd.set("ray_shadow", 1 if shadows else 0)
    cmd.set("surface_quality", 1)
    cmd.set("orthoscopic", 0)


def define_colours(cmd):
    for name, h in (("guideC", GUIDE), ("targetC", TARGET), ("seedC", SEED), ("ago2S", SURF_AGO2),
                    ("bace1S", SURF_BACE1), ("bace1C", AD), ("ligC", LIGAND_C), ("dyadC", SEED)):
        cmd.set_color(name, hex_rgb(h))
    for i, h in enumerate(AD_RAMP):
        cmd.set_color(f"ad{i}", hex_rgb(h))
    for i, h in enumerate(PD_RAMP):
        cmd.set_color(f"pd{i}", hex_rgb(h))


def crop(path, pad=18):
    """EN/PT: trim the transparent margin and keep a small pad | apara a margem transparente."""
    from PIL import Image
    im = Image.open(path).convert("RGBA")
    bb = im.getchannel("A").point(lambda a: 255 if a > 8 else 0).getbbox()
    if bb:
        l, t, r, b = bb
        im = im.crop((max(0, l - pad), max(0, t - pad), min(im.width, r + pad), min(im.height, b + pad)))
    im.save(path)
    return im.size


def render(cmd, name, shadows=False):
    out = f"{FIG_DIR}/{name}.png"
    t0 = time.time()
    cmd.png(out, width=W, height=H, dpi=300, ray=1)
    size = crop(out)
    print(f"  {name}.png  {size[0]} x {size[1]} px  ({time.time() - t0:.0f} s)")
    return out, size


# ----------------------------------------------------------------------------- AGO2
def render_ago2(cmd):
    cmd.delete("all")
    cmd.load(f"{STRUCT_DIR}/6N4O.pdb", "ago2")
    style(cmd)
    define_colours(cmd)
    cmd.hide("everything")
    cmd.show("surface", "ago2 and chain A and polymer.protein")
    cmd.color("ago2S", "ago2 and chain A")
    cmd.set("transparency", 0.62, "ago2 and chain A")
    rna = "ago2 and chain C+D"
    cmd.show("cartoon", rna)
    cmd.show("sticks", rna)
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("stick_radius", 0.22, rna)
    cmd.color("guideC", "ago2 and chain C")
    cmd.color("targetC", "ago2 and chain D")
    cmd.color("seedC", "ago2 and chain C and resi 2-8")
    cmd.orient("ago2 and chain C+D")
    cmd.turn("y", 12)
    cmd.zoom("ago2", 0.5, complete=1)
    return render(cmd, "a2_ago2_guide_target")


# ----------------------------------------------------------------------------- BACE1
def render_bace1(cmd, prov):
    cmd.delete("all")
    cmd.load(f"{STRUCT_DIR}/4D8C.cif", "asu")
    style(cmd)
    define_colours(cmd)
    cmd.create("bace1", "asu and chain A and (polymer.protein or resn BXD)")
    cmd.delete("asu")
    cmd.hide("everything")
    sel = "bace1 and chain A"
    lig = f"{sel} and resn BXD"
    cmd.show("surface", f"{sel} and polymer.protein")
    cmd.color("bace1S", f"{sel} and polymer.protein")
    cmd.set("transparency", 0.78, f"{sel} and polymer.protein")
    cmd.show("cartoon", f"{sel} and polymer.protein")
    cmd.color("bace1C", f"{sel} and polymer.protein")
    cmd.set("cartoon_oval_length", 1.6)
    cmd.show("sticks", lig)
    cmd.show("spheres", lig)
    cmd.set("sphere_scale", 0.36, lig)
    cmd.color("ligC", f"{lig} and elem C")
    cmd.set("stick_radius", 0.30, lig)
    # EN/PT: the catalytic dyad, found by motif and by proximity as in scripts/13, and
    #        required to agree with the dyad that script recorded
    ca = {"l": []}
    cmd.iterate(f"{sel} and polymer.protein and name CA", "l.append((resi, oneletter))", space=ca)
    seq = "".join(x[1] for x in ca["l"])
    by_motif = {}
    for motif in ("DTGS", "DSGT"):
        hits = [m.start() for m in re.finditer(motif, seq)]
        if len(hits) != 1:
            sys.exit(f"EN/PT: BACE1 motif {motif} found {len(hits)} times in chain A")
        by_motif[motif] = int(ca["l"][hits[0]][0])
    if by_motif != prov["4D8C"]["catalytic_dyad_resi_in_this_deposition"]:
        sys.exit(f"EN/PT: BACE1 dyad {by_motif} differs from scripts/13: "
                 f"{prov['4D8C']['catalytic_dyad_resi_in_this_deposition']}")
    cmd.select("dyad", f"{sel} and polymer.protein and resi {'+'.join(str(r) for r in by_motif.values())} and resn ASP")
    cmd.show("sticks", "dyad and not hydro")
    cmd.show("spheres", "dyad and not hydro")
    cmd.set("sphere_scale", 0.34, "dyad")
    cmd.color("dyadC", "dyad and elem C")
    cmd.set("stick_radius", 0.30, "dyad")
    cmd.orient(lig)
    cmd.zoom(sel, 0.5, complete=1)
    v = list(cmd.get_view())
    cam = -v[11]
    # EN/PT: cut-away: the front clipping plane sits CUT_AWAY angstrom in front of the ligand
    #        centre, so the protein between the camera and the active site is not drawn
    v[15] = max(1.0, cam - CUT_AWAY)
    v[16] = cam + 150.0
    cmd.set_view(v)
    out = render(cmd, "a2_bace1_inhibitor")
    cmd.delete("dyad")
    return out, by_motif


# ----------------------------------------------------------------------------- fibrils
def kabsch(P, Q):
    """EN/PT: rotation R and translation t with P @ R + t ~ Q (rows are points) | superposicao de Kabsch."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    R = U @ D @ Vt
    t = qc - pc @ R
    rms = float(np.sqrt(np.mean(np.sum((P @ R + t - Q) ** 2, axis=1))))
    return R, t, rms


def protofilaments(cmd, obj, cutoff=6.0):
    """EN/PT: the same rule as scripts/13: chains whose CA centroids are closer than the cutoff are one stack."""
    chains = cmd.get_chains(obj)
    cen = {c: np.array(cmd.get_coords(f"{obj} and chain {c} and name CA")).mean(axis=0) for c in chains}
    groups, seen = [], set()
    for c in chains:
        if c in seen:
            continue
        comp, todo = [], [c]
        while todo:
            x = todo.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            todo += [y for y in chains if y not in seen and np.linalg.norm(cen[x] - cen[y]) < cutoff]
        groups.append(sorted(comp))
    if len(groups) != 2:
        sys.exit(f"EN/PT: expected two protofilaments in {obj}, found {groups}")
    return groups, cen


def ca_by_resi(cmd, sel):
    d = {"l": []}
    cmd.iterate_state(1, f"{sel} and name CA", "l.append((resi, x, y, z))", space=d)
    return {r: np.array([x, y, z]) for r, x, y, z in d["l"]}


def build_fibril(cmd, pdb_id, name, ramp, prov):
    """EN/PT: load the fibril, order the chains along the axis, add the next layer, colour by layer."""
    cmd.delete("all")
    cmd.load(f"{STRUCT_DIR}/{pdb_id}.cif", "fib")
    style(cmd)
    define_colours(cmd)
    groups, cen = protofilaments(cmd, "fib")
    if [groups[0], groups[1]] != [prov[pdb_id]["chain_split"]["protofilament_1"], prov[pdb_id]["chain_split"]["protofilament_2"]]:
        sys.exit(f"EN/PT: {pdb_id} protofilaments {groups} differ from scripts/13")
    # EN/PT: the fibril axis is the mean of the stack directions of the two protofilaments;
    #        the growing end is the + end of that axis, and chains are ordered along it
    axes = []
    for g in groups:
        C = np.array([cen[c] for c in g])
        u, sv, vt = np.linalg.svd(C - C.mean(axis=0))
        axes.append(vt[0])
    if np.dot(axes[0], axes[1]) < 0:
        axes[1] = -axes[1]
    axis = axes[0] + axes[1]
    axis /= np.linalg.norm(axis)
    ordered = [sorted(g, key=lambda c: float(np.dot(cen[c], axis))) for g in groups]
    # EN/PT: next layer of each protofilament, from the transformation that carries the
    #        second-to-last chain onto the last one
    info = []
    for k, order in enumerate(ordered):
        prev, last = order[-2], order[-1]
        cp, cl = ca_by_resi(cmd, f"fib and chain {prev}"), ca_by_resi(cmd, f"fib and chain {last}")
        common = sorted(set(cp) & set(cl), key=int)
        if len(common) < 20:
            sys.exit(f"EN/PT: {pdb_id} chains {prev}/{last} share only {len(common)} CA positions")
        P = np.array([cp[r] for r in common])
        Q = np.array([cl[r] for r in common])
        R, t, rms = kabsch(P, Q)
        # rise along the axis and twist about it, from the same transformation
        rise = float(np.dot(Q.mean(axis=0) - P.mean(axis=0), axis))
        ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))
        cmd.create(f"ghost{k}", f"fib and chain {last}")
        xyz = np.array(cmd.get_coords(f"ghost{k}"))
        cmd.load_coords(np.ascontiguousarray(xyz @ R + t), f"ghost{k}")
        info.append(dict(protofilament=k + 1, previous_chain=prev, last_chain=last, ca_pairs=len(common),
                         rms_ca_angstrom=round(rms, 2), rise_angstrom=round(rise, 2), rotation_deg=round(ang, 2),
                         chains_along_axis=order))
        if rms > 2.0:
            sys.exit(f"EN/PT: {pdb_id} layer transformation has rms {rms:.2f} A over {len(common)} CA atoms; "
                     "the chains are not equivalent enough to extrapolate one layer")
    cmd.hide("everything")
    cmd.show("cartoon", "fib or ghost0 or ghost1")
    n = max(len(o) for o in ordered)
    for k, order in enumerate(ordered):
        for i, c in enumerate(order):
            level = (i + (n - len(order))) / max(1, n - 1)           # 0 far end .. 1 growing end
            idx = min(len(ramp) - 1, int(round(level * (len(ramp) - 1))))
            if k == 1:
                idx = max(0, idx - 1)                                # second protofilament one step paler
            cmd.color(f"{'ad' if ramp is AD_RAMP else 'pd'}{idx}", f"fib and chain {c}")
        cmd.color(f"{'ad' if ramp is AD_RAMP else 'pd'}{len(ramp) - 1 - k}", f"ghost{k}")
    cmd.show("sticks", "fib and sidechain and not hydro")
    cmd.set("stick_radius", 0.12, "fib")
    cmd.show("sticks", "(ghost0 or ghost1) and sidechain and not hydro")
    cmd.set("stick_radius", 0.12, "ghost0 or ghost1")
    cmd.set("cartoon_transparency", 0.72, "ghost0 or ghost1")
    cmd.set("stick_transparency", 0.72, "ghost0 or ghost1")
    return axis, ordered, info, cen


def set_axis_view(cmd, axis, side_hint, sel="fib or ghost0 or ghost1"):
    """EN/PT: camera with the fibril axis along the screen's y (growing end up) and the two protofilaments side by side."""
    a = axis / np.linalg.norm(axis)
    u = side_hint - np.dot(side_hint, a) * a
    u /= np.linalg.norm(u)
    z = np.cross(u, a)
    R = np.column_stack([u, a, z])                     # columns: screen x, y, z in model space
    cmd.orient(sel)
    v = list(cmd.get_view())
    v[0:9] = R.reshape(9).tolist()
    cmd.set_view(v)
    cmd.zoom(sel, 0.5, complete=1)
    v = list(cmd.get_view())
    cam = -v[11]
    v[15] = max(1.0, cam - 200.0)
    v[16] = cam + 200.0
    cmd.set_view(v)


def render_fibril(cmd, pdb_id, short, ramp, prov):
    axis, ordered, info, cen = build_fibril(cmd, pdb_id, short, ramp, prov)
    c1 = np.mean([cen[c] for c in ordered[0]], axis=0)
    c2 = np.mean([cen[c] for c in ordered[1]], axis=0)
    # side view: axis vertical with the growing end up, the two protofilaments side by side
    set_axis_view(cmd, axis, c2 - c1)
    side = render(cmd, f"a2_{short}_fibril_end")
    # axis view: look down the axis from the growing end, ghost layer removed
    cmd.hide("everything", "ghost0 or ghost1")
    a = axis / np.linalg.norm(axis)
    up = (c2 - c1) - np.dot(c2 - c1, a) * a
    up /= np.linalg.norm(up)
    R = np.column_stack([up, np.cross(a, up), a])      # screen x = protofilament separation, z = + end of the axis (camera on that end)
    cmd.orient("fib")
    v = list(cmd.get_view())
    v[0:9] = R.reshape(9).tolist()
    cmd.set_view(v)
    cmd.zoom("fib", 0.5, complete=1)
    v = list(cmd.get_view())
    cam = -v[11]
    v[15] = max(1.0, cam - 200.0)
    v[16] = cam + 200.0
    cmd.set_view(v)
    axis_img = render(cmd, f"a2_{short}_fibril_axis")
    return side, axis_img, info


# ----------------------------------------------------------------------------- main
def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    s13 = load_module("s13", "13_structure_figures.py")
    if not os.path.exists(PROV13):
        sys.exit(f"EN/PT: missing {PROV13}; run scripts/13 first")
    prov = json.load(open(PROV13, encoding="utf-8"))
    for pdb_id, spec in s13.EXPECT.items():                       # EN/PT: re-read, must agree
        p = s13.read_provenance(pdb_id, os.path.join(STRUCT_DIR, spec["file"]))
        for k in ("title", "method", "resolution_angstrom", "doi", "pmid"):
            if p.get(k) != prov[pdb_id].get(k):
                sys.exit(f"EN/PT: {pdb_id} {k}: file {p.get(k)!r} vs scripts/13 {prov[pdb_id].get(k)!r}")

    import pymol
    pymol.finish_launching(["pymol", "-qc"])
    from pymol import cmd

    print("=" * 78)
    print("EN | Structure renders for Article 2 | PT | Renders de estrutura para o Artigo 2")
    print("=" * 78)
    out = {"_note_en": "Illustration of the molecules at the steps of the models; not a result. The translucent "
                       "layer past the last chain of each fibril is constructed (see the script docstring).",
           "_note_pt": "Ilustracao das moleculas nos passos dos modelos; nao e resultado. A camada translucida "
                       "alem da ultima cadeia de cada fibrila e construida (ver o docstring do script).",
           "image_px": {"width": W, "height": H}, "renders": {}}
    if ONLY and os.path.exists(OUT_JSON):               # EN/PT: partial run keeps the other entries
        out = json.load(open(OUT_JSON, encoding="utf-8"))
        out["image_px"] = {"width": W, "height": H}

    def want(name):
        return not ONLY or name in ONLY

    if want("ago2"):
        img, size = render_ago2(cmd)
        out["renders"]["a2_ago2_guide_target"] = dict(pdb_id="6N4O", crop_px=list(size),
                                                      guide_chain="C", target_chain="D", seed="guide nucleotides 2-8")
    if want("bace1"):
        (img, size), dyad = render_bace1(cmd, prov)
        out["renders"]["a2_bace1_inhibitor"] = dict(pdb_id="4D8C", crop_px=list(size), chain="A",
                                                    catalytic_dyad_resi=dyad)
    for pdb_id, short, ramp in (("5OQV", "abeta42", AD_RAMP), ("6CU7", "asyn", PD_RAMP)):
        if want(short):
            side, axis_img, info = render_fibril(cmd, pdb_id, short, ramp, prov)
            out["renders"][f"a2_{short}_fibril_end"] = dict(pdb_id=pdb_id, crop_px=list(side[1]), layer_transformation=info)
            out["renders"][f"a2_{short}_fibril_axis"] = dict(pdb_id=pdb_id, crop_px=list(axis_img[1]))
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"EN/PT -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
