#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Render the structural figures for the two modelled axes, from deposited
     coordinates archived in data/raw/structures_2026/.
PT | Renderiza as figuras estruturais dos dois eixos modelados, a partir de
     coordenadas depositadas e arquivadas em data/raw/structures_2026/.

EN | What these figures are for, and what they are not. They illustrate the
     molecular machinery behind the miR-29/BACE1/Abeta and miR-7/SNCA/alpha-synuclein
     axes: how a microRNA recognises a target at all, what the enzyme the miR-29
     family represses looks like, and what the aggregated end states of both axes
     look like. They are mechanistic illustration for the introduction and
     discussion. They test no hypothesis and they are not a result. The quantitative
     claims in this project rest on the meta-analysis and on scripts/12.
PT | Para que servem estas figuras, e para que nao servem. Elas ilustram a maquinaria
     molecular por tras dos eixos miR-29/BACE1/Abeta e miR-7/SNCA/alfa-sinucleina:
     como um microRNA reconhece um alvo, qual a cara da enzima que a familia miR-29
     reprime, e qual a cara dos estados agregados finais dos dois eixos. Sao ilustracao
     mecanistica para a introducao e a discussao. Nao testam hipotese e nao sao
     resultado. As afirmacoes quantitativas deste projeto repousam na meta-analise e
     no scripts/12.

EN | Provenance is read from the coordinate files themselves and written to
     results/tables/structure_figure_provenance.json, so the figure legends quote the
     deposited title, method and resolution rather than a description from memory.
PT | A proveniencia e lida dos proprios arquivos de coordenadas e escrita em
     results/tables/structure_figure_provenance.json, para que as legendas citem o
     titulo, o metodo e a resolucao depositados, e nao uma descricao de memoria.

    python scripts/13_structure_figures.py
"""

import json
import os
import re
import sys

STRUCT_DIR = "data/raw/structures_2026"
FIG_DIR = "results/figures/structures"
TAB_DIR = "results/tables"

# EN/PT: what each entry must be, checked against the file before it is drawn
EXPECT = {
    "6N4O": dict(file="6N4O.pdb", must_match=r"argonaute",
                 role_en="How a microRNA recognises its target at all: the guide is held by "
                         "Argonaute2 and the seed pairs with the target. Common to both axes.",
                 role_pt="Como um microRNA reconhece o alvo: o guia e segurado pelo Argonauta2 "
                         "e a semente pareia com o alvo. Comum aos dois eixos."),
    "4D8C": dict(file="4D8C.cif", must_match=r"secretase",
                 role_en="BACE1, the enzyme the miR-29 family represses, with an inhibitor bound "
                         "in the active site.",
                 role_pt="BACE1, a enzima que a familia miR-29 reprime, com inibidor ligado no "
                         "sitio ativo."),
    "6CU7": dict(file="6CU7.cif", must_match=r"synuclein",
                 role_en="The aggregated end state of the miR-7/SNCA axis.",
                 role_pt="O estado agregado final do eixo miR-7/SNCA.",
                 # EN/PT: checked against PubMed (PMID 30190461) on 2026-09-23
                 note_en="The deposition reports 3.5 A (FSC 0.5 cut-off); the article "
                         "abstract reports 3.7 A for both polymorphs. Quote the criterion "
                         "with the number.",
                 note_pt="O deposito informa 3,5 A (corte FSC 0,5); o resumo do artigo "
                         "informa 3,7 A para os dois polimorfos. Citar o criterio junto "
                         "com o numero."),
    "5OQV": dict(file="5OQV.cif", must_match=r"amyloid",
                 role_en="The aggregated end state of the miR-29/BACE1/Abeta axis.",
                 role_pt="O estado agregado final do eixo miR-29/BACE1/Abeta."),
}


def read_provenance(pdb_id, path):
    """
    EN | Pull the deposited title, method, resolution and citation out of the file.
         A figure legend written from memory is a figure legend that can be wrong;
         this makes the legend quote the deposition.
    PT | Extrai titulo, metodo, resolucao e citacao depositados do proprio arquivo.
         Legenda escrita de memoria e legenda que pode estar errada; assim a legenda
         cita o deposito.
    """
    text = open(path, encoding="utf-8", errors="replace").read()
    p = dict(pdb_id=pdb_id, file=os.path.basename(path))

    if path.endswith(".pdb"):
        title = " ".join(re.findall(r"^TITLE\s*\d*\s*(.*?)\s*$", text, re.M))
        p["title"] = re.sub(r"\s+", " ", title).strip()
        m = re.search(r"^EXPDTA\s+(.*?)\s*$", text, re.M)
        p["method"] = m.group(1).strip() if m else ""
        m = re.search(r"^REMARK\s+2\s+RESOLUTION\.\s*([\d.]+)", text, re.M)
        p["resolution_angstrom"] = float(m.group(1)) if m else None
        m = re.search(r"^JRNL\s+DOI\s+(.*?)\s*$", text, re.M)
        p["doi"] = m.group(1).strip().lower() if m else ""
        m = re.search(r"^JRNL\s+PMID\s+(\d+)", text, re.M)
        p["pmid"] = m.group(1) if m else ""
    else:
        m = re.search(r"_struct\.title\s+(?:'([^']*)'|\"([^\"]*)\"|;(.*?)\n;)", text, re.S)
        p["title"] = re.sub(r"\s+", " ", (m.group(1) or m.group(2) or m.group(3) or "")).strip() if m else ""
        m = re.search(r"_exptl\.method\s+(?:'([^']*)'|\"([^\"]*)\"|(\S+))", text)
        p["method"] = ((m.group(1) or m.group(2) or m.group(3)) if m else "").strip()
        res = None
        m = re.search(r"_em_3d_reconstruction\.resolution\s+([\d.]+)", text)
        if m:
            res = float(m.group(1))
        else:
            m = re.search(r"_refine\.ls_d_res_high\s+([\d.]+)", text)
            if m:
                res = float(m.group(1))
        p["resolution_angstrom"] = res
        m = re.search(r"_em_3d_reconstruction\.resolution_method\s+(?:'([^']*)'|(\S+))", text)
        if m:
            p["resolution_method"] = (m.group(1) or m.group(2)).strip()
        cit = parse_cif_citation(text)
        p["doi"] = cit.get("pdbx_database_id_DOI", "")
        p["pmid"] = cit.get("pdbx_database_id_PubMed", "")
    return p


def parse_cif_citation(text):
    """
    EN | Read the primary citation from an mmCIF file, handling both the key-value
         form and the loop_ form, and skipping the PDB deposition DOI.
    PT | Le a citacao primaria de um mmCIF, tratando tanto a forma chave-valor quanto a
         forma loop_, e ignorando o DOI de deposito do PDB.

    EN | Two traps here, both of which produced a wrong DOI in this project before this
         function existed. First, in a loop_ block the tag is a column header, so
         reading the token after it returns the first data row ("primary"). Second,
         matching the first DOI-shaped string in the file returns 10.2210/pdb<id>/pdb -
         the DOI of the DEPOSITION, not of the paper that reports the structure. A
         legend citing the deposition DOI instead of the article looks correct and
         is not.
    PT | Duas armadilhas, e as duas produziram DOI errado neste projeto antes desta
         funcao existir. Primeiro, num bloco loop_ a tag e cabecalho de coluna, entao
         ler o token seguinte devolve a primeira linha de dados ("primary"). Segundo,
         casar a primeira string com cara de DOI no arquivo devolve
         10.2210/pdb<id>/pdb - o DOI do DEPOSITO, nao do artigo que reporta a
         estrutura. Uma legenda citando o DOI do deposito parece certa e nao e.
    """
    def usable(doi):
        d = (doi or "").strip().strip("'\"").lower()
        if not d.startswith("10.") or d in ("?", "."):
            return ""
        return "" if d.startswith("10.2210/pdb") else d.rstrip(".,")

    out = {}
    # --- loop_ form ------------------------------------------------------
    for m in re.finditer(r"^loop_\s*\n((?:\s*_citation\.[^\n]*\n)+)", text, re.M):
        tags = re.findall(r"_citation\.(\S+)", m.group(1))
        body = text[m.end():]
        end = re.search(r"^\s*(?:#|loop_|_\w)", body, re.M)
        body = body[:end.start()] if end else body
        # EN/PT: tokenise, keeping ;multi-line; and quoted values whole
        toks = re.findall(r";(?:[^\n]*\n)*?;|'[^']*'|\"[^\"]*\"|\S+", body)
        if len(toks) >= len(tags):
            row = dict(zip(tags, toks[:len(tags)]))
            for k in ("pdbx_database_id_DOI", "pdbx_database_id_PubMed"):
                v = row.get(k, "").strip().strip("'\"")
                if k.endswith("DOI"):
                    v = usable(v)
                elif not re.fullmatch(r"\d{6,9}", v):
                    v = ""
                if v:
                    out[k] = v
    # --- key-value form --------------------------------------------------
    if "pdbx_database_id_DOI" not in out:
        m = re.search(r"_citation\.pdbx_database_id_DOI\s+(\S+)", text)
        if m:
            v = usable(m.group(1))
            if v:
                out["pdbx_database_id_DOI"] = v
    if "pdbx_database_id_PubMed" not in out:
        m = re.search(r"_citation\.pdbx_database_id_PubMed\s+(\d{6,9})\b", text)
        if m:
            out["pdbx_database_id_PubMed"] = m.group(1)
    return out


def open_clip(cmd, margin=150.0):
    """
    EN | Put the front and rear clipping planes well outside the molecule, set from
         the camera distance rather than nudged with relative moves (whose sign is
         easy to get backwards - doing so once left a thin slab through BACE1).
    PT | Coloca os planos de corte frontal e traseiro bem fora da molecula, definidos a
         partir da distancia da camera e nao por deslocamentos relativos (cujo sinal e
         facil de inverter - isso ja deixou uma fatia fina atravessando a BACE1).
    """
    v = list(cmd.get_view())
    cam = -v[11]
    v[15] = max(1.0, cam - margin)
    v[16] = cam + margin
    cmd.set_view(v)


def style_common(cmd):
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("ambient", 0.22)
    cmd.set("specular", 0.25)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("surface_quality", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)


def render_ago2(cmd, path, out):
    """EN/PT: Ago2 holding the guide, seed paired to the target."""
    cmd.delete("all"); cmd.load(path, "ago2"); style_common(cmd)
    cmd.hide("everything")
    cmd.show("surface", "ago2 and chain A and polymer.protein")
    cmd.color("grey80", "ago2 and chain A")
    cmd.set("transparency", 0.55, "ago2 and chain A")
    # guide (chain C) and target (chain D)
    cmd.show("cartoon", "ago2 and chain C+D")
    cmd.show("sticks", "ago2 and chain C+D")
    cmd.set("cartoon_ring_mode", 3, "ago2 and chain C+D")
    cmd.set("cartoon_ring_finder", 1, "ago2 and chain C+D")
    cmd.color("orange", "ago2 and chain C")
    cmd.color("skyblue", "ago2 and chain D")
    # EN/PT: the seed, guide nucleotides 2-8
    cmd.color("firebrick", "ago2 and chain C and resi 2-8")
    cmd.set("cartoon_nucleic_acid_mode", 4)
    # EN | Orient on the RNA to get a readable view of the duplex, but FRAME on the
    #      whole complex. An earlier version zoomed to the RNA with a small buffer, so
    #      the protein ran off the top of the canvas and left dead space below.
    # PT | Orienta pelo RNA para obter uma vista legivel do duplex, mas ENQUADRA pelo
    #      complexo inteiro. Uma versao anterior dava zoom no RNA com margem pequena, e
    #      a proteina saia pelo topo da tela deixando espaco morto embaixo.
    cmd.orient("ago2 and chain C+D")
    cmd.turn("y", 12)
    cmd.zoom("ago2", 3, complete=1)
    cmd.png(out, width=2000, height=1500, dpi=300, ray=1)


def render_bace1(cmd, path, out):
    """EN/PT: BACE1 with the co-crystallised inhibitor and the catalytic aspartates."""
    cmd.delete("all"); cmd.load(path, "asu"); style_common(cmd)
    # EN | Chain A and its inhibitor go into their own object. Left inside the
    #      three-copy asymmetric unit, the neighbouring copies take part in the
    #      surface calculation and leave a grey cut patch where they touch chain A.
    # PT | A cadeia A e seu inibidor vao para um objeto proprio. Dentro da unidade
    #      assimetrica com tres copias, as copias vizinhas entram no calculo da
    #      superficie e deixam uma placa cinza cortada onde encostam na cadeia A.
    cmd.create("bace1", "asu and chain A and (polymer.protein or resn BXD)")
    cmd.delete("asu")
    cmd.hide("everything")
    sel = "bace1 and chain A"
    # EN | One flat colour for cartoon and surface. An earlier version ran a spectrum
    #      over the CA atoms, and the surface inherited it as a pink/cyan halo that
    #      carried no information at all.
    # PT | Uma cor chapada para cartoon e superficie. Uma versao anterior passava um
    #      espectro nos atomos CA, e a superficie herdava aquilo como um halo rosa/ciano
    #      que nao carregava informacao nenhuma.
    cmd.show("cartoon", f"{sel} and polymer.protein")
    cmd.color("skyblue", f"{sel} and polymer.protein")
    cmd.show("surface", f"{sel} and polymer.protein")
    cmd.set("transparency", 0.72, sel)
    lig = f"{sel} and resn BXD"
    cmd.show("sticks", lig); cmd.color("orange", f"{lig} and elem C")
    cmd.set("stick_radius", 0.22, lig)
    # EN | The catalytic aspartyl dyad, identified TWICE and required to agree: by
    #      sequence motif (the Asp of DTGS and of DSGT, the two catalytic motifs of
    #      pepsin-family aspartic proteases) and by proximity to the inhibitor. Residue
    #      numbers are NOT hard-coded: in this deposition the DSGT aspartate is
    #      numbered 217, and BACE1 depositions do not all share one numbering, so a
    #      legend must quote the motif together with the number read from the file.
    # PT | A diade catalitica de aspartatos, identificada DUAS vezes e exigindo
    #      concordancia: pelo motivo de sequencia (o Asp de DTGS e de DSGT, os dois
    #      motivos cataliticos das aspartil-proteases da familia da pepsina) e pela
    #      proximidade com o inibidor. Os numeros de residuo NAO sao fixos no codigo:
    #      neste deposito o aspartato do DSGT e o 217, e os depositos de BACE1 nao
    #      compartilham uma numeracao unica, entao a legenda deve citar o motivo junto
    #      com o numero lido do arquivo.
    ca = {"l": []}
    cmd.iterate(f"{sel} and polymer.protein and name CA", "l.append((resi, oneletter))", space=ca)
    seq = "".join(x[1] for x in ca["l"])
    by_motif = {}
    for motif in ("DTGS", "DSGT"):
        hits = [m.start() for m in re.finditer(motif, seq)]
        if len(hits) != 1:
            sys.exit(f"EN/PT: BACE1 motif {motif} found {len(hits)} times in chain A")
        by_motif[motif] = ca["l"][hits[0]][0]
    cmd.select("dyad", f"byres ({sel} and resn ASP and polymer.protein within 5 of ({lig}))")
    near = {"l": []}
    cmd.iterate("dyad and name CA", "l.append(resi)", space=near)
    if sorted(near["l"]) != sorted(by_motif.values()):
        sys.exit(f"EN/PT: dyad by motif {by_motif} != dyad by proximity {near['l']}")
    cmd.show("sticks", "dyad and not hydro"); cmd.util.cnc("dyad")
    cmd.color("firebrick", "dyad and elem C")
    cmd.set("stick_radius", 0.22, "dyad")
    # EN | Orient on the pocket, frame the whole domain, and open the clipping planes.
    #      Without the last step the front plane slices the surface and the cut face
    #      renders as a grey slab with coloured edges at the bottom of the figure.
    # PT | Orienta pelo bolso, enquadra o dominio inteiro e abre os planos de corte.
    #      Sem o ultimo passo o plano frontal fatia a superficie e a face cortada
    #      aparece como uma placa cinza de bordas coloridas na base da figura.
    cmd.orient(lig)
    cmd.zoom(sel, 3, complete=1)
    open_clip(cmd)
    cmd.png(out, width=2000, height=1500, dpi=300, ray=1)

    # EN/PT: close-up of the active site, no surface, dyad labelled from the file
    cmd.hide("surface", sel)
    cmd.set("cartoon_transparency", 0.8, sel)
    cmd.set("label_size", 26); cmd.set("label_color", "black")
    cmd.set("label_font_id", 7); cmd.set("float_labels", 1)
    for motif, resi in by_motif.items():
        cmd.label(f"dyad and resi {resi} and name CA", f'"Asp{resi} ({motif})"')
    # EN | Orient on the ligand atoms plus the two dyad CA atoms, so the ligand is seen
    #      along its length and the dyad sits beside it. Orienting on the whole pocket
    #      put the dyad straight behind the ligand, with the two labels on top of each
    #      other; orienting on the ligand centre and the two CAs alone showed the
    #      ligand end-on.
    # PT | Orienta pelos atomos do ligante mais os dois CA da diade, para que o ligante
    #      apareca no comprimento e a diade fique ao lado. Orientar pelo bolso inteiro
    #      punha a diade logo atras do ligante, com os rotulos sobrepostos; orientar so
    #      pelo centro do ligante e pelos dois CA mostrava o ligante de ponta.
    cmd.orient(f"({lig}) or (dyad and name CA)")
    # EN | That view can still stack the two aspartates in depth. Rotate about the
    #      screen's horizontal axis (which keeps the ligand's length in view) to the
    #      angle that puts the two dyad CA atoms farthest apart on screen.
    # PT | Essa vista ainda pode empilhar os dois aspartatos em profundidade. Gira em
    #      torno do eixo horizontal da tela (o que mantem o comprimento do ligante a
    #      vista) ate o angulo que deixa os dois CA da diade mais afastados na tela.
    import numpy as _np

    def _screen_sep():
        v = cmd.get_view()
        R = _np.array(v[:9]).reshape(3, 3)
        xyz = _np.array(cmd.get_coords("dyad and name CA")) - _np.array(v[12:15])
        cam = xyz @ R
        return float(_np.linalg.norm(cam[0, :2] - cam[1, :2]))

    best = max(range(0, 180, 5), key=lambda a: (cmd.turn("x", a), _screen_sep(),
                                                 cmd.turn("x", -a))[1])
    cmd.turn("x", best)
    # EN | The inhibitor sits BETWEEN the two aspartates, as this class of inhibitor
    #      binds, so from any angle it partly covers them. Push each label away from
    #      the ligand along the screen's vertical axis instead of hiding the ligand.
    # PT | O inibidor fica ENTRE os dois aspartatos, como essa classe de inibidor se
    #      liga, entao de qualquer angulo ele os cobre em parte. Cada rotulo e
    #      empurrado para longe do ligante no eixo vertical da tela, em vez de
    #      esconder o ligante.
    v = cmd.get_view()
    R = _np.array(v[:9]).reshape(3, 3)
    lig_y = float(((_np.array(cmd.get_coords(lig)) - _np.array(v[12:15])) @ R)[:, 1].mean())
    for motif, resi in by_motif.items():
        y = float(((_np.array(cmd.get_coords(f"dyad and resi {resi} and name CA"))
                    - _np.array(v[12:15])) @ R)[0, 1])
        cmd.set("label_position", (0.0, 5.0 if y >= lig_y else -5.0, 4.0),
                f"dyad and resi {resi} and name CA")
    cmd.zoom(f"({lig}) or dyad", 5, complete=1)
    open_clip(cmd)
    cmd.png(out.replace(".png", "_active_site.png"), width=2000, height=1500, dpi=300, ray=1)
    cmd.delete("dyad")
    return {m: int(r) for m, r in by_motif.items()}


def render_fibril(cmd, path, out, label):
    """EN/PT: fibril down the axis and from the side, coloured by protofilament."""
    cmd.delete("all"); cmd.load(path, "fib"); style_common(cmd)
    cmd.hide("everything")
    cmd.show("cartoon", "fib")
    cmd.set("cartoon_flat_sheets", 0)
    # EN | Colour by PROTOFILAMENT, not by chain. A rainbow over ten chains hides the
    #      one thing the view down the axis exists to show: that the fibril is two
    #      protofilaments packed against each other, each a stack of identical layers.
    # PT | Colorir por PROTOFILAMENTO, nao por cadeia. Um arco-iris sobre dez cadeias
    #      esconde justamente o que a vista pelo eixo existe para mostrar: que a fibrila
    #      sao dois protofilamentos encaixados, cada um uma pilha de camadas identicas.
    #
    # EN | Protofilaments are found from the structure itself. Consecutive layers of a
    #      cross-beta stack sit ~4.7-4.9 A apart (the inter-strand spacing), while
    #      chains in different protofilaments are tens of angstroms apart. Chains whose
    #      CA centroids are closer than STACK_CUTOFF are linked, and the connected
    #      components are the protofilaments. Two earlier versions were wrong: taking
    #      the first and second half of the chain list assumes a deposition order, and
    #      cutting along the principal axis of the centroids can cut across the stack
    #      instead of between protofilaments when the stack is longer than the gap.
    #      Note that a 9-chain deposition necessarily splits 5/4; that is the model's
    #      length, not a flaw.
    # PT | Os protofilamentos sao obtidos da propria estrutura. Camadas consecutivas de
    #      uma pilha cross-beta ficam a ~4,7-4,9 A (espacamento entre fitas), enquanto
    #      cadeias de protofilamentos diferentes ficam a dezenas de angstroms. Cadeias
    #      com centroides de CA mais proximos que STACK_CUTOFF sao ligadas, e os
    #      componentes conexos sao os protofilamentos. Duas versoes anteriores erravam:
    #      pegar a primeira e a segunda metade da lista supoe uma ordem de deposito, e
    #      cortar pelo eixo principal dos centroides pode cortar a pilha em vez de
    #      separar os protofilamentos quando a pilha e mais longa que a distancia entre
    #      eles. Um deposito com 9 cadeias necessariamente divide 5/4; e o tamanho do
    #      modelo, nao um defeito.
    import numpy as _np
    STACK_CUTOFF = 6.0  # angstrom
    chains = cmd.get_chains("fib")
    cen = {c: _np.array(cmd.get_coords(f"fib and chain {c} and name CA")).mean(axis=0)
           for c in chains}
    dist = {(a, b): float(_np.linalg.norm(cen[a] - cen[b])) for a in chains for b in chains}
    groups, seen = [], set()
    for c in chains:
        if c in seen:
            continue
        comp, todo = [], [c]
        while todo:
            x = todo.pop()
            if x in seen:
                continue
            seen.add(x); comp.append(x)
            todo += [y for y in chains if y not in seen and dist[(x, y)] < STACK_CUTOFF]
        groups.append(sorted(comp))
    if len(groups) != 2:
        sys.exit(f"EN/PT: {label}: expected 2 protofilaments, found {len(groups)}: {groups}")
    pf1, pf2 = groups
    stack = [dist[(a, b)] for g in groups for a in g for b in g
             if a < b and dist[(a, b)] < STACK_CUTOFF]
    between = min(dist[(a, b)] for a in pf1 for b in pf2)
    cmd.color("deepteal", "fib and chain " + "+".join(pf1))
    cmd.color("orange", "fib and chain " + "+".join(pf2))
    cmd.show("sticks", "fib and sidechain and not hydro")
    cmd.set("stick_radius", 0.13)
    cmd.set("cartoon_transparency", 0.1)
    cmd.orient("fib")
    cmd.zoom("fib", 4, complete=1)
    cmd.png(out.replace(".png", "_axis.png"), width=2000, height=1500, dpi=300, ray=1)
    cmd.turn("x", 90)
    cmd.zoom("fib", 4, complete=1)
    cmd.png(out.replace(".png", "_side.png"), width=2000, height=1500, dpi=300, ray=1)
    return dict(n_chains=len(chains), protofilament_1=pf1, protofilament_2=pf2,
                layer_spacing_angstrom=[round(min(stack), 2), round(max(stack), 2)],
                closest_interprotofilament_centroids_angstrom=round(between, 2),
                assignment=f"connected components of chains with CA-centroid distance < {STACK_CUTOFF} A")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(TAB_DIR, exist_ok=True)

    import pymol
    pymol.finish_launching(["pymol", "-qc"])
    from pymol import cmd

    prov = {}
    print("=" * 78)
    print("EN | Structure figures | PT | Figuras de estrutura")
    print("=" * 78)

    for pdb_id, spec in EXPECT.items():
        path = os.path.join(STRUCT_DIR, spec["file"])
        if not os.path.exists(path):
            sys.exit(f"EN/PT: missing {path}")
        p = read_provenance(pdb_id, path)
        # EN | Confirm the file is the molecule the figure claims. A mis-downloaded
        #      structure would otherwise be rendered and captioned as the right one.
        # PT | Confirma que o arquivo e a molecula que a figura diz ser. Uma estrutura
        #      baixada errada seria renderizada e legendada como se fosse a certa.
        if not re.search(spec["must_match"], p["title"], re.I):
            sys.exit(f"EN/PT: {pdb_id} title does not mention "
                     f"'{spec['must_match']}': {p['title']!r}")
        p["role_en"] = spec["role_en"]
        p["role_pt"] = spec["role_pt"]
        for k in ("note_en", "note_pt"):
            if k in spec:
                p[k] = spec[k]
        prov[pdb_id] = p
        res = f"{p['resolution_angstrom']} A" if p["resolution_angstrom"] else "n/a"
        print(f"  {pdb_id}  {p['method']:<20} {res:>8}   {p['title'][:56]}")

    print("\nrendering | renderizando ...")
    render_ago2(cmd, os.path.join(STRUCT_DIR, "6N4O.pdb"), f"{FIG_DIR}/ago2_guide_target.png")
    dyad = render_bace1(cmd, os.path.join(STRUCT_DIR, "4D8C.cif"), f"{FIG_DIR}/bace1_inhibitor.png")
    prov["4D8C"]["catalytic_dyad_resi_in_this_deposition"] = dyad
    prov["4D8C"]["catalytic_dyad_check"] = ("aspartates of DTGS and DSGT motifs == aspartates "
                                            "within 5 A of the inhibitor")
    prov["6CU7"]["chain_split"] = render_fibril(
        cmd, os.path.join(STRUCT_DIR, "6CU7.cif"),
        f"{FIG_DIR}/alpha_synuclein_fibril.png", "alpha-synuclein")
    prov["5OQV"]["chain_split"] = render_fibril(
        cmd, os.path.join(STRUCT_DIR, "5OQV.cif"),
        f"{FIG_DIR}/abeta42_fibril.png", "Abeta42")

    json.dump(prov, open(f"{TAB_DIR}/structure_figure_provenance.json", "w"),
              ensure_ascii=False, indent=1)

    print("\nEN | Figures are mechanistic illustration, not a result.")
    print("PT | As figuras sao ilustracao mecanistica, nao resultado.")
    for f in sorted(os.listdir(FIG_DIR)):
        print(f"  {FIG_DIR}/{f}")
    print(f"  {TAB_DIR}/structure_figure_provenance.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
