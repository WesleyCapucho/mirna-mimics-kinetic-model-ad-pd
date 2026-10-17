#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Sequence-level pairing of the microRNAs of the two axes with the 3'UTR of the mRNA they
     are said to repress: the miR-29 family (miR-29a-3p, -29b-3p, -29c-3p) with BACE1, and
     miR-7-5p and miR-153-3p with SNCA. Nothing here is typed in from memory.

     Inputs (archived in data/raw/sequences_2026/, with their provenance):
       - the mature sequences of the five microRNAs, as written in a copy of miRBase's
         mature.fa (accessions MIMAT...);
       - the two RefSeq transcripts NM_012104 (BACE1) and NM_000345 (SNCA), as rows of the
         UCSC refGene table of assembly hg38;
       - the genomic sequence of their 3'UTRs, read from the hg38 reference by byte range and
         assembled exon by exon (data/extracted/utr3_sequences.fa; each sequence carries its
         SHA-256 in provenance.json).
     The 3'UTR of a transcript is its exonic sequence after the stop codon (refGene cdsEnd on
     the plus strand, cdsStart on the minus strand; both genes lie on the minus strand).

     Method. A site is a stretch of the 3'UTR that is complementary to the seed of the
     microRNA (nucleotides 2 to 8): 8mer (seed match plus an A opposite nucleotide 1),
     7mer-m8 (seed match, no A1) or 7mer-A1 (nucleotides 2 to 7 plus an A1). Every site is
     then laid against the whole microRNA without gaps, nucleotide k of the microRNA facing
     3'UTR position a1 - (k - 1), and each pair is Watson-Crick (A-U, G-C) or wobble (G-U).
     This is a deterministic count of complementarity, not a prediction of binding energy
     and not evidence that the site is functional; whether a site is functional was tested
     in the cited papers by mutating it in a reporter. Positions are 1-based from the first
     nucleotide after the stop codon.

     python scripts/49_mirna_site_pairing.py            # reads the archived inputs, checks their SHA-256
     python scripts/49_mirna_site_pairing.py --fetch    # downloads them again from the public mirrors (needs network)

PT | Pareamento, em nivel de sequencia, dos microRNAs dos dois eixos com o 3'UTR do mRNA que
     eles reprimiriam: a familia miR-29 (miR-29a-3p, -29b-3p, -29c-3p) com o BACE1, e
     miR-7-5p e miR-153-3p com o SNCA. Nada aqui e digitado de memoria.

     Entradas (arquivadas em data/raw/sequences_2026/, com a proveniencia): as sequencias
     maduras dos cinco microRNAs, como escritas numa copia do mature.fa do miRBase
     (acessos MIMAT...); os dois transcritos RefSeq NM_012104 (BACE1) e NM_000345 (SNCA),
     como linhas da tabela refGene do UCSC na montagem hg38; e a sequencia genomica dos
     3'UTR, lida da referencia hg38 por faixa de bytes e montada exon a exon
     (data/extracted/utr3_sequences.fa; cada sequencia tem o SHA-256 em provenance.json).
     O 3'UTR e a sequencia exonica depois do codon de parada (cdsEnd na fita positiva,
     cdsStart na negativa; os dois genes estao na fita negativa).

     Metodo. Um sitio e um trecho do 3'UTR complementar a semente do microRNA (nucleotideos
     2 a 8): 8mer (semente mais um A oposto ao nucleotideo 1), 7mer-m8 (semente, sem A1) ou
     7mer-A1 (nucleotideos 2 a 7 mais A1). Cada sitio e entao posto contra o microRNA inteiro
     sem lacunas, o nucleotideo k do microRNA diante da posicao a1 - (k - 1) do 3'UTR, e cada
     par e Watson-Crick (A-U, G-C) ou oscilante (G-U). E uma contagem deterministica de
     complementaridade, nao uma predicao de energia de ligacao nem prova de que o sitio e
     funcional; se um sitio e funcional foi testado nos artigos citados, mutando-o num
     repórter. As posicoes comecam em 1 no primeiro nucleotideo depois do codon de parada.
"""

import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import urllib.request

RAW = "data/raw/sequences_2026"
EXT = "data/extracted"
TAB = "results/tables"
PROV = f"{RAW}/provenance.json"
UTR_FA = f"{EXT}/utr3_sequences.fa"
MAT_FA = f"{RAW}/mirbase_mature_selected.fa"
REF_TSV = f"{RAW}/refGene_selected.tsv"

URL_MATURE = "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/mature.fa"
URL_REFGENE = "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz"
URL_HG38 = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa"

MIRNAS = ["hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p", "hsa-miR-7-5p", "hsa-miR-153-3p"]
TRANSCRIPTS = {"BACE1": "NM_012104", "SNCA": "NM_000345"}
PAIRS = {("A", "U"): "WC", ("U", "A"): "WC", ("G", "C"): "WC", ("C", "G"): "WC", ("G", "U"): "GU", ("U", "G"): "GU"}
COMP = str.maketrans("ACGU", "UGCA")


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def rc_rna(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, name = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
            out.setdefault("_header_" + name, line[1:])
        elif name:
            out[name] += line.strip()
    return out


def http(url, rng=None, head=False):
    req = urllib.request.Request(url, headers={"User-Agent": "mirna-mimics-kinetic-model/1.0"})
    if rng:
        req.add_header("Range", f"bytes={rng[0]}-{rng[1]}")
    if head:
        req.get_method = lambda: "HEAD"
    with urllib.request.urlopen(req, timeout=120) as r:
        return (r.headers, r.read())


def fetch_all():
    """EN/PT: download the inputs again from the public mirrors and rewrite the archived files | baixa de novo as entradas."""
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(EXT, exist_ok=True)
    prov = {"retrieved_utc": __import__("datetime").datetime.utcnow().strftime("%Y-%m-%d"), "sources": {}, "sha256": {}}
    # --- miRBase mature sequences (a copy of mature.fa kept by nf-core/test-datasets)
    h, body = http(URL_MATURE)
    text = body.decode()
    prov["sources"]["mature"] = {"url": URL_MATURE, "etag": h.get("ETag"), "bytes": len(body), "sha256_whole_file": hashlib.sha256(body).hexdigest(),
                                 "note": "a copy of the miRBase mature.fa kept by nf-core/test-datasets; the release of miRBase is not stated in the file"}
    recs = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = [line, ""]
        elif name:
            recs[name][1] += line.strip()
    with open(MAT_FA, "w", encoding="utf-8") as f:
        for m in MIRNAMES():
            f.write(recs[m][0] + "\n" + recs[m][1] + "\n")
    # --- refGene rows
    h, body = http(URL_REFGENE)
    prov["sources"]["refGene"] = {"url": URL_REFGENE, "last_modified": h.get("Last-Modified"), "etag": h.get("ETag"), "bytes": len(body),
                                  "sha256_whole_file": hashlib.sha256(body).hexdigest(), "table": "UCSC refGene, assembly hg38"}
    rows = {}
    for line in gzip.decompress(body).decode().splitlines():
        f = line.split("\t")
        rows[f[1]] = line
    with open(REF_TSV, "w", encoding="utf-8") as f:
        for g, nm in TRANSCRIPTS.items():
            f.write(rows[nm] + "\n")
    # --- genome by byte range
    h, _ = http(URL_HG38, head=True)
    fai = {}
    _, fb = http(URL_HG38 + ".fai")
    for line in fb.decode().splitlines():
        a = line.split("\t")
        fai[a[0]] = tuple(int(x) for x in a[1:5])
    prov["sources"]["hg38"] = {"url": URL_HG38, "last_modified": h.get("Last-Modified"), "etag": h.get("ETag"), "bytes": int(h.get("Content-Length", 0)),
                               "note": "UCSC hg38 (GRCh38), read by HTTP byte range using the .fai index"}

    def genome(chrom, start, end):
        _, off, lb, lw = fai[chrom]
        b0 = off + start // lb * lw + start % lb
        b1 = off + (end - 1) // lb * lw + (end - 1) % lb
        s = http(URL_HG38, rng=(b0, b1))[1].decode().replace("\n", "").upper()
        assert len(s) == end - start
        return s

    seqs = {}
    for g, nm in TRANSCRIPTS.items():
        f = rows[nm].split("\t")
        chrom, strand, cs, ce = f[2], f[3], int(f[6]), int(f[7])
        es = [int(x) for x in f[9].strip(",").split(",")]
        ee = [int(x) for x in f[10].strip(",").split(",")]
        parts = []
        for s_, e_ in zip(es, ee):
            a, b = (max(s_, ce), e_) if strand == "+" else (s_, min(e_, cs))
            if b > a:
                parts.append((a, b))
        parts.sort()
        dna = "".join(genome(chrom, a, b) for a, b in parts)
        if strand == "-":
            dna = dna[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
        seqs[g] = (dna.replace("T", "U"), nm, chrom, strand, parts)
    with open(UTR_FA, "w", encoding="utf-8") as f:
        for g, (rna, nm, chrom, strand, parts) in seqs.items():
            span = ";".join(f"{chrom}:{a + 1}-{b}" for a, b in parts)
            f.write(f">{g}|{nm}|3UTR|{chrom}({strand})|{span}|hg38\n")
            for i in range(0, len(rna), 70):
                f.write(rna[i:i + 70] + "\n")
            prov["sha256"][g] = sha(rna)
    prov["sha256"]["mirbase_mature_selected.fa"] = hashlib.sha256(open(MAT_FA, "rb").read()).hexdigest()
    prov["sha256"]["refGene_selected.tsv"] = hashlib.sha256(open(REF_TSV, "rb").read()).hexdigest()
    json.dump(prov, open(PROV, "w", encoding="utf-8"), indent=1)
    print("EN/PT: inputs fetched and archived in", RAW, "and", EXT)


def MIRNAMES():
    return MIRNAS


def load_inputs():
    for f in (PROV, UTR_FA, MAT_FA, REF_TSV):
        if not os.path.exists(f):
            sys.exit(f"EN/PT: missing {f}; run once with --fetch")
    prov = json.load(open(PROV, encoding="utf-8"))
    mature = read_fasta(MAT_FA)
    utr = read_fasta(UTR_FA)
    # the archived files must still be the files whose hashes were recorded
    for f, key in ((MAT_FA, "mirbase_mature_selected.fa"), (REF_TSV, "refGene_selected.tsv")):
        if hashlib.sha256(open(f, "rb").read()).hexdigest() != prov["sha256"][key]:
            sys.exit(f"EN/PT: {f} no longer matches the SHA-256 in provenance.json")
    seqs = {}
    for key in list(utr):
        if key.startswith("_header_"):
            continue
        gene = key.split("|")[0]
        if sha(utr[key]) != prov["sha256"][gene]:
            sys.exit(f"EN/PT: the 3'UTR of {gene} no longer matches the SHA-256 in provenance.json")
        seqs[gene] = {"seq": utr[key], "header": utr["_header_" + key]}
    mir = {}
    for key in list(mature):
        if key.startswith("_header_"):
            continue
        mir[key] = {"seq": mature[key], "accession": mature["_header_" + key].split()[1]}
    for m in MIRNAS:
        if m not in mir:
            sys.exit(f"EN/PT: {m} is missing from {MAT_FA}")
    return prov, mir, seqs


def find_sites(utr, mir_seq):
    """EN/PT: every 8mer, 7mer-m8 and 7mer-A1 site, with the position of the nucleotide that faces nucleotide 1 | todos os sitios."""
    seed78 = rc_rna(mir_seq[1:8])        # nucleotides 2 to 8
    seed27 = rc_rna(mir_seq[1:7])        # nucleotides 2 to 7
    out = []
    for i in range(len(utr) - 6):
        if utr[i:i + 7] == seed78:
            a1 = i + 7                    # 0-based position facing nucleotide 1
            kind = "8mer" if a1 < len(utr) and utr[a1] == "A" else "7mer-m8"
            out.append((kind, i, i + 7 if kind == "8mer" else i + 6, a1))
    for i in range(1, len(utr) - 6):
        if utr[i:i + 6] == seed27 and utr[i + 6] == "A" and utr[i - 1:i + 6] != seed78:
            out.append(("7mer-A1", i, i + 6, i + 6))
    return sorted(out, key=lambda x: x[1])


def pair_row(utr, mir_seq, a1):
    """EN/PT: nucleotide k of the microRNA faces 3'UTR position a1 - (k - 1), without gaps | sem lacunas."""
    marks, tgt = [], []
    for k in range(1, len(mir_seq) + 1):
        p = a1 - (k - 1)
        t = utr[p] if 0 <= p < len(utr) else "-"
        tgt.append(t)
        marks.append(PAIRS.get((mir_seq[k - 1], t), ""))
    return marks, tgt


def run():
    prov, mir, seqs = load_inputs()
    rows = []
    for gene, nm in TRANSCRIPTS.items():
        utr = seqs[gene]["seq"]
        for m in MIRNAS:
            ms = mir[m]["seq"]
            for kind, i0, i1, a1 in find_sites(utr, ms):
                marks, tgt = pair_row(utr, ms, a1)
                seed = sum(1 for k in range(2, 9) if marks[k - 1])
                supp = sum(1 for k in range(13, 17) if k <= len(ms) and marks[k - 1])
                tot = sum(1 for x in marks if x)
                wob = sum(1 for x in marks if x == "GU")
                rows.append({
                    "gene": gene, "transcript": nm, "utr3_length_nt": len(utr), "mirna": m, "mirna_accession": mir[m]["accession"],
                    "mirna_length_nt": len(ms), "site_type": kind, "utr_start": i0 + 1, "utr_end": i1 + 1, "a1_position": a1 + 1,
                    "seed_pairs_2_8": seed, "supplementary_pairs_13_16": supp, "pairs_total": tot, "wobble_pairs": wob,
                    "target_5to3": "".join(reversed(tgt)),
                    "pairing_marks": "".join({"WC": "|", "GU": ":", "": " "}[x] for x in reversed(marks)),
                })
    os.makedirs(TAB, exist_ok=True)
    with open(f"{TAB}/mirna_site_pairing.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    summary = {
        "inputs": {"provenance": prov["sources"], "sha256": prov["sha256"]},
        "utr3": {g: {"transcript": TRANSCRIPTS[g], "length_nt": len(seqs[g]["seq"]), "header": seqs[g]["header"]} for g in seqs},
        "mirna": {m: {"accession": mir[m]["accession"], "sequence": mir[m]["seq"], "seed_2_8": mir[m]["seq"][1:8],
                      "site_motif_8mer": rc_rna(mir[m]["seq"][1:8]) + "A"} for m in MIRNAS},
        "n_sites": len(rows),
        "sites_by_pair": {f"{g}|{m}": sum(1 for r in rows if r["gene"] == g and r["mirna"] == m) for g in TRANSCRIPTS for m in MIRNAS},
        "pair_convention": "nucleotide k of the microRNA faces 3'UTR position a1 - (k - 1); Watson-Crick and G-U pairs; no gaps",
        "position_convention": "1-based from the first nucleotide after the stop codon of the RefSeq transcript",
    }
    json.dump(summary, open(f"{TAB}/mirna_site_pairing.json", "w", encoding="utf-8"), indent=1)
    print("=" * 78)
    print("EN | Sequence-level pairing of the microRNAs with the 3'UTR | PT | Pareamento em nivel de sequencia")
    print("=" * 78)
    for g in TRANSCRIPTS:
        print(f"{g} {TRANSCRIPTS[g]} 3'UTR: {len(seqs[g]['seq'])} nt")
    for r in rows:
        print(f"  {r['gene']:5s} {r['mirna']:15s} {r['site_type']:8s} nt {r['utr_start']}-{r['utr_end']}  pairs {r['pairs_total']}/{r['mirna_length_nt']} (seed {r['seed_pairs_2_8']}/7, 13-16: {r['supplementary_pairs_13_16']}/4)")
    print(f"EN/PT -> {TAB}/mirna_site_pairing.csv and .json")
    return 0


if __name__ == "__main__":
    if "--fetch" in sys.argv:
        fetch_all()
    sys.exit(run())
