#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EN | Sequence-level pairing of the microRNAs of the two axes with the 3'UTR of the mRNA they
     are said to repress: the miR-29 family (miR-29a-3p, -29b-3p, -29c-3p) with BACE1, and
     miR-7-5p and miR-153-3p with SNCA. Nothing here is typed in from memory.

     Inputs (archived in data/raw/sequences_2026/, with their provenance):
       - the mature sequences of the five microRNAs, from miRBase (mature.fa, by accession
         MIMAT...);
       - the two RefSeq records NM_000345.4 (SNCA) and NM_012104.6 (BACE1), both MANE Select,
         as GenBank files from NCBI (E-utilities). The 3'UTR of a record is its sequence after
         the end of the CDS feature;
       - as a cross-check, the same two transcripts from the UCSC refGene table of hg38 and the
         hg38 reference (read by byte range): the SNCA 3'UTR is identical, and the BACE1 3'UTR
         of refGene is the RefSeq one plus 12 nt at its 3' end, which lie beyond every site.
     data/extracted/utr3_sequences.fa holds the two 3'UTRs; provenance.json holds the SHA-256 of
     every archived file and of every sequence, and the cross-check.

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
     python scripts/49_mirna_site_pairing.py --fetch    # downloads them again (needs network)

PT | Pareamento, em nivel de sequencia, dos microRNAs dos dois eixos com o 3'UTR do mRNA que
     eles reprimiriam: a familia miR-29 (miR-29a-3p, -29b-3p, -29c-3p) com o BACE1, e
     miR-7-5p e miR-153-3p com o SNCA. Nada aqui e digitado de memoria.

     Entradas (arquivadas em data/raw/sequences_2026/, com a proveniencia): as sequencias
     maduras dos cinco microRNAs, do miRBase (mature.fa, pelo acesso MIMAT...); os dois
     registros RefSeq NM_000345.4 (SNCA) e NM_012104.6 (BACE1), ambos MANE Select, como
     arquivos GenBank do NCBI (E-utilities); o 3'UTR de um registro e a sequencia depois do
     fim da feature CDS; e, como conferencia, os mesmos transcritos da tabela refGene do UCSC
     para hg38 e a referencia hg38 (lida por faixa de bytes): o 3'UTR do SNCA e identico, e o
     do BACE1 no refGene e o do RefSeq mais 12 nt na ponta 3', alem de todos os sitios.
     data/extracted/utr3_sequences.fa tem os dois 3'UTR; provenance.json tem o SHA-256 de cada
     arquivo arquivado e de cada sequencia, e a conferencia.

     Metodo. Um sitio e um trecho do 3'UTR complementar a semente do microRNA (nucleotideos
     2 a 8): 8mer (semente mais um A oposto ao nucleotideo 1), 7mer-m8 (semente, sem A1) ou
     7mer-A1 (nucleotideos 2 a 7 mais A1). Cada sitio e entao posto contra o microRNA inteiro
     sem lacunas, o nucleotideo k do microRNA diante da posicao a1 - (k - 1) do 3'UTR, e cada
     par e Watson-Crick (A-U, G-C) ou oscilante (G-U). E uma contagem deterministica de
     complementaridade, nao uma predicao de energia de ligacao nem prova de que o sitio e
     funcional; se um sitio e funcional foi testado nos artigos citados, mutando-o num
     reporter. As posicoes comecam em 1 no primeiro nucleotideo depois do codon de parada.
"""

import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import sys
import urllib.request

RAW = "data/raw/sequences_2026"
EXT = "data/extracted"
TAB = "results/tables"
PROV = f"{RAW}/provenance.json"
UTR_FA = f"{EXT}/utr3_sequences.fa"
MAT_FA = f"{RAW}/mirbase_mature_selected.fa"
REF_TSV = f"{RAW}/refGene_selected.tsv"

URL_MATURE = "https://mirbase.org/download/mature.fa"
URL_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text"
URL_REFGENE = "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz"
URL_HG38 = "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa"

# label used in the tables and figures -> miRBase accession (the name in the current miRBase file can differ, e.g. hsa-miR-153-1-3p)
MIRNAS = {"hsa-miR-29a-3p": "MIMAT0000086", "hsa-miR-29b-3p": "MIMAT0000100", "hsa-miR-29c-3p": "MIMAT0000681",
          "hsa-miR-7-5p": "MIMAT0000252", "hsa-miR-153-3p": "MIMAT0000439"}
GENES = {"BACE1": "NM_012104.6", "SNCA": "NM_000345.4"}
REFGENE_ID = {"BACE1": "NM_012104", "SNCA": "NM_000345"}
PAIRS = {("A", "U"): "WC", ("U", "A"): "WC", ("G", "C"): "WC", ("C", "G"): "WC", ("G", "U"): "GU", ("U", "G"): "GU"}
COMP = str.maketrans("ACGU", "UGCA")


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def sha_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def rc_rna(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """EN/PT: records keyed by the first word of the header | registros indexados pela primeira palavra."""
    out, name = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = {"header": line[1:], "seq": ""}
        elif name:
            out[name]["seq"] += line.strip()
    return out


def parse_genbank(path):
    """EN/PT: version, MANE keyword, CDS span and sequence of a GenBank record | versao, MANE, CDS e sequencia."""
    txt = open(path, encoding="utf-8").read()
    version = re.search(r"^VERSION\s+(\S+)", txt, re.M).group(1)
    mane = "MANE Select" in re.search(r"^KEYWORDS\s+(.*)$", txt, re.M).group(1)
    cds = re.search(r"^     CDS\s+(\d+)\.\.(\d+)\s*$", txt, re.M)
    origin = txt[txt.index("ORIGIN"):]
    seq = "".join(re.findall(r"[acgtn]+", "".join(l[10:] for l in origin.splitlines()[1:] if l.strip() and not l.startswith("//")))).upper().replace("T", "U")
    return {"version": version, "mane_select": mane, "cds": (int(cds.group(1)), int(cds.group(2))), "seq": seq}


def utr3_of(gb):
    return gb["seq"][gb["cds"][1]:]


def http(url, rng=None, head=False):
    req = urllib.request.Request(url, headers={"User-Agent": "mirna-mimics-kinetic-model/1.0"})
    if rng:
        req.add_header("Range", f"bytes={rng[0]}-{rng[1]}")
    if head:
        req.get_method = lambda: "HEAD"
    with urllib.request.urlopen(req, timeout=180) as r:
        return (r.headers, r.read())


def hg38_utr3(rows, fai, gene):
    """EN/PT: the 3'UTR of a refGene row read from hg38, for the cross-check | 3'UTR de uma linha refGene lido do hg38."""
    f = rows[REFGENE_ID[gene]].split("\t")
    chrom, strand, cs, ce = f[2], f[3], int(f[6]), int(f[7])
    es = [int(x) for x in f[9].strip(",").split(",")]
    ee = [int(x) for x in f[10].strip(",").split(",")]
    parts = sorted((max(s_, ce), e_) if strand == "+" else (s_, min(e_, cs)) for s_, e_ in zip(es, ee))
    parts = [(a, b) for a, b in parts if b > a]
    dna = ""
    for a, b in parts:
        _, off, lb, lw = fai[chrom]
        b0 = off + a // lb * lw + a % lb
        b1 = off + (b - 1) // lb * lw + (b - 1) % lb
        dna += http(URL_HG38, rng=(b0, b1))[1].decode().replace("\n", "").upper()
    if strand == "-":
        dna = dna[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    return dna.replace("T", "U")


def fetch_all():
    """EN/PT: download the inputs again and rewrite the archived files | baixa de novo as entradas e reescreve os arquivos."""
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(EXT, exist_ok=True)
    prov = {"retrieved_utc": datetime.datetime.utcnow().strftime("%Y-%m-%d"), "sources": {}, "sha256": {}}
    # --- miRBase mature sequences, by accession
    h, body = http(URL_MATURE)
    recs, name, acc = {}, None, None
    for line in body.decode().splitlines():
        if line.startswith(">"):
            parts = line[1:].split()
            name, acc = parts[0], parts[1]
            recs[acc] = [line, ""]
        elif acc:
            recs[acc][1] += line.strip()
    prov["sources"]["mature"] = {"url": URL_MATURE, "bytes": len(body), "sha256_whole_file": hashlib.sha256(body).hexdigest(),
                                 "last_modified": h.get("Last-Modified"), "note": "mature.fa from mirbase.org; the file does not state the release"}
    with open(MAT_FA, "w", encoding="utf-8") as f:
        for lab, a in MIRNAS.items():
            f.write(recs[a][0] + "\n" + recs[a][1] + "\n")
    # --- RefSeq records from NCBI
    prov["sources"]["refseq"] = {"url": URL_EFETCH, "records": {}}
    seqs = {}
    for g, acc in GENES.items():
        _, body = http(URL_EFETCH.format(acc=acc))
        path = f"{RAW}/{acc}.gb"
        open(path, "wb").write(body)
        gb = parse_genbank(path)
        assert gb["version"] == acc and gb["mane_select"], (gb["version"], gb["mane_select"])
        seqs[g] = utr3_of(gb)
        prov["sources"]["refseq"]["records"][acc] = {"gene": g, "mane_select": gb["mane_select"], "cds": list(gb["cds"]),
                                                       "transcript_length": len(gb["seq"]), "utr3_length": len(seqs[g])}
        prov["sha256"][f"{acc}.gb"] = sha_file(path)
    with open(UTR_FA, "w", encoding="utf-8") as f:
        for g, rna in seqs.items():
            acc = GENES[g]
            cds = prov["sources"]["refseq"]["records"][acc]["cds"]
            f.write(f">{g}|{acc}|3UTR|{cds[1] + 1}-{cds[1] + len(rna)} of the record|NCBI RefSeq (MANE Select)\n")
            for i in range(0, len(rna), 70):
                f.write(rna[i:i + 70] + "\n")
            prov["sha256"][g] = sha(rna)
    prov["sha256"]["mirbase_mature_selected.fa"] = sha_file(MAT_FA)
    # --- cross-check against UCSC refGene and the hg38 reference
    try:
        hr, body = http(URL_REFGENE)
        rows = {}
        for line in gzip.decompress(body).decode().splitlines():
            rows[line.split("\t")[1]] = line
        with open(REF_TSV, "w", encoding="utf-8") as f:
            for g in GENES:
                f.write(rows[REFGENE_ID[g]] + "\n")
        hh, _ = http(URL_HG38, head=True)
        fai = {}
        for line in http(URL_HG38 + ".fai")[1].decode().splitlines():
            a = line.split("\t")
            fai[a[0]] = tuple(int(x) for x in a[1:5])
        cc = {}
        for g in GENES:
            u = hg38_utr3(rows, fai, g)
            r = seqs[g]
            cc[g] = {"hg38_refgene_utr3_length": len(u), "refseq_utr3_length": len(r), "identical": u == r,
                     "refseq_is_prefix_of_hg38": u.startswith(r), "extra_nt_at_3prime_end_in_hg38": len(u) - len(r) if u.startswith(r) else None}
        prov["crosscheck_hg38"] = {"refGene": {"url": URL_REFGENE, "last_modified": hr.get("Last-Modified"), "bytes": len(body)},
                                   "hg38": {"url": URL_HG38, "last_modified": hh.get("Last-Modified"), "bytes": int(hh.get("Content-Length", 0))},
                                   "result": cc}
        prov["sha256"]["refGene_selected.tsv"] = sha_file(REF_TSV)
    except Exception as exc:       # the cross-check is not needed to run the scan
        print("EN/PT: cross-check against hg38 skipped:", exc)
    json.dump(prov, open(PROV, "w", encoding="utf-8"), indent=1)
    print("EN/PT: inputs fetched and archived in", RAW, "and", EXT)


def load_inputs():
    for f in (PROV, UTR_FA, MAT_FA):
        if not os.path.exists(f):
            sys.exit(f"EN/PT: missing {f}; run once with --fetch")
    prov = json.load(open(PROV, encoding="utf-8"))
    if sha_file(MAT_FA) != prov["sha256"]["mirbase_mature_selected.fa"]:
        sys.exit(f"EN/PT: {MAT_FA} no longer matches the SHA-256 in provenance.json")
    fa = {k.split("|")[0]: v for k, v in read_fasta(UTR_FA).items()}
    seqs = {}
    for g, acc in GENES.items():
        gbp = f"{RAW}/{acc}.gb"
        if sha_file(gbp) != prov["sha256"][f"{acc}.gb"]:
            sys.exit(f"EN/PT: {gbp} no longer matches the SHA-256 in provenance.json")
        gb = parse_genbank(gbp)
        u = utr3_of(gb)
        if u != fa[g]["seq"] or sha(u) != prov["sha256"][g]:
            sys.exit(f"EN/PT: the 3'UTR of {g} in {UTR_FA} is not the sequence after the CDS of {acc}")
        seqs[g] = {"seq": u, "header": fa[g]["header"], "accession": acc, "cds": gb["cds"], "transcript_length": len(gb["seq"])}
    by_acc = {r["header"].split()[1]: r for r in read_fasta(MAT_FA).values()}
    mir = {}
    for lab, a in MIRNAS.items():
        if a not in by_acc:
            sys.exit(f"EN/PT: {a} ({lab}) is missing from {MAT_FA}")
        mir[lab] = {"seq": by_acc[a]["seq"], "accession": a, "mirbase_name": by_acc[a]["header"].split()[0]}
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
    for gene, nm in GENES.items():
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
                    "gene": gene, "transcript": nm, "utr3_length_nt": len(utr), "mirna": m, "mirna_accession": mir[m]["accession"], "mirbase_name": mir[m]["mirbase_name"],
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
        "utr3": {g: {"transcript": GENES[g], "length_nt": len(seqs[g]["seq"]), "cds": list(seqs[g]["cds"]), "transcript_length_nt": seqs[g]["transcript_length"], "header": seqs[g]["header"]} for g in seqs},
        "mirna": {m: {"accession": mir[m]["accession"], "mirbase_name": mir[m]["mirbase_name"], "sequence": mir[m]["seq"], "seed_2_8": mir[m]["seq"][1:8],
                      "site_motif_8mer": rc_rna(mir[m]["seq"][1:8]) + "A"} for m in MIRNAS},
        "n_sites": len(rows),
        "sites_by_pair": {f"{g}|{m}": sum(1 for r in rows if r["gene"] == g and r["mirna"] == m) for g in GENES for m in MIRNAS},
        "pair_convention": "nucleotide k of the microRNA faces 3'UTR position a1 - (k - 1); Watson-Crick and G-U pairs; no gaps",
        "position_convention": "1-based from the first nucleotide after the CDS of the RefSeq record",
    }
    json.dump(summary, open(f"{TAB}/mirna_site_pairing.json", "w", encoding="utf-8"), indent=1)
    print("=" * 78)
    print("EN | Sequence-level pairing of the microRNAs with the 3'UTR | PT | Pareamento em nivel de sequencia")
    print("=" * 78)
    for g in GENES:
        print(f"{g} {GENES[g]} 3'UTR: {len(seqs[g]['seq'])} nt")
    for r in rows:
        print(f"  {r['gene']:5s} {r['mirna']:15s} {r['site_type']:8s} nt {r['utr_start']}-{r['utr_end']}  pairs {r['pairs_total']}/{r['mirna_length_nt']} (seed {r['seed_pairs_2_8']}/7, 13-16: {r['supplementary_pairs_13_16']}/4)")
    print(f"EN/PT -> {TAB}/mirna_site_pairing.csv and .json")
    return 0


if __name__ == "__main__":
    if "--fetch" in sys.argv:
        fetch_all()
    sys.exit(run())
