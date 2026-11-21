"""Propagation des assertions de résistance au niveau ACIDE AMINÉ (codons synonymes).

Brique de conversion AA -> SPDI exacte (gère le brin du gène, CDS H37Rv). Pour un
variant missense `gene_p.{Ref}{N}{Alt}`, génère tous les SPDI mono-base produisant le
même changement d'acide aminé (codons alternatifs codant l'AA muté).

Motivation : un catalogue indexé par UNE représentation SPDI par variant HGVS_p rate
les codons synonymes-d'AA d'un déterminant connu (ex. embB Met306Ile catalogué via
ATG->ATT mais pas ATG->ATC). Cette propagation comble cet angle mort (cf. cahier :
+23 pts de sensibilité sur l'éthambutol).

Validé : katG Ser315Thr (Rv1908c, brin −) -> NC_000962.3:2155167:C:G ;
         embB Met306Ile (brin +) -> :G:T (catalogué) + :G:C + :G:A.
"""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

import paths

_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON2AA = dict(zip(_CODONS, _AAS))
AA2CODONS = defaultdict(list)
for _c, _a in CODON2AA.items():
    AA2CODONS[_a].append(_c)
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PMISSENSE = re.compile(r"_p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def parse_cds_seq():
    """gene_symbol -> (debut1based, fin1based, reverse, coding_seq)."""
    txt = (paths.FASTAS / "NC_000962.3_CDS.fasta").read_text()
    genes, cur, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if cur:
                genes[cur[0]] = (cur[1], cur[2], cur[3], "".join(buf).upper())
            buf = []
            g = re.search(r"\[gene=([^\]]+)\]", line)
            lt = re.search(r"\[locus_tag=([^\]]+)\]", line)
            loc = re.search(r"\[location=([^\]]+)\]", line)
            cur = None
            if loc:
                nums = re.findall(r"\d+", loc.group(1))
                if nums:
                    sym = g.group(1) if g else (lt.group(1) if lt else None)
                    cur = (sym, int(nums[0]), int(nums[-1]), "complement" in loc.group(1))
        else:
            buf.append(line.strip())
    if cur:
        genes[cur[0]] = (cur[1], cur[2], cur[3], "".join(buf).upper())
    return genes


def propagate(gene, n, ref1, alt1, genes):
    """SPDI 0-based mono-base produisant ref1->alt1 au codon n du gène. [] si incohérent."""
    if gene not in genes:
        return []
    deb, fin, rev, seq = genes[gene]
    codon_ref = seq[3 * (n - 1):3 * n]
    if len(codon_ref) != 3 or CODON2AA.get(codon_ref) != ref1:
        return []
    out = []
    for codon_alt in AA2CODONS.get(alt1, []):
        diff = [i for i in range(3) if codon_alt[i] != codon_ref[i]]
        if len(diff) != 1:
            continue
        i = diff[0]
        rn, an = codon_ref[i], codon_alt[i]
        if not rev:
            g1, gref, galt = deb + 3 * (n - 1) + i, rn, an
        else:
            g1, gref, galt = fin - (3 * (n - 1) + i), COMP[rn], COMP[an]
        out.append(f"NC_000962.3:{g1 - 1}:{gref}:{galt}")
    return out


def selftest(genes=None):
    g = genes or parse_cds_seq()
    a = "NC_000962.3:2155167:C:G" in propagate("katG", 315, "S", "T", g)
    e = propagate("embB", 306, "M", "I", g)
    b = "NC_000962.3:4247430:G:T" in e and "NC_000962.3:4247430:G:C" in e
    return a and b


def propagate_catalogue(cat: pd.DataFrame) -> pd.DataFrame:
    """Ajoute au catalogue les SPDI codons-synonymes des variants missense R-associated.

    Les lignes ajoutées héritent du verdict, avec source `<source>_AAprop`.
    """
    genes = parse_cds_seq()
    if not selftest(genes):
        raise RuntimeError("aa_propagation selftest échoué (parsing CDS / brin)")
    existing = set(cat["spdi"])
    rows = []
    for _, r in cat[cat["call"] == "R-associated"].iterrows():
        m = PMISSENSE.search(str(r["variant"]))
        if not m:
            continue
        ref1, n, alt1 = AA3TO1.get(m.group(1)), int(m.group(2)), AA3TO1.get(m.group(3))
        if not ref1 or not alt1 or ref1 == alt1:
            continue
        for spdi in propagate(r["gene"], n, ref1, alt1, genes):
            if spdi in existing or spdi == "":
                continue
            existing.add(spdi)
            nr = dict(r)
            nr.update({"spdi": spdi, "source": f"{r['source']}_AAprop",
                       "notes": f"AA-propagated from {r['variant']}"})
            rows.append(nr)
    if not rows:
        return cat
    return pd.concat([cat, pd.DataFrame(rows)], ignore_index=True)
