#!/usr/bin/env python3
"""Résolveur HGVS -> SPDI, basé sur les coordonnées génomiques officielles de l'OMS.

Le catalogue OMS 2e éd. fournit `WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt` :
    variant  chromosome  position  reference_nucleotide  alternative_nucleotide
qui mappe chaque variant HGVS (ex. `rrs_n.1401A>G`, `katG_p.Ser315Thr`) à sa/ses
coordonnée(s) génomique(s) sur H37Rv (NC_000962.3), en convention **1-based** (VCF).

Le pangénome du projet utilise une convention SPDI **0-based** :
    `NC_000962.3:position:ref:alt`  avec position = position_VCF - 1
(vérifié : rrs_n.1401A>G -> VCF 1473246 -> SPDI projet NC_000962.3:1473245:A:G, présent
dans pan_spdi.pkl). On applique donc systématiquement le décalage -1.

La notation `variant` de l'OMS (`gene_mutation`) est identique à `f"{Gene}_{Mutation}"`
de tb-profiler : ce résolveur sert donc aux deux sources.

Limites : seuls les variants ponctuels (ref et alt d'une seule base) produisent un SPDI
croisable avec le pangénome SNP. Les indels/MNV reçoivent `None` (représentation SPDI
non garantie identique à celle du pangénome) et sont comptés à part.
"""
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

_COORD_FILE = paths.SRC_WHO / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"


@lru_cache(maxsize=1)
def _index():
    """variant -> SPDI 0-based (str), en privilégiant la représentation mono-base.

    Renvoie aussi l'ensemble des variants connus mais non ponctuels (indel/MNV).
    """
    from collections import defaultdict
    df = pd.read_csv(_COORD_FILE, sep="\t", dtype=str)
    snv, allmap = {}, defaultdict(set)
    for v, pos, ref, alt in zip(df["variant"], df["position"],
                                df["reference_nucleotide"],
                                df["alternative_nucleotide"]):
        if not isinstance(v, str) or not isinstance(ref, str) or not isinstance(alt, str):
            continue
        # convention projet : position 0-based = position VCF (1-based) - 1
        spdi = f"NC_000962.3:{int(pos) - 1}:{ref}:{alt}"
        allmap[v].add(spdi)
        if len(ref) == 1 and len(alt) == 1:
            snv.setdefault(v, spdi)
    return snv, allmap


def resolve(variant: str):
    """SPDI 0-based d'un variant HGVS (mono-base préféré, sinon plus court indel/MNV)."""
    snv, allmap = _index()
    if variant in snv:
        return snv[variant]
    cands = allmap.get(variant)
    return min(cands, key=len) if cands else None


def resolve_all(variant: str, max_len: int = 24):
    """Toutes les représentations SPDI courtes d'un variant (SNV + indels/MNV).

    Filtre ref/alt > max_len bases (grosses délétions LoF non croisables avec le
    pangénome). Capture les frameshifts/indels courts (ex. pncA `T:TC`) — correction
    du 4e angle mort de représentation (indels)."""
    _, allmap = _index()
    out = []
    for s in allmap.get(variant, ()):
        _, _, ref, alt = s.split(":")
        if len(ref) <= max_len and len(alt) <= max_len:
            out.append(s)
    return sorted(out)


def coverage_report():
    snv, allmap = _index()
    indels = sum(1 for v in allmap if v not in snv and resolve_all(v))
    print(f"variants ponctuels (SNV) résolus : {len(snv)}")
    print(f"variants non-SNV avec ≥1 indel/MNV court résolu : {indels}")
    for v in ("rrs_n.1401A>G", "katG_p.Ser315Thr", "rpoB_p.Ser450Leu", "pncA_p.Glu173fs"):
        print(f"  {v} -> resolve={resolve(v)} | all={resolve_all(v)[:3]}")


if __name__ == "__main__":
    coverage_report()
