#!/usr/bin/env python3
"""P7.16 — Lignes d'assertion contre chaînes de coordonnées DISTINCTES, par type, dans un catalogue.

Objet      : le manuscrit A1 écrivait « the ingested catalogue holds 52 041 single-base coordinate
             strings and 772 blocks » alors que 52 041 et 772 sont des comptes de LIGNES : un
             catalogue qui ne porte que 37 781 chaînes distinctes ne peut pas en porter 52 041
             (review du 2026-11-21, point 5). Le nombre était exact contre sa source, c'est le mot
             qui le nommait qui était faux, et aucun test d'épinglage ne pouvait voir une erreur
             d'unité. Ce script rend les deux unités côte à côte, et son contrôle est que la somme
             des chaînes distinctes par type redonne le compte de chaînes distinctes du catalogue.
Entrées    : --catalogue (TSV au schéma de `catalogue_consolide.tsv` ; défaut : le catalogue vivant).
Sorties    : --sortie (TSV : type, lignes, chaines_distinctes) ; défaut
             résultats/model/p716_chaines_distinctes.tsv.
Réutilisable : oui (tout catalogue au même schéma, C0 comme C1).
Projet     : Resistance_antibio (P7.16)
Date       : 2026-11-21
Usage      : python3 analyses/p716_chaines_distinctes.py
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

DEFAUT = Path(paths.ROOT) / "résultats" / "model" / "p716_chaines_distinctes.tsv"


def type_de(spdi):
    if not spdi:
        return "sans_spdi"
    p = spdi.split(":")
    ref, alt = p[2], p[3]
    if len(ref) == len(alt) == 1:
        return "snp"
    if len(ref) == len(alt):
        return "mnv"
    return "indel"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", default=str(paths.CATALOGUE_TSV))
    ap.add_argument("--sortie", default=str(DEFAUT))
    a = ap.parse_args()

    lignes, distincts = Counter(), defaultdict(set)
    with open(a.catalogue, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s = (r.get("spdi") or "").strip()
            t = type_de(s)
            lignes[t] += 1
            if s:
                distincts[t].add(s)

    tous = set().union(*distincts.values()) if distincts else set()
    out = Path(a.sortie)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["type", "lignes", "chaines_distinctes"])
        for t in ("snp", "mnv", "indel", "sans_spdi"):
            w.writerow([t, lignes[t], len(distincts[t])])
        w.writerow(["TOTAL", sum(lignes.values()), len(tous)])

    somme = sum(len(distincts[t]) for t in ("snp", "mnv", "indel"))
    for t in ("snp", "mnv", "indel", "sans_spdi"):
        print(f"{t:10s} lignes={lignes[t]:7d}  chaines_distinctes={len(distincts[t]):7d}")
    print(f"{'TOTAL':10s} lignes={sum(lignes.values()):7d}  chaines_distinctes={len(tous):7d}")
    assert somme == len(tous), f"somme par type {somme} != total {len(tous)}"
    print(f"\ncontrôle d'unité OK : {somme} = somme des chaînes distinctes par type\n-> {out}")


if __name__ == "__main__":
    main()
