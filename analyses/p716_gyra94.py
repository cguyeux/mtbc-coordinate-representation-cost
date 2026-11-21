#!/usr/bin/env python3
"""P7.16 — Énumération complète des entrées du codon 94 de gyrA dans le fichier de coordonnées de l'OMS.

Objet      : le manuscrit A1 écrivait, en quatre endroits dont le résumé, que le fichier officiel de
             coordonnées porte CINQ variants protéiques au codon 94 de `gyrA` et qu'il donne une
             entrée combinée avec le voisin Ser95Thr pour DEUX d'entre eux. La review du 2026-11-21
             (point 1) a rouvert le fichier : il en porte HUIT, et CINQ reçoivent une entrée combinée,
             dont deux sous la forme d'un bloc de CINQ bases à la position 7581, c'est-à-dire à la
             position et dans la forme exactes des trois haplotypes que le manuscrit déclare absents.
             La conclusion du manuscrit survit (les trois sans entrée sont bien Asn, His et Tyr) mais
             ses comptes étaient faux, et la preuve devient une identité au lieu d'une analogie.
             Ce script rend le compte depuis la source, pour que le manuscrit l'épingle à un TSV
             rejouable au lieu d'une valeur recopiée d'un registre (leçon B1 de l'état : un énoncé sur
             ce qu'une source externe contient se vérifie en énumérant TOUTES ses lignes pour
             l'ensemble concerné, jamais pour le seul exemple cité).
Entrées    : data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt (positions
             1-based de la convention VCF), --gene et --codon pour rejouer sur un autre codon.
Sorties    : results/model/p716_gyra94.tsv (une ligne par variant protéique : nombre de
             représentations, longueur du plus long bloc, blocs de 4 bases et plus) et
             `<sortie>_synthese.tsv` (les comptes que le manuscrit imprime).
Réutilisable : oui (tout gène et tout codon du fichier de coordonnées de l'OMS).
Projet     : Resistance_antibio (P7.16)
Date       : 2026-11-21
Usage      : python3 analyses/p716_gyra94.py
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

SORTIE = Path(paths.ROOT) / "results" / "model" / "p716_gyra94.tsv"
SOURCE = Path(paths.SRC_WHO) / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
TROIS = ("Asn", "His", "Tyr")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="gyrA")
    ap.add_argument("--codon", type=int, default=94)
    ap.add_argument("--source", default=str(SOURCE))
    ap.add_argument("--sortie", default=str(SORTIE))
    a = ap.parse_args()

    motif = re.compile(rf"^{re.escape(a.gene)}_p\.[A-Z][a-z]{{2}}{a.codon}([A-Z][a-z]{{2}})$")
    par_variant = defaultdict(list)
    with open(a.source, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m = motif.match((r.get("variant") or "").strip())
            if not m:
                continue
            ref = (r.get("reference_nucleotide") or "").strip()
            alt = (r.get("alternative_nucleotide") or "").strip()
            par_variant[m.group(1)].append((int(r["position"]), ref, alt))

    lignes, n_combine, n_cinq = [], 0, 0
    for aa in sorted(par_variant):
        reps = par_variant[aa]
        blocs = [(p, ref, alt) for p, ref, alt in reps if len(ref) == len(alt) >= 4]
        plus_long = max((len(ref) for _, ref, _ in blocs), default=0)
        if blocs:
            n_combine += 1
        if plus_long >= 5:
            n_cinq += 1
        lignes.append({
            "variant": f"{a.gene}_p.Asp{a.codon}{aa}" if a.gene == "gyrA" else f"{a.gene}_p.X{a.codon}{aa}",
            "acide_amine": aa, "n_representations": len(reps),
            "n_blocs_4_ou_plus": len(blocs), "longueur_bloc_max": plus_long,
            "blocs": ";".join(f"{p}:{ref}>{alt}" for p, ref, alt in sorted(blocs)) or "-",
        })

    out = Path(a.sortie)
    out.parent.mkdir(parents=True, exist_ok=True)
    champs = list(lignes[0])
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=champs, delimiter="\t")
        w.writeheader()
        w.writerows(lignes)

    sans = sorted(l["acide_amine"] for l in lignes if l["n_blocs_4_ou_plus"] == 0)
    syn = [
        ("variants_du_codon", len(lignes), "variants protéiques distincts au codon"),
        ("variants_a_entree_combinee", n_combine, "au moins un bloc de 4 bases ou plus"),
        ("variants_a_bloc_cinq_bases", n_cinq, "bloc de 5 bases (position et forme des manquants)"),
        ("variants_sans_entree_combinee", len(lignes) - n_combine, ", ".join(sans)),
        ("lignes_du_codon", sum(l["n_representations"] for l in lignes), "toutes représentations"),
    ]
    syn_path = out.with_name(out.stem + "_synthese.tsv")
    with open(syn_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["indicateur", "valeur", "detail"])
        w.writerows(syn)

    for l in lignes:
        print(f"  {l['variant']:22s} {l['n_representations']:2d} rep, bloc max "
              f"{l['longueur_bloc_max']} : {l['blocs']}")
    print()
    for k, v, d in syn:
        print(f"{k:32s} {v:5d}  {d}")
    assert set(sans) == set(TROIS) if (a.gene, a.codon) == ("gyrA", 94) else True, \
        f"les variants sans entrée combinée devraient être {TROIS}, lus : {sans}"
    print(f"\n-> {out} ; {syn_path}")


if __name__ == "__main__":
    main()
