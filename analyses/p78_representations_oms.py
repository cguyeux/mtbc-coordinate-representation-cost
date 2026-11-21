#!/usr/bin/env python3
"""P7.8 — Comptage des représentations de coordonnées du fichier OMS et de celles que la chaîne retient.

Objet      : donner un TSV aux chiffres du manuscrit recadré (A1) qui portent sur la SOURCE elle-même
             (combien de représentations le fichier de coordonnées de l'OMS fournit par variant) et sur
             la PERTE de l'ingestion (combien le catalogue consolidé en retient). Les chiffres du gain de
             sensibilité, eux, viennent de `p77_C0.tsv` / `p77_C1.tsv` (P7.7) et ne sont pas recalculés ici.
Entrées    : data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt (1-based),
             catalogue/catalogue_consolide.tsv (catalogue vivant, dit C0).
Sorties    : résultats/model/p78_representations_oms.tsv (indicateur, valeur, detail).
Réutilisable : oui pour un autre catalogue HGVS fourni avec un fichier de coordonnées multi-lignes.
Projet     : Resistance_antibio (P7.8)
Date       : 2026-11-16
Usage      : python3 analyses/p78_representations_oms.py
Coût       : lecture de deux TSV (7,8 Mo et 13 Mo), quelques secondes sur un cœur, < 1 Gio.
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
from hgvs_to_spdi import resolve_all

COORD = paths.SRC_WHO / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
SORTIE = paths.ROOT / "résultats" / "model" / "p78_representations_oms.tsv"
MAX_LEN = 24  # filtre de `resolve_all`, donc de ce que la chaîne peut émettre


def nature(ref, alt):
    """SNP, MNV (bloc de même longueur), indel, ou `long` au-delà du filtre de resolve_all."""
    if len(ref) > MAX_LEN or len(alt) > MAX_LEN:
        return "long"
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) == len(alt):
        return "mnv"
    return "indel"


def main():
    df = pd.read_csv(COORD, sep="\t", dtype=str)
    reps = defaultdict(set)          # variant -> {(pos0, ref, alt)}
    par_nature = defaultdict(int)
    for v, pos, ref, alt in zip(df["variant"], df["position"],
                                df["reference_nucleotide"], df["alternative_nucleotide"]):
        if not (isinstance(v, str) and isinstance(ref, str) and isinstance(alt, str)):
            continue
        reps[v].add((int(pos) - 1, ref, alt))
        par_nature[nature(ref, alt)] += 1
    lignes = sum(par_nature.values())
    plusieurs = sum(1 for v, s in reps.items() if len(s) > 1)
    prot = [v for v in reps if "_p." in v]
    prot_plusieurs = sum(1 for v in prot if len(reps[v]) > 1)

    out = [("lignes_fichier_oms", lignes, "lignes à variant, position, ref et alt lisibles"),
           ("variants_fichier_oms", len(reps), ""),
           ("representations_par_variant", round(lignes / len(reps), 2), "moyenne"),
           ("variants_a_plusieurs_representations", plusieurs, f"{100 * plusieurs / len(reps):.1f} %"),
           ("variants_proteiques", len(prot), "notation _p."),
           ("variants_proteiques_a_plusieurs_representations", prot_plusieurs,
            f"{100 * prot_plusieurs / len(prot):.1f} %")]
    for nat in ("snp", "mnv", "indel", "long"):
        out.append((f"lignes_oms_{nat}", par_nature[nat],
                    f"{100 * par_nature[nat] / lignes:.1f} %" if nat == "long" else ""))

    # Perte de l'ingestion : sur les variants du catalogue vivant retrouvés dans le fichier OMS,
    # ce que `resolve_all` rend contre ce que `spdis_for` a retenu (une représentation par variant
    # dès qu'une forme mono-base existe).
    cat = pd.read_csv(paths.CATALOGUE_TSV, sep="\t", dtype=str)
    retenues = defaultdict(set)
    for v, spdi in zip(cat["variant"], cat["spdi"]):
        if isinstance(v, str) and isinstance(spdi, str) and spdi:
            retenues[v].add(spdi)
    communs = [v for v in retenues if v in reps]
    fournies = sum(len(resolve_all(v)) for v in communs)
    gardees = sum(len(retenues[v]) for v in communs)
    perdues_blocs = 0
    perdant = 0
    for v in communs:
        toutes = set(resolve_all(v))
        manquantes = toutes - retenues[v]
        if manquantes:
            perdant += 1
        for s in manquantes:
            _, _, ref, alt = s.split(":")
            if nature(ref, alt) == "mnv":
                perdues_blocs += 1
    # Ce que l'émission complète ajouterait comme déterminants R en bloc : les SPDI de même longueur
    # ref et alt que `resolve_all` rend pour les variants dont une assertion est `R-associated`.
    var_r = {v for v, c in zip(cat["variant"], cat["call"])
             if isinstance(v, str) and c == "R-associated"}
    mnv_r = set()
    for v in var_r:
        for s in resolve_all(v):
            _, _, ref, alt = s.split(":")
            if nature(ref, alt) == "mnv":
                mnv_r.add(s)

    out += [("lignes_catalogue_C0", len(cat), "catalogue consolidé vivant"),
            ("medicaments_du_catalogue_C0", cat["drug"].nunique(), ""),
            ("medicaments_de_la_source_oms", cat[cat["source"] == "WHO_catalogue"]["drug"].nunique(),
             "médicaments distincts des assertions de source OMS"),
            ("variants_R_associes_du_catalogue", len(var_r), "au moins une assertion R-associated"),
            ("spdi_bloc_R_associes_de_l_emission_complete", len(mnv_r),
             "SPDI distincts de même longueur ref et alt"),
            ("variants_du_catalogue_retrouves_dans_oms", len(communs), ""),
            ("representations_fournies_par_oms", fournies, "resolve_all, 24 bases au plus"),
            ("representations_retenues_par_la_chaine", gardees, "catalogue consolidé vivant (C0)"),
            ("representations_perdues", fournies - gardees,
             f"{100 * (fournies - gardees) / fournies:.1f} %"),
            ("representations_perdues_en_bloc", perdues_blocs, "de même longueur ref et alt"),
            ("variants_perdant_au_moins_une_representation", perdant,
             f"{100 * perdant / len(communs):.1f} %")]

    SORTIE.parent.mkdir(parents=True, exist_ok=True)
    with open(SORTIE, "w") as f:
        f.write("indicateur\tvaleur\tdetail\n")
        for k, v, d in out:
            f.write(f"{k}\t{v}\t{d}\n")
    print(f"→ {SORTIE}")
    for k, v, d in out:
        print(f"  {k:52s} {v:>10} {d}")


if __name__ == "__main__":
    main()
