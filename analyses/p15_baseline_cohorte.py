#!/usr/bin/env python3
"""Baseline catalogue à COHORTE CONTRÔLÉE, pour séparer l'effet de la cohorte de celui des features (P1.5).

Objet      : P1.3 a fait passer la cible de 52 496 à 88 017 souches, dont 239 862 paires nouvelles
             venant à 99,9 % de `CRyPTIC_v3`. La sensibilité globale de `phase11` sur la nouvelle
             cible n'est donc PAS comparable à celle d'avant : deux choses ont bougé en même temps,
             la cohorte jugée et les colonnes de la matrice. Ce script rend le baseline du MÊME
             catalogue et de la MÊME matrice sur trois ensembles de souches : celles de l'ancienne
             matrice (`--subset`, la cohorte d'avant), leur complément (les souches ajoutées), et
             tout. Un écart entre les deux premiers mesure l'effet des features ; un écart entre
             les deux derniers mesure l'effet de la cohorte, et lui seul explique une sensibilité
             globale qui bouge sans qu'aucune assertion du catalogue ait changé.
             Même logique que `phase5b_cohorte_egale.py` pour les AUC, transposée au baseline
             déterministe : règle du projet, deux chiffres ne se comparent qu'à cohorte contrôlée.
Entrées    : `results/model/{X.npz, features.txt, rows.tsv}`, `catalogue/catalogue_consolide.tsv`,
             `--subset` (un identifiant de souche par ligne ; défaut : `rows_avant_P1.5.tsv` s'il existe).
Sorties    : `results/model/p15_baseline_cohorte.tsv` (une ligne par médicament et par ensemble,
             avec n_R, n_S, sensibilité, spécificité, sous les trois cibles de `cible_reference`).
Réutilisable : oui pour tout découpage de cohorte d'une même matrice (par étude, par pays, par source
             d'étiquette) : `baseline_sur` prend un masque booléen quelconque.
Projet     : Resistance_antibio (P1.5)
Date       : 2026-11-11
Usage      : python3 analyses/p15_baseline_cohorte.py [--subset ids.txt] [--etiquette-subset ancienne]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
import cible_reference as cr
import phase11_baseline_audit as p11

MODELDIR = paths.RESULTATS / "model"


def baseline_sur(rows, porte, drug, masque_cohorte):
    """n_R, n_S, sensibilité et spécificité du baseline sur les souches de `masque_cohorte`,
    sous chacune des trois cibles ; `porte` est le vecteur 0/1 « porte un déterminant R catalogué »."""
    out = {}
    for nom in cr.CIBLES:
        lab = cr.etiquettes_cible(rows, drug, nom)
        ym = lab.isin(["R", "S"]).values & masque_cohorte
        y = (lab[ym] == "R").astype(int).values
        car = porte[ym]
        tp = int((y & car).sum()); fn = int((y & (1 - car)).sum())
        tn = int(((1 - y) & (1 - car)).sum()); fp = int(((1 - y) & car).sum())
        out[nom] = {"n_R": int(y.sum()), "n_S": int((1 - y).sum()),
                    "sens": round(tp / (tp + fn), 3) if tp + fn else None,
                    "spec": round(tn / (tn + fp), 3) if tn + fp else None}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--subset", help="fichier d'identifiants de souches (défaut : rows_avant_P1.5.tsv)")
    ap.add_argument("--etiquette-subset", default="cohorte_ancienne")
    ap.add_argument("--etiquette-complement", default="souches_ajoutees")
    ap.add_argument("--paires-de", help="ancien rows.tsv : restreint AUSSI aux paires (souche, médicament) "
                                        "qui y portaient déjà une étiquette R/S. Sans cette option, une "
                                        "souche ancienne compte avec ses étiquettes COURANTES, plus "
                                        "nombreuses, et la comparaison reste confondue (mesuré sur la "
                                        "streptomycine : 886 cas R avant contre 2 316 sur les mêmes souches).")
    ap.add_argument("--sortie", default=str(MODELDIR / "p15_baseline_cohorte.tsv"))
    a = ap.parse_args()

    X = sparse.load_npz(MODELDIR / "X.npz").tocsr()
    feats = (MODELDIR / "features.txt").read_text().splitlines()
    rows = pd.read_csv(MODELDIR / "rows.tsv", sep="\t", keep_default_na=False)
    cat = pd.read_csv(paths.CATALOGUE_TSV, sep="\t").fillna("")
    sig = p11.signatures(feats)

    if a.subset:
        ids = set(l.strip() for l in open(a.subset) if l.strip())
    else:
        ancien = MODELDIR / "rows_avant_P1.5.tsv"
        if not ancien.exists():
            raise SystemExit(f"ENTRÉE MANQUANTE : {ancien} (ni --subset donné) ; sans la liste des souches "
                             "de l'ancienne matrice, la cohorte de référence est indéfinissable et un "
                             "baseline « à cohorte contrôlée » n'aurait aucun sens. Fournir --subset.")
        ids = set(pd.read_csv(ancien, sep="\t", keep_default_na=False, usecols=["strain_id"])["strain_id"])
    dans = rows["strain_id"].isin(ids).values
    print(f"souches : {len(rows)} au total, {int(dans.sum())} dans la cohorte de référence, "
          f"{int((~dans).sum())} ajoutées")

    paires_anciennes = None
    if a.paires_de:
        chemin = Path(a.paires_de)
        if not chemin.exists():
            raise SystemExit(f"ENTRÉE MANQUANTE : {chemin} (--paires-de) ; sans l'ancien rows.tsv, les "
                             "paires (souche, médicament) d'alors sont inconnues et la restriction demandée "
                             "serait silencieusement sans effet.")
        anc = pd.read_csv(chemin, sep="\t", keep_default_na=False)
        paires_anciennes = {d: set(anc.loc[anc[d].isin(["R", "S"]), "strain_id"])
                            for d in p11.DRUGS if d in anc.columns}
        print(f"paires de référence lues dans {chemin.name} : "
              + ", ".join(f"{d} {len(s)}" for d, s in sorted(paires_anciennes.items())))

    lignes = []
    for drug in p11.DRUGS:
        if not rows[drug].isin(["R", "S"]).any():
            continue
        porte, n_cols = p11.porteurs(X, feats, sig, cat, drug)
        if paires_anciennes is not None:
            # paire fixée : souche ancienne ET déjà étiquetée sur CE médicament à l'époque
            fixe = rows["strain_id"].isin(paires_anciennes.get(drug, set())).values
            ensembles = [(a.etiquette_subset, dans & fixe), (a.etiquette_complement, ~fixe),
                         ("tout", np.ones(len(rows), dtype=bool))]
        else:
            ensembles = [(a.etiquette_subset, dans), (a.etiquette_complement, ~dans),
                         ("tout", np.ones(len(rows), dtype=bool))]
        for nom_ens, masque in ensembles:
            mes = baseline_sur(rows, porte, drug, masque)
            for nom_cible, v in mes.items():
                lignes.append({"drug": drug, "cohorte": nom_ens, "cible": nom_cible,
                               "det_features": n_cols, **v})
    df = pd.DataFrame(lignes)
    df.to_csv(a.sortie, sep="\t", index=False)
    ref = df[df["cible"] == cr.CIBLE_DEFAUT].pivot(index="drug", columns="cohorte",
                                                   values=["n_R", "sens"])
    print(f"\ncible de référence ({cr.CIBLE_DEFAUT}) :")
    print(ref.to_string())
    print(f"\n-> {a.sortie}")


if __name__ == "__main__":
    main()
