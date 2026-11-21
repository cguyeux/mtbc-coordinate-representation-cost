#!/usr/bin/env python3
"""P5.6 — Cible supervisée de RÉFÉRENCE des rapports de `results/model/`.

Objet      : définir en un seul endroit quelles étiquettes phénotypiques (R/S) entrent dans un
             rapport de résultats, selon la SOURCE retenue au consensus pour chaque paire
             (souche, drogue). Trois cibles, nommées, que `phase5`, `phase10`,
             `phase11_baseline_audit.py` et `p52_multidrogue.py` partagent :
               complete   toutes les étiquettes de `rows.tsv` (la cible d'avant P5.6) ;
               reference  sans `souches_datees` (CIBLE PAR DÉFAUT) : ses étiquettes INH/RIF/EMB/PZA
                          sont dérivées du génotype, donc circulaires face au catalogue (P5.2,
                          cahier 2026-10-18). PATRIC en fait partie depuis P8.5 : son défaut
                          d'encodage est corrigé (P8.4) et sa provenance auditée (dilution en bouillon
                          ou sur gélose, appariement vers SRA à 97 %) ; la règle de conservation d'une
                          source est sa provenance, jamais son accord avec le catalogue ;
               mesuree    sans `souches_datees` ni PATRIC : analyse de SENSIBILITÉ hors PATRIC, pour
                          montrer ce qui dépend de cette source (MXF, AMK, CAP, STM). Le nom date de
                          P5.6, quand PATRIC passait pour suspect ; il est gardé pour la continuité des
                          fichiers (`baseline_audit.tsv`, `p74_angles_morts.tsv`) ; ne pas le lire
                          comme « seules les mesures sont bonnes » ;
             La règle est de SOURCE, jamais de génotype (même principe que `cible_pza.py`) : écarter
             d'après le génotype fabriquerait la circularité qu'on veut retirer.
             REPLI (P8.6) : le consensus ne garde qu'une source par paire (souche, drogue), la plus
             prioritaire. Quand cette source est écartée par la cible, la paire n'est plus perdue :
             son étiquette est celle de la source suivante de la priorité, lue dans la table brute des
             mesures (2 757 paires PATRIC d'INH, RIF, EMB et PZA, masquées jusque-là par `souches_datees`).
             L'étiquette sous une cible peut donc différer de celle de `rows.tsv` : les appelants lisent
             `etiquettes_cible`, jamais `rows[drug]`, et `masque_cible` n'en est que la vue booléenne.
Entrées    : `rows.tsv` de la matrice (colonne `strain_id`, une colonne par drogue valant R/S/vide) ;
             `catalogue/phenotypes_souches.tsv` (toutes les mesures : strain_id, drug, phenotype, source) ;
             ordre de priorité de `phase2_build_phenotypes.PRIORITY` (définition unique du consensus).
Sorties    : étiquettes R/S ou masques booléens alignés sur les lignes de `rows.tsv` (aucun fichier écrit).
Réutilisable : oui — toute table de mesures phénotypiques avec une colonne de source et un ordre de priorité.
Projet     : Resistance_antibio (P5.6, repli de source P8.6)
Date       : 2026-10-18 (docstring révisé en P8.5, repli en P8.6)
Usage      : from cible_reference import etiquettes_cible, masque_cible, CIBLES, CIBLE_DEFAUT
             lab = etiquettes_cible(rows, "isoniazid", "reference")   # R/S/'' sous cette cible
             ym = masque_cible(rows, "isoniazid", "reference")        # lignes étiquetées sous la cible
             python3 analyses/cible_reference.py                      # effectifs par cible et paires récupérées
"""
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
from phase2_build_phenotypes import PRIORITY

CIBLES = {
    "complete": (),
    "reference": ("souches_datees",),
    "mesuree": ("souches_datees", "PATRIC"),
}
CIBLE_DEFAUT = "reference"
LIBELLES = {
    "complete": "complète (toutes sources)",
    "reference": "référence (sans souches_datees)",
    "mesuree": "hors PATRIC (sans souches_datees ni PATRIC, analyse de sensibilité)",
}


def sources_cible(cible):
    """Sources écartées par une cible nommée, ou par une liste séparée par des virgules."""
    if cible in CIBLES:
        return CIBLES[cible]
    sources = tuple(s for s in str(cible).split(",") if s)
    if not sources:
        raise ValueError(f"cible inconnue : {cible!r} (attendu : {', '.join(CIBLES)} ou une liste de sources)")
    return sources


@lru_cache(maxsize=1)
def _mesures():
    """Toutes les mesures R/S de `phenotypes_souches.tsv` (une ligne par souche, drogue et source)."""
    mes = pd.read_csv(paths.PHENOTYPES_TSV, sep="\t", usecols=["strain_id", "drug", "phenotype", "source"],
                      keep_default_na=False)
    return mes[mes["phenotype"].isin(["R", "S"])]


def _calcule_retenu(sources, mesures):
    """{drogue: DataFrame indexé par souche (phenotype, source)} : mesure la plus prioritaire de chaque
    paire une fois `sources` écartées. Même règle que `phase2.calculer_consensus` (tri stable, première
    ligne gardée) ; avec `sources` vide, c'est le consensus."""
    rang = {s: i for i, s in enumerate(PRIORITY)}
    m = mesures[~mesures["source"].isin(sources)]
    m = m.assign(_rang=m["source"].map(lambda s: rang.get(s, 99))).sort_values("_rang", kind="stable")
    m = m.drop_duplicates(["strain_id", "drug"], keep="first").drop(columns="_rang")
    return {d: g.set_index("strain_id")[["phenotype", "source"]] for d, g in m.groupby("drug")}


@lru_cache(maxsize=8)
def _retenu_defaut(sources):
    return _calcule_retenu(sources, _mesures())


def _retenu(sources, mesures=None):
    if mesures is None:
        return _retenu_defaut(tuple(sorted(sources)))
    return _calcule_retenu(sources, mesures)


def source_par_ligne(rows, drug, mesures=None):
    """Source retenue au CONSENSUS pour la drogue, alignée sur `rows` (chaîne vide si absente)."""
    cons = _retenu((), mesures).get(drug)
    if cons is None:
        return pd.Series("", index=rows.index)
    return rows["strain_id"].map(cons["source"]).fillna("")


def etiquettes_cible(rows, drug, cible=CIBLE_DEFAUT, mesures=None, repli=True):
    """Étiquette R/S de `drug` sous `cible`, alignée sur `rows` ('' quand la paire n'a plus d'étiquette).

    Avec `repli` (défaut), une paire dont la source de consensus est écartée prend l'étiquette de la
    source suivante non écartée ; elle n'est perdue que si toutes ses sources le sont. `repli=False` rejoue
    la règle d'avant P8.6 (paire perdue dès que la source de consensus est écartée), pour mesurer ce que le
    repli récupère. Une étiquette absente de `rows` (nettoyage PZA, P5.3) n'est jamais rétablie ; une
    étiquette de `rows` sans aucune mesure dans la table est gardée (source inconnue, donc non écartable)."""
    brut = rows[drug].where(rows[drug].isin(["R", "S"]), "")
    sources = sources_cible(cible)
    if not sources:
        return brut
    sid = rows["strain_id"]
    cons = _retenu((), mesures).get(drug)
    connu = sid.isin(cons.index) if cons is not None else pd.Series(False, index=rows.index)
    if repli:
        ret = _retenu(sources, mesures).get(drug)
        lab = sid.map(ret["phenotype"]).fillna("") if ret is not None else pd.Series("", index=rows.index)
    else:
        src = sid.map(cons["source"]).fillna("") if cons is not None else pd.Series("", index=rows.index)
        lab = brut.where(~src.isin(sources), "")
    return lab.where(connu, brut).where(brut != "", "")


def source_cible_par_ligne(rows, drug, cible=CIBLE_DEFAUT, mesures=None):
    """Source dont vient l'étiquette sous `cible` (chaîne vide si la paire est perdue ou de source inconnue)."""
    sources = sources_cible(cible)
    ret = _retenu(sources, mesures).get(drug)
    lab = etiquettes_cible(rows, drug, cible, mesures)
    if ret is None:
        return pd.Series("", index=rows.index)
    return rows["strain_id"].map(ret["source"]).fillna("").where(lab != "", "")


def masque_cible(rows, drug, cible=CIBLE_DEFAUT, mesures=None, repli=True):
    """Lignes de `rows` dont l'étiquette R/S de `drug` est gardée sous `cible`."""
    return etiquettes_cible(rows, drug, cible, mesures, repli).isin(["R", "S"]).values


def effectifs(rows, drug):
    """{cible: (n_R, n_S)} pour une drogue, sur les trois cibles nommées."""
    out = {}
    for nom in CIBLES:
        lab = etiquettes_cible(rows, drug, nom)
        out[nom] = (int((lab == "R").sum()), int((lab == "S").sum()))
    return out


DROGUES = ["isoniazid", "rifampicin", "ethambutol", "pyrazinamide", "moxifloxacin", "levofloxacin",
           "amikacin", "kanamycin", "streptomycin", "ethionamide", "capreomycin"]


def main():
    rows = pd.read_csv(paths.RESULTATS / "model" / "rows.tsv", sep="\t", keep_default_na=False)
    print("drogue\t" + "\t".join(f"{c}_R\t{c}_S" for c in CIBLES))
    for drug in DROGUES:
        e = effectifs(rows, drug)
        print(drug + "\t" + "\t".join(f"{e[c][0]}\t{e[c][1]}" for c in CIBLES))
    print("\npaires récupérées par le repli de source (P8.6), par cible :")
    for nom in CIBLES:
        gain = {d: int(masque_cible(rows, d, nom).sum() - masque_cible(rows, d, nom, repli=False).sum())
                for d in DROGUES}
        detail = ", ".join(f"{d} {n}" for d, n in gain.items() if n)
        print(f"  {nom}: {detail or 'aucune'}")


if __name__ == "__main__":
    main()
