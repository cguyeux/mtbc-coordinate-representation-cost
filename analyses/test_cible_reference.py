#!/usr/bin/env python3
"""Tests de `cible_reference.py` (P5.6, repli de source P8.6) : règle de SOURCE, jamais de génotype.

Objet      : vérifier sur un jeu synthétique que chaque cible écarte exactement les étiquettes des
             sources prévues, que la cible `complete` ne touche à rien, qu'une paire (souche, drogue)
             sans mesure est gardée, qu'une drogue sans source exclue donne des masques identiques (cas
             que `phase5` et `phase10` ne recalculent pas), et que le REPLI (P8.6) lit la source suivante
             de la priorité au lieu de perdre la paire, sans jamais ressusciter une étiquette retirée de
             `rows` ni dépendre du phénotype.
Entrées    : aucune (jeu synthétique en mémoire).
Sorties    : code de retour 0 si tout passe.
Réutilisable : non (spécifique à `cible_reference`).
Projet     : Resistance_antibio (P5.6)
Date       : 2026-10-18 (repli de source ajouté en P8.6)
Usage      : python3 analyses/test_cible_reference.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cible_reference as cr


def jeu():
    rows = pd.DataFrame({"strain_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
                         "isoniazid": ["R", "S", "R", "S", "", "R"],
                         "amikacin": ["S", "S", "R", "", "S", "S"]})
    cons = pd.DataFrame([
        ("s1", "isoniazid", "R", "souches_datees"), ("s2", "isoniazid", "S", "CRyPTIC"),
        ("s3", "isoniazid", "R", "PATRIC"), ("s4", "isoniazid", "S", "Coll2018"),
        ("s6", "isoniazid", "R", "souches_datees"),                      # s5 : pas d'étiquette
        ("s1", "amikacin", "S", "CRyPTIC"), ("s2", "amikacin", "S", "CRyPTIC"),
        ("s3", "amikacin", "R", "PATRIC"), ("s5", "amikacin", "S", "CRyPTIC"),
        ("s6", "amikacin", "S", "CRyPTIC"),
    ], columns=["strain_id", "drug", "phenotype", "source"])
    return rows, cons


def test_cibles_nommees():
    rows, cons = jeu()
    garde = lambda d, c: list(cr.masque_cible(rows, d, c, mesures=cons))
    assert garde("isoniazid", "complete") == [True, True, True, True, False, True]
    assert garde("isoniazid", "reference") == [False, True, True, True, False, False]
    assert garde("isoniazid", "mesuree") == [False, True, False, True, False, False]
    # amikacine : aucune étiquette souches_datees, la cible de référence ne change rien
    assert garde("amikacin", "reference") == garde("amikacin", "complete")
    assert garde("amikacin", "mesuree") == [True, True, False, False, True, True]


def test_source_libre_et_absente():
    rows, cons = jeu()
    m = cr.masque_cible(rows, "isoniazid", "Coll2018,CRyPTIC", mesures=cons)
    assert list(m) == [True, False, True, False, False, True]
    # étiquette présente dans rows mais absente du consensus : gardée, source vide
    rows2 = rows.copy(); rows2.loc[4, "isoniazid"] = "S"
    assert cr.masque_cible(rows2, "isoniazid", "reference", mesures=cons)[4]
    assert cr.source_par_ligne(rows2, "isoniazid", cons).iloc[4] == ""


def test_regle_de_source_pas_de_genotype():
    """Le masque ne dépend que de (strain_id, drug, source) : changer le phénotype R/S ne le déplace pas."""
    rows, cons = jeu()
    base = cr.masque_cible(rows, "isoniazid", "reference", mesures=cons)
    inv = rows.copy(); inv["isoniazid"] = inv["isoniazid"].map({"R": "S", "S": "R", "": ""})
    assert np.array_equal(base, cr.masque_cible(inv, "isoniazid", "reference", mesures=cons))


def jeu_repli():
    """Deux sources par paire : souches_datees (prioritaire) puis PATRIC, plus un cas à trois sources."""
    rows = pd.DataFrame({"strain_id": ["a", "b", "c", "d", "e", "f"],
                         "isoniazid": ["R", "R", "S", "R", "", "S"]})
    mes = pd.DataFrame([
        ("a", "isoniazid", "R", "souches_datees"), ("a", "isoniazid", "S", "PATRIC"),   # étiquette différente
        ("b", "isoniazid", "R", "souches_datees"), ("b", "isoniazid", "R", "PATRIC"),   # même étiquette
        ("c", "isoniazid", "S", "souches_datees"),                                      # une seule source, écartée
        ("d", "isoniazid", "R", "Coll2018"), ("d", "isoniazid", "S", "souches_datees"),
        ("d", "isoniazid", "S", "PATRIC"),                                              # Coll2018 passe avant
        ("e", "isoniazid", "R", "souches_datees"), ("e", "isoniazid", "S", "PATRIC"),   # retirée de rows (PZA)
        ("f", "isoniazid", "S", "souches_datees"), ("f", "isoniazid", "R", "literature"),
        ("f", "isoniazid", "S", "PATRIC"),
    ], columns=["strain_id", "drug", "phenotype", "source"])
    return rows, mes


def test_repli_source_suivante():
    rows, mes = jeu_repli()
    lab = lambda c, **k: list(cr.etiquettes_cible(rows, "isoniazid", c, mesures=mes, **k))
    # complète : l'étiquette de rows, inchangée
    assert lab("complete") == ["R", "R", "S", "R", "", "S"]
    # référence : a prend l'étiquette de PATRIC (différente), b la même, c est perdue (souche sans autre
    # source), d garde Coll2018, e reste vide (retirée de rows), f prend PATRIC avant literature
    assert lab("reference") == ["S", "R", "", "R", "", "S"]
    # mesurée (sans souches_datees ni PATRIC) : a et b n'ont plus rien, f tombe sur literature
    assert lab("mesuree") == ["", "", "", "R", "", "R"]


def test_repli_faux_rejoue_ancienne_regle():
    rows, mes = jeu_repli()
    ancien = list(cr.etiquettes_cible(rows, "isoniazid", "reference", mesures=mes, repli=False))
    # source de consensus écartée : la paire est perdue entière, comme avant P8.6
    assert ancien == ["", "", "", "R", "", ""]
    assert list(cr.masque_cible(rows, "isoniazid", "reference", mesures=mes, repli=False)) == \
        [False, False, False, True, False, False]


def test_repli_ne_ressuscite_pas_une_etiquette_retiree():
    rows, mes = jeu_repli()
    # la souche e a des mesures (R souches_datees, S PATRIC) mais rows est vide (nettoyage PZA) : rien ne revient
    for cible in ("complete", "reference", "mesuree", "PATRIC"):
        assert cr.etiquettes_cible(rows, "isoniazid", cible, mesures=mes).iloc[4] == ""


def test_repli_regle_de_source_pas_de_genotype():
    """Le repli ne lit que (strain_id, drug, source, priorité) : inverser R et S partout inverse les
    étiquettes sous la cible sans déplacer le masque."""
    rows, mes = jeu_repli()
    sw = {"R": "S", "S": "R", "": ""}
    rows_i = rows.assign(isoniazid=rows["isoniazid"].map(sw))
    mes_i = mes.assign(phenotype=mes["phenotype"].map(sw))
    a = cr.etiquettes_cible(rows, "isoniazid", "reference", mesures=mes)
    b = cr.etiquettes_cible(rows_i, "isoniazid", "reference", mesures=mes_i)
    assert list(a.map(sw)) == list(b)
    assert np.array_equal(cr.masque_cible(rows, "isoniazid", "reference", mesures=mes),
                          cr.masque_cible(rows_i, "isoniazid", "reference", mesures=mes_i))


def test_priorite_unique_et_source_inconnue():
    """L'ordre vient de phase2 (définition unique) ; une source hors liste passe après toutes les autres."""
    import phase2_build_phenotypes as p2
    assert cr.PRIORITY is p2.PRIORITY
    rows = pd.DataFrame({"strain_id": ["x"], "isoniazid": ["R"]})
    mes = pd.DataFrame([("x", "isoniazid", "R", "source_inconnue"), ("x", "isoniazid", "S", "literature")],
                       columns=["strain_id", "drug", "phenotype", "source"])
    assert cr.etiquettes_cible(rows, "isoniazid", "souches_datees", mesures=mes).iloc[0] == "S"


def test_source_cible_par_ligne():
    rows, mes = jeu_repli()
    src = list(cr.source_cible_par_ligne(rows, "isoniazid", "reference", mesures=mes))
    assert src == ["PATRIC", "PATRIC", "", "Coll2018", "", "PATRIC"]
    # source de consensus (avant repli) : souches_datees pour a, b, c, e, f
    assert list(cr.source_par_ligne(rows, "isoniazid", mesures=mes)) == \
        ["souches_datees", "souches_datees", "souches_datees", "Coll2018", "souches_datees", "souches_datees"]


def test_cible_inconnue():
    try:
        cr.sources_cible("")
    except ValueError:
        return
    raise AssertionError("une cible vide doit être refusée")


if __name__ == "__main__":
    for nom, f in sorted(globals().items()):
        if nom.startswith("test_") and callable(f):
            f()
            print("ok", nom)
