#!/usr/bin/env python3
"""P7.12 — Test de non-régression : chaque chiffre épinglé est IMPRIMÉ dans le manuscrit.

Objet      : `test_p78_chiffres_manuscrit.py` vérifie que les valeurs recopiées dans le test sont celles des
             TSV ; il n'ouvre jamais `article/main.tex`. Ce test-ci ferme l'écart : chaque affirmation épinglée
             doit figurer dans `main.tex` ou dans `supplementary.tex`, sauf si elle est déclarée RETIRÉE
             ci-dessous avec son motif. Un chiffre emporté par une coupe fait donc échouer ce test, là où
             l'autre reste vert.
             Les petits entiers (moins de 100) ne sont pas contrôlés : ils se retrouvent partout dans un texte
             et leur présence ne prouve rien (voir `p712_chiffres_dans_texte.py`, § Limites).
             La liste RETIRES est fermée : un nombre qui y figure doit être réellement absent (sinon la ligne
             est périmée et doit en sortir), et tout autre nombre absent fait échouer le test.
Entrées    : analyses/test_p78_chiffres_manuscrit.py (CLAIMS), article/main.tex, article/supplementary.tex.
Sorties    : aucune (code de retour pytest).
Réutilisable : oui, comme patron.
Projet     : Resistance_antibio (P7.12)
Date       : 2026-11-21
Usage      : python3 -m pytest -q analyses/test_p78_chiffres_dans_texte.py
"""
import pytest

import p712_chiffres_dans_texte as ctl

MOTIF_QUATRE_DECIMALES = ("liste de sensibilités à quatre décimales de § 3.4 remplacée par le tableau 2 "
                          "(trois décimales) et un renvoi, réduction de longueur P7.12 (2026-11-21)")
MOTIF_DISCUSSION = ("paragraphe « The cost side favours decomposition » de la Discussion retiré, redite de "
                    "§ 3.5 (2026-11-21, P7.12)")

RETIRES = {
    "pyrazinamide sens. complète C0 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "ethionamide sens. complète C0 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "ethionamide sens. complète C1 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "capreomycin sens. complète C0 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "capreomycin sens. complète C1 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "isoniazid sens. complète C0 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "isoniazid sens. complète C1 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "streptomycin sens. complète C0 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "streptomycin sens. complète C1 (4 déc.)": MOTIF_QUATRE_DECIMALES,
    "facteur des SPDI R-associés distincts": MOTIF_DISCUSSION,
    "colonnes gagnées": MOTIF_DISCUSSION,
}

CLAIMS = ctl.charger_claims()
TEXTE = ctl.lire(ctl.ARTICLE / "main.tex") + "\n" + ctl.lire(ctl.ARTICLE / "supplementary.tex", facultatif=True)
CONTROLES = [c for c in CLAIMS if not ctl.faible(c[1])]


@pytest.mark.parametrize("nom,imprime,_lu,_demi", CONTROLES, ids=[c[0] for c in CONTROLES])
def test_chiffre_present_dans_le_texte(nom, imprime, _lu, _demi):
    if nom in RETIRES:
        pytest.skip(f"déclaré retiré : {RETIRES[nom]}")
    assert ctl.present(imprime, TEXTE), (
        f"« {nom} » ({imprime}) est épinglé mais n'est plus imprimé dans main.tex ni supplementary.tex : "
        "le rétablir, ou le déclarer dans RETIRES avec son motif")


@pytest.mark.parametrize("nom", sorted(RETIRES))
def test_retire_est_bien_absent(nom):
    imprime = next(c[1] for c in CLAIMS if c[0] == nom)
    assert not ctl.present(imprime, TEXTE), (
        f"« {nom} » ({imprime}) est de nouveau imprimé : le retirer de RETIRES")
