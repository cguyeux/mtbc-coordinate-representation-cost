#!/usr/bin/env python3
"""P7.7 — Tests du drapeau `--blocs` (émission des blocs OMS) de `phase1_build_catalogue.spdis_for()`.

Objet      : (1) le défaut reproduit l'ancien comportement (mono-base seule dès qu'elle existe) ; (2) `EMETTRE_BLOCS`
             émet toutes les représentations de `resolve_all`, dont le bloc `ACAG>GCAC` de gyrA Asp94Gly (cahier
             du 2026-11-16, P7.9) ; (3) le défaut est toujours un sous-ensemble du mode blocs, sur tous les variants
             du fichier OMS ; (4) un variant qui n'a qu'une forme (mono-base ou indel) n'est pas modifié.
Entrées    : data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt.
Réutilisable : non (propre à P7.7).
Projet     : Resistance_antibio (P7.7)
Date       : 2026-11-16
Usage      : python3 -m pytest analyses/test_p77.py -q
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import phase1_build_catalogue as p1
from hgvs_to_spdi import _index, resolve_all


@pytest.fixture(autouse=True)
def remettre_defaut():
    p1.EMETTRE_BLOCS = False
    yield
    p1.EMETTRE_BLOCS = False


def est_snp(s):
    _, _, ref, alt = s.split(":")
    return len(ref) == 1 and len(alt) == 1


def test_defaut_ecarte_les_blocs_si_une_forme_mono_existe():
    sp = p1.spdis_for("gyrA_p.Asp94Gly")
    assert sp and all(est_snp(s) for s in sp)
    assert "NC_000962.3:7581:ACAG:GCAC" not in sp


def test_blocs_emet_le_bloc_oms_de_gyrA_asp94gly():
    p1.EMETTRE_BLOCS = True
    sp = p1.spdis_for("gyrA_p.Asp94Gly")
    assert "NC_000962.3:7581:ACAG:GCAC" in sp
    assert set(sp) == set(resolve_all("gyrA_p.Asp94Gly"))


def test_defaut_est_un_sous_ensemble_du_mode_blocs_sur_tout_le_fichier_oms():
    _, allmap = _index()
    defaut, blocs = {}, {}
    for v in allmap:
        p1.EMETTRE_BLOCS = False
        defaut[v] = set(p1.spdis_for(v))
        p1.EMETTRE_BLOCS = True
        blocs[v] = set(p1.spdis_for(v))
    assert all(defaut[v] <= blocs[v] for v in allmap)
    # le mode blocs ajoute bien quelque chose, et beaucoup (le catalogue triple : 64 052 -> 204 522 lignes)
    assert sum(len(blocs[v]) - len(defaut[v]) for v in allmap) > 50_000


def test_variant_a_forme_unique_inchange():
    _, allmap = _index()
    uniques = [v for v in allmap if len(resolve_all(v)) == 1][:200]
    assert uniques
    for v in uniques:
        p1.EMETTRE_BLOCS = False
        a = p1.spdis_for(v)
        p1.EMETTRE_BLOCS = True
        b = p1.spdis_for(v)
        assert a == b


def test_variant_inconnu_rend_la_chaine_vide_dans_les_deux_modes():
    for mode in (False, True):
        p1.EMETTRE_BLOCS = mode
        assert p1.spdis_for("geneInexistant_c.1A>T") == [""]
