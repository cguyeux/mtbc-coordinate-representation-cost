#!/usr/bin/env python3
"""Règles d'interprétation de tbdb 2026 (colonne `comment`) : épistasie Rv0678/mmpL5 et eis.

Objet : porter dans le catalogue consolidé les règles que l'OMS 2e éd. (reprises par tbdb 2026) attachent à
  certains déterminants, et fournir l'application mécanique de ces règles à une baseline.
  Quatre étiquettes, posées dans la colonne `regle` du catalogue :
    LIE_mmpL5     Rv0678, BDQ/CFZ : « ne confère la résistance que si lié à un MmpL5 fonctionnel » ;
    ABROGE_Rv0678 mmpL5 (perte de fonction) : « abroge l'effet des mutations Rv0678 liées » ;
    LIE_eis       eis, AMK/KAN : « ne confère la résistance que si lié à un Eis fonctionnel » ;
    ABROGE_eis    eis (perte de fonction) : « abroge l'effet des mutations eis liées ».
  Règle appliquée (`COUPLES`) : chez une souche qui porte un variant ABROGEANT, les déterminants R du gène soumis
  (Rv0678 pour BDQ/CFZ, eis pour AMK/KAN) ne comptent plus ; les déterminants des autres gènes restent.
  « Lié » est lu comme « porté par la même souche » (génome haploïde clonal) ; un mélange de souches n'est pas distingué.
Entrées : catalogue (DataFrame ou TSV), data/sources/tbprofiler_2026/mutations.csv (colonne `comment`, tous types de lignes).
Sorties : fonctions `annoter_catalogue`, `spdi_regle`, `couple_regle` ; aucune écriture de fichier (le catalogue est écrit par
  `phase1b_tbprofiler2026_v2.py`, la baseline par `phase11_baseline_audit.py`, la mesure de l'effet par `aa_regles_epistasie.py`).
Réutilisable : oui pour toute règle portée par `comment` (ajouter une ligne à REGLES et, si elle neutralise un gène, à COUPLES).
Projet : Resistance_antibio. Date : 2026-11-11 (piste AA1).
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
from vocab import canonical_drug

TBDB = paths.SOURCES / "tbprofiler_2026" / "mutations.csv"

REGLES = [  # (étiquette, motif du commentaire, gène que doit citer le commentaire ou None)
    ("LIE_mmpL5", "functional MmpL5", None),
    ("ABROGE_Rv0678", "Abrogates effect", "Rv0678"),
    ("ABROGE_eis", "Abrogates effect", "eis"),
    ("LIE_eis", "functional Eis", None),
]
# médicament -> (gène dont les déterminants R sont neutralisés, étiquette du variant abrogeant)
COUPLES = {
    "bedaquiline": ("Rv0678", "ABROGE_Rv0678"),
    "clofazimine": ("Rv0678", "ABROGE_Rv0678"),
    "amikacin": ("eis", "ABROGE_eis"),
    "kanamycin": ("eis", "ABROGE_eis"),
}


def regles_tbdb(src=TBDB):
    """{(gène, mutation, médicament canonique): étiquette de règle}, lu dans la colonne `comment`."""
    out = {}
    with open(src, newline="") as f:
        for x in csv.DictReader(f):
            c = x["comment"]
            for nom, motif, gene in REGLES:
                if motif in c and (gene is None or gene in c):
                    out[(x["Gene"], x["Mutation"], canonical_drug(x["drug"])[0])] = nom
                    break
    return out


def annoter_catalogue(cat, src=TBDB):
    """Copie de `cat` avec la colonne `regle` (vide quand aucune règle ne s'applique à la ligne).

    Jointure sur (gène, variante sans le préfixe du gène, médicament) : les lignes propagées aux codons synonymes
    héritent de la règle de leur variant protéique."""
    regles = regles_tbdb(src)
    out = cat.copy()

    def etiq(r):
        v = str(r["variant"])
        if "_" not in v:
            return ""
        return regles.get((r["gene"], v.split("_", 1)[1], r["drug"]), "")

    out["regle"] = out.apply(etiq, axis=1)
    return out


def spdi_regle(cat, etiquette, drug=None):
    """Ensemble des SPDI résolus des lignes portant `etiquette` (colonne `regle`), éventuellement pour un médicament."""
    d = cat[(cat["regle"] == etiquette) & (cat["spdi"] != "")]
    if drug is not None:
        d = d[d["drug"] == drug]
    return set(d["spdi"])


def couple_regle(cat, drug):
    """Pour un médicament soumis à une règle d'abrogation : (SPDI R du gène soumis, SPDI R des autres gènes, SPDI abrogeants).

    None quand le médicament n'a pas de règle. Les SPDI R sont ceux des lignes `R-associated` du médicament."""
    if drug not in COUPLES:
        return None
    gene, abro = COUPLES[drug]
    det = cat[(cat["drug"] == drug) & (cat["call"] == "R-associated") & (cat["spdi"] != "")]
    soumis = set(det[det["gene"] == gene]["spdi"])
    autres = set(det[det["gene"] != gene]["spdi"])
    return soumis, autres, spdi_regle(cat, abro, drug)
