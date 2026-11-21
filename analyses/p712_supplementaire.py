#!/usr/bin/env python3
"""P7.12 — Génère les deux tableaux supplémentaires d'A1 depuis leurs TSV, sans saisie à la main.

Objet      : lors de la réduction du corps d'A1 sous le plafond de la porte de longueur (P7.12), deux
             détails chiffrés ont été déportés VERBATIM du corps vers le supplémentaire, chacun avec la
             question de lecteur à laquelle il répond (plan_narratif.md, §5) :
               Table S1 : « quelles chaînes de coordonnées portent des verdicts contradictoires, et
                           pour quel médicament ? » (les dix-neuf couples du catalogue C1) ;
               Table S2 : « la décomposition des blocs compte-t-elle hors fluoroquinolones ? »
                           (gain de la décomposition, en points de sensibilité, par médicament,
                           sous C0 et sous C1).
             Le script écrit deux CSV et deux fragments LaTeX inclus par `article/supplementary.tex`.
Entrées    : résultats/model/p77_conflits_C1.tsv (séparateur `|`), p77_C0.tsv, p77_C1.tsv
             (cible `reference`, état `mnv`, colonne `gain_sens` et son intervalle).
Sorties    : article/supplementary_materials/TableS1_contradictory_pairs_C1.csv,
             article/supplementary_materials/TableS2_block_decomposition_gain_by_medicine.csv,
             article/supplementary_generated/tableS1.tex et tableS2.tex (tabulars inclus par
             `article/supplementary.tex`).
Réutilisable : non (propre à A1) ; le patron « TSV -> CSV + fragment LaTeX » se recopie.
Projet     : Resistance_antibio (P7.12)
Date       : 2026-11-21
Usage      : python3 analyses/p712_supplementaire.py
"""
from pathlib import Path

import pandas as pd

RACINE = Path(__file__).resolve().parent.parent
M = RACINE / "résultats" / "model"
SUPP = RACINE / "article" / "supplementary_materials"
GEN = RACINE / "article" / "supplementary_generated"

NOMS = {"isoniazid": "Isoniazid", "rifampicin": "Rifampicin", "ethambutol": "Ethambutol",
        "pyrazinamide": "Pyrazinamide", "moxifloxacin": "Moxifloxacin",
        "levofloxacin": "Levofloxacin", "streptomycin": "Streptomycin",
        "ethionamide": "Ethionamide", "amikacin": "Amikacin", "kanamycin": "Kanamycin",
        "capreomycin": "Capreomycin"}
ORDRE = list(NOMS)


def tex(s):
    return str(s).replace("_", r"\_")


def table_s1():
    c = pd.read_csv(M / "p77_conflits_C1.tsv", sep="|")
    c = c.rename(columns={"spdi": "coordinate_string", "drug": "medicine"})
    c["medicine"] = c["medicine"].map(lambda d: NOMS.get(d, d))
    ordre = {NOMS[d]: i for i, d in enumerate(ORDRE)}
    c = c.sort_values(["medicine", "coordinate_string"],
                      key=lambda s: s.map(ordre) if s.name == "medicine" else s).reset_index(drop=True)
    assert len(c) == 19, f"19 couples attendus, {len(c)} lus"
    return c


def table_s2():
    d = None
    for cat, f in [("C0", "p77_C0.tsv"), ("C1", "p77_C1.tsv")]:
        t = pd.read_csv(M / f, sep="\t")
        t = t[(t["cible"] == "reference") & (t["etat"] == "mnv")].set_index("drug")
        t = (100 * t[["gain_sens", "gain_lo", "gain_hi"]]).round(2)
        t.columns = [f"gain_{cat}", f"lo_{cat}", f"hi_{cat}"]
        d = t if d is None else d.join(t)
    d = d.loc[ORDRE].reset_index()
    d["drug"] = d["drug"].map(NOMS)
    return d.rename(columns={"drug": "medicine", "gain_C0": "gain_points_C0", "lo_C0": "ci_low_C0",
                             "hi_C0": "ci_high_C0", "gain_C1": "gain_points_C1",
                             "lo_C1": "ci_low_C1", "hi_C1": "ci_high_C1"})


def main():
    SUPP.mkdir(parents=True, exist_ok=True)
    GEN.mkdir(parents=True, exist_ok=True)
    s1, s2 = table_s1(), table_s2()
    s1.to_csv(SUPP / "TableS1_contradictory_pairs_C1.csv", index=False)
    s2.to_csv(SUPP / "TableS2_block_decomposition_gain_by_medicine.csv", index=False)

    ent = r"% Fragment généré par analyses/p712_supplementaire.py : ne pas éditer à la main."
    l1 = [ent, r"\begin{tabular}{ll}", r"\toprule", r"Coordinate string & Medicine \\", r"\midrule"]
    for chaine, medicament in zip(s1["coordinate_string"], s1["medicine"]):
        l1.append(rf"\texttt{{{tex(chaine)}}} & {medicament} \\")
    l1 += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "tableS1.tex").write_text("\n".join(l1) + "\n", encoding="utf-8")

    l2 = [ent, r"\begin{tabular}{lrrrrrr}", r"\toprule",
          r" & \multicolumn{3}{c}{Ingested catalogue C0} & \multicolumn{3}{c}{Catalogue C1} \\",
          r"Medicine & Gain & 95\,\% low & 95\,\% high & Gain & 95\,\% low & 95\,\% high \\", r"\midrule"]
    for ligne in s2.itertuples(index=False):
        m, g0, a0, b0, g1, a1, b1 = ligne
        l2.append(f"{m} & {g0:.2f} & {a0:.2f} & {b0:.2f} & {g1:.2f} & {a1:.2f} & {b1:.2f} \\\\")
    l2 += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "tableS2.tex").write_text("\n".join(l2) + "\n", encoding="utf-8")
    print(f"Table S1 : {len(s1)} couples ; Table S2 : {len(s2)} médicaments")


if __name__ == "__main__":
    main()
