#!/usr/bin/env python3
"""P7.16 — Table supplémentaire S3 d'A1 : intervalles de Clopper-Pearson de la sensibilité et de la
spécificité, promis par les méthodes et imprimés nulle part (review du 2026-11-21, point 3).

Objet      : les méthodes d'A1 annoncent « Sensitivity and specificity are given per medicine on the
             measured phenotypes, with Clopper-Pearson 95 % intervals » ; les tableaux 1 et 2
             donnaient les points nus. Ces intervalles sont DÉJÀ calculés par
             `p77_baseline_blocs.py` (colonnes sens_lo, sens_hi, spec_lo, spec_hi, via
             `p74.clopper_pearson`) : ce script les met en table plutôt que de retirer la promesse.
             Les effectifs R vont de 991 (capréomycine) à 13 256 (isoniazide), donc la précision des
             onze sensibilités n'est pas comparable, et c'est ce que la table montre.
Entrées    : résultats/model/p77_C0.tsv et p77_C1.tsv (cible `reference`, états `complet` et `mnv`).
Sorties    : article/supplementary_materials/TableS3_clopper_pearson_intervals.csv et
             article/supplementary_generated/tableS3.tex (tabular inclus par `supplementary.tex`).
Réutilisable : non (propre à A1) ; le patron « TSV -> CSV + fragment LaTeX » vient de p712.
Projet     : Resistance_antibio (P7.16)
Date       : 2026-11-21
Usage      : python3 analyses/p716_supplementaire.py
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
FQ = ("moxifloxacin", "levofloxacin")


def lire(fichier, etat):
    t = pd.read_csv(M / fichier, sep="\t")
    return t[(t["cible"] == "reference") & (t["etat"] == etat)].set_index("drug")


def table_s3():
    """Une ligne par médicament et par ingestion, état complet ; plus l'appariement exact seul
    pour les deux fluoroquinolones, qui sont les lignes du tableau 1 du corps."""
    lignes = []
    for etat, nom_etat in (("complet", "complete"), ("mnv", "exact only")):
        for cat in ("C0", "C1"):
            t = lire(f"p77_{cat}.tsv", etat)
            for drug in (ORDRE if etat == "complet" else list(FQ)):
                r = t.loc[drug]
                lignes.append({
                    "medicine": NOMS[drug], "matching": nom_etat, "ingestion": cat,
                    "n_R": int(r["n_R"]), "n_S": int(r["n_S"]),
                    "sensitivity": round(float(r["sens"]), 4),
                    "sens_ci_low": round(float(r["sens_lo"]), 4),
                    "sens_ci_high": round(float(r["sens_hi"]), 4),
                    "specificity": round(float(r["spec"]), 4),
                    "spec_ci_low": round(float(r["spec_lo"]), 4),
                    "spec_ci_high": round(float(r["spec_hi"]), 4),
                })
    d = pd.DataFrame(lignes)
    ordre = {NOMS[k]: i for i, k in enumerate(ORDRE)}
    d["_o"] = d["medicine"].map(ordre)
    d["_m"] = d["matching"].map({"complete": 0, "exact only": 1})
    return d.sort_values(["_m", "_o", "ingestion"]).drop(columns=["_o", "_m"]).reset_index(drop=True)


def main():
    SUPP.mkdir(parents=True, exist_ok=True)
    GEN.mkdir(parents=True, exist_ok=True)
    s3 = table_s3()
    assert len(s3) == 2 * len(ORDRE) + 2 * len(FQ), f"{len(s3)} lignes, attendu {2*len(ORDRE)+2*len(FQ)}"
    s3.to_csv(SUPP / "TableS3_clopper_pearson_intervals.csv", index=False)

    ent = r"% Fragment généré par analyses/p716_supplementaire.py : ne pas éditer à la main."
    l3 = [ent, r"\begin{tabular}{llrrcc}", r"\toprule",
          r"Medicine & Ingestion & $n_R$ & $n_S$ & Sensitivity (95\,\% CI) & Specificity (95\,\% CI) \\",
          r"\midrule"]
    bloc = None
    for r in s3.itertuples(index=False):
        if bloc is not None and r.matching != bloc:
            l3 += [r"\midrule",
                   r"\multicolumn{6}{l}{\textit{Exact matching only (no decomposition of observed"
                   r" blocks)}} \\", r"\midrule"]
        bloc = r.matching
        l3.append(f"{r.medicine} & {r.ingestion} & {r.n_R} & {r.n_S} & "
                  f"{r.sensitivity:.3f} ({r.sens_ci_low:.3f} to {r.sens_ci_high:.3f}) & "
                  f"{r.specificity:.3f} ({r.spec_ci_low:.3f} to {r.spec_ci_high:.3f}) \\\\")
    l3 += [r"\bottomrule", r"\end{tabular}"]
    (GEN / "tableS3.tex").write_text("\n".join(l3) + "\n", encoding="utf-8")
    print(f"Table S3 : {len(s3)} lignes ({len(ORDRE)} médicaments x 2 ingestions, "
          f"plus {len(FQ)} fluoroquinolones x 2 en appariement exact)")
    ic = s3[s3.matching == "complete"]
    large = (ic["sens_ci_high"] - ic["sens_ci_low"]).idxmax()
    etroit = (ic["sens_ci_high"] - ic["sens_ci_low"]).idxmin()
    for nom, i in (("le plus large", large), ("le plus étroit", etroit)):
        r = ic.loc[i]
        print(f"  IC de sensibilité {nom} : {r['medicine']} ({r['ingestion']}), n_R={r['n_R']}, "
              f"{r['sensitivity']:.3f} ({r['sens_ci_low']:.3f} a {r['sens_ci_high']:.3f})")


if __name__ == "__main__":
    main()
