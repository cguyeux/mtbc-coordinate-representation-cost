#!/usr/bin/env python3
"""P7.7 — Comparaison du catalogue vivant (C0) et du catalogue à blocs OMS émis (C1).

Objet      : (1) CONTRÔLE : la chaîne phase1 + phase1b rejouée avec ses défauts redonne-t-elle le
             catalogue v2 vivant ? (égalité d'ensembles de lignes, tous champs) ; (2) mesurer ce que
             change l'émission des blocs : lignes ajoutées par source et par médicament, types de
             représentation (SNP, bloc MNV, indel), SPDI à verdicts contradictoires pour un même
             médicament (le même bloc `ACAG>GCAC` est listé par l'OMS sous gyrA Asp94Gly ET sous
             gyrA Ser95Thr : aparté de P7.9), et part des représentations OMS retrouvées.
Entrées    : --c0 (catalogue vivant), --controle (v2 rejoué avec les défauts), --c1 (v2 à blocs).
Sorties    : un TSV long `--sortie` (indicateur, valeur, détail) et un résumé imprimé.
Réutilisable : oui (toute paire de catalogues au schéma COLUMNS + regle).
Projet     : Resistance_antibio (P7.7)
Date       : 2026-11-16
Usage      : python3 analyses/p77_comparer_catalogues.py --c0 catalogue/catalogue_consolide.tsv \
               --controle catalogue/p77_v2_defaut.tsv --c1 catalogue/p77_v2_blocs.tsv \
               --sortie results/model/p77_catalogues.tsv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import variant_match as vm


def charger(f):
    return pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)


def type_spdi(s):
    if not s:
        return "sans_spdi"
    _, ref, alt = vm.parse_spdi(s)
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    return "mnv" if len(ref) == len(alt) else "indel"


def jeu_lignes(df):
    return set(map(tuple, df.itertuples(index=False, name=None)))


def conflits(df):
    """SPDI portant à la fois un verdict R-associé et un verdict S-associé pour un même médicament."""
    d = df[df["spdi"] != ""]
    g = d.groupby(["spdi", "drug"])["call"].agg(lambda c: set(c))
    mauvais = g[g.map(lambda s: "R-associated" in s and "S-associated" in s)]
    return mauvais


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c0", required=True)
    ap.add_argument("--controle", required=True)
    ap.add_argument("--c1", required=True)
    ap.add_argument("--sortie", required=True)
    a = ap.parse_args()

    c0, ctl, c1 = charger(a.c0), charger(a.controle), charger(a.c1)
    out = []

    def noter(ind, val, det=""):
        out.append({"indicateur": ind, "valeur": val, "detail": det})
        print(f"{ind:<46} {val}   {det}")

    # 1. contrôle : la chaîne rejouée redonne-t-elle le catalogue vivant ?
    j0, jc = jeu_lignes(c0), jeu_lignes(ctl)
    noter("lignes_C0_vivant", len(c0))
    noter("lignes_controle_rejoue", len(ctl))
    noter("lignes_controle_absentes_de_C0", len(jc - j0))
    noter("lignes_C0_absentes_du_controle", len(j0 - jc))
    noter("controle_identique_en_ensemble", int(j0 == jc))

    # 2. effet de l'émission des blocs
    noter("lignes_C1_blocs", len(c1), f"x{len(c1) / len(c0):.2f} de C0")
    j1 = jeu_lignes(c1)
    noter("lignes_C0_absentes_de_C1", len(j0 - j1), "doit valoir 0 : C1 ⊇ C0")
    noter("lignes_ajoutees_par_C1", len(j1 - j0))
    for nom, df in (("C0", c0), ("C1", c1)):
        t = df["spdi"].map(type_spdi).value_counts()
        for k, v in t.items():
            noter(f"types_spdi_{nom}", int(v), k)
        r = df[df["call"] == "R-associated"]["spdi"].map(type_spdi).value_counts()
        for k, v in r.items():
            noter(f"types_spdi_R_associe_{nom}", int(v), k)
        noter(f"spdi_distincts_{nom}", df.loc[df["spdi"] != "", "spdi"].nunique())
        noter(f"spdi_distincts_R_associes_{nom}",
              df.loc[(df["spdi"] != "") & (df["call"] == "R-associated"), "spdi"].nunique())

    nouvelle = [t not in j0 for t in c1.itertuples(index=False, name=None)]
    ajoutees = c1[nouvelle]
    for (src, drug), n in ajoutees.groupby(["source", "drug"]).size().sort_values(ascending=False).head(25).items():
        noter("ajoutees_par_source_et_medicament", int(n), f"{src} | {drug}")
    for k, v in ajoutees["call"].value_counts().items():
        noter("ajoutees_par_verdict", int(v), k)
    for k, v in ajoutees["spdi"].map(type_spdi).value_counts().items():
        noter("ajoutees_par_type", int(v), k)

    # 3. verdicts contradictoires pour un même (SPDI, médicament)
    for nom, df in (("C0", c0), ("C1", c1)):
        cf = conflits(df)
        noter(f"conflits_R_S_meme_spdi_meme_medicament_{nom}", len(cf))
        if nom == "C1":
            pd.Series(cf.index.map(lambda t: f"{t[0]}|{t[1]}")).to_csv(
                Path(a.sortie).with_name("p77_conflits_C1.tsv"), sep="\t", index=False, header=["spdi|drug"])
    # le cas d'école de P7.9 : le bloc ACAG>GCAC de gyrA (0-based 7581)
    for nom, df in (("C0", c0), ("C1", c1)):
        d = df[df["spdi"] == "NC_000962.3:7581:ACAG:GCAC"]
        for _, r in d.iterrows():
            noter(f"bloc_gyrA_7581_ACAG_GCAC_{nom}", r["call"], f"{r['variant']} | {r['drug']} | {r['source']} | {r['source_grade']}")

    pd.DataFrame(out).to_csv(a.sortie, sep="\t", index=False)
    print(f"-> {a.sortie}")


if __name__ == "__main__":
    main()
