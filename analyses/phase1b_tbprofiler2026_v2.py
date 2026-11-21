#!/usr/bin/env python3
"""Phase 1b — Catalogue v2 : source tb-profiler remplacée par tbdb 2026, en parallèle.

N'écrit QUE dans `catalogue_v2/` : le catalogue courant (`catalogue/`) et les chiffres du
manuscrit restent intacts. Le catalogue v2 est le catalogue courant dont on retire les lignes
`tbprofiler` et `tbprofiler_AAprop`, remplacées par les lignes `type == drug_resistance` du
`mutations.csv` de tbdb (jodyphelan/tbdb, commit épinglé dans `data/sources/tbprofiler_2026/`),
puis la propagation acide aminé est rejouée. Le mapping des colonnes, `spdis_for` et
`_call_from_grade` sont ceux de `phase1_build_catalogue.py` (parité voulue : seule la donnée change).

Schéma tbdb 2026 (Gene, Mutation, type, drug, original_mutation, confidence, source, comment)
contre 2022 (Gene, Mutation, Drug, Confers, Interaction, Literature, WHO Confidence).
Les lignes `type == who_confidence` (46 856) dupliquent la source OMS déjà ingérée : exclues.

Confiance vide (644 lignes tbdb « Mutation from literature », 653 lignes développées) : `call=uncertain`
par décision (P4.1, 2026-10-10), comme dans la version 2022. Entrées non gradées par l'OMS, effet d'un
reclassement en R quasi nul ; motifs et chiffres dans `catalogue/SOURCES.md`.

AA1 (2026-11-11) : la colonne `regle` (épistasie Rv0678/mmpL5, eis ; `regles_interpretation.py`) est ajoutée en fin de
schéma, après la propagation. `--bascule` écrit aussi le résultat dans `catalogue/catalogue_consolide.tsv` (sauvegarde
préalable dans `archives/sauvegardes/`) : c'est le geste que fait `update_database.py` après `phase1`, qui reconstruit le
catalogue v1 (tb-profiler 2022) et l'écraserait sans lui.

Usage : .venv-ingest/bin/python analyses/phase1b_tbprofiler2026_v2.py [--bascule]
"""
import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
from vocab import canonical_drug
from phase1_build_catalogue import COLUMNS, _call_from_grade, spdis_for
from regles_interpretation import annoter_catalogue

COLUMNS_V2 = COLUMNS + ["regle"]

SRC = paths.SOURCES / "tbprofiler_2026" / "mutations.csv"
OUT_DIR = paths.ROOT / "catalogue_v2"
OUT = OUT_DIR / "catalogue_consolide_v2.tsv"


def ingest_tbdb2026() -> pd.DataFrame:
    t = pd.read_csv(SRC, dtype=str, keep_default_na=False)
    t = t[t["type"] == "drug_resistance"]
    rows = []
    for _, r in t.iterrows():
        canon, code = canonical_drug(r["drug"])
        variant = f"{r['Gene']}_{r['Mutation']}"
        for spdi in spdis_for(variant):
            rows.append({
                "spdi": spdi, "gene": r["Gene"], "variant": variant,
                "drug": canon, "drug_code": code,
                "call": _call_from_grade(r["confidence"]),
                "source": "tbprofiler",
                "source_grade": r["confidence"],
                "n_R": "", "n_S": "", "metric": "",
                "notes": f"tbdb2026;src={r['source']}"[:120],
            })
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[tbdb2026] {len(t)} lignes drug_resistance -> {len(df)} assertions, "
          f"{df['drug'].nunique()} drogues, {(df['spdi'] != '').sum()} avec SPDI résolu")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bascule", action="store_true",
                    help="écrit aussi catalogue/catalogue_consolide.tsv (sauvegarde préalable)")
    # P7.7 : test falsifiant du recadrage de A1 (blocs OMS émis par spdis_for) ; les défauts
    # reproduisent le comportement historique à l'identique.
    ap.add_argument("--blocs", action="store_true",
                    help="P7.7 : émettre aussi les représentations en bloc de l'OMS (cf. phase1_build_catalogue.EMETTRE_BLOCS)")
    ap.add_argument("--entree", help="catalogue de départ (défaut : le catalogue courant)")
    ap.add_argument("--sortie", help="fichier de sortie (défaut : catalogue_v2/catalogue_consolide_v2.tsv)")
    a = ap.parse_args()
    if a.bascule and (a.entree or a.sortie or a.blocs):
        ap.error("--bascule n'est pas combinable avec --entree, --sortie ni --blocs (test P7.7 isolé)")
    import phase1_build_catalogue as p1
    p1.EMETTRE_BLOCS = a.blocs
    OUT = Path(a.sortie) if a.sortie else globals()["OUT"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(exist_ok=True)
    cur = pd.read_csv(a.entree or paths.CATALOGUE_TSV, sep="\t", dtype=str, keep_default_na=False)
    base = cur[~cur["source"].isin(["tbprofiler", "tbprofiler_AAprop"])]
    v2 = pd.concat([base, ingest_tbdb2026()], ignore_index=True)[COLUMNS]

    from aa_propagation import propagate_catalogue
    n0 = len(v2)
    v2 = propagate_catalogue(v2)[COLUMNS]
    print(f"[AA-prop] +{len(v2) - n0}")
    v2 = annoter_catalogue(v2)[COLUMNS_V2]
    print("[regles] lignes par règle :", v2["regle"].replace("", "(aucune)").value_counts().to_dict())
    v2 = v2.sort_values(["drug", "source", "gene"], kind="stable").reset_index(drop=True)
    v2.to_csv(OUT, sep="\t", index=False)
    if a.bascule:
        sys.path.insert(0, str(Path.home() / ".claude" / "skills" / "pistes"))
        from archives_projet import sauvegarder
        sauv = sauvegarder(paths.CATALOGUE_TSV, "avant_AA1_bascule", racine=paths.ROOT) \
            if paths.CATALOGUE_TSV.exists() else None
        shutil.copyfile(OUT, paths.CATALOGUE_TSV)
        print(f"[bascule] {paths.CATALOGUE_TSV} réécrit (sauvegarde : {sauv})")

    print(f"\n=== v2 -> {OUT} : {len(v2)} assertions (courant : {len(cur)}) ===")
    for name, d in (("courant", cur), ("v2", v2)):
        tp = d[d["source"].str.startswith("tbprofiler")]
        print(f"\n{name} tbprofiler : {len(tp)} assertions, "
              f"{(tp['spdi'] != '').sum()} avec SPDI")
        print(pd.crosstab(tp["drug"], tp["call"]).to_string())


if __name__ == "__main__":
    main()
