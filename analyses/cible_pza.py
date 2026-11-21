#!/usr/bin/env python3
"""Cible supervisée PZA : règle d'exclusion des phénotypes suspects (P5.3).

Objet      : la cible pyrazinamide (PZA) de `rows.tsv` est contaminée par deux blocs que
             le catalogue explique mal et que XGBoost apprend à contresens (cahier
             2026-10-10 13:20, P5.1). (1) *M. bovis* / *M. bovis* BCG : résistance
             INTRINSÈQUE d'espèce (pncA H57D), 831 souches à 85 % R, que GroupKFold met
             dans un pli à part ; (2) le BioProject américain PRJNA1237251, seule étude
             à fournir des phénotypes `TBannotator_biosample` à la matrice (5 893
             *M. tuberculosis*) : 91,7 % de ses R n'ont aucun variant *pncA*, contre
             10,5 % ailleurs. La règle est posée PAR ÉTUDE et PAR ESPÈCE, jamais par
             « R sans variant pncA » (filtrer sur le génotype pour juger le phénotype
             est circulaire et gonflerait l'AUC).
Entrées    : résultats/model/rows.tsv (phase4) ; catalogue/exclusions_pza_bioprojects.tsv
             (strain_id -> BioProject, produit une fois depuis ENA par `--build-registry`).
Sorties    : colonnes `pyrazinamide_brut` (valeur d'origine) et `pyrazinamide_motif`
             (raison d'exclusion ou vide) ajoutées à rows.tsv ; `pyrazinamide` vidée
             pour les souches exclues. Tous les consommateurs de rows.tsv filtrent sur
             `isin(["R","S"])`, donc l'exclusion se propage sans toucher aux phases.
Réutilisable : `appliquer(rows, registre)` pour toute autre drogue au même motif
             (étude à phénotype aberrant, résistance intrinsèque d'espèce).
Projet     : Resistance_antibio (P5.3)
Date       : 2026-10-10
"""
import argparse
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

MODELDIR = paths.RESULTATS / "model"
REGISTRE = paths.CATALOGUE_DIR / "exclusions_pza_bioprojects.tsv"

# Études dont les phénotypes PZA sont écartés (motif mesuré, cf. cahier 2026-10-10).
BIOPROJECTS_PZA_ECARTES = {"PRJNA1237251"}
LIGNEES_BOVIS = {"BOV", "BOV_AFRI"}
DRUG = "pyrazinamide"


def construire_registre(strain_ids, out=REGISTRE):
    """Interroge ENA (POST par lots de 150) : strain_id (run) -> study_accession."""
    ids = list(strain_ids)
    parts = []
    for i in range(0, len(ids), 150):
        q = " OR ".join(f"run_accession={x}" for x in ids[i:i + 150])
        data = urllib.parse.urlencode({
            "result": "read_run", "query": q, "format": "tsv", "limit": 0,
            "fields": "run_accession,study_accession,country,scientific_name"}).encode()
        for _ in range(3):
            try:
                txt = urllib.request.urlopen(
                    "https://www.ebi.ac.uk/ena/portal/api/search", data, timeout=60
                ).read().decode()
                parts.append(pd.read_csv(io.StringIO(txt), sep="\t"))
                break
            except Exception:
                time.sleep(3)
    reg = (pd.concat(parts).rename(columns={"run_accession": "strain_id",
                                            "study_accession": "bioproject"}))
    reg.to_csv(out, sep="\t", index=False)
    return reg


def appliquer(rows, registre):
    """Ajoute `pyrazinamide_brut` / `pyrazinamide_motif` et vide `pyrazinamide` des exclus."""
    rows = rows.copy()
    if f"{DRUG}_brut" not in rows.columns:          # idempotent : on repart de la valeur brute
        rows[f"{DRUG}_brut"] = rows[DRUG]
    rows[DRUG] = rows[f"{DRUG}_brut"]
    bp = rows["strain_id"].map(registre.set_index("strain_id")["bioproject"])
    bovis = rows["lineage_level_1"].isin(LIGNEES_BOVIS)
    etude = bp.isin(BIOPROJECTS_PZA_ECARTES)
    motif = pd.Series("", index=rows.index)
    motif[etude] = "bioproject_" + bp[etude].astype(str)
    motif[bovis] = "bovis_resistance_intrinseque"      # l'espèce prime sur l'étude
    rows[f"{DRUG}_motif"] = motif
    rows.loc[motif != "", DRUG] = ""
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-registry", action="store_true",
                    help="(re)construit le registre strain_id -> BioProject depuis ENA")
    ap.add_argument("--apply", action="store_true",
                    help="applique la règle à résultats/model/rows.tsv (idempotent)")
    args = ap.parse_args()

    rows = pd.read_csv(MODELDIR / "rows.tsv", sep="\t", keep_default_na=False,
                       low_memory=False)
    if args.build_registry:
        cible = rows.loc[rows[f"{DRUG}_brut" if f"{DRUG}_brut" in rows else DRUG]
                         .isin(["R", "S"]), "strain_id"]
        reg = construire_registre(cible)
        print(f"registre : {len(reg)} souches -> {REGISTRE}")
    if args.apply:
        reg = pd.read_csv(REGISTRE, sep="\t")
        out = appliquer(rows, reg)
        out.to_csv(MODELDIR / "rows.tsv", sep="\t", index=False)
        brut = out[f"{DRUG}_brut"].isin(["R", "S"]).sum()
        reste = out[DRUG].isin(["R", "S"]).sum()
        print(f"PZA : {brut} -> {reste} souches phénotypées "
              f"({out[f'{DRUG}_motif'].value_counts().drop('', errors='ignore').to_dict()})")


if __name__ == "__main__":
    main()
