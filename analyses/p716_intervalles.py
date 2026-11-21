#!/usr/bin/env python3
"""P7.16 — Intervalles appariés du gain de la décomposition et du RAPPORT entre les deux ingestions.

Objet      : la review du 2026-11-21 (points 3, 4 et 5 obligatoires) demande (a) les intervalles de
             Clopper-Pearson de la sensibilité et de la spécificité, déjà calculés par
             `p77_baseline_blocs.py` mais jamais imprimés, et (b) un intervalle sur le rapport
             (gain C0 - gain C1) / gain C0, c'est-à-dire la part du gain apparent qui suit de la
             règle d'ingestion (« 77 et 76 % » du résumé). Ce rapport est une fonction de DEUX
             grandeurs rééchantillonnées : son intervalle exige les tirages APPARIÉS, tirage par
             tirage, sous les deux catalogues. `p77_baseline_blocs.py` ne garde que les percentiles
             de chaque gain, pas les tirages, d'où ce script.
             L'appariement est licite parce que les deux bras partagent la matrice, les étiquettes
             et les grappes de sous-lignée, et que `pick` est tiré d'une graine fixe (rng(0)).
Entrées    : --modeldir (X.npz, features.txt, rows.tsv de la matrice à blocs bruts),
             --catalogue-c0 et --catalogue-c1 (TSV au schéma de `catalogue_consolide.tsv`),
             --drugs (défaut : moxifloxacin levofloxacin), --tirages (défaut 2000, comme p77).
Sorties    : --sortie (TSV : une ligne par médicament, gains ponctuels et IC des deux bras, rapport
             et son IC) et `<sortie>_grappes.tsv`, les COMPTES PAR GRAPPE (tp, fn, tn, fp par code
             de sous-lignée, par bras et par état). Ce second fichier est le seul dont dépend tout
             bootstrap ultérieur : il rend tout rejeu LOCAL, sans la matrice ni le mésocentre.
Réutilisable : oui (tout couple de catalogues, toute matrice de `phase4`, tout médicament de p74).
Projet     : Resistance_antibio (P7.16)
Date       : 2026-11-21
Usage      : python3 analyses/p716_intervalles.py --modeldir résultats/model \
               --catalogue-c0 catalogue/p77_v2_defaut.tsv --catalogue-c1 catalogue/p77_v2_blocs.tsv \
               --sortie résultats/model/p716_intervalles.tsv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import variant_match as vm
import cible_reference as cr
import p74_angles_morts as p74
import p77_baseline_blocs as p77

CIBLE = "reference"
ETATS = ("complet", "mnv")


def comptes(X, feats, sig, cat, rows, drug, cible=CIBLE):
    """Comptes par grappe de sous-lignée, pour l'état complet et l'état privé de la décomposition.

    Rend (noms de grappe, {etat: tableau (ncl, 4) de tp, fn, tn, fp}, {etat: sens ponctuelle})."""
    det = vm.determinant_snp_set(cat, drug)
    cols = {"complet": p77.colonnes(feats, sig, det),
            "mnv": p77.colonnes(feats, sig, det, mnv_decompose=False)}
    porte = {k: p74._porte(X, v) for k, v in cols.items()}
    lab = cr.etiquettes_cible(rows, drug, cible)
    ym = lab.isin(["R", "S"]).values
    y = (lab[ym] == "R").astype(int).values
    cl_noms, cl_idx = np.unique(rows["lineage_code"].values, return_inverse=True)
    g = cl_idx[ym]
    cnt = {k: p74.comptes_par_grappe(v[ym], y, g, len(cl_noms)) for k, v in porte.items()}
    sens = {}
    for k, v in cnt.items():
        se, _ = p74.sens_spec_depuis(v.sum(axis=0))
        sens[k] = float(se[0])
    return cl_noms, cnt, sens


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modeldir", required=True)
    ap.add_argument("--catalogue-c0", required=True)
    ap.add_argument("--catalogue-c1", required=True)
    ap.add_argument("--sortie", required=True)
    ap.add_argument("--drugs", nargs="+", default=["moxifloxacin", "levofloxacin"])
    ap.add_argument("--tirages", type=int, default=2000)
    a = ap.parse_args()

    md = Path(a.modeldir)
    X = sparse.load_npz(md / "X.npz").tocsr()
    feats = (md / "features.txt").read_text().splitlines()
    rows = pd.read_csv(md / "rows.tsv", sep="\t", keep_default_na=False)
    assert X.shape[0] == len(rows) and X.shape[1] == len(feats), "X.npz et rows.tsv désalignés"
    print(f"[matrice] {X.shape[0]} souches x {X.shape[1]} colonnes")

    sig = {}
    for s in feats:
        _, ref, alt = vm.parse_spdi(s)
        sig[s] = ("m", frozenset(vm.decompose(s))) if len(ref) == len(alt) and len(ref) > 1 else ("x", s)

    cats = {"C0": pd.read_csv(a.catalogue_c0, sep="\t", dtype=str).fillna(""),
            "C1": pd.read_csv(a.catalogue_c1, sep="\t", dtype=str).fillna("")}
    for k, c in cats.items():
        print(f"[catalogue] {k} : {len(c)} lignes")

    ncl = len(np.unique(rows["lineage_code"].values))
    pick = np.random.default_rng(0).integers(0, ncl, (a.tirages, ncl))

    lignes, grappes = [], []
    for drug in a.drugs:
        res = {}
        for bras, cat in cats.items():
            cl_noms, cnt, sens = comptes(X, feats, sig, cat, rows, drug)
            res[bras] = (cnt, sens)
            for etat in ETATS:
                for i, nom in enumerate(cl_noms):
                    tp, fn, tn, fp = (int(x) for x in cnt[etat][i])
                    if tp + fn + tn + fp == 0:
                        continue
                    grappes.append({"drug": drug, "catalogue": bras, "etat": etat,
                                    "lineage_code": nom, "tp": tp, "fn": fn, "tn": tn, "fp": fp})
        # gains ponctuels et tirages appariés : les deux bras partagent `pick`
        delta = {}
        for bras in ("C0", "C1"):
            cnt, sens = res[bras]
            bc, _ = p74.sens_spec_depuis(cnt["complet"][pick].sum(axis=1))
            bm, _ = p74.sens_spec_depuis(cnt["mnv"][pick].sum(axis=1))
            delta[bras] = bc - bm
        g0, g1 = delta["C0"], delta["C1"]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(g0 > 0, (g0 - g1) / g0, np.nan)
        tot = res["C0"][0]["complet"].sum(axis=0)
        n_R, n_S = int(tot[0] + tot[1]), int(tot[2] + tot[3])
        pt0 = res["C0"][1]["complet"] - res["C0"][1]["mnv"]
        pt1 = res["C1"][1]["complet"] - res["C1"][1]["mnv"]
        lignes.append({
            "drug": drug, "n_R": n_R, "n_S": n_S,
            "gain_C0": 100 * pt0, "gain_C0_lo": 100 * float(np.nanpercentile(g0, 2.5)),
            "gain_C0_hi": 100 * float(np.nanpercentile(g0, 97.5)),
            "gain_C1": 100 * pt1, "gain_C1_lo": 100 * float(np.nanpercentile(g1, 2.5)),
            "gain_C1_hi": 100 * float(np.nanpercentile(g1, 97.5)),
            "part_ingestion": 100 * (pt0 - pt1) / pt0,
            "part_lo": 100 * float(np.nanpercentile(ratio, 2.5)),
            "part_hi": 100 * float(np.nanpercentile(ratio, 97.5)),
            "n_tirages_valides": int(np.isfinite(ratio).sum()), "tirages": a.tirages,
        })
        print(f"[{drug}] gain C0 {100*pt0:.2f} ({100*np.nanpercentile(g0,2.5):.2f} a "
              f"{100*np.nanpercentile(g0,97.5):.2f}) ; gain C1 {100*pt1:.2f} "
              f"({100*np.nanpercentile(g1,2.5):.2f} a {100*np.nanpercentile(g1,97.5):.2f}) ; "
              f"part de l'ingestion {100*(pt0-pt1)/pt0:.1f} % "
              f"({100*np.nanpercentile(ratio,2.5):.1f} a {100*np.nanpercentile(ratio,97.5):.1f})")

    df = pd.DataFrame(lignes)
    df.to_csv(a.sortie, sep="\t", index=False, float_format="%.4f")
    gp = Path(a.sortie).with_name(Path(a.sortie).stem + "_grappes.tsv")
    pd.DataFrame(grappes).to_csv(gp, sep="\t", index=False)
    print(f"\n-> {a.sortie} ; {gp} ({len(grappes)} lignes de comptes par grappe)")


if __name__ == "__main__":
    main()
