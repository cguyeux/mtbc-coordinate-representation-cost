#!/usr/bin/env python3
"""P7.7 — Test falsifiant du recadrage de A1 : baseline et gain des quatre angles morts sous un catalogue qui
émet les blocs OMS, face au catalogue vivant, sur UNE MÊME matrice.

Objet      : le verdict RECADRER (P7.3) affirme que l'angle mort MNV est un artefact de `spdis_for()`, qui écarte
             les représentations en bloc du fichier de coordonnées OMS dès qu'une forme mono-base existe. Prédiction
             falsifiable : si le catalogue émet ces blocs, l'appariement EXACT suffit pour les haplotypes que l'OMS
             énumère, et le gain de la décomposition en SNP (+57,2 MXF, +49,7 LEV, P7.6) s'effondre ; il ne
             subsiste que pour les haplotypes absents du fichier OMS (D94N, D94Y, D94H sur fond S95T, P7.9).
             Rejoue donc `p74_angles_morts.py` par retrait pour UN catalogue donné et rend, en plus,
             le RÉSIDUEL : les colonnes de la matrice qui ne portent un déterminant QUE par décomposition.
Différence avec p74 : l'appariement « complet » d'un bloc MNV est exact OU décomposé (`variant_match.
             carrier_columns`), alors que p74 ne regardait que la décomposition pour un bloc ; sans effet pour le
             catalogue vivant (à vérifier : contrôle (matrice P7.6, C0) contre `p74_angles_morts_blocs_P15.tsv`).
Entrées    : --modeldir (X.npz, features.txt, rows.tsv), --catalogue (TSV), --etiquette (nom de l'état).
Sorties    : --sortie (TSV long : cible, drug, etat, n_R, n_S, sens, spec, IC, gain_sens, gain_lo, gain_hi,
             perte_spec, n_col, catalogue) et `<sortie>_residuel.tsv` (colonnes décomposition-seule).
Réutilisable : oui (tout catalogue au schéma de `catalogue_consolide.tsv`, toute matrice de `phase4`).
Projet     : Resistance_antibio (P7.7)
Date       : 2026-11-16
Usage      : python3 analyses/p77_baseline_blocs.py --modeldir results/model --catalogue catalogue/X.tsv \
               --etiquette C1 --sortie results/model/p77_C1.tsv
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

RETRAITS = p74.RETRAITS


def colonnes(feats, sig, det, mnv_decompose=True, indels=True, freq=None):
    """Indices des colonnes qui portent un déterminant de `det`.

    SNP et indel : appariement exact. Bloc MNV (substitution de plus d'une base) : exact, OU, si
    `mnv_decompose`, décomposition en SNP dont l'un est un déterminant SNP. Sans décomposition, seul
    l'appariement exact reste (retrait « mnv »)."""
    mono = {s for s in det if p74.est_snp(s)}
    cols = []
    for j, s in enumerate(feats):
        typ, val = sig[s]
        if typ == "x":
            ok = s in det and (indels or not p74.est_indel(s))
        else:
            ok = (s in det) or (mnv_decompose and bool(val & mono))
        if ok and (freq is None or freq[j] >= p74.MIN_FREQ):
            cols.append(j)
    return cols


def residuel(X, feats, sig, det, rows, drug, cible="reference"):
    """Colonnes MNV qui portent un déterminant par décomposition SEULEMENT (ni exact, ni indel)."""
    mono = {s for s in det if p74.est_snp(s)}
    lab = cr.etiquettes_cible(rows, drug, cible)
    ym = lab.isin(["R", "S"]).values
    y = (lab[ym] == "R").values
    out = []
    for j, s in enumerate(feats):
        typ, val = sig[s]
        if typ != "m" or s in det or not (val & mono):
            continue
        col = X[:, j].toarray().ravel() > 0
        n = int(col.sum())
        out.append({"drug": drug, "spdi": s, "n_souches": n, "n_etiquetees": int((col & ym).sum()),
                    "n_R": int((col[ym] & y).sum()), "n_S": int((col[ym] & ~y).sum()),
                    "composantes_deterministes": ";".join(sorted(val & mono))})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modeldir", required=True)
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--etiquette", required=True)
    ap.add_argument("--sortie", required=True)
    ap.add_argument("--tirages", type=int, default=2000)
    a = ap.parse_args()

    md = Path(a.modeldir)
    X = sparse.load_npz(md / "X.npz").tocsr()
    feats = (md / "features.txt").read_text().splitlines()
    rows = pd.read_csv(md / "rows.tsv", sep="\t", keep_default_na=False)
    assert X.shape[0] == len(rows) and X.shape[1] == len(feats), "X.npz et rows.tsv désalignés"
    cat = pd.read_csv(a.catalogue, sep="\t", dtype=str).fillna("")
    freq = np.asarray(X.sum(axis=0)).ravel()

    sig = {}
    for s in feats:
        _, ref, alt = vm.parse_spdi(s)
        sig[s] = ("m", frozenset(vm.decompose(s))) if len(ref) == len(alt) and len(ref) > 1 else ("x", s)

    cl_noms, cl_idx = np.unique(rows["lineage_code"].values, return_inverse=True)
    ncl = len(cl_noms)
    pick = np.random.default_rng(0).integers(0, ncl, (a.tirages, ncl))

    long, resid = [], []
    for drug in p74.DRUGS:
        cat_d = cat[cat["drug"] == drug]
        det = vm.determinant_snp_set(cat, drug)
        det_1rep = p74.une_representation(cat_d)
        cols = {
            "complet": colonnes(feats, sig, det),
            "mnv": colonnes(feats, sig, det, mnv_decompose=False),
            "indel": colonnes(feats, sig, det, indels=False),
            "synonymes": colonnes(feats, sig, det_1rep),
            "rares": colonnes(feats, sig, det, freq=freq),
            "naif": colonnes(feats, sig, det_1rep, mnv_decompose=False, indels=False, freq=freq),
        }
        porte = {k: p74._porte(X, v) for k, v in cols.items()}
        resid += residuel(X, feats, sig, det, rows, drug)
        for nom in cr.CIBLES:
            lab = cr.etiquettes_cible(rows, drug, nom)
            ym = lab.isin(["R", "S"]).values
            if not ym.any():
                continue
            y = (lab[ym] == "R").astype(int).values
            g = cl_idx[ym]
            cnt = {k: p74.comptes_par_grappe(v[ym], y, g, ncl) for k, v in porte.items()}
            tot = {k: v.sum(axis=0) for k, v in cnt.items()}
            boot = {k: v[pick].sum(axis=1) for k, v in cnt.items()}
            se_c, sp_c = p74.sens_spec_depuis(tot["complet"])
            bse_c, _ = p74.sens_spec_depuis(boot["complet"])
            for k in ["complet"] + RETRAITS:
                se, sp = p74.sens_spec_depuis(tot[k])
                tp, fn, tn, fp = (int(x) for x in tot[k])
                lo, hi = p74.clopper_pearson(tp, tp + fn)
                slo, shi = p74.clopper_pearson(tn, tn + fp)
                ligne = {"catalogue": a.etiquette, "cible": nom, "drug": drug, "etat": k,
                         "n_col": len(cols[k]), "n_R": tp + fn, "n_S": tn + fp,
                         "sens": float(se[0]), "sens_lo": lo, "sens_hi": hi,
                         "spec": float(sp[0]), "spec_lo": slo, "spec_hi": shi}
                if k != "complet":
                    bse, _ = p74.sens_spec_depuis(boot[k])
                    delta = bse_c - bse
                    ligne.update({"gain_sens": float(se_c[0] - se[0]),
                                  "gain_lo": float(np.nanpercentile(delta, 2.5)),
                                  "gain_hi": float(np.nanpercentile(delta, 97.5)),
                                  "perte_spec": float(sp_c[0] - sp[0])})
                long.append(ligne)
    df = pd.DataFrame(long)
    # P7.16 : %.6f et NON %.4f. Les colonnes de sensibilité et de gain sont des FRACTIONS :
    # arrondies à 4 décimales, elles ne gardent que 2 décimales utiles une fois lues en POINTS
    # (0,117472 -> 0,1175 -> 11,75, qui arrondit à 11,8 quand la valeur vaut 11,7). Un chiffre
    # faux a été imprimé cinq fois dans le manuscrit par ce chemin, et l'épinglage ne pouvait
    # pas le voir puisqu'il comparait l'imprimé à ce TSV déjà arrondi.
    df.to_csv(a.sortie, sep="\t", index=False, float_format="%.6f")
    rs = Path(a.sortie).with_name(Path(a.sortie).stem + "_residuel.tsv")
    pd.DataFrame(resid).sort_values(["drug", "n_souches"], ascending=[True, False]).to_csv(rs, sep="\t", index=False)
    with pd.option_context("display.width", 220, "display.max_rows", 500):
        v = df[df.cible == "reference"]
        print(f"=== {a.etiquette} | {a.modeldir} | cible de référence ===")
        print(v[v.etat == "complet"][["drug", "n_R", "n_S", "sens", "spec", "n_col"]].to_string(index=False))
        print(v[v.etat != "complet"].pivot(index="drug", columns="etat", values="gain_sens")
              .reindex(p74.DRUGS).round(3).to_string())
    print(f"-> {a.sortie} ; {rs}")


if __name__ == "__main__":
    main()
