#!/usr/bin/env python3
"""Phase 11 — Audit du baseline catalogue déterministe par médicament.

Mesure la sensibilité/spécificité du baseline « porte ≥1 déterminant R catalogué »
(réconcilié SNP + MNV + indel via variant_match) sur la matrice souche×SPDI, après
correction des quatre angles morts de représentation. Table finale pour le manuscrit.

Optimisé : signatures de matching pré-calculées une fois pour toutes les features.

P5.6 (2026-10-18) : la sensibilité dépend de la CIBLE (sources de phénotype gardées). Le fichier
écrit est de format LARGE, une ligne par drogue : les colonnes `n_R, n_S, sens, spec` sont lues
sur la cible de RÉFÉRENCE (sans `souches_datees`, étiquettes dérivées du génotype), et les cibles
`complete` (toutes sources, celle d'avant P5.6) et `mesuree` (sans `souches_datees` ni PATRIC)
suivent en colonnes suffixées. Règle de source : `cible_reference.py`. Ancien fichier :
`baseline_audit_avant_P5.6.tsv`. P8.6 : une paire dont la source de consensus est écartée prend
l'étiquette de la source suivante (repli) ; fichier d'avant : `baseline_audit_avant_P8.6.tsv`. Usage : phase11_baseline_audit.py [--cible reference|complete|mesuree]

AA1 (2026-11-11) : les règles d'interprétation de l'OMS portées par la colonne `regle` du catalogue (épistasie
Rv0678/mmpL5 pour BDQ et CFZ, eis pour AMK et KAN ; `regles_interpretation.py`) sont appliquées dans les colonnes
`sens_regles`, `spec_regles` (cible choisie) et `n_neutralises` (souches R du baseline dont aucun déterminant ne survit à
la règle) ; les colonnes `sens`, `spec` restent le baseline SANS règle, comparable aux fichiers antérieurs. Pour un
médicament sans règle, `*_regles` = `*`. BDQ et CFZ entrent dans l'audit. Fichier d'avant : `baseline_audit_avant_AA1.tsv`.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
import variant_match as vm
import cible_reference as cr
from regles_interpretation import couple_regle

MODELDIR = paths.RESULTATS / "model"
DRUGS = ["isoniazid", "rifampicin", "ethambutol", "pyrazinamide", "moxifloxacin",
         "levofloxacin", "amikacin", "kanamycin", "streptomycin", "ethionamide",
         "capreomycin", "bedaquiline", "clofazimine"]


def signatures(feats):
    """Signature de matching par feature (une seule fois) :
    ('x', spdi) pour SNP/indel (match exact) ; ('m', {SNP composants}) pour MNV."""
    sig = {}
    for s in feats:
        p = s.split(":")
        ref, alt = p[2], p[3]
        if len(ref) == len(alt) and len(ref) > 1:
            sig[s] = ("m", frozenset(vm.decompose(s)))
        else:
            sig[s] = ("x", s)
    return sig


def porte_spdi(X, feats, sig, spdis):
    """(vecteur 0/1 « porte au moins un des SPDI », nombre de features) ; SNP/indel exact, MNV par décomposition."""
    mono = {s for s in spdis if len(s.split(":")[2]) == 1 and len(s.split(":")[3]) == 1}
    cols = [j for j, s in enumerate(feats)
            if (sig[s][0] == "x" and sig[s][1] in spdis)
            or (sig[s][0] == "m" and (sig[s][1] & mono))]
    porte = (X[:, cols].sum(axis=1).A1 > 0).astype(int) if cols else np.zeros(X.shape[0], int)
    return porte, len(cols)


def porteurs(X, feats, sig, cat, drug):
    """(vecteur 0/1 « porte au moins un déterminant R catalogué de `drug` », nombre de features)."""
    return porte_spdi(X, feats, sig, vm.determinant_snp_set(cat, drug))


def porteurs_regles(X, feats, sig, cat, drug):
    """Vecteur 0/1 du baseline AVEC les règles d'abrogation (AA1) : le déterminant du gène soumis est neutralisé chez
    une souche qui porte un variant abrogeant, les déterminants des autres gènes restent. Médicament sans règle :
    le vecteur du baseline."""
    couple = couple_regle(cat, drug)
    if couple is None:
        return porteurs(X, feats, sig, cat, drug)[0]
    soumis, autres, abro = couple
    c_soumis = porte_spdi(X, feats, sig, soumis)[0].astype(bool)
    c_autres = porte_spdi(X, feats, sig, autres)[0].astype(bool)
    c_abro = porte_spdi(X, feats, sig, abro)[0].astype(bool)
    return ((c_soumis & ~c_abro) | c_autres).astype(int)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cible", default=cr.CIBLE_DEFAUT, choices=list(cr.CIBLES),
                    help="cible lue dans les colonnes sans suffixe (défaut : référence)")
    cible = ap.parse_args().cible
    X = sparse.load_npz(MODELDIR / "X.npz").tocsr()
    feats = (MODELDIR / "features.txt").read_text().splitlines()
    rows = pd.read_csv(MODELDIR / "rows.tsv", sep="\t", keep_default_na=False)
    cat = pd.read_csv(paths.CATALOGUE_TSV, sep="\t").fillna("")
    sig = signatures(feats)

    out = []
    for drug in DRUGS:
        if not rows[drug].isin(["R", "S"]).any():
            continue
        porte, n_cols = porteurs(X, feats, sig, cat, drug)
        porte_r = porteurs_regles(X, feats, sig, cat, drug)
        ligne = {"drug": drug}
        mesures = {}
        for nom in cr.CIBLES:
            lab = cr.etiquettes_cible(rows, drug, nom)       # étiquette sous la cible, avec repli (P8.6)
            ym = lab.isin(["R", "S"]).values
            y = (lab[ym] == "R").astype(int).values
            m = {"n_R": int(y.sum()), "n_S": int((1 - y).sum())}
            for suf, vec in (("", porte), ("_regles", porte_r)):
                car = vec[ym]
                tp = int((y & car).sum()); fn = int((y & (1 - car)).sum())
                tn = int(((1 - y) & (1 - car)).sum()); fp = int(((1 - y) & car).sum())
                m["sens" + suf] = round(tp / (tp + fn), 3) if tp + fn else None
                m["spec" + suf] = round(tn / (tn + fp), 3) if tn + fp else None
            m["n_neutralises"] = int(((porte == 1) & (porte_r == 0))[ym].sum())
            mesures[nom] = m
        ligne.update(mesures[cible])
        ligne["det_features"] = n_cols
        for nom in cr.CIBLES:
            if nom != cible:
                ligne.update({f"{k}_{nom}": v for k, v in mesures[nom].items()})
        out.append(ligne)
    df = pd.DataFrame(out)
    df.to_csv(MODELDIR / "baseline_audit.tsv", sep="\t", index=False)
    print(df.to_string(index=False))
    print(f"\n-> {MODELDIR/'baseline_audit.tsv'}")


if __name__ == "__main__":
    main()
