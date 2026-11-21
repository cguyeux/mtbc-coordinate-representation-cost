#!/usr/bin/env python3
"""P7.4 — Baseline catalogue et gain des quatre angles morts, rejoués sur chaque cible.

Objet      : recalculer le tableau de sensibilité du baseline déterministe (onze médicaments) sur la
             cible MESURÉE SEULE (sans `souches_datees` ni PATRIC, cf. `cible_reference.py`) et
             rejouer, sur la même cible, le gain de chacune des quatre corrections de représentation
             du manuscrit. Le gain d'une correction est mesuré par RETRAIT : sensibilité du baseline
             complet moins sensibilité du baseline sans cette correction (les trois autres restant
             actives). Un état « naïf » (les quatre retirées) borne l'ensemble. Les gains historiques
             du manuscrit (EMB +23, MXF +55, LEV +51, PZA +7) étaient séquentiels et mesurés sur la
             cible complète d'alors ; les états « avant » n'avaient pas été conservés (claim_check C9).
Retraits   : mnv       les MNV observés ne sont plus décomposés en SNP (appariement exact seul) ;
             indel     aucune feature indel n'est appariée ;
             synonymes une seule représentation nucléotidique par variant protéique du catalogue
                       (mono-base de préférence, sinon la première rencontrée), qui est le
                       comportement documenté du résolveur d'avant la propagation AA ; variante
                       `synonymes_tirage` : représentation tirée au hasard, 20 tirages (moyenne,
                       min, max) pour borner l'arbitraire du choix ;
             rares     les déterminants présents dans moins de MIN_FREQ=20 souches (seuil de
                       `phase4`) sont retirés, comme si le filtre de fréquence les avait élagués.
Incertitude: sensibilité et spécificité : Clopper-Pearson ; gains : rééchantillonnage par grappes de
             sous-lignée (`lineage_code`, 40 grappes), les mêmes tirages pour le complet et le retrait
             (gain apparié), 2 000 tirages, graine fixe.
Entrées    : results/model/{X.npz,features.txt,rows.tsv}, catalogue/catalogue_consolide.tsv.
Sorties    : results/model/p74_angles_morts.tsv (format long), p74_tableau_mesuree.tex (tableau hors
             PATRIC) et, depuis P8.5, p74_tableau_reference.tex (cible de référence, PATRIC corrigé inclus,
             avec la sensibilité hors PATRIC en dernière colonne).
Limite     : X.npz date de la construction de `phase4` ; ses colonnes « déterminants rares » sont
             celles du catalogue de l'époque. Le calcul est donc la même mesure que `phase11`.
Réutilisable : oui (toute cible nommée de `cible_reference.py`).
Projet     : Resistance_antibio (P7.4)
Date       : 2026-10-18
Usage      : python3 analyses/p74_angles_morts.py [--tirages 2000]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
import variant_match as vm
import cible_reference as cr

MODELDIR = paths.RESULTATS / "model"
DRUGS = ["isoniazid", "rifampicin", "ethambutol", "pyrazinamide", "moxifloxacin",
         "levofloxacin", "amikacin", "kanamycin", "streptomycin", "ethionamide",
         "capreomycin"]
MIN_FREQ = 20
RETRAITS = ["mnv", "indel", "synonymes", "rares", "naif"]


def clopper_pearson(k, n, alpha=0.05):
    if n == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return (float(lo), float(hi))


def est_snp(s):
    _, ref, alt = vm.parse_spdi(s)
    return len(ref) == 1 and len(alt) == 1


def est_indel(s):
    _, ref, alt = vm.parse_spdi(s)
    return len(ref) != len(alt)


def une_representation(cat_drug, rng=None):
    """SPDI R-associés gardés quand chaque variant protéique n'a qu'UNE représentation.

    `rng` None : mono-base de préférence puis première rencontrée (déterministe) ;
    sinon : tirage uniforme parmi les représentations du variant protéique.
    Les SPDI qui ne relèvent d'aucun variant protéique (c., empiriques) et tous les indels sont gardés tels quels.
    """
    r = cat_drug[(cat_drug["call"] == "R-associated") & (cat_drug["spdi"] != "")]
    # seules les SUBSTITUTIONS d'un variant protéique sont réduites : les normalisations d'un
    # indel (même variant protéique, coordonnées différentes) relèvent de l'angle mort 4, pas du 1
    indel = r["spdi"].map(est_indel)
    prot = r["variant"].str.contains(r"_p\.", regex=True) & ~indel
    garde = set(r.loc[~prot, "spdi"])
    for _, g in r[prot].groupby("variant", sort=False):
        reps = list(dict.fromkeys(g["spdi"]))
        if len(reps) == 1:
            garde.add(reps[0])
        elif rng is None:
            mono = [s for s in reps if est_snp(s)]
            garde.add((mono or reps)[0])
        else:
            garde.add(reps[rng.integers(len(reps))])
    return garde


def colonnes(feats, sig, det, mnv_decompose=True, indels=True, freq=None):
    mono = {s for s in det if est_snp(s)}
    cols = []
    for j, s in enumerate(feats):
        typ, val = sig[s]
        if typ == "x":
            ok = val in det and (indels or not est_indel(s))
        else:
            ok = (val & mono) if mnv_decompose else (s in det)
            ok = bool(ok)
        if ok and (freq is None or freq[j] >= MIN_FREQ):
            cols.append(j)
    return cols


def comptes_par_grappe(porte, y, grappes, ncl):
    """Tableau (grappes × 4) : tp, fn, tn, fp."""
    out = np.zeros((ncl, 4), dtype=np.int64)
    for k, mask in ((0, (y == 1) & (porte == 1)), (1, (y == 1) & (porte == 0)),
                    (2, (y == 0) & (porte == 0)), (3, (y == 0) & (porte == 1))):
        out[:, k] = np.bincount(grappes[mask], minlength=ncl)
    return out


def sens_spec_depuis(c):
    c = np.atleast_2d(c)
    tp, fn, tn, fp = c[:, 0], c[:, 1], c[:, 2], c[:, 3]
    with np.errstate(divide="ignore", invalid="ignore"):
        return tp / (tp + fn), tn / (tn + fp)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tirages", type=int, default=2000)
    args = ap.parse_args()

    X = sparse.load_npz(MODELDIR / "X.npz").tocsr()
    feats = (MODELDIR / "features.txt").read_text().splitlines()
    rows = pd.read_csv(MODELDIR / "rows.tsv", sep="\t", keep_default_na=False)
    assert X.shape[0] == len(rows) and X.shape[1] == len(feats), "X.npz et rows.tsv désalignés"
    cat = pd.read_csv(paths.CATALOGUE_TSV, sep="\t").fillna("")
    freq = np.asarray(X.sum(axis=0)).ravel()

    sig = {}
    for s in feats:
        _, ref, alt = vm.parse_spdi(s)
        if len(ref) == len(alt) and len(ref) > 1:
            sig[s] = ("m", frozenset(vm.decompose(s)))
        else:
            sig[s] = ("x", s)

    cl_noms, cl_idx = np.unique(rows["lineage_code"].values, return_inverse=True)
    ncl = len(cl_noms)
    rng_boot = np.random.default_rng(0)
    pick = rng_boot.integers(0, ncl, (args.tirages, ncl))

    long = []
    for drug in DRUGS:
        cat_d = cat[cat["drug"] == drug]
        det = vm.determinant_snp_set(cat, drug)
        det_1rep = une_representation(cat_d)
        porte = {}
        porte["complet"] = _porte(X, colonnes(feats, sig, det))
        porte["mnv"] = _porte(X, colonnes(feats, sig, det, mnv_decompose=False))
        porte["indel"] = _porte(X, colonnes(feats, sig, det, indels=False))
        porte["synonymes"] = _porte(X, colonnes(feats, sig, det_1rep))
        porte["rares"] = _porte(X, colonnes(feats, sig, det, freq=freq))
        porte["naif"] = _porte(X, colonnes(feats, sig, det_1rep, mnv_decompose=False,
                                           indels=False, freq=freq))
        tirages = []
        for i in range(20):
            d_i = une_representation(cat_d, np.random.default_rng(100 + i))
            tirages.append(_porte(X, colonnes(feats, sig, d_i)))
        for nom in cr.CIBLES:
            lab = cr.etiquettes_cible(rows, drug, nom)       # étiquette sous la cible, avec repli (P8.6)
            ym = lab.isin(["R", "S"]).values
            if not ym.any():
                continue
            y = (lab[ym] == "R").astype(int).values
            g = cl_idx[ym]
            cnt = {k: comptes_par_grappe(v[ym], y, g, ncl) for k, v in porte.items()}
            tot = {k: v.sum(axis=0) for k, v in cnt.items()}
            boot = {k: v[pick].sum(axis=1) for k, v in cnt.items()}      # (tirages, 4)
            se_c, sp_c = sens_spec_depuis(tot["complet"])
            bse_c, _ = sens_spec_depuis(boot["complet"])
            for k in ["complet"] + RETRAITS:
                se, sp = sens_spec_depuis(tot[k])
                tp, fn, tn, fp = (int(x) for x in tot[k])
                lo, hi = clopper_pearson(tp, tp + fn)
                slo, shi = clopper_pearson(tn, tn + fp)
                ligne = {"cible": nom, "drug": drug, "etat": k, "n_R": tp + fn, "n_S": tn + fp,
                         "sens": float(se[0]), "sens_lo": lo, "sens_hi": hi,
                         "spec": float(sp[0]), "spec_lo": slo, "spec_hi": shi}
                if k != "complet":
                    bse, _ = sens_spec_depuis(boot[k])
                    delta = bse_c - bse
                    ligne.update({"gain_sens": float(se_c[0] - se[0]),
                                  "gain_lo": float(np.nanpercentile(delta, 2.5)),
                                  "gain_hi": float(np.nanpercentile(delta, 97.5)),
                                  "perte_spec": float(sp_c[0] - sp[0])})
                long.append(ligne)
            # variante tirage aléatoire des représentations (borne l'arbitraire du choix)
            gains = []
            for p in tirages:
                c = comptes_par_grappe(p[ym], y, g, ncl).sum(axis=0)
                gains.append(float(se_c[0] - sens_spec_depuis(c)[0][0]))
            long.append({"cible": nom, "drug": drug, "etat": "synonymes_tirage",
                         "n_R": int(tot["complet"][0] + tot["complet"][1]),
                         "n_S": int(tot["complet"][2] + tot["complet"][3]),
                         "gain_sens": float(np.mean(gains)), "gain_lo": float(np.min(gains)),
                         "gain_hi": float(np.max(gains))})
    df = pd.DataFrame(long)
    # P7.16 : %.6f et NON %.4f. Les colonnes de sensibilité et de gain sont des FRACTIONS :
    # arrondies à 4 décimales, elles ne gardent que 2 décimales utiles une fois lues en POINTS
    # (0,117472 -> 0,1175 -> 11,75, qui arrondit à 11,8 quand la valeur vaut 11,7). Un chiffre
    # faux a été imprimé cinq fois dans le manuscrit par ce chemin, et l'épinglage ne pouvait
    # pas le voir puisqu'il comparait l'imprimé à ce TSV déjà arrondi.
    df.to_csv(MODELDIR / "p74_angles_morts.tsv", sep="\t", index=False, float_format="%.6f")

    # effectifs de souches par cible (au moins une étiquette R/S sur les onze médicaments)
    eff = {}
    for nom in cr.CIBLES:
        m = np.zeros(len(rows), bool)
        for d in DRUGS:
            m |= cr.masque_cible(rows, d, nom)
        eff[nom] = int(m.sum())
    print("souches avec au moins une étiquette R/S (onze médicaments) :", eff)
    _ecrire_tableau(df)
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        v = df[(df.cible == "reference")]
        print(v[v.etat == "complet"][["drug", "n_R", "n_S", "sens", "sens_lo", "sens_hi", "spec"]]
              .to_string(index=False))
        print(v[v.etat != "complet"].pivot(index="drug", columns="etat", values="gain_sens")
              .reindex(DRUGS).round(3).to_string())
    print(f"\n-> {MODELDIR/'p74_angles_morts.tsv'}")


def _porte(X, cols):
    return (X[:, cols].sum(axis=1).A1 > 0).astype(np.int8) if cols else np.zeros(X.shape[0], np.int8)


def _ecrire_tableau(df):
    noms = {"isoniazid": "Isoniazid", "rifampicin": "Rifampicin", "ethambutol": "Ethambutol",
            "pyrazinamide": "Pyrazinamide", "moxifloxacin": "Moxifloxacin",
            "levofloxacin": "Levofloxacin", "ethionamide": "Ethionamide",
            "kanamycin": "Kanamycin", "amikacin": "Amikacin", "streptomycin": "Streptomycin",
            "capreomycin": "Capreomycin"}
    ecrire = lambda cible: df[(df.cible == cible) & (df.etat == "complet")].set_index("drug")
    v = ecrire("mesuree")
    ordre = sorted(v.index, key=lambda d: -v.loc[d, "sens"])
    out = []
    for d in ordre:
        r = v.loc[d]
        out.append(f"{noms[d]:<13} & {int(r.n_R):>5} & {int(r.n_S):>5} & "
                   f"{r.sens:.3f} ({r.sens_lo:.2f}--{r.sens_hi:.2f}) & {r.spec:.3f} \\\\")
    (MODELDIR / "p74_tableau_mesuree.tex").write_text("\n".join(out) + "\n")
    # P8.5 : tableau de la cible de référence (PATRIC corrigé inclus) avec, en dernière colonne, la
    # sensibilité hors PATRIC (effectif R entre parenthèses) pour montrer ce qui dépend de cette source.
    ref, hors = ecrire("reference"), v
    ordre = sorted(ref.index, key=lambda d: -ref.loc[d, "sens"])
    out = []
    for d in ordre:
        r = ref.loc[d]
        h = f"{hors.loc[d, 'sens']:.3f} ({int(hors.loc[d, 'n_R'])})" if d in hors.index else "--"
        out.append(f"{noms[d]:<13} & {int(r.n_R):>5} & {int(r.n_S):>5} & "
                   f"{r.sens:.3f} ({r.sens_lo:.2f}--{r.sens_hi:.2f}) & {r.spec:.3f} & {h} \\\\")
    (MODELDIR / "p74_tableau_reference.tex").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
