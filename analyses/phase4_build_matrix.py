#!/usr/bin/env python3
"""Phase 4 — Matrice creuse souche × SPDI (features) + métadonnées alignées.

Construit, pour les souches ayant à la fois un phénotype (consensus) et des features
(pangénome `pan_strains.pkl`), une matrice binaire creuse CSR souche × SPDI, restreinte
aux SPDI de fréquence informative (ni ultra-rares, ni fixés partout). Sauve la matrice
et trois tables alignées sur l'ordre des lignes : souches, phénotypes (large), lignée Coll.

Sorties dans résultats/model/ :
  X.npz                 matrice creuse (CSR), souches × SPDI
  features.txt          liste des SPDI (colonnes), dans l'ordre
  rows.tsv              strain_id, lineage_level_1, lineage_code (lignes), + phénotype par drogue

Option --etiquettes-seules (P8.5) : ne charge ni le pangénome ni bdd/actuelle ; réécrit seulement
rows.tsv à partir du consensus courant, sur les souches et dans l'ordre de l'ancien rows.tsv
(X.npz et features.txt restent ceux de la dernière construction complète).

Option --p88 (P8.8) : reconstruction complète avec la cascade pangénome > bdd/actuelle (extrait
`p88_bdd_extrait.tsv.gz`, bdd n'étant pas sur mp) > tblearn (`--tblearn`, souches de P8.2), une
représentation CANONIQUE des variants pour toutes les sources (blocs MNV découpés en SNP, MNV du
catalogue reconstruits par conjonction de leurs composantes ; `p88_generation.py`) et un masque de
génération (`--masque`) retiré de toutes les souches. `--outdir` écrit ailleurs que dans
résultats/model (témoins A/B/C de la piste) ; `souches_source_spdi.tsv` dit d'où viennent les
variants de chaque souche. Sans `--tblearn`, c'est le témoin « mêmes souches, mêmes sources, nouvelle
représentation » qui isole l'effet de la reconstruction de celui des souches ajoutées.
"""
import pickle
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

MIN_FREQ = 20            # SPDI présent dans au moins 20 souches
MAX_FRAC = 0.99          # et dans au plus 99 % (drop des sites fixés, non informatifs)
OUTDIR = paths.RESULTATS / "model"


import subprocess

BDD = Path.home() / "docs/codes/mtbc/bdd/actuelle"


def build_bdd_index():
    """SRA -> chemin spdi.txt pour les souches présentes sur disque dans bdd/actuelle."""
    out = subprocess.run(["find", str(BDD), "-path", "*/NC_000962.3/spdi.txt"],
                         capture_output=True, text=True).stdout
    idx = {}
    for p in out.splitlines():
        parts = p.split("/")
        # .../actuelle/<lignée>/<SRA>/NC_000962.3/spdi.txt
        if len(parts) >= 3:
            idx[parts[-3]] = p
    return idx


def construire_rows(strains, phe, alias=None):
    """rows.tsv : lignée Coll + une colonne R/S par drogue, alignées sur `strains` (ordre de X.npz).
    `alias` (P8.8) : isolat composite -> runs, dont la lignée Coll sert quand l'isolat n'a pas la sienne."""
    # métadonnées alignées : phénotypes large + lignée
    wide = (phe[phe["strain_id"].isin(strains)]
            .pivot_table(index="strain_id", columns="drug", values="phenotype",
                         aggfunc="first"))
    wide = wide.reindex(strains)
    lin = pd.read_csv(paths.SRC_TBANNOTATOR / "lineage_coll_snapshot.csv")
    lin = lin.drop_duplicates("strain_name").set_index("strain_name")
    rows = pd.DataFrame({"strain_id": strains})
    rows["lineage_level_1"] = rows["strain_id"].map(lin["lineage_level_1"]).fillna("NA").astype(str)
    rows["lineage_code"] = rows["strain_id"].map(lin["lineage_code"]).fillna("NA").astype(str)
    # P1.7 : repli sur bdd/actuelle (taxonomie du groupe, niveau 1 seulement) pour les souches que Coll ne couvre pas
    comp = paths.SRC_TBANNOTATOR / "lineage_bdd_complement.csv"
    if comp.exists():  # entree-facultative : absent avant P1.7
        c = pd.read_csv(comp, dtype=str, keep_default_na=False).drop_duplicates("strain_name").set_index("strain_name")
        m = rows["lineage_level_1"] == "NA"
        rows.loc[m, "lineage_code"] = rows.loc[m, "strain_id"].map(c["lineage_code"]).fillna("NA")
        rows.loc[m, "lineage_level_1"] = rows.loc[m, "strain_id"].map(c["lineage_level_1"]).fillna("NA")
    for iso, runs in (alias or {}).items():
        m = rows["strain_id"] == iso
        for r in runs:
            if m.any() and rows.loc[m, "lineage_level_1"].iloc[0] == "NA" and r in lin.index:
                rows.loc[m, "lineage_level_1"] = str(lin.loc[r, "lineage_level_1"])
                rows.loc[m, "lineage_code"] = str(lin.loc[r, "lineage_code"])
    rows = rows.join(wide.reset_index(drop=True))
    # P5.3 : cible PZA nettoyée (M. bovis + BioProject PRJNA1237251), valeur d'origine gardée
    # dans `pyrazinamide_brut` ; no-op si le registre ENA n'est pas construit.
    from cible_pza import REGISTRE, appliquer
    if REGISTRE.exists():
        rows = appliquer(rows, pd.read_csv(REGISTRE, sep="\t"))
    return rows


def etiquettes_seules():
    """P8.5 : réécrit rows.tsv avec les étiquettes du consensus COURANT, sans toucher à X.npz ni
    à features.txt. Les souches et leur ordre sont ceux de l'ancien rows.tsv (alignés sur X.npz),
    si bien que seule l'étiquette change : une reconstruction complète changerait aussi les
    variables (X.npz et features.txt datent du 31 mai, catalogue d'alors, base d'alors) et
    mélangerait deux effets. L'ancien rows.tsv est gardé en rows_avant_P8.5.tsv (une seule fois)."""
    ancien = OUTDIR / "rows.tsv"
    copie = OUTDIR / "rows_avant_P8.5.tsv"
    if not copie.exists():
        copie.write_bytes(ancien.read_bytes())
        print(f"ancien rows.tsv gardé -> {copie.name}")
    old = pd.read_csv(ancien, sep="\t", keep_default_na=False, usecols=["strain_id"])
    strains = old["strain_id"].tolist()
    phe = pd.read_csv(paths.CATALOGUE_DIR / "phenotypes_souches_consensus.tsv", sep="\t")
    rows = construire_rows(strains, phe)
    rows.to_csv(ancien, sep="\t", index=False)
    drogues = [c for c in rows.columns if c not in ("strain_id", "lineage_level_1", "lineage_code",
                                                    "pyrazinamide_brut")]
    sans = int((~rows[drogues].isin(["R", "S"]).any(axis=1)).sum())
    print(f"rows.tsv réécrit : {len(rows)} souches (ordre de X.npz conservé), {sans} sans aucune étiquette R/S")
    print(f"artefacts -> {OUTDIR} (X.npz et features.txt inchangés)")


def construire_p88(args):
    """P8.8 : matrice sur la cohorte élargie, représentation canonique, masque de génération."""
    import gc
    import p88_generation as g

    out = Path(args.outdir) if args.outdir else OUTDIR
    out.mkdir(parents=True, exist_ok=True)
    phe = pd.read_csv(paths.CATALOGUE_DIR / "phenotypes_souches_consensus.tsv", sep="\t")
    phe_strains = set(phe["strain_id"].dropna().unique())
    conj, cat_r = g.charger_catalogue(paths.CATALOGUE_TSV)
    masque = set()
    if args.masque:
        masque = set(pd.read_csv(args.masque, sep="\t")["variant"]) - cat_r
    print(f"{conj.n} MNV catalogués reconstruits par conjonction ; {len(masque)} variants masqués")

    vocab = {}
    jeux, source = {}, {}

    def ajouter(s, brut, src):
        # --brut : blocs MNV de l'ancien pipeline conservés (mesure de l'angle mort MNV par p74, qui n'existe
        # plus dans une matrice canonique) ; jamais avec --tblearn, dont les SNP scindés ne se recomposent pas
        ens = set(brut) if getattr(args, "brut", False) else g.normaliser(brut, conj)
        if masque:
            ens -= masque
        arr = np.empty(len(ens), dtype=np.int32)
        for k, v in enumerate(ens):
            arr[k] = vocab.setdefault(v, len(vocab))
        jeux[s] = arr
        source[s] = src

    print("chargement pan_strains.pkl…")
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)
    for s in sorted(phe_strains & set(pan)):
        ajouter(s, pan[s], "pangenome")
    del pan
    gc.collect()
    print(f"  pangénome : {len(jeux)} souches")

    n0 = len(jeux)
    for s, src, ens in g.iter_tsv_source(args.bdd_extrait):
        if s in phe_strains and s not in jeux:
            ajouter(s, ens, "bdd")
    print(f"  bdd/actuelle : {len(jeux) - n0} souches")

    if args.tblearn:
        n0 = len(jeux)
        comp = pd.read_csv(args.composites, sep="\t") if args.composites else None
        runs_comp = {}
        if comp is not None:
            for iso, r in zip(comp["strain_id"], comp["run"]):
                if iso in phe_strains and iso not in jeux:
                    runs_comp.setdefault(iso, []).append(r)
        besoin = set(r for rs in runs_comp.values() for r in rs)
        par_run = {}
        for s, ens in g.iter_tsv_spdi(args.tblearn):
            if s in phe_strains and s not in jeux:
                ajouter(s, ens, "tblearn")
            elif s in besoin:
                par_run[s] = ens
        n_simples = len(jeux) - n0
        for iso, ens in g.fusionner_composites(par_run, runs_comp).items():
            ajouter(iso, ens, "tblearn_composite")
        print(f"  tblearn : {n_simples} souches simples + {len(jeux) - n0 - n_simples} isolats composites")

    strains = sorted(jeux)
    n = len(strains)
    print(f"souches utilisables : {n} sur {len(phe_strains)} phénotypées ; sources {dict(Counter(source.values()))}")
    inv = [None] * len(vocab)
    for v, i in vocab.items():
        inv[i] = v
    freq = np.bincount(np.concatenate([jeux[s] for s in strains]), minlength=len(vocab))
    lo, hi = MIN_FREQ, MAX_FRAC * n
    garde = (freq >= lo) & (freq <= hi)
    rare = [vocab[v] for v in cat_r if v in vocab and freq[vocab[v]] > 0 and not garde[vocab[v]]]
    garde[rare] = True
    features = sorted(inv[i] for i in np.flatnonzero(garde))
    col = np.full(len(vocab), -1, dtype=np.int64)
    for j, v in enumerate(features):
        col[vocab[v]] = j
    print(f"features retenues : {len(features)} SPDI (freq∈[{lo},{int(hi)}] + {len(rare)} déterminants R catalogués rares)")

    indptr, indices = [0], []
    for s in strains:
        c = col[jeux[s]]
        indices.append(c[c >= 0].astype(np.int32))
        indptr.append(indptr[-1] + int((c >= 0).sum()))
    indices = np.concatenate(indices)
    X = sparse.csr_matrix((np.ones(len(indices), dtype=np.int8), indices, np.array(indptr, dtype=np.int64)),
                          shape=(n, len(features)))
    sparse.save_npz(out / "X.npz", X)
    (out / "features.txt").write_text("\n".join(features))
    print(f"X : {X.shape[0]} × {X.shape[1]} | nnz={X.nnz} ({100*X.nnz/(X.shape[0]*X.shape[1]):.2f} % densité)")

    pd.DataFrame({"strain_id": strains, "source_spdi": [source[s] for s in strains]}).to_csv(
        out / "souches_source_spdi.tsv", sep="\t", index=False)
    rows = construire_rows(strains, phe, alias=_alias_composites(args.composites, phe_strains))
    rows.to_csv(out / "rows.tsv", sep="\t", index=False)
    cov = (rows["lineage_level_1"] != "NA").sum()
    print(f"couverture lignée Coll : {cov}/{n} ({100*cov/n:.0f} %)")
    print("répartition lignées :", rows["lineage_level_1"].value_counts().to_dict())
    print(f"\nartefacts -> {out}")


def _alias_composites(chemin, phe_strains):
    """isolat CRyPTIC composite -> liste de ses runs (la lignée Coll se lit sur le run, pas sur l'isolat)."""
    if not chemin:
        return {}
    comp = pd.read_csv(chemin, sep="\t")
    alias = {}
    for iso, r in zip(comp["strain_id"], comp["run"]):
        if iso in phe_strains:
            alias.setdefault(iso, []).append(r)
    return alias


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    phe = pd.read_csv(paths.CATALOGUE_DIR / "phenotypes_souches_consensus.tsv", sep="\t")
    phe_strains = set(phe["strain_id"].dropna().unique())

    print("chargement pan_strains.pkl…")
    with open(paths.PAN_STRAINS, "rb") as f:
        pan = pickle.load(f)
    pan_strains = set(pan.keys())

    print("indexation de bdd/actuelle…")
    bdd_idx = build_bdd_index()
    print(f"  {len(bdd_idx)} souches sur disque")

    # cascade : features depuis le pangénome, sinon depuis bdd/actuelle (spdi.txt)
    _bdd_cache = {}

    def spdi_set(s):
        if s in pan:
            return pan[s]
        if s not in _bdd_cache:
            try:
                _bdd_cache[s] = set(Path(bdd_idx[s]).read_text().split())
            except Exception:
                _bdd_cache[s] = set()
        return _bdd_cache[s]

    usable = phe_strains & (pan_strains | set(bdd_idx))
    strains = sorted(usable)
    n = len(strains)
    n_pan = len(phe_strains & pan_strains)
    n_bdd = len(strains) - n_pan
    print(f"souches utilisables : {n}  (pangénome {n_pan} + bdd {n_bdd})")

    # pass 1 : fréquence des SPDI
    print("comptage des fréquences SPDI…")
    freq = Counter()
    for s in strains:
        freq.update(spdi_set(s))
    lo, hi = MIN_FREQ, MAX_FRAC * n
    feat = set(spdi for spdi, c in freq.items() if lo <= c <= hi)
    # toujours inclure les déterminants R catalogués OBSERVÉS, même sous le seuil de
    # fréquence : sinon les déterminants de résistance rares (FQ, injectables, BDQ…)
    # sont exclus et le baseline/modèle les rate (cf. cahier : angle mort MXF).
    cat = pd.read_csv(paths.CATALOGUE_TSV, sep="\t").fillna("")
    catR = {s for s in cat.loc[cat["call"] == "R-associated", "spdi"] if s}
    rare_det = {s for s in catR if s in freq and s not in feat}
    feat |= rare_det
    features = sorted(feat)
    fidx = {spdi: j for j, spdi in enumerate(features)}
    print(f"features retenues : {len(features)} SPDI "
          f"(freq∈[{lo},{int(hi)}] + {len(rare_det)} déterminants R catalogués rares)")

    # pass 2 : matrice creuse
    print("construction de la matrice creuse…")
    indptr = [0]
    indices = []
    for s in strains:
        cols = [fidx[v] for v in spdi_set(s) if v in fidx]
        indices.extend(cols)
        indptr.append(len(indices))
    data = np.ones(len(indices), dtype=np.int8)
    X = sparse.csr_matrix((data, np.array(indices, dtype=np.int32),
                           np.array(indptr, dtype=np.int64)),
                          shape=(n, len(features)))
    sparse.save_npz(OUTDIR / "X.npz", X)
    (OUTDIR / "features.txt").write_text("\n".join(features))
    print(f"X : {X.shape[0]} × {X.shape[1]} | nnz={X.nnz} "
          f"({100*X.nnz/(X.shape[0]*X.shape[1]):.2f} % densité)")

    rows = construire_rows(strains, phe)
    rows.to_csv(OUTDIR / "rows.tsv", sep="\t", index=False)

    cov = (rows["lineage_level_1"] != "NA").sum()
    print(f"couverture lignée Coll sur ces souches : {cov}/{n} ({100*cov/n:.0f} %)")
    print("répartition lignées :", rows["lineage_level_1"].value_counts().to_dict())
    print(f"\nartefacts -> {OUTDIR}")


if __name__ == "__main__":
    if "--etiquettes-seules" in sys.argv[1:]:
        etiquettes_seules()
    elif "--p88" in sys.argv[1:]:
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--p88", action="store_true")
        ap.add_argument("--bdd-extrait", required=True)
        ap.add_argument("--tblearn", help="TSV gzip de p82_fetch_spdi_tblearn.py (souches ajoutées)")
        ap.add_argument("--composites", help="p82_composites.tsv (isolat CRyPTIC -> runs)")
        ap.add_argument("--masque", help="p88_masque_generation.tsv")
        ap.add_argument("--outdir")
        ap.add_argument("--brut", action="store_true", help="sans représentation canonique (blocs MNV gardés)")
        a = ap.parse_args()
        if a.brut and a.tblearn:
            ap.error("--brut est réservé à l'ancienne génération : il exclut --tblearn")
        construire_p88(a)
    else:
        main()
