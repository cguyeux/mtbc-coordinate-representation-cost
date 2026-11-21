#!/usr/bin/env python3
"""Phase 2 — Table consolidée des phénotypes par souche (cible supervisée).

Rassemble les sources souche -> drogue -> phénotype en une table longue unique.
C'est la *cible* du futur modèle prédictif : les *features* (mutations SPDI de
chaque souche) viennent du pangénome (data/pangenome/pan_strains.pkl) ou de
TBannotator, joints sur strain_id.

Sources :
  - CRyPTIC : ~12 287 isolats, 13 drogues, phénotype binaire R/S (+ MIC). Clé ENA_RUN.
  - CRyPTIC_v3 : jeu élargi (Zenodo 22128395 v3.4.1), ~54 206 souches, 21 drogues,
    R/S sans MIC, douze cohortes de DST. Clé ENA_RUN. Branché le 2026-11-11
    (décision de CG, piste P1.3) ; tout chiffre tiré de cette source se donne par
    cohorte, cf. `REVUE_CRITIQUE.md` §2.6.
  - PATRIC  : ~4 359 souches, 20 drogues (résultats/literature_extraction/patric.csv). Clé SRA.
  - Littérature : extraction LLM (résultats/literature_extraction/res.csv). Clé SRA.

Schéma de sortie (catalogue/phenotypes_souches.tsv, TSV) :
  strain_id drug drug_code phenotype mic source
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths
from vocab import canonical_drug

COLUMNS = ["strain_id", "drug", "drug_code", "phenotype", "mic", "source"]


def _norm_pheno(v):
    """Ramène un phénotype hétérogène à 'R' / 'S' / None.

    Gère : 0/1, R/S, resistant/susceptible, Y/N (souches datées),
    et le format `CODE-R`/`CODE-S` de l'antibiogramme TBannotator (ex. 'INH-R').
    """
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip().lower()
    if s.endswith("-r"):
        return "R"
    if s.endswith("-s"):
        return "S"
    if s in ("r", "1", "resistant", "1.0", "true", "y", "yes"):
        return "R"
    if s in ("s", "0", "susceptible", "0.0", "false", "n", "no"):
        return "S"
    return None


def ingest_cryptic() -> pd.DataFrame:
    f = paths.SRC_CRYPTIC / "CRyPTIC_reuse_table_20240917.csv"
    if not f.exists():
        print(f"[cryptic] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    c = pd.read_csv(f, low_memory=False)
    bin_cols = [x for x in c.columns if x.endswith("_BINARY_PHENOTYPE")]
    rows = []
    for col in bin_cols:
        code3 = col.split("_")[0]
        canon, code = canonical_drug(code3)
        mic_col = f"{code3}_MIC"
        sub = c[["ENA_RUN", col] + ([mic_col] if mic_col in c.columns else [])].dropna(subset=[col])
        for _, r in sub.iterrows():
            ph = _norm_pheno(r[col])
            if ph is None:
                continue
            rows.append({
                "strain_id": r["ENA_RUN"], "drug": canon, "drug_code": code,
                "phenotype": ph, "mic": r.get(mic_col, ""), "source": "CRyPTIC",
            })
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[cryptic] {len(df)} mesures | {df['strain_id'].nunique()} souches | "
          f"{df['drug'].nunique()} drogues")
    return df


def ingest_cryptic_v3() -> pd.DataFrame:
    """CRyPTIC élargi (« cryptic-dataset » v3.4.1, Zenodo 22128395), clé ENA_RUN.

    Lit la table normalisée par `ingestion/fetch_cryptic_v3.py` : 490 396 mesures
    R/S sur 54 206 souches et 21 médicaments, contre 12 287 souches et 13
    médicaments pour la reuse table de 2022 lue par `ingest_cryptic`.

    Les deux sources CRyPTIC sont gardées côte à côte, et ce n'est pas une
    redondance : la reuse table porte les MIC quantitatives des plaques UKMYC
    (colonne `mic`), que le jeu élargi ne redonne pas dans ses tables de DST (les
    MIC y vivent dans `UKMYC_PHENOTYPES.parquet`, non rapatrié). D'où la priorité
    de consensus : `CRyPTIC` d'abord, `CRyPTIC_v3` juste après.

    LIMITE À CONNAÎTRE. Le schéma de `phenotypes_souches.tsv` n'a pas de colonne
    de provenance fine, et l'élargir casserait tous les lecteurs en aval
    (`cible_reference.py`, `phase4`, `p85_propagation_patric.py`). La `source`
    vaut donc uniformément `CRyPTIC_v3`, alors que les mesures viennent de DOUZE
    cohortes de DST dont plusieurs ont servi à construire le catalogue OMS. La
    stratification par cohorte, obligatoire pour tout chiffre tiré de cette
    source (cf. `REVUE_CRITIQUE.md` §2.6), se fait en rejoignant
    `data/sources/cryptic/cryptic_v3_dst_by_sra.csv` sur (strain_id, drug_code),
    qui porte `quality`, `dst_source` et `cohort`.
    """
    f = paths.SRC_CRYPTIC / "cryptic_v3_dst_by_sra.csv"
    if not f.exists():
        print(f"[cryptic_v3] table normalisée absente ({f}) — "
              f"lancer ingestion/fetch_cryptic_v3.py")
        return pd.DataFrame(columns=COLUMNS)
    t = pd.read_csv(f, low_memory=False)
    t = t[t["phenotype"].isin(("R", "S"))].copy()
    df = pd.DataFrame({
        "strain_id": t["strain_id"],
        "drug": t["drug"],
        "drug_code": t["drug_code"],
        "phenotype": t["phenotype"],
        "mic": "",
        "source": "CRyPTIC_v3",
    }, columns=COLUMNS)
    print(f"[cryptic_v3] {len(df)} mesures | {df['strain_id'].nunique()} souches | "
          f"{df['drug'].nunique()} drogues")
    return df


def ingest_patric(fichier=None) -> pd.DataFrame:
    """PATRIC / BV-BRC, clé SRA. Lit `patric_corrige.csv` (P8.4, `p84_patric_etiquettes.py`) : seul
    un phénotype Resistant/Susceptible non vide donne une étiquette. L'ancien `patric.csv` codait R
    toute ligne sans phénotype (CMI seule) et laissait la dernière ligne d'un couple écraser les
    autres : il n'est lu qu'en repli, avec un avertissement. `fichier` impose le fichier lu (contrôle
    de P8.5 : rejouer l'ancien encodage)."""
    d = paths.RESULTATS / "literature_extraction"
    f = Path(fichier) if fichier else d / "patric_corrige.csv"
    if fichier is None and not f.exists():
        f = d / "patric.csv"
        print(f"[patric] ATTENTION : patric_corrige.csv absent, lecture de l'ancien {f.name} (étiquettes R fausses)")
    if not f.exists():
        print(f"[patric] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    p = pd.read_csv(f)
    drug_cols = [x for x in p.columns if x not in ("name", "sra")]
    rows = []
    for _, r in p.iterrows():
        sid = r.get("sra") or r.get("name")
        for col in drug_cols:
            ph = _norm_pheno(r[col])
            if ph is None:
                continue
            canon, code = canonical_drug(col)
            rows.append({"strain_id": sid, "drug": canon, "drug_code": code,
                         "phenotype": ph, "mic": "", "source": "PATRIC"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[patric] {len(df)} mesures | {df['strain_id'].nunique()} souches")
    return df


def ingest_literature() -> pd.DataFrame:
    f = paths.RESULTATS / "literature_extraction" / "res.csv"
    if not f.exists():
        print(f"[literature] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    r0 = pd.read_csv(f)
    drug_cols = [x for x in r0.columns if x not in ("SRA", "link")]
    rows = []
    for _, r in r0.iterrows():
        for col in drug_cols:
            ph = _norm_pheno(r[col])
            if ph is None:
                continue
            canon, code = canonical_drug(col)
            rows.append({"strain_id": r["SRA"], "drug": canon, "drug_code": code,
                         "phenotype": ph, "mic": "", "source": "literature"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[literature] {len(df)} mesures | {df['strain_id'].nunique()} souches")
    return df


def ingest_coll() -> pd.DataFrame:
    """Coll 2018 (GWAS Nat Genet) : 2655 isolats, 14 drogues, clé ENA (ERR), binaire 0/1."""
    import glob
    cand = glob.glob(str(paths.SRC_WHO / "Coll*GWAS*.xlsx"))
    if not cand:
        print("[coll] absent")
        return pd.DataFrame(columns=COLUMNS)
    c = pd.read_excel(cand[0])
    drug_cols = [x for x in c.columns if x not in ("Accession", "TOTAL R", "%R")]
    rows = []
    for _, r in c.iterrows():
        sid = r.get("Accession")
        if not isinstance(sid, str):
            continue
        for col in drug_cols:
            ph = _norm_pheno(r[col])
            if ph is None:
                continue
            canon, code = canonical_drug(col)
            rows.append({"strain_id": sid, "drug": canon, "drug_code": code,
                         "phenotype": ph, "mic": "", "source": "Coll2018"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[coll] {len(df)} mesures | {df['strain_id'].nunique()} souches | "
          f"{df['drug'].nunique()} drogues")
    return df


def ingest_souches_datees() -> pd.DataFrame:
    """« Toutes les souches datées » : 11001 isolats, 4 drogues 1re ligne, clé ENA, Y/N."""
    f = paths.METADATA / "Toutes les souches datées.csv"
    if not f.exists():
        print(f"[datees] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    t = pd.read_csv(f)
    drug_cols = [c for c in ("Isoniazid", "Rifampicin", "Ethambutol", "Pyrazinamide")
                 if c in t.columns]
    rows = []
    for _, r in t.iterrows():
        sid = r.get("Id")
        if not isinstance(sid, str):
            continue
        for col in drug_cols:
            ph = _norm_pheno(r[col])
            if ph is None:
                continue
            canon, code = canonical_drug(col)
            rows.append({"strain_id": sid, "drug": canon, "drug_code": code,
                         "phenotype": ph, "mic": "", "source": "souches_datees"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[datees] {len(df)} mesures | {df['strain_id'].nunique()} souches")
    return df


def ingest_tbannotator_dst() -> pd.DataFrame:
    """Antibiogramme TBannotator (DST MESURÉ, BioSample), clé SRA (SRR).

    Lit un snapshot CSV produit par `phase0_fetch_tbannotator_dst.py` :
        strain_name, antibiogram_inh, antibiogram_rif, antibiogram_pza, antibiogram_emb
    Encodage 'INH-R'/'INH-S'… géré par _norm_pheno. ~23 500 souches (4 drogues 1re ligne).
    """
    f = paths.SRC_TBANNOTATOR / "antibiogram_snapshot.csv"
    if not f.exists():
        print(f"[tbann-dst] snapshot absent ({f}) — lancer phase0_fetch_tbannotator_dst.py")
        return pd.DataFrame(columns=COLUMNS)
    t = pd.read_csv(f)
    drug_map = {"antibiogram_inh": "isoniazid", "antibiogram_rif": "rifampicin",
                "antibiogram_pza": "pyrazinamide", "antibiogram_emb": "ethambutol"}
    rows = []
    for _, r in t.iterrows():
        sid = r.get("strain_name")
        if not isinstance(sid, str):
            continue
        for col, canon in drug_map.items():
            ph = _norm_pheno(r.get(col))
            if ph is None:
                continue
            _, code = canonical_drug(canon)
            rows.append({"strain_id": sid, "drug": canon, "drug_code": code,
                         "phenotype": ph, "mic": "", "source": "TBannotator_biosample"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[tbann-dst] {len(df)} mesures | {df['strain_id'].nunique()} souches")
    return df


def ingest_tbportals() -> pd.DataFrame:
    """TB Portals (NIAID) : DST clinique MESURÉ, clé SRA, 1re + 2e ligne.

    Lit la table normalisée produite par `ingestion/fetch_tbportals.py`
    (strain_id, drug, drug_code, phenotype). Absente tant que le bundle DUA
    n'a pas été déposé puis normalisé.
    """
    f = paths.SOURCES / "tbportals" / "tbportals_dst_by_sra.csv"
    if not f.exists():
        print(f"[tbportals] table normalisée absente ({f}) — voir ingestion/fetch_tbportals.py")
        return pd.DataFrame(columns=COLUMNS)
    t = pd.read_csv(f)
    rows = []
    for _, r in t.iterrows():
        ph = _norm_pheno(r.get("phenotype"))
        if ph is None or not isinstance(r.get("strain_id"), str):
            continue
        canon, code = canonical_drug(r.get("drug"))
        rows.append({"strain_id": r["strain_id"], "drug": canon, "drug_code": code,
                     "phenotype": ph, "mic": "", "source": "TB_Portals"})
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[tbportals] {len(df)} mesures | {df['strain_id'].nunique()} souches")
    return df


# Priorité de confiance (la plus fiable d'abord) pour le consensus.
# `CRyPTIC_v3` passe juste après `CRyPTIC` : mêmes plaques UKMYC pour une large
# part, provenance et qualité de mesure documentées cohorte par cohorte, et
# curation par l'équipe qui produit la donnée. Il passe donc devant
# TBannotator_biosample (antibiogramme de BioSample, provenance non documentée),
# devant PATRIC (hétérogène, cf. P8.5) et devant souches_datees (étiquettes
# calculées sur le génotype, donc circulaires, cf. P5.2 et P8.4).
PRIORITY = ["CRyPTIC", "CRyPTIC_v3", "TB_Portals", "TBannotator_biosample",
            "Coll2018", "souches_datees", "PATRIC", "literature"]


def calculer_consensus(phe):
    """Dédoublonnage souche×drogue + détection de conflits inter-sources.

    Une même paire (souche, drogue) peut apparaître dans plusieurs sources ; le consensus garde le
    phénotype de la source la plus prioritaire. Rend (consensus, n_paires, n_conflits)."""
    rank = {s: i for i, s in enumerate(PRIORITY)}
    phe = phe.assign(_rank=phe["source"].map(lambda s: rank.get(s, 99)))
    pairs = phe.groupby(["strain_id", "drug"])["phenotype"].nunique()
    consensus = (phe.sort_values("_rank")
                    .drop_duplicates(["strain_id", "drug"], keep="first")
                    .drop(columns="_rank"))
    return consensus, len(pairs), int((pairs > 1).sum())


def toutes_sources(patric=None):
    """Les huit ingestions, dans l'ordre de `main` ; `patric` impose le fichier PATRIC lu."""
    parts = [ingest_cryptic(), ingest_cryptic_v3(), ingest_patric(patric),
             ingest_literature(), ingest_coll(), ingest_souches_datees(),
             ingest_tbannotator_dst(), ingest_tbportals()]
    return pd.concat([p for p in parts if len(p)], ignore_index=True)[COLUMNS]


def main():
    paths.ensure_dirs()
    phe = toutes_sources()
    phe.to_csv(paths.PHENOTYPES_TSV, sep="\t", index=False)

    print(f"\n=== TABLE PHÉNOTYPES (mesures brutes) -> {paths.PHENOTYPES_TSV} ===")
    print(f"{len(phe)} mesures | {phe['strain_id'].nunique()} souches uniques | "
          f"{phe['drug'].nunique()} drogues")
    print("sources:", phe["source"].value_counts().to_dict())

    consensus, n_pairs, n_conflict = calculer_consensus(phe)
    consensus_path = paths.CATALOGUE_DIR / "phenotypes_souches_consensus.tsv"
    consensus.to_csv(consensus_path, sep="\t", index=False)

    print(f"\n=== CONSENSUS dédoublonné -> {consensus_path} ===")
    print(f"{n_pairs} paires (souche, drogue) uniques | "
          f"{consensus['strain_id'].nunique()} souches")
    print(f"conflits inter-sources (R vs S sur la même paire) : {n_conflict} "
          f"({100*n_conflict/max(n_pairs,1):.2f} %)")
    print("\nmesures par drogue (consensus, top) :")
    print(consensus["drug"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
