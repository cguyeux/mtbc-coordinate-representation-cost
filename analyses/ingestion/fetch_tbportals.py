#!/usr/bin/env python3
"""Fetch/normalise : DST clinique MESURÉ de TB Portals (NIAID), clé SRA.

TB Portals n'expose pas d'API publique anonyme : l'accès bulk passe par le depot
avec accord d'usage (https://tbportals.niaid.nih.gov/download-data). Le flux propre
et répétable est donc :

  1. Télécharger périodiquement le bundle TB Portals (compte + DUA) et déposer les
     CSV bruts dans `data/sources/tbportals/raw/` :
        - un fichier *Genomics*  (lien condition_id -> accession SRA)
        - un fichier *DST*        (condition_id -> drogue -> résultat R/S)
       (ou le fichier *Patient Cases* combiné si c'est ce que fournit l'export).
  2. Lancer ce script : il joint Genomics ⨝ DST sur `condition_id` et produit la
     table normalisée `data/sources/tbportals/tbportals_dst_by_sra.csv` :
        strain_id (SRA), drug, drug_code, phenotype (R/S)
     consommée par phase2 (`ingest_tbportals`).

Détection de colonnes souple (les exports TB Portals varient) : on repère la
colonne d'accession SRA, l'identifiant condition, la drogue et le résultat par
mots-clés, et on journalise ce qui a été trouvé.
"""
import glob
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # racine projet
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analyses/ (pour vocab)
sys.path.insert(0, str(Path(__file__).resolve().parent))       # ingestion/
import paths
from vocab import canonical_drug
import manifest

RAW = paths.SRC_PATRIC.parent / "tbportals" / "raw"
OUT = paths.SRC_PATRIC.parent / "tbportals" / "tbportals_dst_by_sra.csv"


def _find_col(cols, *keywords):
    for c in cols:
        cl = str(c).lower()
        if all(k in cl for k in keywords):
            return c
    return None


def _norm(v):
    if not isinstance(v, str):
        return None
    s = v.strip().lower()
    if s.startswith("r") or s in ("1", "true"):
        return "R"
    if s.startswith("s") or s in ("0", "false"):
        return "S"
    return None


def main():
    paths.ensure_dirs()
    RAW.mkdir(parents=True, exist_ok=True)
    geno = glob.glob(str(RAW / "*[Gg]enomic*.csv"))
    dst = glob.glob(str(RAW / "*DST*.csv")) + glob.glob(str(RAW / "*dst*.csv"))
    if not geno or not dst:
        print(f"[tbportals] bundle absent dans {RAW}")
        print("  → télécharger depuis https://tbportals.niaid.nih.gov/download-data (DUA),")
        print("    déposer les CSV *Genomics* et *DST*, puis relancer ce script.")
        return

    g = pd.read_csv(geno[0])
    sra_col = _find_col(g.columns, "sra") or _find_col(g.columns, "run")
    gid_col = _find_col(g.columns, "condition")
    if not sra_col or not gid_col:
        print(f"[tbportals] colonnes SRA/condition introuvables dans {geno[0]} "
              f"({list(g.columns)[:10]}…)")
        return
    cond2sra = dict(zip(g[gid_col], g[sra_col]))

    d = pd.read_csv(dst[0])
    did = _find_col(d.columns, "condition")
    drug_col = _find_col(d.columns, "drug")
    res_col = (_find_col(d.columns, "result") or _find_col(d.columns, "interpret")
               or _find_col(d.columns, "suscept"))
    if not all([did, drug_col, res_col]):
        print(f"[tbportals] colonnes DST introuvables dans {dst[0]} "
              f"({list(d.columns)[:10]}…)")
        return

    rows = []
    for _, r in d.iterrows():
        sra = cond2sra.get(r[did])
        ph = _norm(r[res_col])
        if not isinstance(sra, str) or ph is None:
            continue
        canon, code = canonical_drug(r[drug_col])
        rows.append({"strain_id": sra, "drug": canon, "drug_code": code,
                     "phenotype": ph})
    out = pd.DataFrame(rows).drop_duplicates()
    out.to_csv(OUT, index=False)
    manifest.record("tbportals", len(out), OUT, note="DST clinique mesuré, 1re+2e ligne")
    print(f"[tbportals] {len(out)} mesures | {out['strain_id'].nunique()} souches -> {OUT}")


if __name__ == "__main__":
    main()
