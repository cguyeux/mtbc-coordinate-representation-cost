#!/usr/bin/env python3
"""Fetch : lignées Coll de TBannotator (pour la stratification CV par lignée).

Récupère strain_name -> lignée Coll (code complet + niveau 1) pour toutes les souches
classées, via le client MCP scriptable. Sert à grouper les souches par lignée lors de
la validation croisée du modèle (GroupKFold), afin de contrôler le biais phylogénétique.

Sortie : data/sources/tbannotator/lineage_coll_snapshot.csv
À exécuter avec le venv :  .venv-ingest/bin/python analyses/ingestion/fetch_tbannotator_lineage.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
from tbannotator_client import query_df
import manifest

SQL = """
SELECT strain_name, lineage_code, lineage_level_1
FROM mv_strain_classification
WHERE system_name = 'Coll' AND strain_name IS NOT NULL
""".strip()

OUT = paths.SRC_TBANNOTATOR / "lineage_coll_snapshot.csv"


def main():
    paths.ensure_dirs()
    print("[lineage] requête MCP (classification Coll)…")
    df = query_df(SQL, max_rows=50000, timeout_s=120)
    df = df.drop_duplicates("strain_name")
    df.to_csv(OUT, index=False)
    manifest.record("tbannotator_lineage_coll", len(df), OUT, note="strain->lignée Coll, pour GroupKFold")
    print(f"[lineage] {len(df)} souches -> {OUT}")


if __name__ == "__main__":
    main()
