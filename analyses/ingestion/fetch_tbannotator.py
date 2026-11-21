#!/usr/bin/env python3
"""Fetch : snapshot de l'antibiogramme MESURÉ de TBannotator (DST, table BioSample).

Interroge le serveur MCP TBannotator (HTTP streamable) et écrit
`data/sources/tbannotator/antibiogram_snapshot.csv` :
    strain_name, ncbi_biosample, antibiogram_inh, antibiogram_rif,
    antibiogram_pza, antibiogram_emb
pour les souches M. tuberculosis ayant au moins un antibiogramme renseigné.

Phénotype MESURÉ (pas une prédiction du catalogue) — vérifié, cf. SOURCES_CANDIDATES.md.
Reproductible : relancer ce script rafraîchit le snapshot et journalise le run.

À exécuter avec le venv d'ingestion :  .venv-ingest/bin/python analyses/ingestion/fetch_tbannotator.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths
from tbannotator_client import query_df
import manifest

SQL = """
SELECT strain_name, ncbi_biosample,
       antibiogram_inh, antibiogram_rif, antibiogram_pza, antibiogram_emb
FROM mv_strain_metadata
WHERE species_group = 'M. tuberculosis'
  AND (antibiogram_inh IS NOT NULL OR antibiogram_rif IS NOT NULL
       OR antibiogram_pza IS NOT NULL OR antibiogram_emb IS NOT NULL)
""".strip()

OUT = paths.SRC_TBANNOTATOR / "antibiogram_snapshot.csv"


def main():
    paths.ensure_dirs()
    print("[tbann] requête MCP en cours…")
    df = query_df(SQL, max_rows=50000, timeout_s=120)
    df.to_csv(OUT, index=False)
    manifest.record("tbannotator_antibiogram", len(df), OUT,
                    note="DST mesuré BioSample, 4 drogues 1re ligne")
    print(f"[tbann] {len(df)} souches -> {OUT}")


if __name__ == "__main__":
    main()
