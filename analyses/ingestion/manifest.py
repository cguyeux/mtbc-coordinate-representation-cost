"""Manifeste d'ingestion : journalise chaque rafraîchissement de source.

Fichier append-only `data/sources/INGESTION_MANIFEST.tsv` :
    timestamp_utc  source  rows  path  note
Permet de savoir quand chaque source a été mise à jour pour la dernière fois et
combien de lignes elle contenait, sans relire les snapshots eux-mêmes.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import paths

MANIFEST = paths.SOURCES / "INGESTION_MANIFEST.tsv"


def record(source: str, rows: int, path, note: str = ""):
    paths.ensure_dirs()
    new = not MANIFEST.exists()
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(MANIFEST, "a") as f:
        if new:
            f.write("timestamp_utc\tsource\trows\tpath\tnote\n")
        f.write(f"{ts}\t{source}\t{rows}\t{path}\t{note}\n")
    print(f"[manifest] {source}: {rows} lignes @ {ts}")
