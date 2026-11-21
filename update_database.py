#!/usr/bin/env python3
"""Orchestrateur de mise à jour de la base mutation/souche ↔ résistance.

Enchaîne, de façon idempotente et journalisée :
  1. Rafraîchissement des sources (fetchers best-effort, chacun peut échouer sans
     bloquer les autres) :
       - antibiogramme TBannotator (MCP)        → data/sources/tbannotator/
       - TB Portals (normalisation d'un bundle) → data/sources/tbportals/
       (WHO, CRyPTIC, PATRIC, Coll, souches datées sont des fichiers statiques
        déjà sur disque ; aucun fetch réseau nécessaire.)
  2. Reconstruction du catalogue mutation→résistance (phase1).
  3. Reconstruction de la table phénotypes + consensus (phase2).

Chaque run est tracé dans `data/sources/INGESTION_MANIFEST.tsv`.

Usage (depuis le venv d'ingestion qui contient mcp + pandas + openpyxl) :
    .venv-ingest/bin/python update_database.py
    .venv-ingest/bin/python update_database.py --no-fetch   # rebuild seul, sans réseau

Pour une mise à jour régulière : planifier cette commande (cron / systemd timer).
"""
import argparse
import runpy
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ANA = ROOT / "analyses"


def _run_module(path: Path, label: str, argv=()):
    print(f"\n{'='*60}\n>>> {label}\n{'='*60}")
    sys.argv = [str(path), *argv]      # le module lit ses propres options, pas celles de ce lanceur
    try:
        runpy.run_path(str(path), run_name="__main__")
    except SystemExit:
        pass
    except Exception as e:  # best-effort : on continue le pipeline
        print(f"[!] {label} a échoué : {type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fetch", action="store_true",
                    help="ne pas rafraîchir les sources réseau, reconstruire seulement")
    args = ap.parse_args()

    if not args.no_fetch:
        _run_module(ANA / "ingestion" / "fetch_tbannotator.py", "FETCH antibiogramme TBannotator (MCP)")
        _run_module(ANA / "ingestion" / "fetch_tbportals.py", "FETCH/normalise TB Portals (bundle DUA)")
        _run_module(ANA / "ingestion" / "fetch_cryptic_v3.py", "FETCH/normalise CRyPTIC élargi (Zenodo)")

    _run_module(ANA / "phase1_build_catalogue.py", "BUILD catalogue mutation→résistance (v1, tb-profiler 2022)")
    # phase1 reconstruit le catalogue v1 : sans cette étape il écraserait le catalogue courant (v2, tbdb 2026, AA1)
    _run_module(ANA / "phase1b_tbprofiler2026_v2.py", "BUILD catalogue v2 (tbdb 2026 + colonne regle) et bascule",
                argv=["--bascule"])
    _run_module(ANA / "phase2_build_phenotypes.py", "BUILD table phénotypes + consensus")
    print("\n=== mise à jour terminée ===")


if __name__ == "__main__":
    main()
