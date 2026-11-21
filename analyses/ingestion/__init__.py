"""Chaîne d'ingestion des sources de résistance (mises à jour répétables).

Chaque source expose une fonction `fetch_*` qui rafraîchit un snapshot brut daté
dans `data/sources/<src>/`, et journalise le run dans le manifeste d'ingestion.
L'orchestrateur `update_database.py` enchaîne les fetchers puis reconstruit le
catalogue (phase1) et la table phénotypes (phase2).
"""
