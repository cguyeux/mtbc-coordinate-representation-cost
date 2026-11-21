#!/usr/bin/env python3
"""Fetch/normalise : CRyPTIC élargi (« cryptic-dataset »), DST MESURÉ, clé SRA.

Objet       : rendre disponible, au format souche→phénotype du projet, le jeu
              CRyPTIC élargi (65 842 échantillons à DST, 53 897 avec génome
              traité) qui remplace la table de 12 287 isolats déjà ingérée.
Entrées     : tables parquet du dépôt Zenodo 22128395 (v3.4.1), téléchargées
              dans `data/sources/cryptic/zenodo_v3.4.1/` ; `analyses/vocab.py`
              pour les noms de médicaments ; `catalogue/phenotypes_souches.tsv`
              pour le seul mode `--rapport`.
Sorties     : `data/sources/cryptic/cryptic_v3_dst_by_sra.csv` (source
              CANDIDATE, non branchée dans phase2) et son
              `cryptic_v3_provenance.json` (DOI, version, sha256, effectifs).
Réutilisable: oui pour tout projet MTBC ayant besoin de DST mesurés appariés à
              des accessions ENA ; la normalisation ne dépend d'aucun objet
              propre à ce projet hors `vocab.py`.
Projet      : Resistance_antibio, piste P1.
Date        : 2026-11-11.

Le jeu CRyPTIC publié en 2022 (`CRyPTIC_reuse_table_20240917.csv`, 12 287
isolats, déjà ingéré par `phase2`) a été remplacé par un jeu bien plus large,
déposé sur Zenodo sous licence CC-BY-4.0 :

    DOI concept  : 10.5281/zenodo.15679730   (toutes versions)
    version lue  : v3.4.1, 12 septembre 2026, enregistrement 22128395
    65 842 échantillons à DST, dont 53 897 avec un génome traité.

Ce script télécharge les seules tables NÉCESSAIRES au volet phénotype (une
vingtaine de Mo ; le dépôt complet pèse 4,5 Go à cause des VCF et des FASTA,
qui n'ont pas leur place sur le portable) et les normalise en une table
    data/sources/cryptic/cryptic_v3_dst_by_sra.csv
        strain_id, drug, drug_code, phenotype, quality, dst_source, cohort

Deux colonnes sont conservées EN PLUS du format des autres fetchers, et ce
n'est pas un ornement :

  * `dst_source` — la cohorte d'où vient la MESURE de DST (CRyPTIC, NEJM2018,
    WHO2019, KMVALIDATION2025, SEQTREAT2020, CLIRES2020…). Le jeu est un
    agrégat de douze cohortes, dont plusieurs ont servi à construire le
    catalogue OMS que ce projet utilise comme baseline. Sans cette colonne, la
    baseline serait évaluée sur ses propres isolats d'entraînement sans qu'on
    puisse le voir. Règle du projet (P8.5) : la provenance, jamais l'accord
    avec le catalogue.
  * `cohort` — `CRyPTIC-v1.0` (jeu remis à FIND/Seq&Treat pour construire la
    1re édition du catalogue OMS), `CRyPTIC-v2.0` (ajouts postérieurs au gel
    de 2020), `CRyPTIC-v3.0` (jeu de validation déclaré). C'est le seul
    découpage entraînement/validation documenté par les auteurs.

Ce script n'écrit RIEN dans `catalogue/` : il produit la table de source. Elle
est lue par `phase2.ingest_cryptic_v3`, branché le 2026-11-11 sur décision de CG
(piste P1.3) ; `update_database.py` appelle ce fetcher avant la reconstruction.
Attention, le branchement a fait passer la cible de 52 496 à 88 017 souches, si
bien que la matrice de features, `baseline_audit.tsv` et les chiffres de
`article/main.tex` décrivent encore l'ancienne cible (cascade P1.5).

Usage
    python3 analyses/ingestion/fetch_cryptic_v3.py            # télécharge + normalise
    python3 analyses/ingestion/fetch_cryptic_v3.py --no-fetch # normalise l'existant
    python3 analyses/ingestion/fetch_cryptic_v3.py --rapport  # + gain net mesuré
"""
import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # racine projet
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # analyses/
sys.path.insert(0, str(Path(__file__).resolve().parent))      # ingestion/
import manifest
import paths
from vocab import canonical_drug

ZENODO_RECORD = "22128395"
ZENODO_VERSION = "v3.4.1"
ZENODO_DOI = "10.5281/zenodo.22128395"
ZENODO_CONCEPT_DOI = "10.5281/zenodo.15679730"

RAW = paths.SRC_CRYPTIC / f"zenodo_{ZENODO_VERSION}"
OUT = paths.SRC_CRYPTIC / "cryptic_v3_dst_by_sra.csv"
PROV = paths.SRC_CRYPTIC / "cryptic_v3_provenance.json"

# tables indispensables au volet phénotype ; les lourdes (MUTATIONS 1,06 Go,
# VARIANTS 1,29 Go, VCF 1,58 Go, FASTA 341 Mo) sont volontairement exclues :
# elles relèvent du mésocentre, pas du portable.
TABLES = ("DST_MEASUREMENTS.parquet", "DST_SAMPLES.parquet",
          "WGS_SAMPLES.parquet", "GENOMES.parquet",
          "DRUG_CODES.csv.gz", "RELEASE_NOTES.md")
OPTIONNELLES = ("PREDICTIONS.parquet", "EFFECTS.parquet")

RANG_QUALITE = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}


def _telecharger(nom: str) -> bool:
    """Télécharge une table si elle manque. Rend True si le fichier est là."""
    cible = RAW / nom
    if cible.exists() and cible.stat().st_size > 0:
        return True
    url = (f"https://zenodo.org/records/{ZENODO_RECORD}/files/"
           f"{nom}?download=1")
    tmp = cible.with_suffix(cible.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as f:
            while True:
                bloc = r.read(1 << 20)
                if not bloc:
                    break
                f.write(bloc)
        tmp.replace(cible)
        print(f"[cryptic_v3] téléchargé {nom} "
              f"({cible.stat().st_size / 1e6:.1f} Mo)")
        return True
    except Exception as exc:                                    # noqa: BLE001
        tmp.unlink(missing_ok=True)
        print(f"[cryptic_v3] échec du téléchargement de {nom} : {exc}")
        return False


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for bloc in iter(lambda: f.read(1 << 20), b""):
            h.update(bloc)
    return h.hexdigest()


def _lire(nom: str, colonnes):
    import pyarrow.parquet as pq
    return pq.read_table(RAW / nom, columns=colonnes).to_pylist()


def normaliser():
    """Rend (lignes, stats) : une ligne par mesure R/S rattachée à un run ENA."""
    wgs = _lire("WGS_SAMPLES.parquet",
                ["UNIQUEID", "run_accession", "dataset", "status",
                 "has_mutations", "has_multiple_ena_run_accessions"])
    gen = _lire("GENOMES.parquet",
                ["UNIQUEID", "SPECIES", "LINEAGE", "TB_COVERAGE", "TB_DEPTH"])
    dstm = _lire("DST_MEASUREMENTS.parquet",
                 ["UNIQUEID", "DRUG", "PHENOTYPE", "QUALITY", "SOURCE"])

    u2run, u2cohorte = {}, {}
    for r in wgs:
        acc = (r.get("run_accession") or "").strip()
        if acc:
            u2run[r["UNIQUEID"]] = acc
            u2cohorte[r["UNIQUEID"]] = r.get("dataset")
    u2gen = {r["UNIQUEID"]: r for r in gen}

    lignes, inconnus, sans_run, non_rs = [], {}, 0, 0
    for r in dstm:
        ph = r.get("PHENOTYPE")
        if ph not in ("R", "S"):          # U (inconnu), I (intermédiaire), vide
            non_rs += 1
            continue
        acc = u2run.get(r["UNIQUEID"])
        if not acc:                       # DST sans génome apparié
            sans_run += 1
            continue
        canon, code = canonical_drug(r.get("DRUG"))
        if code == "?":
            inconnus[r.get("DRUG")] = inconnus.get(r.get("DRUG"), 0) + 1
            continue
        g = u2gen.get(r["UNIQUEID"], {})
        lignes.append({
            "strain_id": acc,
            "drug": canon,
            "drug_code": code,
            "phenotype": ph,
            "quality": r.get("QUALITY") or "",
            "dst_source": r.get("SOURCE") or "",
            "cohort": u2cohorte.get(r["UNIQUEID"]) or "",
            "species": g.get("SPECIES") or "",
            "lineage": g.get("LINEAGE") or "",
        })

    # dédoublonnage exact (même mesure répétée dans la table amont)
    vus, uniques = set(), []
    for d in lignes:
        cle = (d["strain_id"], d["drug_code"], d["phenotype"],
               d["quality"], d["dst_source"])
        if cle in vus:
            continue
        vus.add(cle)
        uniques.append(d)

    stats = {
        "mesures_rs": len(uniques),
        "mesures_non_rs_ecartees": non_rs,
        "mesures_sans_run_ena": sans_run,
        "codes_medicament_inconnus": inconnus,
        "souches": len({d["strain_id"] for d in uniques}),
        "paires": len({(d["strain_id"], d["drug_code"]) for d in uniques}),
    }
    return uniques, stats


def ecrire(lignes, stats, fichiers):
    import csv
    champs = ["strain_id", "drug", "drug_code", "phenotype", "quality",
              "dst_source", "cohort", "species", "lineage"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=champs)
        w.writeheader()
        w.writerows(lignes)
    prov = {
        "source": "CRyPTIC Consortium Dataset (cryptic-dataset)",
        "zenodo_record": ZENODO_RECORD,
        "version": ZENODO_VERSION,
        "doi": ZENODO_DOI,
        "concept_doi": ZENODO_CONCEPT_DOI,
        "licence": "CC-BY-4.0",
        "lu_le": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "fichiers_amont": fichiers,
        "stats": stats,
    }
    PROV.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n")
    return prov


def rapport(lignes):
    """Gain net face à catalogue/phenotypes_souches.tsv, par médicament."""
    import collections
    import csv
    cible = paths.CATALOGUE_DIR / "phenotypes_souches.tsv"
    if not cible.exists():
        # Pas un garde silencieux : sans la cible, « gain net » n'a pas de
        # sens et un rapport vide passerait pour un gain nul.
        raise SystemExit(f"[cryptic_v3] cible absente : {cible} — "
                         f"régénérer par analyses/phase2_build_phenotypes.py "
                         f"avant de demander --rapport")
    projet, souches_projet = set(), set()
    with open(cible) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            projet.add((r["strain_id"], r["drug_code"]))
            souches_projet.add(r["strain_id"])

    paires = collections.defaultdict(list)
    for d in lignes:
        paires[(d["strain_id"], d["drug_code"])].append(d)
    nouvelles = {k for k in paires if k not in projet}
    souches = {s for s, _ in paires}

    print(f"\n[cryptic_v3] GAIN NET face à {cible.name}")
    print(f"  paires du projet, toutes sources : {len(projet):8d}")
    print(f"  paires du jeu élargi             : {len(paires):8d}")
    print(f"  paires NOUVELLES                 : {len(nouvelles):8d} "
          f"(+{100 * len(nouvelles) / len(projet):.1f} %)")
    print(f"  souches nouvelles                : "
          f"{len(souches - souches_projet):8d}")

    par = collections.defaultdict(lambda: [0, 0, 0])
    pj = collections.Counter(d for _, d in projet)
    for k, v in paires.items():
        c = par[k[1]]
        c[0] += 1
        if k not in projet:
            c[1] += 1
            if any(x["phenotype"] == "R" for x in v):
                c[2] += 1
    print(f"\n  {'code':6s} {'élargi':>8s} {'NOUV':>8s} {'NOUV_R':>8s} "
          f"{'projet':>8s} {'gain %':>8s}")
    for d in sorted(par, key=lambda x: -par[x][1]):
        e, n, nr = par[d]
        p = pj.get(d, 0)
        g = f"{100 * n / p:.0f}" if p else "neuf"
        print(f"  {d:6s} {e:8d} {n:8d} {nr:8d} {p:8d} {g:>8s}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-fetch", action="store_true",
                    help="ne rien télécharger, normaliser ce qui est sur disque")
    ap.add_argument("--rapport", action="store_true",
                    help="chiffrer le gain net face à phenotypes_souches.tsv")
    ap.add_argument("--avec-predictions", action="store_true",
                    help="télécharger aussi PREDICTIONS/EFFECTS (catalogue OMS "
                         "appliqué par gnomonicus : comparateur externe)")
    a = ap.parse_args()

    paths.ensure_dirs()
    RAW.mkdir(parents=True, exist_ok=True)

    voulues = list(TABLES) + (list(OPTIONNELLES) if a.avec_predictions else [])
    if a.no_fetch:
        manquantes = [n for n in TABLES if not (RAW / n).exists()]
        if manquantes:
            print(f"[cryptic_v3] tables absentes dans {RAW} : "
                  f"{', '.join(manquantes)}")
            print(f"  → relancer sans --no-fetch, ou télécharger depuis "
                  f"https://doi.org/{ZENODO_DOI}")
            return 1
    else:
        for nom in voulues:
            if not _telecharger(nom) and nom in TABLES:
                print(f"[cryptic_v3] table indispensable manquante : {nom}")
                return 1

    lignes, stats = normaliser()
    # On trace TOUTES les tables présentes, optionnelles comprises : un contrôle
    # qui s'appuie sur PREDICTIONS (comparateur externe, P1.1) doit pouvoir se
    # rejouer depuis la provenance seule.
    fichiers = {n: {"octets": (RAW / n).stat().st_size,
                    "sha256": _sha256(RAW / n)}
                for n in (*TABLES, *OPTIONNELLES) if (RAW / n).exists()}
    ecrire(lignes, stats, fichiers)

    print(f"[cryptic_v3] {stats['mesures_rs']} mesures R/S | "
          f"{stats['souches']} souches | {stats['paires']} paires -> {OUT}")
    print(f"  écartées : {stats['mesures_non_rs_ecartees']} non-R/S, "
          f"{stats['mesures_sans_run_ena']} sans run ENA")
    if stats["codes_medicament_inconnus"]:
        print(f"  codes médicament inconnus de vocab.py : "
              f"{stats['codes_medicament_inconnus']}")
    manifest.record("cryptic_v3", stats["mesures_rs"], OUT,
                    note=f"CRyPTIC élargi {ZENODO_VERSION}, DST mesuré, "
                         f"{stats['souches']} souches ; lue par "
                         f"phase2.ingest_cryptic_v3 (branchée le 2026-11-11, P1.3)")
    if a.rapport:
        rapport(lignes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
