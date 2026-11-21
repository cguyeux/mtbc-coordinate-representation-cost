#!/usr/bin/env python3
"""Phase 1 — Construction du catalogue consolidé mutation -> résistance.

Fusionne les sources de savoir *au niveau mutation* (une mutation, une drogue,
un verdict) en une table longue unique, en conservant la traçabilité de chaque
assertion (colonne `source`).

Deux familles de sources contribuent ici :
  - Savoir expert curé   : OMS (WHO catalogue) et tb-profiler resistance DB.
  - Signal empirique brut : différentiel R-S miné sur CRyPTIC (résultats/resistance_ami.csv).

Les sources qui donnent souche -> phénotype (CRyPTIC complet, PATRIC, littérature)
ne sont PAS ici : elles peuplent phase2 (table phénotypes), parce qu'il faut les
croiser avec les mutations de chaque souche pour produire des couples mutation/résistance.

Schéma de sortie (catalogue/catalogue_consolide.tsv, TSV) :
  spdi gene variant drug drug_code call source source_grade n_R n_S metric notes
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))        # analyses/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # racine projet
import paths
from vocab import canonical_drug
from hgvs_to_spdi import resolve as hgvs_to_spdi
from hgvs_to_spdi import resolve_all


# P7.7 : drapeau du test falsifiant du recadrage de A1. Par défaut (False) le comportement historique,
# qui a produit le catalogue v2 vivant, est conservé à l'octet près ; `--blocs` émet aussi les
# représentations en bloc (MNV) que le fichier de coordonnées de l'OMS fournit pour un variant qui
# a déjà une forme mono-base.
EMETTRE_BLOCS = False


def spdis_for(variant):
    """SPDI à émettre pour un variant : mono-base si SNV, sinon TOUTES les
    représentations indel/MNV courtes (la normalisation varie ; émettre toutes
    maximise le matching exact avec le pangénome). [""] si rien.

    Défaut historique : dès qu'une forme mono-base existe, les autres représentations du même
    variant (blocs MNV) sont écartées (`sp = mono if mono else alls`). Avec `EMETTRE_BLOCS`
    (P7.7), toutes les représentations de `resolve_all` sont émises."""
    alls = resolve_all(variant)
    if EMETTRE_BLOCS:
        return alls or [""]
    mono = [s for s in alls if (lambda p: len(p[2]) == 1 and len(p[3]) == 1)(s.split(":"))]
    sp = mono if mono else alls
    return sp or [""]

COLUMNS = ["spdi", "gene", "variant", "drug", "drug_code", "call",
           "source", "source_grade", "n_R", "n_S", "metric", "notes"]


def _call_from_grade(grade: str) -> str:
    """Normalise un grade de confiance (OMS / WHO / tb-profiler) en verdict."""
    if not isinstance(grade, str):
        return "uncertain"
    g = grade.lower()
    if "not assoc" in g or "5)" in g or "4)" in g:
        return "S-associated"
    if "uncertain" in g or "3)" in g:
        return "uncertain"
    if "assoc w r" in g or "1)" in g or "2)" in g:
        return "R-associated"
    return "uncertain"


def ingest_tbprofiler() -> pd.DataFrame:
    """tb-profiler resistance DB : Gene, Mutation (HGVS), Drug, Confers, WHO Confidence.

    Source experte la plus fiable actuellement (2311 lignes, 21 drogues).
    NB : pas de SPDI — la notation est HGVS (c./p.). La résolution HGVS->SPDI
    est un TODO (cf. SOURCES.md), nécessaire pour croiser avec le pangénome.
    """
    f = paths.SRC_TBPROFILER / "tb-profiler-db.csv"
    if not f.exists():
        print(f"[tbprofiler] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    t = pd.read_csv(f)
    rows = []
    for _, r in t.iterrows():
        canon, code = canonical_drug(r["Drug"])
        variant = f"{r['Gene']}_{r['Mutation']}"
        for spdi in spdis_for(variant):
            rows.append({
                "spdi": spdi,
                "gene": r["Gene"],
                "variant": variant,
                "drug": canon, "drug_code": code,
                "call": _call_from_grade(r.get("WHO Confidence")),
                "source": "tbprofiler",
                "source_grade": r.get("WHO Confidence", ""),
                "n_R": "", "n_S": "", "metric": "",
                "notes": str(r.get("Literature", "") or "")[:120],
            })
    df = pd.DataFrame(rows)
    print(f"[tbprofiler] {len(df)} assertions, {df['drug'].nunique()} drogues, "
          f"{(df['spdi'] != '').sum()} avec SPDI résolu")
    return df


def ingest_who_master() -> pd.DataFrame:
    """Catalogue OMS officiel 2e éd. (WHO-UCN-TB-2023.7), feuille Catalogue_master_file.

    ~48 000 lignes (variant x drogue), 15 drogues, grade de confiance WHO sur 5 niveaux.
    Le SPDI est résolu depuis le variant HGVS via le fichier officiel de coordonnées
    génomiques (hgvs_to_spdi). C'est la source experte de référence.
    """
    import glob
    cand = glob.glob(str(paths.SRC_WHO / "*master*.xlsx")) or \
        glob.glob(str(paths.SRC_WHO / "WHO-UCN-TB-2023*.xlsx"))
    if not cand:
        print("[who] master file absent — télécharger WHO-UCN-TB-2023.7-eng.xlsx")
        return pd.DataFrame(columns=COLUMNS)
    m = pd.read_excel(cand[0], sheet_name="Catalogue_master_file", header=2)
    rows = []
    for _, r in m.iterrows():
        variant = r.get("variant")
        if not isinstance(variant, str):
            continue
        canon, code = canonical_drug(r.get("drug"))
        grade = r.get("FINAL CONFIDENCE GRADING")
        for spdi in spdis_for(variant):
            rows.append({
                "spdi": spdi,
                "gene": r.get("gene", ""),
                "variant": variant,
                "drug": canon, "drug_code": code,
                "call": _call_from_grade(grade),
                "source": "WHO_catalogue",
                "source_grade": grade if isinstance(grade, str) else "",
                "n_R": r.get("Present_R", ""),
                "n_S": r.get("Present_S", ""),
                "metric": f"PPV={r.get('PPV','')}",
                "notes": f"tier={r.get('tier','')};effect={r.get('effect','')}",
            })
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[who] {len(df)} assertions, {df['drug'].nunique()} drogues, "
          f"{(df['spdi'] != '').sum()} avec SPDI résolu")
    return df


def ingest_cryptic_empirical() -> pd.DataFrame:
    """Signal empirique : différentiel R-S miné sur CRyPTIC (résultats/resistance_ami.csv).

    Amikacin uniquement (prototype). Apporte spdi + gene + comptes, mais NON curé :
    confondu par le biais MDR/phylogéographique (cf. cahier de labo). Marqué
    call='empirical-enriched' et source='CRyPTIC_empirical' pour ne pas le confondre
    avec un verdict expert.
    """
    f = paths.RESULTATS / "resistance_ami.csv"
    if not f.exists():
        print(f"[cryptic-emp] absent: {f}")
        return pd.DataFrame(columns=COLUMNS)
    e = pd.read_csv(f)
    rows = []
    for _, r in e.iterrows():
        rows.append({
            "spdi": r.get("spdi", ""),
            "gene": r.get("gene_name", ""),
            "variant": "",
            "drug": "amikacin", "drug_code": "AMK",
            "call": "empirical-enriched",
            "source": "CRyPTIC_empirical",
            "source_grade": "",
            "n_R": r.get("resistant", ""),
            "n_S": r.get("susceptible", ""),
            "metric": f"diff={r.get('diff','')}",
            "notes": f"TBann_total={r.get('total_TBannotator','')}",
        })
    df = pd.DataFrame(rows, columns=COLUMNS)
    print(f"[cryptic-emp] {len(df)} SNPs enrichis (AMK, NON curé)")
    return df


def main():
    global EMETTRE_BLOCS
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blocs", action="store_true",
                    help="P7.7 : émettre aussi les représentations en bloc du fichier de coordonnées OMS "
                         "(test falsifiant du recadrage de A1) ; défaut : comportement historique")
    ap.add_argument("--out", help="fichier de sortie (défaut : le catalogue vivant, paths.CATALOGUE_TSV)")
    args = ap.parse_args()
    EMETTRE_BLOCS = args.blocs
    sortie = Path(args.out) if args.out else paths.CATALOGUE_TSV
    sortie.parent.mkdir(parents=True, exist_ok=True)
    paths.ensure_dirs()
    parts = [ingest_who_master(), ingest_tbprofiler(), ingest_cryptic_empirical()]
    cat = pd.concat([p for p in parts if len(p)], ignore_index=True)[COLUMNS]

    # propagation au niveau acide aminé : ajoute les codons synonymes des
    # déterminants R connus (comble l'angle mort SPMI-indexé, cf. cahier 2026-05-31)
    from aa_propagation import propagate_catalogue
    n0 = len(cat)
    cat = propagate_catalogue(cat)[COLUMNS]
    print(f"[AA-prop] +{len(cat) - n0} SPDI codons-synonymes de déterminants R")

    cat = cat.sort_values(["drug", "source", "gene"]).reset_index(drop=True)
    cat.to_csv(sortie, sep="\t", index=False)

    print(f"\n=== CATALOGUE CONSOLIDÉ -> {sortie} ===")
    print(f"{len(cat)} assertions | {cat['drug'].nunique()} drogues | "
          f"sources: {cat['source'].value_counts().to_dict()}")
    print("\nverdicts par source:")
    print(pd.crosstab(cat["source"], cat["call"]).to_string())
    print("\ncouverture SPDI (mutations croisables avec le pangénome):")
    print(f"  {(cat['spdi'] != '').sum()}/{len(cat)} assertions ont un SPDI résolu")


if __name__ == "__main__":
    main()
