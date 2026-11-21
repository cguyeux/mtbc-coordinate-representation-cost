"""Chemins centralisés du projet Resistance_antibio.

Toute la convention de répertoires est ancrée sur la racine du projet (le
dossier qui contient ce fichier), pas sur le répertoire de travail courant.
Les scripts et notebooks doivent importer ces constantes plutôt que de
hardcoder des chemins relatifs, afin de rester exécutables depuis n'importe où.

    from paths import CATALOGUE_DIR, SRC_CRYPTIC, PANGENOME, ...
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# --- Données ---------------------------------------------------------------
DATA = ROOT / "data"
FASTAS = DATA / "fastas"                      # références génomiques + index BLAST/BWA
LINEAGES = DATA / "lineages"                  # génotypes de lignées publiés
METADATA = DATA / "metadata"                  # dates / lieux des souches

# Sources tierces brutes (read-only, ne pas modifier en place)
SOURCES = DATA / "sources"
SRC_WHO = SOURCES / "who_catalogue"           # catalogue OMS (WHO mutations catalogue)
SRC_CRYPTIC = SOURCES / "cryptic"             # CRyPTIC reuse table
SRC_PATRIC = SOURCES / "patric"               # PATRIC genomes AMR
SRC_TBPROFILER = SOURCES / "tbprofiler"       # tb-profiler resistance DB
SRC_TBANNOTATOR = SOURCES / "tbannotator"     # snapshot antibiogramme TBannotator (DST mesuré)
SRC_LITERATURE = SOURCES / "literature"       # corpus d'articles (chunks .md)

# Pangénome (pickles volumineux, non versionnés)
PANGENOME = DATA / "pangenome"
PAN_SPDI = PANGENOME / "pan_spdi.pkl"         # SPDI -> [souches]
PAN_STRAINS = PANGENOME / "pan_strains.pkl"   # souche -> {SPDI}

# Référence H37Rv
H37RV_FASTA = FASTAS / "NC_000962.3.fasta"
H37RV_CDS = FASTAS / "NC_000962.3_CDS.fasta"

# --- Sorties ---------------------------------------------------------------
CATALOGUE_DIR = ROOT / "catalogue"            # LIVRABLE: catalogue mutation->résistance
CATALOGUE_TSV = CATALOGUE_DIR / "catalogue_consolide.tsv"
PHENOTYPES_TSV = CATALOGUE_DIR / "phenotypes_souches.tsv"  # souche x drogue (long)

RESULTATS = ROOT / "results"                   # sorties dérivées, caches, figures brutes
EXPERIMENTS = ROOT / "experiments"
ARTICLE = ROOT / "article"
LITREVIEW = ROOT / "litterature_review"


def ensure_dirs():
    """Crée les répertoires de sortie s'ils n'existent pas."""
    for d in (CATALOGUE_DIR, RESULTATS, EXPERIMENTS, METADATA,
              SRC_WHO, SRC_CRYPTIC, SRC_PATRIC, SRC_TBPROFILER, SRC_TBANNOTATOR,
              SRC_LITERATURE, PANGENOME):
        d.mkdir(parents=True, exist_ok=True)
