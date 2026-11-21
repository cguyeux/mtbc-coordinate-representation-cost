# Cost of keeping one coordinate representation per variant when ingesting a drug-resistance catalogue of the *Mycobacterium tuberculosis* complex

Ingestion chain, catalogues, measurement scripts and measurement outputs accompanying the
manuscript of the same title (Christophe Guyeux, Université Marie et Louis Pasteur, CNRS,
institut FEMTO-ST, Besançon, France).

## What this repository answers

Catalogues of resistance-conferring mutations of the *M. tuberculosis* complex are published per
variant, in gene-based HGVS notation. A pipeline that compares genomes with such a catalogue must
first express each entry as a genomic coordinate string, and the same change can be written as
several such strings. An ingestion that keeps only one of them per variant makes a silent choice:
an isolate whose call does not use the representation that was kept is counted as carrying no
determinant, and is predicted susceptible.

This repository lets that choice be measured again. It contains two catalogues built from the same
source files and differing by one rule, the scripts that build them, and the scripts that evaluate a
deterministic baseline under each.

## Layout

| Path | Contents |
|---|---|
| `analyses/` | the ingestion chain and the measurement scripts; each file carries a header stating its object, inputs, outputs and date |
| `analyses/hgvs_to_spdi.py` | the HGVS to coordinate-string resolver, applying the one-base offset of the WHO coordinate file |
| `analyses/phase1_build_catalogue.py`, `phase1b_tbprofiler2026_v2.py` | build of the consolidated catalogue; `--blocs` emits every representation (catalogue C1) instead of one per variant (C0) |
| `analyses/variant_match.py` | the reconciled matching: exact single-base, decomposition of an observed block, exact indel |
| `analyses/p77_baseline_blocs.py` | deterministic baseline and contribution of each matching step, by removal, for one catalogue |
| `analyses/p716_intervalles.py` | paired resampling of the two arms, giving the interval of the ratio of two differences |
| `analyses/p716_gyra94.py` | full enumeration of the *gyrA* codon 94 entries of the WHO coordinate file |
| `analyses/p716_chaines_distinctes.py` | assertion rows against distinct coordinate strings, per type |
| `analyses/test_p78_chiffres_manuscrit.py` | every figure printed in the manuscript, pinned to the output file that produces it |
| `catalogue/` | the two catalogues, and `SOURCES.md` on their provenance |
| `results/` | the measurement outputs the manuscript cites |
| `article/` | the manuscript, its supplementary material and the three supplementary data files |

## Reproducing

The measurement needs a strain-by-variant matrix whose variant calls keep multi-nucleotide blocks.
That matrix is built from sequencing data that are not redistributed here; the accessions and the
calling parameters are given in the manuscript and its Supplementary Note S3.

```bash
python3 -m venv .venv-ingest && .venv-ingest/bin/pip install -r requirements-ingest.txt
.venv-ingest/bin/python update_database.py --no-fetch      # rebuild the catalogue offline
python3 analyses/p716_gyra94.py                            # the codon 94 enumeration
python3 analyses/p716_chaines_distinctes.py                # rows against distinct strings
python3 -m pytest -q analyses/test_p78_chiffres_manuscrit.py
```

The last command re-checks every figure of the manuscript against the output files of this
repository. The two scripts before it need only the WHO coordinate file and the catalogue, both
addressed below.

## Third-party sources, and what is not redistributed here

The consolidated catalogue is **derived** from three published sources, and the files of those
sources are not redistributed in this repository. They must be obtained from their publishers:

| Source | Obtain from | Licence |
|---|---|---|
| WHO catalogue of mutations, 2nd edition (2023), and its genomic-coordinates file | World Health Organization, <https://www.who.int/publications/i/item/9789240082410> | CC BY-NC-SA 3.0 IGO |
| TB-Profiler reference database (tbdb) | <https://github.com/jodyphelan/tbdb> | see that repository |
| CRyPTIC compendium, and the expanded release (Zenodo 22128395) | <https://doi.org/10.5281/zenodo.22128395> | CC-BY-4.0 |

Sequencing reads are from the European Nucleotide Archive and the NCBI Sequence Read Archive under
the accessions the manuscript gives. The precomputed pangenome and the per-strain genotype files
are not redistributed: they are the research database of the group, described in the reference the
manuscript cites for it.

## Licence

Two licences apply, because the code and the derived data do not come under the same constraint.

**Code** (`analyses/`, `*.py`, `update_database.py`): MIT, see `LICENSE-CODE`.

**Derived data** (`catalogue/`, `results/`): CC BY-NC-SA 3.0 IGO, see `LICENSE-DATA.md`. These
files are an adaptation of the WHO catalogue of mutations, which WHO publishes under CC BY-NC-SA
3.0 IGO; that licence requires an adaptation to carry the same or an equivalent licence, so the
derived data here are **not** available for commercial use. This adaptation was not created by WHO
and WHO does not endorse it. Original WHO citation: *Catalogue of mutations in Mycobacterium
tuberculosis complex and their association with drug resistance, second edition. Geneva: World
Health Organization; 2023.*

## Citing

Please cite the manuscript. If you use the catalogues or the measurement outputs, cite in addition
the three third-party sources above, and the WHO catalogue in the form WHO gives.
