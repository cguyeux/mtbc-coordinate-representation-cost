# Licence of the derived data

This file covers the data files of this repository, namely the contents of `catalogue/` and of
`results/`, and the three supplementary data files of `article/supplementary_materials/`. The code
is covered by `LICENSE-CODE` (MIT) instead.

## CC BY-NC-SA 3.0 IGO, because the upstream source requires it

The consolidated catalogues `catalogue/catalogue_C0_one_representation.tsv` and
`catalogue/catalogue_C1_every_representation.tsv` are an **adaptation** of the *Catalogue of
mutations in Mycobacterium tuberculosis complex and their association with drug resistance, second
edition* (World Health Organization, 2023), which supplies 58 994 of the 64 052 assertion rows of
the first file, and of the genomic-coordinates file distributed with it.

The World Health Organization publishes that catalogue under the Creative Commons
Attribution-NonCommercial-ShareAlike 3.0 IGO licence (CC BY-NC-SA 3.0 IGO,
<https://creativecommons.org/licenses/by-nc-sa/3.0/igo/>), whose ShareAlike term requires any
adaptation to be licensed under the same or an equivalent licence. The derived data of this
repository are therefore released under **CC BY-NC-SA 3.0 IGO**, and in particular they are **not
available for commercial use**.

## Required notices

**Attribution of the original.** *Catalogue of mutations in Mycobacterium tuberculosis complex and
their association with drug resistance, second edition.* Geneva: World Health Organization; 2023.

**Disclaimer of the adaptation**, in the form the WHO licence prescribes: this is an adaptation of
an original work of the World Health Organization. Views and opinions expressed in the adaptation
are the sole responsibility of the author of the adaptation and are not endorsed by WHO. WHO is not
responsible for the content or accuracy of this adaptation.

**Not for clinical use.** These files are research material. They are not a diagnostic device and
must not be used to guide the treatment of a patient.

## Other upstream sources of the same files

The consolidated catalogue also derives from the reference database of TB-Profiler
(<https://github.com/jodyphelan/tbdb>, see that repository for its terms) and from an empirical
signal computed from the CRyPTIC compendium (CC-BY-4.0). The phenotype labels used by the
measurement, which are **not** redistributed here, come in addition from the expanded CRyPTIC
release (<https://doi.org/10.5281/zenodo.22128395>, CC-BY-4.0), from the NCBI BioSample antibiogram
tag, from the PATRIC antimicrobial-resistance table and from a published association panel; their
provenance is set out in `catalogue/SOURCES.md`.

The measurement outputs of `results/` are our own measurements rather than an extract of a
third-party database, but they are released under the same licence as the catalogues they were
computed from, so that the ShareAlike term is not circumvented.
