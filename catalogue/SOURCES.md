# Catalogue consolidé — provenance et schéma

Construit par `analyses/phase1_build_catalogue.py` (Table 1) et
`analyses/phase2_build_phenotypes.py` (Table 2). Régénérables à tout moment.

## Deux tables, deux natures

Le projet manipule deux objets distincts qu'il ne faut pas confondre.

**Table 1 — `catalogue_consolide.tsv`** : le savoir *au niveau mutation*. Une ligne =
une assertion « telle mutation est (ou non) associée à la résistance à telle drogue »,
avec la source qui la porte. C'est le catalogue mutation→résistance proprement dit.

**Table 2 — `phenotypes_souches.tsv`** : les observations *au niveau souche*. Une ligne =
« telle souche est R/S à telle drogue ». C'est la cible supervisée du futur modèle ;
les *features* (mutations de chaque souche) se récupèrent en joignant `strain_id` au
pangénome (`data/pangenome/pan_strains.pkl`) ou à TBannotator. Les sources souche→phénotype
ne deviennent des couples mutation→résistance qu'**après** cette jointure.

## Table 1 — schéma

| colonne | description |
|---|---|
| `spdi` | `NC_000962.3:pos:ref:alt` quand résolu, sinon vide |
| `gene` | symbole du gène |
| `variant` | notation native de la source (`rrs_n.1401A>G`, `katG_p.Ser315Thr`…) |
| `drug` / `drug_code` | drogue canonique + code 3 lettres (`vocab.py`) |
| `call` | `R-associated` / `S-associated` / `uncertain` / `empirical-enriched` |
| `source` | `WHO_catalogue` / `tbprofiler` / `CRyPTIC_empirical` |
| `source_grade` | grade de confiance original de la source |
| `n_R` / `n_S` | comptes de support si disponibles |
| `metric` | métrique source (PPV OMS, diff empirique…) |
| `notes` | littérature / divers |
| `regle` | règle d'interprétation de l'OMS lue dans la colonne `comment` de tbdb 2026 (AA1, 2026-11-11) : `LIE_mmpL5` (Rv0678 pour BDQ/CFZ, « ne confère la résistance que si lié à un MmpL5 fonctionnel »), `ABROGE_Rv0678` (mmpL5 perte de fonction), `LIE_eis`, `ABROGE_eis` ; vide sinon. Posée par `phase1b_tbprofiler2026_v2.py` (`regles_interpretation.py`) ; la jointure se fait sur (gène, variante, médicament), les lignes propagées aux codons synonymes héritent de la règle de leur variant protéique |

### Sources de la Table 1 (état au 2026-10-10, catalogue v2, tbdb 2026)

| source | fichier | lignes | dont SPDI résolu | notation source |
|---|---|---|---|---|
| `WHO_catalogue` | `data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng.xlsx` (feuille `Catalogue_master_file`) | 58 994 | 58 711 | HGVS `gene_mutation` |
| `tbprofiler` | `data/sources/tbprofiler_2026/` (tbdb commit `618cf0ff`, lignes `drug_resistance`) | 4 558 (lignes développées par `spdis_for`) | 3 800 | HGVS (Gene + c./p.) |
| `CRyPTIC_empirical` | `résultats/resistance_ami.csv` | 500 | 500 | SPDI (AMK, non curé) |
| **total** | | **64 052** | **63 011** | |

18 drogues, 37 781 SPDI uniques, 8 087 lignes R-associated. Aucune ligne `*_AAprop` en v2 (le résolveur HGVS couvre les codons synonymes). Concordance vérifiée (ex. katG Ser315Thr → R-associé isoniazide dans WHO **et** tb-profiler).

**Règles d'interprétation de tbdb 2026 (piste AA, 2026-11-11 ; portées par la colonne `regle` depuis AA1).** La colonne `comment` de `data/sources/tbprofiler_2026/mutations.csv` encode les règles d'interprétation de l'OMS (2e éd.) : « Can only confer resistance if genetically linked to a functional MmpL5 » (120 lignes Rv0678 pour la bédaquiline et la clofazimine, 135 avec 15 lignes qui y ajoutent une mention de biais de recrutement), « Abrogates effect of genetically linked AwR and AwRI Rv0678 mutations » (71 lignes mmpL5, perte de fonction), les deux règles jumelles pour eis (76 et 3 lignes brutes, amikacine et kanamycine), « High-level resistance », « Low-level resistance (… additive …) » et « Confers DLM-PMD cross-resistance ». `phase1_build_catalogue.py` ne recopie pas cette colonne ; depuis AA1, `phase1b_tbprofiler2026_v2.py` en tire la colonne `regle` du catalogue (708 lignes `LIE_mmpL5`, 129 `ABROGE_Rv0678`, 6 `LIE_eis`, 132 `ABROGE_eis`) et `phase11_baseline_audit.py` comme `resistance_profile.py` neutralisent le déterminant du gène soumis (Rv0678 pour BDQ/CFZ, eis pour AMK/KAN) chez une souche qui porte un variant abrogeant (« lié » lu comme « porté par la même souche »). Les colonnes `sens`/`spec` de `baseline_audit.tsv` restent le baseline SANS règle, comparable aux versions antérieures, et `sens_regles`/`spec_regles` portent le baseline avec règles ; `resistance_profile.py --sans-regles` rend la somme brute. `update_database.py` enchaîne `phase1b --bascule` après `phase1` (qui reconstruit le catalogue v1 et l'écraserait). Les mutations « de la littérature » (996 lignes tbdb, dont 15 pour la cyclosérine et 44 pour le PAS, toutes sans SPDI) sont une option de TB-Profiler (supplément S1 de Phelan et al. 2026), laissée `uncertain` par P4.1. Mesure de l'effet sur la matrice P1.5 : `analyses/aa_regles_epistasie.py` → `résultats/model/aa_regles_epistasie.tsv`.

**Confiance vide de tbdb (décision P4.1, 2026-10-10).** 653 lignes `tbprofiler` ont `source_grade` vide : ce sont les entrées tbdb « Mutation from literature » (644 lignes brutes), sans grade OMS. Elles restent `call=uncertain`, comme dans le catalogue v1, donc hors des déterminants R. Motifs : aucune n'est dans le catalogue OMS (0 recoupement sur 653 couples variant/drogue), un `type=drug_resistance` non gradé est une assertion de la littérature et non une validation experte, et l'effet d'un reclassement en R est quasi nul (19 lignes ont un SPDI, dont 2 seulement existent dans les 260 000 souches de TBannotator : `rrs` 1338A>C amikacine, 8 porteurs, aucun phénotypé ; `gyrB` Asn499Tyr moxifloxacine, 4 porteurs, 1 phénotypé R). Les 634 autres lignes sont non résolues (LoF, indels longs) et sans effet.

Le compte est en **lignes** : une ligne = une représentation SPDI d'un triplet (variant, drogue, source).
Une assertion logique unique peut donner plusieurs lignes — un variant protéique se résout en tous ses codons
synonymes, un indel en toutes ses normalisations (corrections des angles morts 1, 3 et 4). D'où la croissance de
51 017 (avant expansion des indels) à 62 019, **sans nouvelle assertion biologique**. Le repère historique
« 50 963 » = somme des 3 sources principales avant propagation AA (état figé dans le manuscrit).

### Codons synonymes (angle mort 1) — `*_AAprop` désormais résiduel

Un déterminant connu via un codon (embB Met306Ile par ATG→ATT) ne doit pas être manqué via un
codon synonyme (ATG→ATC), même AA et même résistance. **Cette couverture est aujourd'hui portée
nativement par le résolveur HGVS** : un variant en notation protéique (`gene_p.X`) est résolu en
**toutes** les représentations codon de l'AA cible (vérifié : `embB_p.Met306Ile` → les 3 SPDI
`4247430:G:A`, `G:C`, `G:T` sous `WHO_catalogue`). La propagation explicite `aa_propagation.py`
(conversion AA→SPDI avec gestion du brin, validée katG S315T brin − / embB M306I brin +,
appliquée en fin de `phase1`) ne complète plus que des cas **résiduels** (2 lignes `*_AAprop`).
**Impact baseline mesuré : éthambutol 0.63→0.86 sensibilité (+23 pts), PZA +5 pts**, spécificité quasi inchangée.

### Le pont HGVS → SPDI (`analyses/hgvs_to_spdi.py`)

Le SPDI est résolu à partir du fichier officiel de coordonnées génomiques de l'OMS,
`WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt` (`variant chromosome position ref alt`,
**1-based**). Le pangénome du projet est **0-based** : on applique `position - 1`
(vérifié : rrs_n.1401A>G → VCF 1473246 → SPDI projet `NC_000962.3:1473245:A:G`, présent
dans `pan_spdi.pkl`, 5413 souches). La notation `variant` de l'OMS étant identique au
`f"{Gene}_{Mutation}"` de tb-profiler, le même résolveur mappe les deux sources.

## Limites résiduelles

1. **Indels / MNV — désormais résolus (angle mort 4 corrigé).** Le résolveur étend chaque
   indel/MNV à **toutes** ses représentations (`resolve_all`), et `variant_match` les apparie
   (SNP exact, décomposition MNV, indel exact). La couverture SPDI est passée de 90 % à **97 %**
   (60 423/62 019). Résiduel non résolu : ~1 600 lignes (3 %), surtout entrées gène-niveau (LoF)
   et indels dont la normalisation est absente du fichier de coordonnées OMS.
2. **tb-profiler partiellement résolu (48 %).** Le reste est constitué d'entrées
   gène-niveau (LoF), de plages de codons ou de variants absents du fichier de
   coordonnées OMS. Le savoir expert reste porté par l'OMS pour ces drogues.
3. **Couche empirique limitée à l'amikacine et biaisée.** Les 500 SNPs
   `CRyPTIC_empirical` viennent du différentiel `n_R − n_S` confondu par le MDR
   (cf. cahier de labo). À régénérer en lineage-matched et étendre aux 13 drogues.
4. **Doublons inter-sources volontaires.** WHO et tb-profiler peuvent porter la même
   mutation : on conserve les deux lignes (traçabilité). Pivoter par `spdi` pour
   consolider un verdict unique côté consommateur.

## Table 2 — schéma

`strain_id`, `drug`, `drug_code`, `phenotype` (R/S), `mic`, `source`.

| source | fichier | mesures |
|---|---|---|
| CRyPTIC | `data/sources/cryptic/CRyPTIC_reuse_table_20240917.csv` | 154 844 |
| TBannotator_biosample | antibiogramme DST mesuré (client MCP, `analyses/ingestion/fetch_tbannotator.py`) | 90 835 |
| souches_datees | « Toutes les souches datées.csv » (export TBannotator de décembre 2022, **étiquettes calculées sur le génotype, non mesurées** : P8.4, `p84_provenance_datees.tsv`) | 44 004 |
| PATRIC | `résultats/literature_extraction/patric_corrige.csv` (P8.4 : seul un phénotype Resistant/Susceptible non vide donne une étiquette ; l'ancien `patric.csv` codait R toute ligne de CMI sans phénotype et laissait la dernière ligne d'un couple gagner, d'où ses 20 746 cases ; **provenance auditée en P8.5** : étiquettes retenues à 99,9 % de dilution en bouillon ou sur gélose, aucune prédiction calculée, appariement vers SRA correct à 97 %) | 19 091 |
| Coll2018 | panel GWAS publié (Coll et al. 2018) | 11 529 |
| littérature | `résultats/literature_extraction/res.csv` | 123 |

Total : **810 822 mesures, 88 017 souches, 22 drogues** (régénéré le 2026-11-11 par P1.3, qui branche la source `CRyPTIC_v3`, jeu CRyPTIC élargi du dépôt Zenodo 22128395 v3.4.1, CC-BY-4.0 : +490 396 mesures, +35 521 souches, +268 200 paires au consensus, soit +84,9 %. Avant : 320 426 mesures et 52 496 souches, copies `archives/sauvegardes/2026-11/*avant_P1.3`). Le branchement est **additif** : aucune paire perdue, et seulement 616 étiquettes préexistantes modifiées sur 315 995 (0,19 %), toutes là où `CRyPTIC_v3` prend la main sur `souches_datees` (547) ou `Coll2018` (69) ; il ne change rien à celles que portait PATRIC (14 954 paires reprises, zéro désaccord). Les conflits inter-sources passent de 0,2 % à 0,76 % (4 435 paires), mécaniquement, le recouvrement entre sources ayant augmenté. **Deux sources CRyPTIC coexistent** : `CRyPTIC` (reuse table 2022) reste prioritaire parce qu'elle seule porte les MIC quantitatives, et `CRyPTIC_v3` la suit. **Les mesures de `CRyPTIC_v3` viennent de douze cohortes de DST dont plusieurs ont servi à construire le catalogue OMS** : tout chiffre tiré de cette source se donne par cohorte, en rejoignant `data/sources/cryptic/cryptic_v3_dst_by_sra.csv` sur (strain_id, drug_code) pour ses colonnes `quality`, `dst_source` et `cohort` (motif complet : `REVUE_CRITIQUE.md` §2.6). Toutes les sources sont des étiquettes de laboratoire **sauf `souches_datees`**, dont les étiquettes se comportent comme une prédiction du catalogue (P5.2, P8.4 : `résultats/model/p52_circularite.tsv`, `p84_provenance_datees.tsv`) ; la cible de référence des rapports l'écarte, PATRIC corrigé y reste (règle de provenance, jamais d'accord avec le catalogue). Les « R » de PATRIC en seconde ligne que l'ancien encodage rendait inexplicables étaient un défaut d'encodage ; après correction, PATRIC reste hétérogène avec les autres sources pour MXF, AMK et CAP (`p85_baseline_par_source.tsv`). Vue dédoublonnée par paire (souche, drogue), avec priorité de source :
`phenotypes_souches_consensus.tsv`.
