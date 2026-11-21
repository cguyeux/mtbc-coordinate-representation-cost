#!/usr/bin/env python3
"""P7.8 — Test de non-régression : chaque chiffre du manuscrit RECADRÉ renvoie à un TSV.

Objet      : même patron que `test_p75_chiffres_manuscrit.py` (P7.5), appliqué au manuscrit réécrit sur
             la thèse recadrée (coût d'une chaîne qui garde une seule représentation de coordonnée par
             variant). Chaque ligne de `CLAIMS` donne le chiffre TEL QU'IMPRIMÉ dans `article/main.tex`,
             le TSV et le filtre qui le produisent, et la demi-unité de l'arrondi imprimé.
             Le régime de mesure est la matrice à BLOCS BRUTS (50 782 souches, P7.6/P7.7) : ces chiffres
             ne se comparent pas à ceux de la matrice canonique (85 509 souches, P7.5), qui sont épinglés
             par `test_p75_chiffres_manuscrit.py` et décrivent le manuscrit d'AVANT.
Entrées    : résultats/model/{p77_C0.tsv, p77_C1.tsv, p77_C0_residuel.tsv, p77_C1_residuel.tsv,
             p77_catalogues.tsv, p77_conflits_C1.tsv, p77_C0_matriceP76.tsv, p74_angles_morts_blocs_P15.tsv,
             p76_effectifs_blocs.tsv, p78_representations_oms.tsv, p75_effectifs_cohorte.tsv}.
Sorties    : aucune (code de retour pytest).
Réutilisable : oui, comme patron.
Projet     : Resistance_antibio (P7.8)
Date       : 2026-11-16
Usage      : python3 -m pytest -q analyses/test_p78_chiffres_manuscrit.py
"""
from pathlib import Path

import pandas as pd
import pytest

M = Path(__file__).resolve().parent.parent / "résultats" / "model"


def lire(nom):
    return pd.read_csv(M / nom, sep="\t")


def etat(fichier, drug, etat_, col, cible="reference"):
    t = lire(fichier)
    r = t[(t["cible"] == cible) & (t["drug"] == drug) & (t["etat"] == etat_)].iloc[0]
    return float(r[col])


def pts(fichier, drug, etat_, col):
    return 100 * etat(fichier, drug, etat_, col)


C0, C1 = "p77_C0.tsv", "p77_C1.tsv"


def _claims():
    c = []

    # ---- Section 3.1 et résumé : ce que la source fournit, ce que l'ingestion garde ----
    p78 = lire("p78_representations_oms.tsv").set_index("indicateur")["valeur"].astype(float)
    for nom, imprime, cle, demi in [
            ("lignes du fichier de coordonnées", 144964, "lignes_fichier_oms", 0.5),
            ("variants du fichier", 41358, "variants_fichier_oms", 0.5),
            ("représentations par variant", 3.5, "representations_par_variant", 0.05),
            ("variants à plusieurs représentations", 25118, "variants_a_plusieurs_representations", 0.5),
            ("variants protéiques", 18732, "variants_proteiques", 0.5),
            ("variants protéiques à plusieurs représentations", 16184,
             "variants_proteiques_a_plusieurs_representations", 0.5),
            ("lignes SNP", 32155, "lignes_oms_snp", 0.5),
            ("lignes en bloc", 97787, "lignes_oms_mnv", 0.5),
            ("lignes indel", 11146, "lignes_oms_indel", 0.5),
            ("lignes de plus de 24 bases", 3876, "lignes_oms_long", 0.5),
            ("médicaments de la source OMS", 15, "medicaments_de_la_source_oms", 0.5),
            ("médicaments du catalogue", 18, "medicaments_du_catalogue_C0", 0.5),
            ("lignes du catalogue ingéré", 64052, "lignes_catalogue_C0", 0.5),
            ("variants communs au catalogue et au fichier", 30534,
             "variants_du_catalogue_retrouves_dans_oms", 0.5),
            ("représentations que le résolveur peut émettre", 120917,
             "representations_fournies_par_oms", 0.5),
            ("représentations retenues", 37309, "representations_retenues_par_la_chaine", 0.5),
            ("représentations écartées", 83608, "representations_perdues", 0.5),
            ("représentations écartées en bloc", 77073, "representations_perdues_en_bloc", 0.5),
            ("variants perdant au moins une représentation", 19154,
             "variants_perdant_au_moins_une_representation", 0.5),
            ("SPDI en bloc R-associés de l'émission complète", 14670,
             "spdi_bloc_R_associes_de_l_emission_complete", 0.5)]:
        c.append((nom, imprime, float(p78[cle]), demi))
    # pourcentages imprimés, recalculés depuis le même TSV
    c += [("part des variants à plusieurs représentations (%)", 60.7,
           100 * p78["variants_a_plusieurs_representations"] / p78["variants_fichier_oms"], 0.05),
          ("part des variants protéiques à plusieurs représentations (%)", 86.4,
           100 * p78["variants_proteiques_a_plusieurs_representations"] / p78["variants_proteiques"], 0.05),
          ("part des lignes de plus de 24 bases (%)", 2.7,
           100 * p78["lignes_oms_long"] / p78["lignes_fichier_oms"], 0.05),
          ("part des représentations écartées (%)", 69.1,
           100 * p78["representations_perdues"] / p78["representations_fournies_par_oms"], 0.05),
          ("part des variants perdant une représentation (%)", 62.7,
           100 * p78["variants_perdant_au_moins_une_representation"]
           / p78["variants_du_catalogue_retrouves_dans_oms"], 0.05),
          # P7.2 : part en BLOC des représentations écartées, qui mesure la distance entre la
          # classe de variants de Verboven 2022 (indels) et celle que mesure A1 (blocs).
          ("part en bloc des représentations écartées (%)", 92.2,
           100 * p78["representations_perdues_en_bloc"] / p78["representations_perdues"], 0.05)]

    cat = lire("p77_catalogues.tsv")
    uniques = cat["indicateur"].value_counts()
    ind = cat[cat["indicateur"].isin(uniques[uniques == 1].index)].set_index("indicateur")["valeur"]
    types = cat[cat["indicateur"] == "types_spdi_C0"].set_index("detail")["valeur"].astype(float)
    # P7.16 : ce sont des LIGNES d'assertion, non des chaînes distinctes (review du 2026-11-21,
    # point 5 : le nombre était exact contre sa source, c'est le mot qui le nommait qui était faux).
    # Les chaînes distinctes par type sont épinglées juste après, et leur somme doit redonner les
    # 37 781 chaînes distinctes que le manuscrit imprime ailleurs : c'est le contrôle d'unité que la
    # review demandait, et qu'un test qui épingle des valeurs ne peut pas déduire tout seul.
    c += [("lignes à chaîne mono-base du catalogue ingéré", 52041, float(types["snp"]), 0.5),
          ("lignes à chaîne en bloc du catalogue ingéré", 772, float(types["mnv"]), 0.5)]
    dist = lire("p716_chaines_distinctes.tsv").set_index("type")["chaines_distinctes"].astype(float)
    c += [("chaînes distinctes mono-base, C0", 32596, float(dist["snp"]), 0.5),
          ("chaînes distinctes en bloc, C0", 543, float(dist["mnv"]), 0.5),
          ("chaînes distinctes d'indel, C0", 4642, float(dist["indel"]), 0.5),
          ("somme des chaînes distinctes par type = chaînes distinctes du catalogue", 37781,
           float(dist["snp"] + dist["mnv"] + dist["indel"]), 0.5)]

    # ---- Codon 94 de gyrA : compte re-dérivé du fichier de coordonnées (P7.16, review point 1) ----
    g94 = lire("p716_gyra94_synthese.tsv").set_index("indicateur")["valeur"].astype(float)
    c += [("variants protéiques du codon 94 de gyrA", 8, float(g94["variants_du_codon"]), 0.5),
          ("ceux qui portent une entrée combinée", 5,
           float(g94["variants_a_entree_combinee"]), 0.5),
          ("ceux dont l'entrée combinée est un bloc de cinq bases", 2,
           float(g94["variants_a_bloc_cinq_bases"]), 0.5),
          ("ceux sans entrée combinée (Asn, His, Tyr)", 3,
           float(g94["variants_sans_entree_combinee"]), 0.5)]

    # ---- Section 3.2 : les deux fluoroquinolones sous les deux ingestions ----
    for d, nR, nS, ex0, ex1, comp, g0, g0lo, g0hi, g1, g1lo, g1hi, ps0, ps1 in [
            ("moxifloxacin", 2916, 14931, 0.267, 0.709, 0.841, 57.4, 29, 63, 13.2, 8, 16, 3.15, 0.58),
            ("levofloxacin", 2358, 10982, 0.300, 0.682, 0.799, 50.0, 24, 57, 11.7, 6, 13, 1.63, 0.30)]:
        c += [(f"{d} n_R", nR, etat(C0, d, "complet", "n_R"), 0.5),
              (f"{d} n_S", nS, etat(C0, d, "complet", "n_S"), 0.5),
              (f"{d} sens. appariement exact seul, C0", ex0, etat(C0, d, "mnv", "sens"), 0.0005),
              (f"{d} sens. appariement exact seul, C1", ex1, etat(C1, d, "mnv", "sens"), 0.0005),
              (f"{d} sens. baseline complet, C0", comp, etat(C0, d, "complet", "sens"), 0.0005),
              (f"{d} sens. baseline complet, C1", comp, etat(C1, d, "complet", "sens"), 0.0005),
              (f"{d} gain de la décomposition, C0", g0, pts(C0, d, "mnv", "gain_sens"), 0.05),
              (f"{d} gain IC bas, C0", g0lo, pts(C0, d, "mnv", "gain_lo"), 0.5),
              (f"{d} gain IC haut, C0", g0hi, pts(C0, d, "mnv", "gain_hi"), 0.5),
              (f"{d} gain de la décomposition, C1", g1, pts(C1, d, "mnv", "gain_sens"), 0.05),
              (f"{d} gain IC bas, C1", g1lo, pts(C1, d, "mnv", "gain_lo"), 0.5),
              (f"{d} gain IC haut, C1", g1hi, pts(C1, d, "mnv", "gain_hi"), 0.5),
              (f"{d} perte de spécificité, C0", ps0, -pts(C0, d, "mnv", "perte_spec"), 0.005),
              (f"{d} perte de spécificité, C1", ps1, -pts(C1, d, "mnv", "perte_spec"), 0.005)]
    # Part du gain qui tient à la règle d'ingestion, et son intervalle (P7.16, review point 4).
    # L'intervalle du RAPPORT vient des 2 000 tirages APPARIÉS des deux bras (p716_intervalles.py) :
    # il ne se déduit pas des bornes des deux gains, qui donneraient 47 % en les combinant à tort,
    # parce que les deux gains montent et descendent ensemble avec la fréquence de Ser95Thr.
    inter = lire("p716_intervalles.tsv").set_index("drug")
    for d, part, plo, phi in [("moxifloxacin", 77.0, 66.2, 79.5),
                              ("levofloxacin", 76.5, 67.7, 79.7)]:
        lu = 100 * (1 - pts(C1, d, "mnv", "gain_sens") / pts(C0, d, "mnv", "gain_sens"))
        r = inter.loc[d]
        c += [(f"{d} part du gain imputable à l'ingestion (%)", part, lu, 0.05),
              (f"{d} part, même chiffre par les tirages appariés (%)", part,
               float(r["part_ingestion"]), 0.05),
              (f"{d} part IC bas (%)", plo, float(r["part_lo"]), 0.05),
              (f"{d} part IC haut (%)", phi, float(r["part_hi"]), 0.05)]
    # Intervalles des gains portés au résumé, à une décimale (le tableau 1 les donne arrondis à
    # l'entier) ; contrôle d'ancrage : le rejeu apparié reproduit les gains publiés par p77.
    for d, g0, g0lo, g0hi, g1, g1lo, g1hi in [
            ("moxifloxacin", 57.4, 29.3, 63.3, 13.2, 7.9, 15.6),
            # P7.16 : 11.7 et non 11.8. Le manuscrit imprimait 11.8 par DOUBLE ARRONDI, la valeur
            # 0.117472 étant écrite 0.1175 dans p77_C1.tsv (float_format %.4f) puis lue 11.75, qui
            # arrondit à 11.8 ; la grandeur vaut 11.7472 et s'arrondit à 11.7. Le test d'avant
            # passait parce qu'il comparait l'imprimé au TSV déjà arrondi, à une demi-unité de 0.05.
            ("levofloxacin", 50.0, 24.1, 57.4, 11.7, 6.3, 13.3)]:
        r = inter.loc[d]
        c += [(f"{d} gain C0, rejeu apparié", g0, float(r["gain_C0"]), 0.05),
              (f"{d} gain C0 IC bas, résumé", g0lo, float(r["gain_C0_lo"]), 0.05),
              (f"{d} gain C0 IC haut, résumé", g0hi, float(r["gain_C0_hi"]), 0.05),
              (f"{d} gain C1, rejeu apparié", g1, float(r["gain_C1"]), 0.05),
              (f"{d} gain C1 IC bas, résumé", g1lo, float(r["gain_C1_lo"]), 0.05),
              (f"{d} gain C1 IC haut, résumé", g1hi, float(r["gain_C1_hi"]), 0.05)]
    # hors fluoroquinolones, le gain de la décomposition tombe à zéro sous C1
    for d, v0 in [("isoniazid", 0.99), ("rifampicin", 0.26), ("streptomycin", 0.15),
                  ("ethionamide", 0.03)]:
        c += [(f"{d} gain décomposition C0", v0, pts(C0, d, "mnv", "gain_sens"), 0.005),
              (f"{d} gain décomposition C1", 0.0, pts(C1, d, "mnv", "gain_sens"), 0.005)]

    # ---- Section 3.3 : le résiduel ----
    r1 = lire("p77_C1_residuel.tsv")
    mxf = r1[r1["drug"] == "moxifloxacin"].set_index("spdi")
    lev = r1[r1["drug"] == "levofloxacin"].set_index("spdi")
    A, T, G = ("NC_000962.3:7580:GACAG:AACAC", "NC_000962.3:7580:GACAG:TACAC",
               "NC_000962.3:7580:GACAG:CACAC")
    c += [("colonnes résiduelles sous C1", 3, float(r1["spdi"].nunique()), 0.5),
          ("souches Asp94Asn", 433, float(mxf.loc[A, "n_souches"]), 0.5),
          ("souches Asp94Tyr", 221, float(mxf.loc[T, "n_souches"]), 0.5),
          ("souches Asp94His", 86, float(mxf.loc[G, "n_souches"]), 0.5),
          ("souches des trois haplotypes", 740,
           float(r1.drop_duplicates("spdi")["n_souches"].sum()), 0.5),
          ("MXF résistantes Asp94Asn", 239, float(mxf.loc[A, "n_R"]), 0.5),
          ("MXF résistantes Asp94Tyr", 111, float(mxf.loc[T, "n_R"]), 0.5),
          ("MXF résistantes Asp94His", 44, float(mxf.loc[G, "n_R"]), 0.5),
          ("MXF sensibles des trois", 86, float(mxf["n_S"].sum()), 0.5),
          ("LEV résistantes Asp94Asn", 172, float(lev.loc[A, "n_R"]), 0.5),
          ("LEV résistantes Asp94Tyr", 70, float(lev.loc[T, "n_R"]), 0.5),
          ("LEV résistantes Asp94His", 41, float(lev.loc[G, "n_R"]), 0.5),
          ("LEV sensibles des trois", 33, float(lev["n_S"].sum()), 0.5),
          ("MXF part des R portée par les trois (pts)", 13.5,
           100 * mxf["n_R"].sum() / etat(C1, "moxifloxacin", "complet", "n_R"), 0.05),
          ("LEV part des R portée par les trois (pts)", 12.0,
           100 * lev["n_R"].sum() / etat(C1, "levofloxacin", "complet", "n_R"), 0.05)]
    # P7.16 (review point 10) : la part des porteurs DÉPASSE la différence mesurée, et le manuscrit
    # dit maintenant de combien, au lieu de surénoncer une clôture exacte.
    for nom, res, d in (("MXF", mxf, "moxifloxacin"), ("LEV", lev, "levofloxacin")):
        ecart = (100 * res["n_R"].sum() / etat(C1, d, "complet", "n_R")
                 - float(inter.loc[d, "gain_C1"]))
        c.append((f"{nom} excès de la part des porteurs sur la différence (pts)", 0.3, ecart, 0.05))
    r0 = lire("p77_C0_residuel.tsv")
    c += [("paires résiduelles sous C0", 49, float(len(r0)), 0.5),
          ("colonnes résiduelles distinctes sous C0", 38, float(r0["spdi"].nunique()), 0.5),
          ("médicaments du résiduel sous C0", 9, float(r0["drug"].nunique()), 0.5)]

    # ---- Section 3.4 : tableau par médicament ----
    for d, nR, nS, co, i0, i1, r_0, r_1, s0, s1 in [
            ("isoniazid", 13256, 25362, 0.17, 0.41, 0.45, 1.05, 1.09, 0.903, 0.904),
            ("rifampicin", 10019, 28665, 0.04, 0.85, 0.85, 1.41, 1.41, 0.931, 0.931),
            ("ethambutol", 5610, 31142, 5.42, 0.00, 0.00, 0.11, 0.11, 0.830, 0.830),
            ("pyrazinamide", 2893, 12962, 1.66, 10.89, 11.55, 17.97, 18.63, 0.733, 0.739),
            ("moxifloxacin", 2916, 14931, 0.00, 0.00, 0.00, 1.27, 1.27, 0.841, 0.841),
            ("levofloxacin", 2358, 10982, 0.00, 0.00, 0.00, 1.27, 1.27, 0.799, 0.799),
            ("streptomycin", 3990, 7642, 0.28, 5.44, 5.46, 2.66, 2.68, 0.764, 0.764),
            ("ethionamide", 3262, 13710, 1.13, 8.28, 8.58, 7.14, 7.45, 0.711, 0.714),
            ("amikacin", 1651, 17795, 0.00, 0.00, 0.00, 0.00, 0.00, 0.738, 0.738),
            ("kanamycin", 1946, 16998, 0.00, 0.05, 0.05, 0.05, 0.05, 0.742, 0.742),
            ("capreomycin", 991, 8806, 0.20, 4.04, 4.24, 4.24, 4.44, 0.692, 0.694)]:
        c += [(f"tableau {d} n_R", nR, etat(C0, d, "complet", "n_R"), 0.5),
              (f"tableau {d} n_S", nS, etat(C0, d, "complet", "n_S"), 0.5),
              (f"tableau {d} codons C0", co, pts(C0, d, "synonymes", "gain_sens"), 0.005),
              (f"tableau {d} codons C1", co, pts(C1, d, "synonymes", "gain_sens"), 0.005),
              (f"tableau {d} indels C0", i0, pts(C0, d, "indel", "gain_sens"), 0.005),
              (f"tableau {d} indels C1", i1, pts(C1, d, "indel", "gain_sens"), 0.005),
              (f"tableau {d} rares C0", r_0, pts(C0, d, "rares", "gain_sens"), 0.005),
              (f"tableau {d} rares C1", r_1, pts(C1, d, "rares", "gain_sens"), 0.005),
              (f"tableau {d} sens. C0", s0, etat(C0, d, "complet", "sens"), 0.0005),
              (f"tableau {d} sens. C1", s1, etat(C1, d, "complet", "sens"), 0.0005)]
    # texte de 3.4 : sensibilités à quatre décimales et hausse maximale
    for d, v0, v1 in [("pyrazinamide", 0.7325, 0.7390), ("ethionamide", 0.7112, 0.7143),
                      ("capreomycin", 0.6922, 0.6942), ("isoniazid", 0.9034, 0.9039),
                      ("streptomycin", 0.7642, 0.7644)]:
        c += [(f"{d} sens. complète C0 (4 déc.)", v0, etat(C0, d, "complet", "sens"), 0.00005),
              (f"{d} sens. complète C1 (4 déc.)", v1, etat(C1, d, "complet", "sens"), 0.00005)]
    c.append(("hausse maximale d'une contribution C0 vers C1 (pts)", 0.66,
              pts(C1, "pyrazinamide", "indel", "gain_sens")
              - pts(C0, "pyrazinamide", "indel", "gain_sens"), 0.005))
    ecart_spec = max(abs(etat(C1, d, "complet", "spec") - etat(C0, d, "complet", "spec"))
                     for d in lire(C0)["drug"].unique())
    c.append(("écart maximal de spécificité C0 vers C1 (pts)", 0.07, 100 * ecart_spec, 0.005))

    # ---- Section 3.5 : coût de l'émission complète ----
    for nom, imprime, cle, demi in [
            ("lignes du catalogue à blocs", 204522, "lignes_C1_blocs", 0.5),
            ("lignes ajoutées", 140470, "lignes_ajoutees_par_C1", 0.5),
            ("SPDI distincts du catalogue ingéré", 37781, "spdi_distincts_C0", 0.5),
            ("SPDI distincts du catalogue à blocs", 115097, "spdi_distincts_C1", 0.5),
            ("SPDI R-associés distincts, C0", 4242, "spdi_distincts_R_associes_C0", 0.5),
            ("SPDI R-associés distincts, C1", 23958, "spdi_distincts_R_associes_C1", 0.5),
            ("couples à verdicts contradictoires, C1", 19,
             "conflits_R_S_meme_spdi_meme_medicament_C1", 0.5)]:
        c.append((nom, imprime, float(ind[cle]), demi))
    ajout_type = cat[cat["indicateur"] == "ajoutees_par_type"].set_index("detail")["valeur"].astype(float)
    ajout_verdict = cat[cat["indicateur"] == "ajoutees_par_verdict"].set_index("detail")["valeur"].astype(float)
    c += [("blocs ajoutés", 129459, float(ajout_type["mnv"]), 0.5),
          ("indels ajoutés", 11011, float(ajout_type["indel"]), 0.5),
          ("verdicts incertains ajoutés", 99103, float(ajout_verdict["uncertain"]), 0.5),
          ("verdicts R ajoutés", 23540, float(ajout_verdict["R-associated"]), 0.5),
          ("verdicts S ajoutés", 17827, float(ajout_verdict["S-associated"]), 0.5),
          ("facteur de taille du catalogue", 3.19,
           float(ind["lignes_C1_blocs"]) / float(ind["lignes_catalogue_C0"])
           if "lignes_catalogue_C0" in ind.index else
           float(ind["lignes_C1_blocs"]) / float(p78["lignes_catalogue_C0"]), 0.005),
          ("facteur des SPDI R-associés distincts", 5.6,
           float(ind["spdi_distincts_R_associes_C1"]) / float(ind["spdi_distincts_R_associes_C0"]), 0.05)]
    colonnes = {n: lire(f).groupby("drug")["n_col"].first().sum()
                for n, f in [("C0", C0), ("C1", C1)]}
    for n, f in [("C0", C0), ("C1", C1)]:
        t = lire(f)
        colonnes[n] = t[(t["cible"] == "reference") & (t["etat"] == "complet")]["n_col"].sum()
    c += [("colonnes à déterminant, C0", 1389, float(colonnes["C0"]), 0.5),
          ("colonnes à déterminant, C1", 1565, float(colonnes["C1"]), 0.5),
          ("colonnes gagnées", 176, float(colonnes["C1"] - colonnes["C0"]), 0.5)]
    conf = pd.read_csv(M / "p77_conflits_C1.tsv", sep="|")
    c += [("couples contradictoires, fichier", 19, float(len(conf)), 0.5),
          ("couples contradictoires de la rifampicine", 11,
           float((conf["drug"] == "rifampicin").sum()), 0.5),
          ("couples contradictoires des fluoroquinolones", 5,
           float(conf["drug"].isin(["levofloxacin", "moxifloxacin"]).sum()), 0.5),
          ("chaînes distinctes contradictoires", 17, float(conf["spdi"].nunique()), 0.5)]

    # ---- Méthodes : effectifs et contrôle de la mesure ----
    ef = lire("p76_effectifs_blocs.tsv").set_index("indicateur")["valeur"]
    c += [("souches à blocs bruts", 50782, float(ef["souches_genotypees"]), 0.5),
          ("souches à blocs bruts phénotypées", 41551,
           float(ef["souches_phenotypees_reference_11_medicaments"]), 0.5),
          ("souches du pangénome", 41484, float(ef["genotypees_source_pangenome"]), 0.5),
          ("souches de la base de souches", 9298, float(ef["genotypees_source_bdd"]), 0.5)]
    ecart_rif = abs(etat("p77_C0_matriceP76.tsv", "rifampicin", "complet", "sens")
                    - etat("p74_angles_morts_blocs_P15.tsv", "rifampicin", "complet", "sens"))
    c.append(("écart du contrôle de la mesure sur la rifampicine (pts)", 0.14,
              100 * ecart_rif, 0.005))
    eff = lire("p75_effectifs_cohorte.tsv").iloc[0]
    c.append(("souches génotypées de la matrice canonique", 85509,
              float(eff["souches_genotypees"]), 0.5))

    # ---- Méthodes : hétérogénéité des sources de test (P7.14, analyses/p714_sources_etudes.py) ----
    sy = lire("p714_synthese.tsv").set_index("indicateur")["valeur"]
    n_paires = float(sy["paires_etiquetees_11_medicaments"])
    n_v3 = float(sy["paires_cryptic_v3"])
    dst = lire("p714_dst_cryptic_v3.tsv")
    v10 = dst[dst["cohorte_cryptic"] == "CRyPTIC-v1.0"]
    n_v10 = float(v10["n_R"].sum() + v10["n_S"].sum())
    bp = lire("p714_etudes_par_medicament.tsv")
    bp = bp[bp["bioproject"] == "PRJNA1237251"].set_index("drug")
    c += [("paires étiquetées, onze médicaments", 237887, n_paires, 0.5),
          ("paires issues de CRyPTIC élargi", 121620, n_v3, 0.5),
          ("part de CRyPTIC élargi dans les paires (%)", 51.1, 100 * n_v3 / n_paires, 0.05),
          ("sources de DST de CRyPTIC élargi présentes dans la cohorte", 11,
           float(sy["dst_sources_retrouvees_dans_les_paires_cryptic_v3"]), 0.5),
          ("paires de la cohorte CRyPTIC-v1.0", 92174, n_v10, 0.5),
          ("part de CRyPTIC-v1.0 dans les paires (%)", 38.7, 100 * n_v10 / n_paires, 0.05)]
    for d, imp in (("isoniazid", 24.3), ("rifampicin", 23.9), ("ethambutol", 22.0)):
        c.append((f"part de PRJNA1237251 dans les S, {d} (%)", imp, 100 * float(bp.loc[d, "part_S"]), 0.05))
    # « under 6 % » : borne imprimée = plafond de la part maximale mesurée (5,64 %) ; 5.6 est réservé au
    # RETIRES de test_p78_chiffres_dans_texte.py (autre grandeur), d'où l'arrondi à l'entier supérieur.
    part_r_max = 100 * float(bp.loc[["isoniazid", "rifampicin", "ethambutol"], "part_R"].max())
    c.append(("borne de la part de PRJNA1237251 dans les R des trois (%)", 6, float(-(-part_r_max // 1)), 0.5))
    return c


CLAIMS = _claims()


@pytest.mark.parametrize("nom,imprime,lu,demi_unite", CLAIMS, ids=[x[0] for x in CLAIMS])
def test_chiffre_imprime(nom, imprime, lu, demi_unite):
    assert abs(imprime - lu) <= demi_unite + 1e-9, f"{nom} : manuscrit {imprime}, TSV {lu}"
