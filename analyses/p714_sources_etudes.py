#!/usr/bin/env python3
"""P7.14 : part de chaque source de test et de chaque étude dans les paires de phénotype de la cohorte d'A1.

Objet      : l'article A1 écrit, dans « Genotypes, phenotypes and the matrix that keeps raw calls », que les
             étiquettes sont hétérogènes : « the twelve testing sources use different methods », « several of
             the cohorts are reported by their authors to have served in building the catalogues » et « one
             surveillance study contributes a large share of the susceptible isolates of three first-line
             medicines ». Aucune pièce chiffrée ne reliait ces énoncés à la cohorte de 50 782 souches. Ce script
             mesure, pour les onze médicaments du tableau 1 et sous la cible de RÉFÉRENCE, la source de chaque
             paire (souche, médicament) étiquetée, la source de DST et la cohorte des paires issues de
             CRyPTIC élargi, le BioProject ENA de la souche, et rend les parts des R et des S.
Cohorte    : souches de `rows.tsv` (matrice P1.5, 85 509 souches) dont la source de SPDI est pangénome ou bdd,
             soit 41 484 + 9 298 = 50 782, la population de la matrice à blocs bruts de P7.6. Les étiquettes ne
             dépendent pas de la représentation des variants (blocs ou SNP), donc la matrice creuse n'est pas lue.
             Contrôles imprimés : 41 551 souches phénotypées (`p76_effectifs_blocs.tsv`) et n_R, n_S par
             médicament égaux à ceux de `p74_angles_morts_blocs_P15.tsv`, cible `reference`.
Entrées    : results/model/{rows.tsv, souches_source_spdi.tsv, p74_angles_morts_blocs_P15.tsv,
             p76_effectifs_blocs.tsv} ; catalogue/{phenotypes_souches.tsv, etudes_ena.tsv} ;
             data/sources/cryptic/cryptic_v3_dst_by_sra.csv.
Sorties    : results/model/p714_sources_par_medicament.tsv (médicament x source de table),
             p714_dst_cryptic_v3.tsv (médicament x source de DST x cohorte, paires issues de CRyPTIC élargi),
             p714_etudes_par_medicament.tsv (médicament x BioProject), p714_synthese.tsv (indicateurs).
Réutilisable : non (effectifs propres à A1) ; la jointure paire -> source -> BioProject l'est.
Projet     : Resistance_antibio (P7.14)
Date       : 2026-11-21
Usage      : python3 analyses/p714_sources_etudes.py [--etude PRJNA1237251]
Limite     : le BioProject est celui du premier run d'un isolat composite `ERRa.ERRb` (même règle que P8.8).
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths  # noqa: E402
import cible_reference as cr  # noqa: E402
from p74_angles_morts import DRUGS  # noqa: E402

M = paths.RESULTATS / "model"
CODES = {"isoniazid": "INH", "rifampicin": "RIF", "ethambutol": "EMB", "pyrazinamide": "PZA",
         "moxifloxacin": "MXF", "levofloxacin": "LEV", "amikacin": "AMK", "kanamycin": "KAN",
         "streptomycin": "STM", "ethionamide": "ETH", "capreomycin": "CAP"}


def cohorte():
    """rows.tsv restreint aux 50 782 souches de la matrice à blocs bruts (SPDI du pangénome ou de bdd)."""
    rows = pd.read_csv(M / "rows.tsv", sep="\t", keep_default_na=False)
    src = pd.read_csv(M / "souches_source_spdi.tsv", sep="\t").set_index("strain_id")["source_spdi"]
    garde = rows["strain_id"].map(src).isin(["pangenome", "bdd"])
    return rows[garde].reset_index(drop=True)


def paires(rows):
    """Une ligne par paire (souche, médicament) étiquetée sous la cible de référence : source et étiquette."""
    out = []
    for d in DRUGS:
        lab = cr.etiquettes_cible(rows, d, "reference")
        src = cr.source_cible_par_ligne(rows, d, "reference")
        ok = lab.isin(["R", "S"])
        out.append(pd.DataFrame({"strain_id": rows.loc[ok, "strain_id"].values, "drug": d,
                                 "etiquette": lab[ok].values, "source": src[ok].values}))
    return pd.concat(out, ignore_index=True)


def parts(df, cle, drug="drug"):
    """Effectifs R et S par (médicament, `cle`) et part de chaque ligne dans les R, les S et l'ensemble."""
    t = df.groupby([drug, cle, "etiquette"]).size().unstack("etiquette", fill_value=0)
    for c in ("R", "S"):
        if c not in t:
            t[c] = 0
    t = t[["R", "S"]].rename(columns={"R": "n_R", "S": "n_S"})
    t["n"] = t["n_R"] + t["n_S"]
    for c, nom in (("n_R", "part_R"), ("n_S", "part_S"), ("n", "part")):
        t[nom] = (t[c] / t.groupby(level=drug)[c].transform("sum")).round(4)
    return t.reset_index().sort_values([drug, "n"], ascending=[True, False])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--etude", default="PRJNA1237251", help="BioProject de l'étude de surveillance citée")
    a = ap.parse_args()
    rows = cohorte()
    p = paires(rows)
    ind = [("souches_cohorte", len(rows)), ("paires_etiquetees_11_medicaments", len(p)),
           ("souches_phenotypees", int(p["strain_id"].nunique()))]

    # 1. source de table
    t_src = parts(p, "source")
    t_src.to_csv(M / "p714_sources_par_medicament.tsv", sep="\t", index=False)
    glob = p.groupby("source").size().sort_values(ascending=False)
    for s, n in glob.items():
        ind.append((f"paires_source_{s}", int(n)))
    ind.append(("sources_de_table_representees", int(glob.size)))

    # 2. source de DST et cohorte pour les paires CRyPTIC élargi
    v3 = pd.read_csv(paths.DATA / "sources" / "cryptic" / "cryptic_v3_dst_by_sra.csv",
                     usecols=["strain_id", "drug_code", "dst_source", "cohort"]).drop_duplicates(["strain_id", "drug_code"])
    pv = p[p["source"] == "CRyPTIC_v3"].copy()
    pv["drug_code"] = pv["drug"].map(CODES)
    pv = pv.merge(v3, on=["strain_id", "drug_code"], how="left")
    pv["dst_source"] = pv["dst_source"].fillna("(non retrouvée)")
    pv["cohorte_cryptic"] = pv["cohort"].fillna("(non retrouvée)")
    t_dst = (pv.groupby(["drug", "dst_source", "cohorte_cryptic", "etiquette"]).size().unstack("etiquette", fill_value=0)
               .rename(columns={"R": "n_R", "S": "n_S"}).reset_index())
    t_dst["n"] = t_dst["n_R"] + t_dst["n_S"]
    t_dst.sort_values(["drug", "n"], ascending=[True, False]).to_csv(M / "p714_dst_cryptic_v3.tsv", sep="\t", index=False)
    ind.append(("paires_cryptic_v3", len(pv)))
    ind.append(("dst_sources_retrouvees_dans_les_paires_cryptic_v3", int((pv["dst_source"] != "(non retrouvée)").groupby(pv["dst_source"]).any().sum())))
    for s, n in pv.groupby("dst_source").size().sort_values(ascending=False).items():
        ind.append((f"paires_cryptic_v3_dst_{s}", int(n)))

    # 3. BioProject
    reg = pd.read_csv(paths.CATALOGUE_DIR / "etudes_ena.tsv", sep="\t", usecols=["strain_id", "bioproject"])
    reg = reg.drop_duplicates("strain_id").set_index("strain_id")["bioproject"]
    p["run"] = p["strain_id"].str.split(".").str[0]
    p["bioproject"] = p["run"].map(reg).fillna("(inconnu)")
    t_bp = parts(p, "bioproject")
    t_bp.to_csv(M / "p714_etudes_par_medicament.tsv", sep="\t", index=False)
    ind.append(("bioprojects_distincts", int(p.loc[p["bioproject"] != "(inconnu)", "bioproject"].nunique())))
    ind.append(("paires_sans_bioproject", int((p["bioproject"] == "(inconnu)").sum())))
    for d in DRUGS:
        sous = t_bp[t_bp["drug"] == d]
        for etat, col in (("S", "part_S"), ("R", "part_R")):
            top = sous.sort_values(col, ascending=False).iloc[0]
            ind.append((f"{CODES[d]}_plus_grande_etude_{etat}", f"{top['bioproject']} ({top[col]:.4f}, n={int(top['n_' + etat])})"))
        cible = sous[sous["bioproject"] == a.etude]
        ind.append((f"{CODES[d]}_part_{a.etude}_S", float(cible["part_S"].iloc[0]) if len(cible) else 0.0))
        ind.append((f"{CODES[d]}_n_{a.etude}_S", int(cible["n_S"].iloc[0]) if len(cible) else 0))
        ind.append((f"{CODES[d]}_part_{a.etude}_R", float(cible["part_R"].iloc[0]) if len(cible) else 0.0))

    # 4. contrôles
    eff = pd.read_csv(M / "p76_effectifs_blocs.tsv", sep="\t").set_index("indicateur")["valeur"]
    att = int(eff["souches_phenotypees_reference_11_medicaments"])
    ind.append(("controle_souches_phenotypees_attendues_P76", att))
    ind.append(("controle_souches_phenotypees_ok", int(att == int(p["strain_id"].nunique()))))
    ref = pd.read_csv(M / "p74_angles_morts_blocs_P15.tsv", sep="\t")
    ref = ref[(ref["cible"] == "reference") & (ref["etat"] == "complet")].set_index("drug")
    nb = p.groupby(["drug", "etiquette"]).size().unstack(fill_value=0)
    ecarts = {d: (int(nb.loc[d, "R"]) - int(ref.loc[d, "n_R"]), int(nb.loc[d, "S"]) - int(ref.loc[d, "n_S"])) for d in DRUGS}
    for d, (er, es) in ecarts.items():
        ind.append((f"controle_ecart_n_R_n_S_{CODES[d]}", f"{er},{es}"))
    ind.append(("controle_effectifs_par_medicament_ok", int(all(v == (0, 0) for v in ecarts.values()))))
    pd.DataFrame(ind, columns=["indicateur", "valeur"]).to_csv(M / "p714_synthese.tsv", sep="\t", index=False)
    print(pd.DataFrame(ind, columns=["indicateur", "valeur"]).to_string(index=False))


if __name__ == "__main__":
    main()
