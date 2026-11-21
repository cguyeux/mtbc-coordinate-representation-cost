#!/usr/bin/env python3
"""Neutraliser la génération du pipeline entre SPDI anciens (pangénome, bdd/actuelle) et SPDI tblearn (P8.8).

Objet      : les SPDI de tblearn ne sont pas ceux des deux sources anciennes. Mesuré en P8.2 puis ici :
             (1) le pipeline ancien fusionne en un seul variant (« MNV », substitution de même longueur
             avec plusieurs positions différentes, ex. 3247850:GCCCGG:ACTCGC) des SNP voisins que tblearn
             garde séparés ; (2) tblearn porte environ 6 % de SNP en plus, surtout sur quelques centaines de
             sites récurrents ; (3) quelques indels complexes n'ont pas d'équivalent. Une matrice qui mélange
             les deux sans précaution rend la GÉNÉRATION prédictive. Ce module (a) ramène les deux
             représentations à un même jeu (fermeture par composantes pour l'ancien, recomposition des MNV
             connus pour tblearn), (b) calcule sur des souches appelées par les deux sources la liste des
             variants dont l'appel dépend de la source (le masque), (c) le valide par un classifieur de
             génération sur des souches tenues hors du calcul du masque.
Entrées    : voir les sous-commandes (`controles`, `masque`, `valider`) ; les appels anciens viennent de
             `p88_bdd_extrait.tsv.gz` (extraction de bdd/actuelle) ou du pangénome, les appels tblearn de
             `p82_fetch_spdi_tblearn.py`.
Sorties    : `results/model/p88_masque_generation.tsv` (un variant par ligne, avec ses effectifs),
             `p88_validation_generation.tsv`.
Réutilisable : les fonctions `est_mnv`, `composantes`, `canon`, `Conjonctions` valent pour toute référence
             bactérienne dont les SPDI sont `<accession>:<pos>:<ref>:<alt>` ; Python 3.8 (tourne sur `mp`).
Projet     : Resistance_antibio (P8.8)
Date       : 2026-10-20
"""
import gzip
import sys
from collections import Counter, defaultdict


def lire_spdi(texte):
    """SPDI séparés par des virgules (tblearn) ou des blancs (spdi.txt) -> ensemble."""
    return set(t for t in texte.replace(",", " ").split() if t)


def decouper(spdi):
    """`acc:pos:ref:alt` -> (acc, pos, ref, alt) ; None si le texte n'a pas cette forme."""
    p = spdi.split(":")
    if len(p) != 4 or not p[1].isdigit():
        return None
    return p[0], int(p[1]), p[2], p[3]


def est_mnv(spdi):
    """Substitution de même longueur (>= 2) portant au moins deux positions différentes."""
    d = decouper(spdi)
    if d is None:
        return False
    _, _, ref, alt = d
    return len(ref) == len(alt) and len(ref) >= 2 and sum(a != b for a, b in zip(ref, alt)) >= 2


def composantes(spdi):
    """SNP composants d'un MNV (positions où ref et alt diffèrent) ; liste vide si ce n'est pas un MNV."""
    d = decouper(spdi)
    if d is None or not est_mnv(spdi):
        return []
    acc, pos, ref, alt = d
    return ["%s:%d:%s:%s" % (acc, pos + i, a, b) for i, (a, b) in enumerate(zip(ref, alt)) if a != b]


def canon(variants):
    """Représentation CANONIQUE d'un jeu de variants, la même pour toute source : chaque bloc MNV est
    remplacé par ses SNP composants. Les deux générations ne découpent pas les blocs aux mêmes bornes
    (mesuré sur 928 souches : tblearn fusionne aussi des SNP que l'ancien pipeline garde séparés), donc
    aucun bloc ne peut servir de variable ; les SNP composants, eux, sont communs."""
    sortie = set()
    for v in variants:
        comp = composantes(v)
        if comp:
            sortie.update(comp)
        else:
            sortie.add(v)
    return sortie


class Conjonctions:
    """MNV du catalogue, reconstruits à partir des SNP d'un jeu canonique.

    Un MNV catalogué est déclaré présent dans une souche quand TOUTES ses composantes y sont : dans un
    génome haploïde, des SNP présents sont en cis. La colonne est donc définie de la même façon pour
    toutes les souches, quelle que soit la façon dont leur pipeline a groupé les SNP voisins."""

    def __init__(self, mnvs=()):
        self.par_tete = defaultdict(list)
        self.n = 0
        for m in set(mnvs):
            comp = composantes(m)
            if comp:
                self.par_tete[comp[0]].append((m, tuple(comp)))
                self.n += 1

    def ajouter(self, variants):
        """Jeu canonique -> jeu canonique + MNV catalogués dont toutes les composantes sont présentes."""
        sortie = set(variants)
        for v in variants:
            for mnv, comp in self.par_tete.get(v, ()):
                if all(c in variants for c in comp):
                    sortie.add(mnv)
        return sortie


def normaliser(variants, conj=None):
    """Jeu brut (ancien ou tblearn) -> jeu canonique, avec les MNV catalogués si `conj` est donné."""
    c = canon(variants)
    return conj.ajouter(c) if conj is not None else c


def lire_tsv_spdi(chemin, col_id=0, col_spdi=5, entete="accession"):
    """TSV gzip d'une ligne par souche (format de `p82_fetch_spdi_tblearn.py`) -> {id: ensemble}."""
    sortie = {}
    with gzip.open(str(chemin), "rt") as f:
        for ligne in f:
            c = ligne.rstrip("\n").split("\t")
            if len(c) <= col_spdi or c[col_id] == entete:
                continue
            sortie[c[col_id]] = lire_spdi(c[col_spdi])
    return sortie


def iter_tsv_source(chemin):
    """Lecture en flux du TSV de `ecrire_tsv_spdi` : (id, source, ensemble), sans tout garder en mémoire."""
    with gzip.open(str(chemin), "rt") as f:
        for ligne in f:
            c = ligne.rstrip("\n").split("\t")
            if len(c) >= 4:
                yield c[0], c[1], lire_spdi(c[3])


def iter_tsv_spdi(chemin, col_id=0, col_spdi=5, entete="accession"):
    """Lecture en flux du TSV de `p82_fetch_spdi_tblearn.py` : (id, ensemble)."""
    with gzip.open(str(chemin), "rt") as f:
        for ligne in f:
            c = ligne.rstrip("\n").split("\t")
            if len(c) > col_spdi and c[col_id] != entete:
                yield c[col_id], lire_spdi(c[col_spdi])


def fusionner_composites(par_run, composites):
    """Isolat CRyPTIC à plusieurs runs : UNION des variants des runs lus (un appel dans un run suffit ;
    l'absence d'appel dans l'autre n'est pas un allèle de référence). `composites` : {isolat: [runs]}."""
    sortie = {}
    for iso, runs in composites.items():
        u = set()
        lus = 0
        for r in runs:
            if r in par_run:
                u |= par_run[r]
                lus += 1
        if lus:
            sortie[iso] = u
    return sortie


def effectifs_discordance(paires):
    """`paires` : itérable de (jeu_ancien_normalise, jeu_tblearn_normalise) pour des MEMES souches.
    Rend trois Counter : variant -> nb de souches où il est appelé par les deux, seulement par l'ancien,
    seulement par tblearn."""
    les_deux, ancien, tbl = Counter(), Counter(), Counter()
    n = 0
    for a, t in paires:
        n += 1
        les_deux.update(a & t)
        ancien.update(a - t)
        tbl.update(t - a)
    return n, les_deux, ancien, tbl


def variants_masques(les_deux, ancien, tbl, part_min=0.02, discordants_min=2):
    """Seuils par défaut (0,02 ; 2) choisis sur une grille de 12 réglages par l'AUC de génération hors échantillon
    (masque sur la moitié A, jugement sur B puis l'inverse : 0,529 et 0,537 contre 0,572 et 0,581 pour 0,10 ; 3).
    Variants dont l'appel dépend de la source : discordants chez au moins `discordants_min` souches ET
    qui forment au moins `part_min` des souches où le variant est appelé par l'une ou l'autre source.
    Le critère est RELATIF (part des appels), pas absolu : un variant rare mais appelé par une seule
    source est un marqueur de génération même s'il ne concerne que 1 % des souches."""
    sortie = {}
    for v in set(ancien) | set(tbl):
        d = ancien.get(v, 0) + tbl.get(v, 0)
        tot = d + les_deux.get(v, 0)
        if d >= discordants_min and d / tot >= part_min:
            sortie[v] = (les_deux.get(v, 0), ancien.get(v, 0), tbl.get(v, 0))
    return sortie


def ecrire_tsv_spdi(chemin, lignes):
    """lignes : itérable de (id, source, ensemble) -> TSV gzip id, source, n, liste."""
    with gzip.open(str(chemin), "wt") as f:
        for i, src, ens in lignes:
            f.write("%s\t%s\t%d\t%s\n" % (i, src, len(ens), ",".join(sorted(ens))))


def lire_tsv_source(chemin):
    """TSV de `ecrire_tsv_spdi` -> ({id: ensemble}, {id: source})."""
    ens, src = {}, {}
    with gzip.open(str(chemin), "rt") as f:
        for ligne in f:
            c = ligne.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            ens[c[0]] = lire_spdi(c[3])
            src[c[0]] = c[1]
    return ens, src


def charger_catalogue(chemin):
    """Catalogue consolidé -> (Conjonctions des MNV catalogués, ensemble des SPDI R-associated)."""
    import pandas as pd
    cat = pd.read_csv(str(chemin), sep="\t").fillna("")
    conj = Conjonctions(s for s in cat["spdi"].unique() if s and est_mnv(s))
    cat_r = set(s for s in cat.loc[cat["call"] == "R-associated", "spdi"] if s)
    return conj, cat_r


def cmd_extraire_bdd(args):
    """Extrait de bdd/actuelle les spdi.txt des souches demandées (local : bdd n'est pas sur mp)."""
    import subprocess
    from pathlib import Path
    import pandas as pd
    racine = Path(args.bdd)
    sortie = subprocess.run(["find", str(racine), "-path", "*/NC_000962.3/spdi.txt"],
                            capture_output=True, text=True).stdout
    idx = {}
    for p in sortie.splitlines():
        parts = p.split("/")
        if len(parts) >= 3:
            idx[parts[-3]] = p
    voulues = set(pd.read_csv(args.phenotypes, sep="\t", usecols=["strain_id"])["strain_id"].dropna())
    for extra in args.en_plus or []:
        voulues |= set(l.strip() for l in open(extra) if l.strip())
    choisies = sorted(s for s in voulues if s in idx)
    ecrire_tsv_spdi(args.sortie, ((s, "bdd", lire_spdi(Path(idx[s]).read_text())) for s in choisies))
    print("bdd/actuelle : %d souches voulues, %d extraites -> %s" % (len(voulues), len(choisies), args.sortie))


def cmd_selectionner(args):
    """Contrôles : échantillon de souches de la matrice (pangénome ou bdd), stratifié par lignée Coll
    (allocation proportionnelle à la racine des effectifs, plancher = tout si moins de 30)."""
    import pickle
    import random
    import pandas as pd
    bdd, _ = lire_tsv_source(args.bdd_extrait)
    with open(args.pan, "rb") as f:
        pan = pickle.load(f)
    deja = set()
    for fic in args.deja or []:
        deja |= set(l.strip() for l in open(fic) if l.strip())
    phe = set(pd.read_csv(args.phenotypes, sep="\t", usecols=["strain_id"])["strain_id"].dropna())
    lin = pd.read_csv(args.lignees, usecols=["strain_name", "lineage_level_1"], dtype=str)
    lin = lin.drop_duplicates("strain_name").set_index("strain_name")["lineage_level_1"]
    pool = sorted(s for s in phe if (s in pan or s in bdd) and s not in deja)
    par_l = defaultdict(list)
    for s in pool:
        par_l[str(lin.get(s, "NA"))].append(s)
    poids = {l: len(v) ** 0.5 for l, v in par_l.items()}
    z = sum(poids.values())
    rng = random.Random(args.graine)
    choix = []
    for l, v in sorted(par_l.items()):
        k = len(v) if len(v) <= 30 else min(len(v), int(round(args.n * poids[l] / z)))
        choix += rng.sample(v, k)
    anc = []
    for s in choix:
        if s in pan:
            anc.append((s, "pan", set(pan[s])))
        else:
            anc.append((s, "bdd", bdd[s]))
    # les contrôles déjà tirés en P8.2 ont leur jeu ancien dans bdd
    for s in sorted(deja):
        if s in bdd:
            anc.append((s, "bdd", bdd[s]))
    ecrire_tsv_spdi(args.sortie_ancien, anc)
    with open(args.sortie_ids, "w") as f:
        f.write("\n".join(choix) + "\n")
    par_src = Counter(a[1] for a in anc)
    # diagnostic : pangénome et bdd/actuelle donnent-ils le même jeu pour une même souche ?
    ecarts = [(len(canon(pan[s]) - canon(bdd[s])), len(canon(bdd[s]) - canon(pan[s])))
              for s in choix if s in pan and s in bdd]
    if ecarts:
        print("pangénome vs bdd (%d souches en commun parmi les contrôles) : pan-seul %.1f, bdd-seul %.1f par souche"
              % (len(ecarts), sum(e[0] for e in ecarts) / len(ecarts), sum(e[1] for e in ecarts) / len(ecarts)))
    print("contrôles : %d tirés (+ %d déjà tirés), sources %s" % (len(choix), len(deja), dict(par_src)))
    print("par lignée :", {l: sum(1 for s in choix if str(lin.get(s, "NA")) == l) for l in par_l})


def charger_controles(args, conj):
    """Appariement ancien/tblearn des mêmes souches -> {souche: (ancien_norm, tbl_norm, source_ancienne)}."""
    anc, src = lire_tsv_source(args.ancien)
    tbl = {}
    for fic in args.tblearn:
        tbl.update(lire_tsv_spdi(fic))
    return {s: (normaliser(anc[s], conj), normaliser(tbl[s], conj), src[s]) for s in anc if s in tbl}


def moitie(ctr, nom):
    """Sous-ensemble déterministe des contrôles : `A` = rangs pairs des souches triées, `B` = rangs impairs,
    `toutes` = tout. Le masque se calcule sur A, se juge sur B (jamais vu par le calcul du masque)."""
    souches = sorted(ctr)
    choix = {"A": souches[0::2], "B": souches[1::2], "toutes": souches}[nom]
    return {s: ctr[s] for s in choix}


def cmd_masque(args):
    """Masque de génération et validation hors échantillon (5 plis par souche)."""
    import random
    conj, cat_r = charger_catalogue(args.catalogue)
    ctr = moitie(charger_controles(args, conj), args.moitie)
    souches = sorted(ctr)
    n, deux, anc, tbl = effectifs_discordance([(ctr[s][0], ctr[s][1]) for s in souches])
    m = variants_masques(deux, anc, tbl, args.part_min, args.discordants_min)
    with open(args.sortie, "w") as f:
        f.write("variant\tdeux\tancien_seul\ttblearn_seul\tcatalogue_R\n")
        for v, (d, a, t) in sorted(m.items(), key=lambda kv: -(kv[1][1] + kv[1][2])):
            f.write("%s\t%d\t%d\t%d\t%d\n" % (v, d, a, t, int(v in cat_r)))
    # déterminants R catalogués : discordance mesurée (jamais masqués)
    det = [(v, deux.get(v, 0), anc.get(v, 0), tbl.get(v, 0)) for v in cat_r
           if deux.get(v, 0) + anc.get(v, 0) + tbl.get(v, 0) > 0]
    with open(args.sortie + ".determinants.tsv", "w") as f:
        f.write("variant\tdeux\tancien_seul\ttblearn_seul\n")
        for v, d, a, t in sorted(det, key=lambda x: -(x[2] + x[3])):
            f.write("%s\t%d\t%d\t%d\n" % (v, d, a, t))
    print("contrôles : %d souches ; masque : %d variants ; déterminants R observés %d dont discordants %d"
          % (n, len(m), len(det), sum(1 for x in det if x[2] + x[3] > 0)))
    # validation hors échantillon : le masque d'un pli est calculé SANS ses souches
    rng = random.Random(args.graine)
    rng.shuffle(souches)
    k = 5
    plis = [souches[i::k] for i in range(k)]
    lignes = []
    brut_ds = {"ancien_seul": 0, "tblearn_seul": 0, "n": 0}
    hors = {"ancien_seul": 0, "tblearn_seul": 0, "n": 0}
    pred = []
    for i in range(k):
        test = set(plis[i])
        ent = [s for s in souches if s not in test]
        n2, d2, a2, t2 = effectifs_discordance([(ctr[s][0], ctr[s][1]) for s in ent])
        m2 = set(variants_masques(d2, a2, t2, args.part_min, args.discordants_min)) - cat_r
        for s in test:
            a, t, _ = ctr[s]
            hors["n"] += 1
            hors["ancien_seul"] += len((a - t) - m2)
            hors["tblearn_seul"] += len((t - a) - m2)
            brut_ds["n"] += 1
            brut_ds["ancien_seul"] += len(a - t)
            brut_ds["tblearn_seul"] += len(t - a)
        pred.append((ent, plis[i], m2))
    lignes.append(("discordance_par_souche_canonique", brut_ds["ancien_seul"] / brut_ds["n"],
                   brut_ds["tblearn_seul"] / brut_ds["n"]))
    lignes.append(("discordance_par_souche_masque_hors_pli", hors["ancien_seul"] / hors["n"],
                   hors["tblearn_seul"] / hors["n"]))
    # classifieur de génération : distingue-t-il l'appel ancien de l'appel tblearn de la MÊME souche ?
    try:
        import numpy as np
        from scipy import sparse
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import roc_auc_score
    except ImportError as e:
        print("classifieur non exécuté :", e)
        np = None
    if np is not None:
        for nom, appliquer in (("auc_generation_canonique", False), ("auc_generation_masque_hors_pli", True)):
            y_all, p_all = [], []
            for ent, test, m2 in pred:
                def jeux(liste):
                    X, y = [], []
                    for s in liste:
                        a, t, _ = ctr[s]
                        for ens, lab in ((a, 0), (t, 1)):
                            X.append(ens - m2 if appliquer else ens)
                            y.append(lab)
                    return X, y
                Xe, ye = jeux(ent)
                Xt, yt = jeux(test)
                voc = {}
                for ens in Xe:
                    for v in ens:
                        voc.setdefault(v, len(voc))
                def mat(X):
                    ind, ptr = [], [0]
                    for ens in X:
                        ind.extend(voc[v] for v in ens if v in voc)
                        ptr.append(len(ind))
                    return sparse.csr_matrix((np.ones(len(ind), dtype=np.float32), np.array(ind), np.array(ptr)),
                                             shape=(len(X), len(voc)))
                clf = LogisticRegression(C=0.5, max_iter=300, solver="liblinear")
                clf.fit(mat(Xe), ye)
                p_all.extend(clf.predict_proba(mat(Xt))[:, 1])
                y_all.extend(yt)
            lignes.append((nom, float(roc_auc_score(y_all, p_all)), float("nan")))
    with open(args.sortie_validation, "w") as f:
        f.write("mesure\tvaleur_a\tvaleur_b\n")
        for nom, a, b in lignes:
            f.write("%s\t%.4f\t%.4f\n" % (nom, a, b))
            print("  %-42s %.4f  %.4f" % (nom, a, b))


def auc_genre_hors_pli(ctr, souches, features, masque, plis=5, graine=88):
    """AUC hors pli d'un classifieur qui distingue l'appel ancien de l'appel tblearn d'une MÊME souche, sur les
    seules colonnes `features` (celles que la matrice porte vraiment), après retrait de `masque`."""
    import random
    import numpy as np
    from scipy import sparse
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    col = {v: j for j, v in enumerate(features)}
    rng = random.Random(graine)
    ordre = list(souches)
    rng.shuffle(ordre)

    def mat(liste):
        ind, ptr, y = [], [0], []
        for s in liste:
            a, t, _ = ctr[s]
            for ens, lab in ((a, 0), (t, 1)):
                ind.extend(col[v] for v in ens if v in col and v not in masque)
                ptr.append(len(ind))
                y.append(lab)
        X = sparse.csr_matrix((np.ones(len(ind), dtype=np.float32), np.array(ind, dtype=np.int64),
                               np.array(ptr, dtype=np.int64)), shape=(len(y), len(features)))
        return X, np.array(y)

    y_all, p_all = [], []
    for i in range(plis):
        test = ordre[i::plis]
        app = [s for s in ordre if s not in set(test)]
        Xe, ye = mat(app)
        Xt, yt = mat(test)
        clf = LogisticRegression(C=0.5, max_iter=300, solver="liblinear")
        clf.fit(Xe, ye)
        p_all.extend(clf.predict_proba(Xt)[:, 1])
        y_all.extend(yt)
    return float(roc_auc_score(y_all, p_all))


def cmd_validation_matrice(args):
    """Neutralisation de la génération au niveau des COLONNES de la matrice : paires appel ancien / appel
    tblearn des mêmes souches, restreintes aux features finales ; modèle nul = AUC 0,5."""
    conj, cat_r = charger_catalogue(args.catalogue)
    ctr = moitie(charger_controles(args, conj), args.moitie)
    features = open(args.features).read().split("\n")
    fset = set(features)
    masque = set()
    if args.masque:
        import pandas as pd
        masque = set(pd.read_csv(args.masque, sep="\t")["variant"]) - cat_r
    souches = sorted(ctr)
    diff = [len(((a ^ t) & fset) - masque) for a, t, _ in (ctr[s] for s in souches)]
    det = Counter()
    for s in souches:
        a, t, _ = ctr[s]
        for v in ((a ^ t) & fset & cat_r):
            det[v] += 1
    auc = auc_genre_hors_pli(ctr, souches, features, masque)
    lignes = [("controles", len(souches)), ("features", len(features)), ("masque_applique", len(masque)),
              ("colonnes_discordantes_par_souche", sum(diff) / len(diff)),
              ("souches_avec_une_colonne_discordante", sum(1 for d in diff if d) / len(diff)),
              ("determinants_R_discordants", len(det)),
              ("part_max_souches_determinant_discordant", (max(det.values()) / len(souches)) if det else 0.0),
              ("auc_generation_sur_features", auc)]
    with open(args.sortie_validation, "w") as f:
        f.write("mesure\tvaleur\n")
        for nom, v in lignes:
            f.write("%s\t%.4f\n" % (nom, v))
            print("  %-42s %.4f" % (nom, v))
    for v, k in det.most_common(8):
        print("  déterminant R discordant", v, k)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extraire-bdd")
    e.add_argument("--bdd", required=True)
    e.add_argument("--phenotypes", required=True)
    e.add_argument("--en-plus", nargs="*")
    e.add_argument("--sortie", required=True)
    e.set_defaults(fn=cmd_extraire_bdd)
    s = sub.add_parser("selectionner")
    s.add_argument("--pan", required=True)
    s.add_argument("--bdd-extrait", required=True)
    s.add_argument("--phenotypes", required=True)
    s.add_argument("--lignees", required=True)
    s.add_argument("--deja", nargs="*")
    s.add_argument("--n", type=int, default=3000)
    s.add_argument("--graine", type=int, default=88)
    s.add_argument("--sortie-ids", required=True)
    s.add_argument("--sortie-ancien", required=True)
    s.set_defaults(fn=cmd_selectionner)
    m = sub.add_parser("masque")
    m.add_argument("--ancien", required=True)
    m.add_argument("--tblearn", nargs="+", required=True)
    m.add_argument("--catalogue", required=True)
    m.add_argument("--part-min", type=float, default=0.02)
    m.add_argument("--discordants-min", type=int, default=2)
    m.add_argument("--graine", type=int, default=88)
    m.add_argument("--moitie", choices=["A", "B", "toutes"], default="toutes")
    m.add_argument("--sortie", required=True)
    m.add_argument("--sortie-validation", required=True)
    m.set_defaults(fn=cmd_masque)
    v = sub.add_parser("validation-matrice")
    v.add_argument("--ancien", required=True)
    v.add_argument("--tblearn", nargs="+", required=True)
    v.add_argument("--catalogue", required=True)
    v.add_argument("--features", required=True)
    v.add_argument("--masque")
    v.add_argument("--moitie", choices=["A", "B", "toutes"], default="toutes")
    v.add_argument("--sortie-validation", required=True)
    v.set_defaults(fn=cmd_validation_matrice)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
