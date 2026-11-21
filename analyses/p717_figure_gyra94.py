#!/usr/bin/env python3
"""P7.17 — Figure 1 d'A1 : anatomie du locus gyrA codons 94 et 95, générée depuis les sources.

Objet      : donner à A1 sa première figure, celle du mécanisme que le manuscrit décrit en prose
             (action 5 du verdict amont P7.10, préoccupations 1 et 6 de la review du 2026-11-21).
             Elle montre, position par position et sur une seule échelle de coordonnées : la
             séquence de référence et les deux codons concernés ; ce que l'appeleur émet pour une
             souche qui porte un variant du codon 94 sur fond Ser95Thr (un bloc unique, de quatre
             ou cinq bases) ; et ce que le fichier de coordonnées de l'OMS liste pour chacun des
             HUIT variants protéiques d'Asp94, cinq recevant une entrée combinée (quatre bases
             pour Gly, Val, Ala ; CINQ bases pour Phe et Leu) et trois n'en recevant aucune
             (Asn, Tyr, His), à l'étendue et dans la forme exactes des entrées de Phe et Leu.
             L'argument de la §3.3 du manuscrit est une IDENTITÉ de position et de forme, non une
             analogie : la figure la rend visible en alignant les barres.
             Le TikZ est GÉNÉRÉ et non écrit à la main, parce que c'est la recopie d'un compte
             d'une passe à l'autre qui avait propagé la faute corrigée par P7.16 : aucune base,
             aucune chaîne et aucun effectif de cette figure n'est saisi ici, tout est lu aux
             fichiers d'entrée.
Entrées    : data/sources/who_catalogue/WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt (positions
             1-based de la convention VCF), data/fastas/NC_000962.3.fasta (séquence de référence),
             results/model/p77_C1_residuel.tsv (effectifs des trois haplotypes du résiduel).
Sorties    : article/figures/fig_gyra94.tex (standalone TikZ, à compiler en PDF vectoriel) et
             results/model/p717_figure_gyra94.tsv (toute valeur imprimée par la figure, pour
             épinglage par test_p78_chiffres_manuscrit.py).
Réutilisable : oui pour un autre codon (--gene, --codon, --voisin) de ce fichier de coordonnées.
Projet     : Resistance_antibio (P7.17)
Date       : 2026-11-21
Usage      : python3 analyses/p717_figure_gyra94.py
             make -C article figures/fig_gyra94.pdf
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

RACINE = Path(paths.ROOT)
SOURCE = Path(paths.SRC_WHO) / "WHO-UCN-TB-2023.7-eng_genomic_coordinates.txt"
FASTA = RACINE / "data" / "fastas" / "NC_000962.3.fasta"
RESIDUEL = RACINE / "results" / "model" / "p77_C1_residuel.tsv"
FIGURE = RACINE / "article" / "figures" / "fig_gyra94.tex"
TSV = RACINE / "results" / "model" / "p717_figure_gyra94.tsv"

# Code à trois lettres -> une lettre, pour l'étiquette des codons de la référence.
AA3 = {"Asp": "D", "Ser": "S", "Asn": "N", "Tyr": "Y", "His": "H",
       "Gly": "G", "Val": "V", "Ala": "A", "Phe": "F", "Leu": "L", "Thr": "T"}

# Palette Okabe-Ito, lisible en vision dichromate et en niveaux de gris.
COULEURS = {
    "ref": "595959", "codon94": "0072B2", "codon95": "7B52AB",
    "appel": "E69F00", "liste": "0072B2", "combine": "009E73", "absent": "D55E00",
}


def sequence_reference(fasta, debut, fin):
    """Rend les bases 1-based [debut, fin] du premier enregistrement du FASTA."""
    morceaux = []
    with open(fasta, encoding="ascii") as f:
        for ligne in f:
            if ligne.startswith(">"):
                if morceaux:
                    break
                continue
            morceaux.append(ligne.strip())
            if sum(len(m) for m in morceaux) >= fin:
                break
    seq = "".join(morceaux)
    if len(seq) < fin:
        raise SystemExit(f"séquence trop courte ({len(seq)} bases) pour atteindre {fin}")
    return seq[debut - 1:fin]


def representations(source, gene, codon):
    """Rend {acide aminé: [(position, ref, alt), ...]} pour gene_p.<X><codon><Y>."""
    motif = re.compile(rf"^{re.escape(gene)}_p\.([A-Z][a-z]{{2}}){codon}([A-Z][a-z]{{2}})$")
    par_aa, origine = defaultdict(list), set()
    with open(source, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            m = motif.match((r.get("variant") or "").strip())
            if not m:
                continue
            origine.add(m.group(1))
            par_aa[m.group(2)].append((int(r["position"]),
                                       (r.get("reference_nucleotide") or "").strip(),
                                       (r.get("alternative_nucleotide") or "").strip()))
    if len(origine) != 1:
        raise SystemExit(f"codon {codon} de {gene} : acide aminé de référence ambigu {origine}")
    return origine.pop(), par_aa


def groupes(reps):
    """Groupe les représentations d'un variant par (position, longueur) et trie du plus long
    au plus court, les alt de chaque groupe étant rendus dans l'ordre du fichier."""
    par_cle = defaultdict(list)
    for position, ref, alt in reps:
        par_cle[(position, ref)].append(alt)
    sortie = []
    for (position, ref), alts in par_cle.items():
        sortie.append({"position": position, "ref": ref, "alts": alts, "longueur": len(ref)})
    sortie.sort(key=lambda g: (-g["longueur"], g["position"]))
    return sortie


def bloc_attendu(g, etendue, locus, debut_locus, pos_fond, alt_fond):
    """Rend (position, ref, alt) du bloc combiné qu'un variant AURAIT avec le marqueur de fond,
    reconstruit depuis sa substitution simple et la séquence de référence. Rend None quand le
    fichier ne donne pas de substitution simple pour ce variant."""
    simple = min((x for x in g if x["longueur"] == 1), key=lambda x: x["position"], default=None)
    if simple is None:
        return None
    debut, largeur = etendue
    i0 = debut - debut_locus
    ref = locus[i0:i0 + largeur]
    alt = list(ref)
    alt[simple["position"] - debut] = simple["alts"][0]
    alt[pos_fond - debut] = alt_fond
    return debut, ref, "".join(alt)


def effectifs_residuel(chemin):
    """Rend {alt du bloc de 5 bases: {souches, par médicament (R, S)}} depuis p77_C1_residuel."""
    par_alt = {}
    with open(chemin, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            alt = r["spdi"].split(":")[-1]
            e = par_alt.setdefault(alt, {"souches": int(r["n_souches"]), "drogues": {}})
            if e["souches"] != int(r["n_souches"]):
                raise SystemExit(f"effectif incohérent pour {alt} entre médicaments")
            e["drogues"][r["drug"]] = (int(r["n_R"]), int(r["n_S"]))
    return par_alt


def muter(ref_locus, debut_locus, position, ref, alt):
    """Rend la chaîne alt alignée sur le locus, et les index des bases qui changent."""
    i = position - debut_locus
    changees = [i + k for k in range(len(ref)) if k < len(alt) and ref[k] != alt[k]]
    if ref_locus[i:i + len(ref)] != ref:
        raise SystemExit(f"allèle de référence {ref} absent du locus en {position}")
    return i, changees


class Dessin:
    """Accumulateur de lignes TikZ, avec une grille d'une base de large par colonne."""

    def __init__(self, debut, fin, dx=0.78, dy=0.40):
        self.debut, self.fin, self.dx, self.dy = debut, fin, dx, dy
        self.n = fin - debut + 1
        self.lignes, self.y = [], 0.0
        self.y_min_cellule = 0.0   # bas de la dernière rangée de cellules dessinée

    def x(self, index):
        return index * self.dx

    def descendre(self, pas=1.0):
        self.y -= pas * self.dy

    def texte(self, index, texte, ancre="center", decalage=0.0, taille="", couleur=None, dy=0.0):
        c = f"text={{{couleur}}}," if couleur else ""
        # L'espace après la commande de taille est obligatoire : « \footnotesizeC » serait lu
        # comme une commande inconnue et le premier mot disparaîtrait silencieusement.
        t = f"{taille} " if taille else ""
        self.lignes.append(
            f"  \\node[anchor={ancre},{c}inner sep=1pt] at "
            f"({self.x(index) + decalage:.3f},{self.y + dy:.3f}) {{{t}{texte}}};")

    def cellules(self, bases, index0, style, couleur, changees=(), taille="\\footnotesize"):
        """Dessine une base par cellule, les bases changées en gras et en couleur."""
        self.y_min_cellule = self.y - self.dy / 2
        for k, base in enumerate(bases):
            i = index0 + k
            chg = i in changees
            fond = f"fill={couleur}!{'22' if chg else '7'}" if couleur else "fill=none"
            self.lignes.append(
                f"  \\draw[{style},{fond}] ({self.x(i) - self.dx / 2:.3f},{self.y - self.dy / 2:.3f}) "
                f"rectangle ({self.x(i) + self.dx / 2:.3f},{self.y + self.dy / 2:.3f});")
            corps = f"\\textbf{{{base}}}" if chg else base
            c = f"text={{{couleur}}}," if (chg and couleur) else ""
            t = f"{taille} " if taille else ""
            self.lignes.append(
                f"  \\node[anchor=center,{c}inner sep=0pt] at "
                f"({self.x(i):.3f},{self.y:.3f}) {{{t}{corps}}};")

    def accolade(self, i0, i1, texte, couleur, haut=True, decalage=0.30):
        d = self.dy / 2 + decalage
        y0 = self.y + d if haut else self.y - d
        sens = "" if haut else "mirror,"
        self.lignes.append(
            f"  \\draw[decorate,decoration={{brace,{sens}amplitude=4pt}},{couleur}] "
            f"({self.x(i0) - self.dx / 2:.3f},{y0:.3f}) -- ({self.x(i1) + self.dx / 2:.3f},{y0:.3f}) "
            f"node[midway,{'above' if haut else 'below'}=3pt,font=\\scriptsize,text={{{couleur}}}] "
            f"{{{texte}}};")

    def rendu(self):
        return "\n".join(self.lignes)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", default="gyrA")
    ap.add_argument("--codon", type=int, default=94)
    ap.add_argument("--voisin", type=int, default=95,
                    help="codon adjacent dont le marqueur de sous-lignée fait le fond")
    ap.add_argument("--debut-codon", type=int, default=7581,
                    help="position 1-based de la première base du codon dans le fichier de l'OMS")
    ap.add_argument("--source", default=str(SOURCE))
    ap.add_argument("--fasta", default=str(FASTA))
    ap.add_argument("--residuel", default=str(RESIDUEL))
    ap.add_argument("--figure", default=str(FIGURE))
    ap.add_argument("--tsv", default=str(TSV))
    a = ap.parse_args()

    # --- Lectures. Rien de ce qui suit n'est saisi à la main. -----------------------------------
    aa_ref, par_aa = representations(a.source, a.gene, a.codon)
    aa_ref_voisin, par_aa_voisin = representations(a.source, a.gene, a.voisin)
    debut_locus, fin_locus = a.debut_codon - 1, a.debut_codon + 6   # une base avant, codon+1 après
    locus = sequence_reference(a.fasta, debut_locus, fin_locus)
    codon = locus[a.debut_codon - debut_locus: a.debut_codon - debut_locus + 3]
    codon_voisin = locus[a.debut_codon - debut_locus + 3: a.debut_codon - debut_locus + 6]
    resid = effectifs_residuel(a.residuel)

    # Le marqueur de fond : le variant du codon voisin dont le fichier donne une substitution
    # d'une base, c'est-à-dire celui que l'appeleur assemble avec le codon étudié.
    fond = None
    for aa, reps in sorted(par_aa_voisin.items()):
        simples = [r for r in reps if len(r[1]) == 1 and len(r[2]) == 1]
        if simples:
            fond = (aa, simples[0])
            break
    if fond is None:
        raise SystemExit(f"aucune substitution simple au codon {a.voisin}")
    aa_fond, (pos_fond, ref_fond, alt_fond) = fond

    # --- Familles : un variant est rangé selon l'étendue de son entrée combinée (ou de celle
    #     qu'il n'a pas), c'est-à-dire selon les bases du codon que le variant change. ----------
    familles = defaultdict(list)
    lignes_tsv = []
    for aa, reps in par_aa.items():
        g = groupes(reps)
        combines = [x for x in g if x["longueur"] >= 4]
        # étendue de l'entrée combinée que ce variant AURAIT : du premier changement du codon
        # jusqu'à la base changée par le marqueur de fond.
        simple = min((x for x in g if x["longueur"] == 1), key=lambda x: x["position"], default=None)
        codon3 = min((x for x in g if x["longueur"] == 3), key=lambda x: x["position"], default=None)
        depart = min(x["position"] for x in g if x["longueur"] <= 3) if (simple or codon3) else None
        if combines:
            depart = min(x["position"] for x in combines)
        etendue = (depart, pos_fond - depart + len(ref_fond))
        familles[(etendue[1], bool(combines))].append((aa, g, combines, etendue))

    # Ordre d'exposition : les entrées combinées présentes d'abord (bloc le plus court en haut),
    # les absentes en dernier, pour que l'identité de forme se lise de haut en bas.
    ordre = sorted(familles, key=lambda k: (not k[1], k[0]))

    d = Dessin(debut_locus, fin_locus)
    gauche = -2.95        # abscisse des étiquettes de gauche, en unités de base
    droite = d.n - 0.5    # abscisse des annotations de droite

    # ---------------- Panneau A : le locus et ce que l'appeleur émet ----------------------------
    d.texte(gauche, "\\textbf{A}", ancre="west", taille="\\normalsize")
    d.texte(gauche + 0.55, f"\\textit{{{a.gene}}} reference, coordinates of the WHO file (1-based)",
            ancre="west", taille="\\footnotesize")
    d.descendre(1.5)
    for i in range(d.n):
        d.texte(i, str(debut_locus + i), taille="\\tiny", couleur="refcol", dy=0.0)
    d.descendre(0.9)
    i_codon = a.debut_codon - debut_locus
    d.cellules(locus, 0, "draw=refcol!55,line width=0.3pt", None)
    d.texte(gauche, "reference", ancre="west", taille="\\footnotesize")
    d.descendre(0.70)
    d.accolade(i_codon, i_codon + 2, f"codon {a.codon} ({aa_ref})", "codon94col",
               haut=False, decalage=0.02)
    d.accolade(i_codon + 3, i_codon + 5, f"codon {a.voisin} ({aa_ref_voisin})",
               "codon95col", haut=False, decalage=0.02)
    d.descendre(2.1)

    d.texte(gauche, f"what the variant caller emits for a strain whose codon {a.codon} "
                    f"carries a change at its", ancre="west", taille="\\footnotesize")
    d.descendre(1.0)
    # Un exemple sans entrée combinée (pris au résiduel, le plus porté) et un avec (le bloc de
    # quatre bases du fichier) : les deux appels que la cohorte contient réellement.
    plus_porte = max(resid, key=lambda alt: resid[alt]["souches"])
    bloc_court = min((x for k in familles for (aa, g, c, e) in familles[k] for x in c),
                     key=lambda x: (len(x["ref"]), x["position"]))
    for alt_bloc, pos_bloc, ref_bloc, legende in (
        (plus_porte, debut_locus + 1, locus[1:1 + len(plus_porte)], "first base"),
        (bloc_court["alts"][0], bloc_court["position"], bloc_court["ref"], "second base"),
    ):
        i0, changees = muter(locus, debut_locus, pos_bloc, ref_bloc, alt_bloc)
        d.cellules(alt_bloc, i0, "draw=appelcol,line width=0.5pt", "appelcol", changees)
        d.texte(gauche, legende, ancre="west", taille="\\scriptsize")
        d.texte(droite, f"\\texttt{{{pos_bloc} {ref_bloc}$>$\\,{alt_bloc}}}, one call of "
                        f"{len(alt_bloc)} bases", ancre="west", taille="\\scriptsize",
                couleur="appelcol")
        d.descendre(1.25)
    d.texte(gauche, f"both on a {a.gene} {aa_ref_voisin}{a.voisin}{aa_fond} background "
                    f"(\\texttt{{{pos_fond} {ref_fond}$>$\\,{alt_fond}}}), a sub-lineage marker: "
                    f"the caller joins the two substitutions into one block",
            ancre="west", taille="\\scriptsize")

    # ---------------- Panneau B : ce que le fichier de coordonnées liste ------------------------
    d.descendre(1.7)
    d.texte(gauche, "\\textbf{B}", ancre="west", taille="\\normalsize")
    d.texte(gauche + 0.55, f"what the WHO coordinate file lists for the "
                           f"{len(par_aa)} protein variants of codon {a.codon}",
            ancre="west", taille="\\footnotesize")
    d.descendre(1.6)
    # La bande d'identité ne couvre que les familles dont l'étendue est celle des lignes
    # absentes : c'est là, et là seulement, que l'identité de position et de forme se lit.
    largeur_cible = max(k[0] for k in ordre if not k[1])
    rang_bande, y_haut_b = None, None

    for cle in ordre:
        largeur, a_combine = cle
        if largeur == largeur_cible and rang_bande is None:
            rang_bande, y_haut_b = len(d.lignes), d.y + d.dy
        def porte(m):
            b = bloc_attendu(m[1], m[3], locus, debut_locus, pos_fond, alt_fond)
            return resid[b[2]]["souches"] if (b and b[2] in resid) else 0

        membres = sorted(familles[cle], key=lambda m: (-porte(m), m[0]))
        depart = min(e[0] for _, _, _, e in membres)
        rang = {0: "first", 1: "second", 2: "third"}
        # Quelles bases du codon le variant change : lu sur sa représentation la plus COURTE,
        # c'est-à-dire la substitution minimale, et non sur un codon entier qui emporte une
        # troisième base synonyme.
        bases_changees = set()
        for _, g, _, _ in membres:
            court = min(g, key=lambda x: x["longueur"])
            for k, (br, ba) in enumerate(zip(court["ref"], court["alts"][0])):
                offset = court["position"] + k - a.debut_codon
                if br != ba and 0 <= offset < 3:
                    bases_changees.add(offset)
        ou = " and ".join(rang.get(k, str(k)) for k in sorted(bases_changees))
        titre = (f"combined entry listed: {largeur} bases from {depart} "
                 f"({ou} base of the codon changed)" if a_combine
                 else f"no combined entry listed: it would take {largeur} bases from {depart}, "
                      f"the form the file already uses above ({ou} base changed)")
        d.texte(gauche, titre, ancre="west", taille="\\scriptsize",
                couleur="combinecol" if a_combine else "absentcol")
        d.descendre(1.05)
        for aa, g, combines, etendue in membres:
            haut_y = d.y
            for x in g:
                i0 = x["position"] - debut_locus
                est_combine = x["longueur"] >= 4
                couleur = "combinecol" if est_combine else "listecol"
                # Pour une entrée combinée, montrer d'abord celle qui combine le variant avec le
                # SEUL marqueur de fond : c'est elle que l'appel de l'appeleur peut égaler. Les
                # autres (ici la combinaison avec un second variant du codon voisin) ne diffèrent
                # de la référence qu'en des positions supplémentaires, et sont reléguées.
                if est_combine:
                    permises = set(range(a.debut_codon, a.debut_codon + 3)) | {pos_fond}
                    x["alts"].sort(key=lambda al: sum(
                        1 for k, (br, ba) in enumerate(zip(x["ref"], al))
                        if br != ba and x["position"] + k not in permises))
                alt0 = x["alts"][0]
                _, changees = muter(locus, debut_locus, x["position"], x["ref"], alt0)
                trait = "line width=0.5pt" if est_combine else "line width=0.3pt"
                d.cellules(alt0, i0, f"draw={couleur},{trait}", couleur, changees,
                           taille="\\scriptsize")
                formes = "/".join(x["alts"])
                d.texte(droite, f"\\texttt{{{x['position']} {x['ref']}$>$\\,{formes}}}",
                        ancre="west", taille="\\scriptsize", couleur=couleur)
                lignes_tsv.append({
                    "famille": f"{largeur}_bases_{'listee' if a_combine else 'absente'}",
                    "variant": f"{a.gene}_p.{aa_ref}{a.codon}{aa}", "position": x["position"],
                    "reference": x["ref"], "alternatives": ";".join(x["alts"]),
                    "longueur": x["longueur"], "entree_combinee": int(est_combine),
                    "n_souches_cohorte": "", "statut": "listee",
                })
                d.descendre(1.0)
            if not combines:
                # La ligne qui manque : même position et même longueur que les entrées listées
                # plus haut, reconstruite depuis la substitution simple et le marqueur de fond.
                attendu = bloc_attendu(g, etendue, locus, debut_locus, pos_fond, alt_fond)
                if attendu is None:
                    raise SystemExit(f"{aa} n'a ni entrée combinée ni substitution simple")
                _, ref_bloc, alt_bloc = attendu
                i0 = etendue[0] - debut_locus
                changees = [i0 + k for k in range(len(ref_bloc)) if ref_bloc[k] != alt_bloc[k]]
                d.cellules(alt_bloc, i0, "draw=absentcol,line width=0.6pt,dash pattern=on 2pt off 1.4pt",
                           "absentcol", changees, taille="\\scriptsize")
                e = resid.get(alt_bloc)
                n = e["souches"] if e else 0
                d.texte(droite, f"\\texttt{{{etendue[0]} {ref_bloc}$>$\\,{alt_bloc}}} "
                                f"\\textbf{{not listed}}\\,\\textemdash\\,{n} strains",
                        ancre="west", taille="\\scriptsize", couleur="absentcol")
                lignes_tsv.append({
                    "famille": f"{etendue[1]}_bases_absente",
                    "variant": f"{a.gene}_p.{aa_ref}{a.codon}{aa}", "position": etendue[0],
                    "reference": ref_bloc, "alternatives": alt_bloc, "longueur": etendue[1],
                    "entree_combinee": 0, "n_souches_cohorte": n, "statut": "absente",
                })
                d.descendre(1.0)
            # Nom aligné sur la PREMIÈRE ligne du variant : centré, il ferait croire que la ligne
            # du dessus appartient au variant précédent.
            d.texte(gauche, f"\\textit{{{a.gene}}} {aa_ref}{a.codon}{aa}", ancre="west",
                    taille="\\footnotesize", dy=haut_y - d.y)
        d.descendre(0.45)

    # Trait d'identité : deux verticales marquent l'étendue que les lignes absentes auraient, et
    # que des lignes listées occupent déjà plus haut. Insérées en tête du panneau B pour passer
    # DERRIÈRE les cellules.
    largeurs_listees = sorted({k[0] for k in ordre if k[1]})
    largeur_absente = sorted({k[0] for k in ordre if not k[1]})
    depart_absent = min(e[0] for k in ordre if not k[1] for _, _, _, e in familles[k])
    i_g = depart_absent - debut_locus
    i_d = i_g + max(largeur_absente) - 1
    d.lignes[rang_bande:rang_bande] = [
        f"  \\draw[absentcol!40,line width=0.4pt,dash pattern=on 1.2pt off 1.6pt] "
        f"({d.x(i) + s * d.dx / 2:.3f},{y_haut_b:.3f}) -- "
        f"({d.x(i) + s * d.dx / 2:.3f},{d.y_min_cellule:.3f});"
        for i, s in ((i_g, -1), (i_d, 1))]
    porteurs = sum(resid[alt]["souches"] for alt in resid)
    residuel = ", ".join(
        f"{d_}: {sum(resid[alt]['drogues'][d_][0] for alt in resid)} of "
        f"{sum(resid[alt]['drogues'][d_][0] + resid[alt]['drogues'][d_][1] for alt in resid)} "
        f"labelled carriers resistant"
        for d_ in sorted({k for alt in resid for k in resid[alt]["drogues"]}))

    entete = "\n".join([
        "% Figure 1 d'A1 — GENEREE par analyses/p717_figure_gyra94.py (P7.17), ne pas editer a la main.",
        f"% Sources : {Path(a.source).name}, {Path(a.fasta).name}, {Path(a.residuel).name}.",
        "\\documentclass[border=3pt]{standalone}",
        "\\usepackage[T1]{fontenc}",
        "\\usepackage[utf8]{inputenc}",
        "\\usepackage{tikz}",
        "\\usetikzlibrary{decorations.pathreplacing}",
        *(f"\\definecolor{{{nom}col}}{{HTML}}{{{val}}}" for nom, val in COULEURS.items()),
        "\\begin{document}",
        f"\\begin{{tikzpicture}}[x=1cm,y=1cm,font=\\sffamily]",
    ])
    pied = "\n".join(["\\end{tikzpicture}", "\\end{document}", ""])
    Path(a.figure).parent.mkdir(parents=True, exist_ok=True)
    Path(a.figure).write_text(entete + "\n" + d.rendu() + "\n" + pied, encoding="utf-8")

    champs = ["famille", "variant", "position", "reference", "alternatives", "longueur",
              "entree_combinee", "n_souches_cohorte", "statut"]
    Path(a.tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(a.tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=champs, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(lignes_tsv, key=lambda r: (r["statut"], r["variant"], -r["longueur"])))

    print(f"codon {a.codon} de {a.gene} : {aa_ref} ({codon}) ; codon {a.voisin} : "
          f"{aa_ref_voisin} ({codon_voisin}) ; locus {debut_locus}-{fin_locus} = {locus}")
    print(f"marqueur de fond : {a.gene} {aa_ref_voisin}{a.voisin}{aa_fond} = "
          f"{pos_fond} {ref_fond}>{alt_fond}")
    for cle in ordre:
        print(f"  famille {cle[0]} bases, entrée combinée {'listée' if cle[1] else 'ABSENTE'} : "
              + ", ".join(f"{aa_ref}{a.codon}{aa}" for aa, *_ in sorted(familles[cle])))
    print(f"largeurs listées {largeurs_listees}, largeur des absents {largeur_absente} ; "
          f"{porteurs} souches porteuses des trois absents ; {residuel}")
    print(f"-> {a.figure}\n-> {a.tsv}")


if __name__ == "__main__":
    main()
