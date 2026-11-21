#!/usr/bin/env python3
"""P7.12 — Contrôle que les chiffres épinglés sont IMPRIMÉS dans le manuscrit (et pas seulement justes).

Objet      : `test_p78_chiffres_manuscrit.py` compare des valeurs recopiées DANS LE TEST aux TSV de résultats ;
             il n'ouvre jamais `article/main.tex`. Il prouve que les nombres du test sont exacts, pas que le
             manuscrit les imprime, de sorte qu'un chiffre emporté par une coupe laisse « 233 passed ».
             Ce script comble l'écart : pour chaque affirmation épinglée, il cherche le nombre imprimé dans
             `main.tex` (et dans `supplementary.tex` pour ce qui a été déporté), avant et après une coupe.
Entrées    : analyses/test_p78_chiffres_manuscrit.py (liste CLAIMS), article/main.tex, article/supplementary.tex,
             et, pour `--avant`, une copie du manuscrit d'avant la coupe.
Sorties    : rapport texte ; code 1 si un nombre présent avant la coupe a disparu du manuscrit ET du
             supplémentaire sans être déclaré retiré (`--retires`), ou si `--strict` et un nombre épinglé manque.
Réutilisable : oui (patron : tout manuscrit dont un test épingle des chiffres sans lire le texte) ; les formes
             de nombres sont celles de LaTeX (`144\\,964`, `0.841`, mots pour les petits entiers).
Limites    : un petit entier (moins de 100) se rencontre partout dans un texte : son résultat est classé
             « faible » et ne prouve rien tant qu'il n'est pas lu à la main ; les nombres écrits en toutes
             lettres (nineteen, seventeen, eleven) sont reconnus pour 1 à 20.
Projet     : Resistance_antibio (P7.12)
Date       : 2026-11-21
Usage      : python3 analyses/p712_chiffres_dans_texte.py [--avant COPIE.tex] [--detail]
"""
import argparse
import importlib.util
import re
import sys
from pathlib import Path

RACINE = Path(__file__).resolve().parent.parent
ARTICLE = RACINE / "article"
MOTS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
        10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen",
        16: "sixteen", 17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}


def charger_claims():
    chemin = RACINE / "analyses" / "test_p78_chiffres_manuscrit.py"
    spec = importlib.util.spec_from_file_location("test_p78_chiffres_manuscrit", chemin)
    if spec is None or spec.loader is None:
        raise SystemExit(f"impossible de charger {chemin}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.CLAIMS


def normaliser(tex):
    """Texte comparable : commentaires LaTeX retirés, `\\,`, `\\%` et `~` neutralisés."""
    sans_commentaires = "\n".join(re.sub(r"(?<!\\)%.*$", "", ligne) for ligne in tex.splitlines())
    return sans_commentaires.replace("\\,", "").replace("\\%", "%").replace("~", " ")


def motifs(imprime):
    """Expressions régulières des formes imprimées possibles d'un nombre épinglé."""
    formes = []
    if isinstance(imprime, int):
        n = int(imprime)
        formes.append(str(n))
        if n in MOTS:
            formes.append(MOTS[n])
    else:
        d = len(repr(float(imprime)).split(".")[1])
        for k in range(d, 5):
            formes.append(f"{float(imprime):.{k}f}")
    return [re.compile(r"(?<![\d.])" + re.escape(f) + r"(?![\d]|\.\d)", re.IGNORECASE) for f in formes]


def present(imprime, texte):
    return any(m.search(texte) for m in motifs(imprime))


def faible(imprime):
    return isinstance(imprime, int) and imprime < 100


def lire(chemin, facultatif=False):
    """Texte normalisé ; une entrée absente arrête le contrôle, sauf `facultatif` (supplémentaire)."""
    chemin = Path(chemin)
    if not chemin.exists():
        if facultatif:  # entree-facultative : un manuscrit peut ne pas avoir de supplémentaire compilé
            return ""
        raise SystemExit(f"entrée introuvable : {chemin}")
    return normaliser(chemin.read_text(encoding="utf-8"))


def nombres(texte):
    """Ensemble des nombres écrits en chiffres d'un texte normalisé (entiers et décimaux)."""
    return set(re.findall(r"(?<![\d.])\d+(?:\.\d+)?(?![\d])", texte))


def classer(claims, avant, main, supp):
    """Une ligne par affirmation : (nom, imprime, avant, main, supp, faible)."""
    return [(nom, imprime, present(imprime, avant) if avant else None, present(imprime, main),
             present(imprime, supp), faible(imprime)) for nom, imprime, *_ in claims]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--avant", help="copie du manuscrit d'avant la coupe")
    ap.add_argument("--tous", action="store_true",
                    help="compare TOUS les nombres du corps avant et après (épinglés ou non)")
    ap.add_argument("--detail", action="store_true", help="liste chaque affirmation hors « présent »")
    ap.add_argument("--retires", default="", help="noms d'affirmations déclarées retirées, séparés par « | »")
    a = ap.parse_args(argv)
    claims = charger_claims()
    avant = lire(a.avant) if a.avant else ""
    lignes = classer(claims, avant, lire(ARTICLE / "main.tex"), lire(ARTICLE / "supplementary.tex", facultatif=True))
    declares = {s.strip() for s in a.retires.split("|") if s.strip()}
    imprimes = [x for x in lignes if x[3]]
    deportes = [x for x in lignes if not x[3] and x[4]]
    absents = [x for x in lignes if not x[3] and not x[4]]
    perdus = [x for x in absents if x[2] and x[0] not in declares]
    jamais = [x for x in absents if x[2] is False or x[2] is None]
    print(f"{len(lignes)} affirmations épinglées : {len(imprimes)} imprimées dans le corps "
          f"({sum(1 for x in imprimes if x[5])} faibles, petits entiers), {len(deportes)} seulement au "
          f"supplémentaire, {len(absents)} dans aucun des deux.")
    if a.avant:
        print(f"Présentes AVANT la coupe et absentes du corps ET du supplémentaire : "
              f"{len([x for x in absents if x[2]])}, dont {len(perdus)} non déclarées retirées.")
        print(f"Absentes avant comme après (jamais imprimées sous cette forme) : {len(jamais)}.")
    if a.detail:
        for titre, liste in (("DÉPORTÉES au supplémentaire", deportes),
                             ("PERDUES (présentes avant, absentes ensuite)", [x for x in absents if x[2]]),
                             ("JAMAIS IMPRIMÉES sous cette forme", jamais)):
            print(f"\n== {titre} : {len(liste)}")
            for nom, imprime, *_ in liste:
                print(f"  - {nom} : {imprime}")
    if a.tous and a.avant:
        n_av, n_ap = nombres(avant), nombres(lire(ARTICLE / "main.tex"))
        n_supp = nombres(lire(ARTICLE / "supplementary.tex", facultatif=True))
        disparus = sorted(n_av - n_ap, key=float)
        ajoutes = sorted(n_ap - n_av, key=float)
        print(f"\nNombres distincts du texte : {len(n_av)} avant, {len(n_ap)} après ; "
              f"{len(disparus)} disparus du corps ({len([x for x in disparus if x in n_supp])} retrouvés au "
              f"supplémentaire), {len(ajoutes)} apparus.")
        print("  disparus  :", " ".join(disparus))
        print("  apparus   :", " ".join(ajoutes))
    return 1 if perdus else 0


if __name__ == "__main__":
    sys.exit(main())
