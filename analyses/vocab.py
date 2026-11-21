"""Vocabulaire contrôlé partagé : normalisation des noms d'antibiotiques.

Toutes les sources (OMS, CRyPTIC, PATRIC, tb-profiler, littérature) utilisent des
graphies différentes (Rifampicin / rifampin / RIF). On les ramène à un nom
canonique + un code à 3 lettres, pour pouvoir joindre les sources entre elles.
"""

# nom canonique -> code 3 lettres
DRUG_CODE = {
    "amikacin": "AMK",
    "bedaquiline": "BDQ",
    "capreomycin": "CAP",
    "ciprofloxacin": "CIP",
    "clofazimine": "CFZ",
    "cycloserine": "CYC",
    "delamanid": "DLM",
    "ethambutol": "EMB",
    "ethionamide": "ETH",
    "isoniazid": "INH",
    "kanamycin": "KAN",
    "levofloxacin": "LEV",
    "linezolid": "LZD",
    "moxifloxacin": "MXF",
    "ofloxacin": "OFX",
    "para-aminosalicylic_acid": "PAS",
    "prothionamide": "PTO",
    "pyrazinamide": "PZA",
    "rifabutin": "RFB",
    "rifampicin": "RIF",
    "streptomycin": "STM",
    # regroupements de classe (tb-profiler)
    "aminoglycosides": "AMG",
    "fluoroquinolones": "FQS",
    # prodrogue PZA listée à part dans PATRIC
    "nicotinamide": "NIC",
}

# toute graphie rencontrée -> nom canonique
_ALIASES = {
    "amk": "amikacin", "ami": "amikacin",
    "bdq": "bedaquiline",
    "cap": "capreomycin",
    "cip": "ciprofloxacin",
    "cfz": "clofazimine",
    "cs": "cycloserine", "cyc": "cycloserine",
    "dlm": "delamanid",
    "emb": "ethambutol",
    "eth": "ethionamide",
    "inh": "isoniazid",
    "kan": "kanamycin", "km": "kanamycin",
    "lev": "levofloxacin", "lfx": "levofloxacin",
    "lzd": "linezolid",
    "mxf": "moxifloxacin",
    "ofx": "ofloxacin", "ofl": "ofloxacin",
    "pas": "para-aminosalicylic_acid",
    "para-aminosalicylic acid": "para-aminosalicylic_acid",
    "para-aminosalicylic_acid": "para-aminosalicylic_acid",
    "para-aminosalisylic_acid": "para-aminosalicylic_acid",  # graphie Coll 2018
    "ethiomide": "ethionamide",                              # graphie Coll 2018
    "pto": "prothionamide", "prothionamide": "prothionamide",
    "pza": "pyrazinamide",
    "rfb": "rifabutin",
    "rif": "rifampicin", "rmp": "rifampicin", "rifampin": "rifampicin",
    "rifampicin": "rifampicin",
    "sm": "streptomycin", "stm": "streptomycin", "str": "streptomycin",
}


def canonical_drug(name):
    """Renvoie (nom_canonique, code) pour une graphie quelconque, ou (name, '?')."""
    if name is None:
        return None, None
    k = str(name).strip().lower()
    canon = _ALIASES.get(k, k)
    return canon, DRUG_CODE.get(canon, "?")
