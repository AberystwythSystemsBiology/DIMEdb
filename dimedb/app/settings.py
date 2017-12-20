import os

MONGO_HOST = os.environ.get('MONGO_HOST', 'localhost')
MONGO_PORT = os.environ.get('MONGO_PORT', 27017)
MONGO_DBNAME = os.environ.get('MONGO_DBNAME', 'dimedb')
MONGO_QUERY_BLACKLIST = []
HATEOAS = False
XML = False
RESOURCE_METHODS = ["GET"]
OPTIMIZE_PAGINATION_FOR_SPEED = True
URL_PREFIX = "api"

mtbl_id_info = {
    "schema": {
        "Name": "string",
        "SMILES": "string",
        "Synonyms": "list",
        "Systematic Name": "string",
        "IUPAC Name": "string",
        "InChI": "string",
        "Molecular Formula": "string"
    }
}

mtbl_phprop = {
    "schema": {
        "Secondary Amines": "float",
        "Ether Oxygens": "float",
        "Heavy Atoms": "float",
        "Rings": "float",
        "Hydrogen Bond Acceptors": "float",
        "Aromatic Rings": "float",
        "Fraction of SP3 Carbon": "float",
        "Carboxylic Acids": "float",
        "Polar Surface Area": "float",
        "Rotatable Bonds": "float",
        "clogP": "float",
        "Hydroxy Groups": "float",
        "Formal Charge": "float",
        "Hydrogen Bond Donors": "float"
    }
}

mtbl_adducts = {
    "type": "list",
    "schema": {
        "Polarity": "string",
        "Adduct": "string",
        "Accurate Mass": "float",
        "Isotopic Distribution": "list"
    }
}


mtbl_extsources = {
    "schema": {
        "MassBank": "string",
        "KEGG Drug": "string",
        "BioCyc": "string",
        "CAS": "string",
        "GOLM": "string",
        "Spektraris": "string",
        "PubChem": "string",
        "KEGG Compound": "string",
        "HMDB Accession": "string",
        "Respect": "string",
        "Wikidata": "string",
        "Chemspider": "string",
        "NMR Shift DB": "string",
        "ChEBI": "string"
    }
}

mtbl_pthwys = {
    "schema": {
        "KEGG": "list",
        "SMPDB": "list",
        "BioCyc": "list"
    }
}

mtbl_taxprop = {
    "schema": {
        "HMDB": {
            "schema": {
                "Origins": "list",
                "Tissue Locations": "list",
                "Biofluid Locations": "list"
            }
        }
    }
}


metabolites = {
    "resource_methods": ["GET"],
    "schema": {
        "Identification Information": mtbl_id_info,
        "Physicochemical Properties": mtbl_phprop,
        "External Sources": mtbl_extsources,
        "Adducts": mtbl_adducts,
        "Pathways": mtbl_pthwys,
        "Taxonomic Properties": mtbl_taxprop
    }
}

pathways = {
    "resource_methods": ["GET"]
}

DOMAIN = {
    "metabolites": metabolites,
    "pathways": pathways
}
