import os

MONGO_HOST = os.environ.get('MONGO_HOST', 'localhost')
MONGO_PORT = os.environ.get('MONGO_PORT', 27017)
MONGO_DBNAME = os.environ.get('MONGO_DBNAME', 'dimedb')
MONGO_QUERY_BLACKLIST = ["$where"]
HATEOAS = False
PAGINATION_DEFAULT = 500
XML = False

URL_PREFIX = "api"

positive_adduct = {
    "schema" : {
        "type" : "string",
        "accurate_mass" : "float",
        "isotopic_distribution" : "list"
    }
}


negative_adduct = {
    "schema" : {
        "type" : "string",
        "accurate_mass" : "float",
        "isotopic_distribution" : "list"
    }
}

neutral_adduct = {
    "schema" : {
        "type" : "string",
        "accurate_mass" : "float",
        "isotopic_distribution" : "list"
    }
}

metabolites = {
    'resource_methods': ["GET"],
    "schema" : {
        "Name" : {
            "type" : "string"
        },
        "Synonyms" : {
            "type" : "list"
        },
        "IUPAC Name" : {
            "type" : "string"
        },
        "Molecular Formula" : {
            "type" : "string"
        },
        "Adducts" : {
            "type" : "dict",
            "schema" : {
                "Negative" : {
                    "type" : "list",
                    "schema" : negative_adduct
                },
                "Neutral" : {
                    "type" : "list",
                    "schema" : neutral_adduct
                },
                "Positive" : {
                    "type" : "list",
                    "schema" : positive_adduct
                }
            }
        },
        "Identifiers" : {
            "type" : "list",
            "schema" : {
                "SMILES" : "string",
                "CAS" : "string",
                "PubChem-compound" : "string",
                "KEGG Compound" : "string",
                "HMDB" : "string",
                "InChI" : "string",
                "Wikidata" : "string",
                "Chemspider" : "string",
                "ChEBI": "string"
            }
        },
        "Properties" : {
            "type" : "list",
            "schema" : {
                "Secondary Amines" : "float",
                "Ether Oxygens" : "float",
                "Heavy Atoms" : "float",
                "Rings" : "float",
                "Hydrogen Bond Acceptors" : "float",
                "Aromatic Rings": "float",
                "Fraction of SP3 Carbon": "float",
                "Carboxylic Acids" : "float",
                "Polar Surface Area" : "float",
                "Rotatable Bonds": "float",
                "logP" : "float",
                "Hydroxy Groups": "float",
                "Formal Charge": "float",
                "Hydrogen Bond Donors" : "float"
            }
        }
    }
}


DOMAIN = {
    "metabolites" : metabolites
}