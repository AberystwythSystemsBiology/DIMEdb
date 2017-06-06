import os

MONGO_HOST = os.environ.get('MONGO_HOST', 'localhost')
MONGO_PORT = os.environ.get('MONGO_PORT', 27017)
MONGO_DBNAME = os.environ.get('MONGO_DBNAME', 'dimedb')
MONGO_QUERY_BLACKLIST = ["$where"]
QUERY_MAX_RESULTS = 1000
PAGINATION_LIMIT = 500
PAGINATION_DEFAULT = 500
XML = False

URL_PREFIX = "api"

adduct = {
"schema" : {
        "type" : "string",
        "accurate_mass" : "float",
        "isotopic_distribution" : "list"
    }
}



sources = {
    "schema" : {
        "kegg_id" : "string",
        "hmdb_id" : "string",
        "chebi_id" : "string",
        "pubchem_id" : "string"
    }
}

metabolites = {
    'resource_methods': ["GET"],
    "schema" : {
        "name" : {
            "type" : "string"
        },
        "other_names" : {
            "type" : "list"
        },
        "chemical_formula" : {
            "type" : "string"
        },
        "accurate_mass" : {
            "type" : "float"
        },
        "num_atoms" : {
            "type" : "integer"
        },
        "inchi" : {
            "type" : "string"
         },
        "smiles" : {
            "type" : "string"
        },
        "origins" : {
            "type" : "list"
        },
        "biofluid_locations" : {
            "type" : "list"
        },
        "tissue_locations" : {
            "type" : "list"
        },
        "pathways" : {
            "type" : "list"
        },
        "sources" : {
          "type" : "dict",
            "schema" : sources
        },
        "adducts" : {
            "type" : "dict",
            "schema" : {
                "positive" : {
                 "type" : "list",
                 "schema" : adduct
                },
                "negative" : {
                 "type" : "list",
                 "schema" : adduct
                },
                "neutral" : {
                 "type" : "list",
                 "schema" : adduct
                }
            }
        }
    }
}


DOMAIN = {
    "metabolites" : metabolites
}