import os

MONGO_HOST = os.environ.get('MONGO_HOST', 'localhost')
MONGO_PORT = os.environ.get('MONGO_PORT', 27017)
MONGO_DBNAME = os.environ.get('MONGO_DBNAME', 'dimedb')
MONGO_QUERY_BLACKLIST = ["$where"]
HATEOAS = False
query_objectid_as_string = False
PAGINATION_DEFAULT = 500
XML = False

URL_PREFIX = "api"

identification_information = {
    "schema" : {
        "Name" : "string",
        "SMILES" : "string",
        "Synonyms" : "list",
        "Systematic Name" : "string",
        "IUPAC Name" : "string",
        "InChI" : "string",
        "Molecular Formula" : "string"
    }
}

physiochemical_properties = {
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
        "clogP" : "float",
        "Hydroxy Groups": "float",
        "Formal Charge": "float",
        "Hydrogen Bond Donors" : "float"
    }
}

adducts = {
    "type" : "list",
    "schema" : {
        "Polarity" : "string",
        "Adduct" : "string",
        "Accurate Mass" : "float",
        "Isotopic Distribution" : "list"
    }
}


external_sources = {
    "schema" : {
        "ChEBI ID" : "string",
        "PubChem ID" :"string",
        "HMDB Accession" : "string",
        "CAS" :  "string",
        "Wikidata" : "string",
        "Chemspider" : "string",
        "BioCyc" : "string"
    }
}

pathways = {
    "schema" : {
        "KEGG" : "list",
        "SMPDB" : "list",
        "BioCyc" : "list"
    }
}

taxonomic_properties = {
    "schema" : {
        "HMDB" : {
            "schema" : {
                "Origins" : "list",
                "Tissue Locations" : "list",
                "Biofluid Locations" : "list"
            }
        }
    }
}



metabolites = {
    "resource_methods": ["GET"],
    "schema" : {
        "Identification Information" : identification_information,
        "Physicochemical Properties" : physiochemical_properties,
        "External Sources" : external_sources,
        "Adducts" : adducts,
        "Pathways" : pathways,
        "Taxonomic Properties" : taxonomic_properties
    }
}


DOMAIN = {
    "metabolites" : metabolites
}