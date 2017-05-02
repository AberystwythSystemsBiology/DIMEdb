class Sources(object):
    def __new__(self):
        return {
    "schema" : {
        "kegg_id" : "string",
        "hmdb_id" : "string",
        "chebi_id" : "string",
        "pubchem_id" : "string"
    }
}

class Adduct(object):
    def __new__(self):
        {
            "schema": {
                "type": "string",
                "accurate_mass": "float",
                "isotopic_distribution": "list"
            }
        }

class Metabolites(object):
    def __new__(self):
        return {
            'resource_methods': ["GET"],
            "schema": {
                "name": {
                    "type": "string"
                },
                "synonyms": {
                    "type": "list"
                },
                "molecular_formula": {
                    "type": "string"
                },
                "accurate_mass": {
                    "type": "float"
                },
                "num_atoms": {
                    "type": "integer"
                },
                "inchi": {
                    "type": "string"
                },
                "smiles": {
                    "type": "string"
                },
                "origins": {
                    "type": "list"
                },
                "biofluid_location": {
                    "type": "list"
                },
                "tissue_locations": {
                    "type": "list"
                },
                "pathways": {
                    "type": "list"
                },
                "sources": {
                    "type": "dict",
                    "schema": Sources
                },
                "adducts": {
                    "type": "dict",
                    "schema": {
                        "positive": {
                            "type": "list",
                            "schema": Adduct
                        },
                        "negative": {
                            "type": "list",
                            "schema": Adduct
                        },
                        "neutral": {
                            "type": "list",
                            "schema": Adduct
                        }
                    }
                }
            }
        }