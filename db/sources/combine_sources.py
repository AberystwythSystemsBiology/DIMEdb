import os, json

def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)

def combine_all(pubchem, hmdb, chebi):
    combined = {}

    for inchikey in chebi.keys():
        combined_compound = {
            "PubChem ID" : None,
            "HMDB Accession" : None,
            "ChEBI ID" : None,
            "InChI" : None
        }
        try:
            combined_compound["PubChem ID"] = pubchem[inchikey]["PubChem ID"]
            del pubchem[inchikey]
        except KeyError:
            pass

        try:
            combined_compound["HMDB Accession"] = hmdb[inchikey]["Accession"]
            del hmdb[inchikey]
        except KeyError:
            pass

        combined_compound["ChEBI ID"] = chebi[inchikey]["ChEBI ID"]
        combined_compound["InChI"] = chebi[inchikey]["InChI"]

        combined[inchikey] = combined_compound
        del chebi[inchikey]

    for inchikey in pubchem.keys():
        combined_compound = {
            "PubChem ID": None,
            "HMDB Accession": None,
            "ChEBI ID": None,
            "InChI": None
        }
        try:
            combined_compound["HMDB Accession"] = hmdb[inchikey]["Accession"]
            del hmdb[inchikey]
        except KeyError:
            pass

        try:
            combined_compound["ChEBI ID"] = chebi[inchikey]["ChEBI ID"]
            del chebi[inchikey]
        except KeyError:
            pass

        combined_compound["PubChem ID"] = pubchem[inchikey]["PubChem ID"]
        combined_compound["InChI"] = pubchem[inchikey]["InChI"]

        combined[inchikey] = combined_compound
        del pubchem[inchikey]

    for inchikey in hmdb.keys():
        combined_compound = {
            "PubChem ID": None,
            "HMDB Accession": None,
            "ChEBI ID": None,
            "InChI": None
        }

        try:
            combined_compound["ChEBI ID"] = chebi[inchikey]["ChEBI ID"]
            del chebi[inchikey]
        except KeyError:
            pass

        try:
            combined_compound["PubChem ID"] = pubchem[inchikey]["PubChem ID"]
            del pubchem[inchikey]
        except KeyError:
            pass

        combined_compound["HMDB Accession"] = hmdb[inchikey]["Accession"]
        combined_compound["InChI"] = hmdb[inchikey]["InChI"]

        combined[inchikey] = combined_compound
        del hmdb[inchikey]

    return combined

if __name__ == "__main__":
    directory = "/home/keo7/.data/dimedb/"

    chebi = load_json(directory+"stripped_chebi.json")
    hmdb = load_json(directory+"stripped_hmdb.json")
    pubchem = load_json(directory+"stripped_pubchem.json")

    combined = combine_all(pubchem, hmdb, chebi)

    print len(combined.keys())

    with open(directory+"/combined_data.json", "wb") as out_file:
        json.dump(combined, out_file, indent=4)