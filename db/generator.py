import json, pubchempy

directory = "/home/keo7/.data/dimedb/"


def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)

combined = load_json(directory + "/combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")

def naming(inchikey):
    naming_info = {
        "Name" : None,
        "IUPAC Name" : None,
        "Systematic Name" : None,
        "Synonyms" : []
    }


    if combined[inchikey]["HMDB Accession"] != None:
        naming_info["Name"] = hmdb[inchikey]["Name"]
        naming_info["Synonyms"] = hmdb[inchikey]["Synonyms"]

    if combined[inchikey]["ChEBI ID"] != None:
        if naming_info["Name"] == None:
            naming_info["Name"] = chebi[inchikey]["Name"]
        if naming_info["IUPAC Name"] == None:
            naming_info["IUPAC Name"] = chebi[inchikey]["IUPAC Name"]
        if naming_info["Synonyms"] == []:
            naming_info["Synonyms"] = chebi[inchikey]["Synonyms"]

    if combined[inchikey]["PubChem ID"] != None:
        compound = pubchempy.get_compounds(combined[inchikey]["PubChem ID"])[0]
        if naming_info["Name"] == None:
            naming_info["Name"] = compound.synonyms[0]
        if naming_info["Synonyms"] == []:
            naming_info["Synonyms"] = compound.synonyms[1:]
        if naming_info["IUPAC Name"] == None:
            naming_info["IUPAC Name"] = compound.iupac_name

    return naming_info

if __name__ == "__main__":
    for inchikey in combined.keys():
        naming_info = naming(inchikey)

        break