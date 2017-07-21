import xmltodict, json

def get_data(compound):
    cpd = xmltodict.parse(compound)["metabolite"]

    name = cpd["name"]

    try:
        synonyms = cpd["synonyms"]["synonym"]
    except TypeError:
        synonyms = None

    if synonyms == None:
        synonyms = []
    elif type(synonyms) == str or type(synonyms) == unicode:
        synonyms = [synonyms]

    id = cpd["accession"]
    inchi = cpd["inchi"]
    inchikey = cpd["inchikey"]

    try:
        biofluids = cpd["biofluid_locations"]["biofluid"]
    except TypeError:
        biofluids = []

    if type(biofluids) == str or type(biofluids) == unicode:
        biofluids = [biofluids]
    elif type(biofluids) == list:
        biofluids = biofluids

    try:
        tissues= cpd["tissue_locations"]["tissue"]
    except TypeError:
        tissues = []

    pathways = []
    try:
        for p in cpd["pathways"]["pathway"]:
            pathways.append(p["smpdb_id"])
    except TypeError:
        pass

    if type(tissues) == str or type(tissues) == unicode:
        tissues = [tissues]
    elif type(tissues) == list:
        tissues = tissues
    try:
        origins = cpd["ontology"]["origins"]["origin"]
    except TypeError:
        origins = []

    if type(origins) == str or type(origins) == unicode:
            origins = [origins]
    elif type(origins) == list:
            origins = origins

    return inchikey, {
        "Accession" : id,
        "InChI" : inchi,
        "Name" : name,
        "Synonyms" : synonyms,
        "Sources" : {
            "Biofluid Locations" : biofluids,
            "Tissue Locations" : tissues,
            "Origins" : origins
        },
        "SMPDB Pathways" : pathways
    }


if __name__ == "__main__":
    directory = "/home/keo7/.data/dimedb/"
    hmdb_xml = open(directory+"hmdb/hmdb_metabolites.xml", "r")

    metabolite_text = ""

    stripped_metabolites = {}

    keep_going = False

    for line in hmdb_xml:
        if line.strip() == "<metabolite>":
            keep_going = True
        if keep_going == True:
            metabolite_text += line
        if line.strip() == "</metabolite>":
            keep_going = False
            inchi, data = get_data(metabolite_text)
            stripped_metabolites[inchi] = data
            metabolite_text = ""

    with open(directory+"/stripped_hmdb.json", "wb") as out_file:
        json.dump(stripped_metabolites, out_file, indent=4)


