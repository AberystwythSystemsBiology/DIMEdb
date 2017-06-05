import os, xmltodict, collections, json
from xml.etree import ElementTree as ET


from joblib import Parallel, delayed

def split_hmdb_file(file_path=None, out_dir=None):
    data = ET.iterparse(file_path, events=("end", ))
    count = 0
    for event, element in data:
        count += 1
        if element.tag == "{http://www.hmdb.ca}metabolite":
            accession = element.find("{http://www.hmdb.ca}accession").text
            with open(out_dir+accession+".xml", "wb") as save_file:
                save_file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
                save_file.write(ET.tostring(element).replace("ns0:", "")) #Forgive me father for I have sinned.



def parse_hmdb(file_name):
    origins_dictonary = {
        "Toxin/Pollutant": "Toxin",
        "Plant": "Plant",
        "Drug metabolite": "Drug",
        "Drug or steroid metabolite": "Drug",
        "Drug": "Drug",
        "Food": "Food",
        "Microbial": "Microbial",
        "Micorbial" : "Microbial",
        "Cosmetic": "Cosmetic",
        "Endogenous": "Endogenous",
        "Exogenous" : "Exogenous"
    }

    with open(file_name, "r") as xml_in:
        metabolite = xmltodict.parse(xml_in.read())["metabolite"]

        try:
            _id = metabolite["inchikey"].replace("InChIKey=", "").replace("-", "_")
        except AttributeError:
            # Not interested in anything that doesn't have a IchiKey.
            return

        name = metabolite["name"]
        try:
            synonyms = metabolite["synonyms"]["synonym"]
            if type(synonyms) != list:
                synonyms = [synonyms]
        except TypeError:
            synonyms = None

        inchi = metabolite["inchi"]
        smiles = metabolite["smiles"]

        try:
            biofluid_locations = metabolite["biofluid_locations"]["biofluid"]
            if type(biofluid_locations) != list:
                biofluid_locations = [biofluid_locations]
        except TypeError:
            biofluid_locations = None

        try:
            tissue_locations = metabolite["tissue_locations"]["tissue"]
            if type(tissue_locations) != list:
                tissue_locations = [tissue_locations]
        except TypeError:
            tissue_locations = None

        try:
            origins = metabolite["ontology"]["origins"]["origin"]
            if type(origins) != list:
                origins = [origins]
            origins = [origins_dictonary[x] for x in origins]
        except TypeError:
            origins = None

        sources = {
            "chebi_id": metabolite["chebi_id"],
            "pubchem_id": metabolite["pubchem_compound_id"],
            "kegg_id": metabolite["kegg_id"],
            "hmdb_id": metabolite["accession"]
        }

        return collections.OrderedDict([("_id" , _id), ("name", name), ("synonyms", synonyms), ("origins", origins),
                                      ("inchi", inchi), ("smiles", smiles), ("biofluid_locations", biofluid_locations),
                                      ("tissue_locations", tissue_locations), ("sources", sources)])


if __name__ == "__main__":
    hmdb_file = "/home/keo7/.data/hmdb/hmdb_metabolites.xml"
    hmdb_directory = "/home/keo7/.data/hmdb/out/"
    output_directory = "./output/"
    files = os.listdir(hmdb_directory)

    if files == []:
        split_hmdb_file(file_path=hmdb_file, out_dir=hmdb_directory)
        files = os.listdir(hmdb_directory)


    processed_list = []

    fr = range(0, len(files), 500)
    for idx, i in enumerate(fr):
        print idx+1,  "/", len(fr)
        processed_data = Parallel(n_jobs=100)(delayed(parse_hmdb)(hmdb_directory+file) for file in files[i:i+500])
        processed_list.extend([x for x in processed_data if x != None])


    dict = {}

    for indx, metabolite in enumerate(processed_list):
        dict[indx] = metabolite

    with open(os.path.dirname(__file__)+"/output/hmdb.json", "wb") as outfile:
        json.dump(dict, outfile, indent=4)

