import os, xmltodict, collections, json

from joblib import Parallel, delayed

from bioservices import KEGG, KEGGParser

def test(file_name):
    with open("/home/keo7/PycharmProjects/DIMEdb/db/dl-files/hmdb/xml_files/"+file_name, "r") as xml_in:
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
        except TypeError:
            origins = None

        try:
            kegg_dict = KEGGParser().parse(KEGG().get(metabolite["kegg_id"]))
            pathways = kegg_dict["PATHWAY"].keys()
        except TypeError:
            pathways = None
        except KeyError:
            pathways = None
        except Exception, err:
            pathways = None

        sources = {
            "chebi_id": metabolite["chebi_id"],
            "pubchem_id": metabolite["pubchem_compound_id"],
            "kegg_id": metabolite["kegg_id"],
            "hmdb_id": metabolite["accession"]
        }

        return collections.OrderedDict([("_id" , _id), ("name", name), ("synonyms", synonyms),
                                      ("inchi", inchi), ("smiles", smiles), ("biofluid_locations", biofluid_locations),
                                      ("tissue_locations", tissue_locations), ("pathways", pathways), ("sources", sources)])




if __name__ == "__main__":
    files = os.listdir("/home/keo7/PycharmProjects/DIMEdb/db/dl-files/hmdb/xml_files/")


    processed_list = []
    for i in range(0, len(files), 1000):
        print i
        processed_list.extend(Parallel(n_jobs=500)(delayed(test)(file) for file in files[i:i+1000]))
        break

    dict = {}

    for indx, metabolite in enumerate(processed_list):
        dict[indx] = metabolite

    with open(os.path.dirname(__file__)+"/output/hmdb_new.json", "w") as outfile:
        json.dump(dict, outfile, indent=4)
