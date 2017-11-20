import xmltodict
import json
import os
import urllib2
import zipfile

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


def extract_from_xml(hmdb_fp, output_fp):
    metabolite_text = ""
    stripped_metabolites = {}
    keep_going = False

    for line in open(hmdb_fp, "r"):
        if line.strip() == "<metabolite>":
            keep_going = True
        if keep_going == True:
            metabolite_text += line
        if line.strip() == "</metabolite>":
            keep_going = False
            inchi, data = get_data(metabolite_text)
            stripped_metabolites[inchi] = data
            metabolite_text = ""

    with open(output_fp, "wb") as outfile:
        json.dump(stripped_metabolites, outfile, indent=4)

def download(data_directory, hmdb_fp, override=False):
    hmdb_url = "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"

    if os.path.isdir(data_directory) != True:
        os.makedirs(data_directory)

    if override != True or os.path.isfile(hmdb_fp) != True:
        hmdb = urllib2.urlopen(hmdb_url)
        with open(hmdb_fp, "wb") as output:
            output.write(hmdb.read())

def unzip(data_directory, hmdb_fp):
    zf = zipfile.ZipFile(hmdb_fp)
    zf.extractall(data_directory)
    zf.close()
