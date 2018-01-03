import xmltodict
import json
import os
import urllib2
import zipfile
import multiprocess
import globals


class hmdb(object):
    pass

def strip_xml(fp):
    text = ""
    metabolites = []

    keep_going = False

    for line in open(fp, "r"):
        if line.strip() == "<metabolite>":
            keep_going = True
        if keep_going == True:
            text += line
        if line.strip() == "</metabolite>":
            keep_going = False
            metabolites.append(text)
            text = ""
    return metabolites



def download(data_dir, fp, override=False):
    url = "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
    if os.path.isdir(data_dir) != True:
        os.makedirs(data_dir)
    if override != True or os.path.isfile(fp) != True:
        hmdb = urllib2.urlopen(url)
        with open(fp, "wb") as output:
            output.write(hmdb.read())

    zf = zipfile.ZipFile(fp)
    zf.extractall(data_dir)
    zf.close()

def process_metabolite(metabolite):
    metabolite = xmltodict.parse(metabolite)["metabolite"]

    try:
        synonyms = metabolite["synonyms"]["synonym"]
        if type(synonyms) == str or type(synonyms) == unicode:
            synonyms = [synonyms]
    except TypeError:
        synonyms = []

    pathways = []
    try:
        for pathway in metabolite["pathways"]["pathway"]:
            pathways.append(pathway["smpdb_id"])
    except TypeError:
        pass

    biofluids = []
    try:
        biofluids = metabolite["biofluid_locations"]["biofluid"]
        if type(biofluids) == str or type(biofluids) == unicode:
            biofluids = [biofluids]
    except TypeError:
        biofluids = []

    tissues = []
    try:
        tissues = metabolite["tissue_locations"]["tissue"]
        if type(tissues) == str or type(tissues) == unicode:
            tissues = [tissues]
    except TypeError:
        tissues = []

    origins = []
    try:
        origins = metabolite["ontology"]["origins"]["origin"]
        if type(origins) == str or type(origins) == unicode:
            origins = [origins]
    except TypeError:
        origins = []

    return [metabolite["inchikey"], {
        "Accession" : metabolite["accession"],
        "InChI" :  metabolite["inchi"],
        "Name" : metabolite["name"],
        "Synonyms" : synonyms,
        "SMPDB Pathways" : pathways,
        "Sources" : {
            "Biofluid Locations" : biofluids,
            "Tissue Locations" : tissues,
            "Origins" : origins
        }
    }]


if __name__ == "__main__":
    data_dir = os.path.join(globals.data_dir, "hmdb/")
    fp = os.path.join(data_dir, "hmdb_metabolites.zip")
    o_fp = os.path.join(globals.data_dir, "hmdb.json")
    #download(data_dir, fp)
    m = strip_xml(fp.replace("zip", "xml"))

    pool = multiprocess.Pool(8)
    s_m = pool.map_async(process_metabolite, m).get()
    pool.close()
    pool.join()

    f_sm = {}

    for inchikey, data in s_m:
        f_sm[inchikey] = data

    with open(o_fp, "wb") as outfile:
        json.dump(f_sm, outfile, indent=4)
