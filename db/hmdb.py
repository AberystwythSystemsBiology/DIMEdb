import xmltodict, os, json, collections, urllib2, zipfile
from tqdm import tqdm
from joblib import Parallel

'''
    NOTE TO SELF: JOBLIB MUST BE INSTALLED USING

    sudo apt-get install python-joblib
'''


def download(fp="./dl-files/hmdb/"):
    if not os.path.exists(fp):
        os.makedirs(fp)
    hmdb_zipfile = fp+"hmdb_metabolites.zip"
    if os.path.isfile(hmdb_zipfile) == True:
        print "File already found, skipping download"
    else:
        download = urllib2.urlopen("http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip")
        save_file = open(hmdb_zipfile, "wb")
        save_file.write(download.read())
        save_file.close()
    if not os.path.exists(fp+"xml_files/"):
        os.makedirs(fp+"xml_files/")
        with zipfile.ZipFile(hmdb_zipfile) as zf:
            zip_list = [x for x in zf.namelist() if x != "hmdb_metabolites.xml"]
            zf.extractall(fp+"xml_files/", members=zip_list)
    else:
        print "You seem to already have it unzipped?"

def parse_hmdb_xml(fd="./dl-files/hmdb/xml_files/"):
    hmdb_dict = collections.defaultdict(dict)
    for file in tqdm(os.listdir(fd)):
        with open(fd+file, "r") as hmdb_xml_f:
            if file.endswith(".xml"):
                try:
                    d = xmltodict.parse(hmdb_xml_f.read())["metabolite"]
                    try:
                        origins = d["ontology"]["origins"]["origin"]
                        if type(origins) != list:
                            origins = [origins]
                    except TypeError:
                        origins = []

                    pathway = []

                    biofluids = []
                    try:
                        for i in d["biofluid_locations"]["biofluid"]:
                            if i == "Cerebrospinal Fluid (CSF)":
                                i = "CSF"
                            biofluids.append(i)
                    except TypeError:
                        pass
                    try:
                        for i in d["pathways"]["pathway"]:
                            try:
                                name = i["name"]
                                try:
                                    kegg_id = i["kegg_map_id"]
                                except TypeError:
                                    kegg_id = None
                                try:
                                    smpdb_id = i["smpdb_id"]
                                except TypeError:
                                    smpdb_id = None
                                pathway.append({
                                    "name": name,
                                    "kegg_id": kegg_id,
                                    "smpdb_id": smpdb_id
                                })
                            except TypeError:
                                pass
                    except TypeError:
                        pass

                    try:
                        synonyms = d["synonyms"]["synonym"]
                        if type(synonyms) != list:
                            synonyms = [synonyms]
                    except TypeError:
                        synonyms = []

                    entity = {
                        "name" : d["name"],
                        "smiles" : d["smiles"],
                        "synonyms": synonyms,
                        "origins" : origins,
                        "pathways" : pathway,
                        "biofluids" : biofluids,
                        "sources" : {
                            "chebi_id" : d["chebi_id"],
                            "pubchem_id" : d["pubchem_compound_id"],
                            "kegg_id" : d["kegg_id"],
                            "hmdb_id" : d["accession"]
                        }

                    }
                    hmdb_dict[d["accession"]] = entity
                except ValueError, err:
                    continue
    return hmdb_dict

def save_hmdb_xml(d, fp="./output/hmdb.json"):
    with open(fp, "w") as hmdb_json:
        json.dump(d, hmdb_json, indent=4)

if __name__ == "__main__":
    download()
    d = parse_hmdb_xml()
    save_hmdb_xml(d)