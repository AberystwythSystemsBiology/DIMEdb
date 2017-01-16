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
            try:
                d = xmltodict.parse(hmdb_xml_f.read())["metabolite"]
                try:
                    origins = d["ontology"]["origins"]["origin"]
                    if type(origins) != list:
                        origins = [origins]
                except TypeError:
                    origins = []
                entity = {
                    "name" : d["name"],
                    "smiles" : d["smiles"],
                    "synonyms": d["synonyms"],
                    "origins" : origins,
                }
                hmdb_dict[d["accession"]] = entity
            except Exception, err:
                continue
    return hmdb_dict

def save_hmdb_xml(d, fp="./output/hmdb.json"):
    with open(fp, "w") as hmdb_json:
        json.dump(d, hmdb_json, indent=4)

if __name__ == "__main__":
    download()
    d = parse_hmdb_xml()
    save_hmdb_xml(d)