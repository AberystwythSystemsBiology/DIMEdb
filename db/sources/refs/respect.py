import urllib2
import zipfile
import os
import requests
import re
import json
import multiprocess

def inchikey_to_inchi(inchikey):
    inchi = None
    response = requests.get("http://www.chemspider.com/InChI.asmx/InChIKeyToInChI?inchi_key=%s" % inchikey)
    if response.status_code == 200:
        inchi = re.sub('<[^<]+?>', '', response.text.splitlines()[1])
    return inchi

def process_metabolite(file_p):
    inchikey, inchi = None, None

    data = {
        "Name" : None,
        "InChI" : None,
        "Accession" : None,
        "Synonyms" : []
    }
    f = open(file_p, 'r')
    lines = f.read().split("\r\n")
    for line in lines:
        try:
            key, value = line.split(": ")
        except:
            pass
        if key.startswith("ACCESSION"):
            data["Accession"] = value

        if key.startswith("CH$NAME"):
            if data["Name"] == None:
                data["Name"] = value
            else:
                data["Synonyms"].append(value)

        if key.startswith("CH$LINK"):
            if value.startswith("CAS "):
                cas_value = value.split("CAS ")[1]
                response = requests.get("http://webservice.bridgedb.org/Human/xrefs/Ca/%s" % cas_value)
                if response.status_code == 200:
                    for line in response.text.splitlines():
                        resource = line.split("\t")
                        if resource[1] == "InChIKey":
                            inchikey = resource[0]

    if inchikey != None:
        inchi = inchikey_to_inchi(inchikey)


    if None not in [inchi, inchikey]:
        data["InChI"] = inchi
    return [inchikey, data]

def process(f_directory, output_fp="/tmp/out.json"):

    stripped_respect = {}

    pool = multiprocess.Pool(16)
    processed_data = pool.map_async(process_metabolite, [os.path.join(f_directory, file) for file in os.listdir(f_directory)]).get()
    pool.close()
    pool.join()

    for inchikey, data in processed_data:
        if inchikey != None and None not in data.values():
            stripped_respect[inchikey] = data

    with open(output_fp, "wb") as outfile:
        json.dump(stripped_respect, outfile, indent=4)

def download(fp, override=False):
    respect_url = "http://spectra.psc.riken.jp/menta.cgi/static/respect/respect.zip"

    if override != True or os.path.isfile(fp) != True:
        respect_zip = urllib2.urlopen(respect_url)
        with open(fp, "wb") as output:
            output.write(respect_zip.read())

def unzip(fp, output_dir):
    zf = zipfile.ZipFile(fp)
    zf.extractall(output_dir)
    zf.close()
