import os
import urllib2
import json

def download(data_directory, override=False):
    urls = [
        "http://gmd.mpimp-golm.mpg.de/download/default.aspx?filename=GMD_20111121_MDN35_ALK_MSL.txt",
        "http://gmd.mpimp-golm.mpg.de/download/default.aspx?filename=GMD_20111121_MDN35_FAME_MSL.txt",
        "http://gmd.mpimp-golm.mpg.de/download/default.aspx?filename=GMD_20111121_VAR5_ALK_MSL.txt",
        "http://gmd.mpimp-golm.mpg.de/download/default.aspx?filename=GMD_20111121_VAR5_FAME_MSL.txt"
        ]

    for url in urls:
        file_name = url.split("filename=")[1]
        if override != True or os.path.join(data_directory, file_name) != True:
            file_url = urllib2.urlopen(url)
            with open(os.path.join(data_directory, file_name), "wb") as output:
                output.write(file_url.read())

def extract(metabolite):
    data = {
        "Name" : None,
        "Synonyms" : [],
        "InChI" : None,
        "Accession" : None
    }


    inchikey = None
    for line in metabolite.split("\r\n"):
        if line.startswith("MST N:"):
            data["Name"] = line.split(": ")[1]
        if line.startswith("SYNON: "):
            data["Synonyms"].append(line.split(": ")[1])
        if line.startswith("MET_INCHI: "):
            data["InChI"] = line.split(": ")[1]
        if line.startswith("GMDLINK: "):
            data["Accession"] = line.split(": ")[1].split("http://gmd.mpimp-golm.mpg.de/Analytes/")[1].split(".aspx")[0]
        if line.startswith("MET_INCHIKEY: "):
            inchikey = line.split(": ")[1]

    if inchikey != None and None not in data.values():
        return inchikey, data
    else:
        return None, None

def convert(data_directory, output_fp):

    golm = {}

    for file in os.listdir(data_directory):
        with open(os.path.join(data_directory, file), "rb") as infile:
            metabolite = ""
            for line in infile.readlines():
                metabolite += line
                if line.startswith("\r\n") == True:
                    inchikey, data = extract(metabolite)
                    if [inchikey, data] != [None, None]:
                        golm[inchikey] = data
                        metabolite = ""

    print len(golm.keys())

    with open(output_fp, "wb") as outfile:
        json.dump(golm, outfile, indent=4)

if __name__ == "__main__":
    data_directory = "/home/keo7/Data/dimedb/golm/"

    if os.path.isdir(data_directory) != True:
        os.makedirs(data_directory)
        download(data_directory)
    convert(data_directory, "/home/keo7/Data/dimedb/golm.json")
