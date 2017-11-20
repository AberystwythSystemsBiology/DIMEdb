import pybel

import os
import json

def download(data_directory):
    os.system("cd "+data_directory+"; svn co http://www.massbank.jp/SVN/OpenData/record/")

def inchi_to_inchikey(inchiString):
    try:
        mol = pybel.readstring("inchi", inchiString)
        return mol.write("inchikey").strip()
    except IOError:
        return None

def convert_data(data_directory, output_fp):

    massbank_data = {}

    # Need to convert InChi to InChi key, and move accession into values.

    for subdir, dir, files in os.walk(data_directory):
        if ".svn" not in subdir:
            for file in files:
                data = {
                    "Name" : None,
                    "Synonyms" : [],
                    "Accession" : None,
                    "InChI" : None
                }
                f = open(os.path.join(subdir, file), 'r')
                lines = f.read().split("\r\n")
                for line in lines:
                    try:
                        record, value = line.split(": ")
                        value = value.encode("utf-8").strip()
                        if record.startswith("ACCESSION"):
                            data["Accession"] = value
                        if record.startswith("CH$NAME"):
                            if data["Name"] == None:
                                data["Name"] = value.title()
                            else:
                                data["Synonyms"].append(value.title())
                        if record.startswith("CH$IUPAC"):
                            data["InChI"] = value.replace("IInChI=", "")
                    except ValueError:
                        pass
                if data["InChI"] != None and data["InChI"] != "N/A":
                    inchikey = inchi_to_inchikey(data["InChI"])
                    if inchikey != None:
                        massbank_data[inchikey] = data

    print len(massbank_data.keys())

    with open(output_fp, "wb") as outfile:
        json.dump(massbank_data, outfile, indent=4)
