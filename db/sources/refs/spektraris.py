import urllib2
import zipfile
import os
import csv
import json

def download(data_directory, spektraris_path, override=False):
    spektraris_url = "http://langelabtools.wsu.edu/nmr/static/downloads/database_flatfiles.zip"

    if os.path.isdir(data_directory) != True:
        os.makedirs(data_directory)

    if override != True or os.path.isfile(spektraris_path) != True:
        spek = urllib2.urlopen(spektraris_url)
        with open(spektraris_path, "wb") as output:
            output.write(spek.read())

def unzip(data_directory, hmdb_fp):
    zf = zipfile.ZipFile(hmdb_fp)
    zf.extractall(data_directory)
    zf.close()

def convert(record_fp, output_fp):
    spektraris_data= {}

    with open(record_fp, "rb") as infile:
        reader = csv.reader(infile, delimiter="\t")
        reader.next()

        for row in reader:
            data = {
                "Name" : None,
                "Synonyms" : [],
                "Accession" : None,
                "InChI" : None
            }
            data["Accession"] = row[0]
            data["InChI"] = row[4]
            if row[5] != "":
                data["Synonyms"].append(unicode(row[5], errors='ignore'))
            data["Name"] = unicode(row[4], errors='ignore')

            spektraris_data[row[2].replace("InChIKey=", "")] = data

    print spektraris_data

    with open(output_fp, "wb") as outfile:
        json.dump(spektraris_data, outfile, indent=4)
