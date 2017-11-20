import urllib2
import os
import pybel
import json
import pubchempy


def download(nmrshiftdb_fp, data_directory, override=False):
    nmrshiftdb_url = "https://sourceforge.net/p/nmrshiftdb2/code/HEAD/tree/trunk/snapshots/nmrshiftdb2.sd?format=raw"

    if os.path.isdir(data_directory) != True:
        os.makedirs(data_directory)

    if override != True or os.path.isfile(nmrshiftdb_fp) != True:
        nmrshiftdb = urllib2.urlopen(nmrshiftdb_url)
        with open(nmrshiftdb_fp, "wb") as output:
            output.write(nmrshiftdb.read())


def cdf(metabolite):
    try:
        mol = pybel.readstring("sdf", metabolite)
        return mol.write("inchikey").strip(), mol.write("inchi").strip()
    except IOError:
        return None, None

def convert(fp, op):
    nmrshift_data = {}

    with open(fp, "rb") as infile:
        lines = infile.readlines()
        # Removal of copyright.
        lines[0] = lines[0].split(" ")[0]+"\n"
        metabolite = ""
        for line in lines:
            metabolite += line
            if "$$$$" in line:
                vals = metabolite.split("\n")
                name = vals[0].title()
                accession = vals[2].split(" ")[1]
                inchikey, inchi = cdf(metabolite)
                if name == "":
                    try:
                        comps = pubchempy.get_compounds(inchikey, "inchikey")
                        if len(comps) > 0:
                            name = comps[0].iupac_name
                    except Exception:
                        pass
                if inchikey != None and name != "":
                    nmrshift_data[inchikey] = {
                        "Accession" : accession,
                        "Name" : name,
                        "InChI" : inchi
                    }

                metabolite = ""

    with open(op, "wb") as outfile:
        json.dump(nmrshift_data, outfile, indent=4)
