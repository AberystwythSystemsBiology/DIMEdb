import json
import os
import operator

def generate_combined_dictonary(files, o_fp):

    combined_dict = {}

    for filename, source in files:
        print filename
        for inchikey, data in source.items():
            if inchikey not in combined_dict.keys():
                combined_data = {
                    "Name" : [],
                    "Synonyms" : [],
                    "InChI" : [],
                    "Source Information" : {
                        "GOLM ID" : None,
                        "HMDB ID" : None,
                        "MassBank ID" : None,
                        "NMR ShiftDB ID" : None,
                        "Respect ID" : None,
                        "Spektraris ID" : None
                        },
                    "HMDB Information" : {
                        "SMPDB Pathways": [],
                        "Sources": {

                        }
                    }
                }
            else:
                combined_data = combined_dict[inchikey]

            combined_data["InChI"].append(data["InChI"])
            combined_data["Name"].append(data["Name"])

            if "Synonyms" in data.keys():
                combined_data["Synonyms"].extend(data["Synonyms"])

            id_dict = {
                "nmrshiftdb2.json" : "NMR ShiftDB ID",
                "hmdb.json" : "HMDB ID",
                "golm.json" : "GOLM ID",
                "massbank.json" : "MassBank ID",
                "respect.json" : "Respect ID",
                "spektraris.json" : "Spektraris ID"
            }

            combined_data["Source Information"][id_dict[filename]] = data["Accession"]

            if filename == "hmdb.json":
                combined_data["HMDB Information"]["SMPDB Pathways"] = data["SMPDB Pathways"]
                combined_data["HMDB Information"]["Sources"] = data["Sources"]

            combined_dict[inchikey] = combined_data


    with open(os.path.join(o_fp), "wb") as out_file:
        json.dump(combined_dict, out_file, indent=4)


def count_and_compare(list_item):
    count_dict = {}
    for i in list_item:
        if type(i) == str:
            i = i.title()
            if i not in count_dict.keys():
                count_dict[i] = 1
            else:
                count_dict[i] += 1

    try:
        return max(count_dict.iteritems(), key=operator.itemgetter(1))[0]
    except ValueError:
        return None

def define_names(combined_fp, output_fp):
    with open(combined_fp, "rb") as in_file:
        combined_data = json.load(in_file)


    for inchikey, data in combined_data.items():
        if len(data["Name"]) > 1:
            data["Name"] = count_and_compare(data["Name"])
        else:
            data["Name"] = data["Name"][0]

        if len(data["InChI"]) > 1:
            data["InChI"] = count_and_compare(data["InChI"])
        else:
            data["InChI"] = data["InChI"][0]

        if data["Name"] == None:
            del combined_data[inchikey]

        data["Synonyms"] = list(set(data["Synonyms"]))

    with open(os.path.join(output_fp), "wb") as out_file:
        json.dump(combined_data, out_file, indent=4)

if __name__ == "__main__":

    d_dir = os.path.join(os.path.expanduser("~"), "Data/dimedb/")
    o_dir = os.path.join(d_dir, "output/")

    combined_fp = os.path.join(o_dir, "stripped_sources_names.json")
    final_fp = os.path.join(o_dir, "stripped_sources_final.json")


    if os.path.isdir(o_dir) != True:
        os.makedirs(o_dir)

    if os.path.isfile(combined_fp) != True:
        json_files = []
        for file in os.listdir(d_dir):
            if ".json" in file:
                with open(os.path.join(d_dir, file), "rb") as infile:
                    json_files.append([file, json.load(infile)])
        generate_combined_dictonary(json_files, combined_fp)

    define_names(combined_fp, final_fp)
