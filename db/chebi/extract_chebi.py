import openbabel, csv, json, collections, os
from tqdm import tqdm

def inchi_to_inchikey(inchi):
    try:
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("inchi", "inchikey")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, inchi)
        return obConversion.WriteString(mol).rstrip()
    except Exception, Warning:
        return None

def read_inchi_tsv(file_name):
    compounds = {}
    with open(file_name) as chebi_inchi_file:
        tsv_reader = csv.reader(chebi_inchi_file, delimiter="\t")
        next(tsv_reader) # SKIP HEADER
        for index, row in enumerate(tsv_reader):
            if index == 4023: #No idea.
                break
            inchi_key = inchi_to_inchikey(row[1])
            if inchi_key != None:
                compounds[inchi_key] = {
                    "InChI": row[1],
                    "ChEBI ID" : row[0]
                }
            else:
                print row[1], "not found?"
    return compounds

def read_inich_names(file_name):
    names = collections.defaultdict(list)
    with open(file_name) as chebi_inchi_file:
        tsv_reader = csv.reader(chebi_inchi_file, delimiter="\t")
        next(tsv_reader)  # SKIP HEADER
        for row in tqdm(tsv_reader):
            if row[1] not in names.keys():
                names[row[1]] = {
                    "Synonyms" : [],
                    "IUPAC Name" : None,
                    "Name" : None
                }
            if row[2] == "SYNONYM":
                names[row[1]]["Synonyms"].append(row[4])
            elif row[2] == "IUPAC NAME":
                names[row[1]]["IUPAC Name"] = row[4]
            elif row[2] == "NAME":
                names[row[1]]["Name"] = row[4]


    return names

if __name__ == "__main__":
    directory = "/home/keo7/.data/dimedb/"

    if os.path.isfile(directory+"chebi/chebi_names.json"):
        with open(directory+"chebi/chebi_names.json", "rb") as name_file:
            names = json.load(name_file)
    else:
        names = read_inich_names(directory+"chebi/names.tsv")
        with open(directory+"chebi/chebi_names.json", "wb") as outfile:
            json.dump(names, outfile, indent=4)

    compounds = read_inchi_tsv(directory+"chebi/chebiId_inchi.tsv")

    for inchikey in compounds.keys():
        try:
            name_info = names[compounds[inchikey]["CHEBI_ID"]]
            compounds[inchikey]["Name"] = name_info["Name"]
            compounds[inchikey]["IUPAC Name"] = name_info["IUPAC Name"]
            compounds[inchikey]["Synonyms"] = name_info["Synonyms"]
        except KeyError:
            compounds[inchikey]["Name"] = None,
            compounds[inchikey]["IUPAC Name"] = None
            compounds[inchikey]["Synonyms"] = []

    with open(directory+"stripped_chebi.json", "wb") as outfile:
        json.dump(compounds, outfile, indent=4)
