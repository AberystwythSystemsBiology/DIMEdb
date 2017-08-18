import json, pybel, requests

from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops, Draw


def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)


directory = "/home/keo7/.data/dimedb/"

combined = load_json(directory + "combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")
pubchem = load_json(directory + "stripped_pubchem.json")


class Metabolite(object):
    def __init__(self, inchikey):
        self.inchikey = inchikey
        self.combined_info = combined[inchikey]

        self.inchi = self.combined_info["InChI"]
        self.smiles = pybel.readstring("inchi", str(self.inchi)).write("smi").rstrip()
        self.rdkit_mol = MolFromSmiles(self.smiles)

    @property
    def identification_information(self):
        id_info = {
            "Name": None,
            "Synonyms": [],
            "IUPAC Name": self._get_iupac_name(),
            "Systematic Name": None,
            "InChI": self.inchi,
            "SMILES": self.smiles,
            "Molecular Formula": rdMolDescriptors.CalcMolFormula(self.rdkit_mol)
        }

        if self.combined_info["HMDB Accession"] != None:
            if combined[self.inchikey]["HMDB Accession"] != None:
                id_info["Name"] = hmdb[self.inchikey]["Name"]
                id_info["Synonyms"].extend(hmdb[self.inchikey]["Synonyms"])

            if combined[inchikey]["PubChem ID"] != None:
                if type(pubchem[self.inchikey]["Name"]) != list:
                    id_info["Name"] = pubchem[self.inchikey]["Name"]
                else:
                    id_info["Name"] = pubchem[self.inchikey]["Name"][0]

                id_info["Synonyms"].extend(pubchem[self.inchikey]["Synonyms"])

            if combined[self.inchikey]["ChEBI ID"] != None:
                if chebi[self.inchikey]["Name"] != None:
                    if type(chebi[self.inchikey]["Name"]) != list:
                        id_info["Name"] = chebi[self.inchikey]["Name"]
                id_info["Synonyms"].extend(chebi[self.inchikey]["Synonyms"])

        return id_info

    def _get_iupac_name(self):
        response = requests.get("http://cactus.nci.nih.gov/chemical/structure/%s/iupac_name" % self.smiles)
        if response.status_code == 200:
            return response.text
        else:
            return None

    @property
    def sources(self):
        pass

    def generate_image(self, fp, s=[250, 250]):
        Draw.MolToFile(self.rdkit_mol, fileName=fp, imageType="svg", size=s)


if __name__ == "__main__":

    for inchikey in combined.keys():
        metabolite = Metabolite(inchikey)
        print metabolite.identification_information
        exit(0)