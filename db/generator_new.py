import json, pybel, requests, re, signal, pyidick, collections, pickle, os
from bson.json_util import dumps as bson_dumps
from tqdm import tqdm
from joblib import Parallel, delayed

from bioservices import KEGG, KEGGParser
from biocyc import biocyc
from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops, Draw


def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)


directory = "/home/keo7/.data/dimedb/"

combined = load_json(directory + "combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")
pubchem = load_json(directory + "stripped_pubchem.json")

class SMILESerror(Exception):
    def __init__(self):
        super(SMILESerror, self).__init__("Unable to generate SMILES")

class MolecularFormulaError(Exception):
    def __init__(self):
        super(MolecularFormulaError, self).__init__("Unable to generate Molecular Formula")


class TimeoutError(Exception):
    pass

class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

class Metabolite(object):
    def __init__(self, inchikey):
        self.inchikey = inchikey
        self.combined_info = combined[inchikey]
        self.inchi = self.combined_info["InChI"]
        try:
            self.smiles = pybel.readstring("inchi", str(self.inchi)).write("smi").rstrip()
        except IOError:
            raise SMILESerror()
        self.rdkit_mol = MolFromSmiles(self.smiles)

        self.identification_information = self.get_identification_information()
        self.external_sources = self.get_external_sources()
        self.physicochemical_properties = self.get_physicochemical_properties()
        self.pathways = self.get_pathways()
        self.taxonomic_properties = self.get_taxonomic_properties()
        self.adduct_information = self.get_adduct_information()

    def get_identification_information(self):
        id_info = {
            "Name": None,
            "Synonyms": [],
            "IUPAC Name": self._get_iupac_name(),
            "Systematic Name": None,
            "InChI": self.inchi,
            "SMILES": self.smiles,
            "Molecular Formula": None
        }

        try:
            id_info["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(self.rdkit_mol)
        except Exception:
            raise MolecularFormulaError()

        if self.combined_info["HMDB Accession"] != None:
            if combined[self.inchikey]["HMDB Accession"] != None:
                id_info["Name"] = hmdb[self.inchikey]["Name"]
                id_info["Synonyms"].extend(hmdb[self.inchikey]["Synonyms"])

            if combined[self.inchikey]["PubChem ID"] != None:
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

    def get_external_sources(self):
        sources = {
            "Wikidata": None,
            "CAS": None,
            "KEGG Compound": None,
            "Chemspider": None,
            "BioCyc": None,
            "ChEBI" : None,
            "PubChem" : None,
            "HMDB Accession" : None
        }

        response = requests.get("http://webservice.bridgedb.org/Human/xrefs/Ik/%s" % self.inchikey)
        if response.status_code == 200:
            for line in response.text.splitlines():
                resource = line.split("\t")
                if resource[1] in sources.keys():
                    sources[resource[1]] = resource[0]

        if sources["KEGG Compound"] != None:
            response = requests.get("https://websvc.biocyc.org/META/foreignid?ids=Kegg:%s" % sources["KEGG Compound"])
            if response.status_code == 200:
                biocyc_info = response.text.split("\t")
                if biocyc_info[1] == "1":
                    sources["BioCyc"] = biocyc_info[2].replace("\n", "")

        sources.update({
            "ChEBI": combined[self.inchikey]["ChEBI ID"],
            "PubChem": combined[self.inchikey]["PubChem ID"],
            "HMDB Accession": combined[self.inchikey]["HMDB Accession"]
        })

        return sources

    def get_physicochemical_properties(self):
        clogP, mr_values = rdMolDescriptors.CalcCrippenDescriptors(self.rdkit_mol)

        return {
            "Molecular Weight": rdMolDescriptors.CalcExactMolWt(self.rdkit_mol),
            "Ether Oxygens": Fragments.fr_ether(self.rdkit_mol),
            "Hydroxy Groups": Fragments.fr_Al_OH(self.rdkit_mol),
            "Carboxylic Acids": Fragments.fr_Al_COO(self.rdkit_mol),
            "Secondary Amines": Fragments.fr_NH2(self.rdkit_mol),
            "Formal Charge": rdmolops.GetFormalCharge(self.rdkit_mol),
            "clogP": clogP,
            "MR Values": mr_values,
            "Fraction of SP3 Carbon": rdMolDescriptors.CalcFractionCSP3(self.rdkit_mol),
            "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(self.rdkit_mol),
            "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(self.rdkit_mol),
            "Hydrogen Bond Acceptors": rdMolDescriptors.CalcNumHBA(self.rdkit_mol),
            "Hydrogen Bond Donors": rdMolDescriptors.CalcNumHBD(self.rdkit_mol),
            "Rings": rdMolDescriptors.CalcNumRings(self.rdkit_mol),
            "Heavy Atoms": self.rdkit_mol.GetNumHeavyAtoms(),
            "Polar Surface Area": sum(MolSurf._pyTPSAContribs(self.rdkit_mol))
        }

    def get_pathways(self):
        pathways = {
            "KEGG": [],
            "SMPDB": [],
            "BioCyc": []
        }

        if self.external_sources["HMDB Accession"] != None:
            pathways["SMPDB"] = hmdb[self.inchikey]["SMPDB Pathways"]

        if self.external_sources["BioCyc"] != None:
            biocyc_object = biocyc.get(str(self.external_sources["BioCyc"]))
            if biocyc_object != None:
                biocyc_pathways = []
                for reaction in biocyc_object.reactions:
                    reaction_pathways = reaction.pathways
                    for pathway in reaction_pathways:
                        biocyc_pathways.append(re.sub('<[^<]+?>', '', pathway.biocyc_link_html))
                pathways["BioCyc"] = list(set(biocyc_pathways))


        if self.external_sources["KEGG Compound"] != None:
            kegg_compound = KEGG().get(self.external_sources["KEGG Compound"])
            if kegg_compound != 404:
                kegg_dict = KEGGParser().parse(kegg_compound)
                if "PATHWAY" in kegg_dict.keys():
                    pathways["KEGG"] = kegg_dict["PATHWAY"].keys()

        return pathways

    def get_taxonomic_properties(self):
        properties = {
            "HMDB": {
                "Origins": [],
                "Biofluid Locations": [],
                "Tissue Locations": []
            }
        }
        if self.combined_info["HMDB Accession"] != None:
            properties["HMDB"] = hmdb[self.inchikey]["Sources"]

        return properties

    def get_adduct_information(self):

        adducts = []
        def calculate(type, polarity, mol, rule_dict=None, electrons=0, charge=0):
            try:
                with timeout(10):
                    iso_dist = mol.isotopic_distribution(rule_dict=rule_dict, electrons=electrons, charge=charge)

                    adducts.append({
                        "Polarity": polarity,
                        "Adduct": type,
                        "Accurate Mass": iso_dist[0][0],
                        "Isotopic Distribution": iso_dist
                    })
            except (RuntimeError, TimeoutError) as e:
                pass


        mol = pyidick.Molecule(self.smiles)

        calculate("[M]", "Neutral", mol)

        if self.physicochemical_properties["Formal Charge"] == -1:
            calculate("[M1-.]1-", "Negative", mol, {"add": {}, "remove": {}}, charge=-1, electrons=-0)
        if self.physicochemical_properties["Hydrogen Bond Donors"] > 0 and self.physicochemical_properties["Formal Charge"] == 0:
            calculate("[3M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}, "multiply": 3}, charge=-1,
                      electrons=1)
            calculate("[2M+Hac-H]1-", "Negative", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}, "multiply": 2},
                      charge=-1, electrons=1)
            calculate("[2M+FA-H]1-", "Negative", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}, "multiply": 2}, charge=-1, electrons=1)
            calculate("[2M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}, "multiply": 2}, charge=-1, electrons=1)
            calculate("[M+TFA-H]1-", "Negative", mol, {"add": {"C": 2, "O": 2, "F": 3}, "remove": {}},
                          charge=-1, electrons=1)
            calculate("[M+Hac-H]1-", "Negative", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}},
                          charge=-1, electrons=1)
            calculate("[M+FA-H]1-", "Negative", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}},
                          charge=-1, electrons=1)
            adducts.append(
                calculate("[M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}}, charge=-1, electrons=1))
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 1:
                adducts.append(
                    calculate("[M-2H]2-", "Negative", mol, {"add": {}, "remove": {"H": 2}}, charge=-2,
                              electrons=2))
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 2:
                adducts.append(
                    calculate("[M-3H]3-", "Negative", mol, {"add": {}, "remove": {"H": 3}}, charge=-3,
                              electrons=3))
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0:
                adducts.append(
                    calculate("[2M+Na-2H]1-", "Negative", mol,
                              {"add": {"Na": 1}, "remove": {"H": 2}, "multiply": 2}, charge=-1,
                              electrons=1))
                adducts.append(
                    calculate("[M+K-2H]1-", "Negative", mol, {"add": {"K": 1}, "remove": {"H": 2}}, charge=-1,
                              electrons=1))
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 1 and self.physicochemical_properties[
                "Hydrogen Bond Acceptors"] > 0:
                adducts.append(
                    calculate("[M+Na-2H]1-", "Negative", mol, {"add": {"Na": 1}, "remove": {"H": 2}}, charge=-1,
                              electrons=1))
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0 and \
                            self.physicochemical_properties["Formal Charge"] == 0:
                adducts.append(
                    calculate("[M+Br]1-", "Negative", mol, {"add": {"Br": 1}, "remove": {}}, charge=-1,
                              electrons=1))
                adducts.append(
                    calculate("[M+Cl]1-", "Negative", mol, {"add": {"Cl": 1}, "remove": {}}, charge=-1,
                              electrons=1))
                # Positive
        if self.physicochemical_properties["Formal Charge"] == 1:
            adducts.append(
                calculate("[M1+.]1+", "Positive", mol, {"add": {}, "remove": {}}, charge=1, electrons=-0))
        if self.physicochemical_properties["Formal Charge"] == 0:
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0:
                adducts.append(
                    calculate("[2M+K]1+", "Positive", mol, {"add": {"K": 1}, "remove": {}, "multiply": 2},
                              charge=1, electrons=-1))
                adducts.append(
                    calculate("[2M+Na]1+", "Positive", mol, {"add": {"Na": 1}, "remove": {}, "multiply": 2},
                              charge=1, electrons=-1))
                adducts.append(
                    calculate("[2M+NH4]1+", "Positive", mol,
                              {"add": {"N": 4, "H": 4}, "remove": {}, "multiply": 2}, charge=1,
                              electrons=-1))
                adducts.append(
                    calculate("[2M+H]1+", "Positive", mol, {"add": {"H": 1}, "remove": {}, "multiply": 2},
                              charge=1, electrons=-1))
                adducts.append(
                    calculate("[M+2K-H]1+", "Positive", mol, {"add": {"K": 2}, "remove": {"H": 1}}, charge=1,
                              electrons=-1))
                adducts.append(
                    calculate("[M+2Na-H]1+", "Positive", mol, {"add": {"Na": 2}, "remove": {"H": 1}}, charge=1,
                              electrons=-1))
                adducts.append(
                    calculate("[M+K]1+", "Positive", mol, {"add": {"K": 1}, "remove": {}}, charge=1,
                              electrons=-1))
                adducts.append(
                    calculate("[M+Na]1+", "Positive", mol, {"add": {"Na": 1}, "remove": {}}, charge=1,
                              electrons=-1))
                adducts.append(
                    calculate("[M+H]1+", "Positive", mol, {"add": {"H": 1}, "remove": {}}, charge=1,
                              electrons=-1))
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 1:
                adducts.append(
                    calculate("[2M+3H2O+2H]2+", "Positive", mol,
                              {"add": {"H": 8, "O": 3}, "remove": {}, "multiply": 2}, charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[2M+3ACN+2H]2+", "Positive", mol,
                              {"add": {"C": 6, "H": 11, "N": 3}, "remove": {}, "multiply": 2},
                              charge=2, electrons=-2))
                adducts.append(
                    calculate("[M+2ACN+2H]2+", "Positive", mol, {"add": {"C": 3, "H": 8, "N": 2}, "remove": {}},
                              charge=2,
                              electrons=-2))

                adducts.append(
                    calculate("[M+2Na]2+", "Positive", mol, {"add": {"Na": 2}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[M+2ACN+2H]2+", "Positive", mol, {"add": {"C": 2, "H": 5, "N": 1}, "remove": {}},
                              charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[M+H+K]2+", "Positive", mol, {"add": {"K": 1, "H": 1}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[M+H+Na]2+", "Positive", mol, {"add": {"Na": 1, "H": 1}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[M+H+NH4]2+", "Positive", mol, {"add": {"N": 5, "H": 5}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts.append(
                    calculate("[M+2H]2+", "Positive", mol, {"add": {"H": 2}, "remove": {}}, charge=2,
                              electrons=-2))
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 2:
                adducts.append(
                    calculate("[M+3Na]3+", "Positive", mol, {"add": {"Na": 3}, "remove": {}}, charge=3,
                              electrons=-3))
                adducts.append(
                    calculate("[M+H+2Na]3+", "Positive", mol, {"add": {"H": 1, "Na": 2}, "remove": {}},
                              charge=3, electrons=-3))
                adducts.append(
                    calculate("[M+3H]3+", "Positive", mol, {"add": {"H": 3}, "remove": {}}, charge=3,
                              electrons=-3))
                adducts.append(
                    calculate("[M+2H+2a]3+", "Positive", mol, {"add": {"H": 2, "Na": 1}, "remove": {}},
                              charge=3, electrons=-3))

        return adducts


    def generate_image(self, fp, s=[250, 250]):
        Draw.MolToFile(self.rdkit_mol, fileName=fp, imageType="svg", size=s)


    def to_dict(self):
        compound = [
            ["_id", self.inchikey],
            ["Identification Information", self.identification_information],
            ["Physicochemical Properties", self.physicochemical_properties],
            ["Taxonomic Properties", self.taxonomic_properties],
            ["External Sources", self.external_sources],
            ["Pathways", self.pathways],
            ["Adducts", self.adduct_information]
        ]

        return collections.OrderedDict(compound)


def handler(inchikey):
    try:
        metabolite = Metabolite(inchikey)
        metabolite.generate_image(directory + "structures/" + metabolite.inchikey + ".svg")
        metabolite_dict = metabolite.to_dict()
    except (SMILESerror, MolecularFormulaError) as e:
        metabolite_dict = None

    return metabolite_dict


if __name__ == "__main__":


    limiter = 1000
    inchikeys = combined.keys()
    slice = range(0, len(inchikeys), limiter)

    for inchikey_index in tqdm(slice):
        processed_data = Parallel(n_jobs=32)(delayed(handler)(id) for id in inchikeys[inchikey_index:inchikey_index + limiter])
        processed_data = [x for x in processed_data if x != None]
        pickle.dump(processed_data, open(directory + "pickles/" + str(inchikey_index) + ".pkl", "wb"))
        break

    db = []

    for file in os.listdir(directory+"pickles/"):
        processed_data = pickle.load(open(directory+"pickles/"+file, "rb"))
        db.extend([metabolite for metabolite in processed_data])

    mongodb_file = json.loads(bson_dumps(db), object_pairs_hook=collections.OrderedDict)
    with open(directory+"dimedb.json", "wb") as outfile:
        json.dump(mongodb_file, outfile, indent=4)
