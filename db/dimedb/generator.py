import json, requests, re, signal, pyidick, collections, urllib, os, warnings, pubchempy
from subprocess import check_output
warnings.filterwarnings("ignore")

from bson.json_util import dumps as bson_dumps
from tqdm import tqdm
from joblib import Parallel, delayed
from bioservices import KEGG, KEGGParser
from biocyc import biocyc
from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops, Draw

VERSION = "1.0NOV2017"

output_dir = os.path.join(os.path.expanduser("~"), "Data/dimedb/output/")

final_dir = os.path.join(output_dir, "dimedb_jsons/")

structures_dir = os.path.join(output_dir, "structures/")

if os.path.isdir(structures_dir) != True:
    os.makedirs(structures_dir)

if os.path.isdir(final_dir) != True:
    os.makedirs(final_dir)

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
    def __init__(self, inchikey, data):
        self.inchikey = inchikey
        try:
            self.combined_info = data
        except KeyError:
            raise SMILESerror()

        self.inchi = self.combined_info["InChI"]
        self.smiles = self.get_smiles()

        self.rdkit_mol = MolFromSmiles(self.smiles)
        self.identification_information = self.get_identification_information()
        self.external_sources = self.get_external_sources()
        self.physicochemical_properties = self.get_physicochemical_properties()
        self.pathways = self.get_pathways()
        self.taxonomic_properties = self.get_taxonomic_properties()
        self.adduct_information = self.get_adduct_information()


    def get_smiles(self):
        with open(os.devnull, 'w') as fp:
            val = check_output('obabel -iinchi -:"%s" -osmi' % self.inchi, shell=True, stderr=open(os.devnull, 'w')).decode().strip()
        if val != "":
            return str(val)
        else:
            raise SMILESerror()


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



        id_info["Name"] = self.combined_info["Name"]
        id_info["Synonyms"] = self.combined_info["Synonyms"]

        if id_info["Name"] == None:
            raise SMILESerror()

        return id_info

    def _get_iupac_name(self):
        response = requests.get("http://cactus.nci.nih.gov/chemical/structure/%s/iupac_name" % urllib.quote_plus(self.smiles))
        if response.status_code == 200:
            return response.text
        else:
            return None

    def get_external_sources(self):
        sources = {
            "Wikidata": None,
            "CAS": None,
            "KEGG Compound": None,
            "KEGG Drug" : None,
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
            "Spektraris": self.combined_info["Source Information"]["Spektraris ID"],
            "NMR Shift DB": self.combined_info["Source Information"]["NMR ShiftDB ID"],
            "MassBank" : self.combined_info["Source Information"]["MassBank ID"],
            "Respect" : self.combined_info["Source Information"]["Respect ID"],
            "GOLM" : self.combined_info["Source Information"]["GOLM ID"],
            "HMDB Accession": self.combined_info["Source Information"]["HMDB ID"]
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
            pathways["SMPDB"] = self.combined_info["HMDB Information"]["SMPDB Pathways"]

        if self.external_sources["BioCyc"] != None:
            biocyc_object = biocyc.get(str(self.external_sources["BioCyc"]))
            if biocyc_object != None:
                biocyc_pathways = []
                for reaction in biocyc_object.reactions:
                    try:
                        reaction_pathways = reaction.pathways
                        for pathway in reaction_pathways:
                            biocyc_pathways.append(re.sub('<[^<]+?>', '', pathway.biocyc_link_html))
                    except AttributeError:
                        biocyc_pathways = []
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
        if self.external_sources["HMDB Accession"] != None:
            properties["HMDB"] = self.combined_info["HMDB Information"]["Sources"]

        return properties

    def get_adduct_information(self):

        adducts = []
        def calculate(type, polarity, mol, rule_dict=None, electrons=0, charge=0):
            try:
                with timeout(30):
                    iso_dist = mol.isotopic_distribution(rule_dict=rule_dict, electrons=electrons, charge=charge)

                    adducts.append({
                        "Polarity": polarity,
                        "Adduct": type,
                        "Accurate Mass": iso_dist[0][0],
                        "Isotopic Distribution": iso_dist
                    })

            except Exception:
                pass

        try:
            mol = pyidick.Molecule(self.smiles)
        except ValueError:
            return []

        calculate("[M]", "Neutral", mol)

        # Negative
        if self.physicochemical_properties["Formal Charge"] == -1:
            calculate("[M1-.]1-", "Negative", mol, {"add": {}, "remove": {}}, charge=-1, electrons=-0)
        if self.physicochemical_properties["Hydrogen Bond Donors"] > 0 and self.physicochemical_properties["Formal Charge"] == 0:
            calculate("[3M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}, "multiply": 3}, charge=-1,
                      electrons=1)
            calculate("[2M+Hac-H]1-", "Negative", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}, "multiply": 2},
                      charge=-1, electrons=1)
            calculate("[2M+FA-H]1-", "Negative", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}, "multiply": 2}, charge=-1, electrons=1)
            calculate("[2M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}, "multiply": 2}, charge=-1, electrons=1)
            calculate("[M+TFA-H]1-", "Negative", mol, {"add": {"C": 2, "O": 2, "F": 3}, "remove": {}}, charge=-1, electrons=1)
            calculate("[M+Hac-H]1-", "Negative", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}}, charge=-1, electrons=1)
            calculate("[M+FA-H]1-", "Negative", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}}, charge=-1, electrons=1)
            calculate("[M-H]1-", "Negative", mol, {"add": {}, "remove": {"H": 1}}, charge=-1, electrons=1)
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 1:
                calculate("[M-2H]2-", "Negative", mol, {"add": {}, "remove": {"H": 2}}, charge=-2, electrons=2)
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 2:
                calculate("[M-3H]3-", "Negative", mol, {"add": {}, "remove": {"H": 3}}, charge=-3, electrons=3)
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0:
                calculate("[2M+Na-2H]1-", "Negative", mol, {"add": {"Na": 1}, "remove": {"H": 2}, "multiply": 2}, charge=-1, electrons=1)
                calculate("[M+K-2H]1-", "Negative", mol, {"add": {"K": 1}, "remove": {"H": 2}}, charge=-1, electrons=1)
            if self.physicochemical_properties["Hydrogen Bond Donors"] > 1 and self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0:
                calculate("[M+Na-2H]1-", "Negative", mol, {"add": {"Na": 1}, "remove": {"H": 2}}, charge=-1, electrons=1)
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0 and self.physicochemical_properties["Formal Charge"] == 0:
                calculate("[M+Br]1-", "Negative", mol, {"add": {"Br": 1}, "remove": {}}, charge=-1, electrons=1)
                calculate("[M+Cl]1-", "Negative", mol, {"add": {"Cl": 1}, "remove": {}}, charge=-1, electrons=1)

        # Positive
        if self.physicochemical_properties["Formal Charge"] == 1:
            calculate("[M1+.]1+", "Positive", mol, {"add": {}, "remove": {}}, charge=1, electrons=-0)
        if self.physicochemical_properties["Formal Charge"] == 0:
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 0:
                calculate("[2M+K]1+", "Positive", mol, {"add": {"K": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1)
                calculate("[2M+Na]1+", "Positive", mol, {"add": {"Na": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1)
                calculate("[2M+NH4]1+", "Positive", mol, {"add": {"N": 4, "H": 4}, "remove": {}, "multiply": 2}, charge=1, electrons=-1)
                calculate("[2M+H]1+", "Positive", mol, {"add": {"H": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1)
                calculate("[M+2K-H]1+", "Positive", mol, {"add": {"K": 2}, "remove": {"H": 1}}, charge=1, electrons=-1)
                calculate("[M+2Na-H]1+", "Positive", mol, {"add": {"Na": 2}, "remove": {"H": 1}}, charge=1, electrons=-1)
                calculate("[M+K]1+", "Positive", mol, {"add": {"K": 1}, "remove": {}}, charge=1, electrons=-1)
                calculate("[M+Na]1+", "Positive", mol, {"add": {"Na": 1}, "remove": {}}, charge=1, electrons=-1)
                calculate("[M+H]1+", "Positive", mol, {"add": {"H": 1}, "remove": {}}, charge=1, electrons=-1)
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 1:
                calculate("[2M+3H2O+2H]2+", "Positive", mol, {"add": {"H": 8, "O": 3}, "remove": {}, "multiply": 2}, charge=2, electrons=-2)
                calculate("[2M+3ACN+2H]2+", "Positive", mol, {"add": {"C": 6, "H": 11, "N": 3}, "remove": {}, "multiply": 2}, charge=2, electrons=-2)
                calculate("[M+2ACN+2H]2+", "Positive", mol, {"add": {"C": 3, "H": 8, "N": 2}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+2Na]2+", "Positive", mol, {"add": {"Na": 2}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+2ACN+2H]2+", "Positive", mol, {"add": {"C": 2, "H": 5, "N": 1}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+H+K]2+", "Positive", mol, {"add": {"K": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+H+Na]2+", "Positive", mol, {"add": {"Na": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+H+NH4]2+", "Positive", mol, {"add": {"N": 5, "H": 5}, "remove": {}}, charge=2, electrons=-2)
                calculate("[M+2H]2+", "Positive", mol, {"add": {"H": 2}, "remove": {}}, charge=2, electrons=-2)
            if self.physicochemical_properties["Hydrogen Bond Acceptors"] > 2:
                calculate("[M+3Na]3+", "Positive", mol, {"add": {"Na": 3}, "remove": {}}, charge=3, electrons=-3)
                calculate("[M+H+2Na]3+", "Positive", mol, {"add": {"H": 1, "Na": 2}, "remove": {}}, charge=3, electrons=-3)
                calculate("[M+3H]3+", "Positive", mol, {"add": {"H": 3}, "remove": {}}, charge=3, electrons=-3)
                calculate("[M+2H+2a]3+", "Positive", mol, {"add": {"H": 2, "Na": 1}, "remove": {}}, charge=3, electrons=-3)

        return adducts


    def generate_image(self, inchikey, s=[250, 250]):
        fp = os.path.join(structures_dir, inchikey+".svg")
        Draw.MolToFile(self.rdkit_mol, fileName=fp, imageType="svg", size=s)


    def to_dict(self):
        compound = [
            ["_id", self.inchikey],
            ["Version", VERSION],
            ["Identification Information", self.identification_information],
            ["Physicochemical Properties", self.physicochemical_properties],
            ["Taxonomic Properties", self.taxonomic_properties],
            ["External Sources", self.external_sources],
            ["Pathways", self.pathways],
            ["Adducts", self.adduct_information]
        ]

        return collections.OrderedDict(compound)


def handler(inchikey, data):
    try:
        metabolite = Metabolite(inchikey, data)
        metabolite.generate_image(inchikey)
        metabolite_dict = metabolite.to_dict()
    except (SMILESerror, MolecularFormulaError) as e:
        metabolite_dict = None

    return metabolite_dict


if __name__ == "__main__":
    limiter = 500

    combined_fp = os.path.join(output_dir, "stripped_sources_final.json")

    with open(combined_fp, "r") as hmdb_file:
        combined = json.load(hmdb_file)

    inchikeys = combined.keys()
    slice = range(0, len(inchikeys), limiter)

    for inchikey_index in tqdm(slice[0:]):
        processed_data = Parallel(n_jobs=6)(delayed(handler)(id, combined[id]) for id in inchikeys[inchikey_index:inchikey_index + limiter])
        processed_data = [x for x in processed_data if x != None]
        mongodb_file = json.loads(bson_dumps(processed_data), object_pairs_hook=collections.OrderedDict)
        with open(os.path.join(final_dir, "dimedb_s"+str(inchikey_index)+".json"), "wb") as outfile:
            json.dump(mongodb_file, outfile, indent=4)
        break
