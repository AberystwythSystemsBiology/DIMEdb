import json, urllib2, collections, pyidick, requests, re, pickle, os
from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops, Draw
from bioservices import KEGG, KEGGParser
from bson.json_util import dumps as bson_dumps
from biocyc import biocyc
import pybel
from joblib import Parallel, delayed
from tqdm import tqdm

import signal, os

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


directory = "/home/keo7/.data/dimedb/"

def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)

combined = load_json(directory + "combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")
pubchem = load_json(directory + "stripped_pubchem.json")

def identification_info(inchikey):
    id_info = {
        "Name" : None,
        "Synonyms": [],
        "IUPAC Name" : None,
        "Systematic Name" : None,
        "InChI" : None,
        "SMILES" : None,
        "Molecular Formula" : None
    }


    try:
        if combined[inchikey]["HMDB Accession"] != None:
            id_info["Name"] = hmdb[inchikey]["Name"]
            id_info["Synonyms"].extend(hmdb[inchikey]["Synonyms"])

        if combined[inchikey]["PubChem ID"] != None:
            if type(pubchem[inchikey]["Name"]) != list:
                id_info["Name"] = pubchem[inchikey]["Name"]
            else:
                id_info["Name"] = pubchem[inchikey]["Name"][0]

            id_info["Synonyms"].extend(pubchem[inchikey]["Synonyms"])

        if combined[inchikey]["ChEBI ID"] != None:
            if chebi[inchikey]["Name"] != None:
                if type(chebi[inchikey]["Name"]) != list:
                    id_info["Name"] = chebi[inchikey]["Name"]
            id_info["Synonyms"].extend(chebi[inchikey]["Synonyms"])

        id_info["InChI"] = str(combined[inchikey]["InChI"])
        try:
            id_info["SMILES"] = pybel.readstring("inchi", id_info["InChI"]).write("smi").rstrip()
        except IOError:
            rdkit_mol = None

        if id_info["SMILES"] != None:
            rdkit_mol = MolFromSmiles(id_info["SMILES"])
            try:
                id_info["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(rdkit_mol)
            except Exception:
                pass
            if id_info["IUPAC Name"] == None:
                try:
                    response = requests.get("http://cactus.nci.nih.gov/chemical/structure/%(smiles)s/iupac_name"  % dict(smiles=id_info["SMILES"]))
                    if response.status_code == 200:
                        id_info["IUPAC Name"] = response.text
                except Exception:
                    pass

        id_info["Synonyms"] = list(set(id_info["Synonyms"]))

        if id_info["Name"] == None:
            id_info["Molecular Formula"] = None

        return id_info, rdkit_mol
    except KeyError:
        return id_info, None

def taxonomic_properties(inchikey):
    properties = {
        "HMDB" : {
            "Origins" : [],
            "Biofluid Locations" : [],
            "Tissue Locations" : []
        }
    }
    if combined[inchikey]["HMDB Accession"] != None:
        properties["HMDB"] = hmdb[inchikey]["Sources"]

    return properties

def generate_sources(inchikey):
    sources = {
        "Wikidata": None,
        "CAS": None,
        "KEGG Compound": None,
        "Chemspider": None,
        "BioCyc" : None
    }

    try:
        response = urllib2.urlopen("http://webservice.bridgedb.org/Human/xrefs/Ik/" + inchikey)
        for line in response.read().splitlines():
            resource = line.split("\t")
            if resource[1] in sources.keys():
                sources[resource[1]] = resource[0]
    except Exception:
        pass

    if sources["KEGG Compound"] != None:
        response = requests.get("https://websvc.biocyc.org/META/foreignid?ids=Kegg:" + sources["KEGG Compound"])
        if response.status_code == 200:
            biocyc_info = response.text.split("\t")
            if biocyc_info[1] == "1":
                sources["BioCyc"] = biocyc_info[2].replace("\n", "")

    sources.update({
        "ChEBI ID" : combined[inchikey]["ChEBI ID"],
        "PubChem ID" : combined[inchikey]["PubChem ID"],
        "HMDB Accession": combined[inchikey]["HMDB Accession"]
    })

    return sources


def generate_pathways(inchikey, sources):
    pathways = {
        "KEGG" : [],
        "SMPDB" : [],
        "BioCyc" : []
    }

    if sources["HMDB Accession"] != None:
        pathways["SMPDB"] = hmdb[inchikey]["SMPDB Pathways"]

    if sources["BioCyc"] != None:
        o = biocyc.get(str(sources["BioCyc"]))
        if o != None:
            try:
                biocyc_pathways = [re.sub('<[^<]+?>', '', r.pathways[0].biocyc_link_html) for r in o.reactions]
                pathways["BioCyc"] = list(set(biocyc_pathways))
            except Exception:
                pass

    if sources["KEGG Compound"] != None:
        try:
            kegg = KEGG().get(sources["KEGG Compound"])
            if kegg != 404:
                kegg_dict = KEGGParser().parse(kegg)
                pathways["KEGG"] = kegg_dict["PATHWAY"].keys()
        except KeyError, AttributeError:
            pass

    return pathways

def physicochemical(rdkit_mol):
    clogP, mr_values = rdMolDescriptors.CalcCrippenDescriptors(rdkit_mol)

    return {
        "Molecular Weight" : rdMolDescriptors.CalcExactMolWt(rdkit_mol),
        "Ether Oxygens": Fragments.fr_ether(rdkit_mol),
        "Hydroxy Groups": Fragments.fr_Al_OH(rdkit_mol),
        "Carboxylic Acids": Fragments.fr_Al_COO(rdkit_mol),
        "Secondary Amines": Fragments.fr_NH2(rdkit_mol),
        "Formal Charge": rdmolops.GetFormalCharge(rdkit_mol),
        "clogP": clogP,
        "MR Values": mr_values,
        "Fraction of SP3 Carbon": rdMolDescriptors.CalcFractionCSP3(rdkit_mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(rdkit_mol),
        "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(rdkit_mol),
        "Hydrogen Bond Acceptors": rdMolDescriptors.CalcNumHBA(rdkit_mol),
        "Hydrogen Bond Donors": rdMolDescriptors.CalcNumHBD(rdkit_mol),
        "Rings": rdMolDescriptors.CalcNumRings(rdkit_mol),
        "Heavy Atoms": rdkit_mol.GetNumHeavyAtoms(),
        "Polar Surface Area": sum(MolSurf._pyTPSAContribs(rdkit_mol))
    }

def adduct_information(smiles, properties):
    adducts = collections.defaultdict(list)
    try:
        mol = pyidick.Molecule(smiles)
    except Exception:
        mol = None

    def calculate(type, mol, rule_dict=None, electrons=0, charge=0):
        iso_dist = mol.isotopic_distribution(rule_dict=rule_dict, electrons=electrons, charge=charge)

        return {
            "Type": type,
            "Accurate Mass": iso_dist[0][0],
            "Isotopic Distribution": iso_dist
        }

    try:
        # Neutral
        adducts["Neutral"].append(calculate("[M]", mol))
        # Negative
        if properties["Formal Charge"] == -1:
            adducts["Negative"].append(calculate("[M1-.]1-", mol, {"add": {}, "remove": {}}, charge=-1, electrons=-0))
        if properties["Hydrogen Bond Donors"] > 0 and properties["Formal Charge"] == 0:
            adducts["Negative"].append(
                calculate("[3M-H]1-", mol, {"add": {}, "remove": {"H": 1}, "multiply": 3}, charge=-1, electrons=1))
            adducts["Negative"].append(
                calculate("[2M+Hac-H]1-", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}, "multiply": 2}, charge=-1,
                          electrons=1))
            adducts["Negative"].append(
                calculate("[2M+FA-H]1-", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}, "multiply": 2}, charge=-1,
                          electrons=1))
            adducts["Negative"].append(
                calculate("[2M-H]1-", mol, {"add": {}, "remove": {"H": 1}, "multiply": 2}, charge=-1, electrons=1))
            adducts["Negative"].append(
                calculate("[M+TFA-H]1-", mol, {"add": {"C": 2, "O": 2, "F": 3}, "remove": {}}, charge=-1, electrons=1))
            adducts["Negative"].append(
                calculate("[M+Hac-H]1-", mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}}, charge=-1, electrons=1))
            adducts["Negative"].append(
                calculate("[M+FA-H]1-", mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}}, charge=-1, electrons=1))
            adducts["Negative"].append(calculate("[M-H]1-", mol, {"add": {}, "remove": {"H": 1}}, charge=-1, electrons=1))
            if properties["Hydrogen Bond Donors"] > 1:
                adducts["Negative"].append(
                    calculate("[M-2H]2-", mol, {"add": {}, "remove": {"H": 2}}, charge=-2, electrons=2))
            if properties["Hydrogen Bond Donors"] > 2:
                adducts["Negative"].append(
                    calculate("[M-3H]3-", mol, {"add": {}, "remove": {"H": 3}}, charge=-3, electrons=3))
            if properties["Hydrogen Bond Acceptors"] > 0:
                adducts["Negative"].append(
                    calculate("[2M+Na-2H]1-", mol, {"add": {"Na": 1}, "remove": {"H": 2}, "multiply": 2}, charge=-1,
                              electrons=1))
                adducts["Negative"].append(
                    calculate("[M+K-2H]1-", mol, {"add": {"K": 1}, "remove": {"H": 2}}, charge=-1, electrons=1))
            if properties["Hydrogen Bond Donors"] > 1 and properties["Hydrogen Bond Acceptors"] > 0:
                adducts["Negative"].append(
                    calculate("[M+Na-2H]1-", mol, {"add": {"Na": 1}, "remove": {"H": 2}}, charge=-1, electrons=1))
            if properties["Hydrogen Bond Acceptors"] > 0 and properties["Formal Charge"] == 0:
                adducts["Negative"].append(
                    calculate("[M+Br]1-", mol, {"add": {"Br": 1}, "remove": {}}, charge=-1, electrons=1))
                adducts["Negative"].append(
                    calculate("[M+Cl]1-", mol, {"add": {"Cl": 1}, "remove": {}}, charge=-1, electrons=1))
        # Positive
        if properties["Formal Charge"] == 1:
            adducts["Positive"].append(calculate("[M1+.]1+", mol, {"add": {}, "remove": {}}, charge=1, electrons=-0))
        if properties["Formal Charge"] == 0:
            if properties["Hydrogen Bond Acceptors"] > 0:
                adducts["Positive"].append(
                    calculate("[2M+K]1+", mol, {"add": {"K": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[2M+Na]1+", mol, {"add": {"Na": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[2M+NH4]1+", mol, {"add": {"N": 4, "H": 4}, "remove": {}, "multiply": 2}, charge=1,
                              electrons=-1))
                adducts["Positive"].append(
                    calculate("[2M+H]1+", mol, {"add": {"H": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[M+2K-H]1+", mol, {"add": {"K": 2}, "remove": {"H": 1}}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[M+2Na-H]1+", mol, {"add": {"Na": 2}, "remove": {"H": 1}}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[M+K]1+", mol, {"add": {"K": 1}, "remove": {}}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[M+Na]1+", mol, {"add": {"Na": 1}, "remove": {}}, charge=1, electrons=-1))
                adducts["Positive"].append(
                    calculate("[M+H]1+", mol, {"add": {"H": 1}, "remove": {}}, charge=1, electrons=-1))
            if properties["Hydrogen Bond Acceptors"] > 1:
                adducts["Positive"].append(
                    calculate("[2M+3H2O+2H]2+", mol, {"add": {"H": 8, "O": 3}, "remove": {}, "multiply": 2}, charge=2,
                              electrons=-2))
                adducts["Positive"].append(
                    calculate("[2M+3ACN+2H]2+", mol, {"add": {"C": 6, "H": 11, "N": 3}, "remove": {}, "multiply": 2},
                              charge=2, electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+2ACN+2H]2+", mol, {"add": {"C": 3, "H": 8, "N": 2}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+2Na]2+", mol, {"add": {"Na": 2}, "remove": {}}, charge=2, electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+2ACN+2H]2+", mol, {"add": {"C": 2, "H": 5, "N": 1}, "remove": {}}, charge=2,
                              electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+H+K]2+", mol, {"add": {"K": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+H+Na]2+", mol, {"add": {"Na": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+H+NH4]2+", mol, {"add": {"N": 5, "H": 5}, "remove": {}}, charge=2, electrons=-2))
                adducts["Positive"].append(
                    calculate("[M+2H]2+", mol, {"add": {"H": 2}, "remove": {}}, charge=2, electrons=-2))
            if properties["Hydrogen Bond Acceptors"] > 2:
                adducts["Positive"].append(
                    calculate("[M+3Na]3+", mol, {"add": {"Na": 3}, "remove": {}}, charge=3, electrons=-3))
                adducts["Positive"].append(
                    calculate("[M+H+2Na]3+", mol, {"add": {"H": 1, "Na": 2}, "remove": {}}, charge=3, electrons=-3))
                adducts["Positive"].append(
                    calculate("[M+3H]3+", mol, {"add": {"H": 3}, "remove": {}}, charge=3, electrons=-3))
                adducts["Positive"].append(
                    calculate("[M+2H+2a]3+", mol, {"add": {"H": 2, "Na": 1}, "remove": {}}, charge=3, electrons=-3))
    except Exception:
        adducts = {"Positive" : [], "Negative" : [], "Neutral" : []}
    return adducts

def process_compound(inchikey):

    id_info, rdkit_mol = identification_info(inchikey)

    if id_info["Molecular Formula"] != None:
        p_properties = physicochemical(rdkit_mol)
        adducts = adduct_information(id_info["SMILES"], p_properties)
        t_properties = taxonomic_properties(inchikey)
        sources = generate_sources(inchikey)
        pathway_info = generate_pathways(inchikey, sources)

        dimedb_compound = [
            ["_id", inchikey],
            ["Identification Information", id_info],
            ["Physicochemical Properties", p_properties],
            ["Taxonomic Properties", t_properties],
            ["External Sources", sources],
            ["Pathways", pathway_info],
            ["Adducts", adducts]
        ]

        return collections.OrderedDict(dimedb_compound), rdkit_mol

    else:
        return None, None

def generate_image(mol, inchikey):
    Draw.MolToFile(mol, fileName=directory+"structures/"+inchikey+".svg", imageType="svg", size=(250, 250))

if __name__ == "__main__":
    limiter = 10
    inchikeys = combined.keys()
    slice = range(0, len(inchikeys), limiter)[1184:]


    failed_slices = [["From", "To"]]

    for inchikey_index in tqdm(slice):
        try:
            with timeout(30):
                processed_data = Parallel(n_jobs=32)(delayed(process_compound)(id) for id in inchikeys[inchikey_index:inchikey_index+limiter])
                processed_data = [[compound, rdkit_mol] for compound, rdkit_mol in processed_data if compound != None]
                [generate_image(rdkit_mol, compound["_id"]) for compound, rdkit_mol in processed_data if compound != None]
                pickle.dump(processed_data, open(directory+"pickles/"+str(inchikey_index)+".pkl", "wb"))
        except TimeoutError:
            print inchikey_index, "failed"
            failed_slices.append([inchikey_index, inchikey_index+limiter])

    for i in failed_slices:
        print i

    db = []
    for file in os.listdir(directory+"pickles/"):
        processed_data = pickle.load(open(directory+"pickles/"+file, "rb"))
        db.extend([compound for compound, rdkit_mol in processed_data])

    mongodb_file = json.loads(bson_dumps(db), object_pairs_hook=collections.OrderedDict)

    with open(directory+"dimedb.json", "wb") as outfile:
        json.dump(mongodb_file, outfile, indent=4)