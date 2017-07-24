import json, pubchempy, urllib2, collections, pyidick
from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops
from bioservices import KEGG, KEGGParser
from bson.json_util import dumps as bson_dumps
import pybel

directory = "/home/keo7/.data/dimedb/"


def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)

combined = load_json(directory + "/combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")

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


    id_info["InChI"] = str(combined[inchikey]["InChI"])

    id_info["SMILES"] = pybel.readstring("inchi", id_info["InChI"]).write("smi").rstrip()

    rdkit_mol = MolFromSmiles(id_info["SMILES"])
    try:
        id_info["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(rdkit_mol)
    except Exception:
        print "Exception caught"

    if combined[inchikey]["HMDB Accession"] != None:
        id_info["Name"] = hmdb[inchikey]["Name"]
        id_info["Synonyms"] = hmdb[inchikey]["Synonyms"]

    if combined[inchikey]["ChEBI ID"] != None:
        if id_info["Name"] == None:
            id_info["Name"] = chebi[inchikey]["Name"]
        if id_info["IUPAC Name"] == None:
            id_info["IUPAC Name"] = chebi[inchikey]["IUPAC Name"]
        if id_info["Synonyms"] == []:
            id_info["Synonyms"] = chebi[inchikey]["Synonyms"]

    if combined[inchikey]["PubChem ID"] != None:
        compound = pubchempy.get_compounds(combined[inchikey]["PubChem ID"])[0]
        if id_info["Name"] == None:
            id_info["Name"] = compound.synonyms[0]
        if id_info["Synonyms"] == []:
            id_info["Synonyms"] = compound.synonyms[1:]
        if id_info["IUPAC Name"] == None:
            id_info["IUPAC Name"] = compound.iupac_name

    return id_info, rdkit_mol

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
        "Chemspider": None
    }


    response = urllib2.urlopen("http://webservice.bridgedb.org/Human/xrefs/Ik/" + inchikey)
    for line in response.read().splitlines():
        resource = line.split("\t")
        if resource[1] in sources.keys():
            sources[resource[1]] = resource[0]

    sources.update({
        "ChEBI ID" : combined[inchikey]["ChEBI ID"],
        "PubChem ID" : combined[inchikey]["PubChem ID"],
        "HMDB Accession": combined[inchikey]["HMDB Accession"]
    })

    return sources


def generate_pathways(inchikey, sources):
    pathways = {
        "KEGG" : [],
        "SMPDB" : []
    }
    if sources["HMDB Accession"] != None:
        pathways["SMPDB"] = hmdb[inchikey]["SMPDB Pathways"]

    if sources["KEGG Compound"] != None:
        try:
            kegg = KEGG().get(sources["KEGG Compound"])
            if kegg != 404:
                kegg_dict = KEGGParser().parse(kegg)
                pathways["KEGG"] = kegg_dict["PATHWAY"].keys()
        except KeyError, AttributeError:
            pass
    return pathways

def physiochemical(rdkit_mol):
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

    mol = pyidick.Molecule(smiles)

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
    except IndexError:
        adducts = {"Positive" : [], "Negative" : [], "Neutral" : []}
    return adducts

def process_compound(inchikey):
    id_info, rdkit_mol = identification_info(inchikey)

    if id_info["Molecular Formula"] != None:
        p_properties = physiochemical(rdkit_mol)
        t_properties = taxonomic_properties(inchikey)
        sources = generate_sources(inchikey)
        pathway_info = generate_pathways(inchikey, sources)
        adducts = adduct_information(id_info["SMILES"], p_properties)

        dimedb_compound = [
            ["_id", inchikey],
            ["Identification Information", id_info],
            ["Physiochemical Properties", p_properties],
            ["Taxonomic Properties", t_properties],
            ["External Sources", sources],
            ["Pathways", pathway_info],
            ["Adducts", adducts]
        ]

        return collections.OrderedDict(dimedb_compound)

    else:
        return None

if __name__ == "__main__":
    db = []

    test_keys = combined.keys()[:100]

    for index, inchikey in enumerate(test_keys):
        print index, "/", len(test_keys)
        dimedb_compound = process_compound(inchikey)
        if dimedb_compound != None:
            db.extend([dimedb_compound])

    mongodb_file = json.loads(bson_dumps(db), object_pairs_hook=collections.OrderedDict)

    with open(directory+"dimedb.json", "wb") as outfile:
        json.dump(mongodb_file, outfile, indent=4)
