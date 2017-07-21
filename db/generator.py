import json, pubchempy, urllib2
from rdkit.Chem import rdMolDescriptors, MolFromSmiles, MolSurf, Fragments, rdmolops
from bioservices import KEGG, KEGGParser

directory = "/home/keo7/.data/dimedb/"


def load_json(fp):
    with open(fp, "r") as hmdb_file:
        return json.load(hmdb_file)

combined = load_json(directory + "/combined_data.json")
chebi = load_json(directory + "stripped_chebi.json")
hmdb = load_json(directory + "stripped_hmdb.json")

def identification_info(inchikey):
    import pybel

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

    id_info["SMILES"] = pybel.readstring("inchi", id_info["InChI"]).write("smi")

    rdkit_mol = MolFromSmiles(id_info["SMILES"])

    id_info["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(rdkit_mol)

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
            kegg_dict = KEGGParser().parse(KEGG().get(sources["KEGG Compound"]))
            pathways["KEGG"] = kegg_dict["PATHWAY"].keys()
        except KeyError:
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

if __name__ == "__main__":
    for inchikey in combined.keys():
        inchikey = "CZMRCDWAGMRECN-UGDNZRGBSA-N"
        #id_info, rdkit_mol = identification_info(inchikey)
        #p_properties = physiochemical(rdkit_mol)
        #t_properties = taxonomic_properties(inchikey)
        #sources = generate_sources(inchikey)
        #pathway_info = generate_pathways(inchikey, sources)
        break

