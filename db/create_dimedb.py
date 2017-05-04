import json, collections
from rdkit.Chem import rdMolDescriptors, Fragments
from bioservices import KEGG, KEGGParser
from pyidick import Molecule
from bson.json_util import dumps
from joblib import Parallel, delayed

def load_files(hmdb_file="./output/hmdb.json"):
    with open(hmdb_file, "r") as hmdb:
        data = json.load(hmdb)
    return data

def get_adducts(py_mol):

    def calculate(type, mol, rule_dict=None, electrons=0, charge=0):

        iso_dist = mol.isotopic_distribution(rule_dict=rule_dict, electrons=electrons, charge=charge)
        accurate_mass = iso_dist[0][0]

        return {
            "type" : type,
            "accurate_mass" : accurate_mass,
            "isotopic_distribution": iso_dist

        }

    nacc = rdMolDescriptors.CalcNumHBA(py_mol._rdkmol)
    ndon = rdMolDescriptors.CalcNumHBD(py_mol._rdkmol)
    noh = sum([Fragments.fr_Al_OH(py_mol._rdkmol), Fragments.fr_Ar_OH(py_mol._rdkmol)])
    nnhh = Fragments.fr_NH2(py_mol._rdkmol)
    ncooh = Fragments.fr_COO(py_mol._rdkmol)
    nch = sum([atom.GetFormalCharge() for atom in py_mol._rdkmol.GetAtoms()])

    adducts = collections.defaultdict(list)
    adducts["neutral"].append(calculate("[M]", py_mol))

    # NEGATIVE
    if nch == -1:
        adducts["negative"].append(calculate("[M1-.]1-", py_mol, {"add": {}, "remove": {}}, charge=-1, electrons=-0))

    if ndon > 0 and nch == 0:
        adducts["negative"].append(
            calculate("[3M-H]1-", py_mol, {"add": {}, "remove": {"H": 1}, "multiply": 3}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[2M+Hac-H]1-", py_mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}, "multiply": 2}, charge=-1,
                      electrons=1))
        adducts["negative"].append(
            calculate("[2M+FA-H]1-", py_mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}, "multiply": 2}, charge=-1,
                      electrons=1))
        adducts["negative"].append(
            calculate("[2M-H]1-", py_mol, {"add": {}, "remove": {"H": 1}, "multiply": 2}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[M+TFA-H]1-", py_mol, {"add": {"C": 2, "O": 2, "F": 3}, "remove": {}}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[M+Hac-H]1-", py_mol, {"add": {"C": 2, "H": 3, "O": 2}, "remove": {}}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[M+FA-H]1-", py_mol, {"add": {"C": 1, "H": 1, "O": 2}, "remove": {}}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[M-H]1-", py_mol, {"add": {}, "remove": {"H": 1}}, charge=-1, electrons=1))
        if ndon > 1:
            adducts["negative"].append(
                calculate("[M-2H]2-", py_mol, {"add": {}, "remove": {"H": 2}}, charge=-2, electrons=2))
        if ndon > 2:
            adducts["negative"].append(
                calculate("[M-3H]3-", py_mol, {"add": {}, "remove": {"H": 3}}, charge=-3, electrons=3))
        if nacc > 0:
            adducts["negative"].append(
                calculate("[2M+Na-2H]1-", py_mol, {"add": {"Na": 1}, "remove": {"H": 2}, "multiply": 2}, charge=-1,
                          electrons=1))
            adducts["negative"].append(
                calculate("[M+K-2H]1-", py_mol, {"add": {"K": 1}, "remove": {"H": 2}}, charge=-1, electrons=1))

        if ndon > 1 and nacc > 0:
            adducts["negative"].append(
                calculate("[M+Na-2H]1-", py_mol, {"add": {"Na": 1}, "remove": {"H": 2}}, charge=-1, electrons=1))
    elif nacc > 0 and nch == 0:
        adducts["negative"].append(
            calculate("[M+Br]1-", py_mol, {"add": {"Br": 1}, "remove": {}}, charge=-1, electrons=1))
        adducts["negative"].append(
            calculate("[M+Cl]1-", py_mol, {"add": {"Cl": 1}, "remove": {}}, charge=-1, electrons=1))

    # POSITIVE

    if nch == 1:
        adducts["positive"].append(calculate("[M1+.]1+", py_mol, {"add": {}, "remove": {}}, charge=1, electrons=-0))

    if nch == 0:
        if nacc > 0:
            adducts["positive"].append(
                calculate("[2M+K]1+", py_mol, {"add": {"K": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[2M+Na]1+", py_mol, {"add": {"Na": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[2M+NH4]1+", py_mol, {"add": {"N": 4, "H": 4}, "remove": {}, "multiply": 2}, charge=1,
                          electrons=-1))
            adducts["positive"].append(
                calculate("[2M+H]1+", py_mol, {"add": {"H": 1}, "remove": {}, "multiply": 2}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+2K-H]1+", py_mol, {"add": {"K": 2}, "remove": {"H": 1}}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+2Na-H]1+", py_mol, {"add": {"Na": 2}, "remove": {"H": 1}}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+2K-H]1+", py_mol, {"add": {"K": 2}, "remove": {"H": 1}}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+K]1+", py_mol, {"add": {"K": 1}, "remove": {}}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+Na]1+", py_mol, {"add": {"Na": 1}, "remove": {}}, charge=1, electrons=-1))
            adducts["positive"].append(
                calculate("[M+H]1+", py_mol, {"add": {"H": 1}, "remove": {}}, charge=1, electrons=-1))
        if nacc > 1:
            adducts["positive"].append(
                calculate("[2M+3H2O+2H]2+", py_mol, {"add": {"H": 8, "O": 3}, "remove": {}, "multiply": 2}, charge=2,
                          electrons=-2))
            adducts["positive"].append(
                calculate("[2M+3ACN+2H]2+", py_mol, {"add": {"C": 6, "H": 11, "N": 3}, "remove": {}, "multiply": 2},
                          charge=2, electrons=-2))
            adducts["positive"].append(
                calculate("[M+2ACN+2H]2+", py_mol, {"add": {"C": 3, "H": 8, "N": 2}, "remove": {}}, charge=2,
                          electrons=-2))
            adducts["positive"].append(
                calculate("[M+2Na]2+", py_mol, {"add": {"Na": 2}, "remove": {}}, charge=2, electrons=-2))
            adducts["positive"].append(
                calculate("[M+2ACN+2H]2+", py_mol, {"add": {"C": 2, "H": 5, "N": 1}, "remove": {}}, charge=2,
                          electrons=-2))
            adducts["positive"].append(
                calculate("[M+H+K]2+", py_mol, {"add": {"K": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2))
            adducts["positive"].append(
                calculate("[M+H+Na]2+", py_mol, {"add": {"Na": 1, "H": 1}, "remove": {}}, charge=2, electrons=-2))
            adducts["positive"].append(
                calculate("[M+H+NH4]2+", py_mol, {"add": {"N": 5, "H": 5}, "remove": {}}, charge=2, electrons=-2))
            adducts["positive"].append(
                calculate("[M+2H]2+", py_mol, {"add": {"H": 2}, "remove": {}}, charge=2, electrons=-2))

        if nacc > 2:
            adducts["positive"].append(
                calculate("[M+3Na]3+", py_mol, {"add": {"Na": 3}, "remove": {}}, charge=3, electrons=-3))
            adducts["positive"].append(
                calculate("[M+H+2Na]3+", py_mol, {"add": {"H": 1, "Na": 2}, "remove": {}}, charge=3, electrons=-3))
            adducts["positive"].append(
                calculate("[M+2H+2a]3+", py_mol, {"add": {"H": 2, "Na": 1}, "remove": {}}, charge=3, electrons=-3))
            adducts["positive"].append(
                calculate("[M+3H]3+", py_mol, {"add": {"H": 3}, "remove": {}}, charge=3, electrons=-3))
    return adducts


def process_metabolite(metabolite):
    try:
        py_mol = Molecule(metabolite["smiles"])

        adducts = get_adducts(py_mol)

        try:
            kegg_dict = KEGGParser().parse(KEGG().get(metabolite["kegg_id"]))
            pathways = kegg_dict["PATHWAY"].keys()
        except Exception, err:
            pathways = None

        metabolite_tuple = [
            ["_id", metabolite["_id"]],
            ["name", metabolite["name"]],
            ["synonyms", metabolite["synonyms"]],
            ["molecular_formula", py_mol.molecular_formula],
            ["accurate_mass", py_mol.accurate_mass],
            ["num_atoms", py_mol.num_atoms],
            ["inchi", metabolite["inchi"]],
            ["smiles", metabolite["smiles"]],
            ["origins", metabolite["origins"]],
            ["biofluid_location", metabolite["biofluid_locations"]],
            ["tissue_locations", metabolite["tissue_locations"]],
            ["pathways", pathways],
            ["sources", metabolite["sources"]],
            ["adducts", adducts]
        ]


        return collections.OrderedDict(metabolite_tuple)
    except Exception, err:
        return None

if __name__ == "__main__":
    data = load_files()


    db = Parallel(n_jobs=300)(delayed(process_metabolite)(data[id]) for id in data)
    db = [x for x in db if x != None]

    fp = "./output/dimedb.json"

    mongodb_file = json.loads(dumps(db), object_pairs_hook=collections.OrderedDict)
    with open(fp, "wb") as output:
        json.dump(mongodb_file, output, indent=4)