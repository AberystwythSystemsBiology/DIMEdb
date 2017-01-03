import os, urllib2, json, collections, re, numpy as np, operator, itertools, tqdm

from bson.json_util import dumps

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Fragments

url="https://gist.githubusercontent.com/KeironO/a2ce6d7fb7e7e10f616a51f511cb27b4/raw/71b0ad0abf92ef101a4833dd9ee0837609939f9d/gistfile1.txt"
periodic_table = json.loads(urllib2.urlopen(url).read())

def load_output(hmdb_fp="./output/hmdb.json"):
    output = collections.defaultdict(dict)
    ds = [hmdb_fp]
    for s in ds:
        if os.path.isfile(s) != True:
            print "File doesn't exist:", s
            # TODO: Interface with modules to create file in this instance.
        else:
            with open(s, "r") as file:
                output.update(json.load(file))
    return output

def formula_splitter(chem_form):
    atom_dict = collections.defaultdict()
    t_atoms = [re.findall('[0-9]+', x) for x in re.findall('[A-Z][a-z]*[0-9]*', chem_form)]
    if int(sum(t_atoms, [])[0]) < 5:
        return False, False
    else:
        for element in re.findall('[A-Z][a-z]*[0-9]*', chem_form):
            atom = re.findall('[A-Z][a-z]*', element)[0]
            atom_count = re.findall('[0-9]+', element)
            if len(atom_count) < 1:
                atom_count = 1
            else:
                atom_count = int(atom_count[0])
            atom_weight = float(
                np.matrix(periodic_table[atom]["isotopic_weight"]) * np.transpose(np.matrix(periodic_table[atom]["isotopic_ratio"])))
            atom_dict[atom] = {
                "Atom Count": atom_count,
                "Atomic Charge": periodic_table[atom]["atomic_charge"],
                "Element Weight": atom_weight,
                "Total Weight": atom_count * atom_weight,
                "Isotopic Ratios": periodic_table[atom]["isotopic_ratio"],
                "Isotopic Weights": periodic_table[atom]["isotopic_weight"]
            }
        return atom_dict

def cartesian(ratios, weights, f_ratios, f_weights, count=1, threshold=0.25):
    n_ratios = []
    n_weights = []
    normalised_ratio = [n / max(f_ratios) for n in f_ratios]

    for i in enumerate(ratios[count]):
        r = ratios[count][i[0]]
        w = weights[count][i[0]]
        for j in enumerate(normalised_ratio):
            current_n_ratio = normalised_ratio[j[0]]*100
            current_s_weight = f_weights[j[0]]
            t_w = current_n_ratio * r
            if t_w > threshold:
                n_ratios+=[current_n_ratio* r]
                n_weights+=[current_s_weight+w]
    count = count + 1
    if count < len(ratios) and len(n_ratios) < 1000:
        n_ratios, n_weights = cartesian(ratios, weights, n_ratios, n_weights, count)
    return n_ratios, n_weights

def isotopes(atoms):
    ratios = []
    weights = []
    for x in atoms:
        for i in range(atoms[x]["Atom Count"]):
            ratios.append(atoms[x]["Isotopic Ratios"])
            weights.append(atoms[x]["Isotopic Weights"])
    return cartesian(ratios, weights, ratios[0], weights[0])

# nom nom nom
def calculate_nom_distribution(weights, ratios):
    paired_w_r = [(weights[index], r) for index, r in enumerate(ratios) if r > 1e-6]
    signals = dict((key, tuple(v for (k, v) in pairs))
                   for (key, pairs) in itertools.groupby(sorted(paired_w_r), operator.itemgetter(0)))
    n_d = {}
    lv = float(max(signals.values())[0])
    for mz, rel_int in signals.items():
        n_d[mz] = float(sum(rel_int)) * 100 / lv
    return sorted(n_d.items(), key=lambda x: x[0])

def get_anion(nacc, ndon, noh, nnhh, ncooh, nch, nominal_distribution):
    anion = {}
    if ndon > 0 and nch == 0:
        anion["[M-H]1-"] = [(x[0] - 1, x[1]) for x in nominal_distribution]
    return anion

def get_canion(nacc, ndon, noh, nnhh, ncooh, nch, nominal_distribution):
    canion = {}
    return canion

def adduct(mol, nominal_distribution):
    nacc = rdMolDescriptors.CalcNumHBA(mol)
    ndon = rdMolDescriptors.CalcNumHBD(mol)
    noh = sum([Fragments.fr_Al_OH(mol), Fragments.fr_Ar_OH(mol)])
    nnhh = Fragments.fr_NH2(mol)
    ncooh = Fragments.fr_COO(mol)
    nch = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    adduct_dict = {}
    adduct_dict["Nominal"] = nominal_distribution
    adduct_dict["Anion"] = get_anion(nacc, ndon, noh, nnhh, ncooh, nch, nominal_distribution)
    adduct_dict["Canion"] = get_canion(nacc, ndon, noh, nnhh, ncooh, nch, nominal_distribution)
    return adduct_dict

def process_entity(entity):
    # TODO: See where these exceptions are...
    try:
        mol = Chem.MolFromSmiles(entity["SMILES"])
        formula = rdMolDescriptors.CalcMolFormula(mol)
        atom_dict = formula_splitter(formula)
        atomic_charge = sum([atom_dict[x]["Atomic Charge"] for x in atom_dict.keys()])
        num_of_atoms = sum([atom_dict[x]["Atom Count"] for x in atom_dict.keys()])
        accurate_mass = sum([atom_dict[x]["Total Weight"] for x in atom_dict.keys()])
        ratios, weights = isotopes(atom_dict)
        nominal_distribution = calculate_nom_distribution(weights, ratios)
        isotopic_weight = nominal_distribution[0][0]
        ad = adduct(mol, nominal_distribution)
        final_d = {
            "name": entity["Name"],
            "synonyms": entity["Synonyms"],
            "accurate_mass": accurate_mass,
            "isotopic_weight": isotopic_weight,
            "molecular_formula": formula,
            "smiles": entity["SMILES"],
            "num_of_atoms": num_of_atoms,
            "atomic_charge": atomic_charge,
            "adducts": ad,
            "origins" : entity["Origins"]
        }
        return final_d
    except Exception, err:
        return None


def generate_db_file(output):
    db = []
    for id in tqdm.tqdm(output):
        pe = process_entity(output[id])
        if pe != None:
            db.append(pe)
    return db

def save_db_file(db, fp="./output/mb-db.json"):
    mongodb_file = json.loads(dumps(db))
    with open(fp, "w") as output:
        json.dump(mongodb_file, output)

if __name__ == "__main__":
    output = load_output(hmdb_fp="./output/hmdb_small.json")
    db = generate_db_file(output)
    save_db_file(db)