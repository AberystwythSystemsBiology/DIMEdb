import os, urllib2, json, collections, re, numpy as np, operator, itertools, tqdm

from joblib import Parallel, delayed

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

def isotopes(elements):
    ratios = []
    weights = []
    for x in elements.keys():
        for i in range(elements[x]["element_count"]):
            ratios.append(elements[x]["isotopic_ratios"])
            weights.append(elements[x]["isotopic_weights"])
    return cartesian(ratios, weights, ratios[0], weights[0])

def calculate_nom_distribution(ratios, weights):
    paired_w_r = [(weights[index], r) for index, r in enumerate(ratios) if r > 1e-6]
    signals = dict((key, tuple(v for (k, v) in pairs))
                   for (key, pairs) in itertools.groupby(sorted(paired_w_r), operator.itemgetter(0)))
    n_d = {}
    lv = float(max(signals.values())[0])
    for mz, rel_int in signals.items():
        n_d[mz] = float(sum(rel_int)) * 100 / lv
    return sorted(n_d.items(), key=lambda x: x[0])

def function_name(element):
    raitos, weights = isotopes(element)
    nd = calculate_nom_distribution(raitos, weights)
    accurate_mass = nd[0][0]
    return nd, accurate_mass

def calculate_element_weight(element):
    iso_weight = periodic_table[element]["isotopic_weight"]
    ratio = periodic_table[element]["isotopic_ratio"]
    return float(np.matrix(ratio) * np.transpose(np.matrix(iso_weight)))

def element_calculator(structure_dict):
    elements = {}
    for e in structure_dict.keys():
        element_weight = calculate_element_weight(e)
        elements[e] = {
            "element_count" : structure_dict[e],
            "atomic_charge" : periodic_table[e]["atomic_charge"],
            "molecular_weight" : int(structure_dict[e]) * element_weight,
            "isotopic_ratios" : periodic_table[e]["isotopic_ratio"],
            "isotopic_weights" :  periodic_table[e]["isotopic_weight"]
        }
    return elements

def gen_rule_dict(t, am, d):
    return {
            "type" : t,
            "accurate_mass" : am,
            "isotopic_distribution" : d
    }

def rule_dict_based(structure_dict, rule_dict):
    for element in rule_dict["remove"]:
        structure_dict[element] = structure_dict[element] - rule_dict["remove"][element]
    for element in rule_dict["add"]:
        try:
            structure_dict[element] = structure_dict[element] + rule_dict["add"][element]
        except KeyError:
            structure_dict[element] = rule_dict["add"][element]
    return structure_dict

def rules(structure_dict, mol):
    nacc = rdMolDescriptors.CalcNumHBA(mol)
    ndon = rdMolDescriptors.CalcNumHBD(mol)
    noh = sum([Fragments.fr_Al_OH(mol), Fragments.fr_Ar_OH(mol)])
    nnhh = Fragments.fr_NH2(mol)
    ncooh = Fragments.fr_COO(mol)
    nch = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])

    adducts = collections.defaultdict(list)
    # M (Neutral)
    nominal_element = element_calculator(structure_dict)
    accurate_mass = sum([x["molecular_weight"] for x in nominal_element.values()])
    d, am = function_name(nominal_element)
    adducts["neutral"].append(gen_rule_dict("M", am, d))

    adducts["negative"] = []
    adducts["positive"] = []

    if ndon > 1 and nacc == 0:
        sd = rule_dict_based(structure_dict, {"remove" : {"H" : 1}, "add" : {}})
        d, am = function_name(sd)
        adducts["negative"].append(gen_rule_dict("M-H", am, d))
    if nacc > 0 and nch == 0:
        sd = rule_dict_based(structure_dict, {"remove" : {}, "add" : {"Na" : 1}})
        element = element_calculator(sd)
        d, am = function_name(element)
        adducts["negative"].append(gen_rule_dict("M+Na", am, d))

    final_adducts = {}
    for ion in adducts.keys():
        data = {
            "count" : len(adducts[ion]),
            "peaks" : [x for x in adducts[ion]]
        }
        final_adducts[ion] = data
    return accurate_mass, final_adducts

def split(formula):
    structure_dict = {}
    elem = re.findall('[A-Z][a-z]*', formula)
    elem_c = [re.findall('[0-9]+', x) for x in re.findall('[A-Z][a-z]*[0-9]*', formula)]
    for idx, e in enumerate(elem):
        if elem_c[idx] == 0:
            elem_c[idx] = 0
        structure_dict[e] = int(elem_c[idx][0])
    return structure_dict

def process_entity(entity):
    # TODO: See where these exceptions are...
    try:
        final_d = {}
        mol = Chem.MolFromSmiles(entity["smiles"])
        formula = rdMolDescriptors.CalcMolFormula(mol)
        structure_dict = split(formula)
        accurate_mass, adducts = rules(structure_dict, mol)

        final_d = {
            "name" : entity["name"],
            "molecular_formula" : formula,
            "smiles" : entity["smiles"],
            "origins" : entity["origins"],
            "accurate_mass" : accurate_mass,
            "adducts" : adducts,
        }
        return final_d
    except Exception, err:
        return None


def generate_db_file(output):
    db = Parallel(n_jobs=4)(delayed(process_entity)(output[id]) for id in output)
    db = [x for x in db if x != None]
    return db

def save_db_file(db, fp="./output/mb-db.json"):
    mongodb_file = json.loads(dumps(db))
    with open(fp, "w") as output:
        json.dump(mongodb_file, output, indent=4)

if __name__ == "__main__":
    output = load_output()
    db = generate_db_file(output)
    save_db_file(db)