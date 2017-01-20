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
    peak_value = float(max(signals.values())[0])
    for mz, rel_int in signals.items():
        mz = round(mz, 6)
        if mz not in n_d:
            n_d[mz] = float(sum(rel_int)) * 100 / peak_value
        else:
            n_d[mz] = n_d[mz] + float(sum(rel_int)) * 100 / peak_value

    return sorted(n_d.items(), key=lambda x: x[0])

def function_name(element):
    ratios, weights = isotopes(element)
    nd = calculate_nom_distribution(ratios, weights)
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

def rule_dict_based(s, rule_dict):
    if "multiply" in rule_dict.keys():
        for element in s:
            s[element] = s[element] * rule_dict["multiply"]
    for element in rule_dict["remove"]:
        s[element] = s[element] - rule_dict["remove"][element]
    for element in rule_dict["add"]:
        try:
            s[element] = s[element] + rule_dict["add"][element]
        except KeyError:
            s[element] = rule_dict["add"][element]
    return s

def adduct_calculator(formula, rule_dict):
    structure_dict = split(formula)
    altered_sd = rule_dict_based(structure_dict, rule_dict)
    calculated_element = element_calculator(altered_sd)
    d, accurate_mass = function_name(calculated_element)
    return accurate_mass, d


def rules(formula, mol):
    nacc = rdMolDescriptors.CalcNumHBA(mol)
    ndon = rdMolDescriptors.CalcNumHBD(mol)
    noh = sum([Fragments.fr_Al_OH(mol), Fragments.fr_Ar_OH(mol)])
    nnhh = Fragments.fr_NH2(mol)
    ncooh = Fragments.fr_COO(mol)
    nch = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])

    structure_dict = split(formula)
    adducts = collections.defaultdict(list)
    # M (Neutral)
    nominal_element = element_calculator(structure_dict)
    accurate_mass = sum([x["molecular_weight"] for x in nominal_element.values()])
    d, am = function_name(nominal_element)
    adducts["neutral"].append(gen_rule_dict("M", am, d))

    adducts["negative"] = []
    adducts["positive"] = []

    # Negative Adducts
    if ndon > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"remove": {"H": 1}, "add": {}})
        adducts["negative"].append(gen_rule_dict("M-H", am, d))

    if ndon > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"remove" : {"H" : 1}, "add" : {}, "multiply" : 2})
        adducts["negative"].append(gen_rule_dict("2M-H", am, d))

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"remove" : {}, "add" : {"Br" : 1}})
        adducts["negative"].append(gen_rule_dict("M+Br", am, d))

    if ndon > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"remove" : {"H" : 1}, "add" : {}, "multiply" : 3})
        adducts["negative"].append(gen_rule_dict("3M-H", am, d))

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"remove": {}, "add": {"Cl" : 1}})
        adducts["negative"].append(gen_rule_dict("M+Cl", am, d))

    # Positive Adducts

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"H": 1}, "remove": {}})
        adducts["positive"].append(gen_rule_dict("M+H", am, d))

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"Na": 1}, "remove": {}})
        adducts["positive"].append(gen_rule_dict("M+Na", am, d))

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"K": 1}, "remove": {}})
        adducts["positive"].append(gen_rule_dict("M+K", am, d))

    if ndon > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"Na": 2}, "remove": {"H":1}})
        adducts["positive"].append(gen_rule_dict("M+2Na-H", am, d))

    if ndon > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"K": 2}, "remove": {"H":1}})
        adducts["positive"].append(gen_rule_dict("M+2K-H", am, d))

    if nacc > 0 and nch == 0:
        am, d = adduct_calculator(formula, {"add": {"H": 1}, "remove": {}, "multiply" : 2})
        adducts["positive"].append(gen_rule_dict("2M+H", am, d))

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
        if len(elem_c[idx]) == 0:
            elem_c[idx] = [0]
        try:
            structure_dict[e] = int(elem_c[idx][0])
        except IndexError:
            print "We have serious problems"
    return structure_dict

def process_entity(entity, id):
    final_d = {}
    try:
        mol = Chem.MolFromSmiles(entity["smiles"])
        formula = rdMolDescriptors.CalcMolFormula(mol)
        accurate_mass, adducts = rules(formula, mol)
        structure_dict = split(formula)

        final_d = {
            "sources" : entity["sources"],
            "synonyms" : entity["synonyms"],
            "name": entity["name"],
            "molecular_formula": formula,
            "smiles": entity["smiles"],
            "origins": entity["origins"],
            "accurate_mass": accurate_mass,
            "pathways" : entity["pathways"],
            "adducts": adducts,
            "num_atoms" : sum([structure_dict[x] for x in structure_dict])-1
        }

    except Exception, err:
        pass
    return final_d

def generate_db_file(output):
    db = Parallel(n_jobs=4)(delayed(process_entity)(output[id], id) for id in output)
    '''
    db = []

    for id in tqdm.tqdm(output):
        db.append(process_entity(output[id], id))
        break
    '''
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