import json, csv, zipfile, os

directory = "/home/keo7/.data/dimedb/"


def load_file(fp):
    with open(fp, "rb") as dimedb_file:
        return json.load(dimedb_file)

def generate_id_tsv(dimedb):
    names = []
    names.append(["InChIKey", "Property", "Entry"])
    for metabolite in dimedb:
        id = metabolite["_id"]
        id_info = metabolite["Identification Information"]
        names.append([id, "Name", id_info["Name"]])

        if id_info["IUPAC Name"] != None:
            names.append([id, "IUPAC Name", id_info["IUPAC Name"]])

        if id_info["Systematic Name"] != None:
            names.append([id, "Systematic Name", id_info["Systematic Name"]])

        for s in id_info["Synonyms"]:
            names.append([id, "Synonym", s])

        names.append([id, "SMILES", id_info["SMILES"]])
        names.append([id, "Molecular Formula", id_info["Molecular Formula"]])
        names.append([id, "InChI", id_info["InChI"]])

    with open(directory+"downloads/dimedb_id_info.tsv", "wb") as names_file:
        writer = csv.writer(names_file, delimiter="\t")
        for line in names:
            writer.writerow(line)

def generate_physiochemical_properties(dimedb):
    pc_p_l = []
    pc_p_l.append(["InChIKey", "Property", "Entry"])
    for metabolite in dimedb:
        id = metabolite["_id"]
        pc_p = metabolite["Physiochemical Properties"]
        for key in pc_p.keys():
            pc_p_l.append([id, key, pc_p[key]])

    with open(directory+"downloads/dimedb_pc_info.tsv", "wb") as names_file:
        writer = csv.writer(names_file, delimiter="\t")
        for line in pc_p_l:
            writer.writerow(line)

def generate_structures():
    zipf = zipfile.ZipFile(directory+"downloads/structures.zip", "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(directory+"structures/"):
        for file in files:
            zipf.write(os.path.join(root, file), arcname=file)
    zipf.close()

def generate_pathways(dimedb):
    pw_l = []
    pw_l.append(["InChIKey", "Pathway", "ID"])
    for metabolite in dimedb:
        id = metabolite["_id"]
        for kegg_id in metabolite["Pathways"]["KEGG"]:
            pw_l.append([id, "KEGG", kegg_id])
        for smpdb_id in metabolite["Pathways"]["SMPDB"]:
            pw_l.append([id, "SMPDB", smpdb_id])

    with open(directory+"downloads/dimedb_pathways.tsv", "wb") as names_file:
        writer = csv.writer(names_file, delimiter="\t")
        for line in pw_l:
            writer.writerow(line)

def generate_sources(dimedb):
    sources_l = []
    sources_l.append(["InChIKey", "Source", "ID"])
    for metabolite in dimedb:
        id = metabolite["_id"]
        for source, sid in metabolite["External Sources"].items():
            if sid != None:
                sources_l.append([id, source, sid])

    with open(directory+"downloads/dimedb_sources.tsv", "wb") as names_file:
        writer = csv.writer(names_file, delimiter="\t")
        for line in sources_l:
            writer.writerow(line)

if __name__ == "__main__":
    dimedb = load_file(directory+"dimedb.json")
    generate_id_tsv(dimedb)
    generate_physiochemical_properties(dimedb)
    generate_structures()
    generate_pathways(dimedb)
    generate_sources(dimedb)