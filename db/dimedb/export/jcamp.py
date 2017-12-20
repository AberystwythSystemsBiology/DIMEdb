VERSION = "1.0NOV2017"

import collections, json, os

output_dir = os.path.join(os.path.expanduser("~"), "Data/dimedb/output/")
final_dir = os.path.join(output_dir, "dimedb_jsons/")
jcamp_dir = os.path.join(output_dir, "download/jcamp/")


if os.path.isdir(jcamp_dir) != True:
    os.makedirs(jcamp_dir)

def convert(data):

    jcamp = collections.OrderedDict([
        ["TITLE", None],
        ["JCAMP", 5.01],
        ["DATA TYPE", "MASS SPECTRUM"],
        ["DATA CLASS", "XYDATA"],
        ["ORIGIN", "DIMEdb neutral DIMS profile."],
        ["OWNER", "Direct Infusion Metabolite Database (DIMEdb) https://dimedb.ibers.aber.ac.uk"],
        ["VERSION", VERSION],
        ["$DIMEDB COMPOUND URL", None],
        ["$SMILES", None],
        ["CAS REGISTRY", None],
        ["MOLFORM", None],
        ["MW", None],
        ["XUNITS", "M/Z"],
        ["YUNITS", "RELATIVE ABUNDANCE"],
        ["XFACTOR", 1],
        ["YFACTOR" , 1],
        ["FIRSTX", None],
        ["FIRSTY", None],
        ["MAXX", None],
        ["MINX", None],
        ["MAXY", None],
        ["MINY", None],
        ["NPOINTS", None],
        ["PEAK TABLE", "(XY..XY)\r\n"],
        ["END", "\r\n"]])


    isotopic_distribution = data["Adducts"]
    for i in isotopic_distribution:
        if i["Polarity"] == "Neutral":
            jcamp["MW"] = i["Accurate Mass"]
            X = [x[0] for x in i["Isotopic Distribution"]]
            Y = [x[1] for x in i["Isotopic Distribution"]]

            jcamp["FIRSTX"] = X[0]
            jcamp["FIRSTY"] = Y[0]

            jcamp["MAXX"] = max(X)
            jcamp["MINX"] = min(X)

            jcamp["MAXY"] = max(Y)
            jcamp["MINY"] = min(Y)

            jcamp["NPOINTS"] = len(X)

            limiter = 0
            for index, mass in enumerate(X):
                vals = {"mass": round(mass, 2), "int": round(Y[index], 2)}
                jcamp["PEAK TABLE"] += ("({mass},{int}) ".format(**vals))
                if limiter == 3 and index != len(X) and len(X) != 3:
                    limiter = 0
                    jcamp["PEAK TABLE"] += "\n"
                else:
                    jcamp["PEAK TABLE"] += "  "
                limiter += 1
            break


    if data["External Sources"]["CAS"] != None:
        jcamp["CAS REGISTRY"] = data["External Sources"]["CAS"]
    else:
        jcamp["CAS REGISTRY"] = "NA"


    jcamp["TITLE"] = data["Identification Information"]["Name"].upper()
    jcamp["$DIMEDB COMPOUND URL"] = "https://dimedb.ibers.aber.ac.uk/view/"+data["_id"]
    jcamp["MOLFORM"] = data["Identification Information"]["Molecular Formula"].upper()

    jcamp["$SMILES"] = data["Identification Information"]["SMILES"].upper()

    jcamp_pt = ""

    for key, value in jcamp.items():
        if type(value) in [str, unicode]:
            value = value.encode('utf-8')
        if key == "END":
            jcamp_pt += "\n"
        jcamp_pt += "##"+key+"="+str(value)
        if key != "PEAK TABLE":
            jcamp_pt += "\n"

    return jcamp_pt

if __name__ == "__main__":
    dimedb = []

    jcamp_str = ""

    for fn in os.listdir(final_dir):
        if fn.endswith(".json"):
            with open(final_dir+fn, "rb") as dimedb_file:
                jfile = json.load(dimedb_file)
                dimedb.append(jfile)

    for i in dimedb:
        for d in i:
            jcamp_str += convert(d)

    with open(os.path.join(jcamp_dir + "dimedb_jcamp.txt"), "wb") as outfile:
        outfile.write(jcamp_str)
        outfile.close()
