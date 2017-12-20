import collections
import os
import json
import dicttoxml
from flask import Blueprint, send_from_directory, render_template
from dimedb import app, abort

view = Blueprint("view", __name__, url_prefix="/view")

source_dict = {
    "HMDB Accession": "http://www.hmdb.ca/metabolites/value",
    "Massbank": "http://www.massbank.jp/jsp/Dispatcher.jsp?type=disp&id=value&site=10",
    "KEGG Compound": "http://www.genome.jp/dbget-bin/www_bget?value",
    "KEGG Drug": "http://www.genome.jp/dbget-bin/www_bget?dr:value",
    "BioCyc": "https://biocyc.org/compound?orgid=META&id=value",
    "CAS": None,
    "GOLM": "http://gmd.mpimp-golm.mpg.de/Analytes/value",
    "Spektraris": "http://langelabtools.wsu.edu/nmr/record/value",
    "PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/value",
    "Wikidata": "https://www.wikidata.org/wiki/value",
    "Respect": "http://spectra.psc.riken.jp/menta.cgi/respect/datail/datail?accession=value",
    "Chemspider": "http://www.chemspider.com/Chemical-Structure.value.html",
    "NMR Shift DB": None,
    "ChEBI": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=value"
}


@view.route("/<inchikey>")
def view_metabolite(inchikey):
    mtbl_driver = app.data.driver.db["metabolites"]

    mtbl = mtbl_driver.find({"_id": inchikey}).limit(1)
    if mtbl.count() != 1:
        abort(404)
    else:
        metabolite = mtbl[0]
        id_info = metabolite["Identification Information"]
        phy_prop = metabolite["Physicochemical Properties"]
        tax_prop = metabolite["Taxonomic Properties"]
        ext_sour = metabolite["External Sources"]
        pathways = metabolite["Pathways"]

        return render_template("view/view.html", id=inchikey, id_info=id_info,
                               phy_prop=phy_prop, tax_prop=tax_prop,
                               ext_sour=ext_sour, pathways=pathways,
                               source_dict=source_dict)


@view.route("/<inchikey>/jcamp/")
def to_jcamp(inchikey):
    mtbl_driver = app.data.driver.db["metabolites"]
    mtbl = mtbl_driver.find({"_id": inchikey}).limit(1)
    if mtbl.count() != 1:
        response = app.response_class(
            status=404,
            mimetype="text/plain"
        )
    else:
        metabolite = mtbl[0]
        jcamp = collections.OrderedDict([
            ["TITLE", None],
            ["JCAMP", 5.01],
            ["DATA TYPE", "MASS SPECTRUM"],
            ["DATA CLASS", "XYDATA"],
            ["ORIGIN", "DIMEdb neutral DIMS profile."],
            ["OWNER", "Direct Infusion Metabolite Database (DIMEdb) https://dimedb.ibers.aber.ac.uk"],
            ["VERSION", "NA"],
            ["$DIMEDB COMPOUND URL", None],
            ["$SMILES", None],
            ["CAS REGISTRY", None],
            ["MOLFORM", None],
            ["MW", None],
            ["XUNITS", "M/Z"],
            ["YUNITS", "RELATIVE ABUNDANCE"],
            ["XFACTOR", 1],
            ["YFACTOR", 1],
            ["FIRSTX", None],
            ["FIRSTY", None],
            ["MAXX", None],
            ["MINX", None],
            ["MAXY", None],
            ["MINY", None],
            ["NPOINTS", None],
            ["PEAK TABLE", "(XY..XY)\r\n"],
            ["END", "\r\n"]])

        isotopic_distribution = metabolite["Adducts"]
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
                    vals = {"mass": round(mass, 2),
                            "int": round(Y[index], 2)}
                    jcamp["PEAK TABLE"] += (
                        "({mass},{int}) ".format(**vals))
                    if limiter == 3 and index != len(X) and len(X) != 3:
                        limiter = 0
                        jcamp["PEAK TABLE"] += "\n"
                    else:
                        jcamp["PEAK TABLE"] += "  "
                    limiter += 1
                break

        if metabolite["External Sources"]["CAS"] is not None:
            jcamp["CAS REGISTRY"] = metabolite["External Sources"]["CAS"]
        else:
            jcamp["CAS REGISTRY"] = "NA"

        jcamp["TITLE"] = metabolite["Identification Information"]["Name"].upper()
        jcamp["$DIMEDB COMPOUND URL"] = "https://dimedb.ibers.aber.ac.uk/view/" + metabolite["_id"]
        jcamp["MOLFORM"] = metabolite["Identification Information"]["Molecular Formula"].upper()

        jcamp["$SMILES"] = metabolite["Identification Information"]["SMILES"].upper()

        jcamp_pt = ""

        for key, value in jcamp.items():
            if type(value) in [str, unicode]:
                value = value.encode('utf-8')
            if key == "END":
                jcamp_pt += "\n"
            jcamp_pt += "##" + key + "=" + str(value)
            if key != "PEAK TABLE":
                jcamp_pt += "\n"

        response = app.response_class(
            status=200,
            response=jcamp_pt,
            mimetype="text/plain",
            headers={"Content-disposition":
                    "attachment; filename="+inchikey+".JCAMP"}
        )
    return response


@view.route("/<inchikey>/json")
def to_json(inchikey):
    mtbl_driver = app.data.driver.db["metabolites"]

    mtbl = mtbl_driver.find({"_id": inchikey}).limit(1)
    if mtbl.count() != 1:
        response = app.response_class(
            status=404,
            mimetype="application/json"
        )
    else:
        response = app.response_class(
            response=json.dumps(mtbl[0]),
            status=200,
            mimetype="application/json",
            headers={"Content-disposition":
                    "attachment; filename="+inchikey+".json"}
        )
    return response


@view.route("/<inchikey>/xml")
def to_xml(inchikey):
    mtbl_driver = app.data.driver.db["metabolites"]

    mtbl = mtbl_driver.find({"_id": inchikey}).limit(1)
    if mtbl.count() != 1:
        response = app.response_class(
            status=404,
            mimetype="application/xml"
        )
    else:
        response = app.response_class(
            response=dicttoxml.dicttoxml(mtbl[0], custom_root="metabolite"),
            status=200,
            mimetype="application/xml",
            headers={"Content-disposition":
                    "attachment; filename="+inchikey+".xml"}
        )
    return response


@view.route("/<inchikey>/structure/")
def get_structure(inchikey):
    d = os.path.expanduser("~/Data/dimedb/output/structures/")
    return send_from_directory(d, inchikey + ".svg")
