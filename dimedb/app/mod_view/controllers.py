import os
import json
from flask import Blueprint, send_from_directory, render_template
from dimedb import app, abort

view = Blueprint("view", __name__, url_prefix="/view")

source_dict = {
    "HMDB Accession" : "http://www.hmdb.ca/metabolites/value",
    "Massbank" : "http://www.massbank.jp/jsp/Dispatcher.jsp?type=disp&id=value&site=10",
    "KEGG Compound": "http://www.genome.jp/dbget-bin/www_bget?value",
    "KEGG Drug": "http://www.genome.jp/dbget-bin/www_bget?dr:value",
    "BioCyc": "https://biocyc.org/compound?orgid=META&id=value",
    "CAS": None,
    "GOLM": "http://gmd.mpimp-golm.mpg.de/Analytes/value",
    "Spektraris": "http://langelabtools.wsu.edu/nmr/record/value",
    "PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/value",
    "Wikidata": "https://www.wikidata.org/wiki/value",
    "Respect": "http://spectra.psc.riken.jp/menta.cgi/respect/datail/datail?accession=value",
    "Chemspider": "www.chemspider.com/Chemical-Structure.value.html",
    "NMR Shift DB": None,
    "ChEBI": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=value"
}

@view.route("/<inchikey>")
def view_metabolite(inchikey):
    metabolite = app.data.driver.db["metabolites"].find({"_id" : inchikey}).limit(1)
    if metabolite.count() != 1:
        abort(404)
    else:
        metabolite = metabolite[0]
        id_info = metabolite["Identification Information"]
        phy_prop = metabolite["Physicochemical Properties"]
        tax_prop = metabolite["Taxonomic Properties"]
        ext_sour = metabolite["External Sources"]
        pathways = metabolite["Pathways"]
        return render_template("view/view.html", id=inchikey, id_info=id_info,
                                phy_prop=phy_prop, tax_prop=tax_prop,
                                ext_sour=ext_sour, pathways=pathways,
                                source_dict=source_dict)

@view.route("/<inchikey>/jcamp")
def to_jcamp(inchikey):
    metabolite = app.data.driver.db["metabolites"].find({"_id" : inchikey}).limit(1)
    if metabolite.count() != 1:
        abort(404)
    else:
        abort(503)

@view.route("/<inchikey>/json")
def to_json(inchikey):
    metabolite = app.data.driver.db["metabolites"].find({"_id" : inchikey}).limit(1)
    if metabolite.count() != 1:
        response = app.response_class(
            status = 404,
            mimetype = "application/json"
        )
    else:
        response = app.response_class(
            response = json.dumps(metabolite[0]),
            status = 200,
            mimetype ="application/json"
        )
    return response

@view.route("/<inchikey>/structure/")
def get_structure(inchikey):
    d = os.path.expanduser("~/Data/dimedb/output/structures/")
    return send_from_directory(d, inchikey + ".svg")
