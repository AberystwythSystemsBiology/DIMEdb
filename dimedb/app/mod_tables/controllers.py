from flask import Blueprint, request, render_template, g, session, redirect, url_for, jsonify
from app import app, db

from app.mod_tables.models import MetaboliteTables, Metabolite
from app.mod_auth.models import User

tables = Blueprint("tables", __name__)

@app.route("/tables")
def tables_home():
    public_tables = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.creation_date,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.public == True).all()

    return render_template("tables/index.html", public_tables=public_tables)

@app.route("/tables/DdbT<id>")
def view_table(id):
    table_info = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.description,
        MetaboliteTables.doi,
        MetaboliteTables.creation_date,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.id == id).first()

    return render_template("tables/view.html", table_info = table_info)


@app.route("/tables/DdbT<id>/get_metabolites")
def get_metabolites(id):
    metabolites = Metabolite.query.filter(Metabolite.table_id == id).all()
    json_dict = {"data" : []}
    for tm in metabolites:
        metabolite_dict = {"Comments": None, "InChIKey": None,
                           "Name": None, "Molecular Formula": None, "Table ID" : None}

        metabolite = app.data.driver.db["metabolites"].find_one({'_id': tm.inchikey})

        metabolite_dict["Comments"] = tm.comment
        metabolite_dict["InChIKey"] = tm.inchikey
        metabolite_dict["Name"] = metabolite["Identification Information"]["Name"]
        metabolite_dict["Molecular Formula"] = metabolite["Identification Information"]["Molecular Formula"]
        metabolite_dict["Table ID"] = tm.id
        json_dict["data"].append(metabolite_dict)

    return jsonify(json_dict)
