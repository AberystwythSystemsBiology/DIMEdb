from flask import Blueprint, flash, render_template, g, redirect, url_for, jsonify
from app import app, db

from app.mod_tables.models import MetaboliteTables, Metabolite
from app.mod_auth.models import User

from app.mod_tables.forms import CreateTableForm


tables = Blueprint("tables", __name__)

@app.route("/tables")
def public_tables():
    public_tables = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.creation_date,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.public == True).all()
    return render_template("tables/index.html", public_tables=public_tables)

@app.route("/tables/new", methods=["GET", "POST"])
def create_table():
    form = CreateTableForm()
    if form.validate_on_submit():
        metabolite_table = MetaboliteTables(
            title=form.title.data,
            description=form.description.data,
            species=form.species.data,
            doi=form.doi.data,
            owner_id=g.user.id
        )
        db.session.add(metabolite_table)
        db.session.commit()
        db.session.flush()
    return render_template("tables/new.html", form=form)

@app.route("/tables/DdbT<id>")
def view_table(id):
    table_info = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.description,
        MetaboliteTables.doi,
        MetaboliteTables.creation_date,
        MetaboliteTables.species,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.id == id).first()

    return render_template("tables/view.html", table_info = table_info)



@app.route("/tables/DdbT<id>/get_metabolites")
def get_metabolites(id):
    metabolites = Metabolite.query.filter(Metabolite.table_id == id).all()
    json_dict = {"data" : []}
    for tm in metabolites:
        metabolite_dict = {"InChIKey": None, "Name": None, "Molecular Formula": None, "Table ID" : None}

        metabolite = app.data.driver.db["metabolites"].find_one({'_id': tm.inchikey})

        metabolite_dict["InChIKey"] = tm.inchikey
        metabolite_dict["Name"] = metabolite["Identification Information"]["Name"]
        metabolite_dict["Molecular Formula"] = metabolite["Identification Information"]["Molecular Formula"]
        metabolite_dict["Table ID"] = tm.id
        json_dict["data"].append(metabolite_dict)

    return jsonify(json_dict)
