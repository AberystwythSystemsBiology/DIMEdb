import json

from flask import Blueprint, flash, render_template, g, abort, redirect, url_for, jsonify, request
from app import app, db

from flask_login import login_required

from app.mod_tables.models import MetaboliteTables, Metabolite
from app.mod_auth.models import User

from app.mod_tables.forms import CreateTableForm, EditTableForm


tables = Blueprint("tables", __name__)

@app.route("/tables")
def public_tables():
    public_tables = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.creation_date,
        MetaboliteTables.doi,
        MetaboliteTables.species,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.public == True).all()
    return render_template("tables/index.html", public_tables=public_tables)

@app.route("/tables/mytables")
@login_required
def mytables():
    mytables = MetaboliteTables.query.filter(MetaboliteTables.owner_id== g.user.id).all()
    return render_template("tables/mytables.html", mytables=mytables)

@app.route("/tables/mytables/api/add_metabolite", methods=["POST"])
@login_required
def add_metabolite():
    payload = request.form
    table = MetaboliteTables.query.filter(MetaboliteTables.owner_id == g.user.id).first()
    if table != None:
        metabolite = Metabolite(
            table_id = payload["Table ID"],
            inchikey = payload["InChI Key"]
        )

        db.session.add(metabolite)
        db.session.commit()

        return json.dumps({'success': True}), 200, {'ContentType': 'application/json'}
    else:
        return json.dumps({'success': False}), 500, {'ContentType': 'application/json'}

@app.route("/tables/mytables/api/get_mytables")
@login_required
def get_mytables():
    mytables = MetaboliteTables.query.filter(MetaboliteTables.owner_id== g.user.id).all()
    json_dict = {"data": []}

    for table in mytables:
        table_dict = {"Title" : table.title, "id" : table.id, "Creation Date" : table.creation_date,
                      "Public" : table.public}
        json_dict["data"].append(table_dict)
    return jsonify(json_dict)


@app.route("/tables/new", methods=["GET", "POST"])
@login_required
def create_table():
    form = CreateTableForm()
    if form.validate_on_submit():
        metabolite_table = MetaboliteTables(
            title=form.title.data,
            description=form.description.data,
            species=form.species.data,
            doi=form.doi.data,
            owner_id=g.user.id,
            public = form.public.data,
        )
        db.session.add(metabolite_table)
        db.session.commit()
        flash("New table created")
        return redirect(url_for("mytables"))

    return render_template("tables/new.html", form=form)

@app.route("/tables/DdbT<id>/edit", methods=["GET", "POST"])
@login_required
def edit_table(id):
    metabolite_table = MetaboliteTables.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id:
        form = EditTableForm()
        if form.validate_on_submit():
            metabolite_table.title = form.title.data,
            metabolite_table.description = form.description.data,
            metabolite_table.species = form.species.data,
            metabolite_table.doi = form.doi.data
            db.session.commit()
            flash("Table edited!")
            return render_template("tables/edit.html", metabolite_table=metabolite_table, form=form)
        return render_template("tables/edit.html", metabolite_table=metabolite_table , form=form)
    else:
        abort(500)

@app.route("/tables/DdbT<id>")
def view_table(id):
    table_info = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.description,
        MetaboliteTables.doi,
        MetaboliteTables.creation_date,
        MetaboliteTables.owner_id,
        MetaboliteTables.species,
        MetaboliteTables.public,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.id == id).first_or_404()

    return render_template("tables/view.html", table_info = table_info)


@app.route("/tables/DdbT<id>/api/get_metabolites")
def get_metabolites(id):
    table = MetaboliteTables.query.filter(id == id).first()
    if table.public == True or table.owner_id == g.user.id:
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
    else:
        abort(500)


