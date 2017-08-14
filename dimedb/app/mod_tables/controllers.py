import json

from flask import Blueprint, flash, render_template, g, abort, redirect, url_for, jsonify, request
from app import app, db

from flask_login import login_required

from app.mod_tables.models import MetaboliteTable, Metabolite, MetaboliteTablePublication

from app.mod_auth.models import User

from app.mod_tables.forms import CreateTableForm, EditTableForm


tables = Blueprint("tables", __name__)

@app.route("/tables")
def public_tables():
    public_tables = User.query.join(MetaboliteTable, User.id == MetaboliteTable.owner_id).add_columns(
        MetaboliteTable.id,
        MetaboliteTable.title,
        MetaboliteTable.creation_date,
        MetaboliteTable.species,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTable.public == True, MetaboliteTable.removed == False).all()
    return render_template("tables/index.html", public_tables=public_tables)

@app.route("/tables/mytables")
@login_required
def mytables():
    mytables = MetaboliteTable.query.filter(MetaboliteTable.owner_id== g.user.id, MetaboliteTable.removed == False).all()
    return render_template("tables/mytables.html", mytables=mytables)

@app.route("/tables/mytables/api/add_metabolite", methods=["POST"])
@login_required
def add_metabolite():
    payload = request.form
    table = MetaboliteTable.query.filter(MetaboliteTable.owner_id == g.user.id).first()
    if table != None:
        metabolite = Metabolite(
            table_id = payload["Table ID"],
            inchikey = payload["InChI Key"],
            comments = payload["Comment"]
        )

        db.session.add(metabolite)
        db.session.commit()

        return json.dumps({'success': True}), 200, {'ContentType': 'application/json'}
    else:
        return json.dumps({'success': False}), 500, {'ContentType': 'application/json'}

@app.route("/tables/mytables/api/get_mytables")
@login_required
def get_mytables():
    mytables = MetaboliteTable.query.filter(MetaboliteTable.owner_id== g.user.id, MetaboliteTable.removed == False).all()
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
        metabolite_table = MetaboliteTable(
            title=form.title.data,
            description=form.description.data,
            species=form.species.data,
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
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.removed != True:
        form = EditTableForm()
        if form.validate_on_submit():
            metabolite_table.title = form.title.data,
            metabolite_table.description = form.description.data,
            metabolite_table.species = form.species.data,
            db.session.commit()
            flash("Table edited!")
            return render_template("tables/edit.html", metabolite_table=metabolite_table, form=form)
        return render_template("tables/edit.html", metabolite_table=metabolite_table , form=form)
    else:
        abort(500)

@app.route("/tables/DdbT<id>")
def view_table(id):
    table_info = User.query.join(MetaboliteTable, User.id == MetaboliteTable.owner_id).add_columns(
        MetaboliteTable.id,
        MetaboliteTable.title,
        MetaboliteTable.description,
        MetaboliteTable.creation_date,
        MetaboliteTable.owner_id,
        MetaboliteTable.species,
        MetaboliteTable.public,
        MetaboliteTable.removed,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTable.id == id).first_or_404()
    if table_info.removed != True:
        return render_template("tables/view.html", table_info = table_info)
    else:
        abort(404)
@app.route("/tables/DdbT<id>/remove")
def delete_table(id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)

    if metabolite_table.owner_id == g.user.id:
        metabolite_table.removed = True
        db.session.commit()
        flash("You have succesfully removed the metabolite table")
        return redirect(url_for("mytables"))
    else:
        abort(500)

@app.route("/tables/DdbT<id>/api/get_metabolites")
def get_metabolites(id):
    table = MetaboliteTable.query.filter(id == id).first()
    if table.public == True or table.owner_id == g.user.id and table.removed != True:
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


