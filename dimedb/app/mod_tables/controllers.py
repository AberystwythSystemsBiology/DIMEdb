import json

from flask import Blueprint, flash, render_template, g, abort, redirect, url_for, jsonify, request
from dimedb import app, db

from flask_login import login_required

from mod_tables.models import MetaboliteTable, Metabolite, MetaboliteTablePublication

from mod_auth.models import User

from mod_tables.forms import CreateTableForm, EditTableForm, NewPublicationForm


tables = Blueprint("tables", __name__, url_prefix="/tables")

@tables.route("/")
def public_tables():
    tables = User.query.join(MetaboliteTable, User.id == MetaboliteTable.owner_id).add_columns(
        MetaboliteTable.id,
        MetaboliteTable.title,
        MetaboliteTable.creation_date,
        MetaboliteTable.species,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTable.published == True, MetaboliteTable.removed == False).all()
    return render_template("tables/index.html", tables=tables)

@tables.route("/mytables")
@login_required
def mytables():
    mytables = MetaboliteTable.query.filter(MetaboliteTable.owner_id== g.user.id, MetaboliteTable.removed == False).all()
    return render_template("tables/mytables.html", mytables=mytables)

@tables.route("/mytables/api/add_metabolite", methods=["POST"])
@login_required
def add_metabolite():
    payload = request.form
    table = MetaboliteTable.query.filter(MetaboliteTable.owner_id == g.user.id).first()
    if table != None:
        metabolite = Metabolite(
            table_id = payload["Table ID"],
            inchikey = payload["InChI Key"],
            comments = payload["Comments"]
        )

        db.session.add(metabolite)
        db.session.commit()

        return json.dumps({'success': True}), 200, {'ContentType': 'application/json'}
    else:
        return json.dumps({'success': False}), 500, {'ContentType': 'application/json'}

@tables.route("/mytables/api/get_mytables")
@login_required
def get_mytables():
    mytables = MetaboliteTable.query.filter(MetaboliteTable.owner_id== g.user.id, MetaboliteTable.removed == False, MetaboliteTable.published == False).all()
    json_dict = {"data": []}

    for table in mytables:
        table_dict = {"Title" : table.title, "id" : table.id, "Creation Date" : table.creation_date,
                      "Published" : table.published}
        json_dict["data"].append(table_dict)
    return jsonify(json_dict)


@tables.route("/new", methods=["GET", "POST"])
@login_required
def create_table():
    form = CreateTableForm()
    if form.validate_on_submit():
        metabolite_table = MetaboliteTable(
            title=form.title.data,
            description=form.description.data,
            species=form.species.data,
            owner_id=g.user.id
        )
        db.session.add(metabolite_table)
        db.session.commit()
        flash("New table created")
        return redirect(url_for("tables.mytables"))

    return render_template("tables/meta/create.html", form=form)

@tables.route("/edit/DDBT<id>", methods=["GET", "POST"])
@login_required
def edit_table(id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.published == False and metabolite_table.owner_id == g.user.id and metabolite_table.removed == False:
        publications = MetaboliteTablePublication.query.filter(MetaboliteTablePublication.table_id == id).all()
        form = EditTableForm()
        if form.validate_on_submit():
            metabolite_table.title = form.title.data,
            metabolite_table.description = form.description.data,
            metabolite_table.species = form.species.data,
            db.session.commit()
            flash("Table edited!")
        return render_template("tables/meta/edit.html", metabolite_table=metabolite_table , form=form, publications=publications)
    return abort(403)

@tables.route("/edit/DDBT<id>/add_publication", methods=["GET", "POST"])
@login_required
def add_publication(id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.removed != True and metabolite_table.published != True:
        form = NewPublicationForm(request.form)
        if form.validate_on_submit():
            publication = MetaboliteTablePublication(
                table_id = id,
                pubmed_id = form.pubmed_id.data,
                doi = form.doi.data,
                author_list = form.author_list.data,
                publication_title = form.publication_title.data
            )
            db.session.add(publication)
            db.session.commit()
            flash("Publication Added")
            return redirect(url_for("tables.edit_table", id=id))
        return render_template("tables/meta/new_publication.html", id=id, form=form)
    else:
        abort(403)

@tables.route("/edit/DDBT<id>/remove_publication/<pub_id>")
@login_required
def remove_publication(id, pub_id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.removed != True and metabolite_table.published != True:
        publication = MetaboliteTablePublication.query.get_or_404(pub_id)
        db.session.delete(publication)
        db.session.commit()
        flash("Publication Successfully Removed")
        return redirect(url_for("tables.edit_table", id=id))
    else:
        abort(403)

@tables.route("/edit/DDBT<id>/remove_metabolite/<m_id>")
@login_required
def remove_metabolite(id, m_id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.removed != True and metabolite_table.published != True:
        metabolite = Metabolite.query.get_or_404(m_id)
        db.session.delete(metabolite)
        db.session.commit()
        flash("Metabolite Successfully Removed")
        return redirect(url_for("tables.edit_table", id=id))
    else:
        abort(403)

@tables.route("/view/DDBT<id>")
def view_table(id):
    table_info = User.query.join(MetaboliteTable, User.id == MetaboliteTable.owner_id).add_columns(
        MetaboliteTable.id,
        MetaboliteTable.title,
        MetaboliteTable.description,
        MetaboliteTable.creation_date,
        MetaboliteTable.owner_id,
        MetaboliteTable.species,
        MetaboliteTable.published,
        MetaboliteTable.removed,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTable.id == id, MetaboliteTable.removed == False).first_or_404()

    publications = MetaboliteTablePublication.query.filter(MetaboliteTablePublication.table_id == id).all()

    return render_template("tables/view.html", table_info = table_info, publications=publications)


@tables.route("/DDBT<id>/edit/publish")
def publish_table(id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.published == False:
        metabolite_table.published = True
        db.session.commit()
        flash("The table has been successfully published")
        return redirect(url_for("tables.mytables"))
    else:
        abort(403)

@tables.route("/DDBT<id>/edit/remove")
def delete_table(id):
    metabolite_table = MetaboliteTable.query.get_or_404(id)
    if metabolite_table.owner_id == g.user.id and metabolite_table.published == False:
        metabolite_table.removed = True
        db.session.commit()
        flash("You have succesfully removed the metabolite table")
        return redirect(url_for("tables.mytables"))
    else:
        abort(403)

@tables.route("/api/DDBT<id>/get_metabolites")
def get_metabolites(id):
    table = MetaboliteTable.query.filter(id == id).first()
    if table.published == True or table.owner_id == g.user.id and table.removed != True:
        metabolites = Metabolite.query.filter(Metabolite.table_id == id).all()
        json_dict = {"data" : []}
        for tm in metabolites:
            metabolite_dict = {"InChIKey": None, "Name": None, "Comments" : None, "Molecular Formula": None, "Table ID" : None}
            metabolite = app.data.driver.db["metabolites"].find_one({'_id': tm.inchikey})
            metabolite_dict["InChIKey"] = tm.inchikey
            metabolite_dict["Comments"] = tm.comments
            metabolite_dict["Name"] = metabolite["Identification Information"]["Name"]
            metabolite_dict["Molecular Formula"] = metabolite["Identification Information"]["Molecular Formula"]
            metabolite_dict["Table ID"] = tm.id
            json_dict["data"].append(metabolite_dict)
        return jsonify(json_dict)
    else:
        abort(500)


