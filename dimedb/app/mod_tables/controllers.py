from flask import Blueprint, request, render_template, g, session, redirect, url_for
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

    metabolite_ids = [str(x.inchikey) for x in Metabolite.query.filter(Metabolite.table_id == id).all()]

    return render_template("tables/view.html", table_info = table_info)

