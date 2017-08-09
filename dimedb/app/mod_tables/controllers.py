from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for
from app import app, db

from app.mod_tables.models import MetaboliteTables, Metabolite
from app.mod_auth.models import User

tables = Blueprint("tables", __name__)

@app.route("/tables")
def tables_home():
    public_tables = MetaboliteTables.query.filter(MetaboliteTables.public == True).all()

    public_tables = User.query.join(MetaboliteTables, User.id == MetaboliteTables.owner_id).add_columns(
        MetaboliteTables.id,
        MetaboliteTables.title,
        MetaboliteTables.creation_date,
        User.first_name,
        User.last_name
    ).filter(MetaboliteTables.public == True).all()

    return render_template("tables/index.html", public_tables=public_tables)
