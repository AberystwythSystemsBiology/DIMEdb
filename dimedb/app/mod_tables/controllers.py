from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for
from app import app, db

from app.mod_tables.models import MetaboliteTables, Metabolite

tables = Blueprint("tables", __name__)

@app.route("/tables")
def tables_home():
    return render_template("tables/index.html")
