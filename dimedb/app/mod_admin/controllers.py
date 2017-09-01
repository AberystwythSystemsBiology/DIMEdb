from flask import Blueprint, abort, render_template, flash, g, session, redirect, url_for
from app import app, db, bcrypt, login_manager
from flask_login import login_required

from app.mod_auth.models import User
from app.mod_tables.models import MetaboliteTable

from decorators import check_admin

admin = Blueprint("admin", __name__, url_prefix="/admin")

@admin.route("/")
@login_required
@check_admin
def home():
    user_count = User.query.count()
    table_count = MetaboliteTable.query.count()
    metabolite_count = app.data.driver.db['metabolites'].count()
    return render_template("admin/index.html", user_count=user_count, table_count=table_count, metabolite_count=metabolite_count)

@admin.route("/users")
@login_required
@check_admin
def users():
    users = User.query.all()
    return render_template("admin/users.html", users=users)

@admin.route("/tables")
@login_required
@check_admin
def tables():
    tables = MetaboliteTable.query.all()
    return render_template("admin/tables.html", tables=tables)