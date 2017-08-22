from flask import Blueprint, abort, render_template, flash, g, session, redirect, url_for
from app import app, db, bcrypt, login_manager
from flask_login import login_required

from app.mod_auth.models import User

from decorators import check_admin

admin = Blueprint("admin", __name__, url_prefix="/admin")

@admin.route("/")
@login_required
@check_admin
def home():
    users = User.query.filter().all()
    return render_template("admin/index.html", users=users)