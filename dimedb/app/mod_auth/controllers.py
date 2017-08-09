from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for

from app import app, db, bcrypt, login_manager
from app.mod_auth.forms import LoginForm, RegistrationForm
from flask.ext.login import login_user, logout_user, current_user, login_required

from app.mod_auth.models import User

authentication = Blueprint("authentication", __name__)

@app.route("/login", methods=["GET", "POST"])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(user_name=form.user_name.data).first()
        if user and user.is_correct_password(form.password.data):
            if user.confirmed == False:
                flash("Still waiting for account approval")
            elif user.banned == True:
                flash("Your account has been banned")
            elif user.disabled == True:
                flash("Your account has been disabled")
            else:
                login_user(user)
                flash("Login successful")
                return redirect(url_for("index"))
        else:
            flash("Incorrect email or password given")
    return render_template("auth/login.html", form=form)

@app.route("/register", methods=["GET", "POST"])
def register():
    pass

@app.route("/logout")
def logout():
    logout_user()
    return redirect(url_for("homepage"))

@login_manager.user_loader
def load_user(id):
    return User.query.get(id)

@app.before_request
def before_request():
    g.user = current_user