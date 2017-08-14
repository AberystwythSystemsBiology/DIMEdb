from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for

from app import app, db, bcrypt, login_manager
from flask.ext.login import login_user, logout_user, current_user, login_required

from app.mod_auth.models import User
from app.mod_auth.forms import LoginForm, RegistrationForm

authentication = Blueprint("authentication", __name__)

@app.route("/login", methods=["GET", "POST"])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(user_name=form.user_name.data).first()
        if user and user.is_correct_password(form.password.data):
            if user.confirmed == False:
                flash("Still waiting for account approval")
            else:
                login_user(user)
                return redirect(url_for("homepage"))
        else:
            flash("Incorrect email or password given")
    return render_template("auth/login.html", form=form)

@login_required
def management():
    return render_template("auth/management.html")

@app.route("/register", methods=["GET", "POST"])
def register():
    form = RegistrationForm(request.form)
    if form.validate_on_submit():
        if User.query.filter_by(email_address=form.email_address.data).first():
            flash("Email address already used!")
        elif User.query.filter_by(user_name=form.user_name.data).first():
            flash("Username already in use!")
        else:
            user = User(
                user_name=form.user_name.data,
                password=form.password.data,
                email_address=form.email_address.data,
                first_name=form.first_name.data,
                address = form.address.data,
                mid_initials = form.mid_initials.data,
                phone = form.phone.data,
                last_name=form.last_name.data,
                affiliation=form.affiliation.data,
                country=form.country.data,
                user_type=form.user_type.data
            )
            db.session.add(user)
            db.session.commit()
            flash("Thank you for registering!")
            return redirect(url_for("login"))
    return render_template("auth/register.html", form=form)

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