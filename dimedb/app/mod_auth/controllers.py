from flask import Blueprint, request, render_template, flash, g, abort, redirect, url_for

from app import app, db, bcrypt, login_manager
from flask_login import login_user, logout_user, current_user, login_required

from app.mod_auth.models import User
from app.mod_auth.forms import LoginForm, RegistrationForm, ResendConfirmationForm, ResetPasswordForm, EditAccountForm

from token import generate_confirmation_token, confirm_token
from email import send_email

authentication = Blueprint("auth", __name__, url_prefix="/auth")

@authentication.before_request
def temporary_disabled():
    return abort(503)

@authentication.route("/login", methods=["GET", "POST"])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(user_name=form.user_name.data).first()
        if user and user.is_correct_password(form.password.data):
            if user.is_confirmed() == False:
                flash("Still awaiting account approval")
            elif user.is_banned():
                flash("Account has been banned")
            else:
                login_user(user)
                return redirect(url_for("homepage"))
        else:
            flash("Incorrect username or password given")
    return render_template("auth/login.html", form=form)

@login_required
@authentication.route("/account_management", methods=["GET", "POST"])
def management():
    form = EditAccountForm(request.form)
    if form.validate_on_submit():
        user = User.query.filter_by(email_address = g.user.email_address).first_or_404()
        user.first_name = form.first_name.data
        user.mid_initials = form.mid_initials.data
        user.last_name = form.last_name.data
        user.phone = form.phone.data
        user.affiliation = form.affiliation.data
        user.address = form.address.data
        user.country = form.country.data
        db.session.commit()
        flash("Account details have been successfuly changed!")
    return render_template("auth/management/management.html", form=form)

@authentication.route("/register", methods=["GET", "POST"])
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
                country=form.country.data
            )
            db.session.add(user)
            db.session.commit()

            token = generate_confirmation_token(user.email_address)
            confirm_url = url_for("auth.confirm_email", token=token, _external=True)
            html = render_template("auth/emails/account_confirmation.html", confirm_url=confirm_url)
            subject = "Please Confirm your email"
            send_email(user.email_address, subject, html)
            flash("Please check your email address for account confirmation")
            return redirect(url_for("auth.login"))
    return render_template("auth/register.html", form=form)

@authentication.route("/logout")
def logout():
    logout_user()
    return redirect(url_for("homepage"))

@authentication.route("/confirm/<token>")
def confirm_email(token):
    email = confirm_token(token)
    if email == False:
        flash("The confirmation link is invalid or has expired.")
        return redirect(url_for("auth.login"))
    else:
        user = User.query.filter_by(email_address = email).first_or_404()
        if user.is_confirmed() == True:
            flash("Already confirmed, feel free to login.")
        else:
            user.confirmed = True
            db.session.add(user)
            db.session.commit()
            flash("Account has been successully confirmed, feel free to login")
        return redirect(url_for("auth.login"))

@authentication.route("/resend_confirmation", methods=["GET", "POST"])
def resend_confirmation():
    form = ResendConfirmationForm(request.form)
    if form.validate_on_submit():
        user = User.query.filter_by(email_address = form.email_address.data, confirmed=False).first()
        if user != None:
            token = generate_confirmation_token(user.email_address)
            confirm_url = url_for("auth.confirm_email", token=token, _external=True)
            html = render_template("auth/emails/account_confirmation.html", confirm_url=confirm_url)
            subject = "DIMEdb: Account Confirmation Required"
            send_email(user.email_address, subject, html)
        flash("Confirmation email resent, please check your email address")
        return redirect(url_for("auth.login"))
    return render_template("auth/management/resend_confirmation.html", form=form)



@authentication.route("/reset_password", methods=["GET", "POST"])
def reset_password_email():
    form = ResendConfirmationForm(request.form)
    if form.validate_on_submit():
        user = User.query.filter_by(email_address = form.email_address.data).first()
        if user != None:
            token = generate_confirmation_token(user.email_address)
            reset_url = url_for("auth.reset_password", token=token, _external=True)
            html = render_template("auth/emails/password_reset.html", reset_url=reset_url)
            subject = "DIMEdb: Account Password Reset"
            send_email(user.email_address, subject, html)
        flash("Password reset email sent, please check your email address")
        return redirect(url_for("auth.login"))
    return render_template("auth/management/password_email.html", form=form)


@authentication.route("/reset_password/<token>", methods=["GET", "POST"])
def reset_password(token):
    form = ResetPasswordForm(request.form)
    email = confirm_token(token)
    if email != None:
        if form.validate_on_submit():
            user = User.query.filter_by(email_address = email).first_or_404()
            user.password = form.password.data
            db.session.add(user)
            db.session.commit()
            flash("Password successfully changed, please sign in.")
            return redirect(url_for("auth.login"))
    return render_template("auth/management/reset_password.html", form=form)

@login_manager.user_loader
def load_user(id):
    return User.query.get(id)

@app.before_request
def before_request():
    g.user = current_user