from flask_wtf import FlaskForm

from wtforms import StringField, PasswordField, RadioField

from wtforms.validators import DataRequired, EqualTo, Length
from wtforms.fields.html5 import  EmailField

class LoginForm(FlaskForm):
    user_name = StringField("User Name", [DataRequired("Forgot your user name?")])
    password = PasswordField("Password", [DataRequired("You must provide a password!")])

class RegistrationForm(FlaskForm):
    user_name = StringField("User Name", [DataRequired("You must provide a username!")])
    password = PasswordField("Password", [Length(min=4),
                                          DataRequired("You must provide a password!"),
                                          EqualTo("password_confirm", message="Passwords must match")])
    password_confirm = PasswordField("Repeat Password")
    email_address = EmailField("Email", [DataRequired("Please provide your email address!")])

    affiliation = StringField("Affiliation", [Length(min=0)])
    address = StringField("Address")
    country = StringField("Country", [Length(min=4)])

    phone = StringField("Telephone Number", [Length(min=0, max=30)])

    first_name = StringField("First Name", [DataRequired("You must provide your first name!")])
    mid_initials = StringField("Middle Initials", [Length(min=0, max=5)])
    last_name = StringField("Last Name", [DataRequired("You must provide your last name!")])

class ResendConfirmationForm(FlaskForm):
    email_address = EmailField("Email", [DataRequired("Please provide your email address!")])

class ResetPasswordForm(FlaskForm):
    password = PasswordField("Password", [Length(min=4),
                                          DataRequired("You must provide a password!"),
                                          EqualTo("password_confirm", message="Passwords must match")])
    password_confirm = PasswordField("Repeat Password")