from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for

from app import app, db, bcrypt, login_manager
from app.mod_auth.forms import LoginForm, RegistrationForm
from flask.ext.login import login_user, logout_user, current_user, login_required

from app.mod_auth.models import User

authentication = Blueprint("auth", __name__, url_prefix="/auth")