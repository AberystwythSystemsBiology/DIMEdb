from flask import Blueprint, request, render_template, flash, g, session, redirect, url_for

from app import app, db, bcrypt, login_manager
from app.modules.auth.forms import LoginForm, RegistrationForm, AccountManagementForm
from flask.ext.login import login_user, logout_user, current_user, login_required

from app.modules.auth.models import User

authentication = Blueprint("auth", __name__, url_prefix="/auth")