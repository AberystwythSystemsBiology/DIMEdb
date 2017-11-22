from eve import Eve

from flask import render_template, request, send_from_directory, abort
from flask_sqlalchemy import SQLAlchemy
from flask_bcrypt import Bcrypt
from flask_login import LoginManager
from hurry.filesize import size, si
from flask_mail import Mail

from config import BaseConfig

import os

app = Eve(__name__)
app.config.from_object(BaseConfig)

db = SQLAlchemy(app)
bcrypt = Bcrypt(app)

mail = Mail(app)

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "auth.login"


@app.before_request
def check_for_maintenance():
    if app.config["MAINTENANCE"]:
        return abort(503)

@app.route("/")
def homepage():
    return render_template("home.html", count='{:,}'.format(app.data.driver.db['metabolites'].count()))

@app.route("/search/mass")
def mass_search():
    return render_template("search/mass.html")

@app.route("/search/text")
def text_search():
    return render_template("search/text.html")

@app.route("/help/")
def help():
    return render_template("help.html", url = request.url_root)

@app.route("/help/downloads")
def downloads_centre():
    d = os.path.expanduser("~/.data/dimedb/downloads/")
    files = [[f, size(os.stat(d+f).st_size, system=si), os.stat(d+f).st_ctime] for f in os.listdir(d)]
    return render_template("misc/downloads.html", url = request.url_root, files = files)

@app.route("/help/downloads/<string:fn>")
def get_file(fn):
    d = os.path.expanduser("~/.data/dimedb/downloads/")
    return send_from_directory(d, fn)


@app.route("/tools/h/")
def isotopic_distribution_calculator():
    return render_template("tools/isotopic_distribution.html")

@app.errorhandler(404)
def page_not_found(e):
    return render_template("./misc/404.html"), 404

@app.errorhandler(403)
def forbidden(e):
    return render_template("./misc/403.html"), 403

@app.errorhandler(503)
def maintenance(e):
    return render_template("./misc/503.html"), 503


@app.route("/test")
def test():
    return render_template("./misc/test.hml")

# Blueprints

from mod_auth.controllers import authentication
from mod_tables.controllers import tables
from mod_admin.controllers import admin
from mod_view.controllers import view


app.register_blueprint(authentication)
app.register_blueprint(tables)
app.register_blueprint(admin)
app.register_blueprint(view)

db.create_all()
