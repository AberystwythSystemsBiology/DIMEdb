from eve import Eve
from flask import render_template, request, send_from_directory, g

from flask_sqlalchemy import SQLAlchemy
from flask.ext.bcrypt import Bcrypt
from flask.ext.login import LoginManager


from hurry.filesize import size, si

import os

app = Eve(__name__)
app.config.update(
    DEBUG = True,
    SQLALCHEMY_DATABASE_URI = "postgresql+psycopg2://keo7:password@localhost/dimedb",
    SQLALCHEMY_TRACK_MODIFICATIONS = False,
    LDAP_LOGIN_VIEW = 'auth.login',
    CSRF_ENABLED = True,
    CSRF_SESSION_KEY = "secret",
    SECRET_KEY = "secret"
)

db = SQLAlchemy(app)
bcrypt = Bcrypt(app)

login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = "login"

@app.route("/")
def homepage():
    return render_template("home.html", count='{:,}'.format(app.data.driver.db['metabolites'].count()))

@app.route("/search/mass")
def mass_search():
    return render_template("search/mass.html")

@app.route("/search/text")
def text_search():
    return render_template("search/text.html")

@app.route("/test")
def test_page():
    return render_template("misc/test.html")

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

@app.route("/view/<string:_id>/")
def view(_id):
        return render_template("view/view.html" , id=_id)

@app.route("/view/structure/<string:id>")
def get_structures_image(id):
    d = os.path.expanduser("~/.data/dimedb/structures/")
    return send_from_directory(d, id+".svg")

@app.route("/tools/isotopic_distribution/")
def isotopic_distribution_calculator():
    return render_template("tools/isotopic_distribution.html")

@app.errorhandler(404)
def page_not_found(e):
    return render_template("./misc/404.html"), 404

# Blueprints

db.create_all()