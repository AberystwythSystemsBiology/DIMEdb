import os

from eve import Eve
from eve_elastic import Elastic

from flask import render_template, request, send_from_directory, abort
from hurry.filesize import size, si
from flask_mail import Mail

from config import BaseConfig


app = Eve(__name__)
app.config.from_object(BaseConfig)


from mod_view.controllers import view

@app.before_request
def check_for_maintenance():
    if app.config["MAINTENANCE"]:
        return abort(503)


@app.route("/")
def homepage():
    count = '{:,}'.format(app.data.driver.db['metabolites'].count())

    metabolites = app.data.driver.db["metabolites"]

    #print metabolites.find_one({"Adducts" : {"$elemMatch" : {"Polarity" : "Positive", "Adduct" : {"$in" : ["[M+H]1+","[M+K]1+","[M+Na]1+"]},"Accurate Mass" : {"$lte" : 362, "$gte" : 360}}}})

    return render_template("home.html", count=count)


@app.route("/search/mass")
def mass_search():
    return render_template("search/mass.html")


@app.route("/search/text")
def text_search():
    return render_template("search/text.html")


@app.route("/help/")
def help():
    return render_template("help.html", url=request.url_root)


@app.route("/help/downloads")
def downloads_centre():
    d = os.path.expanduser("~/.data/dimedb/downloads/")
    files = []
    for f in os.listdir(d):
        rel_path = os.path.join(d, f)
        size = size(os.stat(rel_path).st_size, system=si)
        t = os.stat(rel_path).st_ctime
        files.append([f, size, t])
    url = request.url_root
    return render_template("misc/downloads.html", url=url, files=files)


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


app.register_blueprint(view)
