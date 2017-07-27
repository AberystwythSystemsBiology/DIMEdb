from eve import Eve

from flask import render_template, abort, request, send_from_directory
from hurry.filesize import size, si
import urllib, json, os

app = Eve(__name__)
app.config.update(
    DEBUG = True
)

# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("home.html")

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

@app.errorhandler(404)
def page_not_found(e):
    return render_template("./misc/404.html"), 404
