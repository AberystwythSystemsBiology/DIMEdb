import os
from flask import Blueprint, send_from_directory, render_template

view = Blueprint("view", __name__, url_prefix="/view")


@view.route("/<inchikey>")
def view_metabolite(inchikey):
    return render_template("view/view.html", id=inchikey)

@view.route("/structure/<string:id>")
def get_structures_image(id):
    d = os.path.expanduser("~/.data/dimedb/structures/")
    return send_from_directory(d, id + ".svg")