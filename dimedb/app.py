from flask import Flask, render_template, abort, request, jsonify, url_for

from flask_mongoengine import MongoEngine
from flask_mongorest import MongoRest
from flask_mongorest.views import ResourceView
from flask_mongorest import methods

app = Flask(__name__)

import resources as r, documents as d

app.config.update(
    MONGODB_HOST = "localhost",
    MONGODB_PORT = 27017,
    MONGODB_DB = "dimedb",
    DEBUG = True
)

is_maintenance_mode = False

db = MongoEngine(app)
api = MongoRest(app)

@api.register(name="metabolites_f", url="/api/metabolite/")
class MetaboliteFullView(ResourceView):
    resource = r.MetaboliteFullResource
    methods = [methods.List]

@api.register(name="metabolites", url="/api/metabolites/")
class MetaboliteBasicView(ResourceView):
    resource = r.MetaboliteBasicResource
    methods = [methods.List, methods.Fetch]

@api.register(name="adducts", url="/api/adducts/")
class Adducts(ResourceView):
    resource = r.MetaboliteAdductsResource
    methods = [methods.List, methods.Fetch]

# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("main.html")

@app.route("/search/")
def search():
    return render_template("search.html")


@app.route("/help/")
def help():
    return render_template("help.html", url = request.url_root)

@app.route("/view/<string:_id>/")
def view(_id):
    try:
        url = request.url_root
        return render_template("view.html" , id=_id, base_url=url)
    except Exception, err:
        abort(403)

@app.errorhandler(404)
def page_not_found(e):
    return render_template("./misc/404.html"), 404

@app.before_request
def maintenence_page():
    if is_maintenance_mode:
        return render_template("./misc/503.html"), 503


# Temporary until I find a better solution.
@app.route("/gen_structure/<string:id>/")
def smiles_to_2d(id):
    from rdkit import Chem
    from rdkit.Chem import Draw
    import StringIO
    try:
        smiles = d.MetaboliteBasic.objects.filter(_id=id)[0].smiles
        smiles_image = StringIO.StringIO()
        mol = Chem.MolFromSmiles(smiles)
        Draw.MolToFile(mol, fileName=smiles_image, imageType="png", size=(300, 300))
        return smiles_image.getvalue().encode("base64")
    except Exception, err:
        abort(404)
