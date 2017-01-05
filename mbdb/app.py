from flask import Flask, render_template, abort, request, jsonify

from flask_mongoengine import MongoEngine
from flask_mongorest import MongoRest
from flask_mongorest.views import ResourceView
from flask_mongorest.resources import Resource
from flask_mongorest import operators as ops
from flask_mongorest import methods
from flask_mongorest.operators import Operator

app = Flask(__name__)

app.config.update(
    MONGODB_HOST = "localhost",
    MONGODB_PORT = 27017,
    MONGODB_DB = "mbdb",
    DEBUG = True
)

db = MongoEngine(app)
api = MongoRest(app)


class MetaboliteBasic(db.DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = db.StringField()
    origins = db.ListField(db.StringField())
    molecular_formula = db.StringField()
    accurate_mass = db.FloatField()


class MetaboliteBasicResource(Resource):
    document = MetaboliteBasic
    filters = {
        "name" : [ops.Exact, ops.Startswith, ops.Contains],
        "origins" : [ops.Exact],
        "molecular_formula" : [ops.Exact, ops.Contains],
        "accurate_mass" : [ops.Gt, ops.Lt]
    }


@api.register(name="metabolitesapi", url="/metabolites/")
class MetaboliteBasicView(ResourceView):
    resource = MetaboliteBasicResource
    methods = [methods.List]

class MetaboliteAdduct(db.DynamicDocument):
    meta = {"collection": "metabolites"}
    name = db.StringField()
    molecular_formula = db.StringField()
    adduct_weights = db.StringField()

class MetaboliteAdductResource(Resource):
    document = MetaboliteAdduct
    filters = {}

@api.register(name="adductsapi", url="/adducts/")
class MetaboliteAdductView(ResourceView):
    resource =  MetaboliteAdductResource
    methods = [methods.List]

# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("main.html")

@app.route("/api/")
def api():
    return render_template("api.html", url = request.url)

if __name__ == "__main__":
    app.run()