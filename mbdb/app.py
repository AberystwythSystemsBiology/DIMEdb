from flask import Flask, render_template, abort, request, jsonify

from flask_mongoengine import MongoEngine
from flask_mongorest import MongoRest
from flask_mongorest.views import ResourceView
from flask_mongorest.resources import Resource
from flask_mongorest.operators import Operator
from flask_mongorest import operators as ops
from flask_mongorest import methods
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


class AccurateMassSearch(Operator):
    op = "ppm"
    def prepare_queryset_kwargs(self, field, value, negate=False):
        if value == None:
            value = [0,0]
        else:
            value = [float(x) for x in value.split(',')]

        mz, ppm_threshold = value

        difference = abs(mz * (ppm_threshold * 0.0001))  # PPM to diff.
        return {
            field + '__gt': mz-difference,
            field + '__lt': mz+difference
        }


class MetaboliteBasicResource(Resource):
    document = MetaboliteBasic
    filters = {
        "name" : [ops.Exact, ops.Startswith, ops.Contains],
        "origins" : [ops.Exact],
        "molecular_formula" : [ops.Exact, ops.Contains, ops.IContains],
        "accurate_mass" : [ops.Exact, AccurateMassSearch, ops.Gte, ops.Gt]
    }


@api.register(name="metabolitesapi", url="/metabolites/")
class MetaboliteBasicView(ResourceView):
    resource = MetaboliteBasicResource
    methods = [methods.List, methods.Fetch]

class MetaboliteAdduct(db.DynamicDocument):
    meta = {"collection": "metabolites"}
    name = db.StringField()
    adduct_weights = db.StringField()

class MetaboliteAdductResource(Resource):
    document = MetaboliteAdduct
    filters = {
        "adduct_weights" : [ops.Contains]
    }

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
    app.run(host="0.0.0.0", port=5000)