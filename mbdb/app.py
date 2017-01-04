from flask import Flask, render_template, abort, request

from flask_mongoengine import MongoEngine
from flask_mongorest import MongoRest
from flask_mongorest.views import ResourceView
from flask_mongorest.resources import Resource
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

class Metabolite(db.DynamicDocument):
    meta = {
        "collection" : "metabolites"
    }

    name = db.StringField()
    synonyms = db.StringField()
    smiles = db.StringField()
    origins = db.StringField()
    molecular_formula = db.StringField()
    accurate_mass = db.FloatField()

class MetaboliteResource(Resource):
    document = Metabolite
    filters = {
        "id" : [ops.Exact],
        "name" : [ops.Exact, ops.Startswith, ops.Contains],
        "origins" : [ops.In],
        "molecular_formula" : [ops.Contains, ops.Exact],
        "synonyms" : [ops.Contains, ops.Exact]
    }

@api.register(name="met", url="/met/")
class MetaboliteView(ResourceView):
    resource = MetaboliteResource
    methods = [methods.Fetch, methods.List]

class Ionisation(db.DynamicDocument):
    meta = { "collection" : "metabolites" }
    name = db.StringField()
    origins = db.StringField()
    molecular_formula = db.StringField()
    adducts = db.StringField()



class IonisationResource(Resource):
    document = Ionisation
    filters = {
        "name" : [ops.Exact],
        "origins" : [ops.Contains]
    }

@api.register(name="ion", url="/ion/")
class IonisationView(ResourceView):
    resource = IonisationResource
    methods = [methods.Fetch, methods.List]



@app.route("/")
def homepage():
    return render_template("main.html", n=Metabolite.objects.count(),
                           url = request.url)

@app.route("/metabolites/")
def metabolites():
    return render_template("metabolites.html")

@app.route("/api/")
def api():
    return render_template("api.html")

# TEST
@app.route("/between/")
def between():
    am =  Metabolite.objects(accurate_mass__mod=[300.00, 288.12])
    return str(all_objects)

@app.route("/a/")
def anion():
    am = Ionisation.objects(adducts_Anion_length__gt=0)
    print am
    '''
    replace [M-H]1- with keys.
    db.metabolites.find({"adducts.Anion.[M-H]1-.0.0" : {$mod : [555, 400]}})
    :return:
    '''
    abort(501)

if __name__ == "__main__":
    app.run()