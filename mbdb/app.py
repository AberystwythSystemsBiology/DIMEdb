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
        "origins" : [ops.Contains],
        "molecular_formula" : [ops.Contains, ops.Exact],
        "synonyms" : [ops.Contains, ops.Exact]
    }

@api.register(name="metabolites", url="/metabolites/")
class MetaboliteView(ResourceView):
    resource = MetaboliteResource
    methods = [methods.Fetch, methods.List]

@api.register(name="masses", url="/masses/")
class MassesView(ResourceView):
    pass

@api.register(name="ioinisation", url="/ionisation/")
class IonisationView(ResourceView):
    pass



@app.route("/")
def homepage():
    return render_template("main.html", n=Metabolite.objects.count(),
                           url = request.url)

# TEST
@app.route("/between/")
def between():
    am =  Metabolite.objects(accurate_mass__mod=[300.00, 255.12])
    print len([x["name"] for x in am])
    abort(501)

if __name__ == "__main__":
    app.run()