from flask import Flask, render_template

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

class MetaboliteResource(Resource):
    document = Metabolite
    filters = {
        "name" : [ops.Exact]
    }

@api.register(name="metabolites", url="/metabolites/")
class MetaboliteView(ResourceView):
    resource = MetaboliteResource
    methods = [methods.Fetch, methods.List]

@app.route("/")
def homepage():
    return render_template("main.html", n=len(Metabolite.objects))

if __name__ == "__main__":
    app.run()