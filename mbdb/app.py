from flask import Flask, render_template, abort, request, jsonify

from flask_mongoengine import MongoEngine
from flask_mongorest import MongoRest
from flask_mongorest.views import ResourceView
from flask_mongorest import methods

app = Flask(__name__)

import resources as r

app.config.update(
    MONGODB_HOST = "localhost",
    MONGODB_PORT = 27017,
    MONGODB_DB = "mbdb",
    DEBUG = True
)

db = MongoEngine(app)
api = MongoRest(app)

@api.register(name="metabolites", url="/api/metabolites/")
class MetaboliteBasicView(ResourceView):
    resource = r.MetaboliteBasicResource
    methods = [methods.List, methods.Fetch]

# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("main.html")

@app.route("/api/")
def api():
    return render_template("api.html", url = request.url)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)