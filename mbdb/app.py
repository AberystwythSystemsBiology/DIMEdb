from flask import Flask, render_template, abort, request, jsonify

from flask_mongoengine import MongoEngine
from flask_pymongo import PyMongo
from flask_restful import Api, Resource

app = Flask(__name__)

app.config.update(
    MONGODB_HOST = "localhost",
    MONGODB_PORT = 27017,
    MONGODB_DB = "mbdb",
    DEBUG = True
)

db = MongoEngine(app)
mongo = PyMongo(app, config_prefix='MONGO')

class Metabolite(Resource):
    def get(self, name=None, origins=None):
        data = []
        if name:
            metabolite_info = mongo.db.metabolites.find({"name" : name})
            if metabolite_info:
                return jsonify({"status": "ok", "data": metabolite_info})


# API Blueprints
api = Api(app)
api.add_resource(Metabolite, "/met", endpoint="met")



# Annoying webpage stuff.

@app.route("/")
def homepage():
    return render_template("main.html")

@app.route("/metabolites/")
def metabolites():
    return render_template("metabolites.html", n=0)

@app.route("/api/")
def api():
    return render_template("api.html", url = request.url)

if __name__ == "__main__":
    app.run()