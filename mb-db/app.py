'''
    This is meant to be as lightweight as humanly possible.
'''

import os

from flask import Flask, request

from flask_mongoengine import MongoEngine

app = Flask(__name__)

app.url_map.strict_slashes = False

app.config.update(
    DEBUG = True,
    TESTING = True,
    MONGODB_SETTINGS = {
        "HOST" : "localhost",
        "PORT" : 27017,
        "DB" : "mbdb",

    },
)

db = MongoEngine(app)

class Metabolite():
    pass

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)