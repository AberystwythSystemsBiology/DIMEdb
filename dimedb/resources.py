import documents as d, operators as o
from flask_mongorest import operators as ops
from flask_mongorest.resources import Resource

class AdductResource(Resource):
    document = d.Adduct

class AdductsResource(Resource):
    document = d.Adducts

    related_resources = {
        "positive": AdductResource,
        "negative": AdductResource,
        "neutral": AdductResource
    }

class MetaboliteFullResource(Resource):
    document = d.MetaboliteFull
    max_limit, default_limit = [1000, 1000]

    filters = {
        "id" : [ops.Exact],
        "adducts" : [o.AdductPpm]
    }

    related_resources = {
        "adducts" : AdductsResource
    }