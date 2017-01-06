import documents as d, operators as o
from flask_mongorest import operators as ops
from flask_mongorest.resources import Resource

class MetaboliteBasicResource(Resource):
    document = d.MetaboliteBasic
    filters = {
        "name" : [ops.Exact, ops.Startswith, ops.Contains],
        "origins" : [ops.Exact],
        "molecular_formula" : [ops.Exact, ops.Contains, ops.IContains],
        "accurate_mass" : [ops.Exact, o.AccurateMassSearch, ops.Gte, ops.Gt]
    }

class PositiveResource(Resource):
    document = d.Positive

class NegativeResource(Resource):
    document = d.Negative

class AdductWeightsResource(Resource):
    document = d.AdductWeights
    related_resources = {
        "positive" : PositiveResource,
        "negative" : NegativeResource
    }

class MetaboliteAdductResource(Resource):
    document = d.MetaboliteAdduct
    related_resources = {
        "adduct_weight" : AdductWeightsResource
    }
