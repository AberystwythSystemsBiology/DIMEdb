import documents as d, operators as o
from flask_mongorest import operators as ops
from flask_mongorest.resources import Resource

class MetaboliteBasicResource(Resource):
    document = d.MetaboliteBasic
    filters = {
        "id" : [ops.Exact],
        "name" : [ops.Exact, ops.Startswith, ops.Contains],
        "origins" : [ops.Exact, ops.In(allow_negation=True)],
        "molecular_formula" : [ops.Exact, ops.Contains, ops.IContains],
        "accurate_mass" : [ops.Exact, o.AccurateMassSearch, ops.Gte, ops.Gt]
    }

class NegativeAdductResource(Resource):
    document = d.NegativeAdduct

class PositiveAdductResource(Resource):
    document = d.PositiveAdduct

class AdductWeightsResource(Resource):
    document = d.AdductWeights

    related_resources = {
        "positive" : PositiveAdductResource,
        "negative" : NegativeAdductResource
    }

class MetaboliteAdductResource(Resource):
    document = d.MetaboliteAdduct
    filters = {
        "name" : [ops.Exact]
    }

    related_resources = {
        "adduct_weights" : AdductWeightsResource
    }
