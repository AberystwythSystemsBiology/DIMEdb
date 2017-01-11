import documents as d, operators as o
from flask_mongorest import operators as ops
from flask_mongorest.resources import Resource

class MetaboliteFullResource(Resource):
    document = d.MetaboliteFull
    max_limit, default_limit = [1, 1]
    filters = {
        "id" : [ops.Exact]
    }

class MetaboliteBasicResource(Resource):
    document = d.MetaboliteBasic
    max_limit, default_limit = [1000, 1000]

    filters = {
        "name" : [ops.Exact, ops.Startswith, ops.Contains, ops.IContains],
        "origins" : [ops.Exact, ops.In(allow_negation=True)],
        "molecular_formula" : [ops.Exact, ops.Contains, ops.IContains],
        "accurate_mass" : [ops.Exact, o.AccurateMassSearch, ops.Gte, ops.Gt]
    }

class NegativePeaksResource(Resource):
    document = d.NegativePeaks

class PositivePeaksResource(Resource):
    document = d.PositivePeaks

class NegativeAdductResource(Resource):
    document = d.NegativeAdduct
    related_resources = {
        "peaks" : NegativePeaksResource
    }

class PositiveAdductResource(Resource):
    document = d.PositiveAdduct
    related_resources = {
        "peaks" : PositivePeaksResource
    }

class AdductWeightsResource(Resource):
    document = d.AdductWeights

    related_resources = {
        "positive" : PositiveAdductResource,
        "negative" : NegativeAdductResource
    }

class MetaboliteAdductResource(Resource):
    document = d.MetaboliteAdduct
    filters = {
        "adduct_weights__positive__count": [ops.Gt],
        "adduct_weights": [o.Ionisation, o.IonisationPpm],
        "adduct_weights__negative__peaks" : [ops.Contains],
        "adduct_weights__neutral": [o.AccurateMassSearch]
    }

    related_resources = {
        "adduct_weights" : AdductWeightsResource
    }
