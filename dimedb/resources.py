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
        "molecular_formula" : [ops.Exact, ops.Contains, ops.IContains]    }

class NegativePeaksResource(Resource):
    document = d.NegativePeaks

class NegativeAdductResource(Resource):
    document = d.NegativeAdducts

    related_resources = {
        "peaks" : NegativePeaksResource
    }

class NeutralPeaksResource(Resource):
    document = d.NeutralPeaks

class NeutralAdductResource(Resource):
    document = d.NeutralAdducts

    related_resources = {
        "peaks" : NeutralPeaksResource
    }

class PositivePeaksResource(Resource):
    document = d.PositivePeaks

class PositiveAdductResource(Resource):
    document = d.PositiveAdducts

    related_resources = {
        "peaks" : PositivePeaksResource
    }

'''
class SourcesResource(Resource):
    documents = d.Sources
'''

class AdductsResource(Resource):
    document = d.Adducts

    related_resources = {
        "positive" : PositiveAdductResource,
        "negative" : NegativeAdductResource,
        "neutral" : NeutralAdductResource
    }

class MetaboliteAdductsResource(Resource):
    document = d.MetaboliteAdducts

    filters = {
        "adducts": [o.PiPpm],
        "origins" : [ops.Contains, ops.Exact],
        "biofluids" : [ops.Contains]
    }
    max_limit, default_limit = [200, 200]

    related_resources = {
        "adducts" : AdductsResource
    }

    def get_queryset(self, *args, **kwargs):
        return super(MetaboliteAdductsResource, self).get_queryset(*args, **kwargs).fields(name=1, accurate_mass=1,
                                                                                           molecular_formula=1,
                                                                                           adducts = 1,
                                                                                           adducts_negative_peaks__S = 1
                                                                                           )
