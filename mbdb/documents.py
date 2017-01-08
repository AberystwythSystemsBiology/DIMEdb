from mongoengine import *

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()

class AdductWeights(EmbeddedDocument):
    neutral = FloatField()

class MetaboliteAdduct(Document):
    meta = {" collection" : "metabolites"}
    name = StringField()
    adduct_weights = EmbeddedDocumentField(AdductWeights)
