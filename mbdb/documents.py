from mongoengine import *

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()

class PositiveAdduct(EmbeddedDocument):
    count = IntField()
    peaks = StringField()

class NegativeAdduct(EmbeddedDocument):
    count = IntField()
    peaks = StringField()

class AdductWeights(EmbeddedDocument):
    neutral = FloatField()
    positive = EmbeddedDocumentField(PositiveAdduct)
    negative = EmbeddedDocumentField(NegativeAdduct)

class MetaboliteAdduct(DynamicDocument):
    meta = {"collection" : "metabolites"}
    name = StringField()
    adduct_weights = EmbeddedDocumentField(AdductWeights)
