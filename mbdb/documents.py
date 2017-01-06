from mongoengine import *

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()

class Positive(EmbeddedDocument):
    count = IntField()
    peaks = StringField() # Temporary

class Negative(EmbeddedDocument):
    count = IntField()
    peaks = StringField() # Temporary

class AdductWeights(EmbeddedDocument):
    negative = EmbeddedDocumentField(Negative)
    positive = EmbeddedDocumentField(Positive)
    neutral = FloatField()

class MetaboliteAdduct(Document):
    meta = {"collection" : "metabolites"}
    name = StringField()
    adduct_weights = EmbeddedDocumentField(AdductWeights)