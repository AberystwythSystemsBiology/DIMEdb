from mongoengine import *

class MetaboliteFull(DynamicDocument):
    meta = {"collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    smiles = StringField()
    accurate_mass = FloatField()
    adduct_weights = StringField()
    isotopic_distributions = StringField()

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()

class PositiveAdduct(EmbeddedDocument):
    count = IntField()
    try:
        peaks = ListField(ListField(StringField()))
    except TypeError:
        peaks = StringField()

class NegativeAdduct(EmbeddedDocument):
    count = IntField()
    try:
        peaks = ListField(ListField(StringField()))
    except TypeError:
        peaks = StringField()

class AdductWeights(EmbeddedDocument):
    neutral = FloatField()
    positive = EmbeddedDocumentField(PositiveAdduct)
    negative = EmbeddedDocumentField(NegativeAdduct)

class MetaboliteAdduct(DynamicDocument):
    meta = {"collection" : "metabolites"}
    name = StringField()
    adduct_weights = EmbeddedDocumentField(AdductWeights)
