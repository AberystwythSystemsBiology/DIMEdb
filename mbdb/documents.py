from mongoengine import *

class MetaboliteFull(DynamicDocument):
    meta = {"collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    smiles = StringField()
    accurate_mass = FloatField()
    num_atoms = IntField()
    sources = StringField()
    pathways = StringField()
    synonyms = ListField(StringField())
    adducts = StringField()

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()


class NegativePeaks(EmbeddedDocument):
    type = StringField()
    accurate_mass = FloatField()
    isotopic_distribution = ListField()

class NegativeAdducts(EmbeddedDocument):
    count = IntField()
    peaks = EmbeddedDocumentListField(NegativePeaks)

class PositivePeaks(EmbeddedDocument):
    type = StringField()
    accurate_mass = FloatField()
    isotopic_distribution = ListField()

class PositiveAdducts(EmbeddedDocument):
    count = IntField()
    peaks = EmbeddedDocumentListField(PositivePeaks)

class NeutralPeaks(EmbeddedDocument):
    type = StringField()
    accurate_mass = FloatField()
    isotopic_distribution = ListField()

'''
class Sources(EmbeddedDocument):
    kegg_id = StringField()
    chebi_id = StringField()
    pubchem_id = StringField()
'''

class NeutralAdducts(EmbeddedDocument):
    count = IntField()
    peaks = EmbeddedDocumentListField(PositivePeaks)

class Adducts(EmbeddedDocument):
    neutral = EmbeddedDocumentField(NeutralAdducts)
    positive = EmbeddedDocumentField(PositiveAdducts)
    negative = EmbeddedDocumentField(NegativeAdducts)

class MetaboliteAdducts(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    accurate_mass = FloatField()
    molecular_formula = StringField()
    adducts = EmbeddedDocumentField(Adducts)
    #sources = EmbeddedDocumentField(Sources)