from mongoengine import *


class Adduct(EmbeddedDocument):
    type = StringField()
    accurate_mass = FloatField()
    isotopic_distribution = ListField(ListField(FloatField()))

class Adducts(EmbeddedDocument):
    neutral = EmbeddedDocumentListField(Adduct)
    positive = EmbeddedDocumentListField(Adduct)
    negative = EmbeddedDocumentListField(Adduct)

class MetaboliteFull(DynamicDocument):
    meta = {"collection" : "metabolites"}
    id = StringField(primary_key=True)
    name = StringField()
    synonyms = ListField(StringField())
    origins = ListField(StringField())
    molecular_formula = StringField()
    smiles = StringField()
    inchi = StringField()
    accurate_mass = FloatField()
    num_atoms = IntField()
    sources = StringField()
    biofluid_locations = ListField(StringField())
    tissue_locations = ListField(StringField())
    pathways = ListField(StringField())
    adducts = EmbeddedDocumentField(Adducts)
