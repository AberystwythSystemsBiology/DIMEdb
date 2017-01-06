from mongoengine import *

class MetaboliteBasic(DynamicDocument):
    meta = { "collection" : "metabolites"}
    name = StringField()
    origins = ListField(StringField())
    molecular_formula = StringField()
    accurate_mass = FloatField()

