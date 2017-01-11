from flask_mongorest.operators import Operator

class Ionisation(Operator):
    op = "ion"
    def prepare_queryset_kwargs(self, field, value, negate):
        if value == None:
            value = "positive"
        print {field+"__"+value+"__count__gt" : 0}
        return {
            field+"__"+value+"__count__gt" : 0
        }

class AccurateMassSearch(Operator):
    op = "ppm"
    def prepare_queryset_kwargs(self, field, value, negate=False):
        print field
        if value == None:
            value = [0,0]
        else:
            value = [float(x) for x in value.split(',')]

        mz, ppm_threshold = value
        difference = abs(mz * (ppm_threshold * 0.0001))  # PPM to diff.
        if "adduct_weights__positive" or "adduct_weights__negative" in field:
            return {
                field + "__peaks__peak__gt" : mz - difference,
                field + "__peaks__peak__lt" : mz + difference
            }
        else:
            return {
                field + '__gt': mz-difference,
                field + '__lt': mz+difference
            }