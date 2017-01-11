from flask_mongorest.operators import Operator

class Ionisation(Operator):
    op = "ion"
    def prepare_queryset_kwargs(self, field, value, negate):
        if value == None:
            value = "positive"
        return {
            field + "__"+value+"__count__gt" : 0
        }

class AccurateMassSearch(Operator):
    op = "ppm"
    def prepare_queryset_kwargs(self, field, value, negate=False):
        if value == None:
            value = [0,0]
        else:
            value = [float(x) for x in value.split(',')]

        mz, ppm_threshold = value

        difference = abs(mz * (ppm_threshold * 0.0001))  # PPM to diff.
        return {
            field + '__gt': mz-difference,
            field + '__lt': mz+difference
        }

class IonisationPpm(Operator):
    op = "ppm"
    def prepare_queryset_kwargs(self, field, value, negate=False):
        if value == None:
            value = ["negative", 100, 10]
        else:
            value = [x for x in value.split(",")]

        ionisation = value[0]
        mz, ppm_threshold = [float(x) for x in value[1:]]
        difference = abs(mz * (ppm_threshold * 0.0001))  # PPM to diff.

        print {
            field + "__" + ionisation + "__count__gt": 0,
            field + "__" + ionisation + "__peaks__1__gt": mz - difference,
            field + "__" + ionisation + "__peaks__1__lt": mz + difference
        }
        return {
            field + "__" + ionisation + "__count__gt": 0,
            field + "__" + ionisation + "__peaks__0__1__gt": mz - difference,
            field + "__" + ionisation + "__peaks__0__1__lt": mz + difference
        }
