from flask_mongorest.operators import Operator

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