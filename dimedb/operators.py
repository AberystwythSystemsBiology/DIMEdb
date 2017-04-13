from flask_mongorest.operators import Operator

'''
    Create an dictonary containing the url-suitable URLs?
'''

class AdductPpm(Operator):
    op = "ppm"

    ionisation = None
    mz = None
    ppm = None
    type = None

    def calculate_ppm(self, mz, ppm):
        return abs(float(mz) * (float(ppm) * 0.0001))


    def prepare_queryset_kwargs(self, field, value, negate):
        ionisation, mz, ppm, type = value.split(",")

        self.ionisation = ionisation
        self.mz = mz
        self.ppm = ppm
        self.type = type


        difference = self.calculate_ppm(mz, ppm)

        return {
            "adducts__"+ionisation+"__match": {"accurate_mass__gt": float(mz) - difference,
                "accurate_mass__lt": float(mz) + difference, "type" : "[M-H]1-"},

        }

    def apply(self, queryset, field, value, negate=False):
        kwargs = self.prepare_queryset_kwargs(field, value, negate)
        return queryset.filter(**kwargs).fields(name=1, accruate_mass=1, adducts= 1).only("adducts."+self.ionisation)
