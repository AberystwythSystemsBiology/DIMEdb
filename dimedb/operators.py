from flask_mongorest.operators import Operator
import documents

'''
    Create an dictonary containing the url-suitable Adducts?
'''

class AdductPpm(Operator):
    op = "ppm"

    ionisation = None
    mz = None
    ppm = None
    type = None

    def calculate_ppm(self, mz, ppm):
        return abs(float(mz) * (float(ppm) * 0.0001))

    def apply(self, queryset, field, value, negate=False):
        kwargs = self.prepare_queryset_kwargs(field, value, negate)

        ionisation, mz, ppm = value.split(",")

        difference = self.calculate_ppm(mz, ppm)
        print difference
        return documents.MetaboliteFull.objects(__raw__=(
            {'adducts.'+ionisation: {
                '$elemMatch': {
                    'accurate_mass': {'$gt': float(mz) - difference,
                                      '$lt': float(mz) + difference}
                }
            }}
        )).only(
            "name",
            "adducts."+ionisation+".accurate_mass",
            "adducts."+ionisation+".type")
