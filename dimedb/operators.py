from flask_mongorest.operators import Operator
from documents import *


class AdductPpm(Operator):
    op = "ppm"

    def calculate_ppm(self, mz, ppm):
        tolerance = abs(float(mz) * (float(ppm) * 0.0001))
        return float(mz) - tolerance, float(mz) + tolerance

    def generate_types(self, ionisation, type_list):

        adduct_dict = {
            "negative": {
                0: "[M-H]1-",
                1: "[M+Cl]1-",
                2: "[M+Na-2H]1-"
            },
            "positive": {
                0: "[M+H]1+",
                1: "[M+Na]1+",
                2: "[M+K]1+"
            }
        }

        return [adduct_dict[ionisation][int(x)] for x in type_list]

    def apply(self, queryset, field, value, negate=False):
        values = value.split(",")

        if len(values) == 3:
            ionisation, mz, ppm = values
            types = None
        elif len(values) > 3 and values[0] != "neutral":
            ionisation, mz, ppm = values[:3]
            types = self.generate_types(ionisation, values[3:])

        l_mz, g_mz = self.calculate_ppm(mz, ppm)

        if types == None:
            collection = MetaboliteFull._get_collection().find(
                {"adducts."+str(ionisation): {"$elemMatch": {"accurate_mass": {"$gt": l_mz, "$lt": g_mz}}}},
                {"name" : 1, "adducts."+str(ionisation)+".$": 1}
            )
        else:
            collection = MetaboliteFull._get_collection().find(
                {"adducts." + str(ionisation): {"$elemMatch": {"accurate_mass": {"$gt": l_mz, "$lt": g_mz},
                                                               "type" : {"$in" : types}}}},
                {"name": 1, "adducts." + str(ionisation) + ".$": 1}
            )


        queryset = MetaboliteFull.objects()
        queryset._cursor_obj = collection

        return queryset
