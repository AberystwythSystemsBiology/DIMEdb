from flask_mongorest.operators import Operator

class PiPpm(Operator):
    op = "pi_ppm"
    def prepare_queryset_kwargs(self, field, value, negate):
        ionisation, mz, ppm = value.split(",")
        mz = float(mz)
        ppm = float(ppm)
        difference = abs(mz * (ppm * 0.0001))

        '''
        db.metabolites.find(
            { "adducts.positive.peaks.accurate_mass": {$gt : 1200, $lt : 1220}},
            { "name" : 1, "accurate_mass" :1, "adducts.positive.peaks.$": 1}
        )

        db.metabolites.find(
            {"adducts.negative.peaks": {$elemMatch : {"accurate_mass" : {$gt : 100.00, $lt : 105.00}}}},
            { "name" : 1, "accurate_mass" :1, "adducts.negative.peaks.$": 1}
        )
        '''
        return {
            "adducts__"+ionisation+"__peaks__accurate_mass__lt" : mz+difference,
            "adducts__" + ionisation + "__peaks__accurate_mass__gt": mz-difference,
        }

