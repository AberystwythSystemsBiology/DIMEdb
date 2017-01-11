from flask_mongorest.operators import Operator

class PiPpm(Operator):
    op = "pi_ppm"
    def prepare_queryset_kwargs(self, field, value, negate):
        ionisation, mz, ppm = value.split(",")
        mz = float(mz)
        ppm = float(ppm)
        difference = abs(mz * (ppm * 0.0001))
        return {
            "adducts__"+ionisation+"__peaks__accurate_mass__lt" : mz+difference,
            "adducts__" + ionisation + "__peaks__accurate_mass__gt": mz-difference
        }