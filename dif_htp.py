import numpy as np
from htp import htp


def dif_htp(sg, sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, strength_in, dev, scale, pow, tpl):
    if type(sg) not in set([np.array, np.ndarray, list, tuple]):
        sg = np.asarray([sg])

    sg_m = sg[:] - dev
    sg_p = sg[:] + dev
    absor_minus, absor_im = htp(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, 0.0, sg - dev, Ylm=0.0)
    absor_plus, absor_im = htp(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, 0.0, sg + dev, Ylm=0.0)

    if tpl == 1:
        absor_minus[:] = scale * (1.0 - np.exp(-absor_minus[:] * strength_in)) * (1 + pow * (sg_m[:] - sg0))
        absor_plus[:] = scale * (1.0 - np.exp(-absor_plus[:] * strength_in)) * (1 + pow * (sg_p[:] - sg0))

    elif tpl == 2:
        absor_minus[:] = np.exp(-absor_minus[:] * strength_in) * (scale + pow * (sg_m[:] - sg0))
        absor_plus[:] = np.exp(-absor_plus[:] * strength_in) * (scale + pow * (sg_p[:] - sg0))

    return absor_plus - absor_minus
