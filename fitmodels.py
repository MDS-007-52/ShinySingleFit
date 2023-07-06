import numpy as np
from SDRP import SDRP


def mdl(frq: np.ndarray, params: np.ndarray) -> np.ndarray:
    absor = np.empty_like(frq)
    for ifr in range(len(frq)):
        absor[ifr] = SDRP(frq[ifr], params[0], params[4], params[1], params[2], 0., params[3], params[5])
    absor = absor + params[7] + params[8] * (frq - params[0]) + params[9] * frq ** 2
    return absor


def mdljac(frq: np.ndarray, jac_flag: np.ndarray, params: np.ndarray) -> np.ndarray:
    jac = np.zeros((len(params), len(frq)))
    model1 = mdl(frq, params)
    dpar = 1.E-6
    for ipar in range(0, 5):
        if jac_flag[ipar] == 1.:
            params2 = np.copy(params)
            params2[ipar] = params2[ipar] + dpar
            model2 = mdl(frq, params2)
            jac[ipar] = (model2 - model1) * 1.E6

    if jac_flag[5] == 1.:
        params2 = np.copy(params)
        params2[5] = 1.
        jac[5] = mdl(frq, params2)

    jac[6] = 0.

    if jac_flag[7] == 1.:
        jac[7] = 1.
    if jac_flag[8] == 1.:
        jac[8] = frq - params[0]
    if jac_flag[9] == 1.:
        jac[9] = frq ** 2
    return jac.T
