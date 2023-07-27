import numpy as np
from SDRP import SDRP
from htp import htp


def mdl(frq: np.ndarray, params: np.ndarray, aux_params: np.ndarray) -> np.ndarray:
    # parameters order:
    # 0 - center frq; 1 - strength factor; 2 - width G0; 3 - sd width G2;
    # 4 - sd shift D2; 5 - mixing Y0; 6 - Nu_vc;
    # 7 - power factor; 8 - BL0; 9 - BL1; 10 - BL2; 11 - BL3

    # aux_params order:
    # 0 - absolute Strength (for current Pself and T)
    # 1 - dopler width
    # 2 - deviation value, MHz
    # 3 - cell length
    # 4 - type (basically 0 = CAV, 1 = RAD w deviation, 2 = VID w deviation, 3 = VID natural

    absor = np.empty_like(frq)
    absor1 = np.empty_like(frq)
    if aux_params[-1] == 0:
        # for ifr in range(len(frq)):
        #    absor[ifr] = SDRP(frq[ifr], params[0], params[5], params[2], params[3], 0., params[4], aux_params[0]) \
        #                 * params[1]
        absor, absor1 = htp(params[0], 1.e-6, params[2], params[3], 0.,
                            params[4], 1.e-6, 0., frq, Ylm=params[5])
        absor *= aux_params[0] * params[1]
    absor = absor * (1 + params[7] * (frq - params[0])) \
            + params[8] + params[9] * (frq - params[0]) + params[10] * frq ** 2 \
            + params[11] * (frq - params[0]) ** 3
    return absor


def mdljac(frq: np.ndarray, jac_flag: np.ndarray, params: np.ndarray, aux_params: np.ndarray) -> np.ndarray:
    jac = np.zeros((len(params), len(frq)))
    model1 = mdl(frq, params, aux_params=aux_params)
    dpar = 1.E-6
    params_deriv_numeric = [0, 2, 3, 4, 5, 6]  # for this params index the derivative is calc numerically
    for ipar in params_deriv_numeric:
        if jac_flag[ipar] == 1.:
            params2 = np.copy(params)
            params2[ipar] = params2[ipar] + dpar
            model2 = mdl(frq, params2, aux_params=aux_params)
            jac[ipar] = (model2 - model1) * 1.E6

    # derivative by vertical scale is just profile with params[1] == 1

    if jac_flag[1] == 1.:
        params2 = np.copy(params)
        params2[1] = 1.
        jac[1] = mdl(frq, params2, aux_params=aux_params)

    # derivatives by baseline parameters etc

    if aux_params[-1] == 0:
        jac[7] = 0.

    if jac_flag[8] == 1.:
        jac[8] = 1.
    if jac_flag[9] == 1.:
        jac[9] = frq - params[0]
    if jac_flag[10] == 1.:
        jac[10] = frq ** 2
    if jac_flag[11] == 1.:
        jac[11] = (frq - params[0]) ** 3
    return jac.T
