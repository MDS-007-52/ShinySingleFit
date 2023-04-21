# this is integrated quadratic speed-dependent VVW profile with mixing
# also referred to as qSDVVWLM or qSDRP (quadratic Rosenkranz profile)
# used for elevated pressure recordings (from 100 Torr or where mixing effect is pronounced)
# there is no Dopler or Dicke features as they are unnecessary at these pressures

import math
import scipy.integrate as qdrtr

def FMB(vel: float) -> float:
    """
    Maxwell-Boltzmann distribution versus reduced velocity
    v = V [m/s] / V0 [m/s]
    V0 = sqrt(2 * k_B * T / m)
    :param vel: reduced velocity, unitless
    :return: Maxwell-Boltzmann distribution, unitless
    """
    return 2.256758E0 * vel ** 2 * math.exp(-vel ** 2)


def VVWLM(x: float, frq: float, f0: float, y0: float, g0: float, g2: float, d0: float, d2: float) -> float:
    """
    Van Vleck - Weisskopf profile with line mixing and square speed-dependence of the collisional parameters
    :param x: current reduced velocity (relation of V to V_most_probable), unitless
    :param frq: current frequency, frq units (should be the same for all frequency-based parameters)
    :param f0: central frequency, frq units
    :param y0: 1st order mixing parameter, unitless
    :param g0: Speed-independent part of HWHM, frq units
    :param g2: Speed-dependent part of HWHM, frq units
    :param d0: Speed-independent part of collisional central f0 shift, frq units
    :param d2: Speed-dependent part of collisional central f0 shift, frq units
    :return: spectral function value normalized to area multiplied to the Maxwell-Boltzmann distribution, arbitrary units
    NB: (f/f0)^2 and 1/pi factors are not included!
    """
    gf = g0 + g2 * (x ** 2 - 1.5)  # SD-HWHM at current velocity x
    df = d0 + d2 * (x ** 2 - 1.5)  # SD-shift at current velocity x
    difp = frq - f0 - df
    difm = frq + f0 + df
    VV = (gf + y0 * difp) / (gf ** 2 + difp ** 2) + (gf - y0 * difm) / (gf ** 2 + difm ** 2)
    return VV * FMB(x)


# Speed-dependent Rosenkranz profile
def SDRP(frq: float, f0: float, y0: float, g0: float, g2: float, d0: float, d2: float, strength: float) -> float:
    """
    Returns an absorption coef. calculated using qSDRP (quadratic Speed-Dep-t Rosenkranz profile)
    :param frq: current frequency, frq units (see above)
    :param f0: central frequency of the spectral line, frq units
    :param y0: - 1st-order mixing parameter, unitless
    :param g0: SI-HWHM part, frq units
    :param g2: SD-HWHM part, frq units (setting to 0 turns profile to Rosenkranz one)
    :param d0: Speed-independent part of collisional central f0 shift, frq units
    :param d2: Speed-dependent part of collisional central f0 shift, frq units
    :param strength: integrated intensity of the line profile, should be set out of the function (see https://hitran.org/docs/definitions-and-units/ for the information), 1/cm * frq units / m^-3
    :return : absorption coefficient, units depend on frq units and strength calculation, 1/cm
    """

    (sdshape, tmp_precision) = qdrtr.quad(VVWLM, 0., 99., args=(frq, f0, y0, g0, g2, d0, d2))
    sdshape = sdshape * (frq / f0) ** 2 * strength / 3.1415926
    return sdshape
