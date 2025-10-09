# global kB, clight, pi, k_vs_aem, k_st, c2, hP, m_aem, t_ref
import numpy as np
from typing import Union
from scipy.optimize import curve_fit

# Various functions for T- and P-dependencies

# Single exponent law for T-dependence of broadening coefs
def g0vst(x: float, g0: float, ng: float) -> float:
    return g0 * (t_ref / x)**ng


# Foreign broadening width G_V or G_0 versus pressure with saturation (Rabi frq parameter is "rf")
def g0vsp_frgn(x: float, a: float, b: float, rf: float) -> float:
    return np.sqrt((a + b * x)**2 + rf**2)


# Self broadening width G_V or G_0 versus pressure with saturation (Rabi frq parameter is "rf")
def g0vsp_self(x: float, b: float, rf: float) -> float:
    return np.sqrt((b * x)**2 + rf**2)


# Linear P-dependence with ZERO intercept (for self G_2 or D_2 versus P)
def g0vsp_0(x: float, a: float) -> float:
    return x * a


# Fit routines returning params and corresponding errors

# Fit of the regular linear dependence with slope l_coefs[0] and intercept l_coefs[1] (suits for central frq0 vs P)
def linear_fit(data_x: Union[list, np.ndarray], data_y: Union[list, np.ndarray], data_yerr: Union[list, np.ndarray]) -> list[np.ndarray, np.ndarray]:
    l_coefs, l_cov = np.polyfit(np.copy(data_x)[:], np.copy(data_y)[:], 1, rcond=None, full=False, w = 1. / np.copy(data_yerr)[:]**2, cov=True)
    l_err = [np.sqrt(l_cov[i, i]) for i in range(len(l_coefs))]
    return l_coefs, l_err

# Fit of the linear P-dependence with zero intercept
def linear_fit_0(data_x: Union[list, np.ndarray], data_y: Union[list, np.ndarray], data_yerr: Union[list, np.ndarray]) -> list[np.ndarray, np.ndarray]:
    # coef =  (max(data_y)/max(data_x))
    params, cov = curve_fit(g0vsp_0, np.copy(data_x), np.copy(data_y), sigma=np.copy(data_yerr), absolute_sigma=True, full_output=False)
    return params[0], np.sqrt(cov[0, 0])

# Fit of the saturated P-dependence of the SI-part of collisional width for self-broadening
# order: 0 - broadening, 1 - Rabi frq
def rabi_fit_self(data_x: Union[list, np.ndarray], data_y: Union[list, np.ndarray], data_yerr: Union[list, np.ndarray], rabi0=0.1) -> list[np.ndarray, np.ndarray]:
    slope = np.max(data_y) / np.max(data_x)
    rabi = rabi0
    params0 = np.array([slope, rabi])
    params, cov = curve_fit(g0vsp_self, np.copy(data_x), np.copy(data_y), sigma=np.copy(data_yerr), p0=params0, absolute_sigma=True, full_output=False)
    err = [np.sqrt(cov[i, i]) for i in range(len(params))]
    return params, err

# Fit of the saturated P-dependence of the SI-part of collisional width for foreign-broadening
# order: 0 - intercept, 1 - broadening, 2 - Rabi frq
def rabi_fit_frgn(data_x: Union[list, np.ndarray], data_y: Union[list, np.ndarray], data_yerr: Union[list, np.ndarray]) -> list[np.ndarray, np.ndarray]:
    slope = (np.max(data_y) - np.min(data_y)) / np.max(data_x)
    intercept = np.min(data_y)
    rabi = 0.1
    params0 = np.array([intercept, slope, rabi])
    params, cov = curve_fit(g0vsp_frgn, np.copy(data_x), np.copy(data_y), sigma=np.copy(data_yerr), p0=params0, absolute_sigma=True, full_output=False)
    err = [np.sqrt(cov[i, i]) for i in range(len(params))]
    return params, err

hP = 6.62607004E-34
m_aem = 1.6605402E-27
kB = 1.38064852E-23
k_st = 2894.947919E20  # recalc coef-t for the intensity
pi = 3.14159265359E0
clight = 2.997925E08  # speed of light, m/s
k_vs_aem = 0.120272478907E-3  # atomic mass unit [kg] divided by k_B [J/K]
c2 = 1.4388  # cm*K c2 = hc/k_B for intensity calcs
t_ref = 296.  # reference temperature, K

mhz_in_cm = 29979.25  # recalc factor from 1/cm to MHz and back

npar = 12  # number of common parameters for all the recordings in batch
nauxpar = 5  # number of auxilary params

# n_add_par = 2  # number of individual parameters of a separate recording

n_const_par = 31  # number of common parameters for all recordings in the multifit
n_add_par = 6  # number of individual parameters for each recording in multifit

# f_verbose_fit = True  # this feature was moved to the interface and subr_fit function

single_params_dict = {"frq": "Central freq",
                      "int": "Intensity",
                      "g0": "Gamma 0",
                      "g2": "Gamma 2",
                      "d2": "Delta 2",
                      "y": "Mixing",
                      "nuvc": "Vel.chng.rate",
                      "pow": "Power factor",
                      "bl0": "Baseline 0",
                      "bl1": "Baseline 1",
                      "bl2": "Baseline 2",
                      "bl3": "Baseline 3"}

single_params_adjustable_init = ["frq",
                                 "int",
                                 "g0",
                                 "g2",
                                 "y",
                                 "bl0",
                                 "bl1",
                                 "bl2"]


single_params_preset_si = ["frq",
                            "int",
                            "g0",
                            "pow",                                                        
                            "bl0",
                            ]

single_params_preset_sd = ["frq",
                            "int",
                            "g0",
                            "g2",
                            "d2",
                            "pow",                                                        
                            "bl0",
                            "bl1",
                            "bl2"]


single_params_adjustable_cav = ["frq",
                                 "int",
                                 "g0",
                                 "g2",
                                 "y",                                 
                                 "bl2"]

# preset of the adjustable parameters for LBL fit corresp. to RAD (and VID) recordings
# no mixing, no SD-shift, no Opt.Frq.
single_params_adjustable_rad = ["frq",
                                 "int",
                                 "g0",
                                 "g2",
                                 "pow",                                 
                                 "bl0",
                                 "bl1",
                                 "bl2"]


multi_params_dict = {"mint": "Intensity",
                     "mf0": "Central freq",                      
                     "mg0s": "Gamma 0 self",
                     "mg2s": "Gamma 2 self",
                     "mg0f": "Gamma 0 foreign",
                     "mg2f": "Gamma 2 foreign",                     
                     "md0s": "Delta 0 self",
                     "md2s": "Delta 2 self",
                     "md0f": "Delta 0 foreign",
                     "md2f": "Delta 2 foreign",
                     "my0s": "Mixing self",
                     "my0f": "Mixing foreign",
                     "mnuvcs": "Vel.chng.rate self",
                     "mnuvcf": "Vel.chng.rate foreign",
                     "mcs" : "Continuum self",
                     "mcf" : "Continuum foreign",
                     "mfrab": "Rabi frequency",
                     "mng0s": "T-dependence parameter g0s",
                     "mng2s": "T-dependence parameter g2s",
                     "mng0f": "T-dependence parameter g0f",
                     "mng2f": "T-dependence parameter g2f",
                     "mnd0s": "T-dependence parameter d0s",
                     "mnd2s": "T-dependence parameter d2s",
                     "mnd0f": "T-dependence parameter d0f",
                     "mnd2f": "T-dependence parameter d2f",
                     "mny0s": "T-dependence parameter y0s",
                     "mny0f": "T-dependence parameter y0f",
                     "mnnuvcs": "T-dependence parameter nu_vc self",
                     "mnnuvcf": "T-dependence parameter nu_vc foreign",
                     "mncs": "T-dependence parameter Cs",
                     "mncf": "T-dependence parameter Cf",
                     "mscl": "Intensity scale",
                     "mpow": "Power factor",
                     "mbl0": "Baseline 0",
                     "mbl1": "Baseline 1",
                     "mbl2": "Baseline 2",
                     "mbl3": "Baseline 3"}

multi_params_indx = {"mint": 0,
                     "mf0": 1,
                     "mg0s": 2,
                     "mg2s": 3,
                     "mg0f": 4,
                     "mg2f": 5,                     
                     "md0s": 6,
                     "md2s": 7,
                     "md0f": 8,
                     "md2f": 9,
                     "my0s": 10,
                     "my0f": 11,
                     "mnuvcs": 12,
                     "mnuvcf": 13,
                     "mcs" : 14,
                     "mcf" : 15,
                     "mfrab": 16,
                     "mng0s": 17,
                     "mng2s": 18,
                     "mng0f": 19,
                     "mng2f": 20,
                     "mnd0s": 21,
                     "mnd2s": 22,
                     "mnd0f": 23,
                     "mnd2f": 24,
                     "mny0s": 25,
                     "mny0f": 26,
                     "mnnuvcs": 27,
                     "mnnuvcf": 28,
                     "mncs": 29,
                     "mncf": 30,
                     "mscl": 31,
                     "mpow": 32,
                     "mbl0": 33,
                     "mbl1": 34,
                     "mbl2": 35,
                     "mbl3": 36}

multi_params_numeric_deriv = ("mf0",
                              "mg0s",
                              "mg2s",
                              "mg0f",
                              "mg2f",                     
                              "md0s",
                              "md2s",
                              "md0f",
                              "md2f",
                              "my0s",
                              "my0f",
                              "mnuvcs",
                              "mnuvcf",
                              "mfrab",
                              "mng0s",
                              "mng2s",
                              "mng0f",
                              "mng2f",
                              "mnd0s",
                              "mnd2s",
                              "mnd0f",
                              "mnd2f",
                              "mny0s",
                              "mny0f",
                              "mnnuvcs",
                              "mnnuvcf")

multi_params_adjustable_init = ["mf0" ,
                                "mg0s" ,
                                "mg2s" ,
                                "mg0f" ,
                                "mg2f" ,
                                "md0s" ,
                                "md0f" ,
                                "my0s" ,
                                "my0f" ,
                                "mcs"  ,
                                "mcf"  ,
                                "mfrab" ,
                                "mscl" ,
                                "mbl0" ,
                                "mbl1" ,
                                "mbl2" ,
                                "mbl3"]

multi_params_preset_sis = ["mf0",
                            "mg0s",                                                        
                            "md0s",                            
                            "mscl",
                            "mpow",
                            "mbl0"]

multi_params_preset_sif = ["mg0s",                            
                            "mg0f",                            
                            "md0s",
                            "md0f",
                            "mscl",
                            "mpow",
                            "mbl0"]

multi_params_preset_sds = ["mf0",
                            "mg0s",
                            "mg2s",
                            "md0s",
                            "md2s",
                            "mscl",
                            "mpow",
                            "mbl0", 
                            "mbl1",
                            "mbl2"]

multi_params_preset_sdf = ["mg0s",
                           "mg2s",
                            "mg0f",
                            "mg2f",
                            "md0s",
                            "md2s",
                            "md0f",
                            "md2f",
                            "mscl",
                            "mpow",
                            "mbl0", 
                            "mbl1",
                            "mbl2"]