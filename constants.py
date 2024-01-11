# global kB, clight, pi, k_vs_aem, k_st, c2, hP, m_aem, t_ref

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

multi_params_dict = {"int": "Intensity",
                     "frq": "Central freq",                      
                     "g0s": "Gamma 0 self",
                     "g2s": "Gamma 2 self",
                     "g0f": "Gamma 0 foreign",
                     "g2f": "Gamma 2 foreign",                     
                     "d0s": "Delta 0 self",
                     "d2s": "Delta 2 self",
                     "d0f": "Delta 0 foreign",
                     "d2f": "Delta 2 foreign",
                     "y0s": "Mixing self",
                     "y0f": "Mixing foreign",
                     "nuvcs": "Vel.chng.rate self",
                     "nuvcf": "Vel.chng.rate foreign",
                     "cs" : "Continuum self",
                     "cf" : "Continuum foreign",
                     "frab": "Rabi frequency",
                     "ng0s": "T-dependence parameter g0s",
                     "ng2s": "T-dependence parameter g2s",
                     "ng0f": "T-dependence parameter g0f",
                     "ng2f": "T-dependence parameter g2f",
                     "nd0s": "T-dependence parameter d0s",
                     "nd2s": "T-dependence parameter d2s",
                     "nd0f": "T-dependence parameter d0f",
                     "nd2f": "T-dependence parameter d2f",
                     "ny0s": "T-dependence parameter y0s",
                     "ny0f": "T-dependence parameter y0f",
                     "nnuvcs": "T-dependence parameter nu_vc self",
                     "nnuvcf": "T-dependence parameter nu_vc foreign",
                     "ncs": "T-dependence parameter Cs",
                     "ncf": "T-dependence parameter Cf",
                     "scl": "Intensity scale",
                     "pow": "Power factor",
                     "bl0": "Baseline 0",
                     "bl1": "Baseline 1",
                     "bl2": "Baseline 2",
                     "bl3": "Baseline 3"}
