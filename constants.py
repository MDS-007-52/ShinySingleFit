#global kB, clight, pi, k_vs_aem, k_st, c2, hP, m_aem, t_ref

hP = 6.62607004E-34
m_aem = 1.6605402E-27
kB = 1.38064852E-23
k_st = 2894.947919E20  # recalc coef-t for the intensity
pi = 3.14159265359E0
clight = 2.997925E08  # speed of light, m/s
k_vs_aem = 0.120272478907E-3  # atomic mass unit [kg] divided by k_B [J/K]
c2 = 1.4388  # cm*K c2 = hc/k_B for intensity calcs
t_ref = 296.  # reference temperature, K

n_const_par = 2  # number of common parameters for all the recordings in batch

n_add_par = 2  # number of individual parameters of a separate recording

single_params_dict = {"frq":"Central freq",
                      "int":"Intensity",
                      "g0":"Gamma 0",
                      "g2":"Gamma 2",
                      "d2":"Delta 2",
                      "y":"Mixing",
                      "nuvc":"Vel.chng.rate",
                      "pow":"Power factor",
                      "bl0":"Baseline 0",
                      "bl1":"Baseline 1",
                      "bl2":"Baseline 2",
                      "bl3":"Baseline 3"}