import numpy as np
# from SDRP import SDRP
from htp import htp
from dif_htp import dif_htp
from constants import *


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
        # old code for using hones integration form of the SDRP profile for resonator recordings at high pressure
        # for ifr in range(len(frq)):
        #    absor[ifr] = SDRP(frq[ifr], params[0], params[5], params[2], params[3], 0., params[4], aux_params[0]) \
        #                 * params[1]
        absor, absor1 = htp(params[0], aux_params[1], params[2], params[3], 0.,
                            params[4], params[6], 0., frq, Ylm=params[5])
        absor[:] *= aux_params[0] * params[1] * (frq[:]/params[0])**2
    if aux_params[-1] in [1, 2]:
        absor = dif_htp(frq, params[0], aux_params[1], params[2], params[3], 0.,
                        params[4], params[6], aux_params[0], aux_params[2], params[1], params[7], aux_params[-1])
    absor = absor + params[8] + params[9] * (frq - params[0]) \
            + params[11] * (frq - params[0]) ** 3
    if aux_params[-1] == 0:
        absor += params[10] * frq ** 2
    if aux_params[-1] in [1, 2]:
        absor += params[10] * (frq - params[0])**2
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
            jac[ipar] = (model2 - model1) / dpar

    # derivative by vertical scale is just profile with params[1] == 1

    if jac_flag[1] == 1.:
        params2 = np.copy(params)
        params2[1] = 1.
        jac[1] = mdl(frq, params2, aux_params=aux_params)

    # derivatives by baseline parameters etc

    if aux_params[-1] == 0:
        jac[7] = 0.
    else:
        if jac_flag[7] == 1:
            params2 = np.copy(params)
            params2[7] = params2[7] + dpar
            model2 = mdl(frq, params2, aux_params=aux_params)
            jac[7] = (model2 - model1) / dpar

    if jac_flag[8] == 1.:
        jac[8] = 1.
    if jac_flag[9] == 1.:
        jac[9] = frq - params[0]
    if jac_flag[10] == 1.:
        if aux_params[-1] == 0:
            jac[10] = frq ** 2
        if aux_params[-1] in [1, 2]:
            jac[10] = (frq - params[0]) ** 2
    if jac_flag[11] == 1.:
        jac[11] = (frq - params[0]) ** 3
    return jac.T


def calc_g_simple(p_self: float, p_for: float, temp: float,
                  g_self: float, g_for: float, ng_self: float, ng_for: float) -> float:
    return g_self * p_self * (t_ref / temp) ** ng_self + g_for * p_for * (t_ref / temp) ** ng_for


def mdl_multi(frq: np.ndarray, params: np.ndarray, aux_params: np.ndarray) -> np.ndarray:
    absor = np.empty_like(frq)
    nfil = int(max(aux_params[:, -1]))+1
    # print(params[n_const_par:])
    for ifil in range(nfil):
        ### preparing data for model profile generation
        ### within current recording filtered by the last column of aux_params
        tmp_where = np.where(aux_params[:, -1] == ifil)
        frq_cur = frq[tmp_where]  # frequencies for current recording number ifil
        p_s = aux_params[tmp_where][0, 0]  # self-pressure
        p_f = aux_params[tmp_where][0, 1]  # foreign-pressure
        tmpr = aux_params[tmp_where][0, 2]  # temperature
        dev = aux_params[tmp_where][0, 3]  # temperature
        rtype = int(aux_params[tmp_where][0, 4])  # line index
        clen = aux_params[tmp_where][0, 5]  # temperature        
        gdop = aux_params[tmp_where][0, -3]

        ### Calculation of parameters for the profile

        ### Line strength
        tmpS = aux_params[tmp_where][0, -2]  # line strength
        if rtype == 0:
            tmpS *= 1.E5
        else:
            tmpS *= clen
        ### Line center
        frq0 = params[1]
        ### G0 width (speed-independent HWHM)
        id_s = multi_params_indx["mg0s"]  #    instead of remebering indices of params
        id_f = multi_params_indx["mg0f"]  #    i put them into dictionary
        id_ns = multi_params_indx["mng0s"]  #  with more ore less intuitive keys
        id_nf = multi_params_indx["mng0f"]  #  more code but also more clarity
        
        # print('N = ', ifil, '!!!', p_s, p_f, tmpr, params[id_s], params[id_f], params[id_ns], params[id_nf])
        # print(type(tmpr), type(t_ref), type(params[id_nf]))
        
        tmpG0 = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        if rtype in [1, 2, 3]:          
            id_frab = multi_params_indx["mfrab"]
            frab = params[id_frab]
            tmpG0 = np.sqrt(tmpG0**2 + frab**2)  
            # for RAD and VID spectrometers saturation leads to extra HWHM
        ### G2 width (speed-dependend part)
        id_s = multi_params_indx["mg2s"]
        id_f = multi_params_indx["mg2f"]
        id_ns = multi_params_indx["mng2s"]
        id_nf = multi_params_indx["mng2f"]
        tmpG2 = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        ### SI pressure shift of central frq
        id_s = multi_params_indx["md0s"]
        id_f = multi_params_indx["md0f"]
        id_ns = multi_params_indx["mnd0s"]
        id_nf = multi_params_indx["mnd0f"]
        tmpD0 = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        ### SD part of pressure shift of central frq
        id_s = multi_params_indx["md2s"]
        id_f = multi_params_indx["md2f"]
        id_ns = multi_params_indx["mnd2s"]
        id_nf = multi_params_indx["mnd2f"]
        tmpD2 = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        ### mixing parameter
        id_s = multi_params_indx["my0s"]
        id_f = multi_params_indx["my0f"]
        id_ns = multi_params_indx["mny0s"]
        id_nf = multi_params_indx["mny0f"]
        tmpY0 = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        ### Vel.chng. collisions rate
        id_s = multi_params_indx["mnuvcs"]
        id_f = multi_params_indx["mnuvcf"]
        id_ns = multi_params_indx["mnnuvcs"]
        id_nf = multi_params_indx["mnnuvcf"]
        tmpnu = calc_g_simple(p_s, p_f, tmpr,
                              params[id_s], params[id_f], params[id_ns], params[id_nf])
        ### Continuum
        id_s = multi_params_indx["mcs"]
        id_f = multi_params_indx["mcf"]
        id_ns = multi_params_indx["mncs"]
        id_nf = multi_params_indx["mncf"]
        cs = params[id_s]
        cf = params[id_f]
        ncs = params[id_ns]
        ncf = params[id_nf]
        tmp_conti = p_s * (p_s * cs * (t_ref / tmpr)**ncs + p_f * cf * (t_ref / tmpr)**ncf)

        i_p_rec = n_const_par + ifil * n_add_par  # index of params related to current recording
        if rtype == 0:
            abs_htp, _ =  htp(frq0, gdop, tmpG0, tmpG2, tmpD0, tmpD2, tmpnu, 0., frq_cur, Ylm=tmpY0)       
            abs_htp[:] = abs_htp[:] * tmpS * params[i_p_rec] * (frq_cur[:] / frq0)**2
            abs_htp[:] += tmp_conti * frq_cur[:]**2
            absor[tmp_where] = abs_htp
        else:
            abs_htp = dif_htp(frq_cur, frq0, gdop, tmpG0, tmpG2, tmpD0, tmpD2, tmpnu, 
                              tmpS, dev, params[i_p_rec], params[i_p_rec+1], rtype)            
            abs_htp[:] += params[i_p_rec+2] + params[i_p_rec+3] * (frq_cur[:] - params[1]) \
            + params[i_p_rec+4] * (frq_cur - params[1]) ** 2 + params[i_p_rec+5] * (frq_cur - params[1]) ** 3
            absor[tmp_where] = abs_htp
                
    return absor


def mdljac_multi(frq: np.ndarray, jac_flag: np.ndarray, params: np.ndarray, aux_params: np.ndarray) -> np.ndarray:
    # print('Control output')
    # print(params[n_const_par::n_add_par])
    jac = np.zeros((len(frq), len(params)))    
    model1 = mdl_multi(frq, params, aux_params=aux_params)    
    dpar = 1.E-6
    # list of params requiring numeric derivative
    # generally these are: line center, all width, shift and mixing parameters
    # together with their T-dependence power factors
    params_deriv_numeric = [multi_params_indx[x] for x in multi_params_numeric_deriv]    
    # first take the list of params where numeric derivative is necessary (i put it into special tuple in constants.py)
    for ipar in params_deriv_numeric:
        if jac_flag[ipar] == 1.:  # if parameter number ipar is adjustable (jac_flag=1), calc derivative
            params2 = np.copy(params)  # copy of params
            params2[ipar] = params2[ipar] + dpar  # add small value to parameter number ipar
            model2 = mdl_multi(frq, params2, aux_params=aux_params)  # calculate model with changed parameter            
            jac[:, ipar] = (model2 - model1) / dpar  # the derivative by specified parameter            

    ### Continuum parameters derivatives
    # index of the corresponding parameters are:
    id_s = multi_params_indx["mcs"]
    id_f = multi_params_indx["mcf"]
    id_ns = multi_params_indx["mncs"]
    id_nf = multi_params_indx["mncf"]

    tmp_where = np.where(aux_params[:, -1] == 0)
    cs = params[id_s]
    cf = params[id_f]
    ncs = params[id_ns]
    ncf = params[id_nf]

    if jac_flag[id_s] == 1:
        jac[tmp_where][:, id_s] = aux_params[tmp_where][:, 0]**2 * frq[tmp_where][:]**2 * (t_ref / aux_params[tmp_where][:, 2])**ncs
    if jac_flag[id_f] == 1:
        jac[tmp_where][:, id_f] = aux_params[tmp_where][:, 0] * aux_params[tmp_where][:, 1] * frq[tmp_where][:] * (t_ref / aux_params[tmp_where][:, 2])**ncf
    if jac_flag[id_ns] == 1:
        pass  # currently let's leave T-dependence of continuum as read-only parameter with no adjustment
        #jac[tmp_where][:, id_ns] = aux_params[tmp_where][:, 0]**2 * cs * (t_ref / aux_params[tmp_where][:, 2])**ncs * np.log(ncs)
    if jac_flag[id_nf] == 1:
        pass
        #jac[tmp_where][:, id_nf] = aux_params[tmp_where][:, 0] * aux_params[tmp_where][:, 1] * cs * (t_ref / aux_params[tmp_where][:, 2])**ncs * np.log(ncs)
    
    ### Derivatives by the individual recording parameters
    nfil = int(max(aux_params[:, -1]))+1  # number of recordings

    ### scaling factor
    if jac_flag[n_const_par] == 1:  # only if scaling factor is adjustable (scaling is 0-th additional parameter)
        params2 = np.copy(params)  # first we make copy of params but with scale equal to 1 (deriv by scale is just)        
        for ifil in range(nfil):
            i_start = n_const_par + n_add_par * ifil  # index of first rec-related parameter
            params2[i_start] = 1.
            params2[i_start+2] = 0.
            params2[i_start+3] = 0.
            params2[i_start+4] = 0.
            params2[i_start+5] = 0.
        model2 = mdl_multi(frq, params2, aux_params=aux_params)
        for ifil in range(nfil):        
            i_start = n_const_par + n_add_par * ifil  # index of first rec-related parameter
            tmp_where = np.where(aux_params[:, -1] == ifil)  # filter for the points related to some recording                        
            jac[tmp_where, i_start] = model2[tmp_where]
            # print(jac[tmp_where][:, i_start])

    ### power factor (where applicable)
    if jac_flag[n_const_par+1] == 1:  # only if power factor is adjustable (1-th additional parameter)
        for ifil in range(nfil):
            rtype = int(aux_params[tmp_where][0, 4])  # record type
            if rtype in [1, 2]:  # only RAD and VID recs have power factor in model function
                i_start = n_const_par + n_add_par * ifil  # index of first rec-related parameter
                params2 = np.copy(params)
                params2[i_start+1] += dpar
                model2 = mdl_multi(frq, params2, aux_params=aux_params)
                tmp_where = np.where(aux_params[:, -1] == ifil)  # filter for the points related to some recording
                jac[tmp_where, i_start+1] = (model2[tmp_where] - model1[tmp_where]) / dpar
                #print(jac[tmp_where][:, i_start+1])


    for ifil in range(nfil):
        # tmp_where = np.where(aux_params[:, -1] == ifil)  # filter for the points related to some recording        
        # p_s = aux_params[tmp_where][0, 0]  # self-pressure
        # p_f = aux_params[tmp_where][0, 1]  # foreign-pressure
        # tmpr = aux_params[tmp_where][0, 2]  # temperature
        # dev = aux_params[tmp_where][0, 3]  # temperature
        # rtype = int(aux_params[tmp_where][0, 4])  # record type
        # clen = aux_params[tmp_where][0, 5]  # temperature        
        # gdop = aux_params[tmp_where][0, -3] # dopler width
        # tmpS = aux_params[tmp_where][0, -2]  # line strength
        # if rtype == 0:
        #     tmpS *= 1.E5
        # else:
        #     tmpS *= clen
        i_start = n_const_par + n_add_par * ifil  # index of first rec-related parameter
        if jac_flag[n_const_par+2] == 1:
            jac[tmp_where, i_start+2] = 1.            
        if jac_flag[n_const_par+3] == 1:            
            jac[tmp_where, i_start+3] = frq[tmp_where] - params[1]
        if jac_flag[n_const_par+4] == 1:
            jac[tmp_where, i_start+4] = (frq[tmp_where] - params[1])**2
        if jac_flag[n_const_par+5] == 1:
            jac[tmp_where, i_start+5] = (frq[tmp_where] - params[1])**3
    return jac
