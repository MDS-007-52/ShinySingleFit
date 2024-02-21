import numpy as np
import scipy
import math
from constants import *


def fit_params(x_in: np.ndarray,
               y_in: np.ndarray,
               params0: np.ndarray,
               modelf,
               jacobi,
               jac_flag: np.ndarray,
               fit_precision: float,
               lmstart: float,
               lmstep: float,
               iter_limit: int,
               aux_params: np.ndarray = None,
               f_verbose_fit: bool = None) -> np.ndarray:
    """
    This function fits params of the specially constructed model function modelf to make it fit the best (with minimum
    standard deviation of the difference) to the experimental data (x_in, y_in)
    :param x_in: vector of x values
    :param y_in: vector of y values
    :param params0: vector of initial params
    :param modelf: name of the model function
    :param jacobi: name of the jacobian function (partial derivatives of modelf by each parameter at each x_in value)
    :param jac_flag: list of 0 and 1s showing which parameter is adjustable and which is left constant
    :param fit_precision: relative precision for RMS to stop the iterations (if RMS relative change is less than fit_precision, the fit procedure stops )
    :param lmstart: starting value of lambda in Levenberg-Marquardt method
    :param lmstep: multiplier for lambda in case of too big step of the parameters
    :param iter_limit: max number of iterations to decide that fit does not converge
    :param aux_params: some additional conditions modelf and jacobi may have (e.g., pressure and temperature to calculate line intensity, or list of P, T and other additional things which are used in fit procedure)
    :return: vector of parameters, the same size as params0
    """

    cur_params = np.copy(params0)  # current params
    if aux_params is None:
        model1 = modelf(x_in, cur_params)  # calculated absorption at initial parameters
    else:
        model1 = modelf(x_in, cur_params, aux_params)  # calculated absorption at initial parameters

    rms1 = np.sqrt(np.sum((y_in - model1) ** 2))  # RMS of the measured-minus-calc difference at initial parameters
    # print('Initial rms within fit procedure: ', rms1)
    # print('Initial parameters')
    # print(cur_params)

    k = 0  # outer iteration number
    f_end = False  # flag of fit ending

    while not f_end:
        ap = lmstart  # Leven-Mark lambda initiate
        i = 0  # inner level iteration number
        if aux_params is None:
            jac = jacobi(x_in, jac_flag, cur_params)
        else:
            jac = jacobi(x_in, jac_flag, cur_params, aux_params)
        if f_verbose_fit:
            print('FIT step = ', k, 'substep = ', i, ', jacobi done')
        f_step = False  # flag of parameters step is done
        fit_result_status = ' '
        while not f_step:
            # Levenberg-Markquardt method to calc the step of the params to less residual
            am1 = ap * np.eye(len(params0)) + np.matmul(jac.T, jac)
            resid1 = y_in - model1
            av = np.matmul(jac.T, resid1)
            try:
                am2 = scipy.linalg.inv(am1)  # matrix inversion
            except:
                am2 = scipy.linalg.pinv(am1)  # if matrix is badly conditioned, alternative method is used
            steparams = np.matmul(am2, av)  # step of parameters due to Leven-Mark method
            tmp_params = steparams + cur_params
            if f_verbose_fit:
                # print('Params step:')
                # print(steparams)
                print('Next params:')
                print(tmp_params)
                print('FIT step = ', k, 'substep = ', i, ', calculating new residual')
            if abs(tmp_params[2]) > 1.E10:
                print('Something wrong, fit gone too far')
                fit_result_status = 'Unconverged fit, unreasonably high parameter values. The most latest stable parameters set is returned.'
                f_step = True
                f_end = True
                break
            if aux_params is None:
                model1 = modelf(x_in, tmp_params)  # model calc for new parameters
            else:
                model1 = modelf(x_in, tmp_params, aux_params)  # model calc for new parameters
            rms2 = np.sqrt(np.sum((y_in - model1) ** 2))
            if f_verbose_fit:
                print('old rms = ', rms1, '   new rms = ', rms2, '   lambda = ', ap)
                print('relative difference = ', abs((rms2 - rms1) / rms1))            
            if math.isnan(rms2):
                print('Something wrong, NaN values acquired.')
                fit_result_status = 'Unconverged fit. The most latest stable parameters set is returned.'
                f_step = True
                f_end = True
                break
            if rms2 < rms1:  # case when new residual is smaller
                if abs((rms2 - rms1) / rms1) > fit_precision:  # change of rms is great enough
                    if f_verbose_fit:
                        print('next rms is smaller, STEP FORWARD')
                    f_step = True
                    cur_params += steparams
                    rms1 = rms2
                elif abs((rms2 - rms1) / rms1) <= fit_precision:  # change of rms is smaller then threshold
                    if f_verbose_fit:
                        print('next rms is smaller, but change is less than threshold, fit STOP')
                    fit_result_status = 'Converged successfully in '+str(k)+' steps.'
                    cur_params += steparams
                    rms1 = rms2
                    f_step = True
                    f_end = True
            elif rms2 >= rms1:
                if abs((rms2 - rms1) / rms1) < fit_precision:
                    if f_verbose_fit:
                        print('next rms is greater, but change is less than threshold, fit STOP')
                    fit_result_status = 'Converged successfully in '+str(k)+' steps. Slow convergence.'
                    f_step = True
                    f_end = True
                elif abs((rms2 - rms1) / rms1) > fit_precision:
                    if i <= iter_limit:
                        ap = ap * lmstep
                        if f_verbose_fit:
                            print('next rms is greater, change LAMBDA and repeat')
                    if i > iter_limit:
                        if f_verbose_fit:
                            print('can"t converge to a solution, BREAK')
                        f_end = True
                        f_step = True                        
                        fit_result_status = 'Fit can not converge to a solution. Returning the result of the latest iteration.'
            elif rms2 is None:
                if f_verbose_fit:
                    print('Something went wrong')
                f_end = True
                f_step = True

            i += 1
        k += 1
    return cur_params, fit_result_status


def fit_uncertainties(x_in: np.ndarray,
                      y_in: np.ndarray,
                      params0: np.ndarray,
                      modelf,
                      jacobi,
                      jac_flag: np.ndarray,
                      aux_params: np.ndarray = None) -> np.ndarray:
    """
    This function calculates uncertainties of the parameters defined by fit_params function
    All the params of this function should be the same that used for fit_params, except for the params0
    :param x_in: vector of x values
    :param y_in: vector of y values
    :param params0: vector of params given by fit_params function
    :param modelf: name of the model function
    :param jacobi: name of the jacobian function
    :param jac_flag:  list of 0 and 1s showing which parameter is adjustable and which is left constant
    :param aux_params: some additional conditions modelf and jacobi may have
    :return: vector of the same size as params0 containing uncertainties of the corresponding parameters
    """

    if aux_params is None:
        model1 = modelf(x_in, params0)  # calculated absorption at initial parameters
    else:
        model1 = modelf(x_in, params0, aux_params)  # calculated absorption at initial parameters

    jac_flag_local = jac_flag[0:n_const_par+n_add_par]

    if aux_params is None:
        jac = jacobi(x_in, jac_flag_local, params0)
    else:
        jac = jacobi(x_in, jac_flag_local, params0, aux_params)

    rms = np.sqrt(np.sum((y_in - model1) ** 2) / (len(x_in) - 1))

    am1 = np.matmul(jac.T, jac)

    for iam in range(len(params0)):
        if am1[iam, iam] == 0:
            am1[iam, iam] = 1.E-20

    try: 
        am2 = scipy.linalg.inv(am1)
    except:
        am2 = scipy.linalg.pinv(am1)

    params_err = np.empty_like(params0)

    # print(params0.shape, am1.shape, len(jac_flag_local))

    for ipar in range(len(params_err)):
        try:
            # print(ipar, am2[ipar, ipar])
            params_err[ipar] = rms * np.sqrt(abs(am2[ipar, ipar])) * jac_flag[ipar]
        except:            
            # params_err[ipar] = rms * np.sqrt(abs(am2[ipar, ipar])) * jac_flag[ipar]
            params_err[ipar, ipar] = -1.

    return params_err
