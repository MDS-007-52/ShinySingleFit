import numpy as np
import scipy

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
               aux_params: np.ndarray = None) -> np.ndarray:
    # x_in - vector of x values (np.array)
    # y_in - vector of y values (np.array)
    # params0 - vector of params (np.array)
    # modelf - name of the function
    # jacobi - name of the jacobian function (partial derivatives of modelf at each x_in value)
    # == modelf should have x_in (array) and params (array) as an arguments, returns len(x_in) array
    # == jacobi should have x_in (array), jac_flag (array) and params (array, len(params)==len(jac_flag))
    # as an arguments, returns len(params)*len(x_in) list/array
    # jac_flag - list of 0 and 1s showing which parameter is adjustable and which is left constant
    # fit_precision - relative precision for RMS to stop the iterations (if RMS relative change
    # is less than fit_precision, the fit procedure stops )
    # lmstart - starting value of lambda in Levenberg-Marquardt method
    # lmstep - multiplier for lambda in case of too big step of the parameters
    # iter_limit - max number of iterations to decide that fit does not converge
    # aux_params - some additional conditions modelf and jacobi may have
    # (e.g., pressure and temperature to calculate line intensity, or list of P, T and other additional things
    # which are used in multifit procedure)

    cur_params = np.copy(params0)  # current params
    if aux_params is None:
        model1 = modelf(x_in, cur_params) # calculated absorption at initial parameters
    else:
        model1 = modelf(x_in, cur_params, aux_params) # calculated absorption at initial parameters

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
        # print('FIT step = ', k, 'substep = ', i, ', jacobi done')
        f_step = False  # flag of parameters step is done
        while not f_step:
            # Levenberg-Markquardt method to calc the step of the params to less residual
            am1 = ap * np.eye(len(params0)) + np.matmul(jac.T, jac)
            resid1 = y_in - model1
            av = np.matmul(jac.T, resid1)
            #am2 = np.linalg.inv(am1)
            am2 = scipy.linalg.inv(am1)  #np.linalg.inv(am1)
            steparams = np.matmul(am2, av)  # step of parameters due to Leven-Mark method
            tmp_params = steparams + cur_params
            # print('Params step:')
            # print(steparams)
            # print('Next params:')
            # print(tmp_params)
            # print('FIT step = ', k, 'substep = ', i, ', calcuating new residual')
            if aux_params is None:
                model1 = modelf(x_in, tmp_params)  # model calc for new parameters
            else:
                model1 = modelf(x_in, tmp_params, aux_params)  # model calc for new parameters
            rms2 = np.sqrt(np.sum((y_in - model1) ** 2))
            # print('old rms = ', rms1, '   new rms = ', rms2, '   lambda = ', ap)
            # print('relative difference = ', abs((rms2 - rms1) / rms1))
            if rms2 < rms1:  # case when new residual is smaller
                if abs((rms2 - rms1) / rms1) > fit_precision:  # change of rms is great enough
                    # print('next rms is smaller, STEP FORWARD')
                    f_step = True
                    cur_params += steparams
                    rms1 = rms2
                elif abs((rms2 - rms1) / rms1) <= fit_precision:  # change of rms is smaller then threshold
                    # print('next rms is smaller, but change is less than threshold, fit STOP')
                    cur_params += steparams
                    rms1 = rms2
                    f_step = True
                    f_end = True
            elif rms2 >= rms1:
                if abs((rms2 - rms1) / rms1) < fit_precision:
                    # print('next rms is greater, but change is less than threshold, fit STOP')
                    f_step = True
                    f_end = True
                elif abs((rms2 - rms1) / rms1) > fit_precision:
                    if i <= iter_limit:
                        ap = ap * lmstep
                        # print('next rms is greater, change LAMBDA and repeat')
                    if i > iter_limit:
                        f_end = True
                        f_step = True
                        # print('can"t converge to a solution, BREAK')
            i += 1
        k += 1
    return cur_params


def fit_uncertainties(x_in: np.ndarray,
                      y_in: np.ndarray,
                      params0: np.ndarray,
                      modelf,
                      jacobi,
                      jac_flag: np.ndarray,
                      aux_params: np.ndarray = None) -> np.ndarray:
    # x_in - vector of x values (np.array)
    # y_in - vector of y values (np.array)
    # params0 - vector of params given by fit_params routine (np.array)
    # modelf - name of the function
    # jacobi - name of the jacobian function (partial derivatives of modelf at each x_in value)
    # == modelf should have x_in (array) and params (array) as an arguments, returns len(x_in) array
    # == jacobi should have x_in (array), jac_flag (array) and params (array, len(params)==len(jac_flag))
    # as an arguments, returns len(params)*len(x_in) list/array
    # jac_flag - list of 0 and 1s showing which parameter is adjustable and which is left constant
    # aux_params - some additional conditions modelf and jacobi may have
    # (e.g., pressure and temperature to calculate line intensity, or list of P, T and other additional things
    # which are used in multifit procedure)

    if aux_params is None:
        model1 = modelf(x_in, params0) # calculated absorption at initial parameters
    else:
        model1 = modelf(x_in, params0, aux_params) # calculated absorption at initial parameters

    if aux_params is None:
        jac = jacobi(x_in, jac_flag, params0)
    else:
        jac = jacobi(x_in, jac_flag, params0, aux_params)

    rms = np.sqrt(np.sum((y_in - model1) ** 2) / (len(x_in) - 1))

    am1 = np.matmul(jac.T, jac)

    for iam in range(len(params0)):
        if am1[iam, iam] == 0:
            am1[iam, iam] = 1.E-6

    am1 = np.linalg.inv(am1)

    params_err = np.empty_like(params0)

    for ipar in range(len(params_err)):
        params_err[ipar] = rms * np.sqrt(am1[ipar, ipar]) * jac_flag[ipar]

    return params_err
