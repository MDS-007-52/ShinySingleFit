from numpy import zeros, array, setdiff1d, ndarray, arange
from numpy import place, where, real, polyval
from numpy import complex128, int64, float64
from numpy import sqrt, abs, exp, pi, log, tan
from numpy import flipud
from numpy.fft import fft, fftshift
from numpy import any, minimum, maximum

__ComplexType__ = complex128
__IntegerType__ = int64
__FloatType__ = float64

# ------------------ LINESHAPES -----------------------------------------

# ------------------ complex probability function -----------------------
# define static data
zone = __ComplexType__(1.0e0 + 0.0e0j)
zi = __ComplexType__(0.0e0 + 1.0e0j)
tt = __FloatType__(
    [0.5e0, 1.5e0, 2.5e0, 3.5e0, 4.5e0, 5.5e0, 6.5e0, 7.5e0, 8.5e0, 9.5e0, 10.5e0, 11.5e0, 12.5e0, 13.5e0, 14.5e0])
pipwoeronehalf = __FloatType__(0.564189583547756e0)


# "naive" implementation for benchmarks
def cpf3(X, Y):
    # X,Y,WR,WI - numpy arrays
    if type(X) != ndarray:
        if type(X) not in set([list, tuple]):
            X = array([X])
        else:
            X = array(X)
    if type(Y) != ndarray:
        if type(Y) not in set([list, tuple]):
            Y = array([Y])
        else:
            Y = array(Y)

    zm1 = zone / __ComplexType__(X + zi * Y)  # maybe redundant
    zm2 = zm1 ** 2
    zsum = zone
    zterm = zone

    for tt_i in tt:
        zterm *= zm2 * tt_i
        zsum += zterm

    zsum *= zi * zm1 * pipwoeronehalf

    return zsum.real, zsum.imag


T = __FloatType__([0.314240376e0, 0.947788391e0, 1.59768264e0, 2.27950708e0, 3.02063703e0, 3.8897249e0])
U = __FloatType__([1.01172805e0, -0.75197147e0, 1.2557727e-2, 1.00220082e-2, -2.42068135e-4, 5.00848061e-7])
S = __FloatType__([1.393237e0, 0.231152406e0, -0.155351466e0, 6.21836624e-3, 9.19082986e-5, -6.27525958e-7])


# Complex probability function implementation (Humlicek)
def cpf(X, Y):
    # X,Y,WR,WI - numpy arrays
    if type(X) != ndarray:
        if type(X) not in set([list, tuple]):
            X = array([X])
        else:
            X = array(X)
    if type(Y) != ndarray:
        if type(Y) not in set([list, tuple]):
            Y = array([Y])
        else:
            Y = array(Y)

    # REGION3
    index_REGION3 = where(sqrt(X ** 2 + Y ** 2) > __FloatType__(8.0e0))
    X_REGION3 = X[index_REGION3]
    Y_REGION3 = Y[index_REGION3]
    zm1 = zone / __ComplexType__(X_REGION3 + zi * Y_REGION3)
    zm2 = zm1 ** 2
    zsum_REGION3 = zone
    zterm = zone
    for tt_i in tt:
        zterm *= zm2 * tt_i
        zsum_REGION3 += zterm
    zsum_REGION3 *= zi * zm1 * pipwoeronehalf

    index_REGION12 = setdiff1d(array(arange(len(X))), array(index_REGION3))
    X_REGION12 = X[index_REGION12]
    Y_REGION12 = Y[index_REGION12]

    WR = __FloatType__(0.0e0)
    WI = __FloatType__(0.0e0)

    # REGION12
    Y1_REGION12 = Y_REGION12 + __FloatType__(1.5e0)
    Y2_REGION12 = Y1_REGION12 ** 2

    # REGION2
    subindex_REGION2 = where((Y_REGION12 <= 0.85e0) &
                             (abs(X_REGION12) >= (18.1e0 * Y_REGION12 + 1.65e0)))

    index_REGION2 = index_REGION12[subindex_REGION2]

    X_REGION2 = X[index_REGION2]
    Y_REGION2 = Y[index_REGION2]
    Y1_REGION2 = Y1_REGION12[subindex_REGION2]
    Y2_REGION2 = Y2_REGION12[subindex_REGION2]
    Y3_REGION2 = Y_REGION2 + __FloatType__(3.0e0)

    WR_REGION2 = WR
    WI_REGION2 = WI

    WR_REGION2 = zeros(len(X_REGION2))
    ii = abs(X_REGION2) < __FloatType__(12.0e0)
    WR_REGION2[ii] = exp(-X_REGION2[ii] ** 2)
    WR_REGION2[~ii] = WR

    for I in range(6):
        R_REGION2 = X_REGION2 - T[I]
        R2_REGION2 = R_REGION2 ** 2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D1_REGION2 = Y1_REGION2 * D_REGION2
        D2_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (U[I] * (R_REGION2 * D2_REGION2 - 1.5e0 * D1_REGION2) +
                                               S[I] * Y3_REGION2 * D2_REGION2) / (R2_REGION2 + 2.25e0)
        R_REGION2 = X_REGION2 + T[I]
        R2_REGION2 = R_REGION2 ** 2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D3_REGION2 = Y1_REGION2 * D_REGION2
        D4_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (U[I] * (R_REGION2 * D4_REGION2 - 1.5e0 * D3_REGION2) -
                                               S[I] * Y3_REGION2 * D4_REGION2) / (R2_REGION2 + 2.25e0)
        WI_REGION2 = WI_REGION2 + U[I] * (D2_REGION2 + D4_REGION2) + S[I] * (D1_REGION2 - D3_REGION2)

    # REGION3
    index_REGION1 = setdiff1d(array(index_REGION12), array(index_REGION2))
    X_REGION1 = X[index_REGION1]
    Y_REGION1 = X[index_REGION1]

    subindex_REGION1 = setdiff1d(array(arange(len(index_REGION12))), array(subindex_REGION2))
    Y1_REGION1 = Y1_REGION12[subindex_REGION1]
    Y2_REGION1 = Y2_REGION12[subindex_REGION1]

    WR_REGION1 = WR
    WI_REGION1 = WI

    for I in range(6):
        R_REGION1 = X_REGION1 - T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1 ** 2 + Y2_REGION1)
        D1_REGION1 = Y1_REGION1 * D_REGION1
        D2_REGION1 = R_REGION1 * D_REGION1
        R_REGION1 = X_REGION1 + T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1 ** 2 + Y2_REGION1)
        D3_REGION1 = Y1_REGION1 * D_REGION1
        D4_REGION1 = R_REGION1 * D_REGION1

        WR_REGION1 = WR_REGION1 + U[I] * (D1_REGION1 + D3_REGION1) - S[I] * (D2_REGION1 - D4_REGION1)
        WI_REGION1 = WI_REGION1 + U[I] * (D2_REGION1 + D4_REGION1) + S[I] * (D1_REGION1 - D3_REGION1)

    # total result
    WR_TOTAL = zeros(len(X))
    WI_TOTAL = zeros(len(X))
    # REGION3
    WR_TOTAL[index_REGION3] = zsum_REGION3.real
    WI_TOTAL[index_REGION3] = zsum_REGION3.imag
    # REGION2
    WR_TOTAL[index_REGION2] = WR_REGION2
    WI_TOTAL[index_REGION2] = WI_REGION2
    # REGION1
    WR_TOTAL[index_REGION1] = WR_REGION1
    WI_TOTAL[index_REGION1] = WI_REGION1

    return WR_TOTAL, WI_TOTAL


#hcpf = cpf  # stub for initial cpf

# ------------------ Schreier CPF ------------------------

# "Optimized implementations of rational approximations
#  for the Voigt and complex error function".
# Franz Schreier. JQSRT 112 (2011) 1010-10250
# doi:10.1016/j.jqsrt.2010.12.010

# Enable this if numpy.polyval doesn't perform well.
"""    
def polyval(p, x):
    y = zeros(x.shape, dtype=float)
    for i, v in enumerate(p):
        y *= x
        y += v
    return y
""";


def cef(x, y, N):
    # Computes the function w(z) = exp(-zA2) erfc(-iz) using a rational
    # series with N terms. It is assumed that Im(z) > 0 or Im(z) = 0.
    z = x + 1.0j * y
    M = 2 * N;
    M2 = 2 * M;
    k = arange(-M + 1, M)  # '; # M2 = no. of sampling points.
    L = sqrt(N / sqrt(2));  # Optimal choice of L.
    theta = k * pi / M;
    t = L * tan(theta / 2);  # Variables theta and t.
    # f = exp(-t.A2)*(LA2+t.A2); f = [0; f]; # Function to be transformed.
    f = zeros(len(t) + 1);
    f[0] = 0
    f[1:] = exp(-t ** 2) * (L ** 2 + t ** 2)
    # f = insert(exp(-t**2)*(L**2+t**2),0,0)
    a = real(fft(fftshift(f))) / M2;  # Coefficients of transform.
    a = flipud(a[1:N + 1]);  # Reorder coefficients.
    Z = (L + 1.0j * z) / (L - 1.0j * z);
    p = polyval(a, Z);  # Polynomial evaluation.
    w = 2 * p / (L - 1.0j * z) ** 2 + (1 / sqrt(pi)) / (L - 1.0j * z);  # Evaluate w(z).
    return w


# weideman24 by default
# weideman24 = lambda x,y: cef(x,y,24)
# weideman = lambda x, y, n: cef(x, y, n)


def hum1_wei(x, y, n=24):
    t = y - 1.0j * x
    cerf = 1 / sqrt(pi) * t / (0.5 + t ** 2)
    """
    z = x+1j*y
    cerf = 1j*z/sqrt(pi)/(z**2-0.5)
    """
    mask = abs(x) + y < 15.0
    if any(mask):
        w24 = cef(x[mask], y[mask], n)
        place(cerf, mask, w24)
    return cerf.real, cerf.imag


# VARIABLES['CPF'] = hum1_wei


# VARIABLES['CPF'] = cpf

# ------------------ Hartmann-Tran Profile (HTP) ------------------------
def htp(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, eta, sg, Ylm=0.0):
    """

    :param sg0: Unperturbed line position in cm-1
    :param GamD: Doppler HWHM
    :param Gam0: Speed-averaged line-width
    :param Gam2: Speed dependence of the line-width
    :param Shift0: Speed-averaged line-shift
    :param Shift2: Speed dependence of the line-shift
    :param anuVC: Velocity-changing frequency
    :param eta: Correlation parameter, No unit
    :param sg: Current WaveNumber of the Computation (vector)
    :param Ylm: 1st order (Rosenkranz) line mixing coefficient
    :return: LS_pCqSDHC_R, LS_pCqSDHC_I: Real/Imaginary part of the normalized spectral shape
    """
    # -------------------------------------------------
    #      "pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
    #      Subroutine to Compute the complex normalized spectral shape of an
    #      isolated line by the pCqSDHC model
    #
    #      Reference:
    #      H. Tran, N.H. Ngo, J.-M. Hartmann.
    #      Efficient computation of some speed-dependent isolated line profiles.
    #      JQSRT, Volume 129, November 2013, Pages 199–203
    #      http://dx.doi.org/10.1016/j.jqsrt.2013.06.015
    #
    #      Input/Output Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      T          : Temperature in Kelvin (Input).
    #      amM1       : Molar mass of the absorber in g/mol(Input).
    #      sg0        : Unperturbed line position in cm-1 (Input).
    #      GamD       : Doppler HWHM in cm-1 (Input)
    #      Gam0       : Speed-averaged line-width in cm-1 (Input).
    #      Gam2       : Speed dependence of the line-width in cm-1 (Input).
    #      anuVC      : Velocity-changing frequency in cm-1 (Input).
    #      eta        : Correlation parameter, No unit (Input).
    #      Shift0     : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2     : Speed dependence of the line-shift in cm-1 (Input)
    #      sg         : Current WaveNumber of the Computation in cm-1 (Input).
    #      Ylm        : 1st order (Rosenkranz) line mixing coefficients in cm-1 (Input)
    #
    #      Output Quantities (through Common Statements)
    #      -----------------
    #      LS_pCqSDHC_R: Real part of the normalized spectral shape (cm)
    #      LS_pCqSDHC_I: Imaginary part of the normalized spectral shape (cm)
    #
    #      Called Routines: 'CPF'      (Complex Probability Function)
    #      ---------------  'CPF3'      (Complex Probability Function for the region 3)
    #
    #      Called By: Main Program
    #      ---------
    #
    #     Double Precision Version
    #
    # -------------------------------------------------

    # sg is the only vector argument which is passed to function

    if type(sg) not in set([array, ndarray, list, tuple]):
        sg = array([sg])

    number_of_points = len(sg)
    Aterm_GLOBAL = zeros(number_of_points, dtype=__ComplexType__)
    Bterm_GLOBAL = zeros(number_of_points, dtype=__ComplexType__)

    cte = sqrt(log(2.0e0)) / GamD
    rpi = sqrt(pi)
    iz = __ComplexType__(0.0e0 + 1.0e0j)

    c0 = __ComplexType__(Gam0 + 1.0e0j * Shift0)
    c2 = __ComplexType__(Gam2 + 1.0e0j * Shift2)
    c0t = __ComplexType__((1.0e0 - eta) * (c0 - 1.5e0 * c2) + anuVC)
    c2t = __ComplexType__((1.0e0 - eta) * c2)

    # PART1
    if abs(c2t) == 0.0e0:
        Z1 = (iz * (sg0 - sg) + c0t) * cte
        xZ1 = -Z1.imag
        yZ1 = Z1.real
        WR1, WI1 = hum1_wei(xZ1, yZ1)
        Aterm_GLOBAL = rpi * cte * __ComplexType__(WR1 + 1.0e0j * WI1)
        index_Z1 = abs(Z1) <= 4.0e3
        index_NOT_Z1 = ~index_Z1
        if any(index_Z1):
            Bterm_GLOBAL = rpi * cte * ((1.0e0 - Z1 ** 2) * __ComplexType__(WR1 + 1.0e0j * WI1) + Z1 / rpi)
        if any(index_NOT_Z1):
            Bterm_GLOBAL = cte * (rpi * __ComplexType__(WR1 + 1.0e0j * WI1) + 0.5e0 / Z1 - 0.75e0 / (Z1 ** 3))
    else:
        # PART2, PART3 AND PART4   (PART4 IS A MAIN PART)

        # X - vector, Y - scalar
        X = (iz * (sg0 - sg) + c0t) / c2t
        Y = __ComplexType__(1.0e0 / ((2.0e0 * cte * c2t)) ** 2)
        csqrtY = (Gam2 - iz * Shift2) / (2.0e0 * cte * (1.0e0 - eta) * (Gam2 ** 2 + Shift2 ** 2))

        index_PART2 = abs(X) <= 3.0e-8 * abs(Y)
        index_PART3 = (abs(Y) <= 1.0e-15 * abs(X)) & ~index_PART2
        index_PART4 = ~ (index_PART2 | index_PART3)

        # PART4
        if any(index_PART4):
            X_TMP = X[index_PART4]
            Z1 = sqrt(X_TMP + Y) - csqrtY
            Z2 = Z1 + __FloatType__(2.0e0) * csqrtY
            xZ1 = -Z1.imag
            yZ1 = Z1.real
            xZ2 = -Z2.imag
            yZ2 = Z2.real
            SZ1 = sqrt(xZ1 ** 2 + yZ1 ** 2)
            SZ2 = sqrt(xZ2 ** 2 + yZ2 ** 2)
            DSZ = abs(SZ1 - SZ2)
            SZmx = maximum(SZ1, SZ2)
            SZmn = minimum(SZ1, SZ2)
            length_PART4 = len(index_PART4)
            WR1_PART4 = zeros(length_PART4)
            WI1_PART4 = zeros(length_PART4)
            WR2_PART4 = zeros(length_PART4)
            WI2_PART4 = zeros(length_PART4)
            index_CPF3 = (DSZ <= 1.0e0) & (SZmx > 8.0e0) & (SZmn <= 8.0e0)
            index_CPF = ~index_CPF3  # can be removed
            if any(index_CPF3):
                WR1, WI1 = cpf3(xZ1[index_CPF3], yZ1[index_CPF3])
                WR2, WI2 = cpf3(xZ2[index_CPF3], yZ2[index_CPF3])
                WR1_PART4[index_CPF3] = WR1
                WI1_PART4[index_CPF3] = WI1
                WR2_PART4[index_CPF3] = WR2
                WI2_PART4[index_CPF3] = WI2
            if any(index_CPF):
                WR1, WI1 = hum1_wei(xZ1[index_CPF], yZ1[index_CPF])
                WR2, WI2 = hum1_wei(xZ2[index_CPF], yZ2[index_CPF])
                WR1_PART4[index_CPF] = WR1
                WI1_PART4[index_CPF] = WI1
                WR2_PART4[index_CPF] = WR2
                WI2_PART4[index_CPF] = WI2

            Aterm = rpi * cte * (__ComplexType__(WR1_PART4 + 1.0e0j * WI1_PART4) - __ComplexType__(
                WR2_PART4 + 1.0e0j * WI2_PART4))
            Bterm = (-1.0e0 +
                     rpi / (2.0e0 * csqrtY) * (1.0e0 - Z1 ** 2) * __ComplexType__(WR1_PART4 + 1.0e0j * WI1_PART4) -
                     rpi / (2.0e0 * csqrtY) * (1.0e0 - Z2 ** 2) * __ComplexType__(WR2_PART4 + 1.0e0j * WI2_PART4)) / c2t
            Aterm_GLOBAL[index_PART4] = Aterm
            Bterm_GLOBAL[index_PART4] = Bterm

        # PART2
        if any(index_PART2):
            X_TMP = X[index_PART2]
            Z1 = (iz * (sg0 - sg[index_PART2]) + c0t) * cte
            Z2 = sqrt(X_TMP + Y) + csqrtY
            xZ1 = -Z1.imag
            yZ1 = Z1.real
            xZ2 = -Z2.imag
            yZ2 = Z2.real
            WR1_PART2, WI1_PART2 = hum1_wei(xZ1, yZ1)
            WR2_PART2, WI2_PART2 = hum1_wei(xZ2, yZ2)
            Aterm = rpi * cte * (__ComplexType__(WR1_PART2 + 1.0e0j * WI1_PART2) - __ComplexType__(
                WR2_PART2 + 1.0e0j * WI2_PART2))
            Bterm = (-1.0e0 +
                     rpi / (2.0e0 * csqrtY) * (1.0e0 - Z1 ** 2) * __ComplexType__(WR1_PART2 + 1.0e0j * WI1_PART2) -
                     rpi / (2.0e0 * csqrtY) * (1.0e0 - Z2 ** 2) * __ComplexType__(WR2_PART2 + 1.0e0j * WI2_PART2)) / c2t
            Aterm_GLOBAL[index_PART2] = Aterm
            Bterm_GLOBAL[index_PART2] = Bterm

        # PART3
        if any(index_PART3):
            X_TMP = X[index_PART3]
            xZ1 = -sqrt(X_TMP + Y).imag
            yZ1 = sqrt(X_TMP + Y).real
            WR1_PART3, WI1_PART3 = hum1_wei(xZ1, yZ1)
            index_ABS = abs(sqrt(X_TMP)) <= 4.0e3
            index_NOT_ABS = ~index_ABS
            Aterm = zeros(len(index_PART3), dtype=__ComplexType__)
            Bterm = zeros(len(index_PART3), dtype=__ComplexType__)
            if any(index_ABS):
                xXb = -sqrt(X).imag
                yXb = sqrt(X).real
                WRb, WIb = hum1_wei(xXb, yXb)
                Aterm[index_ABS] = (2.0e0 * rpi / c2t) * (
                            1.0e0 / rpi - sqrt(X_TMP[index_ABS]) * __ComplexType__(WRb + 1.0e0j * WIb))
                Bterm[index_ABS] = (1.0e0 / c2t) * (-1.0e0 +
                                                    2.0e0 * rpi * (1.0e0 - X_TMP[index_ABS] - 2.0e0 * Y) * (
                                                                1.0e0 / rpi - sqrt(X_TMP[index_ABS]) * __ComplexType__(
                                                            WRb + 1.0e0j * WIb)) +
                                                    2.0e0 * rpi * sqrt(X_TMP[index_ABS] + Y) * __ComplexType__(
                            WR1_PART3 + 1.0e0j * WI1_PART3))
            if any(index_NOT_ABS):
                Aterm[index_NOT_ABS] = (1.0e0 / c2t) * (
                            1.0e0 / X_TMP[index_NOT_ABS] - 1.5e0 / (X_TMP[index_NOT_ABS] ** 2))
                Bterm[index_NOT_ABS] = (1.0e0 / c2t) * (-1.0e0 + (1.0e0 - X_TMP[index_NOT_ABS] - 2.0e0 * Y) *
                                                        (1.0e0 / X_TMP[index_NOT_ABS] - 1.5e0 / (
                                                                    X_TMP[index_NOT_ABS] ** 2)) +
                                                        2.0e0 * rpi * sqrt(X_TMP[index_NOT_ABS] + Y) * __ComplexType__(
                            WR1 + 1.0e0j * WI1))
            Aterm_GLOBAL[index_PART3] = Aterm
            Bterm_GLOBAL[index_PART3] = Bterm

    # common part
    # LINE MIXING PART NEEDS FURTHER TESTING!!!
    LS_pCqSDHC = (1.0e0 / pi) * (
                Aterm_GLOBAL / (1.0e0 - (anuVC - eta * (c0 - 1.5e0 * c2)) * Aterm_GLOBAL + eta * c2 * Bterm_GLOBAL))
    return LS_pCqSDHC.real + Ylm * LS_pCqSDHC.imag, LS_pCqSDHC.imag
