import numpy as np

def qpart(Tq: float, Tref: float, tpart_c: np.ndarray, qpart_c: np.ndarray) -> float:
    # function returns value of the partition sum corresponding to the specified T
    # in relation to one at reference T (usually 296K)
    # piecewise linear interpolation is used
    # Tq - temperature where we want Q to be calculated
    # Tref - reference T
    # tpart_c, qpart_c - arrays of T vs Q values precalculated (usually from HITRAN)
    sort_ref = np.argsort(abs(tpart_c - Tref))
    if tpart_c[sort_ref[0]] == Tref:
        Qref = qpart_c[sort_ref[0]]
    else:
        TrefL = tpart_c[sort_ref[0]]  # Temperatures around Tref
        TrefR = tpart_c[sort_ref[1]]
        iL = sort_ref[0]
        iR = sort_ref[1]
        if TrefR < TrefL:
            (TrefL, TrefR) = (TrefR, TrefL)
            (iL, iR) = (iR, iL)
        Qref = qpart_c[iL] + (qpart_c[iR] - qpart_c[iL]) * ((Tref - TrefL) / (TrefR - TrefL))

    sort_q = np.argsort(abs(tpart_c - Tq))  # first two elements of this array point to the closest T in tpart_c
    # print(tpart_c[sort_ref[0]], tpart_c[sort_q[0]])
    if tpart_c[sort_q[0]] == Tq:
        return Qref / qpart_c[sort_q[0]]   # if Tq value is
    else:
        TrefL = tpart_c[sort_q[0]]  # Temperatures around Tq
        TrefR = tpart_c[sort_q[1]]
        iL = sort_q[0]
        iR = sort_q[1]
        if TrefR < TrefL:
            (TrefL, TrefR) = (TrefR, TrefL)
            (iL, iR) = (iR, iL)
        Qq = qpart_c[iL] + (qpart_c[iR] - qpart_c[iL]) * ((Tq - TrefL) / (TrefR - TrefL))
        return Qref / Qq
