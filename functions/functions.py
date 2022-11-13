"""
File containing the functions to calculate the main quantities of the project
"""

from math import pi, log10, floor, sqrt

from l2_massloss_interpolation.L2massloss_fork import constants
from l2_massloss_interpolation.L2massloss_fork.functions.kappa_functions import kap

def f1_the_T_fL2(the, T, fL2, c1, c2):
    return c1*T**4 * the**3 / (1-fL2) - the**2 + c2*T

def f2_the_T_fL2(the, T, fL2, c4, M1dot, alpha_ss, omgK, Rd, c3, PhiRd, GM2, PhiL1, PhiL2, lgRgrid, intp_lgkapgrid):
    x = c4 * (T*the)**3 / (1-fL2)
    U_over_P = (1.5 + x)/(1 + 1./3 * x)
    rho = (1-fL2) * M1dot / (2*pi * alpha_ss * omgK * Rd**3) / the**3
    print(fL2, M1dot, alpha_ss, omgK, Rd, the)
    print(rho)
    return 7./4 - (1.5 * U_over_P + c3*T**4 / kap(rho=rho, T=T, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid) / (1-fL2)**2) * the**2 \
        - PhiRd/(GM2/Rd) + (PhiL1 - fL2*PhiL2)/(GM2/Rd)/(1-fL2)

def T_the_nofL2(the, logthegrid, logTarr, dlogthe):   # under the assumption fL2=0
    # return ((the**2/c2)**(-s) + ((c1*the)**-0.25)**(-s))**(-1./s)  # analytic (not perfect)
    logthe = log10(the)
    if logthe > logthegrid[-2]:  # use analytic extrapolation
        return 10**(logTarr[-2] - 0.25 * (logthe - logthegrid[-2]))
    if logthe < logthegrid[0]:  # analytic extrapolation
        return 10**(logTarr[0] + 2 * (logthe - logthegrid[0]))
    ithe = floor((logthe - logthegrid[0])/dlogthe)
    slope = (logTarr[ithe+1] - logTarr[ithe])/dlogthe
    logT = logTarr[ithe] + (logthe - logthegrid[ithe]) * slope
    return 10**logT

def T_the(the, PhiL2, PhiRd, GM2, Rd, c2):   # only for non-zero fL2
    return (8*the**2 + 1 - 2 * (PhiL2-PhiRd)/(GM2/Rd))/(3*c2)


def fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd):  # only for non-zero fL2
    T = T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2)
    return 1 - c1 * T**4 * the**3/(the**2 - c2*T)

# compute T(Rsph, rho)
def func_T_Rsph_rho(T, Rsph, rho, mug, GM2):
    return constants.arad*T**4 + 1.5*rho*constants.kB*T/mug/constants.mp - GM2*rho/Rsph


# compute fL2 contribution from the inner disk outflow
def func_Rsph(lgRsph, Mdot, GM2, alpha_ss, mug, tol):
    Rsph = 10**lgRsph
    omgK_Rsph = sqrt(GM2/Rsph**3)
    rho = Mdot/(2*pi*alpha_ss*omgK_Rsph*Rsph**3)

    # T = (1./constants.arad * rho * GM2/Rsph)**0.25
    # get the temperature near Rsph by solving aT^4 + 1.5*rho*kB*T/mug/mp - GM2*rho/Rsph = 0
    Tleft, Tright = 0., 1e7   # this range should be wide enough
    fleft = func_T_Rsph_rho(T=Tleft, Rsph=Rsph, rho=rho, mug=mug, GM2=GM2)
    while abs((Tright-Tleft)/Tright) > tol:
        Tmid = 0.5 * (Tleft + Tright)
        fmid = func_T_Rsph_rho(T=Tmid, Rsph=Rsph, rho=rho, mug=mug, GM2=GM2)
        
        if fmid * fleft < 0:
            Tright = Tmid
        else:
            Tleft = Tmid
            fleft = fmid
    T = 0.5 * (Tleft + Tright)

    return Rsph - Mdot*kap(rho, T)/(4*pi*constants.c)
    # return Rsph - Mdot * 0.34 / (4 * pi * constants.c)   # simple Thomson opacity prescription

def solve_Rsph(Mdot, tol, GM2, alpha_ss, mug):
    lgRleft = 7.   # this should be small enough
    lgRright = 15.
    fleft = func_Rsph(lgRsph=lgRleft, Mdot=Mdot, GM2=GM2, alpha_ss=alpha_ss, mug=mug, tol=tol)

    # fleft should be negative
    while abs(lgRleft-lgRright) > tol:
        lgR = (lgRleft + lgRright)/2
        fmid = func_Rsph(lgRsph=lgR, Mdot=Mdot, GM2=GM2, alpha_ss=alpha_ss, mug=mug, tol=tol)
        if fleft * fmid < 0:
            lgRright = lgR
        else:
            lgRleft = lgR
            fleft = fmid

    lgR = (lgRleft + lgRright)/2

    return 10**lgR



def calculate_inner_disk_contribution(NM1dot, logM1dotgrid, Na, logagrid, Rd_over_a, PhiL2_dimless, M2, mu, solu_fL2, solu_fL2_inner, GM2, tol, alpha_ss, mug):
    """
    Function to calculate and set the inner disk L2 mass loss contribution
    """

    for n1 in range(NM1dot):
        M1dot = 10**logM1dotgrid[n1]*constants.msun/constants.yr
        for n2 in range(Na):
            a = 10**logagrid[n2]*constants.rsun
            # auxiliary parameters
            Rd = Rd_over_a * a
            PhiL2 = PhiL2_dimless * constants.G*(M2/mu)/a
            fL2_out = solu_fL2[n1, n2]
            # solve for Rsph
            Rsph = min(Rd, solve_Rsph(Mdot=(1-fL2_out)*M1dot, tol=tol, GM2=GM2, alpha_ss=alpha_ss, mug=mug))
            solu_fL2_inner[n1, n2] = (1-fL2_out) * abs(PhiL2)*Rsph/GM2

    return solu_fL2_inner

def run_thickness_bisection_using_f2_for_fl2ne0(theleft, f2left, theright, tol, PhiL2, PhiRd, GM2, Rd, c2, c1, c4, c3, PhiL1, alpha_ss, M1dot, omgK, lgRgrid, intp_lgkapgrid):
    """
    Function to run the thickness bisection for fl2!=0
    """

    while abs((theleft-theright)/theright) > tol:
        the = (theleft + theright)/2
        f2 = f2_the_T_fL2(the=the, T=T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2), fL2=fL2_the(the=the, c1=c1, c2=c2, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd), c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid)

        if f2 * f2left > 0:
            theleft = the
            f2left = f2
        else:
            theright = the

    return theleft, theright

import numpy as np

def bisect_temperature(c1, c2, Nthe, thegrid, tol):
    """
    Function to do the temperature bisection
    """

    # only T < Tmax is possible
    Tmax = (4./(27*c1**2*c2))**(1./9)

    Tarr = np.zeros(Nthe, dtype=float)
    for i in range(Nthe):
        the = thegrid[i]

        # use bisection method
        Tleft = 0.1 * min(the**2 / c2, Tmax)
        f1left = f1_the_T_fL2(the=the, T=Tleft, fL2=0, c1=c1, c2=c2)
        Tright = Tmax
        while abs((Tleft-Tright)/Tright) > tol:
            T = (Tleft + Tright)/2
            f1 = f1_the_T_fL2(the=the, T=T, fL2=0, c1=c1, c2=c2)
            if f1 * f1left > 0:
                Tleft = T
                f1left = f1
            else:
                Tright = T
        Tarr[i] = (Tleft + Tright)/2

    return Tarr


def calculate_function_constants(alpha_ss, Rd, omgK, M1dot, GM2, mug):
    """
    Function to calculate all the fucntion constant
    """

    c1 = 2*pi * constants.arad * alpha_ss * Rd / (3 * omgK * M1dot)
    c2 = constants.kB * Rd / (GM2 * mug * constants.mp)
    c3 = 8*pi**2 * constants.arad * alpha_ss * constants.c * Rd**2 / (M1dot**2 * omgK)
    c4 = 2*pi * mug * constants.arad * alpha_ss * omgK * constants.mp * Rd**3 / (constants.kB * M1dot)

    return c1, c2, c3, c4

def unnormalise_values(logagrid, n2, Rd_over_a, M2, mu, PhiL1_dimless, PhiL2_dimless, PhiRd_dimless):
    """
    Function to unnormalise some values
    """

    # Un-normalise values
    a = 10**logagrid[n2] * constants.rsun
    Rd = Rd_over_a * a
    Phi_units = constants.G*(M2/mu)/a
    PhiL1 = PhiL1_dimless * Phi_units
    PhiL2 = PhiL2_dimless * Phi_units
    PhiRd = PhiRd_dimless * Phi_units

    return a, Rd, Phi_units, PhiL1, PhiL2, PhiRd


def run_thickness_bisection_using_T_the_mins_T_the_nofL2(theleft, fleft, theright, tol, PhiL2, PhiRd, GM2, Rd, c2, logthegrid, logTarr, dlogthe):
    """
    Function to run the bisection to get the maximum corresponding to fl2=0
    """


    while abs((theleft - theright)/theright) > tol:
        the = (theleft + theright)/2
        f = T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2) - T_the_nofL2(the=the, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe)

        if f * fleft > 0:
            theleft = the
            fleft = f
        else:
            theright = the

    return theleft, theright