"""
File containing the function to run a grid of Mdot and separation configurations for a specific Mdonor and M accretor

TODO: pass settings to all the functions
"""

from math import pi, sqrt, log10
# from math import pi, sqrt, sin, cos, tan, log, log10, floor, ceil

from l2_massloss_interpolation.L2massloss_fork import constants

import numpy as np
from l2_massloss_interpolation.L2massloss_fork.functions.kappa_functions import set_up_kappa_interpolation_table, kap
from l2_massloss_interpolation.L2massloss_fork.functions.output_functions import write_data_to_txt_file, write_data_to_pickle_file
from l2_massloss_interpolation.L2massloss_fork.functions.functions import calculate_inner_disk_contribution, f1_the_T_fL2, f2_the_T_fL2, T_the_nofL2, T_the, fL2_the



def run_fl2_grid_for_gridpoint(settings):
    """
    Function to run a grid of Mdot and separation configurations for a specific Mdonor and M accretor
    """

    ######
    # Read out settings
    savedir = settings['savedir']
    nofL2 = settings['nofL2']

    M2_in_Msun = settings["mass_accretor"]
    q = settings["massratio_accretor_donor"]

    case = settings["kappa_case"] # 1 -- Hrich, solarZ; 2 -- Hpoor, solarZ; 3 -- Hrich, lowZ
    fdir = settings["kappa_fdir"] # kappa data directory
    extrap = settings["kappa_extrap"]   # extrapolation beyond grid may not be accurate [no need to]

    # grid for mass transfer rate and binary separation    # (100, 100) is sufficiently accurate
    NM1dot = settings["log10Mdot_donor_N"]
    logM1dotmin, logM1dotmax = settings["log10Mdot_donor_min"], settings["log10Mdot_donor_max"] # [Msun/yr]

    Na = settings["log10separation_N"]
    logamin, logamax = settings["log10separation_min"], settings["log10separation_max"] # [Rsun]

    # below is for grid search for disk thickness [only used at the beginning]
    Nthe = settings['disk_thickness_search_grid_N']   # ~50 is accurate enough
    thegrid_min, thegrid_max = settings['disk_thickness_search_grid_min'], settings['disk_thickness_search_grid_max']   # 0.1 to 1 is sufficient, we use analytic result below 0.1

    # other parameters (that are not part of the grid)
    alpha_ss = settings["alpha_ss"]   # viscosity
    tol = settings["tol"]   # fractional tolerance for bisection method
    eps_small = settings["eps_small"]   # a very small number
    Tfloor = settings["Tfloor"]   # [K] --- the minimum value for disk temperature solution

    ######
    # Set parameters based on the settings

    # Set savename
    savename = 'fL2grid_M%.1f_q%.1f_case%d' % (M2_in_Msun, q, case)

    # ---caution:
    # for very large binary separation a >> 1 AU and low Mdot << 1e-5 Msun/yr,
    # the disk temperature becomes much less than 1e4 K, then some of the opacities tables
    # are less reliable

    # Set grid of Mdot and separation
    logM1dotgrid = np.linspace(logM1dotmin, logM1dotmax, NM1dot)
    logagrid = np.linspace(logamin, logamax, Na)

    # below is for grid search for disk thickness [only used at the beginning]
    thegrid = np.logspace(log10(thegrid_min), log10(thegrid_max), Nthe, endpoint=True)

    # Set mean molecular weight
    if case in [1, 3]:   # mean molecular weight
        mug = 1.3/2.4   # H-rich
    else:
        mug = 4./3    # H-poor (fully ionized He)

    #
    lgq = log10(q)

    # TODO: make function to get the L1 and L2 positions
    # positions of Lagrangian points (based on analytic fits)
    xL1 = -0.0355 * lgq**2 + 0.251 * abs(lgq) + 0.500   # [a = SMA]
    xL2 = 0.0756 * lgq**2 - 0.424 * abs(lgq) + 1.699  # [a]
    if lgq > 0:   # m2 is more massive
        xL1 = 1 - xL1
        xL2 = 1 - xL2
    mu = q/(1 + q)

    # TODO: make function to get outer disk radius
    # outer disk radius
    Rd_over_a = (1-xL1)**4/mu

    # TODO: make function to get the dimensionless potential energies
    # relavent potential energies
    PhiL1_dimless = -((1-mu)/abs(xL1) + mu/abs(1-xL1) + 0.5*(xL1-mu)**2)   # [G(M1+M2)/a]
    PhiL2_dimless = -((1-mu)/abs(xL2) + mu/abs(1-xL2) + 0.5*(xL2-mu)**2)    # [G(M1+M2)/a]
    PhiRd_dimless = -(1 - mu + mu/Rd_over_a + 0.5 * (1-mu)**2)

    # TODO: consider doing this in the functions
    # recover CGS units
    M2 = M2_in_Msun * constants.msun
    GM2 = constants.G*M2

    # TODO: store in output structure to reduce arguments to all the functions
    # Set up structures for the solutions in 2D grids
    solu_the = np.empty((NM1dot, Na), dtype=float)  # outer disk thickness H/R
    solu_T = np.empty((NM1dot, Na), dtype=float)   # outer disk temperature
    solu_rho = np.empty((NM1dot, Na), dtype=float)
    solu_tau = np.empty((NM1dot, Na), dtype=float)
    solu_fL2 = np.empty((NM1dot, Na), dtype=float)      # outer disk
    solu_fL2_inner = np.empty((NM1dot, Na), dtype=float)      # inner disk
    Lacc_over_LEdd = np.empty((NM1dot, Na), dtype=float)
    QL2Qadv_over_Qrad = np.empty((NM1dot, Na), dtype=float)

    # TODO: put in dict to store all the info
    # Set up interpolation table
    lgRgrid, lgTgrid, lgkapgrid, intp_lgkapgrid = set_up_kappa_interpolation_table(case=case, fdir=fdir, extrap=extrap)

    ###########
    # Main loop
    percent = 5  # progress bar

    # Loop over Mdot
    for n1 in range(NM1dot):

        # Calculate percentage
        if 100*n1/NM1dot > percent:
            print('%d percent' % percent)
            percent += 5

        # Loop over separation
        M1dot = 10**logM1dotgrid[n1] * constants.msun/constants.yr
        for n2 in range(Na):
            # Un-normalise values
            a = 10**logagrid[n2] * constants.rsun
            Rd = Rd_over_a * a
            Phi_units = constants.G*(M2/mu)/a
            PhiL1 = PhiL1_dimless * Phi_units
            PhiL2 = PhiL2_dimless * Phi_units
            PhiRd = PhiRd_dimless * Phi_units

            # Keplerian frequency at Rd
            omgK = sqrt(GM2/Rd**3)

            # TODO: store the constants in a dictionary
            # TODO: make function call for this
            # constants involved in numerical solutions
            c1 = 2*pi * constants.arad * alpha_ss * Rd / (3 * omgK * M1dot)
            c2 = constants.kB * Rd / (GM2 * mug * constants.mp)
            c3 = 8*pi**2 * constants.arad * alpha_ss * constants.c * Rd**2 / (M1dot**2 * omgK)
            c4 = 2*pi * mug * constants.arad * alpha_ss * omgK * constants.mp * Rd**3 / (constants.kB * M1dot)

            # TODO: put the bisection into a function call
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

            # now we have obtained numerical relation between the and T
            logTarr = np.log10(Tarr)
            logthegrid = np.log10(thegrid)
            dlogthe = logthegrid[1] - logthegrid[0]

            # TODO: put in function
            # bisection to find the numerical solution to f2(the, T, fL2=0)=0
            theright = 1.
            f2right = f2_the_T_fL2(the=theright, T=T_the_nofL2(the=theright, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe), fL2=0, c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2)
            separation_factor = 0.95
            theleft = separation_factor * theright
            f2left = f2_the_T_fL2(the=theleft, T=T_the_nofL2(the=theleft, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe), fL2=0, c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2)
            while f2left*f2right > 0:  # need to decrease theleft
                theright = theleft
                f2right = f2left
                theleft *= separation_factor
                f2left = f2_the_T_fL2(the=theleft, T=T_the_nofL2(the=theleft, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe), fL2=0, c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2)

            # now the solution is between theleft and theright
            while abs((theleft-theright)/theright) > tol:
                the = (theleft + theright)/2
                f2 = f2_the_T_fL2(the=the, T=T_the_nofL2(the=the, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe), fL2=0, c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2)
                if f2 * f2left > 0:
                    theleft = the
                    f2left = f2
                else:
                    theright = the

            # solution
            the = (theleft + theright)/2
            T = T_the_nofL2(the=the, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe)

            rho = M1dot / (2*pi * alpha_ss * omgK * Rd**3) / the**3
            themax = sqrt(3./8 * c2 * T + 1./4 * (PhiL2-PhiRd)/(GM2/Rd) - 1./8)

            if the < themax or nofL2:   # this is the correct solution
                # Calculate relevant properties
                tau = rho*kap(rho=rho, T=T, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid)*Rd*the / 2
                Qrad = 2*pi*Rd**2*constants.arad*T**4*constants.c/tau
                U = constants.arad*T**4 + 1.5*rho*constants.kB*T/mug/constants.mp
                P = 1./3*constants.arad*T**4 + rho*constants.kB*T/mug/constants.mp
                tvis = 1./(alpha_ss * the**2 * omgK)
                Md = 2*pi*rho*Rd**3*the
                Qadv = 1.5*U/P * the**2 * GM2/Rd / tvis * Md
                QL2 = 0.

                # Store in result structures
                solu_rho[n1, n2] = rho
                solu_the[n1, n2] = the
                solu_T[n1, n2] = T
                solu_tau[n1, n2] = tau
                solu_fL2[n1, n2] = eps_small
                QL2Qadv_over_Qrad[n1, n2] = (QL2+Qadv)/Qrad
                Lacc_over_LEdd[n1, n2] = M1dot * kap(rho=rho, T=T, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid) / (4*pi * Rd * constants.c)

                continue

            # TODO: put in separate function
            # ---- below is for fL2 \neq 0

            themin = 1./2 * sqrt((PhiL2-PhiRd)/(GM2/Rd) - 1./2)   # corresponding to fL2=1, T=0

            # TODO put in function call
            # need to find the maximum corresponding to fL2=0
            # this is given by the intersection between T_the(the), T_the_nofL2(the)
            theleft = themin
            theright = 1.
            fleft = T_the(the=theleft, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2) - T_the_nofL2(the=theleft, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe)
            while abs((theleft - theright)/theright) > tol:
                the = (theleft + theright)/2
                f = T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2) - T_the_nofL2(the=the, logthegrid=logthegrid, logTarr=logTarr, dlogthe=dlogthe)

                if f * fleft > 0:
                    theleft = the
                    fleft = f
                else:
                    theright = the
            themax = (theleft + theright)/2   # this corresponds to fL2=0

            # --- numerical solution for f2(the, T, fL2)=0 under non-zero fL2

            # -- do not use exactly themin (corresponding to T = 0, bc. kap table breaks down)
            # -- define another themin based on Tfloor (kap table won't be a problem)
            themin = sqrt(3./8*c2*Tfloor + 1./4 * (PhiL2-PhiRd)/(GM2/Rd) - 1./8)
            theleft = themin  #
            f2left = f2_the_T_fL2(the=theleft, T=T_the(the=theleft, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2), fL2=fL2_the(the=theleft, c1=c1, c2=c2, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd), c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid)

            theright = themax / (1 + eps_small)
            # TODO: put in function call
            # bisection again
            while abs((theleft-theright)/theright) > tol:
                the = (theleft + theright)/2
                f2 = f2_the_T_fL2(the=the, T=T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2), fL2=fL2_the(the=the, c1=c1, c2=c2, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd), c4=c4, M1dot=M1dot, alpha_ss=alpha_ss, omgK=omgK, Rd=Rd, c3=c3, PhiRd=PhiRd, GM2=GM2, PhiL1=PhiL1, PhiL2=PhiL2)

                if f2 * f2left > 0:
                    theleft = the
                    f2left = f2
                else:
                    theright = the

            # Determine solution
            the = (theleft + theright)/2
            T = T_the(the=the, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd, c2=c2)
            fL2 = fL2_the(the=the, c1=c1, c2=c2, PhiL2=PhiL2, PhiRd=PhiRd, GM2=GM2, Rd=Rd)
            rho = (1-fL2) * M1dot / (2*pi * alpha_ss * omgK * Rd**3) / the**3

            # Calculate relevant properties
            # TODO: put in function or make function calls for these quantities
            tau = rho*kap(rho=rho, T=T, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid)*Rd*the / 2
            Qrad = 2*pi*Rd**2*constants.arad*T**4*constants.c/tau
            U = constants.arad*T**4 + 1.5*rho*constants.kB*T/mug/constants.mp
            P = 1./3*constants.arad*T**4 + rho*constants.kB*T/mug/constants.mp
            tvis = 1./(alpha_ss * the**2 * omgK)
            Md = 2*pi*rho*Rd**3*the
            Qadv = 1.5*U/P * the**2 * GM2/Rd / tvis * Md
            QL2 = fL2 * M1dot * (PhiL2 - PhiRd - 0.5*GM2/Rd)

            # Set results in output arrays
            # TODO: Put in function 
            QL2Qadv_over_Qrad[n1, n2] = (QL2+Qadv)/Qrad
            solu_rho[n1, n2] = rho
            solu_the[n1, n2] = the
            solu_T[n1, n2] = T
            solu_tau[n1, n2] = tau
            solu_fL2[n1, n2] = fL2
            Lacc_over_LEdd[n1, n2] = M1dot * kap(rho=rho, T=T, lgRgrid=lgRgrid, intp_lgkapgrid=intp_lgkapgrid) / (4*pi * Rd * constants.c)


    #######
    # below is for the contribution to L2 mass by the inner disk
    # fL2_inner is typically not significant and can be ignored for practical purposes
    solu_fL2_inner = calculate_inner_disk_contribution(
        NM1dot=NM1dot,
        logM1dotgrid=logM1dotgrid,
        Na=Na,
        logagrid=logagrid,
        Rd_over_a=Rd_over_a,
        PhiL2_dimless=PhiL2_dimless,
        M2=PhiL2_dimless,
        mu=mu,
        solu_fL2=solu_fL2,
        solu_fL2_inner=solu_fL2_inner,
        GM2=GM2,
        tol=tol,
        alpha_ss=alpha_ss,
        mug=mug
    )

    #######
    # Combine results to a total L2 mass loss fraction
    fL2_tot = solu_fL2 + solu_fL2_inner

    ##############
    # Output section

    # Write data
    write_data_to_pickle_file(
        savedir=savedir,
        savename=savename,
        logM1dotgrid=logM1dotgrid,
        logagrid=logagrid,
        solu_fL2=solu_fL2,
        solu_fL2_inner=solu_fL2_inner,
        solu_the=solu_the,
        solu_T=solu_T,
        Lacc_over_LEdd=Lacc_over_LEdd,
        QL2Qadv_over_Qrad=QL2Qadv_over_Qrad,
        solu_tau=solu_tau,
        solu_rho=solu_rho
    )
    exit()

    write_data_to_txt_file(
        savedir=savedir,
        solu_fL2=solu_fL2,
        solu_fL2_inner=solu_fL2_inner,
        solu_T=solu_T,
        Lacc_over_LEdd=Lacc_over_LEdd,
        QL2Qadv_over_Qrad=QL2Qadv_over_Qrad,
        logM1dotmin=logM1dotmin,
        logM1dotmax=logM1dotmax,
        NM1dot=NM1dot,
        logamin=logamin,
        logamax=logamax,
        M2_in_Msun=M2_in_Msun,
        q=q,
        Na=Na
    )