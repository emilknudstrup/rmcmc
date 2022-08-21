#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import subprocess
import itertools
import os

# =============================================================================
# Keplerian motion 
# =============================================================================
def solveKeplersEq(mean_anomaly, ecc, tolerance=1.e-5):
    '''Solves Kepler's equation.

    Function that solves Kepler's equation:
    .. :math:`M = E - \sin(E)`,

    where :math:`M` is the mean anomaly and :math:`E` the eccentric anomaly.

    This is done following the Newton-Raphson method as described in :cite:t:`Murray2010`.

    :param mean_anomaly: The mean anomaly.
    :type mean_anomaly: array
    :param ecc: Eccentricity.
    :type ecc: float
    :param tolerance: The tolerance for convergene. Defaults to 1.e-5.
    :type tolerance: float, optional

    :return: The new eccentric anomaly.
    :rtype: array 


    '''
    ## Circular orbit
    if ecc == 0: return mean_anomaly 

    new_ecc_anomaly = mean_anomaly
    converged = False

    for ii in range(300):
        old_ecc_anomaly = new_ecc_anomaly

        new_ecc_anomaly = old_ecc_anomaly - (old_ecc_anomaly - ecc*np.sin(old_ecc_anomaly) - mean_anomaly)/(1.0 - ecc*np.cos(old_ecc_anomaly))

        if np.max(np.abs(new_ecc_anomaly - old_ecc_anomaly)/old_ecc_anomaly) < tolerance:
            converged = True
            break

    if not converged:
        print('Calculation of the eccentric anomaly did not converge!')

    return new_ecc_anomaly

def trueAnomaly(time, Tw, ecc, per, ww, T0=True):
    '''Function that returns the true anomaly.

    The approach follows :cite:t:`Murray2010`.


    :param time: Times of observations.
    :type time: array
    :param Tw: Time of periastron.
    :type Tw: float
    :param ecc: Eccentricity.
    :type ecc: float
    :param per: Orbital period.
    :type per: float
    :param ww: Argument of periastron in radians.
    :type ww: float
    :param T0: Whether the reference time is periastron or mid-transit time. Default is mid-transit time.
    :param T0: bool

    :return: cosine, sine of the true anomaly.
    :rtype: (array, array)


    '''

    n = 2.0*np.pi/per
    ## With this you supply the mid-transit time 
    ## and then the time of periastron is calculated
    ## from S. R. Kane et al. (2009), PASP, 121, 886. DOI: 10.1086/648564
    if T0:
        f_Tw = np.pi/2.0 - ww
        E = 2.0*np.arctan(np.tan(f_Tw/2.0)*np.sqrt((1.0 - ecc)/(1.0 + ecc)))
        M = E - ecc*np.sin(E)
        Tc = Tw - M/n
    else:
        Tc = Tw
    mean_anomaly = n*(time-Tc)
    ecc_anomaly = solveKeplersEq(mean_anomaly,ecc)

    cos_E = np.cos(ecc_anomaly)
    sin_E = np.sin(ecc_anomaly)

    ## Cosine and sine of the true anomaly
    cos_f = (cos_E - ecc)/(1.0 - ecc*cos_E)
    sin_f = (np.sqrt(1 - ecc**2)*sin_E)/(1.0 - ecc*cos_E)

    return cos_f, sin_f

# =============================================================================
# Sky projected distance 
# =============================================================================
def projDist(cos_f,sin_f,ww,inc,ar,ecc):
    '''The separation of the centers of the two orbiting objects.

    Function that returns the separation of the centers of the two orbiting objects.
    The approach follows :cite:t:`Kreidberg2015`.


    :param cos_f: cosine of the true anomaly
    :type cos_f: array
    :param sin_f: sine of the true anomaly
    :type sin_f: array            
    :param ww: Argument of periastron in radians.
    :type ww: float
    :param inc: Inclination in radians.
    :type inc: float
    :param ar: Semi-major axis in stellar radii.
    :type ar: float
    :param ecc: Eccentricity.
    :type ecc: float

    :return: separation of centers.
    :rtype: array


    '''

    nn = len(cos_f)
    sep = np.zeros(nn)
    for ii in range(nn):
        ## Huge value for separation to make sure not to model planet passing behind star
        ## NOTE: Expressions like sin(w + f) are expanded to stay clear of arctan
        if np.sin(inc)*(np.sin(ww)*cos_f[ii] + np.cos(ww)*sin_f[ii]) <= 0:
          sep[ii] = 1000.
        else:
            nom = ar*(1.0 - ecc**2)
            nom *= np.sqrt(1.0 - (np.sin(ww)*cos_f[ii] + np.cos(ww)*sin_f[ii])**2*np.sin(inc)**2)
            den = 1.0 + ecc*cos_f[ii]
            sep[ii] = nom/den

    return sep

# =============================================================================
# x,y-position on the stellar disk 
# =============================================================================
def xyPos(cos_f,sin_f,ecc,ww,ar,inc,lam):
    '''Position of planet on stellar disk.

    Function to calculate the position on the stellar disk.
    Stellar disk goes from 0 to 1 in x and y.

    :param cos_f: cosine of the true anomaly
    :type cos_f: array
    :param sin_f: sine of the true anomaly
    :type sin_f: array            
    :param ecc: Eccentricity.
    :type ecc: float
    :param ww: Argument of periastron in radians.
    :type ww: float
    :param ar: Semi-major axis in stellar radii.
    :type ar: float
    :param inc: Inclination in radians.
    :type inc: float
    :param lam: Projected obliquity in radians.
    :type lam: float

    :return: x,y position of planet on stellar disk.
    :rtype: (array, array)


    '''
    r = ar*(1.0 - ecc**2)/(1.0 + ecc*cos_f)
    f = np.arctan2(sin_f,cos_f)

    ## x and y are lists of the positions of the transitting planet on the stellar disk 
    ## normalized to stellar radius (using a/Rs), corresponding to each RV-point
    x_old = -1*r*np.cos(ww + f)
    y_old = -1*r*np.sin(ww + f)*np.cos(inc)

    ## Rotate our coordinate system, such that the projected obliquity becomes the new y-axis
    x = x_old*np.cos(lam) - y_old*np.sin(lam)
    y = x_old*np.sin(lam) + y_old*np.cos(lam)
    return x, y

# =============================================================================
# Rossiter-McLaughlin effect 
# =============================================================================
def getRM(cos_f,sin_f,ww,ecc,ar,inc,rp,c1,c2,lam,vsini,
    xi=3.,gamma=1.,zeta=1.0,alpha=0.,cos_is=0.0,
    mpath='.'):
    '''The Rossiter-McLaughlin effect

    Function to calculate the Rossiter-McLaughlin effect for transiting exoplanets.

    The approach follows :cite:t:`Hirano2011`.

    :param cos_f: cosine of the true anomaly
    :type cos_f: array
    :param sin_f: sine of the true anomaly
    :type sin_f: array            
    :param ww: Argument of periastron in radians.
    :type ww: float
    :param ecc: Eccentricity.
    :type ecc: float
    :param ar: Semi-major axis in stellar radii.
    :type ar: float
    :param inc: Inclination in radians.
    :type inc: float

    :param rp: Planet-to-star radius ratio.
    :type rp: float

    :param c1: Linear limb-darkening coefficient.
    :type c1: float

    :param c2: Quadratic limb-darkening coefficient.
    :type c2: float

    :param lam: Projected obliquity in radians.
    :type lam: float

    :param vsini: Projected stellar rotation in km/s.
    :type vsini: float

    :param xi: Micro-turbulence in km/s. Defaults to 3.
    :type xi: float, optional

    :param zeta: Macro-turbulence in km/s. Defaults to 1.0.
    :type zeta: float, optional

    :param gamma: Coefficient of differential rotation. Defaults to 1.
    :type gamma: float, optional

    :param alpha:  in km/s. Defaults to 0.
    :type alpha: float, optional

    :param cos_is: Cosine of stellar inclination. Defaults to 0.
    :type alpha: float, optional

    :param mpath: Path to the code by [3]. Defaults to './'.
    :type mpath: str, optional

    :return: The RM signal in m/s.
    :rtype: array       

    '''
    x, y = xyPos(cos_f,sin_f,ecc,ww,ar,inc,lam)

    try:
        nn = len(cos_f)
        ## Alternates x and y for Hiranos code
        xy = [str(j) for j in itertools.chain.from_iterable(itertools.zip_longest(x,y))]
    except TypeError:
        nn = 1
        xy = [str(x),str(y)]   

    ## Create list of input to subprocess
    run_input = [mpath+'/Hirano2011_RM.exe']
    pars = [c1,c2,vsini,rp,xi,gamma,zeta,alpha,cos_is,nn]
    for par in pars: run_input.append(str(par))    

    RM = subprocess.check_output(run_input + xy)
    

    RM = [float(k)*1000 for k in RM.split()]

    return RM

# =============================================================================
# Radial velocity curve 
# =============================================================================

def getRV(time, per, T0, ecc, w, K, RVsys, a, Rp, b, vsini, lam, c1=0.3, c2=0.2, zeta=1.0, xi=2.0,  RM=False,mpath=None):
    '''The radial velocity curve

    Function that returns the radial velocity curve of the orbit following :cite:t:`Murray2010`.
    If RM is set to True it will include the RM effect as implemented in :py:func:`get_RM`.

    :param time: Times of observations.
    :type time: array

    :param target: Orbital parameters from :py:class:`Target`.
    :type target: object

    :param RM: Whether to calculate the RM effect or not. Defaults to ``False``.
    :type RM: bool, optional

    :param mpath: Path to the code by :cite:t:`Hirano2011`. Defaults to ``None`` in which case ``os.path.dirname(__file__)``.
    :type mpath: str, optional

    :return: Radial velocity curve.
    :rtype: array

    '''


    ## Convert angle from degree to radians
    w *= np.pi/180.

    ## Get cosine and sine of the true anomaly
    cos_f, sin_f = trueAnomaly(time,T0,ecc,per,w)
    ## Radial velocity
    vr = K*(np.cos(w)*(ecc + cos_f) - np.sin(w)*sin_f)
    

    if RM:
        inc = np.arccos(b/a)
        #inc = np.pi/180.
        ## Convert angle from degree to radians
        lam *= np.pi/180.

        ## Sepration       
        sep = projDist(cos_f,sin_f,w,inc,a,ecc)

        idxs = []
        for idx, dd in enumerate(sep):
            if 1 + Rp < dd: 
                pass
            else: 
                idxs.append(idx)

        if not mpath:
            mpath = os.path.abspath(os.path.dirname(__file__))

        if len(idxs) == 0:
            pass 
        elif len(idxs) == 1:
            cos_f, sin_f = np.array(cos_f[idx]), np.array(sin_f[idx])
            RMs = getRM(cos_f,sin_f,w,ecc,a,inc,Rp,c1,c2,lam,vsini,xi=xi,zeta=zeta,mpath=mpath)
            idx = idxs[0]
            vr[idx] = vr[idx] + RMs
        else:
            RMs = getRM(cos_f,sin_f,w,ecc,a,inc,Rp,c1,c2,lam,vsini,xi=xi,zeta=zeta,mpath=mpath)
            for idx in idxs: 
                vr[idx] = vr[idx] + RMs[idx]

    return vr + RVsys


# def getRV(time, target, RM=False,mpath=None):
#     '''The radial velocity curve

#     Function that returns the radial velocity curve of the orbit following :cite:t:`Murray2010`.
#     If RM is set to True it will include the RM effect as implemented in :py:func:`get_RM`.

#     :param time: Times of observations.
#     :type time: array

#     :param target: Orbital parameters from :py:class:`Target`.
#     :type target: object

#     :param RM: Whether to calculate the RM effect or not. Defaults to ``False``.
#     :type RM: bool, optional

#     :param mpath: Path to the code by :cite:t:`Hirano2011`. Defaults to ``None`` in which case ``os.path.dirname(__file__)``.
#     :type mpath: str, optional

#     :return: Radial velocity curve.
#     :rtype: array

#     '''
#     T0 = target.T0
#     ecc = target.ecc
#     per = target.per
#     w = target.w
#     K = target.K
#     RVsys = target.RVsys

#     ## Convert angle from degree to radians
#     w *= np.pi/180.

#     ## Get cosine and sine of the true anomaly
#     cos_f, sin_f = trueAnomaly(time,T0,ecc,per,w)
#     ## Radial velocity
#     vr = K*(np.cos(w)*(ecc + cos_f) - np.sin(w)*sin_f)
    

#     if RM:
#         a = target.a
#         inc = target.inc
#         ## Convert angle from degree to radians
#         inc *= np.pi/180.

#         Rp = target.rp
#         sep = projDist(cos_f,sin_f,w,inc,a,ecc)

#         ## Projected stellar rotation
#         vsini = target.vsini
#         ## Macroturbulence 
#         zeta = target.zeta
#         ## Microturbulence 
#         xi = target.xi
#         #beta = 3
#         ## The all important obliquity
#         lam = target.lam
#         lam *= np.pi/180.
#         ## LD coefficients
#         c1, c2 = target.c1, target.c2

#         idxs = []
#         for idx, dd in enumerate(sep):
#             if 1 + Rp < dd: 
#                 pass
#             else: 
#                 idxs.append(idx)

#         if not mpath:
#             mpath = os.path.abspath(os.path.dirname(__file__))

#         if len(idxs) == 0:
#             pass 
#         elif len(idxs) == 1:
#             cos_f, sin_f = np.array(cos_f[idx]), np.array(sin_f[idx])
#             RMs = getRM(cos_f,sin_f,w,ecc,a,inc,Rp,c1,c2,lam,vsini,xi=xi,zeta=zeta,mpath=mpath)
#             idx = idxs[0]
#             vr[idx] = vr[idx] + RMs
#         else:
#             RMs = getRM(cos_f,sin_f,w,ecc,a,inc,Rp,c1,c2,lam,vsini,xi=xi,zeta=zeta,mpath=mpath)
#             for idx in idxs: 
#                 vr[idx] = vr[idx] + RMs[idx]

#     return vr + RVsys

def totalDuration(per,rp,ar,inc,ecc,ww):
    '''The total duration of the transit, i.e., :math:`T_{41}`.

    This is the time from the first to the last contact point.

    :param per: Orbital period.
    :type per: float
    :param rp: Planet-to-star radius ratio.
    :type rp: float
    :param ar: Semi-major axis in stellar radii.
    :type ar: float
    :param inc: Inclination in radians.
    :type inc: float
    :param ecc: Eccentricity.
    :type ecc: float
    :param ww: Argument of periastron in radians.
    :type ww: float

    :return: the total duration of the transit
    :rtype: float

    .. note::
        The output will have the same units as the orbital period.


    '''
    b = ar*np.cos(inc)*(1 - ecc**2)/(1 + ecc*np.sin(ww))
    nom = np.arcsin( np.sqrt( ((1 + rp)**2 - b**2))/(np.sin(inc)*ar)  )*np.sqrt(1 - ecc**2)
    den = (1 + ecc*np.sin(ww))
    t41 = per/np.pi * nom/den
    return t41

def fullDuration(per,rp,ar,inc,ecc,ww):
    '''The duration of the transit with the planet completely within the stellar disk, i.e., :math:`T_{32}`.

    This is the time from the second to the third contact point.

    :param per: Orbital period.
    :type per: float
    :param rp: Planet-to-star radius ratio.
    :type rp: float
    :param ar: Semi-major axis in stellar radii.
    :type ar: float
    :param inc: Inclination in radians.
    :type inc: float
    :param ecc: Eccentricity.
    :type ecc: float
    :param ww: Argument of periastron in radians.
    :type ww: float

    :return: the total duration of the transit
    :rtype: float

    .. note::
        The output will have the same units as the orbital period.


    '''
    b = ar*np.cos(inc)*(1 - ecc**2)/(1 + ecc*np.sin(ww))
    nom = np.arcsin( np.sqrt( ((1 - rp)**2 - b**2))/(np.sin(inc)*ar)  )*np.sqrt(1 - ecc**2)
    den = (1 + ecc*np.sin(ww))
    t32 = per/np.pi * nom/den
    return t32