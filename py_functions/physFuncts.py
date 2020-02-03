"""
Functions to calculate physical properties
==============================

"""

import numpy as np
# from __future__ import print_function

def calcAirDensity(temperature, pressure):

    """
    Function to calculate air density from temperature and pressure
    ==============================

    """

        #### EXAMPLE OF USE:
        #### data = calcAirDensity(data['temperature'][:], data['pressure'][:])

    R = 2.8704  #### hPa kg-1 K-1

    print('Calculating air density profile:')
    print('')
    rho = np.zeros([np.size(temperature)])
    for k in range(0,np.size(temperature)):
        rho[k] = pressure[k] / (R * temperature[k])

    # print rho

    return rho


def calcThetaE(temperature, pressure, q, time, height):

    """
    Function to calculate equivalent potential temperature
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd

    eps = Rd/Rv

    # e = rh .* svp(T)  # in hPa
    # m =0.622.*e ./ (p-e); %mixing ratio [g /g]



    print('Calculating theta:')
    theta = np.zeros([len(time),len(height)])
    for k in range(0,len(height)):
        theta[:,k] = temperature[:,k] * np.power(1e5 / pressure[:,k], (Rd/cp))
    print('...')

    print('Calculating theta_e:')
    thetaE = np.zeros([len(time),len(height)])
    for k in range(0,len(height)):
        thetaE[:,k] = theta[:,k] + ((theta[:,k] * L_vap * q[:,k]) / (cp * temperature[:,k]))    ### Stull 1988[4] sect. 13.1 p. 546

    # %% Equation after Bryan 2008
    # % Konstants after Davies-Jones 2009: "On Formulas for equivalent pot. temperature"
    # % Accuracy 0.4 K
    # % T in degC
    # % p in hPa
    # % rh in %
    #
    # p0=1000;
    # Rd=287.04; % dry air J kg^-1 K^-1
    # Rv=461.50;
    #
    # Lvap = 2.555*10^6; %J kg^-1
    # cpd=1005.7; % J kg^1 K^-1
    # kd = Rd/cpd;
    # eps=Rd/Rv;
    #
    # e=rh .* svp(T)/100;  % in hPa
    # m =0.622.*e ./ (p-e); %mixing ratio [g /g]
    #
    # thetad=T.* ((p0./(p-e)).^kd); % in K
    # thetae = thetad .* (rh/100).^(-kd*m/eps) .* exp(L_vap*m./(T.*cpd));


    print('...')
    print('Done!')

    return theta, thetaE

def svp(T):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature = K

    """

    tempC = T - 273.15

    satvappres = 6.112 * np.exp( 17.67*tempC / (tempC + 243.5) ) * 100

    return satvappres

def polysvp(t,type):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature = K

    POLYSVP RETURNED IN UNITS OF PA.
    TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
    COPIED FROM PARCEL MODEL
    """

    ### ! ICE
    if type == 1:
        dt = np.maximum(-80.,t-273.16)
        polysvp = 6.11147274 + dt * (0.503160820 + dt * (0.188439774e-1 + dt * (0.420895665e-3 + dt *
            (0.615021634e-5 + dt * (0.602588177e-7 + dt * (0.385852041e-9 + dt * (0.146898966e-11 +
            0.252751365e-14 * dt)))))))
        polysvp = polysvp*100.

    ### ! LIQUID
    elif type == 0:
        dt = np.maximum(-80.,t-273.16)
        polysvp = 6.11239921 + dt * (0.443987641 + dt * (0.142986287e-1 + dt * (0.264847430e-3 + dt *
            (0.302950461e-5 + dt * (0.206739458e-7 + dt * (0.640689451e-10 + dt * (-0.952447341e-13 +
            -0.976195544e-15 * dt)))))))
        polysvp = polysvp*100.


    return polysvp
