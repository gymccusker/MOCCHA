"""
Function to calculate air density from temperature and pressure
==============================

"""

import numpy as np

def calcAirDensity(temperature, pressure):

        #### EXAMPLE OF USE:
        #### data = calcAirDensity(data['temperature'][:], data['pressure'][:])

    R = 2.8704  #### hPa kg-1 K-1

    print 'Calculating air density profile:'
    rho = np.zeros([np.size(um.variables['temperature'],0),np.size(um.variables['temperature'],1)])
    for k in range(0,np.size(um.variables['temperature'],1)):
        rho[:,k] = pressure[:,k] / (R * temperature[:,k])

    return rho
