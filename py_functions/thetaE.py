"""
Function to calculate equivalent potential temperature
==============================

"""

import numpy as np

def calcThetaE(data, time, height):

    L_vap = 2.5e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K

    print 'Calculating theta:'
    data['theta'] = {}
    data['theta'] = np.zeros([len(time),len(height)])
    for k in range(0,len(height)):
        data['theta'][:,k] = data['temperature'][:,k] * np.power(1e5 / data['pressure'][:,k], 0.2854)

    print 'Calculating theta_e:'
    data['theta_e'] = {}
    data['theta_e'] = np.zeros([len(time),len(height)])
    for k in range(0,len(height)):
        data['theta_e'][:,k] = data['theta'][:,k] + ((data['theta'][:,k] * L_vap * data['q'][:,k]) / (cp * data['temperature'][:,k]))

    return data
