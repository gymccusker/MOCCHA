"""
Functions to convert timestamps
==============================

"""

import time
import datetime
import numpy as np
import pandas as pd

def calcTime_Mat2DOY(cube):

        #### EXAMPLE OF USE:
        #### data = calcThetaE(data_um, time_um, height)

    print 'Converting MATLAB timesteps to DOY:'

    datenums = cube[0].dim_coords[0].points
    timestamps = pd.to_datetime(datenums-719529, unit='D')
    times = timestamps.dayofyear + (timestamps.hour / 24.0) + (timestamps.minute / 1440.0) + (timestamps.second / 86400.0)

    return times
