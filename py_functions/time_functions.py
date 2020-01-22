"""
Functions to convert timestamps
==============================

"""

import time
import datetime
import numpy as np
import pandas as pd

def calcTime_Mat2DOY(matlab_time):

        #### EXAMPLE OF USE:
        #### pytime = calcTime_Mat2DOY(matlab_time)

    print 'Converting MATLAB timesteps to DOY:'

    timestamps = pd.to_datetime(matlab_time-719529, unit='D')
    python_time = timestamps.dayofyear + (timestamps.hour / 24.0) + (timestamps.minute / 1440.0) + (timestamps.second / 86400.0)

    return python_time
