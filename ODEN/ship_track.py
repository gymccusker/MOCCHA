#####################################
##
## Read in and plot ship track
##      -- future release: relate to UM nest
##
#####################################

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime, date




def main():

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
    out_dir = '2_20180801_61DIAGS_TEST/'

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
