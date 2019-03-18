###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###


import iris
import numpy as np
from netCDF4 import Dataset
import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'

filename1 = root_dir + 'umnsaa_pvera000'

cube1 = iris.load(filename1)

print cube1 # lists all diagnostics in file
