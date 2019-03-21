###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

import time
import datetime
import iris
import iris.plot as iplt
import numpy as np
from netCDF4 import Dataset
import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe


def main():

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'

    #### 	FILE NAMES
    #In [15]: ls
    #umnsaa_cb000  umnsaa_cb009  umnsaa_cb018  umnsaa_cb027     umnsaa_pverb000  umnsaa_pverc024
    #umnsaa_cb001  umnsaa_cb010  umnsaa_cb019  umnsaa_cb028     umnsaa_pverb006  umnsaa_pverd000
    #umnsaa_cb002  umnsaa_cb011  umnsaa_cb020  umnsaa_pa000     umnsaa_pverb012  umnsaa_pverd006
    #umnsaa_cb003  umnsaa_cb012  umnsaa_cb021  umnsaa_pb000     umnsaa_pverb018  umnsaa_pverd012
    #umnsaa_cb004  umnsaa_cb013  umnsaa_cb022  umnsaa_pvera000  umnsaa_pverb024  umnsaa_pverd018
    #umnsaa_cb005  umnsaa_cb014  umnsaa_cb023  umnsaa_pvera006  umnsaa_pverc000  umnsaa_pverd024
    #umnsaa_cb006  umnsaa_cb015  umnsaa_cb024  umnsaa_pvera012  umnsaa_pverc006  umnsa.stash
    #umnsaa_cb007  umnsaa_cb016  umnsaa_cb025  umnsaa_pvera018  umnsaa_pverc012  umnsa.xhist
    #umnsaa_cb008  umnsaa_cb017  umnsaa_cb026  umnsaa_pvera024  umnsaa_pverc018

    filename1 = root_dir + 'umnsaa_pb000'

    cube1 = iris.load(filename1)

    print cube1 # lists all diagnostics in file


if __name__ == '__main__':

    main()
