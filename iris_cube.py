###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###
#!/usr/bin/python2.7

import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
# import diags_MOCCHA as diags_ukv

def cart_plot(data1, data2):

    import iris.plot as iplt
    import matplotlib
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt
    import cartopy.crs as crs
    import cartopy.feature as cfe
    import iris.quickplot as qplt

    # Create a figure
    fig = plt.figure(figsize=(8,4))

    ## set axes position
    ax = fig.add_axes([0.1,0.1,0.35,0.8])	# left, bottom, width, height

    ## Allow interactive plotting
    # plt.interactive(True)

    ## Draw the contour with 25 levels.
    # contour = qplt.contourf(data[0,0,:,:], cmap = mpl_cm.Reds)
    # iplt.plot(data[0,:,0,0])
    # iplt.plot(data[0,:,200,200],data.aux_coords[2])
    ## iplt.plot(data[1,0:30,200,200],data.aux_coords[2][0:30])
    for i in range(0,np.size(data1,0)):
        strgi = "%1.f" % (i) # string of timestep
        iplt.plot(data1[i,0:30,200,200],data1.aux_coords[2][0:30],label=strgi)
    plt.legend()
    plt.title(data1.standard_name)

    # iplt.plot(data1[1,0:30,200,200],data1.aux_coords[2][0:30])

    ## set axes position
    ax = fig.add_axes([0.1,0.1,0.35,0.8])	# left, bottom, width, height
    for i in range(0,np.size(data1,0)):
        strgi = "%1.f" % (i) # string of timestep
        iplt.plot(data2[i,0:30,200,200],data2.aux_coords[2][0:30],label=strgi)
    plt.legend()
    plt.title(data2.standard_name)
    plt.show()


def main():

    import iris

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    # root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
    out_dir = '1_20160401_61DIAG_TEST/'

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

    filename1 = root_dir + out_dir + 'umnsaa_pb000'

    ## Set variable constraint (i.e. which variable to load in based on stash code)
    # var_con = iris.AttributeConstraint(STASH='m01s16i222')
    # cube1 = iris.load_cube(filename1, var_con)

    cube1 = iris.load(filename1)

    print cube1 # lists all diagnostics in file

    temperature = cube1[16]    # 3D air temperature, K
    Qvap = cube1[26]    # 3D specific humidity, kg/kg

    plot = cart_plot(temperature, Qvap)

if __name__ == '__main__':

    main()
