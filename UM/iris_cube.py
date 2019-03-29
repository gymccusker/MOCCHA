##
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
# import diags_MOCCHA as diags_ukv

def rotateGrid(data1):
    ##
    ##      ROTATE GRID BACK TO NORMAL
    ##

    rot_lat = data1.coord('grid_latitude').points
    rot_lat = data1.coord('grid_longitude').points

def vertical_profile(data1):

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
    plt.title(data1.standard_name + ', ' + str(data1.units))

    # iplt.plot(data1[1,0:30,200,200],data1.aux_coords[2][0:30])

    ## set axes position
    ax = fig.add_axes([0.5,0.1,0.35,0.8])	# left, bottom, width, height
    for i in range(0,np.size(data1,0)):
        strgi = "%1.f" % (i) # string of timestep
        iplt.plot(data2[i,0:30,200,200],data2.aux_coords[2][0:30],label=strgi)
    plt.legend()
    plt.title(data2.standard_name + ', ' + str(data2.units))

    plt.show()

def cart_map(data1, data2):

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
    ax = fig.add_axes([0.1,0.1,0.35,0.8]), projection=cart_proj)	# left, bottom, width, height
    # ax = plt.axes(projection=cart_proj)

    # Add coastlines
    ax.coastlines('50m', linewidth=0.8)

    # Plot contours
    # plt.contourf(wrf.to_np(lons), wrf.to_np(lats), data2, 10,
    #                 transform=crs.PlateCarree(), cmap = mpl_cm.Reds)

    ## Add a color bar
    # cbar = plt.colorbar(ax=ax, shrink=.62)
    # cbar.set_label(qnwfa2.name[-5:])

    # Set the map limits.  Not really necessary, but used for demonstration.
    # ax.set_xlim(wrf.cartopy_xlim(qnwfa2))
    # ax.set_ylim(wrf.cartopy_ylim(qnwfa2))

    # # Add the gridlines
    # ax.gridlines(color="black", linestyle="dotted")

    # plt.title(qnwfa2.name+'\n'+str(qnwfa2.Time.values))

    plt.show()

def main():

    import iris

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    # root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
    out_dir = '2_20180801_61DIAGS_TEST/'

    ## 1_20160401_61DIAG_TEST/
    ## 2_20180801_61DIAGS_TEST

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

    filename1 = root_dir + out_dir + 'umnsaa_pb024'

    ## Set variable constraint (i.e. which variable to load in based on stash code)
    # var_con = iris.AttributeConstraint(STASH='m01s16i222')
    # cube1 = iris.load_cube(filename1, var_con)

    cube = iris.load(filename1)

    print cube # lists all diagnostics in file

    temperature = cube[15]    # 3D air temperature, K
    Qvap = cube[25]    # 3D specific humidity, kg/kg

    # plot = vertical_profile(temperature, Qvap)
    plot = cart_map(temperature)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
