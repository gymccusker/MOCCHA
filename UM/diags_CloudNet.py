###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import cartopy.crs as ccrs
import iris

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




# def writeNetCDF(cube):


def rotateGrid(data):
    ##
    ##      ROTATE GRID BACK TO STANDARD
    ##

    import iris.analysis.cartography as ircrt

    pole_lat = float(37.5000)       ### from rose suite
    pole_lon = float(177.5000)

    rot_lat = data.coord('grid_latitude').points
    rot_lon = data.coord('grid_longitude').points

    rot_pole = data.coord('grid_latitude').coord_system.as_cartopy_crs()

    # Transform real lat, lon point into rotated coords
    # rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()
    ll = ccrs.Geodetic()
    # rot_lon = np.zeros(np.size(lon,0))
    # rot_lat = np.zeros(np.size(lat,0))
    # for n in range(0,np.size(lon)-1):
    #     target_xy = rot_pole.transform_point(lon[n], lat[n], ll)  # lower left corner
    #     rot_lon[n] = target_xy[0] + 180.
    #     rot_lat[n] = target_xy[1] * float(-1.0)

    rot_lon, rot_lat = np.meshgrid(rot_lon, rot_lat)

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, pole_lon, pole_lat)

    # Print to check sites
    print '******'
    print 'Test to rotate coordinate grid: '
    print 'Rotated lon coord = ', rot_lon[0,0]
    print 'Rotated lat coord = ', rot_lat[0,0]
    print 'Lon = ', lon[0,0]
    print 'Lat = ', lat[0,0]
    print ' '

    # Transform real lat, lon point into rotated coords
    # rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()

    return lon, lat

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    print GlobalStashList
    print GlobalStashList[0]

    return GlobalStashList

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def main():

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

    filename1 = root_dir + out_dir + 'umnsaa_pb000'

    print 'Reading in files: '
    print filename1
    print ' '

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print 'Make stash list for cube read in at ' + time.strftime("%c")
    print ' '
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load

    # -------------------------------------------------------------
    # Load cubes
    # -------------------------------------------------------------
    print 'Cubes read in at ' + time.strftime("%c")
    print ' '


    # cube = iris.load(filenames, global_con, callback)
    # cube = iris.load(filename1, global_con, callback)

    ## Set variable constraint (i.e. which variable to load in based on stash code)
    var_con = iris.AttributeConstraint(STASH='m01s16i222')
    cube1 = iris.load_cube(filename1, var_con)

    # cube = iris.load(filename1)

    print cube1 # lists all diagnostics in file

    rot_pole = cube1.coord('grid_latitude').coord_system.as_cartopy_crs()

    lat, lon = rotateGrid(cube1

    # out = writeNetCDF(cube)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
