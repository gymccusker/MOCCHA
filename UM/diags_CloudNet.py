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
from netCDF4 import Dataset as NetCDFFile
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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

def plot_basemap(ship_data, cube):

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    ##################################################
    ##################################################
    #### 	BASEMAP
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=SMALL_SIZE)
    plt.rc('ytick',labelsize=SMALL_SIZE)
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.rc('figure',titlesize=LARGE_SIZE)

    ## create figure and axes instances
    fig = plt.figure(figsize=(6,8))

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 2550000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=86,lat_0=86,lon_0=0)
    m.drawcoastlines()
    m.bluemarble()

    # define parallels/meridians
    m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,0,0,0],linewidth=0.8,fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],linewidth=0.8,fontsize=10)
    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='lightgrey')
    # m.fillcontinents(color='white')

    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)

    ### MAP ONTO PROJECTION
    x, y = m(ship_data.values[:,6], ship_data.values[:,7])
    x_inIcePeriod, y_inIcePeriod = m(ship_data.values[inIce_index,6],ship_data.values[inIce_index,7])
    x_driftPeriod, y_driftPeriod = m(ship_data.values[drift_index,6],ship_data.values[drift_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, color = 'yellow', linewidth = 2, label = 'Whole')
    plt.plot(x_inIcePeriod, y_inIcePeriod, color = 'darkorange', linewidth = 3, label = 'In Ice')
    plt.plot(x_driftPeriod, y_driftPeriod, color = 'red', linewidth = 4, label = 'Drift')

    ###########################################
    ### PLOT NEST + SWATH FOR INCREASED FREQ DIAGS VIS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629
        #
        # R.P.: original 1.5km nest -> (0, 86.625) @ 500x500

    ### SWATH
    lats = np.arange(80.9998,89.9998,0.09)
    lons = np.arange(3.0,76.0,0.73)
    x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)

    ### NEST (input)
    grx = float(5600)
    gry = float(700)
    centlon = float(39.5)
    centlat = float(85.275)
    latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    ### NEST (output)
    lato, lono = rotateGrid(cube)
    x1o, x2o, x3o, x4o, y1o, y2o, y3o, y4o = gridSetup(lono, lato, m)

    # draw swath
    pols =  Polygon([(x1s,y1s),(x2s,y2s),(x3s,y3s),(x4s,y4s)],\
                  facecolor='none',linestyle='--',edgecolor='g',linewidth=2,label='Swath')
    plt.gca().add_patch(pols)

    # draw nest (input)
    poln =  Polygon([(x1n,y1n),(x2n,y2n),(x3n,y3n),(x4n,y4n)],\
                  facecolor='none',linestyle='-',edgecolor='w',linewidth=2,label='Nest')
    plt.gca().add_patch(poln)

    # draw nest (output)
    polo =  Polygon([(x1o,y1o),(x2o,y2o),(x3o,y3o),(x4o,y4o)],\
                  facecolor='none',linestyle='--',edgecolor='y',linewidth=2,label='Nest (output)')
    plt.gca().add_patch(polo)

    ### ADD LEGEND
    plt.legend()

    plt.show()

def readfile(filename):

    import pandas as pd

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data, values

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = np.arange(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print ''
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print 'Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4])
    print ''

    return drift_index

def inIce(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    Sep_inIce = np.where(np.logical_and(data.values[:,2]<20,data.values[:,1]==9))
    inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    print ''
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ''

    return inIce_index

def gridSetup(lont, latt, m):

    lon, lat = np.meshgrid(lont, latt)

    nx = np.size(lon,0)
    ny = np.size(lat,1)
    x1, y1 = m(lon[nx-1,0],lat[nx-1,0])
    x2, y2 = m(lon[0,0],lat[0,0])
    x3, y3 = m(lon[0,ny-1],lat[0,ny-1])
    x4, y4 = m(lon[nx-1,ny-1],lat[nx-1,ny-1])

    return x1, x2, x3, x4, y1, y2, y3, y4

def rotateGrid(data):
    ##
    ##      ROTATE GRID BACK TO STANDARD
    ##

    import iris.analysis.cartography as ircrt

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 37.5000
    pole_lon = 177.5000

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

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, frst_lon, frst_lat)
    # lon, lat = np.meshgrid(lon, lat)

    # Print to check sites
    print '******'
    print 'Test to rotate coordinate grid: '
    print 'Rotated lon coord = ', rot_lon[0]
    print 'Rotated lat coord = ', rot_lat[0]
    print 'Lon = ', lon[0]
    print 'Lat = ', lat[0]
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

    # rot_pole = cube1.coord('grid_latitude').coord_system.as_cartopy_crs()

    ## LOCATION ON JASMIN
    ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'

    #### CODE NOT WORKING
    # seaice_rootdir = '/nfs/see-fs-02_users/eargy/MOCCHA/parent/data/seaice/AMSR2/'
    # seaice_file = 'asi-AMSR2-n6250-20180801-v5.hdf'
    # seaice_path = seaice_rootdir + seaice_file

    ship_data, values = readfile(ship_filename)

    columns = assignColumns(ship_data)

    # print ''
    # print data.head
    # print ''

    map = plot_basemap(ship_data, cube1)

    # out = writeNetCDF(cube)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
