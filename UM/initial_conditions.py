###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import datetime
from iris.time import PartialDateTime

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')

# STASH_CODE = 'm01s04i118'


def readfile(filename):

    import pandas as pd

    # print '******'
    print ('')
    print ('Reading .txt file with pandas')
    print ('')

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
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print ('Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4]))
    print ('')

    return drift_index

def inIce(data):

    ###################################
    ## DEFINE IN ICE PERIOD
    ###################################
    # Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print ('******')
    print ('')
    print ('CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print ('')
    print ('Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print ('Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print ('Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print ('Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))
    print ('')

    return inIce_index

def readGriddedTrack(grid_filename):

    import pandas as pd

    print ('******')
    print ('')
    print ('Reading ' + grid_filename + ' file with pandas')
    print ('')

    data = pd.read_csv(grid_filename, sep = " ")
    values = data.values

    tim = values[:,1]
    ilon = values[:,2]
    ilat = values[:,3]

    return tim, ilat, ilon

def trackShip(data, date):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################

    print ('Extracting: ' + date)

    month = int(date[5])
    day = int(date[6:8])
    endpt = int(day)+2

    print ('month = ' + str(month))
    print ('start point = ' + str(day))
    print ('end point = ' + str(endpt))

    print (data)

    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==day,data.values[:,1]==month),data.values[:,3]>=12))
    trackShip_end = np.where(np.logical_and(data.values[:,2]==endpt ,data.values[:,1]==month))

    print (trackShip_start[0])
    print (trackShip_end[0])

    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][0])

    print ('******')
    print ('')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print ('Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print ('trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print ('')

    return trackShip_index

def plot_cartmap(ship_data, cube, date_dir):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting cartopy map:')
    print ('')

    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    diag = 0
    if len(cube[diag].dim_coords[-1].points) == 25:
    # if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif len(cube[diag].dim_coords[-1].points) == 56:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif len(cube[diag].dim_coords[-1].points) == 94:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -211
        yoffset = -385
    elif len(cube[diag].dim_coords[-1].points) == 81:          ### 14th and 24th August
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -209
        yoffset = -399
    elif len(cube[diag].dim_coords[-1].points) == 380:         ### needs checked
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -60
        yoffset = -110
    else:
    # elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

    ##################################################
    ##################################################
    #### 	CARTOPY
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


    #################################################################
    ## create figure and axes instances
    #################################################################
    i = 0
    for date in date_dir:
        i = i + 1
        plt.figure(figsize=(5,4))
        plt.subplots_adjust(top = 0.85, bottom = 0.15, right = 0.97, left = 0.1,
                hspace = 0.3, wspace = 0.3)
        # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))


        ### set size
        # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())      ### full track
        # ax.set_extent([20, 50, 88.35, 89.95], crs=ccrs.PlateCarree())     ### SWATH
        ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM

        ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

        #################################################################
        ## add geographic features/guides for reference
        #################################################################
        ax.add_feature(cartopy.feature.OCEAN, zorder=0)
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        # ax.set_global()
        ax.gridlines()

        #################################################################
        ## plot UM data
        #################################################################
        # iplt.pcolormesh(cube[date][0][0,:,:], vmin = 0.8, vmax = 1.0, cmap=mpl_cm.Blues_r)
        # plt.title(cube[date][0].standard_name)
        # plt.colorbar()

        #################################################################
        ## plot UM nest
        #################################################################
        ### draw outline of grid
        qplt.outline(cube[date][0][0,:,:])

        #################################################################
        ## plot ship track
        #################################################################
        ### DEFINE DRIFT + IN_ICE PERIODS
        drift_index = iceDrift(ship_data)
        inIce_index = inIce(ship_data)
        trackShip_index = trackShip(ship_data, date)

        ### Plot tracks as line plot
        plt.plot(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7],
                 color = 'darkorange', linewidth = 3,
                 transform = ccrs.PlateCarree(), label = date,
                 )
        plt.plot(ship_data.values[trackShip_index[0],6], ship_data.values[trackShip_index[0],7],
                 'k^', markerfacecolor = 'darkorange', linewidth = 3,
                 transform = ccrs.PlateCarree(),
                 )
        plt.plot(ship_data.values[trackShip_index[-1],6], ship_data.values[trackShip_index[-1],7],
                 'kv', markerfacecolor = 'darkorange', linewidth = 3,
                 transform = ccrs.PlateCarree(),
                 )

        #################################################################
        ## read in and plot gridded ship track
        #################################################################
        grid_dirname = 'AUX_DATA/'
        if int(date[6:8]) <= 8: grid_filename = grid_dirname + date[:6] + '0' + str(int(date[6:8])+1) + '_ShipTrack_GRIDDED.csv'
        if int(date[6:8]) >= 9: grid_filename = grid_dirname + date[:6] + str(int(date[6:8])+1) + '_ShipTrack_GRIDDED.csv'

        tim, ilat, ilon = readGriddedTrack(grid_filename)

        ### Plot tracks as line plot
        for i in range(0, len(ilon)-1):
            iplt.scatter(cube[date][0].dim_coords[2][int(ilon[i] + xoffset)], cube[date][0].dim_coords[1][int(ilat[i] + yoffset)],color='black')


        plt.legend()

        plt.savefig('../../FIGS/Grid_ZoomedTrack_' + date + '.png')
        plt.close()

    print ('******')
    print ('')
    print ('Finished plotting cartopy map! :)')
    print ('')



def plot_basemap(ship_data, cube):

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plot basemap:')
    print ('')

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
    dim = 2000000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=86,lat_0=86,lon_0=10)
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
    print ('******')
    print ('')
    print ('lat/lon vertices of proposed swath: ', lats[0], lats[-1], lons[0], lons[-1])
    print ('')
    x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)

    ### NEST (input)
    grx = float(500)
    gry = float(500)
    centlon = float(0.0)
    centlat = float(86.625)
    latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    print ('******')
    print ('')
    print ('lat/lon vertices of nest (input): ', latn[0], latn[-1], lonn[0], lonn[-1])
    print ('')
    x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    ### NEST (output)
    lono, lato = unrotateGrid(cube)
    print ('******')
    print ('')
    print ('lat/lon vertices of nest (output): ', lato[0], lato[-1], lono[0], lono[-1])
    print ('')
    x1o, x2o, x3o, x4o, y1o, y2o, y3o, y4o = gridSetup(lono, lato, m)

    ### NEST (required input to cover swath)
    # grx = float(5600)
    # gry = float(700)
    # centlon = float(39.5)
    # centlat = float(85.275)
    # latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    # lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    # print '******'
    # print ''
    # print 'lat/lon vertices of nest (input): ', latn[0], latn[-1], lonn[0], lonn[-1]
    # print ''
    # x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

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
                  facecolor='none',linestyle='-',edgecolor='y',linewidth=2,label='Nest (output)')
    plt.gca().add_patch(polo)

    # draw nest (output)
    # polot =  Polygon([(x1ot,y1ot),(x2ot,y2ot),(x3ot,y3ot),(x4ot,y4ot)],\
    #               facecolor='none',linestyle='--',edgecolor='m',linewidth=2,label='Test nest')
    # plt.gca().add_patch(polot)

    ### ADD LEGEND
    plt.legend()

    plt.show()

def gridSetup(lont, latt, m):

    lon, lat = np.meshgrid(lont, latt)

    nx = np.size(lon,0)
    ny = np.size(lat,1)
    x1, y1 = m(lon[nx-1,0],lat[nx-1,0])
    x2, y2 = m(lon[0,0],lat[0,0])
    x3, y3 = m(lon[0,ny-1],lat[0,ny-1])
    x4, y4 = m(lon[nx-1,ny-1],lat[nx-1,ny-1])

    return x1, x2, x3, x4, y1, y2, y3, y4

def unrotateGrid(data):
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

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, frst_lon, frst_lat)

    # Print to check conversion
    print ('******')
    print ('Test of unrotated coordinate grid: ')
    print ('Rotated lon coord = ', rot_lon[0])
    print ('Rotated lat coord = ', rot_lat[0])
    print ('Lon = ', lon[0])
    print ('Lat = ', lat[0])
    print (' ')

    # ******
    # lat/lon vertices of nest (output):  85.960219407715 80.41973098346767 49.567255645848604 -27.55740381723681
    # ******

    return lon, lat

def rotateGrid(data):
    ##
    ##      RELATE GRID TO ROTATED POLE
    ##

    import iris.analysis.cartography as ircrt

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 37.5000
    pole_lon = 177.5000

    lat = np.arange((centlat-(gry*float(0.5)*dlat)),(centlat+(gry*float(0.5)*dlat)),dlat)
    lon = np.arange((centlon-(grx*float(0.5)*dlon)),(centlon+(grx*float(0.5)*dlon)),dlon)

    rot_lon, rot_lat = ircrt.rotate_pole(lon, lat, frst_lon, frst_lat)

    # Print to check conversion
    print ('******')
    print ('Test of rotated coordinate grid: ')
    print ('Lon = ', lon[0])
    print ('Lat = ', lat[0])
    print ('Rotated lon coord = ', rot_lon[0])
    print ('Rotated lat coord = ', rot_lat[0])
    print (' ')

    # ******
    # lat/lon vertices of nest (output):  85.960219407715 80.41973098346767 49.567255645848604 -27.55740381723681
    # ******

    return lon, lat

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    return GlobalStashList

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def assignTimecoord(cube):
    '''
    assign Time as dimension coordinate if not present in cube
    '''
    if not cube.coords('time', dim_coords=True):
        time = cube.coord('time')
        time.bounds = None
        iris.util.promote_aux_coord_to_dim_coord(cube, time)

    return cube

def fixTimecoord(local_cube_list):

    import cf_units
    import iris.coords as icoords
    from datetime import datetime

    period_1 = 0

    time_ref = datetime(1970,1,1,0,0,0,0)
    time_start = datetime(2018,8,11,12,0,0,0)

    time_unit = cf_units.Unit(
        'hours since 1970-01-01', calendar=cf_units.CALENDAR_GREGORIAN)

    time.gmtime((cube.coord('time')[0].points)*3600)    ## TIME IN UTC

    # build ref Time coord
    LBYR = str(local_cube_list[period_1].attributes['LBYR'])
    if local_cube_list[period_1].attributes['LBMON'] >= 10:
        LBMON = str(local_cube_list[period_1].attributes['LBMON'])
    else:
        LBMON = '0' + \
            str(local_cube_list[period_1].attributes['LBMON'])
    if local_cube_list[period_1].attributes['LBDAT'] >= 10:
        LBDAT = str(local_cube_list[period_1].attributes['LBDAT'])
    else:
        LBDAT = '0' + \
            str(local_cube_list[period_1].attributes['LBDAT'])
    if local_cube_list[period_1].attributes['LBHR'] >= 10:
        LBHR = str(local_cube_list[period_1].attributes['LBHR'])
    else:
        LBHR = '0' + \
            str(local_cube_list[period_1].attributes['LBHR'])
    if local_cube_list[period_1].attributes['LBMIN'] > - 10:
        LBMIN = str(local_cube_list[period_1].attributes['LBMIN'])
    else:
        LBMIN = '0' + \
            str(local_cube_list[period_1].attributes['LBMIN'])

    scanTime = LBYR + '/' + LBMON + '/' + \
        LBDAT + ' ' + LBHR + ':' + LBMIN

    refTime = datetime.datetime.strptime(
        scanTime, "%Y/%m/%d %H:%M")

    refTimeCoord = icoords.AuxCoord(
        time_unit.date2num(refTime),
        standard_name='forecast_reference_time',
        units=time_unit)

    # Data Time coord
    LBYR = str(local_cube_list[period_2].attributes['LBYRD'])
    if local_cube_list[period_2].attributes['LBMOND'] >= 10:
        LBMON = str(local_cube_list[period_2].attributes['LBMOND'])
    else:
        LBMON = '0' + \
            str(local_cube_list[period_2].attributes['LBMOND'])
    if local_cube_list[period_2].attributes['LBDATD'] >= 10:
        LBDAT = str(local_cube_list[period_2].attributes['LBDATD'])
    else:
        LBDAT = '0' + \
            str(local_cube_list[period_2].attributes['LBDATD'])
    if local_cube_list[period_2].attributes['LBHRD'] >= 10:
        LBHR = str(local_cube_list[period_2].attributes['LBHRD'])
    else:
        LBHR = '0' + \
            str(local_cube_list[period_2].attributes['LBHRD'])
    if local_cube_list[period_2].attributes['LBMIND'] > - 10:
        LBMIN = str(local_cube_list[period_2].attributes['LBMIND'])
    else:
        LBMIN = '0' + \
            str(local_cube_list[period_2].attributes['LBMIND'])

    scanTime = LBYR + '/' + LBMON + '/' + \
        LBDAT + ' ' + LBHR + ':' + LBMIN

    dataTime = datetime.datetime.strptime(
        scanTime, "%Y/%m/%d %H:%M")

    dataTimeCoord = icoords.AuxCoord(time_unit.date2num(dataTime),
                                     standard_name='time',
                                     units=time_unit)

def testInput(cube):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Investigating input file:')
    print ('')

    print ('Cube attributes:')
    print (cube.attributes)
    print ('')
    print ('Cube aux coords:')
    print (cube.aux_coords)
    print ('')
    print ('Cube dim coords:')
    print (cube.dim_coords)
    print ('')
    print (np.size(cube.shape))

def write3DNetCDF(cube, outfile):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing NetCDF file:')
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC')

    print ('')
    print (dataset.file_format)
    print ('')

    ###################################
    ## Global Attributes
    ###################################
    desc = 'Test netCDF write out'
    micro = 'Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983)'
    dataset.description = desc
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'N/A'
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    dataset.comment = ''
    dataset.institution = 'University of Leeds.'

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    ###################################
    time = dataset.createDimension('time', np.size(cube.coord('time').points))
    # Z = dataset.createDimension('Z', np.size(icenum1,1))
    lat = dataset.createDimension('grid_latitude', np.size(cube.coord('grid_latitude').points))
    lon = dataset.createDimension('grid_longitude', np.size(cube.coord('grid_longitude').points))

    ###################################
    ## Dimensions variables
    ###################################
    #### Time
    time = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999')
    time.comment = cube.coord('time').standard_name
    time.units = str(cube.coord('time').units)
    time[:] = cube.aux_coords[1].points

    #### Latitude
    lat = dataset.createVariable('grid_latitude', np.float32, ('grid_latitude',),fill_value='-9999')
    lat.comment = cube.coord('grid_latitude').standard_name
    lat.units = str(cube.coord('grid_latitude').units)
    lat[:] = cube.coord('grid_latitude').points

    #### Longitude
    lon = dataset.createVariable('grid_longitude', np.float32, ('grid_longitude',),fill_value='-9999')
    lon.comment = cube.coord('grid_latitude').standard_name
    lon.units = str(cube.coord('grid_longitude').units)
    lon[:] = cube.coord('grid_longitude').points

    ###################################
    ###################################
    ## Create 3-d variables
    ###################################
    ###################################
    data = dataset.createVariable(cube.standard_name, np.float32, ('time','grid_latitude','grid_longitude',),fill_value='-9999')
    data.units = str(cube.units)
    # data.comment = cube.metadata
    data.attributes = str(cube.attributes)
    data[:] = cube.data[:]

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def write4DNetCDF(cube, outfile):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing NetCDF file:')
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC')

    print ('')
    print (dataset.file_format)
    print ('')

    ###################################
    ## Global Attributes
    ###################################
    desc = 'Test netCDF write out'
    micro = 'Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983)'
    dataset.description = desc
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'N/A'
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    dataset.comment = ''
    dataset.institution = 'University of Leeds.'

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    ###################################
    time = dataset.createDimension('time', np.size(cube.coord('time').points))
    # Z = dataset.createDimension('Z', np.size(icenum1,1))
    lat = dataset.createDimension('grid_latitude', np.size(cube.coord('grid_latitude').points))
    lon = dataset.createDimension('grid_longitude', np.size(cube.coord('grid_longitude').points))

    ###################################
    ## Dimensions variables
    ###################################
    #### Time
    time = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999')
    time.comment = cube.coord('time').standard_name
    time.units = str(cube.coord('time').units)
    time[:] = cube.aux_coords[1].points

    #### Latitude
    lat = dataset.createVariable('grid_latitude', np.float32, ('grid_latitude',),fill_value='-9999')
    lat.comment = cube.coord('grid_latitude').standard_name
    lat.units = str(cube.coord('grid_latitude').units)
    lat[:] = cube.coord('grid_latitude').points

    #### Longitude
    lon = dataset.createVariable('grid_longitude', np.float32, ('grid_longitude',),fill_value='-9999')
    lon.comment = cube.coord('grid_latitude').standard_name
    lon.units = str(cube.coord('grid_longitude').units)
    lon[:] = cube.coord('grid_longitude').points

    ###################################
    ###################################
    ## Create 3-d variables
    ###################################
    ###################################
    data = dataset.createVariable(cube.standard_name, np.float32, ('time','grid_latitude','grid_longitude',),fill_value='-9999')
    data.units = str(cube.units)
    # data.comment = cube.metadata
    data.attributes = str(cube.attributes)
    data[:] = cube.data[:]

    ###################################
    ###################################
    ## Create 4-d variables
    ###################################
    ###################################
    # nisg = dataset.createVariable('nisg', np.float32, ('time','Z','X','Y',),fill_value='-9999')
    # nisg.long_name = 'total ice number concentration'
    # nisg.comment = 'Sum of ice, snow, and graupel particles'
    # nisg.units = 'L-1'
    # nisg[:] = ice1[:]

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def readGlobal(cube, ship_data, date_dir):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy

    print ('******')
    print ('')
    print ('Defining longitude and latitude boundaries:')
    print ('')

    if np.ndim(cube[0].data) == 3:      ### 2D data + time
        lats = cube[0].dim_coords[1].points
        lons = cube[0].dim_coords[2].points
    if np.ndim(cube[0].data) == 4:      ### 3D data + time
        lats = cube[0].dim_coords[2].points
        lons = cube[0].dim_coords[3].points

    ###---------------------------------------------------------------------------------
    ### find northern and southern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    nb_lats = lats + ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    sb_lats = lats - ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    print ('sb_lats.shape = ' + str(sb_lats.shape))

    ###---------------------------------------------------------------------------------
    ### find western and eastern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    wb_lons = lons - ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    eb_lons = lons + ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    #################################################################
    ## find date of interest
    #################################################################
    date = date_dir[0:8]
    if date == '20180831':
        date = '20180901'
    else:
        date = date[:6] + str(int(date[-2:]) + 1).zfill(2)
    day_ind = np.array([])
    day_ind = np.where(np.logical_and(ship_data.values[:,2] == float(date[-2:]),ship_data.values[:,1] == float(date[-4:-2])))
    print ('Daily ship track for ' + date + ': ' + str(len(day_ind[0])) + ' pts ')

    #################################################################
    ## save daily lat/lons as temp vars
    #################################################################
    ship_lats = ship_data.values[day_ind[0],7]
    ship_lons = ship_data.values[day_ind[0],6]

    print ('ship_lats.shape = ' + str(ship_lats.shape))

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    ### compare hourly lat-lon with GLM grid
    data = {}
    data['ship_lons'] = np.zeros(24)
    data['ship_lats'] = np.zeros(24)
    data['ship_i'] = np.zeros(24); data['ship_i'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_j'] = np.zeros(24); data['ship_j'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_hour'] = np.zeros(24)
    hours = np.arange(0,24)
    jflag = np.zeros(24)        ### flag for if grid boundary is crossed

    for h in hours:
        print ('')
        print ('hour = ' + str(h))
        for j in range(0,len(sb_lats)):     ### for all latitude points
            if np.logical_and(ship_lats[h] >= sb_lats[j], ship_lats[h] < nb_lats[j]):     ### find where ship lat sits on glm lat grid
                for i in range(0,len(wb_lons)):     ### for all longitude points
                    if np.logical_and(ship_lons[h] >= wb_lons[i], ship_lons[h] < eb_lons[i]):     ### find where ship lon sits on glm lon grid
                        print ('lats and lons match at i = ' + str(i) + ', j = ' + str(j))
                        jflag[h] = jflag[h] + 1
                        data['ship_lons'][h] = lons[i]
                        data['ship_hour'][h] = hours[h]
                        data['ship_lats'][h] = lats[j]
                        data['ship_j'][h] = j         # define grid point indices for use later
                        data['ship_i'][h] = i         # define grid point indices for use later

    # print data['ship_lats']
    # print data['ship_j']
    # print data['ship_lons']
    # print data['ship_i']

    # np.save('working_glm_grid', data)

    ### need to constract an hour, lon index, lat index list like used for the lam
    # for h in hours:
    #     if data['ship_i'][h+1] == data['ship_i'][h]:

    ### arguments to be passed back to pullTrack_CloudNet
    tim = hours
    ilon = data['ship_i']
    ilat = data['ship_j']

    print (tim)
    print (ilon)
    print (ilat)

    # #####--------------------------------------------------------------------------------------------------
    # #####--------------------------------------------------------------------------------------------------
    # ##################################################
    # ##################################################
    # #### 	CARTOPY MAP
    # ##################################################
    # ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=SMALL_SIZE)
    # plt.rc('ytick',labelsize=SMALL_SIZE)
    # plt.rc('legend',fontsize=SMALL_SIZE)
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    #
    # #################################################################
    # ## create figure and axes instances
    # #################################################################
    # plt.figure(figsize=(6,8))
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))
    #
    # ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT
    #
    # #################################################################
    # ## add geographic features/guides for reference
    # #################################################################
    # ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    # ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    # ax.gridlines()
    #
    # #################################################################
    # ## plot global grid outline
    # #################################################################
    # ### draw outline of grid
    # # iplt.pcolormesh(cube[0][0,-10:-2,:-70])      ### covers whole drift
    # iplt.pcolormesh(cube[0][0,-7:-2,70:-70])      ### covers 28Aug - 4Sep subset of drift
    #
    # #################################################################
    # ## plot ship track
    # #################################################################
    # ### Plot tracks as line plot
    # plt.plot(ship_data.values[day_ind[0],6], ship_data.values[day_ind[0],7],
    #          color = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'Ship track',
    #          )
    # plt.plot(data['ship_lons'], data['ship_lats'],
    #          'o', color = 'yellow', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'Ship track',
    #          )
    #
    # plt.show()
    #
    # #####--------------------------------------------------------------------------------------------------

    return tim, ilat, ilon

def loadUMStartDump(filename):

    '''
    Loads in data from UM start dump used to drive global model
    Only loads in specific data relevant for calculating cloud + thermodynamic biases
    '''

    varlist = ['air_potential_temperature', 'sea_ice_thickness', 'mass_fraction_of_cloud_liquid_water_in_air',
                'mass_fraction_of_cloud_ice_in_air', 'specific_humidity']

    cube = iris.load(filename, varlist)

    return cube

def combinePP(root_dir, out_dir, date_dir):

    # -------------------------------------------------------------
    # Define output stream filenames to look at:
    #           start at 012 if 3h dumps (a, b)
    #           start at 011 if 1h dumps (c--e)
    # -------------------------------------------------------------
    print ('Checking if .pp files need to be combined:')
    print ('...')

    dummy = []

    names = ['_pa000','_pb000','_pc000','_pd000','_pe000']

    print ('out_dir is: ' + out_dir)
    print ('...')
    print ('')

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Define time and Stash constraints:
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Combine output pp streams
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    for date in date_dir:
        if date[0:4] == '2018':
        # if date[0:8] == '20180903':
            ## -------------------------------------------------------------------------
            ## Set fixed variable constraint (i.e. which variable to load in based on stash code)
            ## -------------------------------------------------------------------------
            # print '******'
            # print ''
            # print 'Use fixed constraint at ' + time.strftime("%c")
            # print ' '
            # var_con = iris.AttributeConstraint(STASH=STASH_CODE)

            ## -------------------------------------------------------------------------
            ## Time: constrain to take data every hour
            ## -------------------------------------------------------------------------
            # time_con = iris.Constraint(forecast_period=lambda xcell: any(np.isclose(xcell.point % 1, [0, 1./60.])))

            print ('******')
            print ('')
            print ('Naming output files: ')
            print ('')
            pp2_filename = str(int(date[:-6])+1) + '_oden_metum.pp'
            nc_filename = str(int(date[:-6])+1) + '_oden_metum.nc'
            print ('Output will be: ' + nc_filename)

            # print '******'
            print ('')
            print ('Identifying .pp files to read in: ')
            print ('')
            for stream in names:
                ### -------------------------------------------------------------------------
                ### define output filenames
                ### -------------------------------------------------------------------------
                # 20180827T1200Z_glm_pc035.pp
                filename1 = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream + '.pp'
                filename2 = root_dir + out_dir + date + '/' + date + '_glm' + stream + '.pp'
                filenames = [filename1, filename2]
                for filename in filenames:
                    print ('Checking: ' + filename)
                    if os.path.exists(filename):
                        pp_filename = filename[:-3] + '_r0.pp'

                        print ('---')
                        print ('Start files exist, continuing:')
                        print ('')

                        # print (filename[-12:-9])

                        ### define range to loop over
                        if filename[-12:-6] == 'glm_pa':
                            looping = range(0,6)
                            print ('Loop over 6 hours')
                        elif np.logical_or(stream[1:3] == 'pa', stream[1:3] == 'pb'):
                            looping = range(0,12)
                            print ('Loop over 3 hours')
                        else:
                            looping = range(0,36)
                            print ('Loop every hour')

                        for i in looping:
                            if np.size(looping) == 12:
                                res = i*3     # how many hourly dumps in file
                            elif filename[-12:-6] == 'glm_pa':
                                res = i*6
                            else:
                                res = i
                            str_i = "%03d" % res # file number
                            # fileout = root_dir + out_dir + date + stream[:-3] + str_i
                            if filename == filename1: fileout = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream[:-3] + str_i + '.pp'
                            if filename == filename2: fileout = root_dir + out_dir + date + '/' + date + '_glm' + stream[:-3] + str_i + '.pp'
                            # fileout = root_dir + out_dir + date + '/umnsaa_pa' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pb' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pc' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pd' + str_i
                            print (fileout)
                            # # -------------------------------------------------------------
                            # # Load cubes
                            # # -------------------------------------------------------------
                            print ('******')
                            print ('')
                            print ('Begin ' + str_i + ' cube read in at ' + time.strftime("%c"))
                            print (' ')

                            #### LOAD CUBE
                            if 'var_con' in locals():
                                cube = iris.load(fileout, var_con)

                                # -------------------------------------------------------------
                                # Write out data
                                # -------------------------------------------------------------
                                print ('******')
                                print ('')
                                print ('Outputting fixed constraint ' + str_i + ' data:')
                                print ('')
                                iris.save(cube, pp_filename, append=True)

                            elif 'global_con' in locals():
                                cube = iris.load(fileout, global_con, callback)

                                # -------------------------------------------------------------
                                # Write out data
                                # -------------------------------------------------------------
                                print ('******')
                                print ('')
                                print ('Outputting global constraint ' + str_i + ' data at ' + time.strftime("%c"))
                                print (cube)
                                print ('')
                                iris.save(cube, pp_filename, append=True)
                                #### remove file to keep directory tidy
                                print ('Directory clean up: removing ' + fileout)
                                print ('')
                                os.remove(fileout)
                    else:
                        print ('Combined output files already exist, or the directory does not exist')
                        print ('')

    return dummy

def loadPA(root_dir, out_dir, date_dir):

    '''
    Load in PA stream, combined pp file (suffix _r0.pp)
    '''
    model = ['_HighArctic_1p5km_','_glm']
    stream = '_pa000_r0'

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    def hourly_data(cell):
       # return True or False as to whether the cell in question should be kept
       # in this case, should return hourly data
       return cell >= 30


    # time_hourly_with_stash = iris.AttributeConstraint(STASH=lambda stash: str(stash) in GlobalStashList) & iris.Constraint(time=hourly_data)
    # time_hourly = iris.Constraint(time=hourly_data)

    cube = {}
    cubea = {}
    for date in date_dir:
        filename = root_dir + out_dir + date + '/' + date + model[0] + expt + stream + '.pp'

        cubea[date] = iris.load(filename, time_hourly)
        # cubea[date] = date
        #
        # for i in range(0, len(cube[date])):
        #     ### only load swath variables (size n=94)
        #      if np.size(cube[date][i].dim_coords[1],0) <= 100.:
        #          # print (cube[date][i].dim_coords[1])
        #          cubea[date] = cube[date][i]


    return cubea

def loadPD(root_dir, out_dir, date_dir):

    '''
    Load in PD stream, combined pp file (suffix _r0.pp)
    '''
    model = ['_HighArctic_1p5km_','_glm']
    stream = '_pd000_r0'

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    # def hourly_data(cell):
    #    # return True or False as to whether the cell in question should be kept
    #    # in this case, should return hourly data
    #    return cell >= 30


    # time_hourly_with_stash = iris.AttributeConstraint(STASH=lambda stash: str(stash) in GlobalStashList) & iris.Constraint(time=hourly_data)
    # time_hourly = iris.Constraint(time=hourly_data)

    # cube = {}
    cubec = {}
    for date in date_dir:
        filename = root_dir + out_dir + date + '/' + date + model[0] + expt + stream + '.pp'

        cubec[date] = iris.load(filename, 'sea_ice_area_fraction') #global_con, callback)
        # cubea[date] = date
        #
        # for i in range(0, len(cube[date])):
        #     ### only load swath variables (size n=94)
        #      if np.size(cube[date][i].dim_coords[1],0) <= 100.:
        #          # print (cube[date][i].dim_coords[1])
        #          cubea[date] = cube[date][i]


    return cubec

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'JASMIN'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/INITIAL_CONDITIONS_TEST/'
        init_dir = '/gws/nopw/j04/arcticcloud/MOCCHA/UM_STARTFILES/'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'

    ### CHOSEN RUN
    out_dir = '24_u-cc324_RA2T_CON/'
    out_dir_glm = '24_u-cc324_RA2T_CON/'
    date_dir = os.listdir(root_dir + out_dir)

    ## 4_u-bg610_RA2M_CON/              # Wilson and Ballard 1999 uphys
    ## 5_u-bl661_RA1M_CASIM/            # 100/cc accum mode aerosol; ARG + Cooper
    ## 6_u-bm410_RA1M_CASIM/            # 200/cc accum mode aerosol
    ## 7_u-bn068_RA2T_CON/              # RA2T_CON nest + global 4D stash
    ## 8_u-bp738_RA2M_CON/              # ERAI
    ## 10_u-bq791_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Fletcher
    ## 11_u-bq798_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Meyers
    ## 12_u-br210_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507
    ## 13_u-br409_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; passive aerosol processing
    ## 14_u-bu570_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit
    ## 15_u-bu687_RA2M_CON/           # Wilson and Ballard 1999 uphys; new RHcrit
    ## 16_u-bv926_RA2T_CON/              # RA2T_CON nest + global 4D stash + no subgrid mp production
    ## 17_u-bz429_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics
    ## 18_u-ca011_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics; 1A BL Scheme
    ## 19_u-ca012_RA2T_CON/              # RA2T_CON nest + global 4D stash; includes diagnosed turbulent dissipation rate
    ## 20_u-ca362_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; CICE sea ice albedo scheme
    ## 23_u-cc278_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM
    ## 24_u-cc324_RA2T_CON/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included
    ## 25_u-cc568_RA2M_CON/             # Wilson and Ballard 1999 uphys. sea ice albedo and extra BL diags
    ## 26_u-cd847_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. identical suite = u-cd852
    ## 27_u-ce112_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. passive aerosol processing.
    ## 28_u-ce627_RA2T_CON/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included. Mid-level convection switched off in GLM.
    ## 30_u-cg179_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; passive aerosol processing
    ## 31_u-cl349_RA2M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; jules fluxes

    # -------------------------------------------------------------
    # Extract from MASS with:
    # -------------------------------------------------------------
    #### moo select stash_extract_CloudNet.p2 moose:crum/u-bg610/apm.pp 4_u-bg610_RA2M_CON/20180811T1200Z/
    #### moo select stash_extract_CASIM.p5 moose:crum/u-bl661/apm.pp 5_u-bl661_RA1M_CASIM/20180831T1200Z/
    #### moo select stash_extract_CASIM.p5 moose:crum/u-bm410/apm.pp 6_u-bm410_RA1M_CASIM/20180905T1200Z/
    #### moo select stash_extract_all.p4 moose:crum/u-bn068/apm.pp 7_u-bn068_RA2T_CON/20180827T1200Z/
    #### moo select stash_extract_all.p4 moose:crum/u-bp738/apm.pp 8_u-bp738_RA2M_CON/20180904T1200Z/
    #### moo select stash_extract_CASIM_BL.p6 moose:crum/u-br210/apm.pp 12_u-br210_RA1M_CASIM/20180901T1200Z/
    #### moo select stash_extract_CASIM_BL.p6 moose:crum/u-br409/apm.pp 13_u-br409_RA1M_CASIM/20180901T1200Z/

    #### stash CASIM:
    ####    2, 3, 4, 10, 12, 24, 75, 78, 83, 86, 150, 254, 266, 267, 268, 271, 272, 273, 408, 409, 1201,
    ####    2201, 2391, 2392, 3025, 3217, 3234, 3236, 3245, 3247, 3248, 3360, 3361, 3476, 4118, 4203, 4204,
    ####    5216, 9203, 9204, 9205, 16004, 30461

    #### stash CASIM_BL:
    #### 2, 3, 4, 10, 12, 24, 26, 31, 75, 78, 83, 84, 88, 150, 254, 266, 267, 268, 271, 272, 273, 408, 409,
    #### 1201, 2201, 2391, 2392, 3002, 3025, 3208, 3217, 3219, 3220, 3223, 3234, 3236, 3245, 3247, 3248, 3360,
    #### 3361, 3362, 3363, 3460, 3461, 3464, 3465, 3469, 3471, 3473, 3476, 3501, 4118, 4203, 4204, 5216, 9203,
    #### 9204, 9205, 16004, 30461

    ####    RUN SCRIPT IN BACKGROUND (change to executable with chmod +x diags_CloudNet.py)
    #### module load jaspy
    #### nohup python diags_CloudNet.py > nohup_u-bz429_diags_CloudNet.out &

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data, values = readfile(ship_filename)
    columns = assignColumns(ship_data)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Combine PP files
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # dummy = combinePP(root_dir, out_dir, date_dir)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load combined PP files - pa
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    cubed = loadPD(root_dir, out_dir, date_dir)

    print (cubed)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Plot cartopy map
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    figure = plot_cartmap(ship_data, cubed, date_dir)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load in relevant UM start dump data
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------

    # ### list start dumps in UM_STARTFILES/
    # umdumps = os.listdir(init_dir)
    #
    # ### test out with first file
    # startdump = loadUMStartDump(umdumps[0])



    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')


if __name__ == '__main__':

    main()
