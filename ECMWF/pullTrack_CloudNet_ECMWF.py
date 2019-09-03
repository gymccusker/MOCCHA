###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE,
###         PULL SHIP TRACK, AND OUTPUT AS NETCDF FOR CLOUDNET
###

# from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
# import diags_MOCCHA as diags
# import diags_varnames as varnames
# import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os

def readfile(filename):

    import pandas as pd

    # print '******'
    print ''
    print 'Reading .txt file with pandas'
    print ''

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def pullLatLon(filename):

    from netCDF4 import Dataset

    print '*****'
    print 'Extracting lat/lon from ECMWF netCDF file'
    print ''

    nc = Dataset(filename,'r')

    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]
    time = nc.variables['time'][:]

    print 'ECMWF file at: (' + str(lon) + ', ' + str(lat) + ')'

    nc.close()

    return lat, lon, time

def designGrid(data):

    '''
    find eastern and northern boundaries of grid points
    '''

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

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
    # plt.figure(figsize=(12,10))
    # # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))
    #
    # ### set size
    # # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([0, 60, 87.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE


    ### make grid of unique latitude and longitude points
    lats, lons = np.meshgrid(data['ulat'][:], data['ulon'][:])

    # plt.scatter(lons, lats, c = 'lightgrey',#data['pressure'][:,0,0],
    #         label = 'Unique grid',
    #         transform = ccrs.PlateCarree())
    #
    # plt.scatter(data['lons'][29], data['ulat'][15], c = 'purple',#data['pressure'][:,0,0],
    #         label = 'testing',
    #         transform = ccrs.PlateCarree())

    ###---------------------------------------------------------------------------------
    ### find northern boundaries of gridpoints
    nblats = ((data['ulat'][1:] - data['ulat'][0:-1]) / 2.0) + data['ulat'][0:-1]       ## northern bounds for latitude
    data['nb_lats'] = np.zeros([np.size(data['lats'][:])])
    # print 'Northern boundary array has shape: ' + str(np.size(nblats))
    for j in range(0,len(nblats)):
        # print 'j = ' + str(j)
        for i in range(0,len(data['lats'][:])):
            # print 'i = ' + str(i)
            if data['ulat'][j] == data['lats'][i]:
                data['nb_lats'][i] = nblats[j]
    data['nb_lats'][-2:] = 90.0

    ###---------------------------------------------------------------------------------
    ### find southern boundaries of gridpoints
    data['sb_lats'] = np.zeros([np.size(data['lats'][:])])
    # for j in range(0,len(data['nb_lats'])-1):
    data['sb_lats'][:] = data['lats'][:] - (data['nb_lats'][:] - data['lats'][:])
    data['sb_lats'][-2:] = data['nb_lats'][-3]
    data['sb_lats'][26] = data['nb_lats'][19]       ### hard fix for 20180817 and 20180815(?) forecasts
    data['nb_lats'][21] = data['sb_lats'][29]# + 1.0
        # if data['ulat'][j] == data['lats'][i]:
        #     data['nb_lats'][i] = nblats[j]

    ###---------------------------------------------------------------------------------
    ### find eastern boundaries of gridpoints
    rblons = ((data['lons'][1:] - data['lons'][0:-1]) / 2.0) + data['lons'][0:-1]       ## RH bounds for longitude
    print 'rblons = ' + str(rblons.shape)
    print 'lons = ' + str(np.size(data['lons'][:]))
    data['rb_lons'] = np.zeros([np.size(data['lons'][:])])
    data['rb_lons'][0] = data['lons'][1]
    data['rb_lons'][1] = rblons[1]
    data['rb_lons'][2] = data['lons'][3]
    data['rb_lons'][3:9] = rblons[3:9]
    data['rb_lons'][9] = rblons[8]
    data['rb_lons'][10] = data['lons'][9]
    data['rb_lons'][11:12] = rblons[11:12]
    data['rb_lons'][12] = rblons[11]
    data['rb_lons'][13] = rblons[13]
    data['rb_lons'][14] = rblons[13]
    data['rb_lons'][15:17] = rblons[15:17]
    data['rb_lons'][17] = rblons[17] + 1.0
    data['rb_lons'][18] = data['lons'][17] + ((data['lons'][18] - data['lons'][17]) / 2.0)
    data['rb_lons'][19] = rblons[19]
    data['rb_lons'][20] = data['lons'][27]
    data['rb_lons'][21:27] = rblons[21:27]
    data['rb_lons'][27] = data['lons'][27] + (rblons[27] - data['lons'][27])/2.0
    data['rb_lons'][28] = rblons[28]
    data['rb_lons'][29] = rblons[29]
    data['rb_lons'][30:35] = rblons[30:35]
    data['rb_lons'][36] = data['lons'][36] + ((data['lons'][37] - data['lons'][36]) / 2.0) #rblons[34]
    data['rb_lons'][37] = data['lons'][33]#data['lons'][37] + (rblons[35] - data['lons'][3])/2.0

    ###---------------------------------------------------------------------------------
    ### find western boundaries of gridpoints
    data['lb_lons'] = np.zeros([np.size(data['lons'][:])])
    data['lb_lons'][0:1] = data['lons'][0] - (data['rb_lons'][0] - data['lons'][0])
    data['lb_lons'][1] = data['lons'][0]
    data['lb_lons'][2] = data['rb_lons'][1]
    data['lb_lons'][3] = rblons[2]
    data['lb_lons'][4] = data['lons'][3]
    # data['lb_lons'][3:5] = data['lons'][2:4]
    data['lb_lons'][5] = data['lons'][4]; data['rb_lons'][5] = data['lons'][5] + (data['lons'][5] - data['rb_lons'][5])
    data['lb_lons'][6] = data['lons'][6] - (data['rb_lons'][6] - data['lons'][6])
    data['lb_lons'][7] = data['lons'][6]; data['rb_lons'][7] = data['lons'][7] + (data['lons'][7] - data['rb_lons'][7])
    data['lb_lons'][8] = data['lons'][8] - (data['rb_lons'][8] - data['lons'][8])
    data['lb_lons'][9] = data['rb_lons'][9]; data['rb_lons'][9] = data['lons'][9] + (data['lons'][9] - data['rb_lons'][9])
    data['lb_lons'][10] = data['lons'][10] - (data['rb_lons'][10] - data['lons'][10])
    data['lb_lons'][11] = data['lons'][11] - (data['rb_lons'][11] - data['lons'][11])
    data['lb_lons'][12] = data['rb_lons'][12]; data['rb_lons'][12] = data['lons'][12] + (data['lons'][12] - data['rb_lons'][12])
    data['lb_lons'][13] = data['lons'][13] - (data['rb_lons'][13] - data['lons'][13])
    data['lb_lons'][14] = data['rb_lons'][14]; data['rb_lons'][14] = data['lons'][14] + (data['lons'][14] - data['rb_lons'][14])
    data['lb_lons'][15] = data['lons'][15] - (data['rb_lons'][15] - data['lons'][15])
    data['lb_lons'][16] = data['lons'][15]; data['rb_lons'][16] = data['lons'][16] + (data['lons'][16] - data['rb_lons'][16])
    data['lb_lons'][17] = data['lons'][17] - (data['rb_lons'][17] - data['lons'][17])
    data['lb_lons'][18] = data['lons'][17]; data['rb_lons'][18] = data['lons'][18] + (data['lons'][18] - data['rb_lons'][18])
    data['lb_lons'][19] = data['lons'][19] - (data['rb_lons'][19] - data['lons'][19])
    data['lb_lons'][20] = data['rb_lons'][19]
    data['lb_lons'][21] = data['lons'][18]#data['lons'][21] - (data['rb_lons'][21] - data['lons'][21])
    data['lb_lons'][22] = data['rb_lons'][21]; data['rb_lons'][22] = data['lons'][22] + (data['lons'][22] - data['rb_lons'][21])
    data['lb_lons'][23] = data['lons'][23] - (data['rb_lons'][23] - data['lons'][23])
    data['lb_lons'][24:28] = data['rb_lons'][23:27]
    data['lb_lons'][28] = data['lons'][28] - (data['rb_lons'][28] - data['lons'][28])
    data['lb_lons'][29] = data['rb_lons'][28]; data['rb_lons'][29] = data['lons'][29] + (data['lons'][29] - data['rb_lons'][28])
    data['lb_lons'][30] = data['lons'][30] - (data['rb_lons'][30] - data['lons'][30])
    data['lb_lons'][31:35] = data['rb_lons'][30:34]
    data['lb_lons'][35] = data['rb_lons'][34]; data['rb_lons'][35] = data['lons'][35] + (data['lons'][35] - data['rb_lons'][34])
    data['lb_lons'][36] = data['lons'][36] - (data['rb_lons'][36] - data['lons'][36])
    data['lb_lons'][37] = data['rb_lons'][36]#; data['rb_lons'][37] = data['lons'][37] + (data['lons'][37] - data['rb_lons'][36])

    print 'All boundaries loaded :)'

    return data

def checkLatLon(ship_data, date, data):

    print ''
    print 'Finding lat/lon of ship track'
    print '...'

    lats = data['lats'][:]
    lons = data['lons'][:]
    tims = data['tims'][:]

    #################################################################
    ## find date of interest
    #################################################################
    data['day_ind'] = np.array([])
    data['day_ind'] = np.where(np.logical_and(ship_data.values[:,2] == float(date[-2:]),ship_data.values[:,1] == float(date[-4:-2])))
    print 'Daily ship track: ' + str(len(data['day_ind'][0])) + ' pts '
    # print data['day_ind'][0]

    #################################################################
    ## print ship track coordinates
    #################################################################
    print 'Ship start (lon,lat): ' + str(ship_data.values[data['day_ind'][0][0],7]) + ', ' + str(ship_data.values[data['day_ind'][0][0],6])
    print 'Ship end (lon,lat): ' + str(ship_data.values[data['day_ind'][0][-1],7]) + ', ' + str(ship_data.values[data['day_ind'][0][-1],6])

    ship_lats = ship_data.values[data['day_ind'],7]
    ship_lons = ship_data.values[data['day_ind'],6]

    ### compare hourly lat-lon with ECMWF grid
    data['ship_lons'] = np.zeros(24)
    data['ship_lats'] = np.zeros(24)
    data['ship_ind'] = np.zeros(24)
    data['ship_ind'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_hour'] = np.zeros(24)
    hours = np.arange(0,24)
    jflag = np.zeros(24)
    for h in hours:
        # works for hour = 0
        print ''
        print 'hour = ' + str(h)
        for j in range(0,len(data['sb_lats'])):
            if np.logical_and(ship_lats[0][h] >= data['sb_lats'][j], ship_lats[0][h] < data['nb_lats'][j]):
                # print 'j=' + str(j)
                # temp = j
                # templat = lats[j]
                if np.logical_and(ship_lons[0][h] >= data['lb_lons'][j], ship_lons[0][h] < data['rb_lons'][j]):
                    print 'lats and lons match at j = ' + str(j)
                    jflag[h] = jflag[h] + 1
                    # print jflag[h]
                    # templon = lons[j]
                    data['ship_lons'][h] = lons[j]
                    data['ship_hour'][h] = hours[h]
                    data['ship_lats'][h] = lats[j]
                    data['ship_ind'][h] = j         # define grid point indices for use later

        if data['ship_hour'][h] == 0: data['ship_hour'][h] = h

        ### hard code h = 2 and 3 for 20180813 (out-of-grid)
        if data['date'][:] == '20180813':
            if h == 2:
                data['ship_ind'][h] = data['ship_ind'][1]
                data['ship_lons'][h] = data['ship_lons'][1]
                data['ship_lats'][h] = data['ship_lats'][1]
                # jflag[h] = 1
            elif h == 3:
                data['ship_ind'][h] = data['ship_ind'][1]
                data['ship_lons'][h] = data['ship_lons'][1]
                data['ship_lats'][h] = data['ship_lats'][1]
                # jflag[h] = 2

        ### want to increment jflag if h is different gpt than h-1
        if h > 0:
            if data['ship_lons'][h] != data['ship_lons'][h-1]:
                print 'lon changes between h = ' + str(h-1) + ' and h = ' + str(h)
                if jflag[h-1] != 2: jflag[h-1] = 2         # manually increment flag if not already done so
            if data['ship_lons'][h] == data['ship_lons'][h-1]:
                # print 'reset jflag at h = ' + str(h-1)
                if jflag[h-1] != 1: jflag[h-1] = 1         # manually increment flag if not already done so
            if data['ship_lats'][h] != data['ship_lats'][h-1]:
                print 'lat changes between h = ' + str(h-1) + ' and h = ' + str(h)
                if jflag[h-1] != 2: jflag[h-1] = 2         # manually increment flag if not already done so

        data['jflag'] = jflag       ### so jflag is written out for testing
        # print jflag
        # print ''
        ### want to compare two hourly points (i.e. between 0h and 1h where was the ship)
        if h > 0:
            if jflag[h] > jflag[h-1]:
                print 'Grid crossed at h = ' + str(h)
                grid_crossed = h

    ### hard code 20180813 jflags out of loop so as to not confuse algorithm
    if data['date'][:] == '20180813':
        data['jflag'][2] = 1
        data['jflag'][3] = 2

    print ''
    print 'Ship indices are:'
    print data['ship_ind']
    print 'Final jflag is: '
    print data['jflag']

    return data

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=13,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print 'Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4])
    print ''

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
    Aug_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8),data.values[:,3]>=0))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ''
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')'
    # print 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6]))
    # print 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7]))
    print ''

    return inIce_index

def trackShip(data, date):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    if date == '20180831':
        trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==int(date[-2:]),data.values[:,1]==int(date[-4:-2])),data.values[:,3]>=0))
        trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==1,data.values[:,1]==9),data.values[:,3]==1))
    else:
        trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==int(date[-2:]),data.values[:,1]==int(date[-4:-2])),data.values[:,3]>=0))
        trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==(int(date[-2:]) + 1),data.values[:,1]==int(date[-4:-2])),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print '******'
    print ''
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')'
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4])
    print ''

    return trackShip_index

def plot_basemap(ship_data, lats, lons, tim):

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plot basemap:'
    print ''

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
    fig = plt.figure(figsize=(8,10))

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 500000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=89,lat_0=89,lon_0=20)
    m.drawcoastlines()
    # m.bluemarble()

    # define parallels/meridians
    m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,1,0,0],color='k',linewidth=1.,fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,10.),labels=[0,0,0,1],color='k',linewidth=1.,fontsize=10)
    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='aqua')
    # m.fillcontinents(color='coral',lake_color='aqua')

    ### DEFINE DRIFT + IN_ICE PERIODS
    # drift_index = iceDrift(ship_data)
    # inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data)
    edgelats, edgelons = designGrid(lats, lons, tim)

    ### MAP ONTO PROJECTION
    x, y = m(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7])

    # Plot tracks as line plot
    # plt.plot(x, y, color = 'darkorange', linewidth = 2, label = 'Ship track')

    x_ecmwf, y_ecmwf = m(lons, lats)
    # Plot grid box centres as scatter plot
    plt.scatter(x_ecmwf, y_ecmwf, 10,
            color = 'blue', marker = 's',
            edgecolor = 'blue', linewidth = 2,
            label = 'ECMWF')

    x_t, y_t = m(lons, edgelats)
    # Plot grid box centres as scatter plot
    plt.scatter(x_t, y_t, 10,
            color = 'red', marker = '^',
            edgecolor = 'blue', linewidth = 2,
            label = 'ECMWF top edges')

    x_r, y_r = m(edgelons[edgelons>0], lats[edgelons>0])
    # Plot grid box centres as scatter plot
    plt.scatter(x_r, y_r, 10,
            color = 'green', marker = '>',
            edgecolor = 'blue', linewidth = 2,
            label = 'ECMWF right edges')

    ###########################################
    ### PLOT NEST + SWATH FOR INCREASED FREQ DIAGS VIS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629
        #
        # R.P.: original 1.5km nest -> (0, 86.625) @ 500x500

    ### ADD LEGEND
    plt.legend()

    plt.show()

def plot_cartmap(ship_data, data, date): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###################################
    ## CHOOSE DIAGNOSTIC
    ###################################
    # diag = 1
    print ''
    print 'Available diags are: ', data.keys()

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting cartopy map:'
    print ''

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
    # plt.rc('figure',titlesize=LARGE_SIZE)

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.figure(figsize=(12,10))
    # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    ax.set_extent([0, 60, 87.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE

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
    ################################################################
    # if np.size(cube[diag].data.shape) == 4:
    #     iplt.pcolormesh(cube[diag][hour,0,:,:])
    # elif np.size(cube[diag].data.shape) == 3:
    #     iplt.pcolormesh(cube[diag][hour,:,:])
    #     # iplt.pcolormesh(cube[hour,471:495,240:264])
    # elif np.size(cube[diag].data.shape) == 2:
    #     iplt.pcolormesh(cube[diag][:,:])
    # plt.title(cube[diag].standard_name + ', ' + str(cube[diag].units))
    # plt.colorbar()
    # plt.pcolormesh(data['lats'][:],data['lons'][:],data['pressure'][:,0,0])

    #################################################################
    ## design ECMWF grid
    #################################################################

    # data = designGrid(data)

    ###---------------------------------------------------------------------------------
    ### plot grid midpoints
    ###---------------------------------------------------------------------------------
    plt.scatter(data['lons'][:], data['lats'][:], c = 'black',#data['pressure'][:,0,0],
            label = 'Grid mid points',
            transform = ccrs.PlateCarree())

    ###---------------------------------------------------------------------------------
    ### plot latitude boundaries
    ###---------------------------------------------------------------------------------
    plt.scatter(data['lons'][:], data['sb_lats'][:], c = 'pink',
            label = 'southern bounds',
            transform = ccrs.PlateCarree())

    plt.scatter(data['lons'][:], data['nb_lats'][:], c = 'red',
            label = 'northern bounds',
            transform = ccrs.PlateCarree())

    ###---------------------------------------------------------------------------------
    ### plot longitude boundaries
    ###---------------------------------------------------------------------------------
    plt.scatter(data['rb_lons'][:], data['lats'][:], c = 'blue',
            label = 'eastern bounds',
            transform = ccrs.PlateCarree())
    plt.scatter(data['lb_lons'][:], data['lats'][:], c = 'purple',
            label = 'western bounds',
            transform = ccrs.PlateCarree())

    #################################################################
    ## plot gridded ship track
    #################################################################
    plt.scatter(data['ship_lons'][:], data['ship_lats'][:], c = 'yellow',
            label = 'gridded track',
            transform = ccrs.PlateCarree())

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data, date)

    ## Plot tracks as line plot
    # plt.plot(ship_data.values[:,6], ship_data.values[:,7],
    #          color = 'yellow', linewidth = 2,
    #          transform = ccrs.PlateCarree(), label = 'Whole',
    #          )
    # plt.plot(ship_data.values[inIce_index,6], ship_data.values[inIce_index,7],
    #          color = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'In Ice',
    #          )
    # plt.plot(ship_data.values[inIce_index[0],6], ship_data.values[inIce_index[0],7],
    #          'k^', markerfacecolor = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(),
    #          )
    # plt.plot(ship_data.values[inIce_index[-1],6], ship_data.values[inIce_index[-1],7],
    #          'kv', markerfacecolor = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(),
    #          )
    plt.plot(ship_data.values[drift_index,6], ship_data.values[drift_index,7],
             color = 'yellow', linewidth = 2,
             transform = ccrs.PlateCarree(), label = 'Drift',
             )

    ### Plot tracks as line plot
    plt.plot(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7],
             color = 'green', linewidth = 3,
             transform = ccrs.PlateCarree(), label = 'Ship track',
             )
    plt.plot(ship_data.values[trackShip_index[0],6], ship_data.values[trackShip_index[0],7],
             'k^', markerfacecolor = 'green', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )
    plt.plot(ship_data.values[trackShip_index[-1],6], ship_data.values[trackShip_index[-1],7],
             'kv', markerfacecolor = 'green', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )

    plt.legend()

    print '******'
    print ''
    print 'Finished plotting cartopy map! :)'
    print ''

    plt.savefig('FIGS/ECMWF_gridBoundaries_wTRACK.svg')
    plt.show()

    return data

def pullTrack(ship_data, data, date, outfile):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt

    print '******'
    print ''
    #################################################################
    ## design grid boundaries
    #################################################################
    # print '******'
    # print ''
    print 'Designing grid boundaries:'
    print ''
    data = designGrid(data)

    #################################################################
    ## check position of ship track
    #################################################################
    print '******'
    print ''
    print 'Gridding from ship lat/lon:'
    print ''

    data = checkLatLon(ship_data, date, data)
    np.save('working_data', data)

    #################################################################
    ## write out data
    #################################################################
    # print '******'
    # print ''
    # print 'Write out hourly gridded EC IFS data:'
    # print ''
    # out = writeNetCDF(data, date, outfile)

    # #################################################################
    # ## append metadata
    # #################################################################
    # print '******'
    # print ''
    # print 'Appending metadata:'
    # print ''
    # out = appendMetaNetCDF(outfile, date)

    return data

def readCube(name):

    ### LOOP OVER FILENAMES TO EXTRACT DIAGNOSTIC OVER ALL GRIDBOXES

    print 'Filename to load is: ' + name

    diag = 24

    data = {}
    dat = np.zeros([25,137])
    cube = iris.load(name)
    print 'Diag will be ' + cube[diag].var_name
    tims = cube[diag].dim_coords[0].points
    hgts = cube[35].data
    lats = cube[40].data
    lons = cube[41].data
    if np.sum(cube[diag].shape) > 24:        # if 2D diagnostic
        mlevs = cube[diag].dim_coords[1].points
    for t in range(len(tims)):
        dat[t,:] = tims[t]
        for k in range(np.size(hgts,1)):
            dat[:,k] = cube[diag].data[t,k]
    data[cube[diag].var_name] = dat
    data['lats'] = lats
    data['lons'] = lons
    data['tims'] = tims
    data['hgts'] = hgts
    data['mlevs'] = mlevs

    # print data.keys()

    return data, cube, diag

def readDaily(filenames, date):

    from iris.coords import DimCoord
    from iris.cube import Cube

    '''
     function to read in each lat/lon ECMWF IFS (netCDF) file with Iris then
     output required diagnostics for Cloudnet into a new netCDF
    '''

    i = -1
    data = {}
    data['date'] = date
    data['pressure'] = np.zeros([38,25,137])
    data['hgts'] = np.zeros([38,25,137])
    data['tims'] = np.zeros([25])
    data['lats'] = np.zeros([38])
    data['lons'] = np.zeros([38])
    data['mlevs'] = np.zeros([137])
    for name in filenames:
        i = i + 1
        print 'i = ' + str(i)
        dat, cube, diag = readCube(name)
        # print dat
        data['pressure'][i, :, :] = dat['pressure'][:, :]
        data['hgts'][i, :, :] = dat['hgts'][:, :]
        data['lats'][i] = dat['lats']
        data['lons'][i] = dat['lons']
    data['tims'][:] = dat['tims'][:]
    data['mlevs'][:] = dat['mlevs'][:]
    if date == '20180904':
        data['lons'][:] = np.array([36.66999817, 38.08000183, 41.54000092, 43.20000076, 45.        ,
               46.95999908, 45.        , 47.13999939, 45.        , 49.5       ,
               47.36999893, 50.        , 55.        , 52.93999863, 58.24000168,
               56.25      , 60.        , 57.86000061, 62.31000137, 30.        ,
               37.5       , 67.5       , 75.        ,  8.18000031, 16.36000061,
               24.54999924, 32.72999954, 40.90999985, 57.27000046, 65.44999695,
               9.        , 18.        , 27.        , 36.        , 45.        ,
               54.        , 20.        , 30.        ])
        data['lats'][:] = np.array([88.40000153, 88.47000122, 88.47000122, 88.54000092, 88.61000061,
               88.68000031, 88.75      , 88.81999969, 88.88999939, 88.88999939,
               88.95999908, 89.02999878, 89.02999878, 89.09999847, 89.09999847,
               89.16999817, 89.23999786, 89.30999756, 89.37999725, 89.44999695,
               89.44999695, 89.44999695, 89.44999695, 89.52999878, 89.52999878,
               89.52999878, 89.52999878, 89.52999878, 89.52999878, 89.52999878,
               89.59999847, 89.59999847, 89.59999847, 89.59999847, 89.59999847,
               89.59999847, 89.66999817, 89.66999817])

    #################################################################
    ## CREATE EMPTY CUBE
    #################################################################
    ncube = Cube(np.zeros([38,25,137]))

    data['ulat'] = np.zeros([np.size(np.unique(data['lats'][:]))])
    data['ulat'][:] = np.unique(data['lats'][:])
    data['ulon'] = np.zeros([np.size(np.unique(data['lons'][:]))])
    data['ulon'][:] = np.unique(data['lons'][:])
    mlats, mlons = np.meshgrid(data['ulat'][:], data['ulon'][:])

    ntime = DimCoord(data['tims'][:], var_name = 'time', standard_name = 'time', units='hours since ' + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' 00:00:00 +00:00')
    level = DimCoord(data['mlevs'][:], var_name = 'level', standard_name = 'model_level_number', units='m')
    lats = DimCoord(data['ulat'][:], var_name = 'latitude', standard_name = 'latitude', units='degrees_N')
    lons = DimCoord(data['ulon'][:], var_name = 'longitude', standard_name = 'longitude', units='degrees_E')
    # ncube = Cube(data['pressure'][:,:,:],
    #         dim_coords_and_dims=[(lats, 0), (ntime, 1), (level, 2)],
            # standard_name = cube.standard_name,
            # long_name = cube.long_name,
            # units = cube.units,
            # var_name = varname,
            # attributes = cube.attributes,
            # aux_coords_and_dims = None,
            # )

    # nc_outfile = date + '_oden_ecmwf_n38.nc'
    # iris.save(ncube, nc_outfile)

    ### write to combined netCDF file
    # data = writeNetCDF(nc_outfile, data, date, cube)

    ### append metadata to combined netCDF file
    # data = appendMetaNetCDF(nc_outfile, date)

    return data

def writeNetCDF(data, date, outfile):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt
    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    ##################################
    # Open new netCDF file
    ##################################

    dataset = Dataset(outfile, 'w')

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    #### load in a cube to define dimensions
    ###################################
    cube0 = iris.load('DATA/' + date + '_moccha_ecmwf_001.nc')

    ###################################
    ## Data dimensions
    ###################################
    timem = dataset.createDimension('time', np.size(cube0[0].dim_coords[0].points))
    level = dataset.createDimension('model_level_number', np.size(cube0[1].dim_coords[1].points))
    flevel = dataset.createDimension('model_flux_level', np.size(cube0[3].dim_coords[1].points))
    freq = dataset.createDimension('frequency', np.size(cube0[2].dim_coords[0].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('time', np.float64, ('time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since ' + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' 00:00:00 +00:00.'
    timem.units = 'hours'
    timem.long_name = 'hours_UTC'
    timem.standard_name = 'time'
    timem[:] = cube0[0].dim_coords[0].points

    #### model level
    level = dataset.createVariable('model_level_number', np.float64, ('model_level_number',), fill_value='-9999')
    level.scale_factor = float(1)
    level.add_offset = float(0)
    level.comment = ''
    level.units = '1'
    level.long_name = 'model_level'
    level.standard_name = 'model_level_number'
    level.positive = 'down'
    level[:] = cube0[1].dim_coords[1].points

    #### flux model level
    flevel = dataset.createVariable('model_flux_level', np.float64, ('model_flux_level',), fill_value='-9999')
    flevel.scale_factor = float(1)
    flevel.add_offset = float(0)
    flevel.comment = ''
    flevel.units = '1'
    flevel.long_name = 'model_flux_level'
    flevel.positive = 'down'
    flevel[:] = cube0[3].dim_coords[1].points

    #### frequency
    freq = dataset.createVariable('frequency', np.float64, ('frequency',), fill_value='-9999')
    freq.scale_factor = float(1)
    freq.add_offset = float(0)
    freq.comment = ''
    freq.units = 'GHz'
    freq.long_name = 'microwave_frequency'
    freq.missing_value = -999.0
    freq[:] = cube0[2].dim_coords[0].points


    ###################################
    ## Set fluxes to distinguish from model_level diags
    ###################################
    fluxes = ['flx_net_sw','flx_net_lw','flx_ls_snow','flx_ls_rain','flx_turb_mom_v',
    'flx_turb_mom_u','flx_conv_snow','flx_conv_rain','flx_height','flx_turb_moist','flx_down_sens_heat']

    qfields = ['qi','ql']

    ###################################
    ## Loop over diagnostics
    ###################################
    for d in range(0,len(cube0)):
        print ''
        print '**** New diagnostic loop ****'
        for h in range(0,24):
            #### increment ship_ind by 1 to change from python counting (from 0) to 1-38
            file = 'DATA/' + date + '_moccha_ecmwf_' + str(int(data['ship_ind'][h] + 1)).zfill(3) + '.nc'
            print ''
            print 'Ship index = ' + str(data['ship_ind'][h]) + ', h = ' + str(h) + ', and JFLAG = ' + str(data['jflag'][h])
            if data['jflag'][h] == 1:
                print 'So loading ' + file + ' only '
                cube = iris.load(file)
                print 'Writing ' + cube[d].var_name
                if h == 0:
                    if np.size(cube[d].shape) == 0:
                        print 'Diagnostic is a scalar field, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube[d].var_name, np.float64, fill_value='-9999')
                        dat[:] = cube[d].data
                        break
                    elif np.size(cube[d].shape) == 1:
                        print 'Diagnostic is 1D, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube[d].var_name, np.float64, ('time',), fill_value='-9999')
                        dat[h] = cube[d].data[h]
                    elif np.size(cube[d].shape) == 2:
                        print 'Diagnostic is 2D, so writing hour = ' + str(h)
                        if cube[d].var_name in fluxes:
                            print 'Diagnostic is on flux levels.'
                            dat = dataset.createVariable(cube[d].var_name, np.float64, ('time','model_flux_level',), fill_value='-9999')
                        else:
                            print 'Diagnostic is on model levels.'
                            dat = dataset.createVariable(cube[d].var_name, np.float64, ('time','model_level_number',), fill_value='-9999')
                        dat[h,:] = cube[d].data[h,:]
                    elif np.size(cube[d].shape) == 3:
                        print 'Diagnostic is 3D, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube[d].var_name, np.float64, ('frequency','time','model_level_number',), fill_value='-9999')
                        dat[:,h,:] = cube[d].data[:,h,:]
                else:
                    if np.size(cube[d].shape) == 0:
                        print 'Diagnostic is a scalar field, so writing hour = ' + str(h)
                        dat[:] = cube[d].data
                    elif np.size(cube[d].shape) == 1:
                        print 'Diagnostic is 1D, so writing hour = ' + str(h)
                        dat[h] = cube[d].data[h]
                    elif np.size(cube[d].shape) == 2:
                        print 'Diagnostic is 2D, so writing hour = ' + str(h)
                        dat[h,:] = cube[d].data[h,:]
                    elif np.size(cube[d].shape) == 3:
                        print 'Diagnostic is 3D, so writing hour = ' + str(h)
                        dat[:,h,:] = cube[d].data[:,h,:]
            elif data['jflag'][h] == 2:
                print 'h = ' + str(int(data['ship_ind'][h]))
                print 'h+1 = ' + str(int(data['ship_ind'][h+1]))
                file2 = 'DATA/' + date + '_moccha_ecmwf_' + str(int(data['ship_ind'][h+1]) + 1).zfill(3) + '.nc'
                print 'So averaging between ' + file + ' and ' + file2
                print 'Loading ' + file + '...'
                cube1 = iris.load(file)
                if cube1[d].var_name in qfields:
                    print 'Diagnostic is a Q-field, so setting all zeros to NaNs before averaging...'
                    cube1[d].data[cube1[d].data == 0.0] = np.nan
                print 'Loading ' + file2 + '...'
                cube2 = iris.load(file2)
                if cube2[d].var_name in qfields:
                    print 'Diagnostic is a Q-field, so setting all zeros to NaNs before averaging...'
                    cube2[d].data[cube2[d].data == 0.0] = np.nan
                print 'Averaging ' + cube1[d].var_name
                if h == 0:
                    if np.size(cube1[d].shape) == 0:
                        print 'Diagnostic is a scalar field, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube1[d].var_name, np.float64, fill_value='-9999')
                        dat[:] = np.nanmean([cube1[d].data,cube2[d].data])
                        break
                    elif np.size(cube1[d].shape) == 1:
                        print 'Diagnostic is 1D, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube1[d].var_name, np.float64, ('time',), fill_value='-9999')
                        dat[h] = np.nanmean([cube1[d].data[h],cube2[d].data[h]])
                    elif np.size(cube1[d].shape) == 2:
                        print 'Diagnostic is 2D, so writing hour = ' + str(h)
                        if cube1[d].var_name in fluxes:
                            print 'Diagnostic is on flux levels.'
                            dat = dataset.createVariable(cube1[d].var_name, np.float64, ('time','model_flux_level',), fill_value='-9999')
                        else:
                            print 'Diagnostic is on model levels.'
                            dat = dataset.createVariable(cube1[d].var_name, np.float64, ('time','model_level_number',), fill_value='-9999')
                        # dat[h,:] = cube[d].data[h,:]
                        dat[h,:] = np.nanmean([cube1[d].data[h,:],cube2[d].data[h,:]],0)
                    elif np.size(cube1[d].shape) == 3:
                        print 'Diagnostic is 3D, so writing hour = ' + str(h)
                        dat = dataset.createVariable(cube1[d].var_name, np.float64, ('frequency','time','model_level_number',), fill_value='-9999')
                        # dat[:,h,:] = cube1[d].data[:,h,:]
                        dat[:,h,:] = np.nanmean([cube[d].data[:,h,:],cube2[d].data[:,h,:]],1)
                else:
                    if np.size(cube[d].shape) == 0:
                        print 'Diagnostic is a scalar field, so averaging between hour = ' + str(h) + ' and h = ' + str(h+1)
                        # dat[:] = cube[d].data
                        dat[:] = np.nanmean([cube1[d].data,cube2[d].data])
                    elif np.size(cube[d].shape) == 1:
                        print 'Diagnostic is 1D, so averaging between hour = ' + str(h) + ' and h = ' + str(h+1)
                        # dat[h] = cube[d].data[h]
                        dat[h] = np.nanmean([cube1[d].data[h],cube2[d].data[h]])
                    elif np.size(cube[d].shape) == 2:
                        print 'Diagnostic is 2D, so averaging between hour = ' + str(h) + ' and h = ' + str(h+1)
                        # dat[h,:] = cube[d].data[h,:]
                        dat[h,:] = np.nanmean([cube1[d].data[h,:],cube2[d].data[h,:]],0)
                    elif np.size(cube[d].shape) == 3:
                        print 'Diagnostic is 3D, so averaging between hour = ' + str(h) + ' and h = ' + str(h+1)
                        # dat[:,h,:] = cube[d].data[:,h,:]
                        dat[:,h,:] = np.nanmean([cube[d].data[:,h,:],cube2[d].data[:,h,:]],1)
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            dat.units = str(cube[d].units)
            if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
            if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)

    ##################################
    # Loop over hours to load in appropriate file
    ##################################

    # for h in range(0,3):
    #     cube = iris.load('DATA/' + date + '_moccha_ecmwf_' + str(int(data['ship_ind'][h])).zfill(3) + '.nc')
    #
    #     ###################################
    #     ## Write out diagnostics
    #     ###################################
    #     for d in range(0,len(cube)):
    #         print ''
    #         print 'Writing ' + cube[d].var_name
    #         if h == 0:
    #             if np.size(cube[d].shape) == 0:
    #                 print 'Diagnostic is a scalar field, so writing hour = ' + str(h)
    #                 dat = dataset.createVariable(cube[d].var_name, np.float64, fill_value='-9999')
    #                 dat[:] = cube[d].data
    #             elif np.size(cube[d].shape) == 1:
    #                 print 'Diagnostic is 1D, so writing hour = ' + str(h)
    #                 dat = dataset.createVariable(cube[d].var_name, np.float64, ('time',), fill_value='-9999')
    #                 dat[h] = cube[d].data[h]
    #             elif np.size(cube[d].shape) == 2:
    #                 print 'Diagnostic is 2D, so writing hour = ' + str(h)
    #                 if cube[d].var_name in fluxes:
    #                     print 'Diagnostic is on flux levels.'
    #                     dat = dataset.createVariable(cube[d].var_name, np.float64, ('time','model_flux_level',), fill_value='-9999')
    #                 else:
    #                     print 'Diagnostic is on model levels.'
    #                     dat = dataset.createVariable(cube[d].var_name, np.float64, ('time','model_level_number',), fill_value='-9999')
    #                 dat[h,:] = cube[d].data[h,:]
    #             elif np.size(cube[d].shape) == 3:
    #                 print 'Diagnostic is 3D, so writing hour = ' + str(h)
    #                 dat = dataset.createVariable(cube[d].var_name, np.float64, ('frequency','time','model_level_number',), fill_value='-9999')
    #                 dat[:,h,:] = cube[d].data[:,h,:]
    #         else:
    #             if np.size(cube[d].shape) == 0:
    #                 print 'Diagnostic is a scalar field, so writing hour = ' + str(h)
    #                 dat[:] = cube[d].data
    #             elif np.size(cube[d].shape) == 1:
    #                 print 'Diagnostic is 1D, so writing hour = ' + str(h)
    #                 dat[h] = cube[d].data[h]
    #             elif np.size(cube[d].shape) == 2:
    #                 print 'Diagnostic is 2D, so writing hour = ' + str(h)
    #                 dat = dataset.variable[cube[d].var_name][:]
    #                 dat[h,:] = cube[d].data[h,:]
    #             elif np.size(cube[d].shape) == 3:
    #                 print 'Diagnostic is 3D, so writing hour = ' + str(h)
    #                 dat[:,h,:] = cube[d].data[:,h,:]
    #         dat.scale_factor = float(1)
    #         dat.add_offset = float(0)
    #         dat.units = str(cube[d].units)
    #         # dat.attributes = cube[d].attributes
    #         # dat.STASH = str(cube[d].attributes['STASH'])
    #         if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
    #         if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
    #         # dat[:,:] = cube[d].data

    ##################################
    # Write out file
    ##################################
    dataset.close()

    return dataset

def appendMetaNetCDF(outfile, date):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print '******'
    print ''
    print 'Appending metadata to ' + outfile
    print ''

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(outfile, 'a', format ='NETCDF4_CLASSIC')
    # infile = net.Dataset("2015%s%s-160000_0.nc" % (month,day), "a")
    print ''
    print dataset.file_format
    print ''

    ###################################
    ## Global Attributes
    ###################################
    dataset.conventions = 'CF-1.0'
    dataset.title = 'ECMWF Model single-site output during MOCCHA'
    dataset.location = 'MOCCHA'
    # dataset.description = 'Hourly data taken from grid box closest to ship location. Where the ship covers more than one grid box within an hour period, data are averaged from all grid boxes crossed.'
    dataset.description = 'Hourly data combined from n=38 files into 1 daily file.'
    dataset.history = 'ke 5.6.2019 14.09.20 +0300 - NetCDF generated from original data by Ewan O''Connor <ewan.oconnor@fmi.fi> using cnmodel2nc on cloudnet.fmi.fi. Combined from n=38 lat/lon files at ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (Iris and netCDF4).'
    dataset.source = 'ECMWF Integrated Forecast System (IFS)'
    # dataset.references = ''
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    # dataset.comment = ''
    dataset.institution = 'European Centre for Medium-Range Weather Forecasting.'
    # dataset.initialization_time = outfile[0:4] + '-' + outfile[4:6] + '-' + outfile[6:8]) + ' 00:00:00 UTC.'
    dataset.initialization_time = date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' 00:00:00 +00:00'

    dataset.close()

    return dataset

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### JASMIN
    ### LAPTOP
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ECMWF/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '/home/gillian/MOCCHA/ECMWF/DATA/'
        ship_filename = '/home/gillian/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/data/ecmwf_ewan/moccha/ecmwf-all/2018/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    print '******'
    print ''
    print 'Identifying .nc file: '
    print ''

    # # -------------------------------------------------------------
    # # Load data
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin data read in at ' + time.strftime("%c")
    print ' '

    ### -------------------------------------------------------------------------
    ### define input filenames
    ### -------------------------------------------------------------------------
    date = '20180813'
    outfile = date + '_oden_ecmwf.nc'
    print 'Outfile will be: ' + outfile
    base_name = date + '_moccha_ecmwf_'
    names = [None] * 38         ## 'empty' list of 38 elements. can assign index without list.append
    filenames = [None] * 38
    for i in range(0,38):
        id = i+1
        str_i = "%03d" % id
        names[i] = base_name + str_i + '.nc'
        filenames[i] = root_dir + names[i]

    print filenames[0] + ' ... ' + filenames[-1]
    print ''

    # -------------------------------------------------------------
    # Find ECMWF grid
    # -------------------------------------------------------------
    # lats = np.zeros([38])
    # lons = np.zeros([38])
    # tim = np.zeros([24])
    # for i in range(0,38):
    #     lats[i], lons[i], tim = pullLatLon(filenames[i])
    #
    # print 'Lats = ' + str(lats)
    # print 'Lons = ' + str(lons)

    # -------------------------------------------------------------
    # Extract each position file with Iris and write to combined netCDF
    # -------------------------------------------------------------
    data = readDaily(filenames, date)

    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    # map = plot_basemap(ship_data, lats, lons, tim)

    # -------------------------------------------------------------
    # Pull daily gridded ship track from netCDFs
    # -------------------------------------------------------------
    data = pullTrack(ship_data, data, date, outfile)

    ### temporary data save for development/debugging
    # np.save('working_data', data)
    ### load with data = np.load('working_data.npy').item())

    # -------------------------------------------------------------
    # Plot data (cartopy map)
    # -------------------------------------------------------------
    data = plot_cartmap(ship_data, data, date)

    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

    #### DIAGNOSTICS TO CHOOSE FROM:

#     dimensions(sizes): time(25), level(137), flux_level(138), frequency(2)
#     variables(dimensions): float32 latitude(), float32 longitude(),
#           float32 horizontal_resolution(), float32 time(time), float32 forecast_time(time),
#           int16 level(level), int16 flux_level(flux_level), float32 pressure(time,level),
#           float32 uwind(time,level), float32 vwind(time,level), float32 omega(time,level),
#           float32 temperature(time,level), float32 q(time,level), float32 rh(time,level),
#           float32 ql(time,level), float32 qi(time,level), float32 cloud_fraction(time,level),
#           float32 flx_net_sw(time,flux_level), float32 flx_net_lw(time,flux_level),
#           float32 flx_down_sens_heat(time,flux_level), float32 flx_turb_moist(time,flux_level),
#           float32 flx_ls_rain(time,flux_level), float32 flx_ls_snow(time,flux_level),
#           float32 flx_conv_rain(time,flux_level), float32 flx_conv_snow(time,flux_level),
#           float32 flx_turb_mom_u(time,flux_level), float32 flx_turb_mom_v(time,flux_level),
#           float32 sfc_pressure(time), float32 sfc_net_sw(time), float32 sfc_net_lw(time),
#           float32 sfc_down_sw(time), float32 sfc_down_lw(time), float32 sfc_cs_down_sw(time),
#           float32 sfc_cs_down_lw(time), float32 sfc_down_lat_heat_flx(time),
#           float32 sfc_down_sens_heat_flx(time), float32 sfc_ls_rain(time),
#           float32 sfc_conv_rain(time), float32 sfc_ls_snow(time), float32 sfc_conv_snow(time),
#           float32 sfc_ls_precip_fraction(time), float32 sfc_cloud_fraction(time),
#           float32 sfc_bl_height(time), float32 sfc_albedo(time), float32 sfc_temp_2m(time),
#           float32 sfc_q_2m(time), float32 sfc_rough_mom(time), float32 sfc_rough_heat(time),
#           float32 sfc_skin_temp(time), float32 sfc_wind_u_10m(time), float32 sfc_wind_v_10m(time),
#           float32 sfc_geopotential(time), float32 height(time,level), float32 sfc_height_amsl(time),
#           float32 flx_height(time,flux_level), float32 wwind(time,level), float32 frequency(frequency),
#           float32 gas_atten(frequency,time,level), float32 specific_gas_atten(frequency,time,level),
#           float32 specific_saturated_gas_atten(frequency,time,level),
#           float32 specific_dry_gas_atten(frequency,time,level), float32 K2(frequency,time,level),
#           float32 specific_liquid_atten(frequency,time,level)

if __name__ == '__main__':

    main()
