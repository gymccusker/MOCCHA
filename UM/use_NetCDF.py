###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE
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

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
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
    # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=13,data.values[:,1]==8),data.values[:,3]>=0))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ''
    print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')'
    print 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6]))
    print 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7]))
    print ''

    return inIce_index

def writeoutGrid(tim, lat, lon, date):

    import pandas as pd

    # ******
    # write to csv file
    # ******

    print '******'
    print 'Writing daily grid to file:'
    print ''
    dat = np.zeros([len(tim), 3])
    dat[:,0] = tim
    dat[:,1] = lon
    dat[:,2] = lat
    df = pd.DataFrame(dat)
    filename = 'AUX_DATA/' + date + '_ShipTrack_GRIDDED.csv'
    df.to_csv(filename,  sep = " ")
    print '... finished!'
    print ''
    print '******'

def gridShipTrack(cube, xoffset, yoffset):

    import iris.plot as iplt
    import pandas as pd
    # cube.dim_coords[1].coord_system

    ###---------------------------------
    ### 12th August 2018
    ###---------------------------------

    date = '20180812'

    ### HOURS 0 TO 19 IN REGION WHERE BCS CAUSE PROBLEMS

    # ### box pick 0-1h
    lon0 = np.array([253,253,252,252,252,251])
    lat0 = np.array([489,490,490,491,492,492])
    tim0 = np.zeros([np.size(lon0)])
    # for i in range(0,np.size(lon0)-1):
    #     iplt.scatter(cube.dim_coords[2][lon0[i] + xoffset], cube.dim_coords[1][lat0[i] + yoffset],color='yellow')
    #
    # ### box pick 1-2h
    lon1 = np.array([251,251])
    lat1 = np.array([492,493])
    tim1 = np.zeros([np.size(lon1)])
    tim1[:] = 1.0
    tim_128 = np.append(tim0, tim1)
    lat_128 = np.append(lat0, lat1)
    lon_128 = np.append(lon0, lon1)
    # for i in range(0,np.size(lon1)-1):
    #     iplt.scatter(cube.dim_coords[2][lon1[i] + xoffset], cube.dim_coords[1][lat1[i] + yoffset],color='blue')
    #
    # ### box pick 2-3h
    lon2 = np.array([251,251])
    lat2 = np.array([493,492])
    tim2 = np.zeros([np.size(lon2)])
    tim2[:] = 2.0
    tim_128 = np.append(tim_128, tim2)
    lat_128 = np.append(lat_128, lat2)
    lon_128 = np.append(lon_128, lon2)
    # for i in range(0,np.size(lon2)-1):
    #     iplt.scatter(cube.dim_coords[2][lon2[i] + xoffset], cube.dim_coords[1][lat2[i] + yoffset],color='green')
    #
    # ### box pick 3-18h
    tim3 = np.arange(3,18)
    lon3 = np.zeros([np.size(tim3)])
    lon3[:] = 251
    lat3 = np.zeros([np.size(tim3)])
    lat3[:] = 492
    tim_128 = np.append(tim_128, tim3)
    lat_128 = np.append(lat_128, lat3)
    lon_128 = np.append(lon_128, lon3)
    # for i in range(0,np.size(lon3)-1):
    #     iplt.scatter(cube.dim_coords[2][int(lon3[i] + xoffset)], cube.dim_coords[1][int(lat3[i] + yoffset)],color='black')
    #
    # ### box pick 18-19h
    lon18 = np.array([251,251,251])
    lat18 = np.array([492,491,490])
    tim18 = np.zeros([np.size(lon18)])
    tim18[:] = 18.0
    tim_128 = np.append(tim_128, tim18)
    lat_128 = np.append(lat_128, lat18)
    lon_128 = np.append(lon_128, lon18)
    # for i in range(0,np.size(lon18)-1):
    #     iplt.scatter(cube.dim_coords[2][lon18[i] + xoffset], cube.dim_coords[1][lat18[i] + yoffset],color='red')

    ### box pick 19-20h
    lon19 = np.array([251,251,251,252,252])
    lat19 = np.array([490,489,488,488,487])
    tim19 = np.zeros([np.size(lon19)])
    tim19[:] = 19.0
    tim_128 = np.append(tim_128, tim19)
    lat_128 = np.append(lat_128, lat19)
    lon_128 = np.append(lon_128, lon19)
    for i in range(0,np.size(lon19)-1):
        iplt.scatter(cube.dim_coords[2][lon19[i] + xoffset], cube.dim_coords[1][lat19[i] + yoffset],color='blue')

    ### box pick 20-21h
    lon20 = np.array([252,252,252,252,252])
    lat20 = np.array([487,486,485,484,483])
    tim20 = np.zeros([np.size(lon20)])
    tim20[:] = 20.0
    tim_128 = np.append(tim_128, tim20)
    lat_128 = np.append(lat_128, lat20)
    lon_128 = np.append(lon_128, lon20)
    for i in range(0,np.size(lon20)-1):
        iplt.scatter(cube.dim_coords[2][lon20[i] + xoffset], cube.dim_coords[1][lat20[i] + yoffset],color='green')

    ### box pick 21-22h
    lon21 = np.array([252,251,250,249,248])
    lat21 = np.array([483,483,483,483,483])
    tim21 = np.zeros([np.size(lon21)])
    tim21[:] = 21.0
    tim_128 = np.append(tim_128, tim21)
    lat_128 = np.append(lat_128, lat21)
    lon_128 = np.append(lon_128, lon21)
    for i in range(0,np.size(lon21)-1):
        iplt.scatter(cube.dim_coords[2][lon21[i] + xoffset], cube.dim_coords[1][lat21[i] + yoffset],color='black')

    ### box pick 22-23h
    lon22 = np.array([248,248,248,248,248])
    lat22 = np.array([483,482,481,480,479])
    tim22 = np.zeros([np.size(lon22)])
    tim22[:] = 22.0
    tim_128 = np.append(tim_128, tim22)
    lat_128 = np.append(lat_128, lat22)
    lon_128 = np.append(lon_128, lon22)
    for i in range(0,np.size(lon22)-1):
        iplt.scatter(cube.dim_coords[2][lon22[i] + xoffset], cube.dim_coords[1][lat22[i] + yoffset],color='red')

    ### box pick 23-00h
    lon23 = np.array([248,249,249,250,250])
    lat23 = np.array([479,479,478,478,477])
    tim23 = np.zeros([np.size(lon23)])
    tim23[:] = 23.0
    tim_128 = np.append(tim_128, tim23)
    lat_128 = np.append(lat_128, lat23)
    lon_128 = np.append(lon_128, lon23)
    for i in range(0,np.size(lon23)-1):
        iplt.scatter(cube.dim_coords[2][lon23[i] + xoffset], cube.dim_coords[1][lat23[i] + yoffset],color='blue')

    # ******
    # write out index arrays
    # ******

    out = writeoutGrid(tim_128, lat_128, lon_128, date)

    ###---------------------------------
    ### 13th August 2018
    ###---------------------------------

    date = '20180813'

    ### box pick 0-1h
    lon0 = np.array([250,251,252,252])
    lat0 = np.array([477,477,477,476])
    for i in range(0,np.size(lon0)-1):
        iplt.scatter(cube.dim_coords[2][lon0[i] + xoffset], cube.dim_coords[1][lat0[i] + yoffset],color='yellow')

    ### box pick 1-2h    (13th aug)
    lon1 = np.array([252,253])
    lat1 = np.array([476,476])
    for i in range(0,np.size(lon1)-1):
        iplt.scatter(cube.dim_coords[2][lon1[i] + xoffset], cube.dim_coords[1][lat1[i] + yoffset],color='black')

    ### box pick 2-3h    (13th aug)
    lon2 = np.array([253,253,254])
    lat2 = np.array([476,475,475])
    for i in range(0,np.size(lon2)-1):
        iplt.scatter(cube.dim_coords[2][lon2[i] + xoffset], cube.dim_coords[1][lat2[i] + yoffset],color='red')

    ### box pick 3-4h    (13th aug)
    lon3 = np.array([254,254,255,255])
    lat3 = np.array([475,474,474,473])
    for i in range(0,np.size(lon3)-1):
        iplt.scatter(cube.dim_coords[2][lon3[i] + xoffset], cube.dim_coords[1][lat3[i] + yoffset],color='blue')

    ### box pick 4-5h    (13th aug)
    lon4 = np.array([255,256])
    lat4 = np.array([473,473])
    for i in range(0,np.size(lon4)-1):
        iplt.scatter(cube.dim_coords[2][lon4[i] + xoffset], cube.dim_coords[1][lat4[i] + yoffset],color='green')

    ### box pick 5-6h
    lon5 = np.array([256,256,255,255])
    lat5 = np.array([473,472,472,471])
    for i in range(0,np.size(lon5)-1):
        iplt.scatter(cube.dim_coords[2][lon5[i] + xoffset], cube.dim_coords[1][lat5[i] + yoffset],color='black')

    ### box pick 6-17h
    lon6 = np.array([255])
    lat6 = np.array([471])
    for i in range(0,np.size(lon6)-1):
        iplt.scatter(cube.dim_coords[2][lon6[i] + xoffset], cube.dim_coords[1][lat6[i] + yoffset],color='red')

    ### box pick 17-18h
    lon17 = np.array([255,255])
    lat17 = np.array([471,472])
    for i in range(0,np.size(lon17)-1):
        iplt.scatter(cube.dim_coords[2][lon17[i] + xoffset], cube.dim_coords[1][lat17[i] + yoffset],color='blue')

    ### box pick 18-23h
    lon18 = np.array([255])
    lat18 = np.array([472])
    for i in range(0,np.size(lon18)-1):
        iplt.scatter(cube.dim_coords[2][lon18[i] + xoffset], cube.dim_coords[1][lat18[i] + yoffset],color='green')

    ### box pick 23-00h
    lon23 = np.array([255,254])
    lat23 = np.array([472,472])
    for i in range(0,np.size(lon23)-1):
        iplt.scatter(cube.dim_coords[2][lon23[i] + xoffset], cube.dim_coords[1][lat23[i] + yoffset],color='black')

    # print ''
    # print 'tim_128 = ' + str(tim_128.shape)
    # print 'lat_128 = ' + str(lat_128.shape)
    # print 'lon_128 = ' + str(lon_128.shape)
    # print ''

    # out = writeoutGrid(tim_138, lat_138, lon_138, date)


def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==13,data.values[:,1]==8),data.values[:,3]>=23))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]==1))
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

def findLatLon(ship_data, cube, hour):

    print ''
    print 'Finding lat/lon of ship track'
    print '...'

    #################################################################
    ## find ship track coordinates
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)

    # -------------------------------------------------------------
    # Define unrotated coordinate model grid
    # -------------------------------------------------------------
    #### the following uses iris to unrotate the coordinate grid.
    ####    this only works with square domains (i.e. paXXX files)
    ####    only needs to be performed once -- saved grid as .csv file
    lon, lat = unrotateGrid(cube)

    print 'findLatLon testing:'
    print 'Ship (lon,lat): ' + str(ship_data.values[drift_index,7][0]) + ', ' + str(ship_data.values[drift_index,6][0])
    print 'Model unrotated [max, median], (lat,lon): ' + str(np.max(lat)) + ', ' + str(np.median(lon))


    ship_index = np.where(np.logical_and(np.greater_equal(lat[:],ship_data.values[drift_index,7][0]), np.less_equal(lat[:],ship_data.values[drift_index,7][1])))
    print 'Ship index test'
    print ship_index
    # print lat[ship_index[0]


    print 'test complete!'

### Plot tracks as line plot
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
# plt.plot(ship_data.values[drift_index,6], ship_data.values[drift_index,7],
#          color = 'red', linewidth = 4,
#          transform = ccrs.PlateCarree(), label = 'Drift',
#          )

    return lat, lon

def plot_cartmap(ship_data, cube, hour): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

    # print 'xoffset = ', xoffset
    # print 'yoffset = ', yoffset

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
    plt.figure(figsize=(5,4))
    # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    ax.set_extent([0, 60, 89.5, 90], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([0, 60, 87.75, 90], crs=ccrs.PlateCarree())     ### SWATH
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
    # if np.size(cube.shape) == 4:
    #     iplt.pcolormesh(cube[hour,0,:,:])
    # elif np.size(cube.shape) == 3:
    #     iplt.pcolormesh(cube[hour,:,:])
    #     # iplt.pcolormesh(cube[hour,471:495,240:264])
    # elif np.size(cube.shape) == 2:
    #     iplt.pcolormesh(cube[:,:])
    # plt.title(cube.standard_name + ', ' + str(cube.units))
    # plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    qplt.outline(cube[hour,:,:])

    gridship = gridShipTrack(cube, xoffset, yoffset)

            #### MID POINT: (433, 258)

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    # drift_index = iceDrift(ship_data)
    # inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data)

    ### Plot tracks as line plot
    plt.plot(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7],
             color = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(), label = 'Ship track',
             )
    plt.plot(ship_data.values[trackShip_index[0],6], ship_data.values[trackShip_index[0],7],
             'k^', markerfacecolor = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )
    plt.plot(ship_data.values[trackShip_index[-1],6], ship_data.values[trackShip_index[-1],7],
             'kv', markerfacecolor = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )

    ### Plot tracks as line plot
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
    # plt.plot(ship_data.values[drift_index,6], ship_data.values[drift_index,7],
    #          color = 'red', linewidth = 4,
    #          transform = ccrs.PlateCarree(), label = 'Drift',
    #          )

    #### test plotting of unrotated grid
    # lon, lat = unrotateGrid(cube)

    # plt.plot(np.nanmin(lon),np.nanmin(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmin(lon),np.nanmax(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmax(lon),np.nanmin(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmax(lon),np.nanmax(lat),
    #         color='black',transform = ccrs.PlateCarree())

    plt.legend()

    print '******'
    print ''
    print 'Finished plotting cartopy map! :)'
    print ''

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

def unrotateGrid(cube):
    ##
    ##      ROTATE GRID BACK TO STANDARD
    ##

    import iris.analysis.cartography as ircrt
    import pandas as pd

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 3.375 # UM SUITE: 37.5000
    pole_lon = 210.0 # UM SUITE: 177.5000

    rot_lat = cube.coord('grid_latitude').points
    rot_lon = cube.coord('grid_longitude').points

    rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, pole_lon, pole_lat)

    # Print to check conversion
    print '******'
    print 'Test of unrotated coordinate grid: '
    print 'Rotated lon coord = ', rot_lon[250]
    print 'Rotated lat coord = ', rot_lat[250]
    print 'Lon = ', lon[250]
    print 'Lat = ', lat[250]
    print ' '

    # ******
    # write to csv file
    # ******

    # print '******'
    # print 'Writing unrotated coordinate grid to file:'
    # print ''
    # lonp = pd.DataFrame(lon)
    # latp = pd.DataFrame(lat)
    # dat = {'Latitude': lat, 'Longitude': lon}
    # df = pd.DataFrame(dat,columns=['Latitude','Longitude'])
    # df.to_csv('POSITION_UNROTATED.csv',  sep = " ")
    # print '... finished!'
    # print ''
    # print '******'

    return lon, lat

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'JASMIN'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
        position_filename = 'POSITION_UNROTATED.csv'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir = '3_12AUG_SWATH_2FCSTS/'

    ## 1_20160401_61DIAG_TEST/
    ## 2_20180801_61DIAGS_TEST/2_30_86.625/
    ## 3_12AUG_SWATH_2FCSTS/

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

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    filename1 = root_dir + out_dir + 'umnsaa_pa012_r0.nc'
    print filename1
    print ''

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin cube read in at ' + time.strftime("%c")
    print ' '
    var = 'surface_net_downward_shortwave_flux'
    cube = iris.load_cube(filename1, var)
    # data = Dataset(filename1,'r')

    print cube
    print ''

    ### paXXX
    # <iris 'Cube' of air_pressure_at_sea_level / (Pa) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of air_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of dew_point_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of relative_humidity / (%) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of specific_humidity / (1) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_air_pressure / (Pa) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_net_downward_longwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_net_downward_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of toa_incoming_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of toa_outgoing_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 8; grid_latitude: 501; grid_longitude: 500)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 8; grid_latitude: 501; grid_longitude: 500)>]

    #### 12 AUG ONLY - NO FULL NEST DIAGNOSTICS
    # <iris 'Cube' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_longwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_incoming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>]


    ### pbXXX
    # [<iris 'Cube' of cloud_area_fraction_assuming_maximum_random_overlap / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of cloud_area_fraction_assuming_random_overlap / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of m01s03i241 / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of Turbulent mixing height after boundary layer / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of wet_bulb_freezing_level_altitude / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of air_pressure_at_sea_level / (Pa) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of air_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of atmosphere_boundary_layer_thickness / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of dew_point_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of high_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of low_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of medium_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of relative_humidity / (%) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of specific_humidity / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of stratiform_rainfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of stratiform_snowfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_air_pressure / (Pa) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_upward_latent_heat_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_upward_sensible_heat_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>]


    ### pcXXX
    # <iris 'Cube' of cloud_volume_fraction_in_atmosphere_layer / (1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of m01s04i118 / (1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of air_pressure / (Pa) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of air_temperature / (K) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of specific_humidity / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of upward_air_velocity / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>]

    # 0: total_radar_reflectivity / (unknown) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 1: air_pressure / (Pa)                 (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 2: air_temperature / (K)               (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 3: eastward_wind / (m s-1)             (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 4: large_scale_cloud_area_fraction / (1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 5: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 6: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 7: northward_wind / (m s-1)            (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 8: specific_humidity / (kg kg-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 9: upward_air_velocity / (m s-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)

    # FORECAST_PERIOD = cube1.aux_coords[1][:]

    # -------------------------------------------------------------
    # Define unrotated coordinate grid
    # -------------------------------------------------------------
    #### the following uses iris to unrotate the coordinate grid.
    ####    this only works with square domains (i.e. paXXX files)
    ####    only needs to be performed once -- saved grid as .csv file
    # lon, lat = unrotateGrid(cube)

    # hour = 0
    # test = findLatLon(ship_data, cube, hour)

    ############## DOESN'T WORK
    #### read in saved unrotated coordinate grid
    # position_data = readfile(position_filename)
    # lon = position_data.values[:,2]     ### unrotated longitude
    # lat = position_data.values[:,1]     ### unrotated latitude
    # lon, lat = np.meshgrid(lon,lat)     ### mesh for use with model diags

    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    ### select hour to plot
    hour = 12
    map = plot_cartmap(ship_data, cube, hour)#, lon, lat)



    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''


if __name__ == '__main__':

    main()
