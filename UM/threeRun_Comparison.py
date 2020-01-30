###
###
### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
###
###

# from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import seaborn as sns

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY
from readMAT import readMatlabStruct

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
    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    Sep_inIce = np.where(np.logical_and(data.values[:,2]<20,data.values[:,1]==9))
    inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    # Aug_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8),data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=13,data.values[:,1]==8),data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    # inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

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

def trackShip(data, date):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==25,data.values[:,1]==8),data.values[:,3]==1))
    # trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==int(date[-2:]),data.values[:,1]==int(date[-4:-2])),data.values[:,3]>=0))
    # trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==(int(date[-2:]) + 1),data.values[:,1]==int(date[-4:-2])),data.values[:,3]==1))
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

def plot_cartmap(ship_data, cube, hour, date): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.path as mpath
        # from matplotlib.patches import Polygon

    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print 'What grid are we looking at?'
    if len(cube[0].dim_coords[-1].points) == 25:
    # if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif len(cube[0].dim_coords[-1].points) == 56:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 95:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    else:
    # elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

    print 'Because cube shape = ', str(len(cube[0].dim_coords[-1].points))
    print 'xoffset = ', xoffset
    print 'yoffset = ', yoffset

    ###################################
    ## CHOOSE DIAGNOSTIC
    ###################################
    diag = 2
    print ''
    print 'Diag is: ', cube[diag].long_name
    ### pcXXX
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
    plt.figure(figsize=(6,8))#, dpi=300)
    ax = plt.axes(projection=ccrs.Orthographic(30, 70))    # NP Stereo
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    ax.set_extent([0, 60, 86.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE
    # ax.set_extent([-180, 180, 70, 90], crs=ccrs.PlateCarree())    ### V LARGE
    # ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())    ### POSTER
    # ax.set_global()

    ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

    #################################################################
    ## add geographic features/guides for reference
    #################################################################
    ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=0)
    ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=0, edgecolor='black')
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.set_global()
    ax.gridlines()

    # # Compute a circle in axes coordinates, which we can use as a boundary
    # # for the map. We can pan/zoom as much as we like - the boundary will be
    # # permanently circular.
    # theta = np.linspace(0, 2*np.pi, 100)
    # center, radius = [0.5, 0.5], 0.5
    # verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    # circle = mpath.Path(verts * radius + center)
    #
    # ax.set_boundary(circle, transform=ax.transAxes)

    #################################################################
    ## plot UM data
    ################################################################
    # if np.size(cube[diag].data.shape) == 4:
    #     iplt.pcolormesh(cube[diag][hour,0,:,:])
    # elif np.size(cube[diag].data.shape) == 3:
    #     iplt.pcolormesh(cube[diag][hour,:,:])
    #     # iplt.pcolormesh(cube[hour,471:495,240:264])
    # elif np.size(cube[diag].data.shape) == 2:
    iplt.pcolormesh(cube[diag][290:370,150:230])
    # # plt.title(cube[diag].standard_name + ', ' + str(cube[diag].units))
    # plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[diag][hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[diag][hour,386:495,211:305])          ### misc
    # qplt.outline(cube[diag][290:370,150:230])



    # gridship = gridShipTrack(cube[diag], xoffset, yoffset)

            #### MID POINT: (433, 258)

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

    ### Plot full track as line plot
    # plt.plot(ship_data.values[:,6], ship_data.values[:,7], '--',
    #          color = 'pink', linewidth = 2,
    #          transform = ccrs.PlateCarree(), label = 'Whole',
    #          )
    # plt.plot(ship_data.values[inIce_index,6], ship_data.values[inIce_index,7],
    #          color = 'palevioletred', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'In Ice',
    #          )
    # plt.plot(ship_data.values[inIce_index[0],6], ship_data.values[inIce_index[0],7],
    #          'k^', markerfacecolor = 'palevioletred', linewidth = 3,
    #          transform = ccrs.PlateCarree(),
    #          )
    # plt.plot(ship_data.values[inIce_index[-1],6], ship_data.values[inIce_index[-1],7],
    #          'kv', markerfacecolor = 'palevioletred', linewidth = 3,
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

    # plt.savefig('FIGS/HighArctic_vPOSTER.svg', dpi=100)
    plt.show()

def plot_contour_TS(cube, filename): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###################################
    ## CHOOSE DIAGNOSTIC
    ###################################
    diag = 2
    print ''
    print 'Diag is: '
    print cube[diag]
    ### pcXXX
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

    ###################################
    ## DEFINE DIMENSIONS COORDS DEPENDING ON DIAG
    ###################################

    time = cube[diag].dim_coords[0].points
    height = cube[diag].dim_coords[1].points

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting contour timeseries:'
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
    plt.figure(figsize=(8,6))
    ax = plt.gca()

    # plt.plot(cube[diag].dim_coords[0].points,cube[diag][:,0].data)        # line plot
    # plt.contourf(cube[0].data)
    # plt.plot(cube[2][0,:].data,height);plt.show()
    #################################################################
    ## plot contour timeseries
    ################################################################
    plt.contourf(time,height,np.transpose(cube[diag].data))
    # plt.pcolormesh(time,height,np.transpose(cube[2].data))
    plt.title(cube[diag].standard_name + ', ' + str(cube[diag].units))
    plt.colorbar()
    ax.set_ylim([0, 3000])

    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

def plot_multicontour_TS(cube, filename, out_dir): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting contour timeseries:'
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
    plt.figure(figsize=(12,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    l = -1
    for i in range(0,len(cube)):
        ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
        if np.sum(cube[i].data.shape) > 24:

                ###################################
                ## CHOOSE DIAGNOSTIC
                ###################################
                diag = i
                print ''
                print 'Diag is: '
                print cube[diag]

                ### define empty array for cube data
                data = []

                ### pcXXX
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

                ###################################
                ## DEFINE DIMENSIONS COORDS DEPENDING ON DIAG
                ###################################

                time = cube[diag].dim_coords[0].points
                height = cube[diag].dim_coords[1].points

                ### if mass mixing ratio, *1e3 to change to g/kg
                # if cube[diag].var_name[0] == 'q':
                #     data = np.transpose(np.squeeze(cube[diag].data[:,ind]*1e3))
                # elif cube[diag].var_name == 'pressure':
                #     data = np.transpose(np.squeeze(cube[diag].data[:,ind]/1e2))
                # else:
                #     data = np.transpose(np.squeeze(cube[diag].data[:,ind]))

                #################################################################
                ## data corrections
                #################################################################
                ### set height limit to consider
                ind = np.where(height<5000)

                if cube[diag].var_name == 'temperature':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'
                elif cube[diag].var_name == 'qice':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]*1e3))
                    title = cube[diag].var_name + ' [g/kg]'
                elif cube[diag].var_name == 'qliq':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]*1e3))
                    title = cube[diag].var_name + ' [g/kg]'
                elif cube[diag].var_name == 'q':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]*1e3))
                    title = cube[diag].var_name + ' [g/kg]'
                elif cube[diag].var_name == 'pressure':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]/1e2))
                    title = cube[diag].var_name + ' [hPa]'
                elif cube[diag].var_name == 'uwind':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'
                elif cube[diag].var_name == 'wwind':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'
                elif cube[diag].var_name == 'radr_refl':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'
                elif cube[diag].var_name == 'cloud_fraction':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'
                elif cube[diag].var_name == 'vwind':
                    data = np.transpose(np.squeeze(cube[diag].data[:,ind]))
                    title = cube[diag].var_name + ' [' + str(cube[diag].units) + ']'

                #################################################################
                ## create figure and axes instances
                #################################################################
                if len(data) > 0:
                    l = l + 1 ## increment index for positive data association

                    print 'l = ' + str(l)
                    print title

                    plt.subplot(5,2,l+1)
                    ax = plt.gca()

                    #################################################################
                    ## plot timeseries
                    #################################################################
                    # plt.contourf(time,height,np.transpose(cube[diag].data))
                    if cube[diag].var_name == 'temperature':
                        plt.pcolormesh(time, height[ind], data, vmin = 250, vmax = np.nanmax(data))
                    elif cube[diag].var_name == 'uwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -20, vmax = 20)
                    elif cube[diag].var_name == 'vwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -20, vmax = 20)
                    elif cube[diag].var_name == 'wwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -0.1, vmax = 0.1)
                    else:
                        plt.pcolormesh(time, height[ind], data, vmin = np.nanmin(data), vmax = np.nanmax(data))

                    #################################################################
                    ## set plot properties
                    #################################################################
                    ### colormaps:
                    if cube[diag].var_name == 'wwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif cube[diag].var_name == 'uwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif cube[diag].var_name == 'vwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif cube[diag].var_name[0] == 'q':
                        plt.set_cmap(mpl_cm.Blues)
                    else:
                        plt.set_cmap(mpl_cm.viridis)

                    ### title and axes properties
                    plt.title(title)
                    plt.colorbar()
                    ax.set_ylim([0, 5000])

                    ### global plot properties
                    plt.subplot(5,2,9)
                    plt.xlabel('Time [UTC]')
                    plt.ylabel('Z [m]')
                    plt.subplot(5,2,10)
                    plt.xlabel('Time [UTC]')
                    plt.subplot(5,2,1)
                    plt.ylabel('Z [m]')
                    plt.subplot(5,2,3)
                    plt.ylabel('Z [m]')
                    plt.subplot(5,2,5)
                    plt.ylabel('Z [m]')
                    plt.subplot(5,2,7)
                    plt.ylabel('Z [m]')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if out_dir[:18] == '5_u-bl616_RA2M_CAS':
        fileout = 'FIGS/' + out_dir[:21] + filename[-22:-3] + '.png'
    elif out_dir[:18] == '4_u-bg610_RA2M_CON':
        fileout = 'FIGS/' + out_dir[:19] + filename[-22:-3] + '.png'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_multicontour_multidate_TS(timem, data, cube, month_flag, missing_files, out_dir): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined contour timeseries:'
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
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    print data.keys()

    l = -1

    # for i in range(0,len(data.keys())):
    for diag in range(0,len(data.keys())):
        ###################################
        ## CHOOSE DIAGNOSTIC
        ###################################
        # diag = i
        print ''
        print 'Diag is: '
        print data.keys()[diag]

        ### pcXXX - CUBE
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

        ###################################
        ## DEFINE DIMENSIONS COORDS DEPENDING ON DIAG
        ###################################

        ### the following works for now, but could do with finding an easier
        ###         way to index cube by string
        for j in range(0,len(cube)):
            if cube[j].var_name == 'qliq': height = cube[j].dim_coords[1].points

        dat = []
        print 'dat = ' + str(len(dat))

        ### set diag-specific titles
        if str(data.keys()[diag]) == "radr_refl":
            title = str(data.keys()[diag]) + ' [dBz]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "pressure":
            title = str(data.keys()[diag]) + ' [hPa]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data/1e2))
        elif str(data.keys()[diag]) == "temperature":
            title = str(data.keys()[diag]) + ' [K]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "uwind":
            title = str(data.keys()[diag]) + ' [m/s]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "wwind":
            title = str(data.keys()[diag]) + ' [m/s]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "vwind":
            title = str(data.keys()[diag]) + ' [m/s]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "cloud_fraction":
            title = str(data.keys()[diag]) + ' []'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qice":
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qliq":
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qnrain":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qnsnow":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qnice":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qprisice":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qnliq":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qngraup":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qrain":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "qgraup":
            casim_flag = 1
            dat = []
        elif str(data.keys()[diag]) == "q":
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))

        print 'dat = ' + str(len(dat))
        print ''

        if len(dat) > 0:
            print 'Plotting since len(dat)>0'
            print ''
            print title

            l = l + 1 ## increment index for positive data association
            print 'l = ' +  str(l)
            #################################################################
            ## create figure and axes instances
            #################################################################
            plt.subplot(5,2,l+1)
            ax = plt.gca()

            #################################################################
            ## plot timeseries
            #################################################################
            # plt.contourf(time,height,np.transpose(cube[diag].data))
            if str(data.keys()[diag]) == 'temperature':
                plt.pcolormesh(timem, height, dat, vmin = 250, vmax = np.nanmax(dat))
            elif str(data.keys()[diag]) == 'pressure':
                plt.pcolormesh(timem, height, dat, vmin = 500, vmax = np.nanmax(dat))
            elif str(data.keys()[diag]) == 'uwind':
                plt.pcolormesh(timem, height, dat, vmin = -20, vmax = 20)
            elif str(data.keys()[diag]) == 'vwind':
                plt.pcolormesh(timem, height, dat, vmin = -20, vmax = 20)
            elif str(data.keys()[diag]) == 'wwind':
                plt.pcolormesh(timem, height, dat, vmin = -0.1, vmax = 0.1)
            elif str(data.keys()[diag]) == 'qice':
                plt.pcolormesh(timem, height, dat, vmin = 0, vmax = 0.05)
            else:
                plt.pcolormesh(timem, height, dat, vmin = np.nanmin(dat), vmax = np.nanmax(dat))
            #################################################################
            ## set plot properties
            #################################################################
            ### colormaps:
            if str(data.keys()[diag]) == 'wwind':
                plt.set_cmap(mpl_cm.RdBu_r)
            elif str(data.keys()[diag]) == 'uwind':
                plt.set_cmap(mpl_cm.RdBu_r)
            elif str(data.keys()[diag]) == 'vwind':
                plt.set_cmap(mpl_cm.RdBu_r)
            elif str(data.keys()[diag][0]) == 'q':
                plt.set_cmap(mpl_cm.Blues)
            else:
                plt.set_cmap(mpl_cm.viridis)

            plt.title(title)
            plt.colorbar()
            ax.set_ylim([0, 5000])
            if month_flag == 8: ax.set_xlim([13.0, 32.0])
            if month_flag == 9: ax.set_xlim([1.0, 15.0])
            if month_flag == -1: ax.set_xlim([225.0, 258.0])

            print ''
            print 'Zero out any data from missing files:'
            print ''
            for mfile in missing_files:
                # mtime = float(mfile[6:8]) + ((cube[0].dim_coords[0].points)/24.0)
                # nans = np.zeros([len(height),len(mtime)])
                # # nans[nans == 0] = np.nan
                # plt.pcolormesh(mtime, height, nans)
                mtime = float(mfile[6:8]) + ((cube[0].dim_coords[0].points)/24.0)
                nans = ax.get_ylim()
                ax.fill_between(mtime, nans[0], nans[-1], facecolor = 'lightgrey', zorder = 3)

    ### global plot properties
    plt.subplot(5,2,9)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.ylabel('Z [m]')
    plt.subplot(5,2,10)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(5,2,1)
    plt.ylabel('Z [m]')
    plt.subplot(5,2,3)
    plt.ylabel('Z [m]')
    plt.subplot(5,2,5)
    plt.ylabel('Z [m]')
    plt.subplot(5,2,7)
    plt.ylabel('Z [m]')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201808_oden_metum.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201808_oden_metum.png'
    if month_flag == 9:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201809_oden_metum.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201809_oden_metum.png'
    if month_flag == -1:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum.png'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined 1d timeseries:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    # datenums_obs['foremast'] = obs['foremast'].variables['time'][:] ### obs['foremast'] data on different timestep
    # time_obs['foremast'] = calcTime_Mat2DOY(datenums_obs['foremast'])

    # datenums_obs['deck7th'] = obs['deck7th'].variables['time'][:] ### 7th deck data on different timestep
    # time_obs['deck7th'] = calcTime_Mat2DOY(datenums_obs['deck7th'])

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)
    # ax = plt.gca()
    # plt.plot(data1['time'], data1['sfc_pressure'].data/1e2, color = 'steelblue', label = label1)
    # plt.plot(data2['time'], data2['sfc_pressure'].data/1e2, color = 'forestgreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'darkorange', label = label3)
    # else:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'darkorange',label = label3)
    # plt.plot(obs['foremast'].variables['doy'][:],obs['foremast'].variables['Psurf'][:], 'k')
    # plt.title('sfc_pressure [hPa]')
    # plt.ylim([980, 1010])
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'steelblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'forestgreen')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    # data1['surface_net_SW_radiation'].data[data1['surface_net_SW_radiation'].data == 0] = np.nan
    # data2['surface_net_SW_radiation'].data[data2['surface_net_SW_radiation'].data == 0] = np.nan
    # if out_dir4 != 'OUT_25H/': data3['surface_net_SW_radiation'].data[data3['surface_net_SW_radiation'].data == 0] = np.nan

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'forestgreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    # plt.subplot(3,3,6)
    # ax = plt.gca()
    # plt.plot(time_um, data['snowfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_snow'].data)
    # plt.title('sfc_snow_amount [kg/m2]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,7)
    # ax = plt.gca()
    # plt.plot(time_um, data['rainfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_rain'].data)
    # plt.title('sfc_rain_amount [kg/m2]')
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station']['mday'],obs['ice_station']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    # if month_flag == -1:
        # if out_dir1[:20] == '5_u-bl661_RA1M_CASIM':
            # if out_dir2[:20] == '6_u-bm410_RA1M_CASIM':
            #     if out_dir4 == 'OUT_25H/':
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] + '_oden_metum_ifs_casim_TSa.svg'
            #     else:
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_casim_TSa.png'
            # elif out_dir2[:9] == '4_u-bg410':
            #     if out_dir4 == 'OUT_25H/':
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_ifs_casim_TSa.svg'
            # else:
            #     fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_oden_metum_casim_TSa.svg'
        # if out_dir2[:20] == '5_u-bl661_RA1M_CASIM':
        #     if out_dir4[:20] == '6_u-bm410_RA1M_CASIM':
        #         fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_' + out_dir4[:9] + '_oden_metum_casim-100_200_TSa.png'
        #     elif out_dir4 == 'OUT_25H/':
        #         fileout = '../FIGS/comparisons/' + out_dir2[:20] + '_oden_metum_ifs_casim-100_TSa.svg'
        # elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
        #     fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] +'_oden_metum_casim-100_TSa.svg'
    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_ra2t_ifs_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined 1d timeseries:'
    print ''

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 18

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(8,9))
    plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.9, left = 0.1,
            hspace = 0.4, wspace = 0.15)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    # time_um1 = data1['time'][:]
    # time_um2 = data2['time'][:]
    # time_um3 = data3['time'][:]

    #################################################################
    ## create figure and axes instances
    #################################################################

    plt.subplot(211)
    ax = plt.gca()
    plt.plot(data1['time'][:], data1['temp_1.5m'].data - 273.15, color = 'steelblue', label = label1)
    plt.plot(data2['time'][:], data2['temp_1.5m'].data - 273.15, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_temp_2m'].data - 273.15, color = 'darkorange', label =  label3)
    else:
        plt.plot(data3['time'][:], data3['temp_1.5m'].data - 273.15, color = 'darkorange')#, label = '2m')
    plt.plot(time_temp,obs['obs_temp'].variables['Tice'][:] - 273.15, color = 'black', label = 'Observations')
    plt.legend()
    plt.title('Temperature [$^{o}C$]')
    plt.ylim([263 - 273,275 - 273])
    # plt.grid('on')
    if month_flag == 8:
        ax.set_xlim([13.0, 31.0])
        plt.xlabel('Day of month [Aug]')
    if month_flag == 9:
        ax.set_xlim([1.0, 15.0])
        plt.xlabel('Day of month [Sep]')
    if month_flag == -1:
        ax.set_xlim([doy[0],doy[-1]])
        # plt.xlabel('Day of year')

    plt.subplot(2,1,2)
    ax = plt.gca()
    plt.plot(data1['time'][:], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'][:], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'][:], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.title('Net SW radiation [W/m2]')
    plt.ylim([0,100])
    # plt.grid('on')
    if month_flag == 8:
        ax.set_xlim([13.0, 31.0])
        plt.xlabel('Day of month [Aug]')
    if month_flag == 9:
        ax.set_xlim([1.0, 15.0])
        plt.xlabel('Day of month [Sep]')
    if month_flag == -1:
        ax.set_xlim([doy[0],doy[-1]])
        plt.xlabel('Day of year')

    # plt.subplot(3,1,3)
    # ax = plt.gca()
    # # data['surface_net_LW_radiation'].data[data['surface_net_LW_radiation'].data == 0] = np.nan
    # plt.plot(data1['time'][:], data1['surface_net_LW_radiation'].data, color = 'steelblue', label = label1)
    # plt.plot(data2['time'][:], data2['surface_net_LW_radiation'].data, color = 'forestgreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'][:], data3['sfc_net_lw'].data, color = 'darkorange')
    # else:
    #     plt.plot(data3['time'][:], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    # plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'Observations')
    # # plt.legend()
    # plt.title('Net LW radiation [W/m2]')
    # # plt.ylim([260,275])
    # # plt.grid('on')
    # if month_flag == 8:
    #     ax.set_xlim([13.0, 31.0])
    #     plt.xlabel('Day of month [Aug]')
    # if month_flag == 9:
    #     ax.set_xlim([1.0, 15.0])
    #     plt.xlabel('Day of month [Sep]')
    # if month_flag == -1:
    #     ax.set_xlim([doy[0],doy[-1]])
    #     plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        if out_dir2[0:20] == '6_u-bm410_RA1M_CASIM':
            fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-200_tempoC_SW.png'
        elif out_dir2[0:20] == '5_u-bl661_RA1M_CASIM':
            if ifs_flag == True:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
            else:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-100_tempoC_SW.png'
        elif out_dir2[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_tempoC_SW.png'
        else:
            fileout = '../FIGS/comparisons/' + out_dir2[:18] + '_oden_metum_tempoC_SW.png'
    # print 'Saving as: ' + fileout
    fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
    plt.savefig(fileout)#, dpi=400)
    plt.show()

def plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined 1d timeseries:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'steelblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'firebrick')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'firebrick')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'firebrick')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station']['mday'],obs['ice_station']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'firebrick')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    fileout = '../FIGS/comparisons/CASIM-NiceTest_C86-F62-M92_DOY243-249.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined 1d timeseries:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'steelblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'purple')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'purple')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'purple')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station']['mday'],obs['ice_station']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'purple')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_ra2t_ifs_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_ERAI_GLM(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined 1d timeseries:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'steelblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'firebrick')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'firebrick')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'firebrick')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station']['mday'],obs['ice_station']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'firebrick')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_erai-glm_ifs_TSa.svg'
    # plt.savefig(fileout, dpi=300)
    plt.show()

def plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting timeseries of turbulent fluxes:'
    print ''

    ##################################################
    ##################################################
    #### 	SET AXES PROPERTIES
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 18
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### set diagnostic naming flags for if IFS being used
    ### -------------------------------
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### -------------------------------
    ### for reference in figures
    ### -------------------------------
    zeros = np.zeros(len(data2['time']))

    ### -------------------------------
    ### change matlab time for obs
    ### -------------------------------
    time_iceStation = calcTime_Mat2DOY(np.squeeze(obs['ice_station']['mday'][:]))

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12))

    ax  = fig.add_axes([0.07,0.55,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1],
        obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    plt.plot(time_iceStation[np.squeeze(obs['ice_station']['taflag'][:]==1)],
        np.squeeze(obs['ice_station']['taflux'][obs['ice_station']['taflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 6, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')
    plt.ylim([-20, 40])
    plt.title('sensible_heat_flux [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.1,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1],
        obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    # index = np.logical_and(obs['ice_station']['lrflux']>=-30, obs['ice_station']['lrflux']<=70)
    plt.plot(time_iceStation[np.squeeze(obs['ice_station']['lrflag'][:]==1)],
        np.squeeze(obs['ice_station']['lrflux'][obs['ice_station']['lrflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 6, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.title('latent_heat_flux [W/m2]')
    plt.ylim([-20, 60])
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)
    ax  = fig.add_axes([0.7,0.55,0.27,0.35])   # left, bottom, width, height
    # zerosC = np.zeros(len(data2['time']))
    yCmax = 0.16
    plt.plot([0,0],[0,yCmax],'--', color='lightgrey')
    ##---
    indextaum = np.logical_and(data1['sensible_heat_flux'].data>=-50, data1['sensible_heat_flux'].data<=50)
    sns.distplot(data1['sensible_heat_flux'][indextaum].data, hist=False, color="steelblue", kde_kws={"shade": True}, label = label1)
    ##---
    indextaifs = np.logical_and(data3['sfc_down_sens_heat_flx'].data * -1.0 >=-50, data3['sfc_down_sens_heat_flx'].data * -1.0 <=50)
    sns.distplot(data3['sfc_down_sens_heat_flx'].data * -1.0, hist=False, color="darkorange", kde_kws={"shade": True}, label = label3)
    ##---
    indextacasim = np.logical_and(data2['sensible_heat_flux'].data>=-50, data2['sensible_heat_flux'].data<=50)
    sns.distplot(data2['sensible_heat_flux'][indextacasim].data, hist=False, color="forestgreen", kde_kws={"shade": True}, label = label2)
    ##---
    fmst_taflux = obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1]
    indextafmst = np.logical_and(fmst_taflux>=-50, fmst_taflux<=50)
    sns.distplot(fmst_taflux[indextafmst], hist=False, color="black", label = 'Foremast')#, kde_kws={"shade": True}, label = 'Foremast')
    ##---
    taflux = np.squeeze(obs['ice_station']['taflux'][obs['ice_station']['taflag'][:] == 1])
    indexta = np.logical_and(taflux>=-50, taflux<=50)
    sns.distplot(taflux[indexta], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3}, label = 'Ice_station')
    plt.title('sensible_heat_flux [W/m2]')
    plt.legend()
    plt.xlim([-20,50])
    plt.ylim([0,yCmax])

    ax  = fig.add_axes([0.7,0.1,0.27,0.35])   # left, bottom, width, height
    yDmax = 0.12
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    ##---
    indexlrum = np.logical_and(data1['latent_heat_flux'].data>=-50, data1['latent_heat_flux'].data<=50)
    sns.distplot(data1['latent_heat_flux'][indexlrum].data, hist=False, color="steelblue", kde_kws={"shade": True})
    ##---
    indexlrifs = np.logical_and(data3['sfc_down_lat_heat_flx'].data * -1.0>=-50, data3['sfc_down_lat_heat_flx'].data * -1.0<=50)
    sns.distplot(data3['sfc_down_lat_heat_flx'][indexlrifs].data * -1.0, hist=False, color="darkorange", kde_kws={"shade": True})
    ##---
    indexlrcasim = np.logical_and(data2['latent_heat_flux'].data>=-50, data2['latent_heat_flux'].data<=50)
    sns.distplot(data2['latent_heat_flux'][indexlrcasim].data, hist=False, color="forestgreen", kde_kws={"shade": True})
    ##---
    fmst_lrflux = obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1]
    indexlrfmst = np.logical_and(fmst_lrflux>=-50, fmst_lrflux<=50)
    sns.distplot(fmst_lrflux[indexlrfmst], hist=False, color="black")#, kde_kws={"shade": True})
    ##---
    lrflux = np.squeeze(obs['ice_station']['lrflux'][obs['ice_station']['lrflag'][:] == 1])
    indexlr = np.logical_and(lrflux>=-50, lrflux<=50)
    sns.distplot(lrflux[indexlr], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3})
    plt.title('latent_heat_flux [W/m2]')
    plt.xlim([-20,50])
    plt.ylim([0,yDmax])

    fileout = '../FIGS/comparisons/SHF_LHF_line+PDFS_oden_foremast+iceStationQC_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_paperRadiation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined timeseries and PDFs of radiation terms:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(16,12))

    ax  = fig.add_axes([0.07,0.7,0.56,0.22])   # left, bottom, width, height
    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.4,0.56,0.22])   # left, bottom, width, height
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    # plt.legend()
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.1,0.56,0.22])   # left, bottom, width, height
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'forestgreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)
    ax  = fig.add_axes([0.7,0.7,0.27,0.22])   # left, bottom, width, height
    yDmax = 0.05
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    sns.distplot(data1['surface_net_SW_radiation'].data + data1['surface_net_LW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(data2['surface_net_SW_radiation'].data + data2['surface_net_LW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW + netSW, hist=False, color="black")
    plt.title('CRF [W/m2]')
    plt.xlim([-50,80])
    plt.ylim([0,yDmax])

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.4,0.27,0.22])   # left, bottom, width, height
    yEmax = 0.06
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(data1['surface_net_SW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True}, label = label1)
    sns.distplot(data3['sfc_net_sw'].data, hist=False, color="darkorange", kde_kws={"shade": True}, label = label3)
    sns.distplot(data2['surface_net_SW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True}, label = label2)
    sns.distplot(netSW, hist=False, color="black", label = 'Ice_station')
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.1,0.27,0.22])   # left, bottom, width, height
    yFmax = 0.12
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    sns.distplot(data1['surface_net_LW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(data3['sfc_net_lw'].data, hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(data2['surface_net_LW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW, hist=False, color="black")
    plt.title('surface_net_LW_radiation [W/m2]')
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    fileout = '../FIGS/comparisons/CRF_netSW_netLW_line+PDFS_oden_iceStation_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_Precipitation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting combined timeseries and PDFs of precipitation terms:'
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
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    data1['rainfall_flux'][data1['rainfall_flux'] < 0] = np.nan
    data2['rainfall_flux'][data2['rainfall_flux'] < 0] = np.nan
    # flx_ls_rain = np.nansum(data3['flx_ls_rain'],1)
    # flx_conv_rain = np.nansum(data3['flx_conv_rain'],1)
    # flx_conv_rain[flx_conv_rain < 0] = np.nan
    # flx_ls_snow = np.nansum(data3['flx_ls_snow'],1)

    flx_ls_rain = data3['flx_ls_rain'][:,0]         #### just take surface value
    flx_ls_snow = data3['flx_ls_snow'][:,0]             #### assumes all precip which forms at altitudes
                                                        #### above evaporates/sublimes before it reaches
                                                        #### the surface
    #### remove flagged values
    flx_ls_rain[flx_ls_rain < 0] = np.nan
    flx_ls_snow[flx_ls_snow < 0] = np.nan

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(16,12))

    ax  = fig.add_axes([0.07,0.7,0.56,0.22])   # left, bottom, width, height
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # # plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['rainfall_flux'].data*3600, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['rainfall_flux'].data*3600, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], (flx_ls_rain)*3600, color = 'darkorange', label = label3)
        # plt.plot(data3['time'], (flx_conv_rain)*3600, color = 'k', label = label3)
    else:
        plt.plot(data3['time'], data3['rainfall_flux'].data, color = 'darkorange', label = label3)
    plt.title('Rainfall flux [mm/hr]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.4,0.56,0.22])   # left, bottom, width, height
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # # plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['snowfall_flux'].data*3600, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['snowfall_flux'].data*3600, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], (flx_ls_snow)*3600, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['snowfall_flux'].data, color = 'darkorange', label = label3)
    plt.title('Snowfall flux [mm/hr]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.1,0.56,0.22])   # left, bottom, width, height
    ax = plt.gca()
    # plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    # plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    # plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'forestgreen')
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    # else:
    #     plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    # plt.title('surface_net_LW_radiation [W/m2]')
    # ax.set_xlim([doy[0],doy[-1]])
    # plt.xlabel('Day of year')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    ax  = fig.add_axes([0.7,0.7,0.27,0.22])   # left, bottom, width, height
    # yDmax = 0.05
    # plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    # sns.distplot(data1['surface_net_SW_radiation'].data + data1['surface_net_LW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True})
    # sns.distplot(data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, hist=False, color="darkorange", kde_kws={"shade": True})
    # sns.distplot(data2['surface_net_SW_radiation'].data + data2['surface_net_LW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True})
    # sns.distplot(netLW + netSW, hist=False, color="black")
    # plt.title('CRF [W/m2]')
    # plt.xlim([-50,80])
    # plt.ylim([0,yDmax])


    ax  = fig.add_axes([0.7,0.4,0.27,0.22])   # left, bottom, width, height
    # yEmax = 0.06
    # plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    # sns.distplot(data1['surface_net_SW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True}, label = label1)
    # sns.distplot(data3['sfc_net_sw'].data, hist=False, color="darkorange", kde_kws={"shade": True}, label = label3)
    # sns.distplot(data2['surface_net_SW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True}, label = label2)
    # sns.distplot(netSW, hist=False, color="black", label = 'Ice_station')
    # plt.title('surface_net_SW_radiation [W/m2]')
    # plt.legend()
    # plt.xlim([-10,110])
    # plt.ylim([0,yEmax])


    ax  = fig.add_axes([0.7,0.1,0.27,0.22])   # left, bottom, width, height
    # yFmax = 0.12
    # plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    # sns.distplot(data1['surface_net_LW_radiation'].data, hist=False, color="steelblue", kde_kws={"shade": True})
    # sns.distplot(data3['sfc_net_lw'].data, hist=False, color="darkorange", kde_kws={"shade": True})
    # sns.distplot(data2['surface_net_LW_radiation'].data, hist=False, color="forestgreen", kde_kws={"shade": True})
    # sns.distplot(netLW, hist=False, color="black")
    # plt.title('surface_net_LW_radiation [W/m2]')
    # plt.xlim([-80,20])
    # plt.ylim([0,yFmax])

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    fileout = '../FIGS/comparisons/CRF_netSW_netLW_line+PDFS_oden_iceStation_metum_ifs_casim-100.svg'
    # plt.savefig(fileout)
    plt.show()

def plot_Radiosondes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting radiosonde and model temperature profiles:'
    print ''

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

    ###
    temp = np.zeros([len(data1['time'])])
    for i in range(0, len(temp)-1):
        if data1['time'][i] == data1['time'][i+1]:
            continue
        else:
            temp[i] = data1['time'][i]
    ii = np.where(temp != 0.0)      ### picks out where data are non-zero

    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_hrly'] = temp[ii]
    data2['time_hrly'] = temp[ii]
    data3['time_hrly'] = temp[ii]

    #### ---------------------------------------------------------------
    #### save hourly temperature model profiles (using the ii index defined by the time indices)
    #### ---------------------------------------------------------------
    data1['temp_hrly'] = np.squeeze(data1['temperature'][ii,:])
    data2['temp_hrly'] = np.squeeze(data2['temperature'][ii,:])
    data3['temp_hrly'] = np.squeeze(data3['temperature'][ii,:])

    #### ---------------------------------------------------------------
    #### index to only look at altitudes <10km
    #### ---------------------------------------------------------------
    iTim = 0
    iObs = np.where(obs['sondes']['gpsaltitude'][:,iTim] <= 10000)
    iUM = np.where(data1['height'] <= 10000)
    iIFS = np.where(data3['height'][iTim,:] <= 10000)

    # print obs['sondes']['gpsaltitude'][iObs,iTim].shape
    # print data1['height'][iUM].shape
    # print data3['height'][iTim,iIFS].shape
    # (1, 1001)
    # (53,)
    # (1, 58)
        #### so, will interpolate on to um grid (lowest resolution)


    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    print ''
    print 'Defining IFS temperature profile as a function:'
    data3['temp_hrly_UM'] = np.zeros([np.size(data3['time_hrly'][::6],0),len(data1['height'][iUM[0][2:]])])
    for iTim in range(0,np.size(data3['time_hrly'][::6],0)):
        if np.squeeze(data3['height'][iTim,iIFS[0][0]]) <= 0.0:     ### if height data flagged or not present, skip
            continue
        else:
            fnct_IFS = interp1d(np.squeeze(data3['height'][iTim,iIFS]), np.squeeze(data3['temp_hrly'][iTim,iIFS]))
        data3['temp_hrly_UM'][iTim,:] = fnct_IFS(data1['height'][iUM[0][2:]].data)
    print '...'
    print 'IFS(UM Grid) function worked!'
    print '*****'

    #### INTERPOLATION TESTING:
    print data3['temp_hrly_UM'].shape
    print data3['time_hrly'][::6].shape
    plt.plot(data3['temp_hrly_UM'][10,:],data1['height'][iUM[0][2:]])
    plt.plot(np.squeeze(data3['temp_hrly'][10,iIFS]),np.squeeze(data3['height'][10,iIFS]))
    plt.show()

    print ''
    print 'Defining Radiosonde temperature profile as a function:'
    fnct_Obs = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes']['temperature'][iObs,iTim]))
    print '...'
    print 'Applying to UM vertical grid:'
    obs['sondes']['temp_hrly_UM'] = fnct_Obs(data1['height'][iUM[0][2:]].data)
    print '...'
    print 'Sonde(UM Grid) function worked!'
    print '*****'

    #### INTERPOLATION TESTING:
    # print obs['sondes']['temp_hrly_UM']
    # plt.plot(obs['sondes']['temp_hrly_UM'],data1['height'][iUM[0][2:]])
    # plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,iTim]),np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]))
    # plt.show()

    print 'Starting radiosonde figure (quite slow!)...:'
    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    #
    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(12,10))
    #
    # ax  = fig.add_axes([0.1,0.78,0.9,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy'],data1['height'][iUM[0][2:]],obs['sondes']['temp_hrly_UM'], vmin = -25, vmax = 5)
    # plt.ylim([0,4000])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title('Radiosondes, T[degC]')
    #
    # ax  = fig.add_axes([0.1,0.54,0.9,0.17])   # left, bottom, width, height
    # plt.pcolor(data3['time_hrly'][::6],data1['height'][iUM[0][2:]],np.transpose(data3['temp_hrly'][::6,:])-273.15, vmin = -25, vmax = 5)
    # plt.ylim([0,4000])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label3 + ', T[degC]')
    #
    # ax  = fig.add_axes([0.1,0.3,0.9,0.17])   # left, bottom, width, height
    # plt.pcolor(data1['time_hrly'][::6],data1['height'][2:],np.transpose(data1['temp_hrly'][::6,2:])-273.15, vmin = -25, vmax = 5)
    # plt.ylim([0,4000])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label1 + ', T[degC]')
    #
    # ax  = fig.add_axes([0.1,0.06,0.9,0.17])   # left, bottom, width, height
    # plt.pcolor(data2['time_hrly'][::6],data2['height'][2:],np.transpose(data2['temp_hrly'][::6,2:])-273.15, vmin = -25, vmax = 5)
    # plt.ylim([0,4000])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    # plt.title(label2 + ', T[degC]')
    #
    #
    # print '******'
    # print ''
    # print 'Finished plotting! :)'
    # print ''
    #
    # fileout = '../FIGS/comparisons/TemperatureProfiles_sondes_metum_ifs_casim-100.svg'
    # # plt.savefig(fileout)
    # plt.show()

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    # print GlobalStashList
    # print GlobalStashList[0]

    return GlobalStashList

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### only works on laptop for now

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        misc_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        # misc_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/'
    out_dir2 = '5_u-bl661_RA1M_CASIM/OUT_R0/'
    # out_dir3 = 'MET_DATA/'
    out_dir4 = 'OUT_25H/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R0/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param

    print '******'
    print ''
    print 'Identifying .nc file: '
    print ''

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    print 'Loading observations:'
            # -------------------------------------------------------------
            # Which file does what?
            # -------------------------------------------------------------
            #### ice station: net LW / net SW
                    #### obs['ice_station']/mast_radiation_30min_v2.3.mat
            #### obs['foremast']:
                    #### obs['foremast']/ACAS_AO2018_obs['foremast']_30min_v2_0.nc
            #### 7th deck: temperature, surface temperature, RH, downwelling SW, downwelling LW
                    #### 7thDeck/ACAS_AO2018_WX_30min_v2_0.nc

    obs = {}

    print 'Load temporary ice station data from Jutta...'
    obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_wTemp1p5m.nc','r')

    print 'Load ice station data from Jutta...'
    obs['ice_station'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')
            #### mast_radiation_30min_v2.3.mat
            #### flux30_trhwxrel.mat

    print 'Load radiosonde data from Jutta...'
    obs['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print 'Load foremast data from John...'
    obs['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print 'Load 7th deck weather station data from John...'
    obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print '...'

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print '******'
    print ''
    print 'Make stash list for cube read in at ' + time.strftime("%c")
    print ' '
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin cube read in at ' + time.strftime("%c")
    print ' '

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    # tempnames = ['umnsaa_pa012_r0.nc','umnsaa_pb012_r0.nc','umnsaa_pc011_r0.nc','umnsaa_pd011_r0.nc','20180812_oden_metum.nc']
    Aug_names = ['20180813_oden_','20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_']

    Sep_names = ['20180901_oden_','20180902_oden_','20180903_oden_','20180904_oden_',
            '20180905_oden_','20180906_oden_','20180907_oden_','20180908_oden_',
            '20180909_oden_','20180910_oden_','20180911_oden_','20180912_oden_',
            '20180913_oden_','20180914_oden_']

    moccha_names = ['20180813_oden_','20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = []

    doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)          ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):
        filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
        filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
        if out_dir4 == 'OUT_25H/':
            print '***IFS being compared***'
            ifs_flag = True
            filename_um3 = misc_root_dir + out_dir4 + names[i] + 'ecmwf.nc'
        else:
            print '***IFS NOT being compared***'
            filename_um3 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            ifs_flag = False
        print filename_um1
        print filename_um2
        print filename_um3
        print ''

        #### LOAD CUBE
        print 'Loading first run diagnostics:'
        nc1 = Dataset(filename_um1,'r')
        print '...'
        print 'Loading second run diagnostics:'
        nc2 = Dataset(filename_um2,'r')
        print '...'
        print 'Loading third run diagnostics:'
        nc3 = Dataset(filename_um3,'r')
        print '...'
        # -------------------------------------------------------------
        # print 'i = ' + str(i)
        print ''

        #### LOAD IN SPECIFIC DIAGNOSTICS
        # if out_dir == '4_u-bg610_RA2M_CON/OUT_R1/':
        var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux','latent_heat_flux',
            'temp_1.5m', 'rainfall_flux','snowfall_flux']
        var_list2 = var_list1
        if ifs_flag: var_list3 = ['height','temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
            'sfc_temp_2m','flx_ls_rain','flx_conv_rain','flx_ls_snow']
        if not ifs_flag: var_list3 = var_list1

        if i == 0:
            ## ------------------
            #### UM
            ## ------------------
            data1 = {}
            data2 = {}
            data3 = {}
            if month_flag == -1:
                time_um1 = doy[i] + (nc1.variables['forecast_time'][:]/24.0)
                time_um2 = doy[i] + (nc2.variables['forecast_time'][:]/24.0)
                if ifs_flag: time_um3 = doy[i] + (nc3.variables['time'][:]/24.0)
                if not ifs_flag: time_um3 = doy[i] + (nc3.variables['forecast_time'][:]/24.0)
            else:
                time_um1 = float(filename_um1[-16:-14]) + (nc1.variables['forecast_time'][:]/24.0)
                time_um2 = float(filename_um2[-16:-14]) + (nc2.variables['forecast_time'][:]/24.0)
                if ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['time'][:]/24.0)
                if not ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['forecast_time'][:]/24.0)

            ### define height arrays explicitly
            data1['height'] = nc1.variables['height'][:]
            data2['height'] = nc2.variables['height'][:]
            if not ifs_flag: data3['height'] = nc3.variables['height'][:]

            for j in range(0,len(var_list1)):
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) >= 1:
                    data1[var_list1[j]] = nc1.variables[var_list1[j]][:]
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            for j in range(0,len(var_list2)):
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
                    data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
            nc2.close()
            ## ------------------
            #### um3
            ## ------------------
            for j in range(0,len(var_list3)):
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) >= 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
            nc3.close()
        else:
            if month_flag == -1:
                time_um1 = np.append(time_um1, doy[i] + (nc1.variables['forecast_time'][:]/24.0))
                time_um2 = np.append(time_um2, doy[i] + (nc2.variables['forecast_time'][:]/24.0))
                if ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['time'][:]/24.0))
                if not ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['forecast_time'][:]/24.0))
            ## ------------------
            #### UM
            ## ------------------
            for j in range(0,len(var_list1)):
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data1[var_list1[j]] = np.append(data1[var_list1[j]].data,nc1.variables[var_list1[j]][:])
                elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                    data1[var_list1[j]] = np.append(data1[var_list1[j]].data,nc1.variables[var_list1[j]][:],0)
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            for j in range(0,len(var_list2)):
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) == 1:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]].data,nc2.variables[var_list2[j]][:])
                elif np.ndim(nc2.variables[var_list2[j]]) == 2:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]].data,nc2.variables[var_list2[j]][:],0)
            nc2.close()
            ## ------------------
            #### um3 / ifs
            ## ------------------
            for j in range(0,len(var_list3)):
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]].data,nc3.variables[var_list3[j]][:])
                elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]].data,nc3.variables[var_list3[j]][:],0)
            nc3.close()

    #################################################################
    ## save time to dictionary now we're not looping over all diags anymore
    #################################################################
    data1['time'] = time_um1
    data2['time'] = time_um2
    data3['time'] = time_um3

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T'
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'

    label2 = 'undefined_label'
    if out_dir2[:10] == '11_u-bq798': label2 = 'UM_CASIM-100_Meyers'
    if out_dir2[:10] == '10_u-bq791': label2 = 'UM_CASIM-100_Fletcher'
    if out_dir2[:9] == '8_u-bp738': label2 = 'UM_ERAI-GLM'
    if out_dir2[:9] == '7_u-bn068': label2 = 'UM_RA2T'
    if out_dir2[:9] == '6_u-bm410': label2 = 'UM_CASIM-200'
    if out_dir2[:9] == '5_u-bl661': label2 = 'UM_CASIM-100'
    if out_dir2[:9] == '4_u-bg610': label2 = 'UM_RA2M'

    label3 = 'undefined_label'
    if out_dir4 == 'OUT_25H/': label3 = 'ECMWF_IFS'
    if out_dir4[:10] == '11_u-bq798': label3 = 'UM_CASIM-100_Meyers'
    if out_dir4[:10] == '10_u-bq791': label3 = 'UM_CASIM-100_Fletcher'
    if out_dir4[:9] == '8_u-bp738': label3 = 'UM_ERAI-GLM'
    if out_dir4[:9] == '7_u-bn068': label3 = 'UM_RA2T'
    if out_dir4[:9] == '6_u-bm410': label3 = 'UM_CASIM-200'
    if out_dir4[:9] == '5_u-bl661': label3 = 'UM_CASIM-100'
    if out_dir4[:9] == '4_u-bg610': label3 = 'UM_RA2M'

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['sondes'])

    # -------------------------------------------------------------
    # Plot combined column data (5x2 timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_multidate_TS(timem, data, cube, month_flag, missing_files, out_dir)
                ### doesn't matter which cube, just needed for dim_coords

    # -------------------------------------------------------------
    # Plot combined CASIM column data (4x3 timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_multidate_casim_TS(timem, data, cube, month_flag, missing_files, out_dir)
                ### doesn't matter which cube, just needed for dim_coords

    # -------------------------------------------------------------
    # Plot combined timeseries as lineplot
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # figure = plot_line_BLDepth(time_um1, time_um2, data1, data2, cube_um1, cube_um2, month_flag,
    #             missing_files, out_dir1, obs, doy)

    # figure = plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3,
    #     month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # figure = plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # Plot paper figures
    # -------------------------------------------------------------
    # figure = plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_paperRadiation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_Precipitation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    figure = plot_Radiosondes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_ERAI_GLM(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # Plot data (5x2 monthly timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_TS(cube, filename, out_dir)


    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

    #### DIAGNOSTICS TO CHOOSE FROM:

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. sfc_temperature -> sfc_temperature
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

if __name__ == '__main__':

    main()
