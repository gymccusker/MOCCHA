###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE
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

def plot_multicontour_multidate_casim_TS(timem, data, cube, month_flag, missing_files, out_dir): #, lon, lat):

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
        print 'Diag = '  + str(diag)
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
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "pressure":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "temperature":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "uwind":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "wwind":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "vwind":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "cloud_fraction":
            casim_flag = 0
            dat = []
        elif str(data.keys()[diag]) == "qice":
            casim_flag = 0
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qliq":
            casim_flag = 0
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qnrain":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qnsnow":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qnice":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qprisice":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qnliq":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qngraup":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
        elif str(data.keys()[diag]) == "qrain":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "qgraup":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
        elif str(data.keys()[diag]) == "q":
            casim_flag = 1
            title = str(data.keys()[diag]) + ' [g/kg]'
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))

        print 'dat = ' + str(len(dat))
        # print ''

        if len(dat) > 0:
            print 'Plotting since len(dat)>0'
            # print ''
            print title

            l = l + 1 ## increment index for positive data association
            print 'l = ' +  str(l)
            #################################################################
            ## create figure and axes instances
            #################################################################
            plt.subplot(4,3,l+1)
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
            if month_flag == 8: ax.set_xlim([13.0, 31.0])
            if month_flag == 9: ax.set_xlim([1.0, 15.0])

            # print ''
            print 'Zero out any data from missing files:'
            # print ''
            for mfile in missing_files:
                # mtime = float(mfile[6:8]) + ((cube[0].dim_coords[0].points)/24.0)
                # nans = np.zeros([len(height),len(mtime)])
                # # nans[nans == 0] = np.nan
                # plt.pcolormesh(mtime, height, nans)
                mtime = float(mfile[6:8]) + ((cube[0].dim_coords[0].points)/24.0)
                nans = ax.get_ylim()
                ax.fill_between(mtime, nans[0], nans[-1], facecolor = 'lightgrey', zorder = 3)

    ### global plot properties
    plt.subplot(4,3,9)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    plt.subplot(4,3,10)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    plt.ylabel('Z [m]')
    plt.subplot(4,3,11)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    plt.subplot(4,3,1)
    plt.ylabel('Z [m]')
    plt.subplot(4,3,4)
    plt.ylabel('Z [m]')
    plt.subplot(4,3,7)
    plt.ylabel('Z [m]')
    plt.subplot(4,3,11)
    plt.ylabel('Z [m]')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201808_oden_metum_casim.png'
            plt.savefig(fileout, dpi=300)
    if month_flag == 9:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201809_oden_metum_casim.png'
            plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_TSa(time_um1, time_um2, time_um3, data1d_um1, data1d_um2, data1d_um3, cube_um1, cube_um2, cube_um3, month_flag,
            missing_files, out_dir1, out_dir2, out_dir4, cube_obs, doy): #, lon, lat):

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
            hspace = 0.4, wspace = 0.15)

    #################################################################
    ## sort out observations' timestamp
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

    datenums_temp = cube_obs[0].dim_coords[0].points
    timestamps_temp = pd.to_datetime(datenums_temp-719529, unit='D')
    time_temp = timestamps_temp.dayofyear + (timestamps_temp.hour / 24.0) + (timestamps_temp.minute / 1440.0) + (timestamps_temp.second / 86400.0)

    datenums_tice = cube_obs[4].dim_coords[0].points
    timestamps_tice = pd.to_datetime(datenums_tice-719529, unit='D')
    time_tice = timestamps_tice.dayofyear + (timestamps_tice.hour / 24.0) + (timestamps_tice.minute / 1440.0) + (timestamps_tice.second / 86400.0)

    datenums_radice = cube_obs[1].dim_coords[0].points
    timestamps_radice = pd.to_datetime(datenums_radice-719529, unit='D')
    time_radice = timestamps_radice.dayofyear + (timestamps_radice.hour / 24.0) + (timestamps_radice.minute / 1440.0) + (timestamps_radice.second / 86400.0)

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT2/':
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
    zeros = np.zeros(len(time_um2))
    # t1 = 12
    # # t2 = 46
    # t2 = 66
    # # t3 = 61
    # # t1 = 28
    # # t2 = 72
    # # t3 = 84
    # t3 = 86

    t1 = 11
    # t2 = 13
    t3 = 45
    # t4 = 61
    # t5 = 72    ###good
    # t6 = 83
    # t7 = 86


    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)
    ax = plt.gca()
    plt.plot(time_um1, data1d_um1['sfc_pressure'].data/1e2, label = 'CASIM-100')
    plt.plot(time_um2, data1d_um2['sfc_pressure'].data/1e2, label = 'CASIM-200')
    if ifs_flag == True:
        plt.plot(time_um3, data1d_um3['sfc_pressure'].data/1e2, color = 'grey', label = 'IFS')
    else:
        plt.plot(time_um3, data1d_um3['sfc_pressure'].data/1e2, label = 'CASIM-200')
    plt.title('sfc_pressure [hPa]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,cube_obs[4].data + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(time_um1, data1d_um1['temp_1.5m'].data, label = '1.5m')
    ax1.plot(time_um2, data1d_um2['temp_1.5m'].data)#, label = '2m')
    if ifs_flag == True:
        ax1.plot(time_um3, data1d_um3['sfc_temp_2m'].data, color = 'grey', label = '2m')
    else:
        ax1.plot(time_um3, data1d_um3['temp_1.5m'].data)#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    lims = plt.ylim()
    plt.plot([time_um1[t1],time_um1[t1]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t2],time_um1[t2]],[lims[0],lims[-1]],'--',color = 'purple')
    plt.plot([time_um1[t3],time_um1[t3]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t4],time_um1[t4]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t5],time_um1[t5]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t6],time_um1[t6]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t7],time_um1[t7]],[lims[0],lims[-1]],'--',color = 'purple')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    data1d_um1['surface_net_SW_radiation'].data[data1d_um1['surface_net_SW_radiation'].data == 0] = np.nan
    data1d_um2['surface_net_SW_radiation'].data[data1d_um2['surface_net_SW_radiation'].data == 0] = np.nan
    if out_dir4 != 'OUT2/': data1d_um3['surface_net_SW_radiation'].data[data1d_um3['surface_net_SW_radiation'].data == 0] = np.nan

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(time_um2, zeros,'--', color='lightgrey')
    plt.plot(time_radice,(cube_obs[7].data - cube_obs[8].data), color = 'black', label = 'obs: ice')
    plt.plot(time_um1, data1d_um1['surface_net_SW_radiation'].data)
    plt.plot(time_um2, data1d_um2['surface_net_SW_radiation'].data)
    if ifs_flag == True:
        plt.plot(time_um3, data1d_um3['sfc_net_sw'].data, color = 'grey')
    else:
        plt.plot(time_um3, data1d_um3['surface_net_SW_radiation'].data)
    plt.title('surface_net_SW_radiation [W/m2]')
    lims = plt.ylim()
    plt.plot([time_um1[t1],time_um1[t1]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t2],time_um1[t2]],[lims[0],lims[-1]],'--',color = 'purple')
    plt.plot([time_um1[t3],time_um1[t3]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t4],time_um1[t4]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t5],time_um1[t5]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t6],time_um1[t6]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t7],time_um1[t7]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])


    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(time_um2, zeros,'--', color='lightgrey')
    plt.plot(time_radice,(cube_obs[1].data - cube_obs[2].data), color = 'black', label = 'obs: ice')
    plt.plot(time_um1, data1d_um1['surface_net_LW_radiation'].data)
    plt.plot(time_um2, data1d_um2['surface_net_LW_radiation'].data)
    if ifs_flag == True:
        plt.plot(time_um3, data1d_um3['sfc_net_lw'].data, color = 'grey')
    else:
        plt.plot(time_um3, data1d_um3['surface_net_LW_radiation'].data)
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    # plt.subplot(3,3,6)
    # ax = plt.gca()
    # plt.plot(time_um, data1d_um['snowfall_flux'].data)
    # plt.plot(time_um2, data1d_um2['sfc_ls_snow'].data)
    # plt.title('sfc_snow_amount [kg/m2]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,7)
    # ax = plt.gca()
    # plt.plot(time_um, data1d_um['rainfall_flux'].data)
    # plt.plot(time_um2, data1d_um2['sfc_ls_rain'].data)
    # plt.title('sfc_rain_amount [kg/m2]')
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(time_um2, zeros,'--', color='lightgrey')
    plt.plot(time_um1, data1d_um1['sensible_heat_flux'].data)
    plt.plot(time_um2, data1d_um2['sensible_heat_flux'].data)# * -1.0)
    if ifs_flag == True:
        plt.plot(time_um3, data1d_um3['sfc_down_sens_heat_flx'].data * -1.0, color = 'grey')
    else:
        plt.plot(time_um3, data1d_um3['sensible_heat_flux'].data)# * -1.0)
    plt.title('sensible_heat_flux [W/m2]')
    lims = plt.ylim()
    plt.plot([time_um1[t1],time_um1[t1]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t2],time_um1[t2]],[lims[0],lims[-1]],'--',color = 'purple')
    plt.plot([time_um1[t3],time_um1[t3]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t4],time_um1[t4]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t5],time_um1[t5]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t6],time_um1[t6]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t7],time_um1[t7]],[lims[0],lims[-1]],'--',color = 'purple')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(time_um2, zeros,'--', color='lightgrey')
    plt.plot(time_um1, data1d_um1['latent_heat_flux'].data)
    plt.plot(time_um2, data1d_um2['latent_heat_flux'].data)# * -1.0)
    if ifs_flag == True:
        plt.plot(time_um3, data1d_um3['sfc_down_lat_heat_flx'].data * -1.0, color = 'grey')
    else:
        plt.plot(time_um3, data1d_um3['latent_heat_flux'].data)# * -1.0)
    plt.title('latent_heat_flux [W/m2]')
    lims = plt.ylim()
    plt.plot([time_um1[t1],time_um1[t1]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t2],time_um1[t2]],[lims[0],lims[-1]],'--',color = 'purple')
    plt.plot([time_um1[t3],time_um1[t3]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t4],time_um1[t4]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t5],time_um1[t5]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t6],time_um1[t6]],[lims[0],lims[-1]],'--',color = 'purple')
    # plt.plot([time_um1[t7],time_um1[t7]],[lims[0],lims[-1]],'--',color = 'purple')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    # print ''
    # print 'Zero out any data from missing files:'
    # print ''
    # for mfile in missing_files:
    #     mtime = float(mfile[6:8]) + ((cube_um[0].dim_coords[0].points)/24.0)
    #     nans = ax.get_ylim()
    #     ax.fill_between(mtime, nans[0], nans[-1], facecolor = 'lightgrey', zorder = 3)


    ### global plot properties
    # plt.subplot(3,2,4)
    # if month_flag == 8: plt.xlabel('Day of month [Aug]')
    # if month_flag == 9: plt.xlabel('Day of month [Sep]')
    # if month_flag == -1: plt.xlabel('Day of year')
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

    # if month_flag == 8:
    #     if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
    #         if out_dir4 in locals():
    #             fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir4[:9] + '201808_oden_metum_TS.png'
    #         else:
    #             fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201808_oden_metum_TS.png'
    #     elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
    #         fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201808_oden_metum_TS.png'
    # if month_flag == 9:
    #     if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
    #         fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201809_oden_metum_TS.png'
    #     elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
    #         fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201809_oden_metum_TS.png'
    # if month_flag == -1:
    #     if out_dir1[:20] == '5_u-bl661_RA1M_CASIM':
    #         if out_dir2[:20] == '6_u-bm410_RA1M_CASIM':
    #             if out_dir4 == 'OUT2/':
    #                 fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] + '_oden_metum_ifs_casim_TSa.svg'
    #             else:
    #                 fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_casim_TSa.png'
    #         elif out_dir2[:9] == '4_u-bg410':
    #             if out_dir4 == 'OUT2/':
    #                 fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_ifs_casim_TSa.svg'
    #         else:
    #             fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_oden_metum_casim_TSa.svg'
    #     if out_dir2[:20] == '5_u-bl661_RA1M_CASIM':
    #         if out_dir4[:20] == '6_u-bm410_RA1M_CASIM':
    #             fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_' + out_dir4[:9] + '_oden_metum_casim-100_200_TSa.png'
    #         elif out_dir4 == 'OUT2/':
    #             fileout[:20] = '../FIGS/comparisons/' + out_dir2[:20] + '_oden_metum_ifs_casim-100_TSa.png'
        # elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
        #     fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] +'_oden_metum_casim_TSa.png'
    fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] + '_oden_metum_ifs_casim_TSa.svg'
    plt.savefig(fileout, dpi=400)
    plt.show()

def plot_BL_profiles(time_um1, time_um2, time_um3, data1d_um1, data1d_um2, data1d_um3, cube_um1, cube_um2, cube_um3, month_flag,
            missing_files, out_dir1, out_dir2, out_dir4, cube_obs, doy): #, lon, lat):

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
    print 'Plotting BL profile comparisons:'
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
            hspace = 0.4, wspace = 0.15)

    #################################################################
    ## sort out observations' timestamp
    #################################################################
    datenums_temp = cube_obs[0].dim_coords[0].points
    timestamps_temp = pd.to_datetime(datenums_temp-719529, unit='D')
    time_temp = timestamps_temp.dayofyear + (timestamps_temp.hour / 24.0) + (timestamps_temp.minute / 1440.0) + (timestamps_temp.second / 86400.0)

    datenums_tice = cube_obs[4].dim_coords[0].points
    timestamps_tice = pd.to_datetime(datenums_tice-719529, unit='D')
    time_tice = timestamps_tice.dayofyear + (timestamps_tice.hour / 24.0) + (timestamps_tice.minute / 1440.0) + (timestamps_tice.second / 86400.0)

    datenums_radice = cube_obs[1].dim_coords[0].points
    timestamps_radice = pd.to_datetime(datenums_radice-719529, unit='D')
    time_radice = timestamps_radice.dayofyear + (timestamps_radice.hour / 24.0) + (timestamps_radice.minute / 1440.0) + (timestamps_radice.second / 86400.0)


    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
            if out_dir4 in locals():
                fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir4[:9] + '201808_oden_metum_TS.png'
            else:
                fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201808_oden_metum_TS.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201808_oden_metum_TS.png'
    if month_flag == 9:
        if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201809_oden_metum_TS.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201809_oden_metum_TS.png'
    if month_flag == -1:
        if out_dir2[:20] == '5_u-bl661_RA1M_CASIM':
            if out_dir4[:20] == '6_u-bm410_RA1M_CASIM':
                fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_' + out_dir4[:9] + '_oden_metum_casim-100_200_BLprofiles.png'
            else:
                fileout = '../FIGS/comparisons/' + out_dir2[:20] + '_oden_metum_casim-200_BLprofiles.png'
        # elif out_dir2[:18] == '4_u-bg610_RA2M_CON':
        #     fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_casim_TS.png'
    plt.savefig(fileout, dpi=400)
    plt.show()




def plot_line_BLDepth(time_um1, time_um2, data1d_um1, data1d_um2, cube_um1, cube_um2, month_flag,
        missing_files, out_dir1, cube_obs, doy): #, lon, lat):

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
    plt.figure(figsize=(8,5))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
    #         hspace = 0.4, wspace = 0.15)

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
    zeros = np.zeros(len(time_um2))

    #################################################################
    ## create figure and axes instances
    #################################################################

    ax = plt.gca()
    plt.plot(time_um1, data1d_um1['bl_depth'].data, label = 'UM')
    plt.plot(time_um2, data1d_um1['bl_depth'].data, label = 'CASIM-100')
    plt.legend()
    plt.title('BL_depth [m]')
    if month_flag == 8:
        ax.set_xlim([13.0, 31.0])
        plt.xlabel('Day of month [Aug]')
    if month_flag == 9:
        ax.set_xlim([1.0, 15.0])
        plt.xlabel('Day of month [Sep]')
    if month_flag == -1:
        ax.set_xlim([doy[0],doy[-1]])
        plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201808_oden_metum_ecmwf_BLDepth.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201808_oden_metum_ecmwf_BLDepth.png'
    if month_flag == 9:
        if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201809_oden_metum_ecmwf_BLDepth.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201809_oden_metum_ecmwf_BLDepth.png'
    if month_flag == -1:
        if out_dir1[:18] == '6_u-bm410_RA1M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:20] + '_oden_metum_casim-200_BLDepth.png'
        if out_dir1[:18] == '5_u-bl661_RA1M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:20] + '_oden_metum_casim-100_BLDepth.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_casim_BLDepth.png'
    plt.savefig(fileout, dpi=400)
    plt.show()

def plot_line_RAD(time_um1, time_um2, data1d_um1, data1d_um2, cube_um1, cube_um2,
    month_flag, missing_files, out_dir1, out_dir2, cube_obs, doy): #, lon, lat):

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
    ## sort out observations' timestamp
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

    datenums_temp = cube_obs[0].dim_coords[0].points
    timestamps_temp = pd.to_datetime(datenums_temp-719529, unit='D')
    time_temp = timestamps_temp.dayofyear + (timestamps_temp.hour / 24.0) + (timestamps_temp.minute / 1440.0) + (timestamps_temp.second / 86400.0)

    datenums_radice = cube_obs[1].dim_coords[0].points
    timestamps_radice = pd.to_datetime(datenums_radice-719529, unit='D')
    time_radice = timestamps_radice.dayofyear + (timestamps_radice.hour / 24.0) + (timestamps_radice.minute / 1440.0) + (timestamps_radice.second / 86400.0)

    #################################################################
    ## create figure and axes instances
    #################################################################

    plt.subplot(211)
    ax = plt.gca()
    plt.plot(time_um1, data1d_um1['temp_1.5m'].data - 273.15, color = 'r', label = 'Oper')
    plt.plot(time_um2, data1d_um2['temp_1.5m'].data - 273.15, color = 'b', label = 'CASIM-100')
    plt.plot(time_temp,cube_obs[0].data - 273.15, color = 'black', label = 'Observations')
    plt.legend()
    plt.title('Temperature [$^{o}C$]')
    plt.ylim([260 - 273,275 - 273])
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
    data1d_um1['surface_net_SW_radiation'].data[data1d_um1['surface_net_SW_radiation'].data == 0] = np.nan
    data1d_um2['surface_net_SW_radiation'].data[data1d_um2['surface_net_SW_radiation'].data == 0] = np.nan
    plt.plot(time_um1, data1d_um1['surface_net_SW_radiation'].data, color = 'r', label = 'Oper')
    plt.plot(time_um2, data1d_um2['surface_net_SW_radiation'].data, color = 'b', label = 'CASIM-100')
    plt.plot(time_radice,(cube_obs[7].data - cube_obs[8].data), color = 'black', label = 'Observations')
    # plt.legend()
    plt.title('Net SW radiation [W/m2]')
    # plt.ylim([260,275])
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
    # data1d_um['surface_net_LW_radiation'].data[data1d_um['surface_net_LW_radiation'].data == 0] = np.nan
    # plt.plot(time_um, data1d_um['surface_net_LW_radiation'].data, color = 'r', label = 'MetUM')
    # plt.plot(time_radice,(cube_obs[1].data - cube_obs[2].data), color = 'black', label = 'Observations')
    # # plt.legend()
    # plt.title('Net SW radiation [W/m2]')
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

    # if month_flag == 8:
    #     if out_dir[:18] == '5_u-bl616_RA2M_CAS':
    #         fileout = '../FIGS/UM/' + out_dir[:21] + '201808_oden_metum_temp.png'
    #     elif out_dir[:18] == '4_u-bg610_RA2M_CON':
    #         fileout = '../FIGS/UM/' + out_dir[:19] + '201808_oden_metum_temp.png'
    # if month_flag == 9:
    #     if out_dir[:18] == '5_u-bl616_RA2M_CAS':
    #         fileout = '../FIGS/UM/' + out_dir[:21] + '201809_oden_metum_temp.png'
    #     elif out_dir[:18] == '4_u-bg610_RA2M_CON':
    #         fileout = '../FIGS/UM/' + out_dir[:19] + '201809_oden_metum_temp.png'
    if month_flag == -1:
        if out_dir2[0:20] == '6_u-bm410_RA1M_CASIM':
            fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-200_tempoC_SW.png'
        if out_dir2[0:20] == '5_u-bl661_RA1M_CASIM':
            fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-100_tempoC_SW.png'
        elif out_dir2[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_tempoC_SW.png'
    print 'Saving as: ' + fileout
    plt.savefig(fileout, dpi=400)
    plt.show()

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
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir1 = '5_u-bl661_RA1M_CASIM/OUT/'
    out_dir2 = '6_u-bm410_RA1M_CASIM/OUT/'
    out_dir3 = 'MET_DATA/'
    out_dir4 = 'OUT2/'

    ### TESTING/domain_tests/umnsaa_pa000
    ### 4_u-bg610_RA2M_CON/OUT_R1/papbpc_combined/
    ### 5_u-bl661_RA1M_CASIM/OUT/
    ### 6_u-bm410_RA1M_CASIM/OUT/

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
    filename_obs = obs_root_dir + out_dir3 + 'MetData_Gillian_wTemp1p5m.nc'
    cube_obs = iris.load(filename_obs)#, global_con, callback)
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

    moccha_names = [#'20180813_oden_','20180814_oden_','20180815_oden_','20180816_oden_',
            # '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            # '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            # '20180825_oden_','20180826_oden_','20180827_oden_',
            '20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_']#,'20180903_oden_','20180904_oden_','20180905_oden_',
            # '20180906_oden_','20180907_oden_']#,'20180908_oden_','20180909_oden_',
            # '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = []

    # doy = np.arange(225,259)        ## set DOY for full moccha figures
    doy = np.arange(240,247)        ## set DOY for subset of moccha figures
    # doy = np.arange(240,251)        ## set DOY for subset of moccha figures

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    if combine == 0:
        name = '20180813_oden_'
        filename_um1 = um_root_dir + out_dir1 + name + 'metum.nc'
        filename_um2 = um_root_dir + out_dir2 + name + 'metum.nc'
        if out_dir4 == 'OUT2/':
            print '***IFS being compared***'
            filename_um3 = misc_root_dir + out_dir4 + name + 'ecmwf.nc'
        else:
            filename_um3 = misc_root_dir + out_dir4 + name + 'metum.nc'
        print filename_um1
        print filename_um2
        print filename_um3
        print ''

        #### LOAD CUBE
        print 'Loading first run diagnostics:'
        cube_um1 = iris.load(filename_um1)#, global_con, callback)
        print '...'
        print 'Loading second run diagnostics:'
        cube_um2 = iris.load(filename_um2)#, global_con, callback)
        print '...'
        print 'Loading third run diagnostics:'
        cube_um3 = iris.load(filename_um3)#, global_con, callback)
        # -------------------------------------------------------------
        print '...'

        print cube_um1
        print ''
        print cube_um2
        print ''
        print cube_um3
        print ''


    else:
        for i in range(0,len(names)):
            filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
            filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
            if out_dir4 == 'OUT2/':
                print '***IFS being compared***'
                filename_um3 = misc_root_dir + out_dir4 + names[i] + 'ecmwf.nc'
            else:
                filename_um3 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            print filename_um1
            print filename_um2
            print filename_um3
            print ''

            #### LOAD CUBE
            print 'Loading first run diagnostics:'
            cube_um1 = iris.load(filename_um1)#, global_con, callback)
            print '...'
            print 'Loading second run diagnostics:'
            cube_um2 = iris.load(filename_um2)#, global_con, callback)
            print '...'
            print 'Loading third run diagnostics:'
            cube_um3 = iris.load(filename_um3)#, global_con, callback)
            print '...'
            # -------------------------------------------------------------
            # print 'i = ' + str(i)
            print ''

            if i == 0:
                ## ------------------
                #### UM
                ## ------------------
                data_um1 = {}
                data1d_um1 = {}
                data_um2 = {}
                data1d_um2 = {}
                data_um3 = {}
                data1d_um3 = {}
                if month_flag == -1:
                    time_um1 = doy[i] + ((cube_um1[0].dim_coords[0].points)/24.0)
                    time_um2 = doy[i] + ((cube_um2[0].dim_coords[0].points)/24.0)
                    time_um3 = doy[i] + ((cube_um3[0].dim_coords[0].points)/24.0)
                else:
                    time_um1 = float(filename_um1[-16:-14]) + ((cube_um1[0].dim_coords[0].points)/24.0)
                    time_um2 = float(filename_um2[-16:-14]) + ((cube_um2[0].dim_coords[0].points)/24.0)
                    time_um3 = float(filename_um3[-16:-14]) + ((cube_um3[0].dim_coords[0].points)/24.0)
                for j in range(0,len(cube_um1)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um1[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um1[j].data.shape) == 24:  # 1d timeseries only
                        data1d_um1[cube_um1[j].var_name] = cube_um1[j].data
                    else:                                   # 2d column data
                        data_um1[cube_um1[j].var_name] = cube_um1[j].data
                ## ------------------
                #### um2
                ## ------------------
                for j in range(0,len(cube_um2)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um2[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um2[j].data.shape) == 24:  # 1d timeseries only
                        data1d_um2[cube_um2[j].var_name] = cube_um2[j].data
                    else:                                   # 2d column data
                        data_um2[cube_um2[j].var_name] = cube_um2[j].data
                ## ------------------
                #### um3
                ## ------------------
                for j in range(0,len(cube_um3)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um3[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um3[j].data.shape) == 24:  # 1d timeseries only
                        data1d_um3[cube_um3[j].var_name] = cube_um3[j].data
                    else:                                   # 2d column data
                        data_um3[cube_um3[j].var_name] = cube_um3[j].data

            else:
                if month_flag == -1:
                    time_um1 = np.append(time_um1, doy[i] + ((cube_um1[0].dim_coords[0].points)/24.0))
                    time_um2 = np.append(time_um2, doy[i] + ((cube_um2[0].dim_coords[0].points)/24.0))
                    time_um3 = np.append(time_um3, doy[i] + ((cube_um3[0].dim_coords[0].points)/24.0))
                else:
                    time_um1 = np.append(time_um1,float(filename_um1[-16:-14]) + ((cube_um1[0].dim_coords[0].points)/24.0))
                    time_um2 = np.append(time_um2,float(filename_um2[-16:-14]) + ((cube_um2[0].dim_coords[0].points)/24.0))
                    time_um3 = np.append(time_um3,float(filename_um3[-16:-14]) + ((cube_um3[0].dim_coords[0].points)/24.0))
                ## ------------------
                #### UM
                ## ------------------
                for j in range(0,len(cube_um1)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um1[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um1[j].data.shape) == 24:
                        data1d_um1[cube_um1[j].var_name] = np.append(data1d_um1[cube_um1[j].var_name].data,cube_um1[j].data)
                    else:
                        data_um1[cube_um1[j].var_name] = np.append(data_um1[cube_um1[j].var_name].data,cube_um1[j].data,0)
                ## ------------------
                #### um2
                ## ------------------
                for j in range(0,len(cube_um2)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um2[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um2[j].data.shape) == 24:
                        data1d_um2[cube_um2[j].var_name] = np.append(data1d_um2[cube_um2[j].var_name].data,cube_um2[j].data)
                    else:
                        data_um2[cube_um2[j].var_name] = np.append(data_um2[cube_um2[j].var_name].data,cube_um2[j].data,0)
                ## ------------------
                #### um3
                ## ------------------
                for j in range(0,len(cube_um3)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um3[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um3[j].data.shape) == 24:
                        data1d_um3[cube_um3[j].var_name] = np.append(data1d_um3[cube_um3[j].var_name].data,cube_um3[j].data)
                    else:
                        data_um3[cube_um3[j].var_name] = np.append(data_um3[cube_um3[j].var_name].data,cube_um3[j].data,0)


        # -------------------------------------------------------------
        # Plot combined column data (5x2 timeseries)
        # -------------------------------------------------------------
        # np.save('working_data', data)
        # figure = plot_multicontour_multidate_TS(timem, data, cube, month_flag, missing_files, out_dir)
                    ### doesn't matter which cube, just needed for dim_coords

        # -------------------------------------------------------------
        # Plot combined CASIM column data (4x3 timeseries)
        # -------------------------------------------------------------
        # np.save('working_data', data)
        # figure = plot_multicontour_multidate_casim_TS(timem, data, cube, month_flag, missing_files, out_dir)
                    ### doesn't matter which cube, just needed for dim_coords

        # -------------------------------------------------------------
        # Plot combined timeseries as lineplot
        # -------------------------------------------------------------
        figure = plot_line_TSa(time_um1, time_um2, time_um3, data1d_um1, data1d_um2, data1d_um3, cube_um1, cube_um2, cube_um3, month_flag,
                    missing_files, out_dir1, out_dir2, out_dir4, cube_obs, doy)

        # figure = plot_line_BLDepth(time_um1, time_um2, data1d_um1, data1d_um2, cube_um1, cube_um2, month_flag,
        #             missing_files, out_dir1, cube_obs, doy)

        # figure = plot_line_RAD(time_um1, time_um2, data1d_um1, data1d_um2, cube_um1, cube_um2,
        #     month_flag, missing_files, out_dir1, out_dir2, cube_obs, doy)

        # figure = plot_BL_profiles(time_um1, time_um2, time_um3, data1d_um1, data1d_um2, data1d_um3, cube_um1, cube_um2, cube_um3, month_flag,
        #             missing_files, out_dir1, out_dir2, out_dir4, cube_obs, doy)

        data_um1['time'] = time_um1
        data_um2['time'] = time_um2
        data_um3['time'] = time_um3
        np.save('working_data1', data_um1)
        np.save('working_data2', data_um2)
        np.save('working_data3', data_um3)

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
