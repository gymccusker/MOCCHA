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

def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==12,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==20,data.values[:,1]==9),data.values[:,3]==1))
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

def plot_cartmap(ship_data, cube, hour, grid_filename): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
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
    plt.figure(figsize=(12,10))
    ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([0, 60, 86.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE
    # ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())    ### V LARGE

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
    if np.size(cube[diag].data.shape) == 4:
        iplt.pcolormesh(cube[diag][hour,0,:,:])
    elif np.size(cube[diag].data.shape) == 3:
        iplt.pcolormesh(cube[diag][hour,:,:])
        # iplt.pcolormesh(cube[hour,471:495,240:264])
    elif np.size(cube[diag].data.shape) == 2:
        iplt.pcolormesh(cube[diag][:,:])
    # plt.title(cube[diag].standard_name + ', ' + str(cube[diag].units))
    plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[diag][hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[diag][hour,386:495,211:305])          ### misc
    # qplt.outline(cube[diag][hour,:,:])

    # gridship = gridShipTrack(cube[diag], xoffset, yoffset)

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

def plot_multicontour_TS(cube, filename): #, lon, lat):

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

    for i in range(0,len(cube)):

        ###################################
        ## CHOOSE DIAGNOSTIC
        ###################################
        diag = i
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

        #################################################################
        ## create figure and axes instances
        #################################################################
        plt.subplot(5,2,i+1)
        ax = plt.gca()

        #################################################################
        ## data corrections
        #################################################################
        ### set height limit to consider
        ind = np.where(height<3000)

        ### if mass mixing ratio, *1e3 to change to g/kg
        if cube[diag].var_name[0] == 'q':
            data = np.transpose(np.squeeze(cube[diag].data[:,ind]*1e3))
        else:
            data = np.transpose(np.squeeze(cube[diag].data[:,ind]))

        #################################################################
        ## plot timeseries
        #################################################################
        # plt.contourf(time,height,np.transpose(cube[diag].data))
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
        else:
            plt.set_cmap(mpl_cm.viridis)

        ### title and axes properties
        if cube[diag].var_name[0] == 'q':
            plt.title(cube[diag].var_name + ' [g/kg]')
        else:
            plt.title(cube[diag].var_name + ' [' + str(cube[diag].units) + ']')
        plt.colorbar()
        ax.set_ylim([0, 3000])

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

    fileout = 'FIGS/' + filename[-22:-3] + '.png'
    plt.savefig(fileout)
    plt.show()

def plot_multicontour_multidate_TS(timem, data, cube, month_flag): #, lon, lat):

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
    plt.figure(figsize=(12,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    print data.keys()

    for i in range(0,len(data.keys())):

        ###################################
        ## CHOOSE DIAGNOSTIC
        ###################################
        diag = i
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

        if data.keys()[diag] == 'radr_refl':
            height = cube[0].dim_coords[1].points
            title = cube[0].var_name + ' [dBz]'
        elif data.keys()[diag] == 'pressure':
            height = cube[1].dim_coords[1].points
            title = cube[1].var_name + ' [' + str(cube[1].units) + ']'
        elif data.keys()[diag] == 'temperature':
            height = cube[2].dim_coords[1].points
            title = cube[2].var_name + ' [' + str(cube[2].units) + ']'
        elif data.keys()[diag] == 'uwind':
            height = cube[3].dim_coords[1].points
            title = cube[3].var_name + ' [' + str(cube[3].units) + ']'
        elif data.keys()[diag] == 'cloud_fraction':
            height = cube[4].dim_coords[1].points
            title = cube[4].var_name + ' [' + str(cube[4].units) + ']'
        elif data.keys()[diag] == 'qice':
            height = cube[5].dim_coords[1].points
            title = cube[5].var_name + ' [g/kg]'
        elif data.keys()[diag] == 'qliq':
            height = cube[6].dim_coords[1].points
            title = cube[6].var_name + ' [g/kg]'
        elif data.keys()[diag] == 'vwind':
            height = cube[7].dim_coords[1].points
            title = cube[7].var_name + ' [' + str(cube[7].units) + ']'
        elif data.keys()[diag] == 'q':
            height = cube[8].dim_coords[1].points
            title = cube[8].var_name + ' [g/kg]'
        elif data.keys()[diag] == 'wwind':
            height = cube[9].dim_coords[1].points
            title = cube[9].var_name + ' [' + str(cube[9].units) + ']'

        #################################################################
        ## create figure and axes instances
        #################################################################
        plt.subplot(5,2,i+1)
        ax = plt.gca()

        #################################################################
        ## data corrections
        #################################################################
        ### set height limit to consider
        # ind = np.where(height<5000)

        ### if mass mixing ratio, *1e3 to change to g/kg
        if data.keys()[diag][0] == 'q':
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data*1e3))
            # print 'transposed for plotting' + str(dat.shape)
        else:
            dat = np.transpose(np.squeeze(data[data.keys()[diag]].data))
            # print 'transposed for plotting' + str(dat.shape)

        #################################################################
        ## plot timeseries
        #################################################################
        # plt.contourf(time,height,np.transpose(cube[diag].data))
        plt.pcolormesh(timem, height, dat, vmin = np.nanmin(dat), vmax = np.nanmax(dat))

        #################################################################
        ## set plot properties
        #################################################################
        ### colormaps:
        if data.keys()[diag] == 'wwind':
            plt.set_cmap(mpl_cm.RdBu_r)
        elif data.keys()[diag] == 'uwind':
            plt.set_cmap(mpl_cm.RdBu_r)
        elif data.keys()[diag] == 'vwind':
            plt.set_cmap(mpl_cm.RdBu_r)
        else:
            plt.set_cmap(mpl_cm.viridis)

        plt.title(title)
        plt.colorbar()
        ax.set_ylim([0, 5000])
        ax.set_xlim([15.0, 18.0])

    ### global plot properties
    plt.subplot(5,2,9)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    plt.ylabel('Z [m]')
    plt.subplot(5,2,10)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
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

    if month_flag == 8: fileout = 'FIGS/201808_oden_metum.png'
    plt.savefig(fileout)
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

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
        ### CHOSEN RUN
        out_dir = '4_RA2M_CON/OUT/'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        ### CHOSEN RUN
        out_dir = '4_RA2M_CON/'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ## Flag for individual file or monthly:
    combine = 1

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
    tempnames = ['umnsaa_pa012_r0.nc','umnsaa_pb012_r0.nc','umnsaa_pc011_r0.nc','umnsaa_pd011_r0.nc','20180812_oden_metum.nc']
    names = ['20180812_oden_metum.nc','20180815_oden_metum.nc','20180816_oden_metum.nc','20180817_oden_metum.nc']
    # names = ['umnsaa_pa000.nc','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    if combine == 0:
        filename1 = root_dir + out_dir + names[3]
        print filename1
        print ''

        # date = names[1][0:8]

        #### LOAD CUBE
        if 'var_con' in locals():
            print 'Loading single diagnostic:'
            print var_con
            cube1 = iris.load_cube(filename1, var_con, callback)
        elif 'global_con' in locals():
            print 'Loading multiple diagnostics:'
            # cube = iris.load_cubes(filename1, global_con)
            cube = iris.load(filename1, global_con, callback)

            # -------------------------------------------------------------

        print cube
        print ''

        # -------------------------------------------------------------
        # Plot data (5x2 timeseries)
        # -------------------------------------------------------------
        figure = plot_multicontour_TS(cube, filename1)

    else:
        for i in range(1,len(names)):
            month_flag = 8
            filename = root_dir + out_dir + names[i]
            print filename
            print ''

            print 'Loading multiple diagnostics:'
            cube = iris.load(filename, global_con, callback)

            print cube
            print ''

            if i == 1:
                data = {}
                # data['time'] = []
                # data['time'] = float(filename[-16:-14]) + ((cube[0].dim_coords[0].points)/24.0)
                timem = float(filename[-16:-14]) + ((cube[0].dim_coords[0].points)/24.0)
                for j in range(0,len(cube)):
                    data[cube[j].var_name] = cube[j].data
                # print data[cube[0].var_name]
            else:
                # data['time'] = np.append(data['time'],float(filename[-16:-14]) + ((cube[0].dim_coords[0].points)/24.0))
                timem = np.append(timem,float(filename[-16:-14]) + ((cube[0].dim_coords[0].points)/24.0))
                # print data
                for j in range(0,len(cube)):
                    data[cube[j].var_name] = np.append(data[cube[j].var_name].data,cube[j].data,0)

            print 'Data dict = ' + str(data['radr_refl'].shape)

        # -------------------------------------------------------------
        # Plot combined data (5x2 timeseries)
        # -------------------------------------------------------------
        figure = plot_multicontour_multidate_TS(timem, data, cube, month_flag)
                    ### doesn't matter which cube, just needed for dim_coords

        # -------------------------------------------------------------
        # Plot data (5x2 monthly timeseries)
        # -------------------------------------------------------------
        # figure = plot_multicontour_TS(cube, filename)


    # -------------------------------------------------------------
    # Plot data (single timeseries)
    # -------------------------------------------------------------
    # figure = plot_contour_TS(cube, filename1)
    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    ### select hour to plot
    # hour = 0
    # figure = plot_cartmap(ship_data, cube, hour, grid_filename)#, lon, lat)

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

    #### DIAGNOSTICS TO CHOOSE FROM:

    ### paXXX
    # 0: northward_wind_at_10m / (m s-1)     (time: 8; grid_latitude: 501; grid_longitude: 500)
    # 1: eastward_wind_at_10m / (m s-1)      (time: 8; grid_latitude: 501; grid_longitude: 500)
    # 2: surface_downwelling_SW_radiation / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 3: surface_net_LW_radiation / (W m-2)  (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 4: surface_net_SW_radiation / (W m-2)  (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 5: relative_humidity_at_1.5m / (%)     (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 6: surface_downwelling_LW_radiation / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 7: specific_humidity_at_1.5m / (1)     (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 8: surface_net_SW_radiation / (W m-2)  (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 9: air_temperature_at_1.5m / (K)       (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 10: dew_point_temperature_at_1.5m / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 11: surface_net_LW_radiation / (W m-2)  (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 12: air_pressure_at_sea_level / (Pa)    (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 13: surface_air_pressure / (Pa)         (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 14: surface_temperature / (K)           (time: 8; grid_latitude: 500; grid_longitude: 500)
    # 15: toa_incoming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 16: toa_outgoing_longwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 17: toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)

    #### 12 AUG ONLY - NO FULL NEST DIAGNOSTICS
    # <iris 'Cube' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_longwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_incoming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>]

    ### pbXXX
    # 0: specific_humidity_at_1.5m / (1)     (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 1: cloud_area_fraction_assuming_maximum_random_overlap / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 2: cloud_area_fraction_assuming_random_overlap / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 3: air_temperature_at_1.5m / (K)       (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 4: northward_wind_at_10m / (m s-1)     (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 5: eastward_wind_at_10m / (m s-1)      (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 6: height_of_stratocumulus_cloud_base / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 7: dew_point_temperature_at_1.5m / (K) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 8: total_column_q / (1)                (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 9: turbulent mixing height after boundary layer / (m) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 10: height_of_decoupled_layer_base / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 11: combined_boundary_layer_type / (1)  (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 12: wet_bulb_freezing_level_altitude / (m) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 13: relative_humidity_at_1.5m / (%)     (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 14: large_scale_ice_water_path / (1)    (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 15: large_scale_liquid_water_path / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 16: air_pressure_at_sea_level / (Pa)    (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 17: atmosphere_boundary_layer_thickness / (m) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 18: high_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 19: low_type_cloud_area_fraction / (1)  (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 20: medium_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 21: stratiform_rainfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 22: stratiform_snowfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 23: surface_air_pressure / (Pa)         (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 24: surface_temperature / (K)           (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 25: surface_upward_latent_heat_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 26: surface_upward_sensible_heat_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 27: water_evaporation_amount / (1)      (time: 24; grid_latitude: 94; grid_longitude: 95)


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


    ### pdXXX
    # 0: entrainment_rate_for_surface_mixed_layer / (unknown) (grid_latitude: 25; grid_longitude: 25)
    # 1: entrainment_rate_for_boundary_layer / (unknown) (grid_latitude: 25; grid_longitude: 25)
    # 2: obukhov_length / (unknown)          (grid_latitude: 25; grid_longitude: 25)
    # 3: atmosphere_downward_eastward_stress / (Pa) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 4: atmosphere_downward_northward_stress / (Pa) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 5: turbulent_kinetic_energy / (unknown) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 6: air_pressure / (Pa)                 (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 7: surface_downward_eastward_stress / (Pa) (grid_latitude: 25; grid_longitude: 25)
    # 8: surface_downward_northward_stress / (Pa) (grid_latitude: 25; grid_longitude: 25)
    # 9: surface_upward_water_flux / (kg m-2 s-1) (grid_latitude: 25; grid_longitude: 25)

if __name__ == '__main__':

    main()
