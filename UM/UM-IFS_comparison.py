###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE
###
###

# from __future__ import print_function
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

def plot_line_TSa(time_um, time_ifs, data1d_um, data1d_ifs, cube_um, cube_ifs, month_flag, missing_files, out_dir): #, lon, lat):

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

    # print data1d_um.keys()

    # UM -> IFS comparisons:
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
    zeros = np.zeros(len(time_ifs))

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)
    ax = plt.gca()
    plt.plot(time_um, data1d_um['sfc_pressure'].data, label = 'UM')
    plt.plot(time_ifs, data1d_ifs['sfc_pressure'].data, label = 'IFS')
    plt.title('sfc_pressure [Pa]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,2)
    ax = plt.gca()
    plt.plot(time_um, data1d_um['bl_depth'].data)
    plt.plot(time_ifs, data1d_ifs['sfc_bl_height'].data)
    plt.title('BL_depth [m]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(time_ifs, zeros,'k--')
    plt.plot(time_um, data1d_um['surface_net_SW_radiation'].data)
    plt.plot(time_ifs, data1d_ifs['sfc_net_sw'].data)
    plt.title('surface_net_SW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(time_ifs, zeros,'k--')
    plt.plot(time_um, data1d_um['surface_net_LW_radiation'].data)
    plt.plot(time_ifs, data1d_ifs['sfc_net_lw'].data)
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,5)
    # ax = plt.gca()
    # plt.plot(time_um, data1d_um['temp_1.5m'].data, label = '1.5m')
    # plt.plot(time_ifs, data1d_ifs['sfc_temp_2m'].data, label = '2m')
    # plt.title('near-sfc_temperature [K]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,6)
    # ax = plt.gca()
    # plt.plot(time_um, data1d_um['snowfall_flux'].data)
    # plt.plot(time_ifs, data1d_ifs['sfc_ls_snow'].data)
    # plt.title('sfc_snow_amount [kg/m2]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,7)
    # ax = plt.gca()
    # plt.plot(time_um, data1d_um['rainfall_flux'].data)
    # plt.plot(time_ifs, data1d_ifs['sfc_ls_rain'].data)
    # plt.title('sfc_rain_amount [kg/m2]')
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(time_ifs, zeros,'k--')
    plt.plot(time_um, data1d_um['sensible_heat_flux'].data)
    plt.plot(time_ifs, data1d_ifs['sfc_down_sens_heat_flx'].data * -1.0)
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(time_ifs, zeros,'k--')
    plt.plot(time_um, data1d_um['latent_heat_flux'].data)
    plt.plot(time_ifs, data1d_ifs['sfc_down_lat_heat_flx'].data * -1.0)
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([225.0, 258.0])

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

    if month_flag == 8:
        if out_dir1[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201808_oden_metum_ecmwf_1Da.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201808_oden_metum_ecmwf_1Da.png'
    if month_flag == 9:
        if out_dir1[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:21] + '201809_oden_metum_ecmwf_1Da.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:19] + '201809_oden_metum_ecmwf_1Da.png'
    if month_flag == -1:
        if out_dir1[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/comparisons/' + out_dir1[:20] + '_oden_metum_ecmwf_1Da.png'
        elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_ecmwf_1Da.png'
    plt.savefig(fileout, dpi=400)
    plt.show()

def plot_line_TSb(timem, data, cube, month_flag, missing_files, out_dir): #, lon, lat):

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
    plt.rc('xtick',labelsize=SMALL_SIZE)
    plt.rc('ytick',labelsize=SMALL_SIZE)
    plt.rc('legend',fontsize=SMALL_SIZE)
    plt.figure(figsize=(16,12))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.15)

    print data.keys()

    l = -1
    for diag in range(0,len(data.keys())):

        ###################################
        ## CHOOSE DIAGNOSTIC
        ###################################
        print ''
        print 'Diag is: '
        print data.keys()[diag]

                                            # [u'vis_1.5m',
                                            #  u'surface_net_LW_radiation',
        #  u'bl_type',
                                            #  u'fogfrac_1.5m',
                                            #  u'surface_net_SW_radiation',
                                            #  u'temp_1.5m',
        #  u'snowfall_flux',
        #  u'h_sc_cloud_base',
        #  u'h_decoupled_layer_base',
        #  u'latent_heat_flux',
        #  u'high_cloud',
                                            #  u'sfc_temperature',
        #  u'IWP',
        #  u'total_column_q',
        #  u'bl_depth',
        #  u'LWP',
                                            #  u'rh_1.5m',
        #  u'medium_cloud',
        #  u'sensible_heat_flux',
                                            #  u'sfc_pressure',
        #  u'rainfall_flux',
        #  u'low_cloud']

        ###################################
        ## DEFINE DIMENSIONS COORDS/UNITS DEPENDING ON DIAG
        ###################################

        dat = []

        if str(data.keys()[diag]) == 'LWP':
            dat = data[data.keys()[diag]].data*1e3
            title = str(data.keys()[diag]) + ' [g/m2]'
        elif str(data.keys()[diag]) == 'IWP':
            dat = data[data.keys()[diag]].data*1e3
            title = str(data.keys()[diag]) + ' [g/m2]'
        elif str(data.keys()[diag]) == 'total_column_q':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [kg/m2]'

        elif str(data.keys()[diag]) == 'rainfall_flux':
            dat = data[data.keys()[diag]].data*3600
            title = str(data.keys()[diag]) + ' [mm/hr]'
        elif str(data.keys()[diag]) == 'snowfall_flux':
            dat = data[data.keys()[diag]].data*3600
            title = str(data.keys()[diag]) + ' [mm/hr]'

        elif str(data.keys()[diag]) == 'bl_type':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'
        elif str(data.keys()[diag]) == 'bl_depth':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'

        elif str(data.keys()[diag]) == 'latent_heat_flux':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [W/m2]'
        elif str(data.keys()[diag]) == 'sensible_heat_flux':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [W/m2]'

        elif str(data.keys()[diag]) == 'high_cloud':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'
        elif str(data.keys()[diag]) == 'medium_cloud':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'
        elif str(data.keys()[diag]) == 'low_cloud':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'

        elif str(data.keys()[diag]) == 'h_sc_cloud_base':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [m]'
        elif str(data.keys()[diag]) == 'h_decoupled_layer_base':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [m]'

        #################################################################
        ## create figure and axes instances
        #################################################################
        if len(dat) > 0:
            l = l + 1 ## increment index for positive data association
            plt.subplot(5,3,l+1)
            # print 'l = ' + str(l)
            print title
            ax = plt.gca()

            #################################################################
            ## plot timeseries
            #################################################################

            plt.plot(timem, dat)
            if np.logical_or(np.nanmin(dat) != np.nan, np.nanmax(dat) != np.nan):
                plt.ylim([np.nanmin(dat),np.nanmax(dat)])
            plt.title(title)

            if month_flag == 8: ax.set_xlim([13.0, 31.0])
            if month_flag == 9: ax.set_xlim([1.0, 15.0])
            if month_flag == -1: ax.set_xlim([225.0, 258.0])

            print ''
            print 'Zero out any data from missing files:'
            print ''
            for mfile in missing_files:
                mtime = float(mfile[6:8]) + ((cube[0].dim_coords[0].points)/24.0)
                nans = ax.get_ylim()
                ax.fill_between(mtime, nans[0], nans[-1], facecolor = 'lightgrey', zorder = 3)

    ### global plot properties
    plt.subplot(5,3,12)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(5,3,13)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(5,3,14)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201808_oden_metum_1Db.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201808_oden_metum_1Db.png'
    if month_flag == 9:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201809_oden_metum_1Db.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201809_oden_metum_1Db.png'
    if month_flag == -1:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum_1Db.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum_1Db.png'
    plt.savefig(fileout, dpi=300)
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
        ifs_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/papbpc_combined/'
    out_dir2 = 'OUT2/'

    ### TESTING/domain_tests/umnsaa_pa000
    ### 4_u-bg610_RA2M_CON/OUT_R1/papbpc_combined/
    ### 5_u-bl616_RA2M_CASIM/OUT_R0/

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

    doy = np.arange(225,259)        ## set DOY for full moccha figures

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    if combine == 0:
        name = '20180813_oden_'
        filename_um = um_root_dir + out_dir1 + name + 'metum.nc'
        filename_ifs = ifs_root_dir + out_dir2 + name+ 'ecmwf.nc'
        print filename_um
        print filename_ifs
        print ''

        #### LOAD CUBE
        print 'Loading UM diagnostics:'
        cube_um = iris.load(filename_um)#, global_con, callback)
        print '...'
        print 'Loading IFS diagnostics:'
        cube_ifs = iris.load(filename_ifs)#, global_con, callback)
        # -------------------------------------------------------------

        print cube_um
        print ''
        print cube_ifs
        print ''

    else:
        for i in range(0,len(names)):
            filename_um = um_root_dir + out_dir1 + names[i] + 'metum.nc'
            filename_ifs = ifs_root_dir + out_dir2 + names[i] + 'ecmwf.nc'
            print filename_um
            print filename_ifs
            print ''

            #### LOAD CUBE
            print 'Loading UM diagnostics:'
            cube_um = iris.load(filename_um)#, global_con, callback)
            print '...'
            print 'Loading IFS diagnostics:'
            cube_ifs = iris.load(filename_ifs)#, global_con, callback)
            print '...'
            # -------------------------------------------------------------
            # print 'i = ' + str(i)
            print ''

            if i == 0:
                ## ------------------
                #### UM
                ## ------------------
                data_um = {}
                data1d_um = {}
                if month_flag == -1:
                    time_um = doy[i] + ((cube_um[0].dim_coords[0].points)/24.0)
                else:
                    time_um = float(filename_um[-16:-14]) + ((cube_um[0].dim_coords[0].points)/24.0)
                for j in range(0,len(cube_um)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um[j].data.shape) == 24:  # 1d timeseries only
                        data1d_um[cube_um[j].var_name] = cube_um[j].data
                    else:                                   # 2d column data
                        data_um[cube_um[j].var_name] = cube_um[j].data
                ## ------------------
                #### IFS
                ## ------------------
                data_ifs = {}
                data1d_ifs = {}
                if month_flag == -1:
                    time_ifs = doy[i] + ((cube_ifs[0].dim_coords[0].points)/24.0)
                else:
                    time_ifs = float(filename_ifs[-16:-14]) + ((cube_ifs[0].dim_coords[0].points)/24.0)
                for j in range(0,len(cube_ifs)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_ifs[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_ifs[j].data.shape) == 24:  # 1d timeseries only
                        data1d_ifs[cube_ifs[j].var_name] = cube_ifs[j].data
                    else:                                   # 2d column data
                        data_ifs[cube_ifs[j].var_name] = cube_ifs[j].data
            else:
                if month_flag == -1:
                    time_um = np.append(time_um, doy[i] + ((cube_um[0].dim_coords[0].points)/24.0))
                    time_ifs = np.append(time_ifs, doy[i] + ((cube_ifs[0].dim_coords[0].points)/24.0))
                else:
                    time_um = np.append(time_um,float(filename_um[-16:-14]) + ((cube_um[0].dim_coords[0].points)/24.0))
                    time_ifs = np.append(time_ifs,float(filename_ifs[-16:-14]) + ((cube_ifs[0].dim_coords[0].points)/24.0))
                ## ------------------
                #### UM
                ## ------------------
                for j in range(0,len(cube_um)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_um[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_um[j].data.shape) == 24:
                        data1d_um[cube_um[j].var_name] = np.append(data1d_um[cube_um[j].var_name].data,cube_um[j].data)
                    else:
                        data_um[cube_um[j].var_name] = np.append(data_um[cube_um[j].var_name].data,cube_um[j].data,0)
                ## ------------------
                #### IFS
                ## ------------------
                for j in range(0,len(cube_ifs)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(cube_ifs[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(cube_ifs[j].data.shape) == 24:
                        data1d_ifs[cube_ifs[j].var_name] = np.append(data1d_ifs[cube_ifs[j].var_name].data,cube_ifs[j].data)
                    else:
                        data_ifs[cube_ifs[j].var_name] = np.append(data_ifs[cube_ifs[j].var_name].data,cube_ifs[j].data,0)

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
        figure = plot_line_TSa(time_um, time_ifs, data1d_um, data1d_ifs, cube_um, cube_ifs, month_flag, missing_files, out_dir1)
                    ### doesn't matter which cube, just needed for dim_coords + cube structure

        # figure = plot_line_TSb(timem, data1d, cube, month_flag, missing_files, out_dir)
                    ### doesn't matter which cube, just needed for dim_coords + cube structure

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

    # UM -> IFS comparisons:
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
