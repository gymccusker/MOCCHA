###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE,
###         PULL SHIP TRACK, AND OUTPUT AS NETCDF FOR CLOUDNET
###


import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
# import diags_MOCCHA as diags
# import diags_varnames as varnames
import cartopy.crs as ccrs
# import iris
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

    return lat, lon

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
    print 'Writing ' + date + ' grid to file:'
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

def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==30,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==31,data.values[:,1]==8),data.values[:,3]==1))
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

    # import iris.plot as iplt
    # import iris.quickplot as qplt
    # import iris.analysis.cartography
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
    diag = 1
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
    # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
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
    # if np.size(cube[diag].data.shape) == 4:
    #     iplt.pcolormesh(cube[diag][hour,0,:,:])
    # elif np.size(cube[diag].data.shape) == 3:
    #     iplt.pcolormesh(cube[diag][hour,:,:])
    #     # iplt.pcolormesh(cube[hour,471:495,240:264])
    # elif np.size(cube[diag].data.shape) == 2:
    #     iplt.pcolormesh(cube[diag][:,:])
    # plt.title(cube[diag].standard_name + ', ' + str(cube[diag].units))
    # plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[diag][hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[diag][hour,386:495,211:305])          ### misc
    qplt.outline(cube[diag][hour,:,:])

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

    #################################################################
    ## read in and plot gridded ship track
    #################################################################
    tim, ilat, ilon = readGriddedTrack(grid_filename)

    ### Plot tracks as line plot
    for i in range(0, len(ilon)-1):
        iplt.scatter(cube[diag].dim_coords[2][int(ilon[i] + xoffset)], cube[diag].dim_coords[1][int(ilat[i] + yoffset)],color='black')


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

def pullTrack(cube, grid_filename, con):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt

    print '******'
    print ''

    #################################################################
    ## load gridded ship track
    #################################################################
    # print '******'
    print ''
    print 'Pulling gridded track from netCDF:'
    print ''

    tim, ilat, ilon = readGriddedTrack(grid_filename)

    #################################################################
    ## fix time index
    #################################################################

    if np.size(cube)>1:
        print ''
        print 'More than one variable constraint. Proceeding...'
        print ''

        cubetime = np.round(cube[0].coord('forecast_period').points - 12.0)      ### forecast period (ignore first 12h)
        print ''
        print 'Cube times relative to forecast start:', cubetime[:-1]
        print ''

        #################################################################
        ## CREATE EMPTY CUBE
        #################################################################
        ncube = Cube(np.zeros([np.size(cube),70,len(cubetime)-1]))

        #################################################################
        ## POPULATE NP ARRAY WITH DATA
        #################################################################
        ### populate 0th dimension with time field
        # data[:,0] = cubetime[:,:-1]

        for k in range(0,np.size(cube)):            ### loop over number of variables
            print ''
            print 'k = ', k, ###', so processing', con[k]   # doesn't work with global_con
            print ''
            #################################################################
            ## PROBE VARIABLE
            #################################################################
            ### do we want to average exluding zeros?
            stash_flag, stash = excludeZeros(cube[k])

            ### do we need to re-grid?  -- DOESN'T WORK LIKE WRF, GRID NOT SPACED SAME WAY
            # cube[k], wind_stash = checkWind(cube[k])

            #################################################################
            ## CHECK DIMENSIONS
            #################################################################
            if np.logical_and(np.size(cube[k].data,1) >= 69, np.size(cube[k].data,1) < 71):
                print 'Variable is 4D:'
                print ''
                #### create empty arrays to be filled
                data = np.zeros([len(cube[k].coord('model_level_number').points),len(cubetime)-1])
                ### make empty cube
                dim_flag = 1        ### for next loops
                print 'data.shape = ', str(data.shape)
                print ''
            else:
                print 'Variable is 3D:'
                print ''
                #### create empty arrays to be filled
                data = np.zeros([len(cubetime)-1])
                dim_flag = 0       ### for next loops
                print 'data.shape = ', str(data.shape)
                print ''

            #################################################################
            ## LOOP OVER TIME INDEX, DECOMPOSE ONTO 24H TIMESERIES
            #################################################################
            for j in range(0,len(cubetime)-1):              ### loop over time
                if j < len(cubetime[:-1]):
                    itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
                else:
                    ### end point (23h)
                    itime = np.where(tim >= cubetime[-1])
                print ''
                print 'For ', str(j), 'h, itime = ', itime
                if dim_flag == 1: dat = np.zeros([len(cube[k].coord('model_level_number').points),len(itime[0])])
                if dim_flag == 0: dat = np.zeros([len(itime[0])])
                for i in range(0, len(itime[0])):                   ### loop over time gridded by ship track
                    if np.size(itime) > 1:
                        # print 'Processing i = ', str(itime[0][i])
                        # print '...'
                        if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                        if dim_flag == 0: temp = cube[k][j,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                    else:
                        # print 'Processing i = ', str(itime[i])
                        # print '...'
                        if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                        if dim_flag == 0: temp = cube[k][j,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                    if dim_flag == 1: dat[:,i] = np.squeeze(temp.data)
                    if dim_flag == 0: dat[i] = np.squeeze(temp.data)
                    if np.size(itime) > 1:
                        if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                        if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                        if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                        # print 'averaging over itime ...'
                        # print ''
                    else:
                        if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                        if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                        # print 'no averaging, itime = 1 ...'
                        print ''
                # print data
        # print 'data.shape = ', data.shape

        #################################################################
        ## FIGURES TO TEST OUTPUT
        #################################################################
        ### timeseries of lowest model level
        # plt.figure(figsize=(7,5))
        # plt.plot(cubetime[:-1],data[0:10,:])
        # plt.show()

        ### vertical profile of 1st timestep
        # plt.figure(figsize=(7,5))
        # plt.plot(data[:,0],cube.coord('model_level_number').points)
        # plt.show()

        ### pcolormesh of timeseries
        # plt.figure(figsize=(7,5))
        # plt.pcolormesh(cubetime[:-1], cube.coord('model_level_number').points, data)
        # plt.colorbar()
        # plt.show()

        #################################################################
        ## CREATE CUBE
        #################################################################
        ### ECMWF FIELD NAMES
        # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
        #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
        #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};
            varname = varnames.findfieldName(stash)
            print 'standard_name = ', cube[k].standard_name
            print 'long name = ', cube[k].long_name
            print 'varname = ', varname
            print ''

            ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
            if dim_flag == 1:         ### 4D VARIABLE
                model_height = DimCoord(cube[k].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
                ncube = Cube(np.transpose(data),
                        dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
                        standard_name = cube[k].standard_name,
                        long_name = cube[k].long_name,
                        units = cube[k].units,
                        var_name = varname,
                        attributes = cube[k].attributes,
                        aux_coords_and_dims = None,
                        )
            elif dim_flag == 0:         ### 3D VARIABLE
                ncube = Cube(np.transpose(data),
                        dim_coords_and_dims=[(ntime, 0)],
                        standard_name = cube[k].standard_name,
                        long_name = cube[k].long_name,
                        units = cube[k].units,
                        var_name = varname,
                        attributes = cube[k].attributes,
                        aux_coords_and_dims = None,
                        )
            # ncube.attributes = cube[k].attributes
            # iris.save(ncube, pp_outfile, append=True)
            if k == 0:
                print 'Assigning fcube'
                print ''
                fcube = [ncube]
            else:
                print 'Appending to fcube'
                print ''
                fcube.append(ncube)

        # print fcube

    else:
        print ''
        print 'Only one variable constraint. Proceeding...'
        print ''

        cubetime = np.round(cube.coord('forecast_period').points - 12.0)      ### forecast period (ignore first 12h)
        print ''
        print 'Cube times relative to forecast start:', cubetime[:-1]
        print ''

        #################################################################
        ## CREATE EMPTY CUBE
        #################################################################
        ncube = Cube(np.zeros([len(cube.coord('model_level_number').points),len(cubetime)-1]))

        #################################################################
        ## PROBE VARIABLE
        #################################################################
        ### do we want to average exluding zeros?
        stash_flag, stash = excludeZeros(cube)

        #################################################################
        ## FIND ARRAY SIZE AND CREATE EMPTY NP ARRAY
        #################################################################
        if np.logical_and(np.size(cube.data,1) >= 69, np.size(cube.data,1) < 71):
            print 'Variable is 4D:'
            print ''
            #### create empty arrays to be filled
            data = np.zeros([len(cube.coord('model_level_number').points),len(cubetime)-1])
            dim_flag = 1        ### for next loops
            print 'data.shape = ', str(data.shape)
            print ''
        else:
            print 'Variable is 3D:'
            print ''
            #### create empty arrays to be filled
            data = np.zeros([len(cubetime)-1])
            dim_flag = 0       ### for next loops
            print 'data.shape = ', str(data.shape)
            print ''

        #################################################################
        ## POPULATE NP ARRAY WITH DATA
        #################################################################
        ### populate 0th dimension with time field
        # data[:,0] = cubetime[:,:-1]

        for j in range(0,len(cubetime)-1):
            if j < len(cubetime[:-1]):
                itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
            else:
                ### end point (23h)
                itime = np.where(tim >= cubetime[-1])
            print 'For ', str(j), 'h, itime = ', itime
            if dim_flag == 1: dat = np.zeros([len(cube.coord('model_level_number').points),len(itime[0])])
            if dim_flag == 0: dat = np.zeros([len(itime[0])])
            for i in range(0, len(itime[0])):
                if np.size(itime) > 1:
                    # print 'Processing i = ', str(itime[0][i])
                    if dim_flag == 1: temp = cube[j,:,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                    if dim_flag == 0: temp = cube[j,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                else:
                    # print 'Processing i = ', str(itime[i])
                    if dim_flag == 1: temp = cube[j,:,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                    if dim_flag == 0: temp = cube[j,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                if dim_flag == 1: dat[:,i] = temp.data
                if dim_flag == 0: dat[i] = temp.data
                if np.size(itime) > 1:
                    if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                    if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                    if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                    # print 'averaging over itime...'
                else:
                    if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                    if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                    # print 'no averaging, itime = 1...'
        # print data
        # print 'data.shape = ', data.shape

        #################################################################
        ## FIGURES TO TEST OUTPUT
        #################################################################
        ### timeseries of lowest model level
        # plt.figure(figsize=(7,5))
        # plt.plot(cubetime[:-1],data[0:10,:])
        # plt.show()

        ### vertical profile of 1st timestep
        # plt.figure(figsize=(7,5))
        # plt.plot(data[:,0],cube.coord('model_level_number').points)
        # plt.show()

        ### pcolormesh of timeseries
        # plt.figure(figsize=(7,5))
        # plt.pcolormesh(cubetime[:-1], cube.coord('model_level_number').points, data)
        # plt.colorbar()
        # plt.show()

        #################################################################
        ## CREATE CUBE
        #################################################################
        ### ECMWF FIELD NAMES
        # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
        #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
        #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};

        varname = varnames.findfieldName(stash)
        print 'standard_name = ', cube.standard_name
        print 'long name = ', cube.long_name
        print 'varname = ', varname
        print ''

        ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
        if dim_flag == 1:             ### 4D VARIABLE
            model_height = DimCoord(cube.aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
            ncube = Cube(np.transpose(data),
                    dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
                    standard_name = cube.standard_name,
                    long_name = cube.long_name,
                    units = cube.units,
                    var_name = varname,
                    )
        elif dim_flag == 0:             ### 3D VARIABLE
            ncube = Cube(np.transpose(data),
                    dim_coords_and_dims=[(ntime, 0)],
                    standard_name = cube.standard_name,
                    long_name = cube.long_name,
                    units = cube.units,
                    var_name = varname,
                    )
        ncube.attributes = cube.attributes
        ### for consistency with multi-diag option
        fcube = ncube

    #################################################################
    ## CREATE NETCDF
    #################################################################

    #################################################################
    ## define output filename
    #################################################################
    print '******'
    print 'Define outfile:'
    # pp_outfile = out_dir + grid_filename[9:17] + '_oden_metum.pp'
    # nc_outfile = out_dir + grid_filename[9:17] + '_oden_metum.nc'
    # pp_outfile = grid_filename[9:17] + '_oden_metum.pp'
    nc_outfile = grid_filename[9:17] + '_oden_metum.nc'
    print 'Outfile = ', nc_outfile

    ### save cube to netcdf file
    print ''
    print 'Writing fcube to NetCDF file:'
    print ''
    # iris.save(fcube, nc_outfile)
    print fcube

    return fcube, nc_outfile

def appendNetCDF(outfile):

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
    dataset =  Dataset(outfile, 'a', format ='NETCDF4_CLASSIC')
    # infile = net.Dataset("2015%s%s-160000_0.nc" % (month,day), "a")
    print ''
    print dataset.file_format
    print ''

    ###################################
    ## Global Attributes
    ###################################
    dataset.title = 'ECMWF Model single-site output during MOCCHA'
    dataset.description = 'Hourly data taken from grid box closest to ship location. Where the ship covers more than one grid box within an hour period, data are averaged from all grid boxes crossed.'
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> from NetCDF generated from original data by Ewan O''Connor <ewan.oconnor@fmi.fi> using cnmodel2nc on cloudnet.fmi.fi.'
    dataset.source = 'ECMWF Integrated Forecast System (IFS)'
    # dataset.references = 'N/A'
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    dataset.comment = micro + wind
    dataset.institution = 'University of Leeds/FMI.'
    # dataset.initialization_time = outfile[0:4] + '-' + outfile[4:6] + '-' + outfile[6:8]) + ' 00:00:00 UTC.'
    dataset.initialization_time = outfile[0:4] + '-' + outfile[4:6] + '-' + str(int(outfile[6:8]) - 1) + ' 12:00:00 UTC.'

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'DESKTOP'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
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
    ## 4_OPER/20180830T0000Z_TRIAL/

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    grid_dirname = 'AUX_DATA/'
    date = '20180812'
    grid_filename = grid_dirname + date + '_ShipTrack_GRIDDED.csv'

    print '******'
    print ''
    print 'Identifying .nc file: '
    print ''

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
    # var_con = 'specific_humidity'
    # cube = iris.load_cube(filename1, var_con)
    # global_con = ['atmosphere_downward_eastward_stress','atmosphere_downward_northward_stress']

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    names = ['umnsaa_pa012_r0.nc','umnsaa_pb012_r0.nc','umnsaa_pc011_r0.nc','umnsaa_pd011_r0.nc']
    filename1 = root_dir + out_dir + names[2]
    print filename1
    print ''

    #### LOAD CUBE
    if 'var_con' in locals():
        print 'Loading single diagnostic:'
        print var_con
        cube1 = iris.load_cube(filename1, var_con, callback)
        con_flag = 0            # constraint flag
    elif 'global_con' in locals():
        print 'Loading multiple diagnostics:'
        # cube = iris.load_cubes(filename1, global_con)
        cube = iris.load(filename1, global_con, callback)
        con_flag = 1            # constraint flag

        # -------------------------------------------------------------

    print cube
    print ''

    # -------------------------------------------------------------
    # Pull gridded ship track from cube
    # -------------------------------------------------------------

    #### LOAD CUBE
    if con_flag == 0: fcube, outfile = pullTrack(cube, grid_filename, var_con)
    if con_flag == 1: fcube, outfile = pullTrack(cube, grid_filename, global_con)
    ## Update netCDF comments
    out = appendNetCDF(outfile)
    # final_outfile = out_dir + grid_filename[9:17] + '_oden_metum.nc'
    # os.rename(outfile, final_outfile)

    # print outfile

    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    ### select hour to plot
    # hour = 0
    # map = plot_cartmap(ship_data, cube, hour, grid_filename)#, lon, lat)


    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

    #### DIAGNOSTICS TO CHOOSE FROM:

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
    # 0: large_scale_ice_water_path / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 1: large_scale_liquid_water_path / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 2: eastward_wind_at_10m / (m s-1)      (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 3: northward_wind_at_10m / (m s-1)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 4: air_temperature_at_1.5m / (K)       (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 5: specific_humidity_at_1.5m / (1)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 6: relative_humidity_at_1.5m / (%)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 7: dew_point_temperature_at_1.5m / (K) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 8: turbulent mixing height after boundary layer / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 9: height_of_decoupled_layer_base / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 10: height_of_stratocumulus_cloud_base / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 11: combined_boundary_layer_type / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 12: cloud_area_fraction_assuming_random_overlap / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 13: cloud_area_fraction_assuming_maximum_random_overlap / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 14: wet_bulb_freezing_level_altitude / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 15: total_column_q / (unknown)          (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 16: air_pressure_at_sea_level / (Pa)    (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 17: atmosphere_boundary_layer_thickness / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 18: high_type_cloud_area_fraction / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 19: low_type_cloud_area_fraction / (1)  (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 20: medium_type_cloud_area_fraction / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 21: stratiform_rainfall_flux / (kg m-2 s-1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 22: stratiform_snowfall_flux / (kg m-2 s-1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 23: surface_air_pressure / (Pa)         (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 24: surface_temperature / (K)           (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 25: surface_upward_latent_heat_flux / (W m-2) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 26: surface_upward_sensible_heat_flux / (W m-2) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 27: water_evaporation_amount / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)

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
