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
# import cartopy.crs as ccrs
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

def pullLatLon(filename):

    from netCDF4 import Dataset

    print '*****'
    print 'Extracting lat/lon from ECMWF netCDF file'
    print ''

    nc = Dataset(filename,'r')

    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]

    print 'ECMWF file at: (' + str(lon) + ', ' + str(lat) + ')'

    nc.close()

    return lat, lon

def checkLatLon(ship_data, filenames):

    print ''
    print 'Finding lat/lon of ship track'
    print '...'

    #################################################################
    ## find ship track coordinates
    #################################################################
    ship_index = trackShip(ship_data)

    print 'findLatLon testing:'
    print 'Ship (lon,lat): ' + str(ship_data.values[ship_index,7][0]) + ', ' + str(ship_data.values[ship_index,6][0])

    # lat, lon = pullLatLon(filenames[0])

    # ship_index = np.where(np.logical_and(np.greater_equal(lat[:],ship_data.values[drift_index,7][0]), np.less_equal(lat[:],ship_data.values[drift_index,7][1])))
    # print 'Ship index test'
    # print ship_index
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

def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==12,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==9),data.values[:,3]==1))
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

def plot_basemap(ship_data, lats, lons):

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
    dim = 300000

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

    ### MAP ONTO PROJECTION
    x, y = m(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, color = 'darkorange', linewidth = 2, label = 'Ship track')

    # lat, lon = np.meshgrid(lats, lons)
    x_ecmwf, y_ecmwf = m(lons, lats)
    # Plot grid box centres as scatter plot
    plt.scatter(x_ecmwf, y_ecmwf, 400,
            color = 'white', marker = 's',
            edgecolor = 'blue', linewidth = 2,
            label = 'ECMWF')

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
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
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
    date = '20180812'
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
    # Pull gridded ship track from data
    # -------------------------------------------------------------
    lats = np.zeros([38])        ## 'empty' list of 38 elements. can assign index without list.append
    lons = np.zeros([38])
    for i in range(0,38):
        lats[i], lons[i] = pullLatLon(filenames[i])

    print 'Lats = ' + str(lats)
    print 'Lons = ' + str(lons)

    # -------------------------------------------------------------
    # Pull gridded ship track from cube
    # -------------------------------------------------------------

    #### LOAD CUBE
    # if con_flag == 0: fcube, outfile = pullTrack(cube, grid_filename, var_con)
    # if con_flag == 1: fcube, outfile = pullTrack(cube, grid_filename, global_con)
    # ## Update netCDF comments
    # out = appendNetCDF(outfile)
    # # final_outfile = out_dir + grid_filename[9:17] + '_oden_metum.nc'
    # # os.rename(outfile, final_outfile)

    # print outfile

    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    ### select hour to plot
    # hour = 0
    map = plot_basemap(ship_data, lats, lons)

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
