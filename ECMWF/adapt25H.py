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
import pandas as pd
# import diags_MOCCHA as diags
# import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os

def readfile(filename):

    import pandas as pd

    # print '******'
    print ('')
    print ('Reading .txt file with pandas')
    print ('')

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

    # print '******'
    # print ''
    # # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    # print 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    # print ''
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    # print 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')'
    # print 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')'
    # print 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6]))
    # print 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7]))
    # print ''

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

    # print '******'
    # print ''
    # # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    # print 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')'
    # print 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')'
    # # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    # print 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4])
    # print ''

    return trackShip_index

def combineNC(nc1, nc2, filename1, filename2, date):

    '''
    Load in two netCDF files at a time, then join to make nc1 25h long
    for compatibility with Cloudnet
    '''

    ###################################
    ## Set fluxes to distinguish from model_level diags
    ###################################
    fluxes = ['flx_net_sw','flx_net_lw','flx_ls_snow','flx_ls_rain','flx_turb_mom_v',
    'flx_turb_mom_u','flx_conv_snow','flx_conv_rain','flx_height','flx_turb_moist','flx_down_sens_heat']

    qfields = ['qi','ql']

    #################################################################
    ## CREATE NEW NETCDF
    #################################################################
    nc = Dataset(filename1[-22:], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc.file_format)
    print ('')

    ###################################
    #### load in a cube to define dimensions
    ###################################
    cube0 = iris.load('DATA/' + date + '_moccha_ecmwf_001.nc')

    ###################################
    ## Data dimensions
    ###################################
    timem = nc.createDimension('time', 25)
    level = nc.createDimension('model_level_number', np.size(cube0[9].dim_coords[1].points))
    flevel = nc.createDimension('model_flux_level', np.size(cube0[0].dim_coords[1].points))
    freq = nc.createDimension('frequency', np.size(cube0[10].dim_coords[0].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = nc.createVariable('time', np.float64, ('time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since ' + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' 00:00:00 +00:00.'
    timem.units = 'hours'
    timem.long_name = 'hours_UTC'
    timem.standard_name = 'time'
    timem[:] = cube0[0].dim_coords[0].points

    #### model level
    level = nc.createVariable('model_level_number', np.float64, ('model_level_number',), fill_value='-9999')
    level.scale_factor = float(1)
    level.add_offset = float(0)
    level.comment = ''
    level.units = '1'
    level.long_name = 'model_level'
    level.standard_name = 'model_level_number'
    level.positive = 'down'
    level[:] = cube0[9].dim_coords[1].points

    #### flux model level
    flevel = nc.createVariable('model_flux_level', np.float64, ('model_flux_level',), fill_value='-9999')
    flevel.scale_factor = float(1)
    flevel.add_offset = float(0)
    flevel.comment = ''
    flevel.units = '1'
    flevel.long_name = 'model_flux_level'
    flevel.positive = 'down'
    flevel[:] = cube0[0].dim_coords[1].points

    #### frequency
    freq = nc.createVariable('frequency', np.float64, ('frequency',), fill_value='-9999')
    freq.scale_factor = float(1)
    freq.add_offset = float(0)
    freq.comment = ''
    freq.units = 'GHz'
    freq.long_name = 'microwave_frequency'
    freq.missing_value = -999.0
    freq[:] = cube0[10].dim_coords[0].points

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    for diag in nc1.variables:
    # diag = 'sfc_pressure'
    # if diag == 'sfc_pressure':
        print ('Writing ' + diag)
        print ('')
        ### 1Dimension
        if np.size(np.shape(nc1.variables[diag])) == 1:
            if diag == 'time':
                print ('Diagnostic is forecast_time which is already defined... skipping.')
                continue
            if diag == 'model_level_number':
                print ('Diagnostic is model_level_number which is already defined... skipping.')
                continue
            if diag == 'model_flux_level':
                print ('Diagnostic is model_flux_level which is already defined... skipping.')
                continue
            if diag == 'frequency':
                print ('Diagnostic is frequency which is already defined... skipping.')
                continue
            dat = nc.createVariable(diag, np.float64, ('time',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            if 'units' in nc1.variables[diag].ncattrs(): dat.units = nc1.variables[diag].units
            if 'STASH' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].STASH
            if 'um_stash_source' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].um_stash_source
            if 'standard_name' in nc1.variables[diag].ncattrs(): dat.standard_name = nc1.variables[diag].standard_name
            if 'long_name' in nc1.variables[diag].ncattrs(): dat.long_name = nc1.variables[diag].long_name
            dat[0:24] = nc1.variables[diag][0:]
            dat[24] = nc2.variables[diag][0]

        ### 2Dimensions:         time: 24; model_level_number: 137 / time: 24; model_flux_level: 138
        elif np.size(np.shape(nc1.variables[diag])) == 2:
            if diag in fluxes:
                print ('Diagnostic is on flux levels.')
                dat = nc.createVariable(diag, np.float64, ('time','model_flux_level',), fill_value='-9999')
            else:
                print ('Diagnostic is on model levels.')
                dat = nc.createVariable(diag, np.float64, ('time','model_level_number',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            if 'units' in nc1.variables[diag].ncattrs(): dat.units = nc1.variables[diag].units
            if 'STASH' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].STASH
            if 'um_stash_source' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].um_stash_source
            if 'standard_name' in nc1.variables[diag].ncattrs(): dat.standard_name = nc1.variables[diag].standard_name
            if 'long_name' in nc1.variables[diag].ncattrs(): dat.long_name = nc1.variables[diag].long_name
            dat[0:24,:] = nc1.variables[diag][0:,:]
            dat[24,:] = nc2.variables[diag][0,:]

        ### 3Dimensions:         microwave_frequency: 2; time: 24; model_level_number: 137
        elif np.size(np.shape(nc1.variables[diag])) == 3:
                dat = nc.createVariable(diag, np.float64, ('frequency','time','model_level_number',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                if 'units' in nc1.variables[diag].ncattrs(): dat.units = nc1.variables[diag].units
                if 'STASH' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].STASH
                if 'um_stash_source' in nc1.variables[diag].ncattrs(): dat.um_stash_source = nc1.variables[diag].um_stash_source
                if 'standard_name' in nc1.variables[diag].ncattrs(): dat.standard_name = nc1.variables[diag].standard_name
                if 'long_name' in nc1.variables[diag].ncattrs(): dat.long_name = nc1.variables[diag].long_name
                dat[:,0:24,:] = nc1.variables[diag][:,0:,:]
                dat[:,24,:] = nc2.variables[diag][:,0,:]

        ### 0Dimensions
        else:
            if diag == 'horizontal_resolution':
                print ('Diagnostic is horizontal_resolution which needs to be defined separately...')
                dat = nc.createVariable('horizontal_resolution', np.float32, fill_value='-9999')
                dat.comment = 'Horizontal grid size.'
                dat.units = 'km'
                dat[:] = nc1.variables['horizontal_resolution'][:]
                continue
            elif diag == 'latitude':
                print ('Diagnostic is latitude which needs to be defined separately...')
                dat = nc.createVariable('latitude', np.float32, fill_value='-9999')
                dat.standard_name = 'latitude'
                dat.units = 'degrees_N'
                dat[:] = nc1.variables['latitude'][:]
                continue
            elif diag == 'longitude':
                print ('Diagnostic is longitude which needs to be defined separately...')
                dat = nc.createVariable('longitude', np.float32, fill_value='-9999')
                dat.standard_name = 'longitude'
                dat.units = 'degrees_E'
                dat[:] = nc1.variables['longitude'][:]
                continue

    ###################################
    ## Add Global Attributes
    ###################################
    if 'title' in nc1.ncattrs():
        nc.title = nc1.title
    elif 'title' in nc2.ncattrs():
        nc.title = nc2.title
    if 'description' in nc1.ncattrs():
        nc.description = nc1.description
    elif 'description' in nc2.ncattrs():
        nc.description = nc2.description
    if 'history' in nc1.ncattrs():
        nc.history = nc1.history
    elif 'history' in nc2.ncattrs():
        nc.history = nc2.history
    if 'source' in nc1.ncattrs():
        nc.source = nc1.source
    elif 'source' in nc2.ncattrs():
        nc.source = nc2.source
    if 'project' in nc1.ncattrs():
        nc.project = nc1.project
    elif 'project' in nc2.ncattrs():
        nc.project = nc2.project
    if 'institution' in nc1.ncattrs():
        nc.institution = nc1.institution
    elif 'institution' in nc2.ncattrs():
        nc.institution = nc2.institution
    if 'initialization_time' in nc1.ncattrs():
        nc.initialization_time = nc1.initialization_time
    if date == '20180904':
        nc.initialization_time = '2018-09-04 00:00:00 +00:00'

    nc.close()

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### JASMIN
    ### LAPTOP
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ECMWF/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '/home/gillian/MOCCHA/ECMWF/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/'
        ship_filename = '/home/gillian/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/data/ecmwf_ewan/moccha/ecmwf-all/2018/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'

    ### CHOSEN RUN
    out_dir = 'OUT2/'
    out_dir3 = 'MET_DATA/'

    ### TESTING/domain_tests/umnsaa_pa000
    ### 4_u-bg610_RA2M_CON/OUT_R1/papbpc_combined/
    ### 5_u-bl661_RA1M_CASIM/OUT/
    ### 6_u-bm410_RA1M_CASIM/OUT/

    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    # print 'Loading observations:'
    # filename_obs = obs_root_dir + out_dir3 + 'MetData_Gillian_wTemp1p5m.nc'
    # cube_obs = iris.load(filename_obs)#, global_con, callback)
    # print '...'

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin cube read in at ' + time.strftime("%c"))
    print (' ')

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    # tempnames = ['umnsaa_pa012_r0.nc','umnsaa_pb012_r0.nc','umnsaa_pc011_r0.nc','umnsaa_pd011_r0.nc','20180812_oden_ecmwf.nc']
    Aug_names = ['20180813_oden_ecmwf.nc','20180814_oden_ecmwf.nc','20180815_oden_ecmwf.nc','20180816_oden_ecmwf.nc',
            '20180817_oden_ecmwf.nc','20180818_oden_ecmwf.nc','20180819_oden_ecmwf.nc','20180820_oden_ecmwf.nc',
            '20180821_oden_ecmwf.nc','20180822_oden_ecmwf.nc','20180823_oden_ecmwf.nc','20180824_oden_ecmwf.nc',
            '20180825_oden_ecmwf.nc','20180826_oden_ecmwf.nc','20180827_oden_ecmwf.nc','20180828_oden_ecmwf.nc',
            '20180829_oden_ecmwf.nc','20180830_oden_ecmwf.nc','20180831_oden_ecmwf.nc']

    Sep_names = ['20180901_oden_ecmwf.nc','20180902_oden_ecmwf.nc','20180903_oden_ecmwf.nc','20180904_oden_ecmwf.nc',
            '20180905_oden_ecmwf.nc','20180906_oden_ecmwf.nc','20180907_oden_ecmwf.nc','20180908_oden_ecmwf.nc',
            '20180909_oden_ecmwf.nc','20180910_oden_ecmwf.nc','20180911_oden_ecmwf.nc','20180912_oden_ecmwf.nc',
            '20180913_oden_ecmwf.nc','20180914_oden_ecmwf.nc']

    moccha_names = ['20180813_oden_ecmwf.nc','20180814_oden_ecmwf.nc','20180815_oden_ecmwf.nc','20180816_oden_ecmwf.nc',
            '20180817_oden_ecmwf.nc','20180818_oden_ecmwf.nc','20180819_oden_ecmwf.nc','20180820_oden_ecmwf.nc',
            '20180821_oden_ecmwf.nc','20180822_oden_ecmwf.nc','20180823_oden_ecmwf.nc','20180824_oden_ecmwf.nc',
            '20180825_oden_ecmwf.nc','20180826_oden_ecmwf.nc','20180827_oden_ecmwf.nc','20180828_oden_ecmwf.nc',
            '20180829_oden_ecmwf.nc','20180830_oden_ecmwf.nc','20180831_oden_ecmwf.nc','20180901_oden_ecmwf.nc',
            '20180902_oden_ecmwf.nc','20180903_oden_ecmwf.nc','20180904_oden_ecmwf.nc','20180905_oden_ecmwf.nc',
            '20180906_oden_ecmwf.nc','20180907_oden_ecmwf.nc','20180908_oden_ecmwf.nc','20180909_oden_ecmwf.nc',
            '20180910_oden_ecmwf.nc','20180911_oden_ecmwf.nc','20180912_oden_ecmwf.nc','20180913_oden_ecmwf.nc',
            '20180914_oden_ecmwf.nc']

    Aug_missing_files = ['20180812_oden_ecmwf.nc']

    Sep_missing_files = []

    moccha_missing_files = []

    doy = np.arange(225,258)        ## set DOY for full moccha figures

    # date = '20180813'

    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    # i = 0
    for i in range(0,len(moccha_names) - 1):
        filename1 = root_dir + out_dir + names[i]
        filename2 = root_dir + out_dir + names[i+1]
        print (filename1)
        print (filename2)
        print ('')
        date = names[i][0:8]
        print ('Date = ' + date)

        #### -------------------------------------------------------------
        #### LOAD NETCDF FILES
        #### -------------------------------------------------------------
        # cube1 = iris.load(filename1)
        nc1 = Dataset(filename1,'r')
        print (nc1)
        print ('')

        # cube2 = iris.load(filename1)
        nc2 = Dataset(filename2,'r')
        print (nc2)
        print ('')

        #### -------------------------------------------------------------
        #### COMBINE NETCDF FILES
        #### -------------------------------------------------------------
        out = combineNC(nc1, nc2, filename1, filename2, date)

        #### -------------------------------------------------------------
        #### CLOSE ORIGINAL NETCDF FILES
        #### -------------------------------------------------------------
        nc1.close()
        nc2.close()

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')

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
