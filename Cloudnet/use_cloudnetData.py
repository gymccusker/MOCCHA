###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS nc
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


def plot_contour_TS(nc, filename): #, lon, lat):

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
    print nc[diag]
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

    time = nc[diag].dim_coords[0].points
    height = nc[diag].dim_coords[1].points

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

    # plt.plot(nc[diag].dim_coords[0].points,nc[diag][:,0].data)        # line plot
    # plt.contourf(nc[0].data)
    # plt.plot(nc[2][0,:].data,height);plt.show()
    #################################################################
    ## plot contour timeseries
    ################################################################
    plt.contourf(time,height,np.transpose(nc[diag].data))
    # plt.pcolormesh(time,height,np.transpose(nc[2].data))
    plt.title(nc[diag].standard_name + ', ' + str(nc[diag].units))
    plt.colorbar()
    ax.set_ylim([0, 3000])

    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

def plot_multicontour_TS(nc, filename, out_dir): #, lon, lat):

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
    for i in range(0,len(nc)):
        ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
        if np.sum(nc[i].data.shape) > 24:

                ###################################
                ## CHOOSE DIAGNOSTIC
                ###################################
                diag = i
                print ''
                print 'Diag is: '
                print nc[diag]

                ### define empty array for nc data
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

                time = nc[diag].dim_coords[0].points
                height = nc[diag].dim_coords[1].points

                ### if mass mixing ratio, *1e3 to change to g/kg
                # if nc[diag].var_name[0] == 'q':
                #     data = np.transpose(np.squeeze(nc[diag].data[:,ind]*1e3))
                # elif nc[diag].var_name == 'pressure':
                #     data = np.transpose(np.squeeze(nc[diag].data[:,ind]/1e2))
                # else:
                #     data = np.transpose(np.squeeze(nc[diag].data[:,ind]))

                #################################################################
                ## data corrections
                #################################################################
                ### set height limit to consider
                ind = np.where(height<5000)

                if nc[diag].var_name == 'temperature':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'
                elif nc[diag].var_name == 'qice':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]*1e3))
                    title = nc[diag].var_name + ' [g/kg]'
                elif nc[diag].var_name == 'qliq':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]*1e3))
                    title = nc[diag].var_name + ' [g/kg]'
                elif nc[diag].var_name == 'q':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]*1e3))
                    title = nc[diag].var_name + ' [g/kg]'
                elif nc[diag].var_name == 'pressure':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]/1e2))
                    title = nc[diag].var_name + ' [hPa]'
                elif nc[diag].var_name == 'uwind':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'
                elif nc[diag].var_name == 'wwind':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'
                elif nc[diag].var_name == 'radr_refl':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'
                elif nc[diag].var_name == 'cloud_fraction':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'
                elif nc[diag].var_name == 'vwind':
                    data = np.transpose(np.squeeze(nc[diag].data[:,ind]))
                    title = nc[diag].var_name + ' [' + str(nc[diag].units) + ']'

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
                    # plt.contourf(time,height,np.transpose(nc[diag].data))
                    if nc[diag].var_name == 'temperature':
                        plt.pcolormesh(time, height[ind], data, vmin = 250, vmax = np.nanmax(data))
                    elif nc[diag].var_name == 'uwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -20, vmax = 20)
                    elif nc[diag].var_name == 'vwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -20, vmax = 20)
                    elif nc[diag].var_name == 'wwind':
                        plt.pcolormesh(time, height[ind], data, vmin = -0.1, vmax = 0.1)
                    else:
                        plt.pcolormesh(time, height[ind], data, vmin = np.nanmin(data), vmax = np.nanmax(data))

                    #################################################################
                    ## set plot properties
                    #################################################################
                    ### colormaps:
                    if nc[diag].var_name == 'wwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc[diag].var_name == 'uwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc[diag].var_name == 'vwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc[diag].var_name[0] == 'q':
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

def plot_multicontour_multidate_TS(timem, data, nc, month_flag, missing_files, out_dir, doy): #, lon, lat):

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

        ### pcXXX - nc
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
        ###         way to index nc by string
        for j in range(0,len(nc)):
            if nc[j].var_name == 'qliq': height = nc[j].dim_coords[1].points

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
            # plt.contourf(time,height,np.transpose(nc[diag].data))
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
            if month_flag == -1: ax.set_xlim([doy[0], doy[-1]])

            print ''
            print 'Zero out any data from missing files:'
            print ''
            for mfile in missing_files:
                # mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
                # nans = np.zeros([len(height),len(mtime)])
                # # nans[nans == 0] = np.nan
                # plt.pcolormesh(mtime, height, nans)
                mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
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
        if out_dir[:18] == '7_u-bn068_RA2M_PC2':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum.png'
        if out_dir[:18] == '6_u-bm410_RA1M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum.png'
        if out_dir[:18] == '5_u-bl661_RA1M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum.png'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_multicontour_multidate_casim_TS(timem, data, nc, month_flag, missing_files, out_dir): #, lon, lat):

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

        ### pcXXX - nc
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
        ###         way to index nc by string
        for j in range(0,len(nc)):
            if nc[j].var_name == 'qliq': height = nc[j].dim_coords[1].points

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
            # plt.contourf(time,height,np.transpose(nc[diag].data))
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
                # mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
                # nans = np.zeros([len(height),len(mtime)])
                # # nans[nans == 0] = np.nan
                # plt.pcolormesh(mtime, height, nans)
                mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
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

def plot_line_TSa(timem, data, nc, month_flag, missing_files, out_dir): #, lon, lat):

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
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.9, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    print data.keys()
    # [u'LWP', u'sfc_temperature', u'IWP', u'sfc_pressure', u'rainfall_flux', u'snowfall_flux']

    l = -1
    jindex = 0

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

        # if str(data.keys()[diag]) == 'LWP':
        #     dat = data[data.keys()[diag]].data*1e3
        #     title = str(data.keys()[diag]) + ' [g/m2]'
        # elif str(data.keys()[diag]) == 'IWP':
        #     dat = data[data.keys()[diag]].data*1e3
        #     title = str(data.keys()[diag]) + ' [g/m2]'
        # elif str(data.keys()[diag]) == 'rainfall_flux':
        #     dat = data[data.keys()[diag]].data*3600
        #     title = str(data.keys()[diag]) + ' [mm/hr]'
        # elif str(data.keys()[diag]) == 'snowfall_flux':
        #     dat = data[data.keys()[diag]].data*3600
        #     title = str(data.keys()[diag]) + ' [mm/hr]'

        if str(data.keys()[diag]) == 'sfc_temperature':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [K]'
        elif str(data.keys()[diag]) == 'sfc_pressure':
            dat = data[data.keys()[diag]].data/1e2
            title = str(data.keys()[diag]) + ' [hPa]'
        elif str(data.keys()[diag]) == 'rh_1.5m':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [%]'
        elif str(data.keys()[diag]) == 'temp_1.5m':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [K]'
        elif str(data.keys()[diag]) == 'surface_net_SW_radiation':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [W/m2]'
        elif str(data.keys()[diag]) == 'surface_net_LW_radiation':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' [W/m2]'
        elif str(data.keys()[diag]) == 'fogfrac_1.5m':
            dat = data[data.keys()[diag]].data
            title = str(data.keys()[diag]) + ' []'
        elif str(data.keys()[diag]) == 'vis_1.5m':
            dat = data[data.keys()[diag]].data/1.0e3
            title = str(data.keys()[diag]) + ' [km]'

        # str(data.keys()[0][-4:])

        #################################################################
        ## create figure and axes instances
        #################################################################
        if len(dat) > 0:
            l = l + 1 ## increment index for positive data association
            plt.subplot(4,2,l+1)
            print 'l = ' + str(l)
            print title
            ax = plt.gca()

            #################################################################
            ## plot timeseries
            #################################################################

            plt.plot(timem, dat)
            plt.ylim([np.nanmin(dat),np.nanmax(dat)])
            plt.title(title)

            if month_flag == 8: ax.set_xlim([13.0, 31.0])
            if month_flag == 9: ax.set_xlim([1.0, 15.0])
            if month_flag == -1: ax.set_xlim([225.0, 258.0])

            print ''
            print 'Zero out any data from missing files:'
            print ''
            for mfile in missing_files:
                mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
                nans = ax.get_ylim()
                ax.fill_between(mtime, nans[0], nans[-1], facecolor = 'lightgrey', zorder = 3)

    ### global plot properties
    plt.subplot(4,2,7)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(4,2,8)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201808_oden_metum_1Da.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201808_oden_metum_1Da.png'
    if month_flag == 9:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:21] + '201809_oden_metum_1Da.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:19] + '201809_oden_metum_1Da.png'
    if month_flag == -1:
        if out_dir[:18] == '7_u-bn068_RA2M_PC2':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum_1Da.png'
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum_1Da.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum_1Da.png'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_TSb(timem, data, nc, month_flag, missing_files, out_dir): #, lon, lat):

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
                mtime = float(mfile[6:8]) + ((nc[0].dim_coords[0].points)/24.0)
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
        if out_dir[:18] == '7_u-bn068_RA2M_PC2':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum_1Db.png'
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = 'FIGS/' + out_dir[:20] + '_oden_metum_1Db.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = 'FIGS/' + out_dir[:18] + '_oden_metum_1Db.png'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_TEMP(time_um, data1d_um, nc_um, month_flag, missing_files, out_dir, nc_obs, doy): #, lon, lat):

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

    datenums_temp = nc_obs[0].dim_coords[0].points
    timestamps_temp = pd.to_datetime(datenums_temp-719529, unit='D')
    time_temp = timestamps_temp.dayofyear + (timestamps_temp.hour / 24.0) + (timestamps_temp.minute / 1440.0) + (timestamps_temp.second / 86400.0)

    #################################################################
    ## create figure and axes instances
    #################################################################

    ax = plt.gca()
    plt.plot(time_um, data1d_um['temp_1.5m'].data - 273.15, color = 'r', label = 'MetUM')
    plt.plot(time_temp,nc_obs[0].data - 273.15, color = 'black', label = 'Observations')
    plt.legend()
    plt.title('Temperature_at_1.5m [$^{o}C$]')
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
        plt.xlabel('Day of year')

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == 8:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/UM/' + out_dir[:21] + '201808_oden_metum_temp.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/UM/' + out_dir[:19] + '201808_oden_metum_temp.png'
    if month_flag == 9:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/UM/' + out_dir[:21] + '201809_oden_metum_temp.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/UM/' + out_dir[:19] + '201809_oden_metum_temp.png'
    if month_flag == -1:
        if out_dir[:18] == '5_u-bl616_RA2M_CAS':
            fileout = '../FIGS/UM/' + out_dir[:20] + '_oden_metum_temp.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/UM/' + out_dir[:18] + '_oden_metum_temp_degC.png'
    plt.savefig(fileout, dpi=600)
    plt.show()

def plot_line_RAD(time_um, data1d_um, nc_um, month_flag, missing_files, out_dir, nc_obs, doy): #, lon, lat):

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

    datenums_temp = nc_obs[0].dim_coords[0].points
    timestamps_temp = pd.to_datetime(datenums_temp-719529, unit='D')
    time_temp = timestamps_temp.dayofyear + (timestamps_temp.hour / 24.0) + (timestamps_temp.minute / 1440.0) + (timestamps_temp.second / 86400.0)

    datenums_radice = nc_obs[1].dim_coords[0].points
    timestamps_radice = pd.to_datetime(datenums_radice-719529, unit='D')
    time_radice = timestamps_radice.dayofyear + (timestamps_radice.hour / 24.0) + (timestamps_radice.minute / 1440.0) + (timestamps_radice.second / 86400.0)

    #################################################################
    ## create figure and axes instances
    #################################################################

    plt.subplot(211)
    ax = plt.gca()
    plt.plot(time_um, data1d_um['temp_1.5m'].data - 273.15, color = 'r', label = 'MetUM')
    plt.plot(time_temp,nc_obs[0].data - 273.15, color = 'black', label = 'Observations')
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
    data1d_um['surface_net_SW_radiation'].data[data1d_um['surface_net_SW_radiation'].data == 0] = np.nan
    plt.plot(time_um, data1d_um['surface_net_SW_radiation'].data, color = 'r', label = 'MetUM')
    plt.plot(time_radice,(nc_obs[7].data - nc_obs[8].data), color = 'black', label = 'Observations')
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
    # plt.plot(time_radice,(nc_obs[1].data - nc_obs[2].data), color = 'black', label = 'Observations')
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
        if out_dir[:20] == '5_u-bl661_RA1M_CASIM':
            fileout = '../FIGS/UM/' + out_dir[:20] + '_oden_metum_SW+LW.png'
        elif out_dir[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/UM/' + out_dir[:18] + '_oden_metum_tempdegC_SW.png'
    plt.savefig(fileout, dpi=400)
    plt.show()

def callback(nc, field, filename):
    '''
    rename nc diagnostics per list of wanted stash diags
    '''

    iStash = nc.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if nc.name() != diags.findfieldName(iStash):
            nc.rename(diags.findfieldName(iStash))

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
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/Cloudnet_data/UM/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        # position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir = 'cloud-fraction-metum-grid/2018/'
    # out_dir3 = 'MET_DATA/'

    ### lwc-adiabatic-metum-grid/2018/20180814_oden_lwc-adiabatic-metum-grid.nc
    ###             -> liquid water content derived using measurements averaged on to model grid
    ### cloud-fraction-metum-grid/2018/20180814_oden_cloud-fraction-metum-grid.nc
    ###             -> cloud fraction both from a forecast model and derived from the high-resolution observations on the grid of that model.

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
    # print 'Loading observations:'
    # filename_obs = obs_root_dir + out_dir3 + 'MetData_Gillian_wTemp1p5m.nc'
    # nc_obs = iris.load(filename_obs)#, global_con, callback)
    # print '...'

    # # -------------------------------------------------------------
    # # Load nc
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin nc read in at ' + time.strftime("%c")
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

    doy = np.arange(225,258)        ## set DOY for full moccha figures
    # doy = np.arange(240,244)        ## set DOY for subset of moccha figures

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    if combine == 0:
        # name = '20180902_oden_metum.nc'
        filename = root_dir + out_dir + names[0] + out_dir[:-6] + '.nc'
        print filename
        print ''

        #### LOAD nc
        if 'var_con' in locals():
            print 'Loading single diagnostic:'
            print var_con
            nc1 = iris.load_nc(filename)#, var_con, callback)
        elif 'global_con' in locals():
            print 'Loading multiple diagnostics:'
            # nc = iris.load_ncs(filename1, global_con)
            nc = iris.load(filename)#, global_con, callback)

            # -------------------------------------------------------------

        print nc
        print ''

        # -------------------------------------------------------------
        # Plot data (5x2 timeseries)
        # -------------------------------------------------------------
        # figure = plot_multicontour_TS(nc, filename, out_dir)

        # -------------------------------------------------------------
        # Plot data (map)
        # -------------------------------------------------------------
        ### select hour to plot
        # date = '20180814'
        # hour = 0
        # figure = plot_cartmap(ship_data, nc, hour, date)#, lon, lat)

    else:
        for i in range(0,len(names)):
            filename = root_dir + out_dir + names[i] + out_dir[:-6] + '.nc'
            print filename
            print ''

            print 'Loading multiple diagnostics:'
            nc = Dataset(filename,'r')#, global_con, callback)

            # print 'i = ' + str(i)
            print ''

            if i == 0:
                data = {}
                data1d = {}
                # data['time'] = []
                # data['time'] = float(filename[-16:-14]) + ((nc[0].dim_coords[0].points)/24.0)
                # timem = float(filename[-16:-14]) + ((nc[0].dim_coords[0].points)/24.0)
                if month_flag == -1:
                    timem = doy[i] + ((nc.variables['time'][:])/24.0)
                else:
                    timem = float(names[i][6:8]) + ((nc.variables['time'][:])/24.0)
                for j in range(0,len(nc.variables.keys())):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    if np.sum(nc.variables[nc.variables.keys()[j]].shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(nc.variables[nc.variables.keys()[j]]) >= 24:  # 1d timeseries only
                        data1d[nc.variables.keys()[j]] = nc.variables[nc.variables.keys()[j]][:]
                    else:                                   # 2d column data
                        data[nc.variables.keys()[j]] = nc.variables[nc.variables.keys()[j]][:]
                # print data[nc[0].var_name]
            else:
                # data['time'] = np.append(data['time'],float(filename[-16:-14]) + ((nc[0].dim_coords[0].points)/24.0))
                if month_flag == -1:
                    timem = np.append(timem, doy[i] + ((nc[0].dim_coords[0].points)/24.0))
                else:
                    timem = np.append(timem,float(filename[-16:-14]) + ((nc[0].dim_coords[0].points)/24.0))
                # print data
                for j in range(0,len(nc.variables.keys())):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    # print 'j = ' + str(j)
                    if np.sum(nc[j].data.shape) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.sum(nc[j].data.shape) >= 24:
                        data1d[nc[j].var_name] = np.append(data1d[nc[j].var_name].data,nc[j].data)
                    else:
                        data[nc[j].var_name] = np.append(data[nc[j].var_name].data,nc[j].data,0)

            # print 'Data dict = ' + str(data['radr_refl'].shape)

        # -------------------------------------------------------------
        # Plot combined column data (5x2 timeseries)
        # -------------------------------------------------------------
        # np.save('working_data', data)
        # figure = plot_multicontour_multidate_TS(timem, data, nc, month_flag, missing_files, out_dir, doy)
                    ### doesn't matter which nc, just needed for dim_coords

        # -------------------------------------------------------------
        # Plot combined CASIM column data (4x3 timeseries)
        # -------------------------------------------------------------
        # np.save('working_data', data)
        # figure = plot_multicontour_multidate_casim_TS(timem, data, nc, month_flag, missing_files, out_dir)
                    ### doesn't matter which nc, just needed for dim_coords

        # -------------------------------------------------------------
        # Plot combined timeseries as lineplot
        # -------------------------------------------------------------
        # figure = plot_line_TS(timem, data1d, nc, month_flag, missing_files, out_dir)
                    ### doesn't matter which nc, just needed for dim_coords + nc structure

        # -------------------------------------------------------------
        # Plot combined timeseries as lineplot
        # -------------------------------------------------------------
        # figure = plot_line_TEMP(timem, data1d, nc, month_flag, missing_files, out_dir, nc_obs, doy)
        # figure = plot_line_RAD(timem, data1d, nc, month_flag, missing_files, out_dir, nc_obs, doy)



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
    # <iris 'nc' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc' of surface_net_downward_longwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc' of surface_net_downward_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc' of toa_incoming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc' of toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>]

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
