###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS nc1
###
###

# from __future__ import print_func1tion
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

def readfile(filename_um):

    import pandas as pd

    # print '******'
    print ''
    print 'Reading .txt file with pandas'
    print ''

    um_data = pd.read_csv(filename_um, sep = " ")
    values = um_data.values

    return um_data

def assignColumns(um_data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def iceDrift(um_data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(um_data.values[:,2]>=14,um_data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(um_data.values[:,2]<=14,um_data.values[:,1]==9),um_data.values[:,3]<=22))
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(um_data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(um_data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(um_data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(um_data.values[Sep_drift_index[0][-1],0:3])
    print 'Whole drift: ' + str(um_data.values[drift_index[0],0:4]) + ' - ' + str(um_data.values[drift_index[-1],0:4])
    print ''

    return drift_index

def inIce(um_data):

    ###################################
    ## DEFINE IN ICE PERIOD
    ###################################
    Aug_inIce = np.where(np.logical_and(um_data.values[:,2]>=3,um_data.values[:,1]==8))
    Sep_inIce = np.where(np.logical_and(um_data.values[:,2]<20,um_data.values[:,1]==9))
    inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    # Aug_inIce = np.where(np.logical_and(np.logical_and(um_data.values[:,2]>=12,um_data.values[:,1]==8),um_data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(um_data.values[:,2]>=13,um_data.values[:,1]==8),um_data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(um_data.values[:,2]<=20,um_data.values[:,1]==9))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(um_data.values[:,2]<=20,um_data.values[:,1]==9),um_data.values[:,3]<=1))
    # inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(um_data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(um_data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(um_data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(um_data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(um_data.values[inIce_index[0],0:4]) + ' - ' + str(um_data.values[inIce_index[-1],0:4])
    print 'CloudNET: ' + str(um_data.values[inIce_index[0],0:4]) + ' - ' + str(um_data.values[inIce_index[-1],0:4])
    print ''
    print 'Mean lon/lat of ship track: (' + str(np.nanmedian(um_data.values[inIce_index,6])) + ', ' + str(np.nanmedian(um_data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(um_data.values[inIce_index[0],6]) + ', ' + str(um_data.values[inIce_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(um_data.values[inIce_index[-1],6]) + ', ' + str(um_data.values[inIce_index[-1],7]) + ')'
    print 'Min/max longitude: ' + str(np.nanmin(um_data.values[inIce_index,6])) + ', ' + str(np.nanmax(um_data.values[inIce_index,6]))
    print 'Min/max latitude: ' + str(np.nanmin(um_data.values[inIce_index,7])) + ', ' + str(np.nanmax(um_data.values[inIce_index,7]))
    print ''

    return inIce_index

def trackShip(um_data, date):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(um_data.values[:,2]==14,um_data.values[:,1]==8),um_data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(um_data.values[:,2]==25,um_data.values[:,1]==8),um_data.values[:,3]==1))
    # trackShip_start = np.where(np.logical_and(np.logical_and(um_data.values[:,2]==int(date[-2:]),um_data.values[:,1]==int(date[-4:-2])),um_data.values[:,3]>=0))
    # trackShip_end = np.where(np.logical_and(np.logical_and(um_data.values[:,2]==(int(date[-2:]) + 1),um_data.values[:,1]==int(date[-4:-2])),um_data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print '******'
    print ''
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(um_data.values[inIce_index,6])) + ', ' + str(np.nanmedian(um_data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(um_data.values[trackShip_index[0],6]) + ', ' + str(um_data.values[trackShip_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(um_data.values[trackShip_index[-1],6]) + ', ' + str(um_data.values[trackShip_index[-1],7]) + ')'
    # print 'Start: ' + str(um_data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(um_data.values[trackShip_end[0][-1],0:4])
    print 'trackShip: ' + str(um_data.values[trackShip_index[0],0:4]) + ' - ' + str(um_data.values[trackShip_index[-1],0:4])
    print ''

    return trackShip_index


def plot_contour_TS(nc1, filename_um): #, lon, lat):

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
    print nc1[diag]
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

    time = nc1[diag].dim_coords[0].points
    height = nc1[diag].dim_coords[1].points

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
    ## create figure and axes instanc1es
    #################################################################
    plt.figure(figsize=(8,6))
    ax = plt.gca()

    # plt.plot(nc1[diag].dim_coords[0].points,nc1[diag][:,0].um_data)        # line plot
    # plt.contourf(nc1[0].um_data)
    # plt.plot(nc1[2][0,:].um_data,height);plt.show()
    #################################################################
    ## plot contour timeseries
    ################################################################
    plt.contourf(time,height,np.transpose(nc1[diag].um_data))
    # plt.pcolormesh(time,height,np.transpose(nc1[2].um_data))
    plt.title(nc1[diag].standard_name + ', ' + str(nc1[diag].units))
    plt.colorbar()
    ax.set_ylim([0, 3000])

    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

def plot_multicontour_TS(nc1, filename_um, um_out_dir): #, lon, lat):

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
    for i in range(0,len(nc1)):
        ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
        if np.sum(nc1[i].um_data.shape) > 24:

                ###################################
                ## CHOOSE DIAGNOSTIC
                ###################################
                diag = i
                print ''
                print 'Diag is: '
                print nc1[diag]

                ### define empty array for nc1 um_data
                um_data = []

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

                time = nc1[diag].dim_coords[0].points
                height = nc1[diag].dim_coords[1].points

                ### if mass mixing ratio, *1e3 to change to g/kg
                # if nc1[diag].var_name[0] == 'q':
                #     um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]*1e3))
                # elif nc1[diag].var_name == 'pressure':
                #     um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]/1e2))
                # else:
                #     um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))

                #################################################################
                ## um_data corrections
                #################################################################
                ### set height limit to consider
                ind = np.where(height<5000)

                if nc1[diag].var_name == 'temperature':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'
                elif nc1[diag].var_name == 'qice':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]*1e3))
                    title = nc1[diag].var_name + ' [g/kg]'
                elif nc1[diag].var_name == 'qliq':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]*1e3))
                    title = nc1[diag].var_name + ' [g/kg]'
                elif nc1[diag].var_name == 'q':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]*1e3))
                    title = nc1[diag].var_name + ' [g/kg]'
                elif nc1[diag].var_name == 'pressure':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]/1e2))
                    title = nc1[diag].var_name + ' [hPa]'
                elif nc1[diag].var_name == 'uwind':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'
                elif nc1[diag].var_name == 'wwind':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'
                elif nc1[diag].var_name == 'radr_refl':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'
                elif nc1[diag].var_name == 'cloud_fraction':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'
                elif nc1[diag].var_name == 'vwind':
                    um_data = np.transpose(np.squeeze(nc1[diag].um_data[:,ind]))
                    title = nc1[diag].var_name + ' [' + str(nc1[diag].units) + ']'

                #################################################################
                ## create figure and axes instanc1es
                #################################################################
                if len(um_data) > 0:
                    l = l + 1 ## inc1rement index for positive um_data association

                    print 'l = ' + str(l)
                    print title

                    plt.subplot(5,2,l+1)
                    ax = plt.gca()

                    #################################################################
                    ## plot timeseries
                    #################################################################
                    # plt.contourf(time,height,np.transpose(nc1[diag].um_data))
                    if nc1[diag].var_name == 'temperature':
                        plt.pcolormesh(time, height[ind], um_data, vmin = 250, vmax = np.nanmax(um_data))
                    elif nc1[diag].var_name == 'uwind':
                        plt.pcolormesh(time, height[ind], um_data, vmin = -20, vmax = 20)
                    elif nc1[diag].var_name == 'vwind':
                        plt.pcolormesh(time, height[ind], um_data, vmin = -20, vmax = 20)
                    elif nc1[diag].var_name == 'wwind':
                        plt.pcolormesh(time, height[ind], um_data, vmin = -0.1, vmax = 0.1)
                    else:
                        plt.pcolormesh(time, height[ind], um_data, vmin = np.nanmin(um_data), vmax = np.nanmax(um_data))

                    #################################################################
                    ## set plot properties
                    #################################################################
                    ### colormaps:
                    if nc1[diag].var_name == 'wwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc1[diag].var_name == 'uwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc1[diag].var_name == 'vwind':
                        plt.set_cmap(mpl_cm.RdBu_r)
                    elif nc1[diag].var_name[0] == 'q':
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

    if um_out_dir[:18] == '5_u-bl616_RA2M_CAS':
        fileout = 'FIGS/' + um_out_dir[:21] + filename_um[-22:-3] + '.png'
    elif um_out_dir[:18] == '4_u-bg610_RA2M_CON':
        fileout = 'FIGS/' + um_out_dir[:19] + filename_um[-22:-3] + '.png'
    plt.savefig(fileout, dpi=300)
    plt.show()


def plot_CvProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print 'Plotting Cv statistics based on melt/freeze up periods:'
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
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(9,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.3, wspace = 0.1)

    # print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    plt.plot(np.nanmean(np.squeeze(um_data['Cv'][melt,:]),0),np.nanmean(np.squeeze(um_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.legend()

    plt.subplot(122)
    plt.plot(np.nanmean(np.squeeze(um_data['Cv'][freeze,:]),0),np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    plt.xlabel('Cloud Fraction')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,10000])
    # plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv_splitSeason.svg'
    plt.savefig(fileout)
    plt.show()

def plot_CvProfiles(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print 'Plotting Cv statistics for whole drift period:'
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
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(6,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan

    plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv.svg'
    plt.savefig(fileout)
    plt.show()

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
        um_dir = '/gws/nopw/j04/nc1as_weather/gyoung/MOCCHA/UM/'
        ship_filename_um = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        um_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'
        ifs_dir = '/home/gillian/MOCCHA/Cloudnet/IFS_DATA/'
        obs_um_dir = '/home/gillian/MOCCHA/ODEN/'
        ship_filename_um = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        um_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        um_dir = '/nfs/a96/MOCCHA/working/gillian/Cloudnet_data/UM/'
        ship_filename_um = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        # position_filename_um = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    um_out_dir = 'cloud-fraction-metum-grid/2018/'
    ifs_out_dir = 'cloud-fraction-ecmwf-grid/2018/'
    # out_dir3 = 'MET_DATA/'

    ### lwc-adiabatic-metum-grid/2018/20180814_oden_lwc-adiabatic-metum-grid.nc1
    ###             -> liquid water content derived using measurements averaged on to model grid
    ### cloud-fraction-metum-grid/2018/20180814_oden_cloud-fraction-metum-grid.nc1
    ###             -> cloud fraction both from a forecast model and derived from the high-resolution observations on the grid of that model.

    print '******'
    print ''
    print 'Identifying .nc1 file: '
    print ''

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data = readfile(ship_filename_um)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    # print 'Loading observations:'
    # filename_um_obs = obs_um_dir + um_out_dir3 + 'MetData_Gillian_wTemp1p5m.nc1'
    # nc1_obs = iris.load(filename_um_obs)#, global_con, callback)
    # print '...'

    # # -------------------------------------------------------------
    # # Load nc1
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin nc1 read in at ' + time.strftime("%c")
    print ' '

    ### -------------------------------------------------------------------------
    ### define input filename_um
    ### -------------------------------------------------------------------------
    # tempnames = ['umnsaa_pa012_r0.nc1','umnsaa_pb012_r0.nc1','umnsaa_pc011_r0.nc1','umnsaa_pd011_r0.nc1','20180812_oden_metum.nc1']
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
            '20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180910_oden_']   ### cloud radar not working

    doy = np.arange(225,258)        ## set DOY for full moccha figures
    # doy = np.arange(240,244)        ## set DOY for subset of moccha figures

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):
        filename_um = um_dir + um_out_dir + names[i] + um_out_dir[:-6] + '.nc'
        filename_ifs = ifs_dir + ifs_out_dir + names[i] + ifs_out_dir[:-6] + '.nc'
        print filename_um
        print filename_ifs
        print ''

        print 'Loading multiple diagnostics:'
        nc1 = Dataset(filename_um,'r')
        nc2 = Dataset(filename_ifs,'r')

        # print 'i = ' + str(i)
        print ''

        #### LOAD IN SPECIFIC DIAGNOSTICS
        if um_out_dir[:-6] == 'cloud-fraction-metum-grid':
            var_list = ['height','Cv','model_iwc','model_lwc','model_temperature','model_Cv_filtered']   ### time always read in separately

        ###     LOAD IN UM DATA FIRST
        if i == 0:
            um_data = {}
            um_data1d = {}
            if month_flag == -1:
                time_um = doy[i] + ((nc1.variables['time'][:])/24.0)
            else:
                time_um = float(names[i][6:8]) + ((nc1.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.sum(nc1.variables[var_list[j]].shape) == 24:  # 1d timeseries only
                    um_data1d[var_list[j]] = nc1.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    um_data[var_list[j]] = nc1.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                time_um = np.append(time_um, doy[i] + ((nc1.variables['time'][:])/24.0))
            else:
                time_um = np.append(time_um,float(filename_um[-16:-14]) + ((nc1.variables['time'][:])/24.0))
            print um_data
            for j in range(0,len(var_list)):
                ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                # print 'j = ' + str(j)
                if np.sum(nc1.variables[var_list[j]].shape) == 24:
                    um_data1d[var_list[j]] = np.append(um_data1d[var_list[j]].data,nc1.variables[var_list[j]][:])
                else:
                    um_data[var_list[j]] = np.append(um_data[var_list[j]].data,nc1.variables[var_list[j]][:],0)
        nc1.close()

        if ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
            var_list = ['height','Cv','model_iwc','model_lwc','model_temperature','model_snow_Cv_filtered']   ### time always read in separately

        ###     LOAD IN IFS DATA
        if i == 0:
            ifs_data = {}
            ifs_data1d = {}
            if month_flag == -1:
                time_ifs = doy[i] + ((nc2.variables['time'][:])/24.0)
            else:
                time_ifs = float(names[i][6:8]) + ((nc2.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.sum(nc2.variables[var_list[j]].shape) == 24:  # 1d timeseries only
                    ifs_data1d[var_list[j]] = nc2.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    ifs_data[var_list[j]] = nc2.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                time_ifs = np.append(time_ifs, doy[i] + ((nc2.variables['time'][:])/24.0))
            else:
                time_ifs = np.append(time_ifs,float(filename_ifs[-16:-14]) + ((nc2.variables['time'][:])/24.0))
            print ifs_data
            for j in range(0,len(var_list)):
                ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                # print 'j = ' + str(j)
                if np.sum(nc2.variables[var_list[j]].shape) == 24:
                    ifs_data1d[var_list[j]] = np.append(ifs_data1d[var_list[j]].data,nc2.variables[var_list[j]][:])
                else:
                    ifs_data[var_list[j]] = np.append(ifs_data[var_list[j]].data,nc2.variables[var_list[j]][:],0)
        nc2.close()

        ### PUT TIME INTO DATA DICTIONARIES FOR EASE
        ifs_data['time'] = time_ifs
        um_data['time'] = time_um


        ######  LOAD ALL DIAGNOSTICS
        # if i == 0:
        #     um_data = {}
        #     um_data1d = {}
        #     # um_data['time'] = []
        #     # um_data['time'] = float(filename_um[-16:-14]) + ((nc1[0].dim_coords[0].points)/24.0)
        #     # time_um = float(filename_um[-16:-14]) + ((nc1[0].dim_coords[0].points)/24.0)
        #     if month_flag == -1:
        #         time_um = doy[i] + ((nc1.variables['time'][:])/24.0)
        #     else:
        #         time_um = float(names[i][6:8]) + ((nc1.variables['time'][:])/24.0)
        #     for j in range(0,len(nc1.variables.keys())):
        #         ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
        #         if np.sum(nc1.variables[nc1.variables.keys()[j]].shape) == 0:     # ignore horizontal_resolution
        #             continue
        #         elif nc1.variables.keys()[j] == 'forecast_time':     # ignore forecast_time
        #             continue
        #         elif nc1.variables.keys()[j] == 'time':     # ignore forecast_time
        #             continue
        #         elif np.sum(nc1.variables[nc1.variables.keys()[j]].shape) == 24:  # 1d timeseries only
        #             um_data1d[nc1.variables.keys()[j]] = nc1.variables[nc1.variables.keys()[j]][:]
        #         else:                                   # 2d column um_data
        #             um_data[nc1.variables.keys()[j]] = nc1.variables[nc1.variables.keys()[j]][:]
        #     # nc1.close()
        #     # np.save('working_um_data', um_data)
        #     # np.save('working_um_data1d', um_data1d)
        # else:
        #     if month_flag == -1:
        #         time_um = np.append(time_um, doy[i] + ((nc1.variables['time'][:])/24.0))
        #     else:
        #         time_um = np.append(time_um,float(filename_um[-16:-14]) + ((nc1.variables['time'][:])/24.0))
        #     print um_data
        #     for j in range(0,len(nc1.variables.keys())):
        #         ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
        #         print 'j = ' + str(j)
        #         if np.sum(nc1.variables[nc1.variables.keys()[j]].shape) == 0:     # ignore horizontal_resolution
        #             continue
        #         elif nc1.variables.keys()[j] == 'forecast_time':     # ignore forecast_time
        #             continue
        #         elif nc1.variables.keys()[j] == 'time':     # ignore time, already defined
        #             continue
        #         elif np.sum(nc1.variables[nc1.variables.keys()[j]].shape) == 24:
        #             um_data1d[nc1.variables.keys()[j]] = np.append(um_data1d[nc1.variables.keys()[j]].um_data,nc1.variables[nc1.variables.keys()[j]][:])
        #         else:
        #             um_data[nc1.variables.keys()[j]] = np.append(um_data[nc1.variables.keys()[j]].um_data,nc1.variables[nc1.variables.keys()[j]][:])
        # nc1.close()



    # -------------------------------------------------------------
    # Save working data for debugging
    # -------------------------------------------------------------
    np.save('working_um_data', um_data)
    np.save('working_ifs_data', ifs_data)
    #### um_data = np.load('working_um_data.npy').item()

    # -------------------------------------------------------------
    # Plot Cv statistics from drift period
    # -------------------------------------------------------------
    figure = plot_CvProfiles(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot Cv statistics based on melt/freeze up
    # -------------------------------------------------------------
    figure = plot_CvProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot combined timeseries as lineplot
    # -------------------------------------------------------------
    # figure = plot_line_TS(time_um, um_data1d, month_flag, missing_files, um_out_dir)
                ### doesn't matter which nc1, just needed for dim_coords + nc1 structure

    # -------------------------------------------------------------
    # Plot combined timeseries as lineplot
    # -------------------------------------------------------------
    # figure = plot_line_TEMP(time_um, um_data1d, nc1, month_flag, missing_files, um_out_dir, nc1_obs, doy)
    # figure = plot_line_RAD(time_um, um_data1d, nc1, month_flag, missing_files, um_out_dir, nc1_obs, doy)



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
    # 15: toa_inc1oming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 16: toa_outgoing_longwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)
    # 17: toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 94; grid_longitude: 95)

    #### 12 AUG ONLY - NO FULL NEST DIAGNOSTICS
    # <iris 'nc1' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc1' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc1' of surface_net_downward_longwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc1' of surface_net_downward_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc1' of toa_inc1oming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'nc1' of toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>]

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
