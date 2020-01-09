###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS nc1
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

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.15,
            hspace = 0.22, wspace = 0.1)

    # print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['Cv'][melt,:]),0),np.nanmean(np.squeeze(um_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['Cv'][melt,:]),0) - np.nanstd(np.squeeze(um_data['Cv'][melt,:]),0),
        np.nanmean(np.squeeze(um_data['Cv'][melt,:]),0) + np.nanstd(np.squeeze(um_data['Cv'][melt,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0) - np.nanstd(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0),
        np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0) + np.nanstd(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0) - np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0),
        np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0) + np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0), color = 'navajowhite', alpha = 0.35)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['Cv'][freeze,:]),0),np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['Cv'][freeze,:]),0) - np.nanstd(np.squeeze(um_data['Cv'][freeze,:]),0),
        np.nanmean(np.squeeze(um_data['Cv'][freeze,:]),0) + np.nanstd(np.squeeze(um_data['Cv'][freeze,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0) - np.nanstd(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0),
        np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0) + np.nanstd(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0) - np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0),
        np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0) + np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0), color = 'navajowhite', alpha = 0.35)
    plt.xlabel('Cloud Fraction')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,10000])
    plt.xlim([0,1])
    # plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv_splitSeason.svg'
    # plt.savefig(fileout)
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

    ### define axis instance
    ax = plt.gca()

    print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan

    plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
        np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
        np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
        np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv_240-250DOY.svg'
    # plt.savefig(fileout)
    plt.show()

def plot_CvProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    ### define axis instance
    ax = plt.gca()

    print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan

    plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
        np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
        np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
        np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv_240-250DOY.svg'
    # plt.savefig(fileout)
    plt.show()

def plot_lwcProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print 'Plotting LWC cloudnet statistics based on melt/freeze up periods:'
    print ''

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.12,
            hspace = 0.22, wspace = 0.14)

    # print um_data.keys()

    #### set flagged um_data to nans
    um_data['lwc'][um_data['lwc'] == -999] = np.nan
    um_data['lwc'][um_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['lwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['lwc'][melt,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)

    plt.xlabel('Liquid water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,6000])
    plt.xlim([0,0.2])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['lwc'][freeze,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.xlabel('Liquid water content [g/m3]')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,6000])
    plt.xlim([0,0.2])
    # plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_LWC_splitSeason.png'
    # plt.savefig(fileout, dpi=300)
    plt.show()

def plot_iwcProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print 'Plotting IWC cloudnet statistics based on melt/freeze up periods:'
    print ''

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.12,
            hspace = 0.22, wspace = 0.14)

    # print um_data.keys()

    #### set flagged um_data to nans
    um_data['iwc'][um_data['iwc'] == -999] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] == -999.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] == -999.0] = np.nan
    ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] == -999.0] = np.nan

    um_data['iwc'][um_data['iwc'] == 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['iwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['iwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['iwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['iwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['iwc'][melt,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)

    plt.xlabel('Ice water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,8000])
    plt.xlim([0,0.05])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(um_data['iwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['iwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['iwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['iwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['iwc'][freeze,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3,
    #     np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.xlabel('Ice water content [g/m3]')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,8000])
    plt.xlim([0,0.05])
    # plt.legend()

    print '******'
    print ''
    print 'Finished plotting! :)'
    print ''

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_IWC_splitSeason.png'
    # plt.savefig(fileout, dpi=300)
    plt.show()

def main():

    import sys
    sys.path.insert(1, '../py_functions/')
        ### include py function in path

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
        misc_dir = '/home/gillian/MOCCHA/UM/DATA/'
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
    misc_out_dir = '5_u-bl661_RA1M_CASIM/OUT_R0/'
    # out_dir3 = 'MET_DATA/'

    ### lwc-adiabatic-metum-grid/2018/
    ###             -> liquid water content derived using measurements averaged on to model grid
    ### cloud-fraction-metum-grid/2018/ + cloud-fraction-ecmwf-grid/2018/
    ###             -> cloud fraction both from a forecast model and derived from the high-resolution observations on the grid of that model.
    ### lwc-scaled-metum-grid/2018/ + lwc-scaled-ecmwf-grid/2018/
    ###             -> dataset contains liquid water content derived using radar/lidar cloud boundaries and liquid water path from dual-wavelength
    ###                 microwave radiometers, averaged on to the grid of a forecast model.
    ###                 It also contains the liquid water content and liquid water path from that model, so may be used to calculate statistics
    ###                 quantifying the model performance.
    ### iwc-Z-T-metum-grid/2018/ + iwc-Z-T-ecmwf-grid/2018/
    ###             -> dataset contains ice water content derived using radar reflectivity and temperature, averaged on to the grid of a forecast
    ###                 model. It also contains the ice water content from that model, so may be used to calculate statistics quantifying the
    ###                 model performance.

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

    moccha_names = [#'20180813_oden_','20180814_oden_','20180815_oden_','20180816_oden_',
            # '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            # '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            # '20180825_oden_','20180826_oden_','20180827_oden_',
            '20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_']#,'20180908_oden_','20180909_oden_',
            # '20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180910_oden_']   ### cloud radar not working

    # doy = np.arange(225,258)        ## set DOY for full moccha figures
    doy = np.arange(240,251)        ## set DOY for subset of moccha figures

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1
    misc_flag = 1       ## flag to compare non-cloudnet model data

    for i in range(0,len(names)):
        filename_um = um_dir + um_out_dir + names[i] + um_out_dir[:-6] + '.nc'
        filename_ifs = ifs_dir + ifs_out_dir + names[i] + ifs_out_dir[:-6] + '.nc'
        if misc_flag == 1: filename_misc = misc_dir + misc_out_dir + names[i] + 'metum.nc'
        print filename_um
        print filename_ifs
        if misc_flag == 1: print filename_misc
        print ''

        print 'Loading multiple diagnostics:'
        nc1 = Dataset(filename_um,'r')
        nc2 = Dataset(filename_ifs,'r')
        if misc_flag == 1: nc3 = Dataset(filename_misc,'r')

        # print 'i = ' + str(i)
        print ''

        #### LOAD IN SPECIFIC DIAGNOSTICS
        if um_out_dir[:-6] == 'cloud-fraction-metum-grid':
            var_list = ['height','Cv','model_Cv_filtered']   ### time always read in separately
        elif um_out_dir[:-6] == 'lwc-scaled-metum-grid':
            var_list = ['height','lwc','model_lwc']   ### time always read in separately
        elif um_out_dir[:-6] == 'iwc-Z-T-metum-grid':
            var_list = ['height','iwc','model_iwc_filtered']   ### time always read in separately

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
            var_list = ['height','Cv','model_snow_Cv_filtered']   ### time always read in separately
        elif ifs_out_dir[:-6] == 'lwc-scaled-metum-grid':
            var_list = ['height','lwc','model_lwc']   ### time always read in separately
        elif ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
            var_list = ['height','iwc','model_snow_iwc_filtered','model_iwc_filtered']   ### time always read in separately

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

        ### -------------------------------------------------------------------------
        ###     LOAD IN MISC DATA IF COMPARING
        ###             Only load in what variables are needed based on IFS file chosen
        ### -------------------------------------------------------------------------
        if misc_flag == 1:
            if ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
                var_list = ['height','cloud_fraction']   ### time always read in separately
            elif ifs_out_dir[:-6] == 'lwc-scaled-metum-grid':
                var_list = ['height','qliq']   ### time always read in separately
            elif ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
                var_list = ['height','qice']   ### time always read in separately

            if i == 0:
                misc_data = {}
                misc_data1d = {}
                if month_flag == -1:
                    time_misc = doy[i] + ((nc3.variables['forecast_time'][:])/24.0)
                else:
                    time_misc = float(names[i][6:8]) + ((nc3.variables['forecast_time'][:])/24.0)
                for j in range(0,len(var_list)):
                    if np.sum(nc3.variables[var_list[j]].shape) == 24:  # 1d timeseries only
                        misc_data1d[var_list[j]] = nc3.variables[var_list[j]][:]
                    else:                                   # 2d column um_data
                        misc_data[var_list[j]] = nc3.variables[var_list[j]][:]
            else:
                if month_flag == -1:
                    time_misc = np.append(time_misc, doy[i] + ((nc3.variables['forecast_time'][:])/24.0))
                else:
                    time_misc = np.append(time_misc,float(filename_misc[-16:-14]) + ((nc3.variables['forecast_time'][:])/24.0))
                print misc_data
                for j in range(0,len(var_list)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    # print 'j = ' + str(j)
                    if np.sum(nc3.variables[var_list[j]].shape) == 24:
                        misc_data1d[var_list[j]] = np.append(misc_data1d[var_list[j]].data,nc3.variables[var_list[j]][:])
                    else:
                        misc_data[var_list[j]] = np.append(misc_data[var_list[j]].data,nc3.variables[var_list[j]][:],0)
            nc3.close()

        ### -------------------------------------------------------------------------
        ### PUT TIME INTO DATA DICTIONARIES FOR EASE
        ### -------------------------------------------------------------------------
        ifs_data['time'] = time_ifs
        um_data['time'] = time_um
        if misc_flag == 1: misc_data['time'] = time_misc

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
    if misc_flag == 1: np.save('working_misc_data', misc_data)
    #### um_data = np.load('working_um_data.npy').item()

    # -------------------------------------------------------------
    # Plot Cv statistics from drift period
    # -------------------------------------------------------------
    # figure = plot_CvProfiles(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot Cv statistics from drift period with a 3rd dataset (not run through cloudnet)
    # -------------------------------------------------------------
    figure = plot_CvProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot Cv statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_CvProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot LWC statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_lwcProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot IWC statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_iwcProfiles_SplitSeason(um_data, ifs_data, month_flag, missing_files, um_out_dir, doy)

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
