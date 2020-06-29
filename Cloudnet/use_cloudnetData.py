###
###
### SCRIPT TO READ IN CLOUDNET DATA
###
###

from __future__ import print_function
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

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY
from readMAT import readMatlabStruct

def readfile(filename_um):

    import pandas as pd

    print ('')
    print ('Reading .txt file with pandas')
    print ('')

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

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(um_data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(um_data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(um_data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(um_data.values[Sep_drift_index[0][-1],0:3])
    print ('Whole drift: ' + str(um_data.values[drift_index[0],0:4]) + ' - ' + str(um_data.values[drift_index[-1],0:4]))
    print ('')

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

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(um_data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(um_data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(um_data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(um_data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(um_data.values[inIce_index[0],0:4]) + ' - ' + str(um_data.values[inIce_index[-1],0:4])
    print ('CloudNET: ' + str(um_data.values[inIce_index[0],0:4]) + ' - ' + str(um_data.values[inIce_index[-1],0:4]))
    print ('')
    print ('Mean lon/lat of ship track: (' + str(np.nanmedian(um_data.values[inIce_index,6])) + ', ' + str(np.nanmedian(um_data.values[inIce_index,7])) + ')')
    print ('Lon/lat of start point: (' + str(um_data.values[inIce_index[0],6]) + ', ' + str(um_data.values[inIce_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(um_data.values[inIce_index[-1],6]) + ', ' + str(um_data.values[inIce_index[-1],7]) + ')')
    print ('Min/max longitude: ' + str(np.nanmin(um_data.values[inIce_index,6])) + ', ' + str(np.nanmax(um_data.values[inIce_index,6])))
    print ('Min/max latitude: ' + str(np.nanmin(um_data.values[inIce_index,7])) + ', ' + str(np.nanmax(um_data.values[inIce_index,7])))
    print ('')

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

    print ('******')
    print ('')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(um_data.values[inIce_index,6])) + ', ' + str(np.nanmedian(um_data.values[inIce_index,7])) + ')'
    print ('Lon/lat of start point: (' + str(um_data.values[trackShip_index[0],6]) + ', ' + str(um_data.values[trackShip_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(um_data.values[trackShip_index[-1],6]) + ', ' + str(um_data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(um_data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(um_data.values[trackShip_end[0][-1],0:4])
    print ('trackShip: ' + str(um_data.values[trackShip_index[0],0:4]) + ' - ' + str(um_data.values[trackShip_index[-1],0:4]))
    print ('')

    return trackShip_index

def plot_LWP(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    ###################################
    ## PLOT TIMESERIES
    ###################################

    print ('******')
    print ('')
    print ('Plotting LWP timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	SET UP FIGURE PROPERTIES
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
    plt.figure(figsize=(8,4.5))
    plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged and bad data to nans
    um_data['model_lwp'][um_data['model_lwp'] < 0] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] < 0] = np.nan
    misc_data['model_lwp'][misc_data['model_lwp'] < 0] = np.nan
    # um_data['model_lwp'][um_data['model_lwp'] >= 1000] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] >= 1.0] = np.nan
    # misc_data['model_lwp'][misc_data['model_lwp'] >= 1000] = np.nan
    obs_data['lwp'][obs_data['lwp'][:,0] < 0, 0] = np.nan     ### index 0 is mean
    obs_data['lwp'][obs_data['lwp'][:,0] > 0.8, 0] = np.nan    ### >0.8 == >800g/m2

    # plt.plot(obs_data['time'][:],obs_data['lwp'][:,0]*1e3, 'k', label = 'Obs')
    plt.plot(obs_data['deck7th']['doy'][:],obs_data['deck7th']['lwp'][:], 'k', label = 'Obs_HATPRO')
    plt.plot(um_data['time'][::3],um_data['model_lwp'][::3]*1e3,
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = 'UM_RA2M')
    plt.plot(ifs_data['time'][::3],ifs_data['model_lwp'][::3]*1e3,
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown', label = 'ECMWF_IFS')
    plt.plot(misc_data['time'][::3],misc_data['model_lwp'][::3]*1e3,
        'v', color = 'forestgreen', markeredgecolor = 'darkgreen', label = 'UM_CASIM-100')
    plt.xlabel('Day of Year')
    plt.ylabel('LWP [g/m2]')
    # plt.ylim([0,10000])
    plt.xlim([doy[0],doy[-1]])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-HATPRO_UM_IFS_CASIM-100_LWP_226-257DOY.svg'
    plt.savefig(fileout)
    plt.show()

def plot_CvProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting Cv statistics based on melt/freeze up periods:')
    print ('')

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
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.15,
            hspace = 0.22, wspace = 0.1)

    # print um_data.keys()

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0),np.nanmean(np.squeeze(obs_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][melt,:]),0),np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][melt,:]),0),
        np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][melt,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0) - np.nanstd(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0),
        np.nanmean(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0) + np.nanstd(np.squeeze(um_data['model_Cv_filtered'][melt,:]),0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0) - np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0),
        np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0) + np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][melt,:]),0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][melt,:]),0),np.nanmean(np.squeeze(misc_data['height'][melt,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax1.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][melt,:]),0),np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][melt,:]),0) - np.nanstd(np.squeeze(misc_data['model_Cv_filtered'][melt,:]),0),
        np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][melt,:]),0) + np.nanstd(np.squeeze(misc_data['model_Cv_filtered'][melt,:]),0), color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][freeze,:]),0),
        np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][freeze,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0) - np.nanstd(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0),
        np.nanmean(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0) + np.nanstd(np.squeeze(um_data['model_Cv_filtered'][freeze,:]),0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0) - np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0),
        np.nanmean(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0) + np.nanstd(np.squeeze(ifs_data['model_snow_Cv_filtered'][freeze,:]),0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][freeze,:]),0),np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0),np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][freeze,:]),0) - np.nanstd(np.squeeze(misc_data['model_Cv_filtered'][freeze,:]),0),
        np.nanmean(np.squeeze(misc_data['model_Cv_filtered'][freeze,:]),0) + np.nanstd(np.squeeze(misc_data['model_Cv_filtered'][freeze,:]),0), color = 'mediumaquamarine', alpha = 0.15)
    plt.xlabel('Cloud Fraction')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,10000])
    plt.xlim([0,1])
    # plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_Cv_splitSeason_226-257DOY.svg'
    plt.savefig(fileout)
    plt.show()

def plot_CvProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting Cv statistics for whole drift period:')
    print ('')

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
    plt.figure(figsize=(6,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    plt.plot(np.nanmean(obs_data['Cv'],0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
        np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
        np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
        np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0),np.nanmean(misc_data['height'],0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_Cv_filtered'],0) - np.nanstd(misc_data['model_Cv_filtered'],0),
        np.nanmean(misc_data['model_Cv_filtered'],0) + np.nanstd(misc_data['model_Cv_filtered'],0), color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_Cv_226-257DOY.svg'
    plt.savefig(fileout)
    plt.show()

def plot_CvProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT PROFILE
    ###################################

    print ('******')
    print ('')
    print ('Plotting Cv statistics for whole drift period:')
    print ('')

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
    plt.figure(figsize=(6,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    # misc_data['cloud_fraction'][misc_data['cloud_fraction'] < 0.0] = np.nan

    plt.plot(np.nanmean(obs_data['Cv'],0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
        np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
        np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
        np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['cloud_fraction'],0),misc_data['height'], color = 'forestgreen', linewidth = 3, label = 'UM_GLM')
    ax.fill_betweenx(misc_data['height'],np.nanmean(misc_data['cloud_fraction'],0) - np.nanstd(misc_data['cloud_fraction'],0),
        np.nanmean(misc_data['cloud_fraction'],0) + np.nanstd(misc_data['cloud_fraction'],0), color = 'mediumaquamarine', alpha = 0.15)

    # plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'Obs_UM')
    # ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
    #     np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightblue', alpha = 0.5)
    # plt.plot(np.nanmean(ifs_data['Cv'],0),np.nanmean(ifs_data['height'],0), color='darkorange', linewidth = 3, label = 'Obs_IFS')
    # ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['Cv'],0) - np.nanstd(ifs_data['Cv'],0),
    #     np.nanmean(ifs_data['Cv'],0) + np.nanstd(ifs_data['Cv'],0), color = 'navajowhite', alpha = 0.5)
    # plt.plot(np.nanmean(obs_data['Cv'],0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    # ax.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
    #     np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM-RA2M_IFS_UM-RA2T-LAM_Cv.png'
    plt.savefig(fileout)
    plt.show()

def plot_lwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting LWC cloudnet statistics based on melt/freeze up periods:')
    print ('')

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
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.12,
            hspace = 0.22, wspace = 0.14)

    # print um_data.keys()

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
    misc_data['model_lwc'][misc_data['model_lwc'] <= 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][melt,:]),0),np.nanmean(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][melt,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax1.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][melt,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Liquid water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.xlim([0,0.2])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][freeze,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Liquid water content [g/m3]')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,10000])
    plt.xlim([0,0.2])
    # plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_LWC_splitSeason.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_lwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting LWC statistics for whole drift period:')
    print ('')

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
    plt.figure(figsize=(6,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    obs_data['lwc'][obs_data['lwc'] >= 1e-3] = np.nan       ## exclude >1g/m3
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
    misc_data['model_lwc'][misc_data['model_lwc'] <= 0] = np.nan
    # um_data['lwc'][um_data['lwc'] <= 0.0] = np.nan
    # ifs_data['lwc'][ifs_data['lwc'] <= 0.0] = np.nan
    # ifs_data['lwc'][ifs_data['lwc'] >= 20.0] = np.nan
    # misc_data['lwc'][misc_data['lwc'] <= 0] = np.nan

    print (np.nanmean(obs_data['lwc'],0).shape)
    print (obs_data['height'].shape)

    # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
        np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_lwc'],0)*1e3,np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc']*1e3,0),
        np.nanmean(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_lwc'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3,
        np.nanmean(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['model_lwc'],0)*1e3,np.nanmean(misc_data['height'],0), color = 'forestgreen', linewidth = 3, label = 'CASIM-100')
    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_lwc'],0)*1e3 - np.nanstd(misc_data['model_lwc']*1e3,0),
        np.nanmean(misc_data['model_lwc'],0)*1e3 + np.nanstd(misc_data['model_lwc'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    # plt.plot(np.nanmean(um_data['lwc'],0)*1e3,um_data['height'], color = 'steelblue', linewidth = 3, label = 'UM')
    # ax1.fill_betweenx(um_data['height'],np.nanmean(um_data['lwc'],0)*1e3 - np.nanstd(um_data['lwc']*1e3,0),
    #     np.nanmean(um_data['lwc'],0)*1e3 + np.nanstd(um_data['lwc'],0)*1e3, color = 'lightblue', alpha = 0.4)
    # plt.plot(np.nanmean(ifs_data['lwc'],0)*1e3,ifs_data['height'], color = 'darkorange', linewidth = 3, label = 'IFS')
    # ax1.fill_betweenx(ifs_data['height'],np.nanmean(ifs_data['lwc'],0)*1e3 - np.nanstd(ifs_data['lwc'],0)*1e3,
    #     np.nanmean(ifs_data['lwc'],0)*1e3 + np.nanstd(ifs_data['lwc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(misc_data['lwc'],0)*1e3,misc_data['height'], color = 'forestgreen', linewidth = 3, label = 'CASIM-100')
    # ax1.fill_betweenx(misc_data['height'],np.nanmean(misc_data['lwc'],0)*1e3 - np.nanstd(misc_data['lwc']*1e3,0),
    #     np.nanmean(misc_data['lwc'],0)*1e3 + np.nanstd(misc_data['lwc'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Liquid water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,1.0])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_LWC-InstRes_226-257DOY.png'
    # plt.savefig(fileout)
    plt.show()

def plot_lwcProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting LWC statistics for whole drift period:')
    print ('')

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
    ax1 = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
    misc_data['qliq'][misc_data['qliq'] <= 0] = np.nan

    plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
        np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_lwc'],0)*1e3,np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc']*1e3,0),
        np.nanmean(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_lwc'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3,
        np.nanmean(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['qliq'],0)*1e3,misc_data['height'], color = 'forestgreen', linewidth = 3, label = 'UM_RA2T')
    ax1.fill_betweenx(misc_data['height'],(np.nanmean(misc_data['qliq'],0) - np.nanstd(misc_data['qliq'],0))*1e3,
        (np.nanmean(misc_data['qliq'],0) + np.nanstd(misc_data['qliq'],0))*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Liquid water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,0.2])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM-RA2M_IFS_LWC_UM-RA2T_qliq.png'
    plt.savefig(fileout)
    plt.show()

def plot_iwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting IWC cloudnet statistics based on melt/freeze up periods:')
    print ('')

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
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.12,
            hspace = 0.22, wspace = 0.14)

    # print um_data.keys()

    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] == -999.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] == -999.0] = np.nan
    ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] == -999.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] == -999.0] = np.nan

    obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    # ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] <= 0.0] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][melt,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['iwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][melt,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][melt,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][melt,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][melt,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax1.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][melt,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][melt,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][melt,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Ice water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.xlim([0,0.05])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs')
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['iwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][freeze,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][freeze,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][freeze,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Ice water content [g/m3]')
    plt.title('Freeze up')
    plt.yticks([])
    plt.ylim([0,10000])
    plt.xlim([0,0.05])
    # plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_IWC_splitSeason.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_iwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting IWC statistics for whole drift period:')
    print ('')

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
    plt.figure(figsize=(6,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
    obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] <= 0.0] = np.nan

    # plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
    #     np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
        np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_iwc_filtered'],0)*1e3,np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered']*1e3,0),
        np.nanmean(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3,
        np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered'],0)*1e3,np.nanmean(misc_data['height'],0), color = 'forestgreen', linewidth = 3, label = 'CASIM-100')
    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 - np.nanstd(misc_data['model_iwc_filtered']*1e3,0),
        np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 + np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Ice water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,0.05])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_IWC-InstRes_226-257DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_iwcProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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

    print ('******')
    print ('')
    print ('Plotting IWC statistics for whole drift period:')
    print ('')

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
    ax1 = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
    obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
    misc_data['qice'][misc_data['qice'] <= 0] = np.nan

    plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
        np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_iwc_filtered'],0)*1e3,np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered']*1e3,0),
        np.nanmean(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3,
        np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['qice'],0)*1e3,misc_data['height'], color = 'forestgreen', linewidth = 3, label = 'UM_RA2T')
    ax1.fill_betweenx(misc_data['height'],(np.nanmean(misc_data['qice'],0) - np.nanstd(misc_data['qice'],0))*1e3,
        (np.nanmean(misc_data['qice'],0) + np.nanstd(misc_data['qice'],0))*1e3, color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Ice water content [g/m3]')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    plt.xlim([0,0.05])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM-RA2M_IFS_IWC_UM-RA2T_qice.png'
    plt.savefig(fileout)
    plt.show()

def plot_TempProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT PROFILE
    ###################################

    print ('******')
    print ('')
    print ('Plotting Cv statistics for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	SET FONT SIZE AND AXES INSTANCE
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

    print (um_data.keys())
    #### set flagged um_data to nans
    um_data['model_temperature'][um_data['model_temperature'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    # um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_temperature'][ifs_data['model_temperature'] == -999] = np.nan
    misc_data['temperature'][misc_data['temperature'] == -9999] = np.nan

    # plt.plot(np.nanmean(um_data['model_temperature'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    # ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_temperature'],0) - np.nanstd(um_data['Cv'],0),
    #     np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(um_data['model_temperature'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_temperature'],0) - np.nanstd(um_data['model_temperature'],0),
        np.nanmean(um_data['model_temperature'],0) + np.nanstd(um_data['model_temperature'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_temperature'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'IFS')
    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_temperature'],0) - np.nanstd(ifs_data['model_temperature'],0),
        np.nanmean(ifs_data['model_temperature'],0) + np.nanstd(ifs_data['model_temperature'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(misc_data['temperature'],0),misc_data['height'], color = 'forestgreen', linewidth = 3, label = 'CASIM-AeroProf')
    ax.fill_betweenx(misc_data['height'],np.nanmean(misc_data['temperature'],0) - np.nanstd(misc_data['temperature'],0),
        np.nanmean(misc_data['temperature'],0) + np.nanstd(misc_data['temperature'],0), color = 'mediumaquamarine', alpha = 0.15)

    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    # plt.xlim([0,1])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_Cv_CASIM-AeroProf_temp_229-250DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_CvTimeseries_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Cv timeseries for whole drift period:')
    print ('')

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
    plt.figure(figsize=(14,12))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['cloud_fraction'][misc_data['cloud_fraction'] < 0.0] = np.nan

    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr
    newcmp = ListedColormap(newcolors)

    plt.subplot(411)
    plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['Cv']),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('Measured cloud fraction by volume, 1 hour sampling')
    plt.colorbar()

    plt.subplot(412)
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_Cv_filtered']),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_RA2M; modelled cloud fraction')
    plt.colorbar()

    plt.subplot(413)
    plt.contourf(misc_data['time'], misc_data['height'][:], np.transpose(misc_data['cloud_fraction']),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_CASIM-100; modelled cloud fraction')
    plt.colorbar()

    plt.subplot(414)
    plt.contourf(ifs_data['time'], np.squeeze(ifs_data['height'][0,:]), np.transpose(ifs_data['model_snow_Cv_filtered']),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlabel('DOY')
    plt.title('ECMWF_IFS; modelled cloud fraction (including snow)')
    plt.colorbar()


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_AP_CvTimeseries_226-257DOY.svg'
    plt.savefig(fileout)
    plt.show()

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
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        um_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename_um = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        um_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'
        ifs_dir = '/home/gillian/MOCCHA/Cloudnet/IFS_DATA/'
        misc_dir = '/home/gillian/MOCCHA/UM/DATA/'                ### FOR NON-CLOUDNET UM DATA
        # misc_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'        ### FOR CLOUDNET UM DATA
        obs_dir = '/home/gillian/MOCCHA/Cloudnet/OBS_DATA/QF10_ecmwf/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename_um = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        um_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        um_dir = '/nfs/a96/MOCCHA/working/gillian/Cloudnet_data/UM/'
        ship_filename_um = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        # position_filename_um = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    um_out_dir = '4_u-bg610_RA2M_CON/cloud-fraction-metum-grid/2018/'
    ifs_out_dir = 'cloud-fraction-ecmwf-grid/2018/'
    obs_out_dir = ifs_out_dir
    if misc_dir == '/home/gillian/MOCCHA/Cloudnet/UM_DATA/':
        misc_out_dir = '5_u-bl661_RA1M_CASIM/lwc-scaled-metum-grid/2018/'
        misc_flag = 0       ## flag to compare cloudnet model data
    elif misc_dir == '/home/gillian/MOCCHA/UM/DATA/':
        misc_out_dir = '7_u-bn068_RA2T_CON/OUT_R1/'
        misc_flag = 1       ## flag to compare non-cloudnet model data

    print ('Misc_flag = ' + str(misc_flag) + '... so third simulation for comparison is:')
    if misc_flag == 0: print ('Cloudnet-ed data!')
    if misc_flag == 1: print ('standard model output!')

    ######## lwc-adiabatic-metum-grid/2018/
    ########             -> liquid water content derived using measurements averaged on to model grid
    ### cloud-fraction-metum-grid/2018/ + cloud-fraction-ecmwf-grid/2018/
    ###             -> cloud fraction both from a forecast model and derived from the high-resolution observations on the grid of that model.
    ### lwc-scaled-metum-grid/2018/ + lwc-scaled-ecmwf-grid/2018/ + lwc-scaled-adiabatic/2018/
    ###             -> dataset contains liquid water content derived using radar/lidar cloud boundaries and liquid water path from dual-wavelength
    ###                 microwave radiometers, averaged on to the grid of a forecast model.
    ###                 It also contains the liquid water content and liquid water path from that model, so may be used to calculate statistics
    ###                 quantifying the model performance.
    ### iwc-Z-T-metum-grid/2018/ + iwc-Z-T-ecmwf-grid/2018/ + iwc-Z-T-method/2018/
    ###             -> dataset contains ice water content derived using radar reflectivity and temperature, averaged on to the grid of a forecast
    ###                 model. It also contains the ice water content from that model, so may be used to calculate statistics quantifying the
    ###                 model performance.

    print ('******')
    print ('')
    print ('Identifying .nc1 file: ')
    print ('')

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
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
    print ('******')
    print ('')
    print ('Begin nc1 read in at ' + time.strftime("%c"))
    print (' ')

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

    moccha_names = [#'20180814_oden_','20180815_oden_','20180816_oden_',
            # '20180817_oden_','20180819_oden_','20180820_oden_',
            # '20180821_oden_','20180822_oden_','20180823_oden_',
            # '20180825_oden_','20180826_oden_','20180827_oden_',
            '20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_',
            '20180901_oden_','20180902_oden_','20180903_oden_','20180904_oden_']#,'20180905_oden_',
            # '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            # '20180911_oden_','20180912_oden_','20180913_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180813_oden_','20180818_oden_','20180824_oden_','20180910_oden_','20180914_oden_']   ### cloud radar not working

    # doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)        ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)
    # doy = np.arange(229,259)        ## set DOY for CASIM-AeroProf (17th Aug to 14th Sep)
    # doy = np.arange(226,259)        ## set DOY for CASIM-100_AP (1st Sep to 9th Sep)

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):

        ### --------------------------------------------------------------------
        ###     DEFINE FILENAMES
        ### --------------------------------------------------------------------
        if um_out_dir[-31:-6] == 'cloud-fraction-metum-grid':
            out_dir = 'cloud-fraction-metum-grid'
        elif um_out_dir[-27:-6] == 'lwc-scaled-metum-grid':
            out_dir = 'lwc-scaled-metum-grid'
        elif um_out_dir[-26:-6] == 'lwc-scaled-adiabatic':
            out_dir = 'lwc-scaled-adiabatic'
        elif um_out_dir[-24:-6] == 'iwc-Z-T-metum-grid':
            out_dir = 'iwc-Z-T-metum-grid'
        filename_um = um_dir + um_out_dir + names[i] + out_dir + '.nc'
        filename_ifs = ifs_dir + ifs_out_dir + names[i] + ifs_out_dir[:-6] + '.nc'
        filename_obs = obs_dir + obs_out_dir + names[i] + obs_out_dir[:-6] + '.nc'
        if misc_flag == 1: filename_misc = misc_dir + misc_out_dir + names[i] + 'metum.nc'
        if misc_flag == 0: filename_misc = misc_dir + misc_out_dir + names[i] + out_dir + '.nc'
        print (filename_um)
        print (filename_ifs)
        if misc_flag == 1: print (filename_misc)
        print ('')

        ### --------------------------------------------------------------------
        ###     READ IN ALL FILES
        ### --------------------------------------------------------------------

        print ('Loading multiple diagnostics:')
        nc1 = Dataset(filename_um,'r')
        nc2 = Dataset(filename_ifs,'r')
        if misc_flag != -1: nc3 = Dataset(filename_misc,'r')
        nc4 = Dataset(filename_obs,'r')

        # print 'i = ' + str(i)
        print ('')

        ### --------------------------------------------------------------------
        ###     LOAD UM DATA INTO DICTIONARY
        ### --------------------------------------------------------------------
        #### LOAD IN SPECIFIC DIAGNOSTICS
        if out_dir == 'cloud-fraction-metum-grid':
            var_list = ['height','Cv','model_Cv_filtered','model_temperature']   ### time always read in separately
        elif out_dir == 'lwc-scaled-metum-grid':
            var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
        elif out_dir == 'lwc-scaled-adiabatic':
            var_list = ['height','lwc']   ### time always read in separately
        elif out_dir == 'iwc-Z-T-metum-grid':
            var_list = ['height','iwc','model_iwc_filtered']   ### time always read in separately

        ###     LOOP OVER TIME DUMP
        if i == 0:
            um_data = {}
            # um_data1d = {}
            if month_flag == -1:
                time_um = doy[i] + ((nc1.variables['time'][:])/24.0)
            else:
                time_um = float(names[i][6:8]) + ((nc1.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.ndim(nc1.variables[var_list[j]]) == 1:  # 1d timeseries only
                    um_data[var_list[j]] = nc1.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    um_data[var_list[j]] = nc1.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                time_um = np.append(time_um, doy[i] + ((nc1.variables['time'][:])/24.0))
            else:
                time_um = np.append(time_um,float(filename_um[-16:-14]) + ((nc1.variables['time'][:])/24.0))
            print (um_data)
            for j in range(0,len(var_list)):
                # print 'j = ' + str(j)
                if np.ndim(nc1.variables[var_list[j]]) == 1:
                    um_data[var_list[j]] = np.append(um_data[var_list[j]].data,nc1.variables[var_list[j]][:])
                else:
                    um_data[var_list[j]] = np.append(um_data[var_list[j]].data,nc1.variables[var_list[j]][:],0)
        nc1.close()

        ### --------------------------------------------------------------------
        ### LOAD IN IFS DATA INTO DICTIONARY
        ### --------------------------------------------------------------------
        if ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
            var_list = ['height','Cv','model_snow_Cv_filtered','model_temperature']   ### time always read in separately
        elif ifs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
            var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
        elif out_dir == 'lwc-scaled-adiabatic':
            var_list = ['height','lwc']   ### time always read in separately
        elif ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
            var_list = ['height','iwc','model_snow_iwc_filtered','model_iwc_filtered']   ### time always read in separately

        ###     LOOP OVER TIME DUMP
        if i == 0:
            ifs_data = {}
            # ifs_data1d = {}
            if month_flag == -1:
                time_ifs = doy[i] + ((nc2.variables['time'][:])/24.0)
            else:
                time_ifs = float(names[i][6:8]) + ((nc2.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.ndim(nc2.variables[var_list[j]]) == 1:  # 1d timeseries only
                    ifs_data[var_list[j]] = nc2.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    ifs_data[var_list[j]] = nc2.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                time_ifs = np.append(time_ifs, doy[i] + ((nc2.variables['time'][:])/24.0))
            else:
                time_ifs = np.append(time_ifs,float(filename_ifs[-16:-14]) + ((nc2.variables['time'][:])/24.0))
            print (ifs_data)
            for j in range(0,len(var_list)):
                ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                # print 'j = ' + str(j)
                if np.ndim(nc2.variables[var_list[j]]) == 1:
                    ifs_data[var_list[j]] = np.append(ifs_data[var_list[j]].data,nc2.variables[var_list[j]][:])
                else:
                    ifs_data[var_list[j]] = np.append(ifs_data[var_list[j]].data,nc2.variables[var_list[j]][:],0)
        nc2.close()

        ### -------------------------------------------------------------------------
        ###     LOAD IN MISC DATA INTO DICTIONARY IF COMPARING
        ###             Only load in what variables are needed based on IFS file chosen
        ### -------------------------------------------------------------------------
        if misc_flag == -1:
            continue
        elif misc_flag == 1:
            if ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
                var_list = ['cloud_fraction','temperature']   ### time always read in separately
            elif ifs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
                var_list = ['qliq']   ### time always read in separately
            elif ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
                var_list = ['qice']   ### time always read in separately
        elif misc_flag == 0:
            if out_dir == 'cloud-fraction-metum-grid':
                var_list = ['height','Cv','model_Cv_filtered','model_temperature']   ### time always read in separately
            elif out_dir == 'lwc-scaled-metum-grid':
                var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
            elif out_dir == 'iwc-Z-T-metum-grid':
                var_list = ['height','iwc','model_iwc_filtered']   ### time always read in separately

        print ('')
        print ('misc file variable list is:')
        print (var_list)
        print ('')

        if i == 0:
            misc_data = {}
            # misc_data1d = {}
            if month_flag == -1:
                if misc_flag == 1:
                    time_misc = doy[i] + ((nc3.variables['forecast_time'][:])/24.0)
                    misc_data['height'] = nc3.variables['height'][:]
                if misc_flag == 0: time_misc = doy[i] + ((nc3.variables['time'][:])/24.0)
            else:
                if misc_flag == 1: time_misc = float(names[i][6:8]) + ((nc3.variables['forecast_time'][:])/24.0)
                if misc_flag == 0: time_misc = float(names[i][6:8]) + ((nc3.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.ndim(nc3.variables[var_list[j]]) == 1:  # 1d timeseries only
                    misc_data[var_list[j]] = nc3.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    misc_data[var_list[j]] = nc3.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                if misc_flag == 1: time_misc = np.append(time_misc, doy[i] + ((nc3.variables['forecast_time'][:])/24.0))
                if misc_flag == 0: time_misc = np.append(time_misc, doy[i] + ((nc3.variables['time'][:])/24.0))
            else:
                if misc_flag == 1: time_misc = np.append(time_misc,float(filename_misc[-16:-14]) + ((nc3.variables['forecast_time'][:])/24.0))
                if misc_flag == 0: time_misc = np.append(time_misc,float(filename_misc[-16:-14]) + ((nc3.variables['time'][:])/24.0))
            print (misc_data)
            for j in range(0,len(var_list)):
                # print 'j = ' + str(j)
                if np.ndim(nc3.variables[var_list[j]]) == 1:
                    misc_data[var_list[j]] = np.append(misc_data[var_list[j]].data,nc3.variables[var_list[j]][:])
                # elif var_list[j] == 'height':#np.sum(nc3.variables[var_list[j]].shape) == 71:
                #     continue
                else:
                    misc_data[var_list[j]] = np.append(misc_data[var_list[j]].data,nc3.variables[var_list[j]][:],0)
        nc3.close()

        ### -------------------------------------------------------------------------
        ###     LOAD IN OBS DATA
        ###             Only load in what variables are needed based on IFS file chosen
        ### -------------------------------------------------------------------------
        if obs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
            var_list = ['height','Cv']   ### time always read in separately
        elif obs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
            var_list = ['height','lwc','lwp']   ### time always read in separately
        elif obs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
            var_list = ['height','iwc']   ### time always read in separately
        elif obs_out_dir == 'lwc-scaled-adiabatic/2018/':
            var_list = ['height','lwc','lwp']   ### time always read in separately
        elif obs_out_dir == 'iwc-Z-T-method/2018/':
            var_list = ['height','iwc']   ### time always read in separately

        if i == 0:
            obs_data = {}
            # misc_data1d = {}
            if month_flag == -1:
                time_obs = doy[i] + ((nc4.variables['time'][:])/24.0)
            else:
                time_obs = float(names[i][6:8]) + ((nc4.variables['time'][:])/24.0)
            for j in range(0,len(var_list)):
                if np.ndim(nc4.variables[var_list[j]]) == 1:  # 1d timeseries only
                    obs_data[var_list[j]] = nc4.variables[var_list[j]][:]
                else:                                   # 2d column um_data
                    obs_data[var_list[j]] = nc4.variables[var_list[j]][:]
        else:
            if month_flag == -1:
                time_obs = np.append(time_obs, doy[i] + ((nc4.variables['time'][:])/24.0))
            else:
                time_obs = np.append(time_obs,float(filename_obs[-16:-14]) + ((nc4.variables['time'][:])/24.0))
            print (obs_data)
            for j in range(0,len(var_list)):
                # print 'j = ' + str(j)
                if np.ndim(nc4.variables[var_list[j]]) == 1:
                    obs_data[var_list[j]] = np.append(obs_data[var_list[j]].data,nc4.variables[var_list[j]][:])
                elif np.sum(nc4.variables[var_list[j]].shape) == 71:
                    continue
                else:
                    obs_data[var_list[j]] = np.append(obs_data[var_list[j]].data,nc4.variables[var_list[j]][:],0)
        nc4.close()

        ### -------------------------------------------------------------------------
        ### PUT TIME INTO DATA DICTIONARIES FOR EASE
        ### -------------------------------------------------------------------------
        ifs_data['time'] = time_ifs
        um_data['time'] = time_um
        if misc_flag != -1: misc_data['time'] = time_misc
        obs_data['time'] = time_obs

    ### -------------------------------------------------------------------------
    ### Load in other measurement data
    ### -------------------------------------------------------------------------

    print ('Load temporary ice station data from Jutta...')
    obs_data['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_wTemp1p5m.nc','r')

    print ('Load ice station data from Jutta...')
    obs_data['ice_station'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')
            #### mast_radiation_30min_v2.3.mat
            #### flux30_trhwxrel.mat

    print ('Load radiosonde data from Jutta...')
    obs_data['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print ('Load observations inversion height data from Jutta...')
    obs_data['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/InversionHeights_RSh05int_final_V03.mat')

    print ('Load foremast data from John...')
    obs_data['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print ('Load 7th deck weather station data from John...')
    obs_data['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')


    # -------------------------------------------------------------
    # Save working data for debugging
    # -------------------------------------------------------------
    np.save('working_um_data', um_data)
    np.save('working_ifs_data', ifs_data)
    if misc_flag != -1: np.save('working_misc_data', misc_data)
    # np.save('working_obs_data', obs_data)
    #### um_data = np.load('working_um_data.npy').item()

    # -------------------------------------------------------------
    # Plot timeseries of 1D diagnostics from drift period
    # -------------------------------------------------------------
    # figure = plot_LWP(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot Cv statistics from drift period
    # -------------------------------------------------------------
    # figure = plot_CvProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_lwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_iwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot statistics from drift period with a 3rd dataset (not run through cloudnet)
    # -------------------------------------------------------------
    figure = plot_CvProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_lwcProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_iwcProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_TempProfiles_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)
    # figure = plot_CvTimeseries_3rdNoCloudnet(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot Cv statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_CvProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot LWC statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_lwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

    # -------------------------------------------------------------
    # Plot IWC statistics based on melt/freeze up
    # -------------------------------------------------------------
    # figure = plot_iwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy)

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
