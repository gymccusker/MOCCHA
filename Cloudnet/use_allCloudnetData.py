###
###
### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
###
###     ### test change

from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import numpy as np
# import diags_MOCCHA as diags
# import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import seaborn as sns

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY
from readMAT import readMatlabStruct
from physFuncts import calcThetaE, calcThetaVL
# from conversionFuncts import reGrid_Sondes

def readfile(filename):

    import pandas as pd

    # print '******'
    print( '')
    print( 'Reading .txt file with pandas')
    print( '')

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

    print( '******')
    print( '')
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print( 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print( '')
    print( 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print( 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print( 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print( 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))
    print( '')

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

    print( '******')
    print( '')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print( 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print( 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print( '')

    return trackShip_index

def plot_CvProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs, obs_switch): #, lon, lat):

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
    plt.figure(figsize=(6,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.2,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    obs_data['Cv_adv'][obs_data['Cv_adv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    if obs_switch == 'RADAR':
        # plt.plot(np.nanmean(obs_data['Cv'],0),obs_data['height'], 'k--', linewidth = 3, label = 'Obs-on-model-grid')
        # ax.fill_betweenx(obs_data['height'],np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
        #     np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
        # plt.plot(np.nanmean(obs_data['Cv_adv'],0),np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs')
        # plt.plot(np.nanmean(obs['cloudfractions']['cloudfraction_total'],0),np.squeeze(obs['cloudfractions']['height']),
        #     'k--', linewidth = 3, label = 'Obs-on-artificial-grid')
        ax.fill_betweenx(np.squeeze(obs['cloudfractions']['height']),np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) - np.nanstd(obs['cloudfractions']['cloudfraction_total'],0),
            np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) + np.nanstd(obs['cloudfractions']['cloudfraction_total'],0), color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
        ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
            np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
        plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
        ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
            np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)
        plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0),np.nanmean(misc_data['height'],0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
        ax.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_Cv_filtered'],0) - np.nanstd(misc_data['model_Cv_filtered'],0),
            np.nanmean(misc_data['model_Cv_filtered'],0) + np.nanstd(misc_data['model_Cv_filtered'],0), color = 'mediumaquamarine', alpha = 0.15)
    else:
        plt.plot(np.nanmean(obs_data['Cv_adv'],0),np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs Cv_adv')
        plt.plot(np.nanmean(obs_data['Cv_adv'],0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs Cv')
        ax.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv_adv'],0) - np.nanstd(obs_data['Cv'],0),
            np.nanmean(obs_data['Cv_adv'],0) + np.nanstd(obs_data['Cv_adv'],0), color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(obs['cloudfractions']['cloudfraction_total'],0),np.squeeze(obs['cloudfractions']['height']),
            '--', color = 'grey', linewidth = 3, label = 'Obs-on-artificial-grid')
        # ax.fill_betweenx(np.squeeze(obs['cloudfractions']['height']),np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) - np.nanstd(obs['cloudfractions']['cloudfraction_total'],0),
        #     np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) + np.nanstd(obs['cloudfractions']['cloudfraction_total'],0), color = 'lightgrey', alpha = 0.5)
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
        fileout = 'FIGS/Obs-JVObsGridded_UM_IFS_CASIM-100_Cv_226-257DOY.svg'
    # plt.savefig(fileout)
    plt.show()

def plot_CvTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    obs_data['Cv_adv'][obs_data['Cv_adv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr
    newcmp = ListedColormap(newcolors)

    plt.subplot(411)
    plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['Cv_adv']),
        np.arange(0,1.1,0.1),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.title('Measured cloud fraction by volume, 1 hour sampling')
    plt.title('Measured cloud fraction by volume, 1.5km sampling')
    plt.colorbar()

    plt.subplot(412)
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_Cv_filtered']),
        np.arange(0,1.1,0.1),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_RA2M; modelled cloud fraction')
    plt.colorbar()

    plt.subplot(413)
    plt.contourf(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_Cv_filtered']),
        np.arange(0,1.1,0.1),
        cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_CASIM-100; modelled cloud fraction')
    plt.colorbar()

    plt.subplot(414)
    plt.contourf(ifs_data['time'], np.squeeze(ifs_data['height'][0,:]), np.transpose(ifs_data['model_snow_Cv_filtered']),
        np.arange(0,1.1,0.1),
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
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_CvTimeseries_226-257DOY.png'
    # plt.savefig(fileout)
    plt.show()

def plot_lwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

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


    if obs_switch == 'RADAR':
        plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs')
        ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
            np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.xlim([0,1.0])
    else:
        plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
            np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.xlim([0,0.2])
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
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_LWC-InstRes_226-257DOY.png'
    # plt.savefig(fileout)
    plt.show()

def plot_LWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print ('Plotting LWC timeseries for whole drift period:')
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
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    misc_data['model_lwc'][misc_data['model_lwc'] <= 0.0] = np.nan

    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    plt.subplot(411)
    plt.pcolor(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['lwc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        # cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('Obs')
    plt.colorbar()

    plt.subplot(412)
    plt.pcolor(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_lwc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        # cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_RA2M')
    plt.colorbar()

    plt.subplot(413)
    plt.pcolor(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_lwc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        # cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_CASIM-100')
    plt.colorbar()

    plt.subplot(414)
    plt.pcolor(ifs_data['time'], np.squeeze(ifs_data['height'][0,:]), np.transpose(ifs_data['model_lwc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        # cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlabel('DOY')
    plt.title('ECMWF_IFS')
    plt.colorbar()


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_LWCTimeseries_226-257DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_IWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print ('Plotting IWC timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 12
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

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

    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    plt.subplot(411)
    plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['iwc'])*1e3,
        )#cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('Obs')
    plt.colorbar()

    plt.subplot(412)
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_iwc_filtered'])*1e3,
        )#cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_RA2M')
    plt.colorbar()

    plt.subplot(413)
    plt.contourf(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_iwc_filtered'])*1e3,
        )#cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_CASIM-100')
    plt.colorbar()

    plt.subplot(414)
    plt.contourf(ifs_data['time'], np.squeeze(ifs_data['height'][0,:]), np.transpose(ifs_data['model_snow_iwc_filtered'])*1e3,
        )#cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlabel('DOY')
    plt.title('ECMWF_IFS')
    plt.colorbar()


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_IWCTimeseries_226-257DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_TWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    import matplotlib.colors as colors
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting IWC timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 12
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(8,9))
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.2)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] <= 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    # ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 1.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] <= 0.0] = np.nan

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    misc_data['model_lwc'][misc_data['model_lwc'] <= 0.0] = np.nan

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']

    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    twc0 = np.transpose(obs_data['twc'])*1e3

    plt.subplot(411)
    plt.pcolormesh(obs_data['time'], np.squeeze(obs_data['height'][0,:]), twc0,
        # norm=colors.LogNorm(vmin=0.0, vmax=0.5))
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('Obs')
    plt.colorbar()

    plt.subplot(412)
    plt.pcolormesh(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_twc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_RA2M')
    plt.colorbar()

    plt.subplot(413)
    plt.pcolormesh(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_twc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.title('UM_CASIM-100')
    plt.colorbar()

    plt.subplot(414)
    plt.pcolormesh(ifs_data['time'], np.squeeze(ifs_data['height'][0,:]), np.transpose(ifs_data['model_twc'])*1e3,
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlabel('DOY')
    plt.title('ECMWF_IFS')
    plt.colorbar()


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_UM_IFS_CASIM-100_TWCTimeseries_226-257DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_TWCTesting(um_data, ifs_data, misc_data, obs_data, data1, data2, data3, obs, month_flag, missing_files, doy):

    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    um_data['model_iwc'][um_data['model_iwc'] < 0.0] = 0.0
    # ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 1.0] = 0.0
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = 0.0
    misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content - Cloudnet
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']

    ###----------------------------------------------------------------
    ###         Calculate total water content - Model
    ###----------------------------------------------------------------
    data1['qliq'][data1['qliq'] < 0.0] = 0.0
    data1['qice'][data1['qice'] < 0.0] = 0.0
    data1['qtot'] = data1['qliq'] + data1['qice']
    data2['qliq'][data2['qliq'] < 0.0] = 0.0
    data2['qice'][data2['qice'] < 0.0] = 0.0
    data2['qtot'] = data2['qliq'] + data2['qice']
    data3['ql'][data3['ql'] < 0.0] = 0.0
    data3['qi'][data3['qi'] < 0.0] = 0.0
    data3['qtot'] = data3['ql'] + data3['qi']


    ###----------------------------------------------------------------
    ###         New figure - timeseries
    ###----------------------------------------------------------------

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    plt.subplot(211)
    plt.pcolormesh(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_twc'])*1e3,
        # )
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlim([226,257])
    plt.title('TWC [g/m3]; Cloudnet')
    plt.colorbar()

    plt.subplot(212)
    plt.pcolormesh(data1['time'], data1['height'], np.transpose(data1['qtot'])*1e3,
        # )
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlim([226,257])
    plt.title('TWC [g/kg]; Model')
    plt.colorbar()

    plt.show()

    ###----------------------------------------------------------------
    ###         New figure - profiles
    ###----------------------------------------------------------------
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.45, wspace = 0.25)
    ### define axis instance
    ax = plt.gca()

    plt.subplot(331)
    plt.plot(np.nanmean(um_data['model_lwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data1['qliq']*1e3,0), data1['height'], label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: LWC')

    plt.subplot(332)
    plt.plot(np.nanmean(um_data['model_iwc_filtered']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data1['qice']*1e3,0), data1['height'], label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: IWC')

    plt.subplot(333)
    plt.plot(np.nanmean(um_data['model_twc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data1['qtot']*1e3,0), data1['height'], label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: TWC')

    plt.subplot(334)
    plt.plot(np.nanmean(misc_data['model_lwc']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data2['qliq']*1e3,0), data2['height'], label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: LWC')

    plt.subplot(335)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data2['qice']*1e3,0), data2['height'], label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: IWC')

    plt.subplot(336)
    plt.plot(np.nanmean(misc_data['model_twc']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data2['qtot']*1e3,0), data2['height'], label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: TWC')

    plt.subplot(337)
    plt.plot(np.nanmean(ifs_data['model_lwc']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data3['ql']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('ECMWF_IFS: LWC')

    plt.subplot(338)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data3['qi']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    plt.xlim([0.0,0.011])
    plt.title('ECMWF: IWC')

    plt.subplot(339)
    plt.plot(np.nanmean(ifs_data['model_twc']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data3['qtot']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    plt.xlim([0.0,0.12])
    plt.title('ECMWF_IFS: TWC')

    plt.savefig('../FIGS/comparisons/CN-ModelComparison_LWC-IWC-TWC_profiles.png')
    plt.show()

def plot_LWP(um_data, ifs_data, misc_data, obs_data, obs, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

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
    print (obs_data.keys())

    #### set flagged and bad data to nans
    um_data['model_lwp'][um_data['model_lwp'] < 0] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] < 0] = np.nan
    misc_data['model_lwp'][misc_data['model_lwp'] < 0] = np.nan
    # um_data['model_lwp'][um_data['model_lwp'] >= 1000] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] >= 1.0] = np.nan
    # misc_data['model_lwp'][misc_data['model_lwp'] >= 1000] = np.nan
    if obs_switch == 'RADAR':
        obs_data['lwp'][obs_data['lwp'] < 0] = np.nan     ### index 0 is mean
        obs_data['lwp'][obs_data['lwp'] > 0.8] = np.nan    ### >0.8 == >800g/m2
    else:
        obs_data['lwp'][obs_data['lwp'][:,0] < 0, 0] = np.nan     ### index 0 is mean
        obs_data['lwp'][obs_data['lwp'][:,0] > 0.8, 0] = np.nan    ### >0.8 == >800g/m2
        obs_data['lwp'][obs_data['lwp'][:,1] < 0, 1] = np.nan     ### index 1 is min
        obs_data['lwp'][obs_data['lwp'][:,1] > 0.8, 1] = np.nan    ### >0.8 == >800g/m2
        obs_data['lwp'][obs_data['lwp'][:,2] < 0, 2] = np.nan     ### index 2 is max
        obs_data['lwp'][obs_data['lwp'][:,2] > 0.8, 2] = np.nan    ### >0.8 == >800g/m2

    if obs_switch == 'RADAR':
        plt.plot(obs_data['time'][:],obs_data['lwp'][:]*1e3, color = 'purple', label = 'Obs_' + obs_switch + 'grid')
    else:
        plt.plot(obs_data['time'][:],obs_data['lwp'][:,0]*1e3, color = 'purple', label = 'Obs_' + obs_switch + 'grid')
        # plt.plot(obs_data['time'][:],obs_data['lwp'][:,1]*1e3, 'k--')
        # plt.plot(obs_data['time'][:],obs_data['lwp'][:,2]*1e3, 'k--')
        ax.fill_between(obs_data['time'][:], obs_data['lwp'][:,1]*1e3, obs_data['lwp'][:,2]*1e3, color = 'thistle', alpha = 0.5)
    plt.plot(obs['deck7th']['doy'][:],obs['deck7th']['lwp'][:], color = 'k', label = 'Obs_HATPRO')
    plt.plot(um_data['time'][::3],um_data['model_lwp'][::3]*1e3,
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = 'UM_RA2M')
    plt.plot(ifs_data['time'][::3],ifs_data['model_lwp'][::3]*1e3,
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown', label = 'ECMWF_IFS')
    plt.plot(misc_data['time'][::3],misc_data['model_lwp'][::3]*1e3,
        'v', color = 'forestgreen', markeredgecolor = 'darkgreen', label = 'UM_CASIM-100')
    plt.xlabel('Day of Year')
    plt.ylabel('LWP [g/m2]')
    plt.ylim([0,800])
    plt.xlim([doy[0],doy[-1]])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid_UM_IFS_CASIM-100_LWP_226-257DOY_wMissingFiles.svg'
    plt.savefig(fileout)
    plt.show()

def plot_ObsGridComparison(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

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
    print ('Comparing observed Cv based on grid results interpolated on to:')
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
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan

    melt = np.where(um_data['time'] < 240.0)
    freeze = np.where(um_data['time'] >= 240.0)

    plt.subplot(121)
    ax1 = plt.gca()
    ax1.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][melt,:]),0),np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][melt,:]),0),
        np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][melt,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['Cv'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'Obs_IFS')
    ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['Cv'][melt,:]),0) - np.nanstd(np.squeeze(ifs_data['Cv'][melt,:]),0),
        np.nanmean(np.squeeze(ifs_data['Cv'][melt,:]),0) + np.nanstd(np.squeeze(ifs_data['Cv'][melt,:]),0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][melt,:]),0),np.nanmean(np.squeeze(obs_data['height'][melt,:]),0), 'k--', linewidth = 3, label = 'Obs_UM')

    plt.xlabel('Cloud Fraction')
    plt.ylabel('Height [m]')
    plt.title('Melt')
    plt.ylim([0,10000])
    plt.xlim([0,1])
    plt.legend()

    plt.subplot(122)
    ax2 = plt.gca()
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][freeze,:]),0),
        np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][freeze,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['Cv'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'Obs_IFS')
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['Cv'][freeze,:]),0) - np.nanstd(np.squeeze(ifs_data['Cv'][freeze,:]),0),
        np.nanmean(np.squeeze(ifs_data['Cv'][freeze,:]),0) + np.nanstd(np.squeeze(ifs_data['Cv'][freeze,:]),0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][freeze,:]),0),np.nanmean(np.squeeze(obs_data['height'][freeze,:]),0), 'k--', linewidth = 3, label = 'Obs_UM')
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
        fileout = 'FIGS/Obs_GridComparison_Cv_splitSeason_226-257DOY_wMissingFiles.svg'
    # plt.savefig(fileout)
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
        fileout = 'FIGS/Obs-UMgrid_UM_IFS_CASIM-100_Cv_splitSeason_226-257DOY_wMissingFiles.svg'
    # plt.savefig(fileout)
    plt.show()

def plot_lwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

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

    ### make bespoke obs indices
    meltobs = np.where(obs_data['time'] < 240.0)
    freezeobs = np.where(obs_data['time'] >= 240.0)

    if obs_switch == 'RADAR':
        plt.subplot(121)
        ax1 = plt.gca()
        plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3,obs_data['height'], 'k--', linewidth = 3, label = 'Obs')
        ax1.fill_betweenx(obs_data['height'],np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3,
            np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
        ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][melt,:])*1e3,0),
            np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
        plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
        ax1.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,
            np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
        plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][melt,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
        ax1.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][melt,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][melt,:])*1e3,0),
            np.nanmean(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][melt,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

        plt.xlabel('Liquid water content [g/m3]')
        plt.ylabel('Height [m]')
        plt.title('Melt')
        plt.ylim([0,10000])
        plt.xlim([0,0.4])
        plt.legend()

        plt.subplot(122)
        ax2 = plt.gca()
        plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3,obs_data['height'], 'k--', linewidth = 3, label = 'Obs')
        ax2.fill_betweenx(obs_data['height'],np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3,
            np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
        ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,
            np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
        plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
        ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,
            np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
        plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
        ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][freeze,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][freeze,:])*1e3,0),
            np.nanmean(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][freeze,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)

        plt.xlabel('Liquid water content [g/m3]')
        plt.title('Freeze up')
        plt.yticks([])
        plt.ylim([0,10000])
        plt.xlim([0,0.4])
        # plt.legend()

        print ('******')
        print ('')
        print ('Finished plotting! :)')
        print ('')

        if month_flag == -1:
            fileout = 'FIGS/Obs-RADARgrid_UM_IFS_CASIM-100_LWC_splitSeason_wMissingFiles.svg'
        # plt.savefig(fileout, dpi=300)
        plt.show()

    else:
        plt.subplot(121)
        ax1 = plt.gca()
        plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][meltobs,:]),0), 'k--', linewidth = 3, label = 'Obs')
        ax1.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][meltobs,:]),0),np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][melt,:]),0)*1e3,
            np.nanmean(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][meltobs,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][melt,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
        ax1.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][melt,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][melt,:])*1e3,0),
            np.nanmean(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][melt,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
        plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][melt,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][melt,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
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
        plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][freezeobs,:]),0), 'k--', linewidth = 3, label = 'Obs')
        ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][freezeobs,:]),0),np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc'][freeze,:]),0)*1e3,
            np.nanmean(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc'][freezeobs,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
        ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3,
            np.nanmean(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
        plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
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
            fileout = 'FIGS/Obs-UMgrid_UM_IFS_CASIM-100_LWC_splitSeason_wMissingFiles.svg'
        # plt.savefig(fileout, dpi=300)
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
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][freeze,:]),0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][freeze,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][freeze,:]),0)*1e3, color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][freeze,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][freeze,:]),0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
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
        fileout = 'FIGS/Obs-UMgrid_UM_IFS_CASIM-100_IWC_splitSeason_wMissingFiles.svg'
    # plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined 1d timeseries:')
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
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
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

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    # datenums_obs['foremast'] = obs['foremast'].variables['time'][:] ### obs['foremast'] data on different timestep
    # time_obs['foremast'] = calcTime_Mat2DOY(datenums_obs['foremast'])

    # datenums_obs['deck7th'] = obs['deck7th'].variables['time'][:] ### 7th deck data on different timestep
    # time_obs['deck7th'] = calcTime_Mat2DOY(datenums_obs['deck7th'])

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/',out_dir4 == 'ECMWF_IFS/'):
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
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)
    # ax = plt.gca()
    # plt.plot(data1['time'], data1['sfc_pressure'].data/1e2, color = 'steelblue', label = label1)
    # plt.plot(data2['time'], data2['sfc_pressure'].data/1e2, color = 'forestgreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'darkorange', label = label3)
    # else:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'darkorange',label = label3)
    # plt.plot(obs['foremast'].variables['doy'][:],obs['foremast'].variables['Psurf'][:], 'k')
    # plt.title('sfc_pressure [hPa]')
    # plt.ylim([980, 1010])
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'steelblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'forestgreen')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    # data1['surface_net_SW_radiation'].data[data1['surface_net_SW_radiation'].data == 0] = np.nan
    # data2['surface_net_SW_radiation'].data[data2['surface_net_SW_radiation'].data == 0] = np.nan
    # if out_dir4 != 'OUT_25H/': data3['surface_net_SW_radiation'].data[data3['surface_net_SW_radiation'].data == 0] = np.nan

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'forestgreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    # plt.subplot(3,3,6)
    # ax = plt.gca()
    # plt.plot(time_um, data['snowfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_snow'].data)
    # plt.title('sfc_snow_amount [kg/m2]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,7)
    # ax = plt.gca()
    # plt.plot(time_um, data['rainfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_rain'].data)
    # plt.title('sfc_rain_amount [kg/m2]')
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station_fluxes']['mday'],obs['ice_station_fluxes']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/244-256_oden_metum_ifs_casim-aeroprof_TSa.svg'
    # plt.savefig(fileout, dpi=300)
    plt.show()

def plot_scaledBLCv_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Scaling cloudnet data by identified thetaE inversions:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # #################################################################
    # ## save data into temp variables to allow subsampling
    # #################################################################
    # bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    # bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    # if ifs_flag == True:
    #     bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    # else:
    #     bldepth3 = data3['bl_depth'][data3['hrly_flag']]

    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    # #### ---------------------------------------------------------------
    # #### ONLY LOOK AT SONDES FROM THE DRIFT
    # #### ---------------------------------------------------------------
    # drift = np.where(np.logical_and(obs['inversions']['thetaE']['time'] >= 225.9, obs['inversions']['thetaE']['time'] <= 258.0))

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = obs['inversions']['thetaE']['invbase']
    obsmlh = obs['inversions']['thetaE']['sfmlheight']

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    temp_time = obs['inversions']['thetaE']['time']         #### temporary time array for indexing
    temp_time2 = np.zeros(len(temp_time))
    temp_time2[:] = np.nan
    temp_time2[:nanindices[0]] = temp_time[:nanindices[0]]
    temp_time2[nanindices[3]+1:nanindices[4]] = temp_time[nanindices[3]+1:nanindices[4]]
    temp_time2[nanindices[7]+1:nanindices[8]] = temp_time[nanindices[7]+1:nanindices[8]]
    temp_inv = np.zeros(len(obsinv))
    temp_inv[:] = np.nan
    temp_inv[:nanindices[0]] = obsinv[:nanindices[0]]
    temp_inv[nanindices[3]+1:nanindices[4]] = obsinv[nanindices[3]+1:nanindices[4]]
    temp_inv[nanindices[7]+1:nanindices[8]] = obsinv[nanindices[7]+1:nanindices[8]]
    temp_sfml = np.zeros(len(obsmlh))
    temp_sfml[:] = np.nan
    temp_sfml[:nanindices[0]] = obsmlh[:nanindices[0]]
    temp_sfml[nanindices[3]+1:nanindices[4]] = obsmlh[nanindices[3]+1:nanindices[4]]
    temp_sfml[nanindices[7]+1:nanindices[8]] = obsmlh[nanindices[7]+1:nanindices[8]]

    ### reset time array with new one accounting for missing FILES
    obs['inversions']['thetaE']['time'] = temp_time2
    obsinv = temp_inv
    obsmlh = temp_sfml

    ### save non-nan values to dictionary
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['thetaE']['time'][~np.isnan(obs['inversions']['thetaE']['time'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data['Cv_6hrly'] = obs_data['Cv'][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); obssfml[:] = np.nan

    ### look for altitudes < invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    plt.figure()
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),obs['inversions']['thetaE']['invbase'])
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaledCv'] = {}
    obs['inversions']['scaledCv']['binned'] = {}
    obs['inversions']['scaledCv']['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['mean'][:] = np.nan
    obs['inversions']['scaledCv']['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['blCv'] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['blCv'][:] = np.nan

    ### fill arrays with cloudnet data below each invbase
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        ###         +1 includes invbase in array
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['blCv'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)] = obs_data['Cv_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['blCv'][i,tempvar]
            if np.size(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaledCv']['mean'][i,k] = np.nanmean(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaledCv']['stdev'][i,k] = np.nanstd(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    ####        data from thetaE algorithm is on radiosonde (6 hourly) timesteps already
    #### ---------------------------------------------------------------

    ### need to make sure times where we don't have an inversion (invbase == nan) in the IFS doesn't mess up the algorithm
    ###     set these cases to zero for now
    data3['inversions']['invbase'][np.isnan(data3['inversions']['invbase'])] = 0.0

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    ####        time
    tim1 = np.zeros(len(data1['inversions']['time']))
    tim1[:] = np.nan
    tim1[:nanindices[0]] = data1['inversions']['time'][:nanindices[0]]
    tim1[nanindices[3]+1:nanindices[4]] = data1['inversions']['time'][nanindices[3]+1:nanindices[4]]
    tim1[nanindices[7]+1:nanindices[8]] = data1['inversions']['time'][nanindices[7]+1:nanindices[8]]
    tim2 = tim1
    tim3 = tim1
    ####        inversions
    inv1 = np.zeros(len(data1['inversions']['invbase']))
    inv1[:] = np.nan
    inv1[:nanindices[0]] = data1['inversions']['invbase'][:nanindices[0]]
    inv1[nanindices[3]+1:nanindices[4]] = data1['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv1[nanindices[7]+1:nanindices[8]] = data1['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv2 = np.zeros(len(data2['inversions']['invbase']))
    inv2[:] = np.nan
    inv2[:nanindices[0]] = data2['inversions']['invbase'][:nanindices[0]]
    inv2[nanindices[3]+1:nanindices[4]] = data2['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv2[nanindices[7]+1:nanindices[8]] = data2['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv3 = np.zeros(len(data3['inversions']['invbase']))
    inv3[:] = np.nan
    inv3[:nanindices[0]] = data3['inversions']['invbase'][:nanindices[0]]
    inv3[nanindices[3]+1:nanindices[4]] = data3['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv3[nanindices[7]+1:nanindices[8]] = data3['inversions']['invbase'][nanindices[7]+1:nanindices[8]]

    ### save non-nan values
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    tim1 = tim1[~np.isnan(tim1)]
    tim2 = tim2[~np.isnan(tim2)]
    tim3 = tim3[~np.isnan(tim3)]
    inv1 = inv1[~np.isnan(inv1)]
    inv2 = inv2[~np.isnan(inv2)]
    inv3 = inv3[~np.isnan(inv3)]

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    # #### ---------------------------------------------------------------
    # #### remove flagged IFS heights
    # #### ---------------------------------------------------------------
    # data3['height'][data3['height'] == -9999] = 0.0
    #         #### set all heights to zero if flagged. setting to nan caused problems further on
    # data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    # #### ---------------------------------------------------------------
    # #### Define meteorological periods from Jutta's paper
    # #### ---------------------------------------------------------------
    #
    # ## Meteorological period definitions from Jutta's paper:
    # ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    # ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    # ##     radar measurements were not possible during heavy ice breaking because of excessive
    # ##     vibration, cloud characteristics and fog heights are not available during period 2.
    # ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    # ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    # ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    # ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    # ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    # ##     out to the ice edge. "
    # #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...
    #
    # ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    # p3 = np.where(um_data['time'] < 230.0)
    # p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    # p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    # p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    # p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    # p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time'].data[::6]
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time'].data[::6]
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time'].data[::6]

    #### define empty arrays of nans to fill with scaled data
    data1['scaledCv'] = {}
    data1['scaledCv']['binned'] = {}
    data1['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaledCv']['mean'][:] = np.nan
    data1['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaledCv']['stdev'][:] = np.nan
    data1['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(um_data['height'],1)]); data1['blCv'][:] = np.nan

    data2['scaledCv'] = {}
    data2['scaledCv']['binned'] = {}
    data2['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaledCv']['mean'][:] = np.nan
    data2['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaledCv']['stdev'][:] = np.nan
    data2['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(misc_data['height'],1)]); data2['blCv'][:] = np.nan

    data3['scaledCv'] = {}
    data3['scaledCv']['binned'] = {}
    data3['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaledCv']['mean'][:] = np.nan
    data3['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaledCv']['stdev'][:] = np.nan
    data3['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(ifs_data['height'],1)]); data3['blCv'][:] = np.nan

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data1['scaledCv']['inversion_Tindex'][:] = np.nan
    data1['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data1['scaledCv']['inversionForCloudnet'][:] = np.nan
    data2['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data2['scaledCv']['inversion_Tindex'][:] = np.nan
    data2['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data2['scaledCv']['inversionForCloudnet'][:] = np.nan
    data3['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data3['scaledCv']['inversion_Tindex'][:] = np.nan
    data3['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data3['scaledCv']['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.size(np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))[0][0]
        if np.size(np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))[0][0]
        if np.size(np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data1['scaledCv']['inversionForCloudnet'][int(data1['scaledCv']['inversion_Tindex'][i])] = inv1[i]
        if data2['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data2['scaledCv']['inversionForCloudnet'][int(data2['scaledCv']['inversion_Tindex'][i])] = inv2[i]
        if data3['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data3['scaledCv']['inversionForCloudnet'][int(data3['scaledCv']['inversion_Tindex'][i])] = inv3[i]
    #
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    #
    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(data1['scaledTime'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(data2['scaledTime'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(data3['scaledTime'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(data1['scaledTime'])):        ### all can go in this loop, data1['scaledTime'] == 6-hourly data

        ### main inversion base assignments
        ###         (1) find where UM height array matches the invbase index for index i
        ###         (2) find where UM height array matches the invbase index for index i
        ###         (3) find where IFS height array is less than or equal to the UM-gridded invbase index for index i
        if np.size(np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(misc_data['height'][i,:].data == data2['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(misc_data['height'][i,:].data == data2['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(ifs_data['height'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            temp = ifs_data['height'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    plt.figure()
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(data1['scaledTime'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1,inv1)
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(data2['scaledTime'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2, inv2)
    plt.ylim([0,3e3])
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(data3['scaledTime'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3, inv3)
    plt.ylim([0,3e3])
    plt.show()

    ### set 6 hourly cloudnet Cv arrays as tempvars
    ra2m_Cv = um_data['model_Cv_filtered'][::6,:]
    casim_Cv = misc_data['model_Cv_filtered'][::6,:]
    ifs_Cv = ifs_data['model_snow_Cv_filtered'][::6,:]

    ### find all Cv data below identified inversion
    for i in range(0,np.size(data1['scaledTime'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaledCv']['binned']['t' + str(i)] = {}
        data2['scaledCv']['binned']['t' + str(i)] = {}
        data3['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['blCv'][i,:int(data1['inversions']['invbase_kIndex'][i]+1)] = ra2m_Cv[i,:int(data1['inversions']['invbase_kIndex'][i]+1)]
        data2['blCv'][i,:int(data2['inversions']['invbase_kIndex'][i]+1)] = casim_Cv[i,:int(data2['inversions']['invbase_kIndex'][i]+1)]
        data3['blCv'][i,:int(data3['inversions']['invbase_kIndex'][i]+1)] = ifs_Cv[i,:int(data3['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            data1['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data1['blCv'][i,tempvar1]
            if np.size(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaledCv']['mean'][i,k] = np.nanmean(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data1['scaledCv']['stdev'][i,k] = np.nanstd(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data2['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data2['blCv'][i,tempvar2]
            if np.size(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaledCv']['mean'][i,k] = np.nanmean(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data2['scaledCv']['stdev'][i,k] = np.nanstd(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data3['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data3['blCv'][i,tempvar3]
            if np.size(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaledCv']['mean'][i,k] = np.nanmean(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data3['scaledCv']['stdev'][i,k] = np.nanstd(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    ### save working data for debug
    np.save('working_data1', data1)

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.show()

    ### obs
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['blCv'])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),np.squeeze(obs['inversions']['thetaE']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_ra2m
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],um_data['height'][0,:].data,np.transpose(data1['blCv'])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['time'],data1['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_casim-100
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],misc_data['height'][0,:].data,np.transpose(data2['blCv'])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['time'],data2['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### ecmwf_ifs
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],ifs_data['height'][0,:].data,np.transpose(data3['blCv'])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['time'],data3['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### profiles
    ###         loop through and set all zeros to nans
    obsmean = obs['inversions']['scaledCv']['mean']
    ra2mmean = data1['scaledCv']['mean']
    casimmean = data2['scaledCv']['mean']
    ifsmean = data3['scaledCv']['mean']
    # for i in range(0, len(data1['scaledTime'])):
    #     obsmean[i,obs['inversions']['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     ra2mmean[i,data1['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     casimmean[i,data2['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     ifsmean[i,data3['scaledCv']['mean'][i,:] == 0.0] = np.nan
    plt.plot(np.nanmean(obsmean,0),obs['inversions']['scaledZ'], '--', color = 'k', linewidth = 2, label = 'Obs')
    plt.plot(np.nanmean(ra2mmean,0),data1['scaledZ'], '^-', color = 'steelblue', linewidth = 2, label = label1)
    plt.plot(np.nanmean(casimmean,0),data2['scaledZ'], 'v-', color = 'forestgreen', linewidth = 2, label = label2)
    plt.plot(np.nanmean(ifsmean,0),data3['scaledZ'], 'd-', color = 'darkorange', linewidth = 2, label = label3)
    plt.legend()
    plt.show()

def plot_scaledBLCv_JVInv(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Plotting BL depth timeseries:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ## save data into temp variables to allow subsampling
    #################################################################
    bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    if ifs_flag == True:
        bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    else:
        bldepth3 = data3['bl_depth'][data3['hrly_flag']]


    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    #################################################################
    ## convert model inversion timesteps
    #################################################################
    data1['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data1['inversions']['mday']))
    data2['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data2['inversions']['mday']))
    data3['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data3['inversions']['mday']))

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['inversions']['doy'] >= 225.9, obs['inversions']['doy'] <= 258.0))

    ### save in dict for ease
    obs['inversions']['doy_drift'] = obs['inversions']['doy'][drift]

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = np.squeeze(obs['inversions']['invbase'][drift])
    obsmlh = np.squeeze(obs['inversions']['sfmlheight'][drift])

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])
    obs['inversions']['doy_drift'][nanindices] = np.nan
    obsinv[nanindices] = np.nan
    obsmlh[nanindices] = np.nan

    ### save non-nan values to dictionary
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['doy_drift'][~np.isnan(obs['inversions']['doy_drift'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data['Cv_6hrly'] = obs_data['Cv'][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); obssfml[:] = np.nan

    ### look for altitudes < invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    plt.figure()
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']))
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaledCv'] = {}
    obs['inversions']['scaledCv']['binned'] = {}
    obs['inversions']['scaledCv']['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['mean'][:] = np.nan
    obs['inversions']['scaledCv']['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['blCv'] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['blCv'][:] = np.nan

    ###
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['blCv'][i,:int(obs['inversions']['invbase_kIndex'][i])] = obs_data['Cv_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i])]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['blCv'][i,tempvar]
            if np.size(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaledCv']['mean'][i,k] = np.nanmean(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaledCv']['stdev'][i,k] = np.nanstd(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    #### ---------------------------------------------------------------
    #### ---------------------------------------------------------------
    #### ONLY LOOK AT DATA FROM THE DRIFT
    #### ---------------------------------------------------------------
    driftmod = np.where(np.logical_and(data1['inversions']['doy'] >= 226.0, data1['inversions']['doy'] <= 258.0))
    data1['inversions']['doy_drift'] = data1['inversions']['doy'][driftmod[0][1:-2]]
    data1['inversions']['invbase_drift'] = data1['inversions']['invbase'][driftmod[0][1:-2]]
    data2['inversions']['doy_drift'] = data2['inversions']['doy'][driftmod[0][1:-2]]
    data2['inversions']['invbase_drift'] = data2['inversions']['invbase'][driftmod[0][1:-2]]
    data3['inversions']['doy_drift'] = data3['inversions']['doy'][driftmod[0][1:-2]]
    data3['inversions']['invbase_drift'] = data3['inversions']['invbase'][driftmod[0][1:-2]]

    #### make inversion tempvars to allow for easy subsampling
    tim1 = np.squeeze(data1['inversions']['doy_drift'][data1['hrly_flag']])
    tim2 = np.squeeze(data2['inversions']['doy_drift'][data2['hrly_flag']])
    tim3 = np.squeeze(data3['inversions']['doy_drift'][data3['hrly_flag']])
    inv1 = np.squeeze(data1['inversions']['invbase_drift'][data1['hrly_flag'],0])
    inv2 = np.squeeze(data2['inversions']['invbase_drift'][data2['hrly_flag'],0])
    inv3 = np.squeeze(data3['inversions']['invbase_drift'][data3['hrly_flag'],0])

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    #### ---------------------------------------------------------------
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems further on
    data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### Define meteorological periods from Jutta's paper
    #### ---------------------------------------------------------------

    ## Meteorological period definitions from Jutta's paper:
    ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    ##     radar measurements were not possible during heavy ice breaking because of excessive
    ##     vibration, cloud characteristics and fog heights are not available during period 2.
    ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    ##     out to the ice edge. "
    #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...

    ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    p3 = np.where(um_data['time'] < 230.0)
    p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    #### define empty arrays of nans to fill with scaled data
    data1['scaledCv'] = {}
    data1['scaledCv']['binned'] = {}
    data1['scaledCv']['mean'] = np.zeros([np.size(um_data['height'],0),len(Zpts)]); data1['scaledCv']['mean'][:] = np.nan
    data1['scaledCv']['stdev'] = np.zeros([np.size(um_data['height'],0),len(Zpts)]); data1['scaledCv']['stdev'][:] = np.nan
    data1['blCv'] = np.zeros([np.size(um_data['height'],0),np.size(um_data['height'],1)]); data1['blCv'][:] = np.nan

    data2['scaledCv'] = {}
    data2['scaledCv']['binned'] = {}
    data2['scaledCv']['mean'] = np.zeros([np.size(misc_data['height'],0),len(Zpts)]); data2['scaledCv']['mean'][:] = np.nan
    data2['scaledCv']['stdev'] = np.zeros([np.size(misc_data['height'],0),len(Zpts)]); data2['scaledCv']['stdev'][:] = np.nan
    data2['blCv'] = np.zeros([np.size(misc_data['height'],0),np.size(misc_data['height'],1)]); data2['blCv'][:] = np.nan

    data3['scaledCv'] = {}
    data3['scaledCv']['binned'] = {}
    data3['scaledCv']['mean'] = np.zeros([np.size(ifs_data['height'],0),len(Zpts)]); data3['scaledCv']['mean'][:] = np.nan
    data3['scaledCv']['stdev'] = np.zeros([np.size(ifs_data['height'],0),len(Zpts)]); data3['scaledCv']['stdev'][:] = np.nan
    data3['blCv'] = np.zeros([np.size(ifs_data['height'],0),np.size(ifs_data['height'],1)]); data3['blCv'][:] = np.nan

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time']
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time']
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time']

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim1)); data1['scaledCv']['inversion_Tindex'][:] = np.nan
    data1['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim1)); data1['scaledCv']['inversionForCloudnet'][:] = np.nan
    data2['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim2)); data2['scaledCv']['inversion_Tindex'][:] = np.nan
    data2['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim2)); data2['scaledCv']['inversionForCloudnet'][:] = np.nan
    data3['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim3)); data3['scaledCv']['inversion_Tindex'][:] = np.nan
    data3['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim3)); data3['scaledCv']['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.floor(tim1[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv1[i] = np.nan
            continue
        elif np.size(np.where(np.round(um_data['time'].data,3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaledCv']['inversion_Tindex'][i] = np.where(np.round(um_data['time'].data,3) == np.round(tim1[i],3))[0][0]
        if np.floor(tim2[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv2[i] = np.nan
            continue
        elif np.size(np.where(np.round(misc_data['time'].data,3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaledCv']['inversion_Tindex'][i] = np.where(np.round(misc_data['time'].data,3) == np.round(tim2[i],3))[0][0]
        if np.floor(tim3[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv3[i] = np.nan
            continue
        elif np.size(np.where(np.round(ifs_data['time'].data,3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaledCv']['inversion_Tindex'][i] = np.where(np.round(ifs_data['time'].data,3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data1['scaledCv']['inversionForCloudnet'][int(data1['scaledCv']['inversion_Tindex'][i])] = inv1[i]
        if data2['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data2['scaledCv']['inversionForCloudnet'][int(data2['scaledCv']['inversion_Tindex'][i])] = inv2[i]
        if data3['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data3['scaledCv']['inversionForCloudnet'][int(data3['scaledCv']['inversion_Tindex'][i])] = inv3[i]

    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)

    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(um_data['time'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(misc_data['time'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(ifs_data['time'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(um_data['time'])):        ### all can go in this loop, um_data['time'] == hourly data

        ### main inversion base assignments
        if np.size(np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(data2['height'][1:].data == data2['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(data2['height'][1:].data == data2['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(data3['height_hrly'][i].data <= data3['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            temp = data3['height_hrly'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    plt.figure()
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(um_data['time'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1[:-22],inv1[:-22])
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(misc_data['time'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2[:-22], inv2[:-22])
    plt.ylim([0,3e3])
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(ifs_data['time'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3[:-22], inv3[:-22])
    plt.ylim([0,3e3])
    plt.show()

    # print (zind3)
    #### re-check inversion algorithm success rate to make sure no !0 values got dropped
    zzind1 = np.where(data1['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind1rate = np.size(zzind1) / np.float(np.size(inv1)) * 100.0
    zzind2 = np.where(data2['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind2rate = np.size(zzind2) / np.float(np.size(inv2)) * 100.0
    zzind3 = np.where(data3['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind3rate = np.size(zzind3) / np.float(np.size(inv3)) * 100.0

    ### try i = 0 first to see if it works
    ### this will go into a loop once tested
    # i = 110
    for i in range(0,np.size(um_data['time'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaledCv']['binned']['t' + str(i)] = {}
        data2['scaledCv']['binned']['t' + str(i)] = {}
        data3['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['blCv'][i,:int(data1['inversions']['invbase_kIndex'][i])] = um_data['model_Cv_filtered'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        data2['blCv'][i,:int(data2['inversions']['invbase_kIndex'][i])] = misc_data['model_Cv_filtered'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        data3['blCv'][i,:int(data3['inversions']['invbase_kIndex'][i])] = ifs_data['model_snow_Cv_filtered'][i,:int(data3['inversions']['invbase_kIndex'][i])]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            data1['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data1['blCv'][i,tempvar1]
            if np.size(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaledCv']['mean'][i,k] = np.nanmean(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data1['scaledCv']['stdev'][i,k] = np.nanstd(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data2['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data2['blCv'][i,tempvar2]
            if np.size(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaledCv']['mean'][i,k] = np.nanmean(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data2['scaledCv']['stdev'][i,k] = np.nanstd(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data3['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data3['blCv'][i,tempvar3]
            if np.size(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaledCv']['mean'][i,k] = np.nanmean(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data3['scaledCv']['stdev'][i,k] = np.nanstd(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.show()

    ### obs
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['blCv'])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['doy_drift']),np.squeeze(obs['inversions']['invbase'][drift]),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_ra2m
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(um_data['time'].data,um_data['height'][0,:].data,np.transpose(data1['blCv'])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['doy_drift'],np.squeeze(data1['inversions']['invbase_drift']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_casim-100
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(misc_data['time'].data,misc_data['height'][0,:].data,np.transpose(data2['blCv'])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['doy'],np.squeeze(data2['inversions']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### ecmwf_ifs
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(ifs_data['time'].data,ifs_data['height'][0,:].data,np.transpose(data3['blCv'])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['doy'],np.squeeze(data3['inversions']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### profiles
    plt.plot(np.nanmean(obs['inversions']['scaledCv']['mean'],0),obs['inversions']['scaledZ'], color = 'k', linewidth = 2, label = 'Obs')
    plt.plot(np.nanmean(data1['scaledCv']['mean'][::6,:],0),data1['scaledZ'], color = 'steelblue', linewidth = 2, label = label1)
    plt.plot(np.nanmean(data2['scaledCv']['mean'][::6,:],0),data2['scaledZ'], color = 'forestgreen', linewidth = 2, label = label2)
    plt.plot(np.nanmean(data3['scaledCv']['mean'][::6,:],0),data3['scaledZ'], color = 'darkorange', linewidth = 2, label = label3)
    plt.legend()
    plt.show()

def plot_scaledBL_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3, var):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Scaling cloudnet data by identified thetaE inversions:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ### interpolate obs cloudnet data for continuous array
    #################################################################
    obs_data = interpCloudnet(obs_data, month_flag, missing_files, doy)

    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    if var == 'Cv':
        obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
        um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
        ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
        misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan
    elif var == 'lwc':
        obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
        obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
        um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
        ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
        ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
        misc_data['model_lwc'][misc_data['model_lwc'] <= 0.0] = np.nan
        #### change units to g/cm3
        obs_data['lwc'] = obs_data['lwc'] * 1e3
        um_data['model_lwc'] = um_data['model_lwc'] * 1e3
        misc_data['model_lwc'] = misc_data['model_lwc'] * 1e3
        ifs_data['model_lwc'] = ifs_data['model_lwc'] * 1e3
    elif var == 'iwc':
        obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
        obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
        um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
        ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] <= 0.0] = np.nan
        # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
        misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] <= 0.0] = np.nan
        #### change units to g/cm3
        obs_data['iwc'] = obs_data['iwc'] * 1e3
        um_data['model_iwc_filtered'] = um_data['model_iwc_filtered'] * 1e3
        misc_data['model_iwc_filtered'] = misc_data['model_iwc_filtered'] * 1e3
        ifs_data['model_snow_iwc_filtered'] = ifs_data['model_snow_iwc_filtered'] * 1e3

    # #### ---------------------------------------------------------------
    # #### ONLY LOOK AT SONDES FROM THE DRIFT
    # #### ---------------------------------------------------------------
    # drift = np.where(np.logical_and(obs['inversions']['thetaE']['time'] >= 225.9, obs['inversions']['thetaE']['time'] <= 258.0))

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = obs['inversions']['thetaE']['invbase']
    obsmlh = obs['inversions']['thetaE']['sfmlheight']

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    temp_time = obs['inversions']['thetaE']['time']         #### temporary time array for indexing
    temp_time2 = np.zeros(len(temp_time))
    temp_time2[:] = np.nan
    temp_time2[:nanindices[0]] = temp_time[:nanindices[0]]
    temp_time2[nanindices[3]+1:nanindices[4]] = temp_time[nanindices[3]+1:nanindices[4]]
    temp_time2[nanindices[7]+1:nanindices[8]] = temp_time[nanindices[7]+1:nanindices[8]]
    temp_inv = np.zeros(len(obsinv))
    temp_inv[:] = np.nan
    temp_inv[:nanindices[0]] = obsinv[:nanindices[0]]
    temp_inv[nanindices[3]+1:nanindices[4]] = obsinv[nanindices[3]+1:nanindices[4]]
    temp_inv[nanindices[7]+1:nanindices[8]] = obsinv[nanindices[7]+1:nanindices[8]]
    temp_sfml = np.zeros(len(obsmlh))
    temp_sfml[:] = np.nan
    temp_sfml[:nanindices[0]] = obsmlh[:nanindices[0]]
    temp_sfml[nanindices[3]+1:nanindices[4]] = obsmlh[nanindices[3]+1:nanindices[4]]
    temp_sfml[nanindices[7]+1:nanindices[8]] = obsmlh[nanindices[7]+1:nanindices[8]]

    ### reset time array with new one accounting for missing FILES
    obs['inversions']['thetaE']['time'] = temp_time2
    obsinv = temp_inv
    obsmlh = temp_sfml

    ### save non-nan values to dictionary
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['thetaE']['time'][~np.isnan(obs['inversions']['thetaE']['time'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data[var + '_6hrly'] = obs_data[var][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); obssfml[:] = np.nan

    ### look for altitudes <= invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    fig = plt.figure(figsize=(8,6))
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),obs['inversions']['thetaE']['invbase'])
    plt.xlabel('DOY')
    plt.ylabel('Z [m]')
    plt.savefig('FIGS/' + var + '_obs_inversionDetection_timeseries.png')
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaled' + var] = {}
    obs['inversions']['scaled' + var]['binned'] = {}
    obs['inversions']['scaled' + var]['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaled' + var]['mean'][:] = np.nan
    obs['inversions']['scaled' + var]['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaled' + var]['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['bl' + var] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['bl' + var][:] = np.nan

    ### fill arrays with cloudnet data below each invbase
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaled' + var]['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        ###         +1 includes invbase in array
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['bl' + var][i,:int(obs['inversions']['invbase_kIndex'][i]+1)] = obs_data[var + '_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['bl' + var][i,tempvar]
            if np.size(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaled' + var]['mean'][i,k] = np.nanmean(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaled' + var]['stdev'][i,k] = np.nanstd(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    ####        data from thetaE algorithm is on radiosonde (6 hourly) timesteps already
    #### ---------------------------------------------------------------

    ### need to make sure times where we don't have an inversion (invbase == nan) in the IFS doesn't mess up the algorithm
    ###     set these cases to zero for now
    data3['inversions']['invbase'][np.isnan(data3['inversions']['invbase'])] = 0.0

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    ####        time
    tim1 = np.zeros(len(data1['inversions']['time']))
    tim1[:] = np.nan
    tim1[:nanindices[0]] = data1['inversions']['time'][:nanindices[0]]
    tim1[nanindices[3]+1:nanindices[4]] = data1['inversions']['time'][nanindices[3]+1:nanindices[4]]
    tim1[nanindices[7]+1:nanindices[8]] = data1['inversions']['time'][nanindices[7]+1:nanindices[8]]
    tim2 = tim1
    tim3 = tim1
    ####        inversions
    inv1 = np.zeros(len(data1['inversions']['invbase']))
    inv1[:] = np.nan
    inv1[:nanindices[0]] = data1['inversions']['invbase'][:nanindices[0]]
    inv1[nanindices[3]+1:nanindices[4]] = data1['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv1[nanindices[7]+1:nanindices[8]] = data1['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv2 = np.zeros(len(data2['inversions']['invbase']))
    inv2[:] = np.nan
    inv2[:nanindices[0]] = data2['inversions']['invbase'][:nanindices[0]]
    inv2[nanindices[3]+1:nanindices[4]] = data2['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv2[nanindices[7]+1:nanindices[8]] = data2['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv3 = np.zeros(len(data3['inversions']['invbase']))
    inv3[:] = np.nan
    inv3[:nanindices[0]] = data3['inversions']['invbase'][:nanindices[0]]
    inv3[nanindices[3]+1:nanindices[4]] = data3['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv3[nanindices[7]+1:nanindices[8]] = data3['inversions']['invbase'][nanindices[7]+1:nanindices[8]]

    ### save non-nan values
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    tim1 = tim1[~np.isnan(tim1)]
    tim2 = tim2[~np.isnan(tim2)]
    tim3 = tim3[~np.isnan(tim3)]
    inv1 = inv1[~np.isnan(inv1)]
    inv2 = inv2[~np.isnan(inv2)]
    inv3 = inv3[~np.isnan(inv3)]

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    # #### ---------------------------------------------------------------
    # #### remove flagged IFS heights
    # #### ---------------------------------------------------------------
    # data3['height'][data3['height'] == -9999] = 0.0
    #         #### set all heights to zero if flagged. setting to nan caused problems further on
    # data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    # #### ---------------------------------------------------------------
    # #### Define meteorological periods from Jutta's paper
    # #### ---------------------------------------------------------------
    #
    # ## Meteorological period definitions from Jutta's paper:
    # ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    # ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    # ##     radar measurements were not possible during heavy ice breaking because of excessive
    # ##     vibration, cloud characteristics and fog heights are not available during period 2.
    # ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    # ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    # ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    # ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    # ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    # ##     out to the ice edge. "
    # #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...
    #
    # ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    # p3 = np.where(um_data['time'] < 230.0)
    # p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    # p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    # p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    # p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    # p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time'].data[::6]
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time'].data[::6]
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time'].data[::6]

    #### define empty arrays of nans to fill with scaled data
    data1['scaled' + var] = {}
    data1['scaled' + var]['binned'] = {}
    data1['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaled' + var]['mean'][:] = np.nan
    data1['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaled' + var]['stdev'][:] = np.nan
    data1['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(um_data['height'],1)]); data1['bl' + var][:] = np.nan

    data2['scaled' + var] = {}
    data2['scaled' + var]['binned'] = {}
    data2['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaled' + var]['mean'][:] = np.nan
    data2['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaled' + var]['stdev'][:] = np.nan
    data2['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(misc_data['height'],1)]); data2['bl' + var][:] = np.nan

    data3['scaled' + var] = {}
    data3['scaled' + var]['binned'] = {}
    data3['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaled' + var]['mean'][:] = np.nan
    data3['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaled' + var]['stdev'][:] = np.nan
    data3['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(ifs_data['height'],1)]); data3['bl' + var][:] = np.nan

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data1['scaled' + var]['inversion_Tindex'][:] = np.nan
    data1['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data1['scaled' + var]['inversionForCloudnet'][:] = np.nan
    data2['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data2['scaled' + var]['inversion_Tindex'][:] = np.nan
    data2['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data2['scaled' + var]['inversionForCloudnet'][:] = np.nan
    data3['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data3['scaled' + var]['inversion_Tindex'][:] = np.nan
    data3['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data3['scaled' + var]['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.size(np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))[0][0]
        if np.size(np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))[0][0]
        if np.size(np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data1['scaled' + var]['inversionForCloudnet'][int(data1['scaled' + var]['inversion_Tindex'][i])] = inv1[i]
        if data2['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data2['scaled' + var]['inversionForCloudnet'][int(data2['scaled' + var]['inversion_Tindex'][i])] = inv2[i]
        if data3['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data3['scaled' + var]['inversionForCloudnet'][int(data3['scaled' + var]['inversion_Tindex'][i])] = inv3[i]
    #
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    #
    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(data1['scaledTime'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(data2['scaledTime'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(data3['scaledTime'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(data1['scaledTime'])):        ### all can go in this loop, data1['scaledTime'] == 6-hourly data

        ### main inversion base assignments
        ###         (1) find where UM height array matches the invbase index for index i
        ###         (2) find where UM height array matches the invbase index for index i
        ###         (3) find where IFS height array is less than or equal to the UM-gridded invbase index for index i
        if np.size(np.where(um_data['height'][i,:].data == data1['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaled' + var]['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(misc_data['height'][i,:].data == data2['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(misc_data['height'][i,:].data == data2['scaled' + var]['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(ifs_data['height'][i,:].data <= data3['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            temp = ifs_data['height'][i,:].data <= data3['scaled' + var]['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    fig = plt.figure(figsize=(9,10))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(data1['scaledTime'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1,inv1)
    plt.ylabel('Z [m]')
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(data2['scaledTime'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2, inv2)
    plt.ylim([0,3e3])
    plt.ylabel('Z [m]')
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(data3['scaledTime'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3, inv3)
    plt.ylim([0,3e3])
    plt.ylabel('Z [m]')
    plt.xlabel('DOY')
    plt.savefig('FIGS/' + var + '_model_inversionDetection_timeseries.png')
    plt.show()

    ### set 6 hourly cloudnet Cv arrays as tempvars
    if var == 'Cv':
        ra2m_var = um_data['model_Cv_filtered'][::6,:]
        casim_var = misc_data['model_Cv_filtered'][::6,:]
        ifs_var = ifs_data['model_snow_Cv_filtered'][::6,:]
    elif var == 'lwc':
        ra2m_var = um_data['model_lwc'][::6,:]
        casim_var = misc_data['model_lwc'][::6,:]
        ifs_var = ifs_data['model_lwc'][::6,:]
    elif var == 'iwc':
        ra2m_var = um_data['model_iwc_filtered'][::6,:]
        casim_var = misc_data['model_iwc_filtered'][::6,:]
        ifs_var = ifs_data['model_snow_iwc_filtered'][::6,:]

    ### find all Cv data below identified inversion
    for i in range(0,np.size(data1['scaledTime'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaled' + var]['binned']['t' + str(i)] = {}
        data2['scaled' + var]['binned']['t' + str(i)] = {}
        data3['scaled' + var]['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['bl' + var][i,:int(data1['inversions']['invbase_kIndex'][i]+1)] = ra2m_var[i,:int(data1['inversions']['invbase_kIndex'][i]+1)]
        data2['bl' + var][i,:int(data2['inversions']['invbase_kIndex'][i]+1)] = casim_var[i,:int(data2['inversions']['invbase_kIndex'][i]+1)]
        data3['bl' + var][i,:int(data3['inversions']['invbase_kIndex'][i]+1)] = ifs_var[i,:int(data3['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            ### find mean and stdev of points within given Zpts range
            data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data1['bl' + var][i,tempvar1]
            if np.size(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaled' + var]['mean'][i,k] = np.nanmean(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data1['scaled' + var]['stdev'][i,k] = np.nanstd(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

            data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data2['bl' + var][i,tempvar2]
            if np.size(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaled' + var]['mean'][i,k] = np.nanmean(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data2['scaled' + var]['stdev'][i,k] = np.nanstd(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

            data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data3['bl' + var][i,tempvar3]
            if np.size(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaled' + var]['mean'][i,k] = np.nanmean(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data3['scaled' + var]['stdev'][i,k] = np.nanstd(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

    ### save working data for debug
    # np.save('working_data1', data1)

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    fig = plt.figure(figsize=(8,12))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.xlabel('DOY')
    plt.savefig('FIGS/' + var + '_ALL_scaledZ_timeseries.png')
    plt.show()

    ### obs
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['bl' + var])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),np.squeeze(obs['inversions']['thetaE']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_obs_scaledZ_timeseries.png')
    plt.show()

    ### um_ra2m
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],um_data['height'][0,:].data,np.transpose(data1['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['time'],data1['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_RA2M_scaledZ_timeseries.png')
    plt.show()

    ### um_casim-100
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],misc_data['height'][0,:].data,np.transpose(data2['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['time'],data2['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_CASIM-100_scaledZ_timeseries.png')
    plt.show()

    ### ecmwf_ifs
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],ifs_data['height'][0,:].data,np.transpose(data3['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['time'],data3['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_IFS_scaledZ_timeseries.png')
    plt.show()

    ####================================================================
    ### profiles
    ###         loop through and set all zeros to nans
    obsmean = obs['inversions']['scaled' + var]['mean']
    ra2mmean = data1['scaled' + var]['mean']
    casimmean = data2['scaled' + var]['mean']
    ifsmean = data3['scaled' + var]['mean']
    # for i in range(0, len(data1['scaledTime'])):
    #     obsmean[i,obs['inversions']['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     ra2mmean[i,data1['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     casimmean[i,data2['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     ifsmean[i,data3['scaled' + var]['mean'][i,:] == 0.0] = np.nan

    ##################################################
    ##################################################
    #### create figure and axes instances
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

    ### -------------------------------
    ### Build figure (profiles)
    ### -------------------------------
    fig = plt.figure(figsize=(5,6))
    ax1  = fig.add_axes([0.2,0.1,0.7,0.8])   # left, bottom, width, height
    plt.plot(np.nanmean(obsmean,0),obs['inversions']['scaledZ'], '--', color = 'k', linewidth = 2, label = 'Obs')
    ax1.fill_betweenx(data1['scaledZ'],np.nanmean(obsmean,0) - np.nanstd(obsmean,0),
        np.nanmean(obsmean,0) + np.nanstd(obsmean,0), color = 'lightgrey', alpha = 0.5)
    # ax1.fill_betweenx(data1['scaledZ'],np.nanmean(obsmean,0) - np.nanstd(obs['inversions']['scaled' + var]['stdev'],0),
    #     np.nanmean(obsmean,0) + np.nanstd(obs['inversions']['scaled' + var]['stdev'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(ra2mmean,0),data1['scaledZ'], '^-',
        color = 'steelblue', markeredgecolor = 'midnightblue', linewidth = 2, label = label1)
    ax1.fill_betweenx(data1['scaledZ'],np.nanmean(ra2mmean,0) - np.nanstd(ra2mmean,0),
        np.nanmean(ra2mmean,0) + np.nanstd(ra2mmean,0), color = 'lightblue', alpha = 0.4)
    # ax1.fill_betweenx(data1['scaledZ'],np.nanmean(ra2mmean,0) - np.nanstd(data1['scaled' + var]['stdev'],0),
    #     np.nanmean(ra2mmean,0) + np.nanstd(data1['scaled' + var]['stdev'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifsmean,0),data3['scaledZ'], 'd-',
        color = 'darkorange', markeredgecolor = 'saddlebrown', linewidth = 2, label = label3)
    ax1.fill_betweenx(data3['scaledZ'],np.nanmean(ifsmean,0) - np.nanstd(ifsmean,0),
        np.nanmean(ifsmean,0) + np.nanstd(ifsmean,0), color = 'navajowhite', alpha = 0.35)
    # ax1.fill_betweenx(data3['scaledZ'],np.nanmean(ifsmean,0) - np.nanstd(data3['scaled' + var]['stdev'],0),
    #     np.nanmean(ifsmean,0) + np.nanstd(data3['scaled' + var]['stdev'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(casimmean,0),data2['scaledZ'], 'v-',
        color = 'forestgreen', markeredgecolor = 'darkslategrey', linewidth = 2, label = label2)
    ax1.fill_betweenx(data2['scaledZ'],np.nanmean(casimmean,0) - np.nanstd(casimmean,0),
        np.nanmean(casimmean,0) + np.nanstd(casimmean,0), color = 'mediumaquamarine', alpha = 0.15)
    # ax1.fill_betweenx(data2['scaledZ'],np.nanmean(casimmean,0) - np.nanstd(data2['scaled' + var]['stdev'],0),
    #     np.nanmean(casimmean,0) + np.nanstd(data2['scaled' + var]['stdev'],0), color = 'mediumaquamarine', alpha = 0.15)
    plt.xlim([0,1])
    plt.ylim([0,1])
    if var == 'Cv': plt.xlabel('Cv')
    plt.ylabel('scaled Z \n (0 = lowest level; 1 = inversion base height)')
    plt.legend()
    plt.savefig('FIGS/' + var + '_scaledZ.svg')
    plt.show()

def interpCloudnet(obs_data, month_flag, missing_files, doy):

    from scipy.interpolate import interp1d

    print ('*******')
    print ('Interpolate obs cloudnet field for continuous array:')
    print ('*******')
    print ('')

    varlist = ['Cv', 'lwc', 'iwc']

    for var in varlist:
        ### remove bad and flagged data
        obs_data[var][obs_data[var] < 0.0] = np.nan

        ### save relevant fields as tempvars for ease
        cv = np.copy(obs_data[var].data)
        times = np.copy(obs_data['time'].data)
        height = np.copy(obs_data['height'][0,:])        ### height array constant in time, so just take first column

        ### check if the mean of the column is NaN
        # np.isnan(np.nanmean(cv[0,:]))

        ### need 3 points for interp, so start for loop at i = 2 (remember to finish at i-1!)
        ### check if the column mean == nan but next timestep is non-nan:
        for i in range(2,len(times)-1):
            mn = np.nanmean(cv[i,:])
            # print(str(mn))
            if np.isnan(np.nanmean(cv[i,:])) == True:
                # print ('column i = ' + str(i) + ' needs to be fixed:')
                if np.isnan(np.nanmean(cv[i+1,:])) == False:
                    if np.isnan(np.nanmean(cv[i-1,:])) == False:        ### if the timestep before is non-nan
                        # print (str(i-1) + ' and ' + str(i+1) + ' = yes')
                        cv[i,:] = (cv[i-1,:] + cv[i+1,:]) / 2.0
                        mncv = np.nanmean(cv[i,:])
                        # print ('new mean for i = ' + str(i) + ' is: ' + str(mncv))
                    else:
                        # print ('need to find last non-nan instance...')
                        if np.isnan(np.nanmean(cv[i-2,:])) == False:        ### if the timestep before is non-nan
                            # print (str(i-2) + ' and ' + str(i+1) + ' = yes')
                            cv[i,:] = (cv[i-2,:] + cv[i+1,:]) / 2.0
                            mncv = np.nanmean(cv[i,:])
                            # print ('new mean for i = ' + str(i) + ' is: ' + str(mncv))

        for i in range(2,len(times)-1):
            mn = np.nanmean(cv[i,:])
            # print(str(mn))
            if np.isnan(np.nanmean(cv[i,:])) == True:
                # print ('column i = ' + str(i) + ' needs to be fixed:')
                if np.isnan(np.nanmean(cv[i+1,:])) == False:
                    if np.isnan(np.nanmean(cv[i-1,:])) == False:        ### if the timestep before is non-nan
                        # print (str(i-1) + ' and ' + str(i+1) + ' = yes')
                        cv[i,:] = (cv[i-1,:] + cv[i+1,:]) / 2.0
                        mncv = np.nanmean(cv[i,:])
                        # print ('new mean for i = ' + str(i) + ' is: ' + str(mncv))
                    else:
                        # print ('need to find last non-nan instance...')
                        if np.isnan(np.nanmean(cv[i-2,:])) == False:        ### if the timestep before is non-nan
                            # print (str(i-2) + ' and ' + str(i+1) + ' = yes')
                            cv[i,:] = (cv[i-2,:] + cv[i+1,:]) / 2.0
                            mncv = np.nanmean(cv[i,:])
                            # print ('new mean for i = ' + str(i) + ' is: ' + str(mncv))

        # fig = plt.figure(figsize=(12,6))
        # plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.12,
        #         hspace = 0.3, wspace = 0.3)
        # plt.subplot(211)
        # plt.pcolor(obs_data['time'].data,obs_data['height'][0,:].data,np.transpose(obs_data['Cv'].data), vmin = 0, vmax = 1)
        # plt.xlim([226,258])
        # plt.ylim([0,3000])
        # plt.ylabel('Z [m]')
        # plt.title('original')
        # plt.subplot(212)
        # plt.pcolor(times,height,np.transpose(cv), vmin = 0, vmax = 1)
        # plt.xlim([226,258])
        # plt.ylim([0,3000])
        # plt.xlabel('DOY')
        # plt.ylabel('Z [m]')
        # plt.title('interpolated')
        # plt.savefig('FIGS/Cv_TS_orig_interpd.png')
        # plt.show()
        #
        #
        # fig = plt.figure(figsize=(12,6))
        # plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.12,
        #         hspace = 0.3, wspace = 0.3)
        # plt.subplot(211)
        # plt.pcolor(obs_data['time'][::6].data,obs_data['height'][0,:].data,np.transpose(obs_data['Cv'][::6,:].data), vmin = 0, vmax = 1)
        # plt.xlim([226,258])
        # plt.ylim([0,3000])
        # plt.ylabel('Z [m]')
        # plt.title('original, 6hrly')
        # plt.subplot(212)
        # plt.pcolor(times[::6],height,np.transpose(cv[::6,:]), vmin = 0, vmax = 1)
        # plt.xlim([226,258])
        # plt.ylim([0,3000])
        # plt.xlabel('DOY')
        # plt.ylabel('Z [m]')
        # plt.title('interpolated, 6hrly')
        # plt.savefig('FIGS/Cv_TS_6hrly_orig_interpd.png')
        # plt.show()


        ### save back to dictionary after completion of updates
        obs_data[var] = cv

    return obs_data

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### Choose observations vertical gridding used in Cloudnet processing (UM/IFS/RADAR)
    obs_switch = 'RADAR'

    ### only works on laptop for now

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        um_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        ifs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        obs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/'
        cn_um_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/Cloudnet/UM_RA2M/'
        cn_ifs_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/Cloudnet/ECMWF_IFS/'
        # cn_misc_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/DATA/'; cn_misc_flag = 1              ### FOR NON-CLOUDNET UM DATA
        cn_misc_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/Cloudnet/UM_CASIM-100/'; cn_misc_flag = 0  ### FOR CLOUDNET UM DATA
        cn_obs_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/Cloudnet/Observations/'
    if platform == 'LAPTOP':
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        ifs_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        cn_um_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'
        cn_ifs_dir = '/home/gillian/MOCCHA/Cloudnet/IFS_DATA/'
        # cn_ifs_dir = '/home/gillian/MOCCHA/Cloudnet/OBS_DATA/'
        # cn_misc_dir = '/home/gillian/MOCCHA/UM/DATA/'; cn_misc_flag = 1              ### FOR NON-CLOUDNET UM DATA
        cn_misc_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'; cn_misc_flag = 0  ### FOR CLOUDNET UM DATA
        if obs_switch == 'UM':
            cn_obs_dir = '/home/gillian/MOCCHA/Cloudnet/OBS_DATA/QF30_metum/'
        elif obs_switch == 'IFS':
            cn_obs_dir = '/home/gillian/MOCCHA/Cloudnet/OBS_DATA/QF10_ecmwf/'
        else:
            cn_obs_dir = '/home/gillian/MOCCHA/Cloudnet/OBS_DATA/'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### -----------------------------------------------------------------
    ### CHOSEN RUN - MODEL DATA
    if platform == 'LAPTOP':
        out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/'
        out_dir2 = '5_u-bl661_RA1M_CASIM/OUT_R0/'
        out_dir4 = 'OUT_25H/'
    elif platform == 'JASMIN':
        out_dir1 = 'UM_RA2M/'
        out_dir2 = 'UM_CASIM-100/'
        out_dir4 = 'ECMWF_IFS/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R0/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param
    ### 12_u-br210_RA1M_CASIM/OUT_R0/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507

    ### -----------------------------------------------------------------
    ### CHOSEN RUN - CLOUDNET DATA
    if platform == 'LAPTOP':
        cn_um_out_dir = ['4_u-bg610_RA2M_CON/cloud-fraction-metum-grid/2018/',
                        '4_u-bg610_RA2M_CON/lwc-scaled-metum-grid/2018/',
                        '4_u-bg610_RA2M_CON/iwc-Z-T-metum-grid/2018/']
        cn_ifs_out_dir = ['cloud-fraction-ecmwf-grid/2018/',
                    'lwc-scaled-ecmwf-grid/2018/',
                    'iwc-Z-T-ecmwf-grid/2018/']
        if obs_switch == 'IFS':
            cn_obs_out_dir = cn_ifs_out_dir
        elif obs_switch == 'UM':
            cn_obs_out_dir = ['cloud-fraction-metum-grid/2018/',
                        'lwc-scaled-metum-grid/2018/',
                        'iwc-Z-T-metum-grid/2018/']
        elif obs_switch == 'RADAR':
            cn_obs_out_dir = ['cloud-fraction-ecmwf-grid/2018/',
                        'lwc-scaled-adiabatic/2018/',
                        'iwc-Z-T-method/2018/']
        if cn_misc_flag == 0:       ## flag to compare cloudnet model data
            cn_misc_out_dir = ['5_u-bl661_RA1M_CASIM/cloud-fraction-metum-grid/2018/',
                            '5_u-bl661_RA1M_CASIM/lwc-scaled-metum-grid/2018/',
                            '5_u-bl661_RA1M_CASIM/iwc-Z-T-metum-grid/2018/']
        elif cn_misc_flag == 1:       ## flag to compare non-cloudnet model data
            cn_misc_out_dir = '12_u-br210_RA1M_CASIM/OUT_R0/'
    elif platform == 'JASMIN':
        cn_um_out_dir = 'cloud-fraction-metum-grid/2018/'
        cn_ifs_out_dir = 'cloud-fraction-ecmwf-grid/2018/'
        if obs_switch == 'IFS':
            cn_obs_out_dir = cn_ifs_out_dir
        if cn_misc_flag == 0:       ## flag to compare cloudnet model data
            cn_misc_out_dir = cn_um_out_dir
        elif cn_misc_flag == 1:       ## flag to compare non-cloudnet model data
            cn_misc_out_dir = '12_u-br210_RA1M_CASIM/OUT_R0/'

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

    print ('Misc_flag = ' + str(cn_misc_flag) + '... so third simulation for Cloudnet comparison is:')
    if cn_misc_flag == 0: print ('Cloudnet-ed data!')
    if cn_misc_flag == 1: print ('standard model output!')

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

###################################################################################################################
###################################################################################################################

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    print ('Loading observations:')
            # -------------------------------------------------------------
            # Which file does what?
            # -------------------------------------------------------------
            #### ice station: net LW / net SW
                    #### obs['ice_station_fluxes']/mast_radiation_30min_v2.3.mat
            #### obs['foremast']:
                    #### obs['foremast']/ACAS_AO2018_obs['foremast']_30min_v2_0.nc
            #### 7th deck: temperature, surface temperature, RH, downwelling SW, downwelling LW
                    #### 7thDeck/ACAS_AO2018_WX_30min_v2_0.nc

    obs = {}

    if platform == 'LAPTOP':
        print ('Load temporary ice station data from Jutta...')
        obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_wTemp1p5m.nc','r')

        print ('Load ice station flux data from Jutta...')
        obs['ice_station_fluxes'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')

        print ('Load cloud fraction data from Jutta/Peggy...')
        obs['cloudfractions'] = readMatlabStruct(cn_obs_dir + '../Gillian_cloudfraction.mat')

    # print ('Load ice station radiation data from Jutta...')
    # obs['ice_station_radiation'] = readMatlabStruct(obs_root_dir + 'ice_station/mast_radiation_30min_v2.3.mat')

    print ('Load radiosonde data from Jutta...')
    obs['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print ('Load observations inversion height data from Jutta...')
    obs['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/InversionHeights_RSh05int_final_V03.mat')

    print ('Load foremast data from John...')
    obs['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print ('Load 7th deck weather station data from John...')
    obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print ('Load weather sensor data from John...')
    obs['pws'] = readMatlabStruct(obs_root_dir + '7thDeck/ACAS_AO2018_PWD_30min_v1_0.mat')

    print ('...')

###################################################################################################################
###################################################################################################################

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin nc read in at ' + time.strftime("%c"))
    print (' ')

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

    moccha_names = ['20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_',
            '20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180813_oden_','20180818_oden_','20180910_oden_','20180914_oden_']   ### cloud radar not working
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)        ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)
    # doy = np.arange(244,256)        ## set DOY for CASIM-AeroProf (1st Sep to 11th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Choose month:
    names = moccha_names
    # missing_files = moccha_missing_files
    month_flag = -1


###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

    for i in range(0,len(names)):

        ### --------------------------------------------------------------------
        #### DEFINE FILENAMES TO BE READ IN
        ### --------------------------------------------------------------------
        print ('Load raw model data first: ')
        filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
        filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
        if np.logical_or(out_dir4 == 'OUT_25H/',out_dir4 == 'ECMWF_IFS/'):
            print( '***IFS being compared***')
            ifs_flag = True
            filename_um3 = ifs_root_dir + out_dir4 + names[i] + 'ecmwf.nc'
        else:
            print ('***IFS NOT being compared***')
            print (filename_um1)
            filename_um3 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            ifs_flag = False
        print (filename_um2)
        print (filename_um3)
        print ('')

        print ('Now load cloudnet data:')
        # if cn_um_out_dir[-31:-6] == 'cloud-fraction-metum-grid':
        #     cn_out_dir = 'cloud-fraction-metum-grid'
        # elif cn_um_out_dir[-27:-6] == 'lwc-scaled-metum-grid':
        #     cn_out_dir = 'lwc-scaled-metum-grid'
        # elif cn_um_out_dir[-24:-6] == 'iwc-Z-T-metum-grid':
        #     cn_out_dir = 'iwc-Z-T-metum-grid'
        cn_filename_um = [cn_um_dir + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc']
        cn_filename_ifs = [cn_ifs_dir + cn_ifs_out_dir[0] + names[i] + cn_ifs_out_dir[0][:-6] + '.nc',
                        cn_ifs_dir + cn_ifs_out_dir[1] + names[i] + cn_ifs_out_dir[1][:-6] + '.nc',
                        cn_ifs_dir + cn_ifs_out_dir[2] + names[i] + cn_ifs_out_dir[2][:-6] + '.nc']
        cn_filename_obs = [cn_obs_dir + cn_obs_out_dir[0] + names[i] + cn_obs_out_dir[0][:-6] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[1] + names[i] + cn_obs_out_dir[1][:-6] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[2] + names[i] + cn_obs_out_dir[2][:-6] + '.nc']
        if cn_misc_flag == 1: cn_filename_misc = cn_misc_dir + cn_misc_out_dir[0] + names[i] + 'metum.nc'
        if cn_misc_flag == 0:
            cn_filename_misc = [cn_misc_dir + cn_misc_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                            cn_misc_dir + cn_misc_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                            cn_misc_dir + cn_misc_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc']
        print (cn_filename_um)
        print (cn_filename_ifs)
        if cn_misc_flag != 1: print (cn_filename_misc)
        print ('')

        ### --------------------------------------------------------------------
        #### CHOOSE MODEL DIAGNOSTICS FIRST
        ### --------------------------------------------------------------------
        var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux','latent_heat_flux',
            'temp_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice']
        var_list2 = var_list1
        if ifs_flag: var_list3 = ['height', 'flx_height', 'temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
            'sfc_temp_2m','flx_ls_rain','flx_conv_rain','flx_ls_snow','q','pressure','sfc_bl_height','ql','qi']
        if not ifs_flag: var_list3 = var_list1

        if names[i] in moccha_missing_files:        ### NOTE THIS WON'T WORK IF IT'S THE FIRST FILE THAT'S MISSING!!
            print ('File not available...')
            print ('***Filling arrays with nans***')
            print (str(doy[i]))
            nanarray = np.zeros(24)
            nanarray[:] = np.nan
            timarray = np.arange(0.041666666666666664,1.01,0.041666666666666664)
            time_um1 = np.append(time_um1, doy[i] + timarray)
            time_um2 = np.append(time_um2, doy[i] + timarray)
            if ifs_flag:
                time_um3 = np.append(time_um3, doy[i] + timarray)

            for j in range(0,len(var_list1)):
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                    nanarray = np.zeros(24)
                    nanarray[:] = np.nan
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nanarray)
                    data2[var_list2[j]] = np.append(data2[var_list2[j]],nanarray)
                elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                    nanarray = np.zeros([24,71])
                    nanarray[:] = np.nan
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nanarray,0)
                    data2[var_list2[j]] = np.append(data2[var_list2[j]],nanarray,0)
            for j in range(0,len(var_list3)):
                print (j)
                print (var_list3[j])
                # np.save('testing', data3)
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                    nanarray = np.zeros(24)
                    nanarray[:] = np.nan
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nanarray)
                elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                    if var_list3[j][:3] == 'flx':
                        nanarray = np.zeros([24,138])
                    else:
                        nanarray = np.zeros([24,137])
                    nanarray[:] = np.nan
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nanarray,0)

        else:
            ### --------------------------------------------------------------------
            ###     READ IN ALL MODEL FILES
            ### --------------------------------------------------------------------
            print( 'Loading first run diagnostics:')
            nc1 = Dataset(filename_um1,'r')
            print ('...')
            print( 'Loading second run diagnostics:')
            nc2 = Dataset(filename_um2,'r')
            print ('...')
            print ('Loading third run diagnostics:')
            nc3 = Dataset(filename_um3,'r')
            print ('...')

            ### --------------------------------------------------------------------
            ###     READ IN ALL CLOUDNET FILES
            ### --------------------------------------------------------------------
            print ('Loading multiple diagnostics:')
            cn_nc1 = {}
            cn_nc2 = {}
            cn_nc3 = {}
            cn_nc0 = {}
            for c in range(0,3):
                cn_nc1[c] = Dataset(cn_filename_um[c],'r')
                cn_nc2[c] = Dataset(cn_filename_ifs[c],'r')
                if cn_misc_flag != -1: cn_nc3[c] = Dataset(cn_filename_misc[c],'r')
                cn_nc0[c] = Dataset(cn_filename_obs[c],'r')

            # -------------------------------------------------------------
            print ('')

    ###################################################################################################################
    ### ---------------------------------------------------------------------------------------------------------------

            if i == 0:      ### first file
                ## ------------------
                #### UM
                ## ------------------
                data1 = {}
                data2 = {}
                data3 = {}
                if month_flag == -1:
                    time_um1 = doy[i] + (nc1.variables['forecast_time'][:]/24.0)
                    time_um2 = doy[i] + (nc2.variables['forecast_time'][:]/24.0)
                    if ifs_flag: time_um3 = doy[i] + (nc3.variables['time'][:]/24.0)
                    if not ifs_flag: time_um3 = doy[i] + (nc3.variables['forecast_time'][:]/24.0)
                else:
                    time_um1 = float(filename_um1[-16:-14]) + (nc1.variables['forecast_time'][:]/24.0)
                    time_um2 = float(filename_um2[-16:-14]) + (nc2.variables['forecast_time'][:]/24.0)
                    if ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['time'][:]/24.0)
                    if not ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['forecast_time'][:]/24.0)

                ### define height arrays explicitly
                data1['height'] = nc1.variables['height'][:]
                data2['height'] = nc2.variables['height'][:]
                if not ifs_flag: data3['height'] = nc3.variables['height'][:]
                # data3['height'] = nc3.variables['height'][:]
                # data3['flx_height'] = nc3.variables['flx_height'][:]

                for j in range(0,len(var_list1)):
                    if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc1.variables[var_list1[j]]) >= 1:
                        data1[var_list1[j]] = nc1.variables[var_list1[j]][:]
                nc1.close()
                ## ------------------
                #### um2
                ## ------------------
                for j in range(0,len(var_list2)):
                    if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
                        data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
                nc2.close()
                ## ------------------
                #### um3
                ## ------------------
                for j in range(0,len(var_list3)):
                    if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc3.variables[var_list3[j]]) >= 1:
                        data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
                nc3.close()
            else:
                # if doy[i] in missing_files:
                #     print (str(doy[i]))
                #     nanarray = np.zeros(24)
                #     nanarray[:] = np.nan
                #     time_um1 = np.append(time_um1, nanarray)
                #     time_um2 = np.append(time_um2, nanarray)
                #     if ifs_flag: time_um3 = np.append(time_um3, nanarray)
                #     if not ifs_flag: time_um3 = np.append(time_um3, nanarray)
                #     for j in range(0,len(var_list1)):
                #         if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                #             continue
                #         elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                #             nanarray = np.zeros(24)
                #             nanarray[:] = np.nan
                #             data1[var_list1[j]] = np.append(data1[var_list1[j]],nanarray)
                #             data2[var_list2[j]] = np.append(data2[var_list2[j]],nanarray)
                #         elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                #             nanarray = np.zeros([24,71])
                #             nanarray[:] = np.nan
                #             data1[var_list1[j]] = np.append(data1[var_list1[j]],nanarray,0)
                #             data2[var_list2[j]] = np.append(data2[var_list2[j]],nanarray,0)
                #     for j in range(0,len(var_list3)):
                #         print (j)
                #         print (var_list3[j])
                #         np.save('testing', data3)
                #         if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                #             continue
                #         elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                #             nanarray = np.zeros(24)
                #             nanarray[:] = np.nan
                #             data3[var_list3[j]] = np.append(data3[var_list3[j]],nanarray)
                #         elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                #             if var_list3[j][:3] == 'flx':
                #                 nanarray = np.zeros([24,138])
                #             else:
                #                 nanarray = np.zeros([24,137])
                #             nanarray[:] = np.nan
                #             data3[var_list3[j]] = np.append(data3[var_list3[j]],nanarray,0)


                time_um1 = np.append(time_um1, doy[i] + (nc1.variables['forecast_time'][:]/24.0))
                time_um2 = np.append(time_um2, doy[i] + (nc2.variables['forecast_time'][:]/24.0))
                if ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['time'][:]/24.0))
                if not ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['forecast_time'][:]/24.0))

                ## ------------------
                #### UM
                ## ------------------
                for j in range(0,len(var_list1)):
                    if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                        # data1[cube_um1[j].var_name] = cube_um1[j].data
                        data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:])
                    elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                        data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:],0)
                nc1.close()
                ## ------------------
                #### um2
                ## ------------------
                for j in range(0,len(var_list2)):
                    if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc2.variables[var_list2[j]]) == 1:
                        data2[var_list2[j]] = np.append(data2[var_list2[j]],nc2.variables[var_list2[j]][:])
                    elif np.ndim(nc2.variables[var_list2[j]]) == 2:
                        data2[var_list2[j]] = np.append(data2[var_list2[j]],nc2.variables[var_list2[j]][:],0)
                nc2.close()
                ## ------------------
                #### um3 / ifs
                ## ------------------
                for j in range(0,len(var_list3)):
                    if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                        data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:])
                    elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                        data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:],0)
                nc3.close()

###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

            ### --------------------------------------------------------------------
            ###     LOAD UM CLOUDNET DIAGS INTO DICTIONARY
            ### --------------------------------------------------------------------
            #### LOAD IN SPECIFIC DIAGNOSTICS
            um_var_list = [['Cv','model_Cv_filtered','model_temperature'],
                    ['lwc','lwp','model_lwc','model_lwp'],
                    ['height','iwc','model_iwc','model_iwc_filtered']]   ### time always read in separately

            ###     LOOP OVER TIME DUMP
            if i == 0:
                um_data = {}
                # um_data1d = {}
                if month_flag == -1:
                    time_um = doy[i] + ((cn_nc1[0].variables['time'][:])/24.0)
                else:
                    time_um = float(names[i][6:8]) + ((cn_nc1[0].variables['time'][:])/24.0)
                for c in range(0,3):
                    for j in range(0,len(um_var_list[c])):
                        if np.ndim(cn_nc1[c].variables[um_var_list[c][j]]) == 1:  # 1d timeseries only
                            um_data[um_var_list[c][j]] = cn_nc1[c].variables[um_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            um_data[um_var_list[c][j]] = cn_nc1[c].variables[um_var_list[c][j]][:]
            else:
                if month_flag == -1:
                    time_um = np.append(time_um, doy[i] + ((cn_nc1[0].variables['time'][:])/24.0))
                else:
                    time_um = np.append(time_um,float(cn_filename_um[-16:-14]) + ((cn_nc1[0].variables['time'][:])/24.0))
                print (um_data)
                for c in range(0,3):
                    for j in range(0,len(um_var_list[c])):
                        # print 'j = ' + str(j)
                        if np.ndim(cn_nc1[c].variables[um_var_list[c][j]]) == 1:
                            um_data[um_var_list[c][j]] = np.append(um_data[um_var_list[c][j]],cn_nc1[c].variables[um_var_list[c][j]][:])
                        else:
                            um_data[um_var_list[c][j]] = np.append(um_data[um_var_list[c][j]],cn_nc1[c].variables[um_var_list[c][j]][:],0)
            for c in range(0,3): cn_nc1[c].close()

            ### --------------------------------------------------------------------
            ### LOAD IN IFS DATA INTO DICTIONARY
            ### --------------------------------------------------------------------
            ifs_var_list = [['Cv','model_snow_Cv_filtered','model_temperature'],
                    ['lwc','lwp','model_lwc','model_lwp'],
                    ['height','iwc','model_iwc','model_snow_iwc_filtered','model_iwc_filtered']]   ### time always read in separately

            ###     LOOP OVER TIME DUMP
            if i == 0:
                ifs_data = {}
                # ifs_data1d = {}
                if month_flag == -1:
                    time_ifs = doy[i] + ((cn_nc2[0].variables['time'][:])/24.0)
                else:
                    time_ifs = float(names[i][6:8]) + ((cn_nc2[0].variables['time'][:])/24.0)
                for c in range(0,3):
                    for j in range(0,len(ifs_var_list[c])):
                        if np.ndim(cn_nc2[c].variables[ifs_var_list[c][j]]) == 1:  # 1d timeseries only
                            ifs_data[ifs_var_list[c][j]] = cn_nc2[c].variables[ifs_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            ifs_data[ifs_var_list[c][j]] = cn_nc2[c].variables[ifs_var_list[c][j]][:]
            else:
                if month_flag == -1:
                    time_ifs = np.append(time_ifs, doy[i] + ((cn_nc2[0].variables['time'][:])/24.0))
                else:
                    time_ifs = np.append(time_ifs,float(cn_filename_ifs[-16:-14]) + ((cn_nc2[0].variables['time'][:])/24.0))
                print (ifs_data)
                for c in range(0,3):
                    for j in range(0,len(ifs_var_list[c])):
                        ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                        # print 'j = ' + str(j)
                        if np.ndim(cn_nc2[c].variables[ifs_var_list[c][j]]) == 1:
                            ifs_data[ifs_var_list[c][j]] = np.append(ifs_data[ifs_var_list[c][j]],cn_nc2[c].variables[ifs_var_list[c][j]][:])
                        else:
                            ifs_data[ifs_var_list[c][j]] = np.append(ifs_data[ifs_var_list[c][j]],cn_nc2[c].variables[ifs_var_list[c][j]][:],0)
            for c in range(0,3): cn_nc2[c].close()

            ### -------------------------------------------------------------------------
            ###     LOAD IN MISC DATA INTO DICTIONARY IF COMPARING
            ###             Only load in what variables are needed based on IFS file chosen
            ### -------------------------------------------------------------------------
            if cn_misc_flag == -1:
                continue
            elif cn_misc_flag == 1:
                misc_var_list = [['cloud_fraction','temperature'],
                        ['qliq'],
                        ['qice']]   ### time always read in separately
            elif cn_misc_flag == 0:
                misc_var_list = [['Cv','model_Cv_filtered','model_temperature'],
                        ['lwc','lwp','model_lwc','model_lwp'],
                        ['height','iwc','model_iwc','model_iwc_filtered']]   ### time always read in separately

            print ('')
            print ('misc file variable list is:')
            print (misc_var_list)
            print ('')

            if i == 0:
                misc_data = {}
                # misc_data1d = {}
                if month_flag == -1:
                    if cn_misc_flag == 1:
                        time_misc = doy[i] + ((cn_nc3[0].variables['forecast_time'][:])/24.0)
                        misc_data['height'] = cn_nc3[0].variables['height'][:]
                    if cn_misc_flag == 0: time_misc = doy[i] + ((cn_nc3[0].variables['time'][:])/24.0)
                else:
                    if cn_misc_flag == 1: time_misc = float(names[i][6:8]) + ((cn_nc3[0].variables['forecast_time'][:])/24.0)
                    if cn_misc_flag == 0: time_misc = float(names[i][6:8]) + ((cn_nc3[0].variables['time'][:])/24.0)
                for c in range(0,3):
                    for j in range(0,len(misc_var_list[c])):
                        if np.ndim(cn_nc3[c].variables[misc_var_list[c][j]]) == 1:  # 1d timeseries only
                            misc_data[misc_var_list[c][j]] = cn_nc3[c].variables[misc_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            misc_data[misc_var_list[c][j]] = cn_nc3[c].variables[misc_var_list[c][j]][:]
            else:
                if month_flag == -1:
                    if cn_misc_flag == 1: time_misc = np.append(time_misc, doy[i] + ((cn_nc3[0].variables['forecast_time'][:])/24.0))
                    if cn_misc_flag == 0: time_misc = np.append(time_misc, doy[i] + ((cn_nc3[0].variables['time'][:])/24.0))
                else:
                    if cn_misc_flag == 1: time_misc = np.append(time_misc,float(cn_filename_misc[-16:-14]) + ((cn_nc3[0].variables['forecast_time'][:])/24.0))
                    if cn_misc_flag == 0: time_misc = np.append(time_misc,float(cn_filename_misc[-16:-14]) + ((cn_nc3[0].variables['time'][:])/24.0))
                print (misc_data)
                for c in range(0,3):
                    for j in range(0,len(misc_var_list[c])):
                        # print 'j = ' + str(j)
                        if np.ndim(cn_nc3[c].variables[misc_var_list[c][j]]) == 1:
                            misc_data[misc_var_list[c][j]] = np.append(misc_data[misc_var_list[c][j]],cn_nc3[c].variables[misc_var_list[c][j]][:])
                        # elif var_list[j] == 'height':#np.sum(nc3.variables[var_list[j]].shape) == 71:
                        #     continue
                        else:
                            misc_data[misc_var_list[c][j]] = np.append(misc_data[misc_var_list[c][j]],cn_nc3[c].variables[misc_var_list[c][j]][:],0)
            for c in range(0,3): cn_nc3[c].close()

            ### -------------------------------------------------------------------------
            ###     LOAD IN OBS DATA
            ###             Only load in what variables are needed based on IFS file chosen
            ### -------------------------------------------------------------------------
            obs_var_list = [['Cv', 'Cv_adv'],
                        ['lwc','lwp'],
                        ['height','iwc']]

            if obs_switch == 'RADAR':
                if i == 0:
                    obs_data = {}
                    # misc_data1d = {}
                    if month_flag == -1:
                        time_obs = doy[i] + ((cn_nc0[1].variables['time'][:])/24.0)
                    else:
                        time_obs = float(names[i][6:8]) + ((cn_nc0[1].variables['time'][:])/24.0)
                    for c in range(0,3):
                        for j in range(0,len(obs_var_list[c])):
                            if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:  # 1d timeseries only
                                obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                            else:                                   # 2d column um_data
                                obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                else:
                    if month_flag == -1:
                        time_obs = np.append(time_obs, doy[i] + ((cn_nc0[1].variables['time'][:])/24.0))
                    else:
                        time_obs = np.append(time_obs,float(cn_filename_obs[-16:-14]) + ((cn_nc0[1].variables['time'][:])/24.0))
                    print (obs_data)
                    for c in range(0,3):
                        for j in range(0,len(obs_var_list[c])):
                            # print 'j = ' + str(j)
                            if obs_var_list[c][j] == 'height':
                                continue
                            elif np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:
                                obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:])
                            elif np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 2:
                                obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:],0)
                for c in range(0,3): cn_nc0[c].close()

            else:
                if i == 0:
                    obs_data = {}
                    # misc_data1d = {}
                    if month_flag == -1:
                        time_obs = doy[i] + ((cn_nc0[0].variables['time'][:])/24.0)
                    else:
                        time_obs = float(names[i][6:8]) + ((cn_nc0[0].variables['time'][:])/24.0)
                    for c in range(0,3):
                        for j in range(0,len(obs_var_list[c])):
                            if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:  # 1d timeseries only
                                obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                            else:                                   # 2d column um_data
                                obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                else:
                    if month_flag == -1:
                        time_obs = np.append(time_obs, doy[i] + ((cn_nc0[0].variables['time'][:])/24.0))
                    else:
                        time_obs = np.append(time_obs,float(cn_filename_obs[-16:-14]) + ((cn_nc0[0].variables['time'][:])/24.0))
                    print (obs_data)
                    for c in range(0,3):
                        for j in range(0,len(obs_var_list[c])):
                            # print 'j = ' + str(j)
                            if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:
                                obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:])
                            elif np.sum(cn_nc0[c].variables[obs_var_list[c][j]].shape) == 71:
                                continue
                            else:
                                obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:],0)
                for c in range(0,3): cn_nc0[c].close()


    #################################################################
    ## save time to dictionaries now we're not looping over all diags anymore
    #################################################################
    ### models
    data1['time'] = time_um1
    data2['time'] = time_um2
    data3['time'] = time_um3

    ### cloudnet
    ifs_data['time'] = time_ifs
    um_data['time'] = time_um
    if cn_misc_flag != -1: misc_data['time'] = time_misc
    obs_data['time'] = time_obs

    #################################################################
    ### stop double counting of 0000 and 2400 from model data
    #################################################################
    temp = np.zeros([len(data1['time'])])
    for i in range(0, len(temp)-1):
        if data1['time'][i] == data1['time'][i+1]:
            continue
        else:
            temp[i] = data1['time'][i]
    ii = np.where(temp != 0.0)      ### picks out where data are non-zero

    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_hrly'] = temp[ii]
    data2['time_hrly'] = temp[ii]
    data3['time_hrly'] = temp[ii]
    data1['hrly_flag'] = ii
    data2['hrly_flag'] = ii
    data3['hrly_flag'] = ii

    #################################################################
    ### if RADAR flag used, force 2D height array for function compatibility
    #################################################################
    # if obs_switch == 'RADAR':
    #     tmp_height = np.zeros([np.size(obs_data['time']), np.size(obs_data['height'])])
    #     for t in range(0,len(obs_data['time'])): tmp_height[t,:] = obs_data['height'][:]
    #     obs_data['height'] = tmp_height

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '13_u-br409': label1 = 'UM_CASIM-100_AP'
    if out_dir1[:10] == '12_u-br210': label1 = 'UM_CASIM-AeroProf'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T'
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'

    label2 = 'undefined_label'
    if out_dir2[:10] == '13_u-br409': label2 = 'UM_CASIM-100_AP'
    if out_dir2[:10] == '12_u-br210': label2 = 'UM_CASIM-AeroProf'
    if out_dir2[:10] == '11_u-bq798': label2 = 'UM_CASIM-100_Meyers'
    if out_dir2[:10] == '10_u-bq791': label2 = 'UM_CASIM-100_Fletcher'
    if out_dir2[:9] == '8_u-bp738': label2 = 'UM_ERAI-GLM'
    if out_dir2[:9] == '7_u-bn068': label2 = 'UM_RA2T'
    if out_dir2[:9] == '6_u-bm410': label2 = 'UM_CASIM-200'
    if out_dir2[:9] == '5_u-bl661': label2 = 'UM_CASIM-100'
    if out_dir2[:9] == '4_u-bg610': label2 = 'UM_RA2M'

    label3 = 'undefined_label'
    if np.logical_or(out_dir4 == 'OUT_25H/',out_dir4 == 'ECMWF_IFS/'): label3 = 'ECMWF_IFS'
    if out_dir4[:10] == '13_u-br409': label3 = 'UM_CASIM-100_AP'
    if out_dir4[:10] == '12_u-br210': label3 = 'UM_CASIM-AeroProf'
    if out_dir4[:10] == '11_u-bq798': label3 = 'UM_CASIM-100_Meyers'
    if out_dir4[:10] == '10_u-bq791': label3 = 'UM_CASIM-100_Fletcher'
    if out_dir4[:9] == '8_u-bp738': label3 = 'UM_ERAI-GLM'
    if out_dir4[:9] == '7_u-bn068': label3 = 'UM_RA2T'
    if out_dir4[:9] == '6_u-bm410': label3 = 'UM_CASIM-200'
    if out_dir4[:9] == '5_u-bl661': label3 = 'UM_CASIM-100'
    if out_dir4[:9] == '4_u-bg610': label3 = 'UM_RA2M'

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    ### model/measurement data
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    # np.save('working_dataObs', obs['sondes'])
    #
    # ### cloudnet
    # np.save('working_um_data', um_data)
    # np.save('working_ifs_data', ifs_data)
    # if cn_misc_flag != -1: np.save('working_misc_data', misc_data)

###################################################################################################################
###################################################################################################################
################################################ FIGURES ##########################################################
###################################################################################################################
###################################################################################################################

    # -------------------------------------------------------------
    # Cloudnet plot: Plot Cv statistics from drift period
    # -------------------------------------------------------------
    # figure = plot_CvProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs, obs_switch)
    # figure = plot_lwcProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)

    # -------------------------------------------------------------
    # Cloudnet plot: Plot contour timeseries
    # -------------------------------------------------------------

    # obs_data = interpCloudnet(obs_data, month_flag, missing_files, doy)
    # figure = plot_CvTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_LWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_IWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_TWCTimeseries(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_TWCTesting(um_data, ifs_data, misc_data, obs_data, data1, data2, data3, obs, month_flag, missing_files, doy)

    # -------------------------------------------------------------
    # Model plots
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # plot LWP timeseries with missing files accounted for
    # -------------------------------------------------------------
    figure = plot_LWP(um_data, ifs_data, misc_data, obs_data, obs, month_flag, missing_files, cn_um_out_dir, doy, obs_switch) #, lon, lat):

    # -------------------------------------------------------------
    # make obs comparison fig between um and ifs grids
    # -------------------------------------------------------------
    # figure = plot_ObsGridComparison(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)

    # -------------------------------------------------------------
    # plot cloudnet split season figures with missing files accounted for
    # -------------------------------------------------------------
    # figure = plot_CvProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_lwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)
    # figure = plot_iwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)

    # -------------------------------------------------------------
    # cloud properties scaled by BL depth
    # -------------------------------------------------------------
    # -------------------------------------------------------------
    # Identify inversions in model / radiosonde data
    # -------------------------------------------------------------
    #################################################################
    ## load calculated model inversion heights
    #################################################################
    # print ('Load calculated model inversion heights (JV algorithm)...')
    # data1['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_RA2M_inversion_results.mat')
    # data2['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_CASIM-100_inversion_results.mat')
    # data3['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/ECMWF_IFS_inversion_results.mat')

    # print ('Load calculated model inversion heights (GY algorithm)...')
    # obs['inversions']['thetaE'] = np.load(um_root_dir[:-5] + 'obs_inversions_v2.npy').item()
    # data1['inversions'] = np.load(um_root_dir[:-5] + 'um_ra2m_inversions_v2.npy').item()
    # data2['inversions'] = np.load(um_root_dir[:-5] + 'um_casim-100_inversions_v2.npy').item()
    # data3['inversions'] = np.load(um_root_dir[:-5] + 'ecmwf_ifs_inversions_v2.npy').item()

    # -------------------------------------------------------------
    ### use IFS named directory to allocate variable to plot
    # -------------------------------------------------------------
    # if cn_ifs_out_dir[0] == 'cloud-fraction-ecmwf-grid/2018/': var = 'Cv'
    # if cn_ifs_out_dir[0] == 'lwc-scaled-ecmwf-grid/2018/': var = 'lwc'
    # if cn_ifs_out_dir[0] == 'iwc-Z-T-ecmwf-grid/2018/': var = 'iwc'

    # obs_data = interpCloudnet(obs_data, month_flag, missing_files, doy)
    # figure = plot_scaledBL_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3, var)
    # figure = plot_scaledBLCv_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_scaledBLCv_JVInv(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_scaledBLlwc(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)


    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    ### model/measurement data
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['inversions'])

    ### cloudnet
    np.save('working_um_data', um_data)
    np.save('working_ifs_data', ifs_data)
    if cn_misc_flag != -1: np.save('working_misc_data', misc_data)
    np.save('working_obs_data', obs_data)

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')

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
