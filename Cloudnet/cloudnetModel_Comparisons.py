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
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
        np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
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
    # plt.savefig(fileout)
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

def plot_scaledBL(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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

    #### split into melt and freeze in scatter plots:
    ####        all model data share a timestamp
    melt = np.where(np.logical_and(data1['time_hrly'] >= obs['inversions']['doy_drift'][0], data1['time_hrly'] < 240.0))
    freeze = np.where(data1['time_hrly'] >= 240.0)
    #### allow some leeway in radiosonde timesteps
    obsmelt = np.where(np.logical_and(obs['inversions']['doy'] >= data1['time_hrly'][0]-0.2, obs['inversions']['doy'] < 240.0))
    obsfreeze = np.where(np.logical_and(obs['inversions']['doy'] >= 240.0, obs['inversions']['doy'] <= data1['time_hrly'][-1]))

    # print ('obsmelt.shape = ', obsmelt.shape)

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    #### ---------------------------------------------------------------
    #### make inversion tempvars to allow for easy subsampling
    inv1 = np.squeeze(data1['inversions']['invbase'][data1['hrly_flag'],0])
    inv2 = np.squeeze(data2['inversions']['invbase'][data2['hrly_flag'],0])
    inv3 = np.squeeze(data3['inversions']['invbase'][data3['hrly_flag'],0])

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
    #### Look at data below main inversion base only - model data first
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(inv1))
    zind2 = np.zeros(np.size(inv2))
    zind3 = np.zeros(np.size(inv3))#,np.size(data3['height'],1))

    #### fill arrays with height index of main inversion base
    for i in range(0, np.size(inv1)):        ### all can go in this loop, inv1 == hourly data
        # print (i)
        if np.size(np.where(data1['height'][1:].data == inv1[i])) > 0.0:
            zind1[i] = np.where(data1['height'][1:].data == inv1[i])[0][0]
        else:
            zind1[i] = np.nan
        if np.size(np.where(data2['height'][1:].data == inv2[i])) > 0.0:
            zind2[i] = np.where(data2['height'][1:].data == inv2[i])[0][0]
        else:
            zind2[i] = np.nan
        if np.size(np.where(data3['height_hrly'][i].data <= inv3[i])) > 0.0:
            temp = data3['height_hrly'][i,:].data <= inv3[i]
            zind3[i] = np.where(temp == True)[0][-1]
            # zind3[i] = np.where(data3['height_hrly'][i,:].data == inv3[i])[0][0]
        else:
            zind3[i] = np.nan

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    # print (zind3)
    #### re-check inversion algorithm success rate to make sure no !0 values got dropped
    zzind1 = np.where(data1['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind1rate = np.size(zzind1) / np.float(np.size(inv1)) * 100.0
    zzind2 = np.where(data2['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind2rate = np.size(zzind2) / np.float(np.size(inv2)) * 100.0
    zzind3 = np.where(data3['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind3rate = np.size(zzind3) / np.float(np.size(inv3)) * 100.0

    # print (label1 + ' zind1 success rate = ' + str(zind1rate))
    # print (label2 + ' zind2 success rate = ' + str(zind2rate))
    # print (label3 + ' zind3 success rate = ' + str(zind3rate))
    # print ('****')

    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged um_data to nans
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------
    # um['scaled_height'] =
    ### try i = 0 first to see if it works
    i = 0
    hgts1 = data1['height'][0:int(data1['inversions']['invbase_kIndex'][i]+1)]
    scaled_hgts1 = hgts1 / data1['height'][int(data1['inversions']['invbase_kIndex'][i])]
    cv1 = um['Cv'][i,0:int(data1['inversions']['invbase_kIndex'][i]+1)]

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
    p3 = np.where(um['time'] < 230.0)
    p4 = np.where(np.logical_and(um['time'] >= 230.0, um['time'] < 240.0))
    p5 = np.where(np.logical_and(um['time'] >= 240.0, um['time'] < 247.0))
    p6 = np.where(np.logical_and(um['time'] >= 247.0, um['time'] < 251.0))
    p7 = np.where(np.logical_and(um['time'] >= 251.0, um['time'] < 255.5))
    p8 = np.where(um['time'] >= 255.5)

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
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(5,7))
    ax = plt.gca()

    # plt.plot(np.nanmean(um_data['Cv'],0),np.nanmean(um_data['height'],0), 'k--', linewidth = 3, label = 'Obs')
    # ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['Cv'],0) - np.nanstd(um_data['Cv'],0),
    #     np.nanmean(um_data['Cv'],0) + np.nanstd(um_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2M')
    # ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
    #     np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.4)
    # plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'darkorange', linewidth = 3, label = 'ECMWF_IFS')
    # ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
    #     np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0),np.nanmean(misc_data['height'],0), color = 'forestgreen', linewidth = 3, label = 'UM_CASIM-100')
    # ax.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_Cv_filtered'],0) - np.nanstd(misc_data['model_Cv_filtered'],0),
    #     np.nanmean(misc_data['model_Cv_filtered'],0) + np.nanstd(misc_data['model_Cv_filtered'],0), color = 'mediumaquamarine', alpha = 0.15)
    # plt.xlabel('Cloud Fraction')
    # plt.ylabel('Height [m]')
    # plt.ylim([0,10000])
    # plt.xlim([0,1])
    # plt.legend()

    plt.plot(data1['inversions']['invbase_kIndex'])
    plt.plot(data2['inversions']['invbase_kIndex'])
    plt.plot(data3['inversions']['invbase_kIndex'])




    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_oden_metum_ifs_casim-100.svg'
    # fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_metum.svg'
    # plt.savefig(fileout)
    plt.show()

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

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
        # cn_misc_dir = '/home/gillian/MOCCHA/UM/DATA/'; cn_misc_flag = 1              ### FOR NON-CLOUDNET UM DATA
        cn_misc_dir = '/home/gillian/MOCCHA/Cloudnet/UM_DATA/'; cn_misc_flag = 0  ### FOR CLOUDNET UM DATA
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
        cn_um_out_dir = '4_u-bg610_RA2M_CON/cloud-fraction-metum-grid/2018/'
        cn_ifs_out_dir = 'cloud-fraction-ecmwf-grid/2018/'
        cn_obs_out_dir = cn_ifs_out_dir
        if cn_misc_flag == 0:       ## flag to compare cloudnet model data
            cn_misc_out_dir = '5_u-bl661_RA1M_CASIM/cloud-fraction-metum-grid/2018/'
        elif cn_misc_flag == 1:       ## flag to compare non-cloudnet model data
            cn_misc_out_dir = '12_u-br210_RA1M_CASIM/OUT_R0/'
    elif platform == 'JASMIN':
        cn_um_out_dir = 'cloud-fraction-metum-grid/2018/'
        cn_ifs_out_dir = 'cloud-fraction-ecmwf-grid/2018/'
        cn_obs_out_dir = cn_ifs_out_dir
        if cn_misc_flag == 0:       ## flag to compare cloudnet model data
            cn_misc_out_dir = cn_um_out_dir
        elif cn_misc_flag == 1:       ## flag to compare non-cloudnet model data
            cn_misc_out_dir = '12_u-br210_RA1M_CASIM/OUT_R0/'

    ######## lwc-adiabatic-metum-grid/2018/
    ########             -> liquid water content derived using measurements averaged on to model grid
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
            filename_um3 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            ifs_flag = False
        print (filename_um1)
        print (filename_um2)
        print (filename_um3)
        print ('')

        print ('Now load cloudnet data:')
        if cn_um_out_dir[-31:-6] == 'cloud-fraction-metum-grid':
            cn_out_dir = 'cloud-fraction-metum-grid'
        elif cn_um_out_dir[-27:-6] == 'lwc-scaled-metum-grid':
            cn_out_dir = 'lwc-scaled-metum-grid'
        elif cn_um_out_dir[-24:-6] == 'iwc-Z-T-metum-grid':
            out_dir = 'iwc-Z-T-metum-grid'
        cn_filename_um = cn_um_dir + cn_um_out_dir + names[i] + cn_out_dir + '.nc'
        cn_filename_ifs = cn_ifs_dir + cn_ifs_out_dir + names[i] + cn_ifs_out_dir[:-6] + '.nc'
        cn_filename_obs = cn_obs_dir + cn_obs_out_dir + names[i] + cn_obs_out_dir[:-6] + '.nc'
        if cn_misc_flag == 1: cn_filename_misc = cn_misc_dir + cn_misc_out_dir + names[i] + 'metum.nc'
        if cn_misc_flag == 0: cn_filename_misc = cn_misc_dir + cn_misc_out_dir + names[i] + cn_out_dir + '.nc'
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
        if ifs_flag: var_list3 = ['temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
            'sfc_temp_2m','flx_ls_rain','flx_conv_rain','flx_ls_snow','q','pressure','sfc_bl_height']
        if not ifs_flag: var_list3 = var_list1

        if names[i] in moccha_missing_files:        ### NOTE THIS WON'T WORK IF IT'S THE FIRST FILE THAT'S MISSING!!
            print ('File not available...')
            print ('***Filling arrays with nans***')
            print (str(doy[i]))
            nanarray = np.zeros(24)
            nanarray[:] = np.nan
            time_um1 = np.append(time_um1, nanarray)
            time_um2 = np.append(time_um2, nanarray)
            if ifs_flag: time_um3 = np.append(time_um3, nanarray)
            if not ifs_flag: time_um3 = np.append(time_um3, nanarray)
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
                np.save('testing', data3)
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
            cn_nc1 = Dataset(cn_filename_um,'r')
            cn_nc2 = Dataset(cn_filename_ifs,'r')
            if cn_misc_flag != -1: cn_nc3 = Dataset(cn_filename_misc,'r')
            cn_nc4 = Dataset(cn_filename_obs,'r')

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
                data3['height'] = nc3.variables['height'][:]
                data3['flx_height'] = nc3.variables['flx_height'][:]

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
    ### ----------------------------------------------------------------------------------------------------------------
            ### --------------------------------------------------------------------
            ###     LOAD UM CLOUDNET DIAGS INTO DICTIONARY
            ### --------------------------------------------------------------------
            #### LOAD IN SPECIFIC DIAGNOSTICS
            if cn_out_dir == 'cloud-fraction-metum-grid':
                cn_var_list = ['height','Cv','model_Cv_filtered','model_temperature']   ### time always read in separately
            elif cn_out_dir == 'lwc-scaled-metum-grid':
                cn_var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
            elif cn_out_dir == 'iwc-Z-T-metum-grid':
                cn_var_list = ['height','iwc','model_iwc_filtered']   ### time always read in separately

            ###     LOOP OVER TIME DUMP
            if i == 0:
                um_data = {}
                # um_data1d = {}
                if month_flag == -1:
                    time_um = doy[i] + ((cn_nc1.variables['time'][:])/24.0)
                else:
                    time_um = float(names[i][6:8]) + ((cn_nc1.variables['time'][:])/24.0)
                for j in range(0,len(cn_var_list)):
                    if np.ndim(cn_nc1.variables[cn_var_list[j]]) == 1:  # 1d timeseries only
                        um_data[cn_var_list[j]] = cn_nc1.variables[cn_var_list[j]][:]
                    else:                                   # 2d column um_data
                        um_data[cn_var_list[j]] = cn_nc1.variables[cn_var_list[j]][:]
            else:
                if month_flag == -1:
                    time_um = np.append(time_um, doy[i] + ((cn_nc1.variables['time'][:])/24.0))
                else:
                    time_um = np.append(time_um,float(cn_filename_um[-16:-14]) + ((cn_nc1.variables['time'][:])/24.0))
                print (um_data)
                for j in range(0,len(cn_var_list)):
                    # print 'j = ' + str(j)
                    if np.ndim(cn_nc1.variables[cn_var_list[j]]) == 1:
                        um_data[cn_var_list[j]] = np.append(um_data[cn_var_list[j]],cn_nc1.variables[cn_var_list[j]][:])
                    else:
                        um_data[cn_var_list[j]] = np.append(um_data[cn_var_list[j]],cn_nc1.variables[cn_var_list[j]][:],0)
            cn_nc1.close()

            ### --------------------------------------------------------------------
            ### LOAD IN IFS DATA INTO DICTIONARY
            ### --------------------------------------------------------------------
            if cn_ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
                cn_var_list = ['height','Cv','model_snow_Cv_filtered','model_temperature']   ### time always read in separately
            elif cn_ifs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
                cn_var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
            elif cn_ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
                cn_var_list = ['height','iwc','model_snow_iwc_filtered','model_iwc_filtered']   ### time always read in separately

            ###     LOOP OVER TIME DUMP
            if i == 0:
                ifs_data = {}
                # ifs_data1d = {}
                if month_flag == -1:
                    time_ifs = doy[i] + ((cn_nc2.variables['time'][:])/24.0)
                else:
                    time_ifs = float(names[i][6:8]) + ((cn_nc2.variables['time'][:])/24.0)
                for j in range(0,len(cn_var_list)):
                    if np.ndim(cn_nc2.variables[cn_var_list[j]]) == 1:  # 1d timeseries only
                        ifs_data[cn_var_list[j]] = cn_nc2.variables[cn_var_list[j]][:]
                    else:                                   # 2d column um_data
                        ifs_data[cn_var_list[j]] = cn_nc2.variables[cn_var_list[j]][:]
            else:
                if month_flag == -1:
                    time_ifs = np.append(time_ifs, doy[i] + ((cn_nc2.variables['time'][:])/24.0))
                else:
                    time_ifs = np.append(time_ifs,float(cn_filename_ifs[-16:-14]) + ((cn_nc2.variables['time'][:])/24.0))
                print (ifs_data)
                for j in range(0,len(cn_var_list)):
                    ## ONLY WANT COLUMN VARIABLES - IGNORE TIMESERIES FOR NOW
                    # print 'j = ' + str(j)
                    if np.ndim(cn_nc2.variables[cn_var_list[j]]) == 1:
                        ifs_data[cn_var_list[j]] = np.append(ifs_data[cn_var_list[j]],cn_nc2.variables[cn_var_list[j]][:])
                    else:
                        ifs_data[cn_var_list[j]] = np.append(ifs_data[cn_var_list[j]],cn_nc2.variables[cn_var_list[j]][:],0)
            cn_nc2.close()

            ### -------------------------------------------------------------------------
            ###     LOAD IN MISC DATA INTO DICTIONARY IF COMPARING
            ###             Only load in what variables are needed based on IFS file chosen
            ### -------------------------------------------------------------------------
            if cn_misc_flag == -1:
                continue
            elif cn_misc_flag == 1:
                if cn_ifs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
                    cn_var_list = ['cloud_fraction','temperature']   ### time always read in separately
                elif cn_ifs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
                    cn_var_list = ['qliq']   ### time always read in separately
                elif cn_ifs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
                    cn_var_list = ['qice']   ### time always read in separately
            elif cn_misc_flag == 0:
                if cn_out_dir == 'cloud-fraction-metum-grid':
                    cn_var_list = ['height','Cv','model_Cv_filtered','model_temperature']   ### time always read in separately
                elif cn_out_dir == 'lwc-scaled-metum-grid':
                    cn_var_list = ['height','lwc','model_lwc','model_lwp']   ### time always read in separately
                elif cn_out_dir == 'iwc-Z-T-metum-grid':
                    cn_var_list = ['height','iwc','model_iwc_filtered']   ### time always read in separately

            print ('')
            print ('misc file variable list is:')
            print (cn_var_list)
            print ('')

            if i == 0:
                misc_data = {}
                # misc_data1d = {}
                if month_flag == -1:
                    if cn_misc_flag == 1:
                        time_misc = doy[i] + ((cn_nc3.variables['forecast_time'][:])/24.0)
                        misc_data['height'] = cn_nc3.variables['height'][:]
                    if cn_misc_flag == 0: time_misc = doy[i] + ((cn_nc3.variables['time'][:])/24.0)
                else:
                    if cn_misc_flag == 1: time_misc = float(names[i][6:8]) + ((cn_nc3.variables['forecast_time'][:])/24.0)
                    if cn_misc_flag == 0: time_misc = float(names[i][6:8]) + ((cn_nc3.variables['time'][:])/24.0)
                for j in range(0,len(cn_var_list)):
                    if np.ndim(cn_nc3.variables[cn_var_list[j]]) == 1:  # 1d timeseries only
                        misc_data[cn_var_list[j]] = cn_nc3.variables[cn_var_list[j]][:]
                    else:                                   # 2d column um_data
                        misc_data[cn_var_list[j]] = cn_nc3.variables[cn_var_list[j]][:]
            else:
                if month_flag == -1:
                    if cn_misc_flag == 1: time_misc = np.append(time_misc, doy[i] + ((cn_nc3.variables['forecast_time'][:])/24.0))
                    if cn_misc_flag == 0: time_misc = np.append(time_misc, doy[i] + ((cn_nc3.variables['time'][:])/24.0))
                else:
                    if cn_misc_flag == 1: time_misc = np.append(time_misc,float(cn_filename_misc[-16:-14]) + ((cn_nc3.variables['forecast_time'][:])/24.0))
                    if cn_misc_flag == 0: time_misc = np.append(time_misc,float(cn_filename_misc[-16:-14]) + ((cn_nc3.variables['time'][:])/24.0))
                print (misc_data)
                for j in range(0,len(cn_var_list)):
                    # print 'j = ' + str(j)
                    if np.ndim(cn_nc3.variables[cn_var_list[j]]) == 1:
                        misc_data[cn_var_list[j]] = np.append(misc_data[cn_var_list[j]],cn_nc3.variables[cn_var_list[j]][:])
                    # elif var_list[j] == 'height':#np.sum(nc3.variables[var_list[j]].shape) == 71:
                    #     continue
                    else:
                        misc_data[cn_var_list[j]] = np.append(misc_data[cn_var_list[j]],cn_nc3.variables[cn_var_list[j]][:],0)
            cn_nc3.close()

            ### -------------------------------------------------------------------------
            ###     LOAD IN OBS DATA
            ###             Only load in what variables are needed based on IFS file chosen
            ### -------------------------------------------------------------------------
            if cn_obs_out_dir[:-6] == 'cloud-fraction-ecmwf-grid':
                cn_var_list = ['height','Cv']   ### time always read in separately
            elif cn_obs_out_dir[:-6] == 'lwc-scaled-ecmwf-grid':
                cn_var_list = ['height','lwc','lwp']   ### time always read in separately
            elif cn_obs_out_dir[:-6] == 'iwc-Z-T-ecmwf-grid':
                cn_var_list = ['height','iwc']   ### time always read in separately

            if i == 0:
                obs_data = {}
                # misc_data1d = {}
                if month_flag == -1:
                    time_obs = doy[i] + ((cn_nc4.variables['time'][:])/24.0)
                else:
                    time_obs = float(names[i][6:8]) + ((cn_nc4.variables['time'][:])/24.0)
                for j in range(0,len(cn_var_list)):
                    if np.ndim(cn_nc4.variables[cn_var_list[j]]) == 1:  # 1d timeseries only
                        obs_data[cn_var_list[j]] = cn_nc4.variables[cn_var_list[j]][:]
                    else:                                   # 2d column um_data
                        obs_data[cn_var_list[j]] = cn_nc4.variables[cn_var_list[j]][:]
            else:
                if month_flag == -1:
                    time_obs = np.append(time_obs, doy[i] + ((cn_nc4.variables['time'][:])/24.0))
                else:
                    time_obs = np.append(time_obs,float(cn_filename_obs[-16:-14]) + ((cn_nc4.variables['time'][:])/24.0))
                print (obs_data)
                for j in range(0,len(cn_var_list)):
                    # print 'j = ' + str(j)
                    if np.ndim(cn_nc4.variables[cn_var_list[j]]) == 1:
                        obs_data[cn_var_list[j]] = np.append(obs_data[cn_var_list[j]],cn_nc4.variables[cn_var_list[j]][:])
                    elif np.sum(cn_nc4.variables[cn_var_list[j]].shape) == 71:
                        continue
                    else:
                        obs_data[cn_var_list[j]] = np.append(obs_data[cn_var_list[j]],cn_nc4.variables[cn_var_list[j]][:],0)
            cn_nc4.close()

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
    ## load calculated model inversion heights
    #################################################################
    print ('Load calculated model inversion heights...')
    data1['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_RA2M_inversion_results.mat')
    data2['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_CASIM-100_inversion_results.mat')
    data3['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/ECMWF_IFS_inversion_results.mat')

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '12_u-br210': label1 = 'UM_CASIM-AeroProf'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T'
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'

    label2 = 'undefined_label'
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
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['sondes'])

    ### cloudnet
    np.save('working_um_data', um_data)
    np.save('working_ifs_data', ifs_data)
    if cn_misc_flag != -1: np.save('working_misc_data', misc_data)


###################################################################################################################
###################################################################################################################
################################################ FIGURES ##########################################################
###################################################################################################################
###################################################################################################################

    # -------------------------------------------------------------
    # Test cloudnet plot: Plot Cv statistics from drift period
    # -------------------------------------------------------------
    # figure = plot_CvProfiles(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)

    # -------------------------------------------------------------
    # Test model plot: flux timeseries
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # cloud properties scaled by BL depth
    # -------------------------------------------------------------
    figure = plot_scaledBL(data1, data2, data3,um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    ### model/measurement data
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['sondes'])

    ### cloudnet
    np.save('working_um_data', um_data)
    np.save('working_ifs_data', ifs_data)
    if cn_misc_flag != -1: np.save('working_misc_data', misc_data)

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
