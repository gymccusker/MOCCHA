###
###
### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
###
###

from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import diags_varnames as varnames
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
    if out_dir4 == 'OUT_25H/':
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
    plt.plot(obs['ice_station']['mday'],obs['ice_station']['tafluxB'], 'r.')
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

    # if month_flag == -1:
        # if out_dir1[:20] == '5_u-bl661_RA1M_CASIM':
            # if out_dir2[:20] == '6_u-bm410_RA1M_CASIM':
            #     if out_dir4 == 'OUT_25H/':
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] + '_oden_metum_ifs_casim_TSa.svg'
            #     else:
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_casim_TSa.png'
            # elif out_dir2[:9] == '4_u-bg410':
            #     if out_dir4 == 'OUT_25H/':
            #         fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:20] + '_oden_metum_ifs_casim_TSa.svg'
            # else:
            #     fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_oden_metum_casim_TSa.svg'
        # if out_dir2[:20] == '5_u-bl661_RA1M_CASIM':
        #     if out_dir4[:20] == '6_u-bm410_RA1M_CASIM':
        #         fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_' + out_dir4[:9] + '_oden_metum_casim-100_200_TSa.png'
        #     elif out_dir4 == 'OUT_25H/':
        #         fileout = '../FIGS/comparisons/' + out_dir2[:20] + '_oden_metum_ifs_casim-100_TSa.svg'
        # elif out_dir1[:18] == '4_u-bg610_RA2M_CON':
        #     fileout = '../FIGS/comparisons/' + out_dir1[:9] + '_' + out_dir2[:9] +'_oden_metum_casim-100_TSa.svg'
    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_ra2t_ifs_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Plotting BL depth timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	set up figure properties
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
    plt.figure(figsize=(10,5))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.94, left = 0.12,
            hspace = 0.4, wspace = 0.15)

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ## create figure and axes instances
    #################################################################

    ax = plt.gca()
    plt.plot(data1['time'], data1['bl_depth'], color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['bl_depth'], color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_bl_height'], color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['bl_depth'], color = 'darkorange', label = label3)
    plt.legend()
    plt.title('BL_depth [m]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')
    plt.ylabel('Z [m]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLDepth_timeseries_oden_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_BLType(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    ###################################
    ## Plot occurrences of BL type from UM
    ###################################

    print ('******')
    print ('')
    print ('Plotting BL type histograms:')
    print ('')

    ##################################################
    ##################################################
    #### 	set up figure properties
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
    plt.figure(figsize=(9,8))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.5, left = 0.08,
            hspace = 0.4, wspace = 0.10)

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ## split data into relative contributions from each BL type
    #################################################################

    data1['histogram_n'], data1['histogram_bins'] = np.histogram(data1['bl_type'], bins = np.arange(1,9))
    data2['histogram_n'], data2['histogram_bins'] = np.histogram(data2['bl_type'], bins = np.arange(1,9))

    ### calculate total number of occurrences for normalised charts
    total1 = float(np.sum(data1['histogram_n']))
    total2 = float(np.sum(data2['histogram_n']))

    ### create arranys of instances for e.g. type I, type II etc. between each simulation
    types = {}
    for i in range(0,7): types[i+1] = [data1['histogram_n'][i]/total1, data2['histogram_n'][i]/total2]

    #### list of BL types from UM documentation
    doc = ['1: Stable BL', '2: Sc over stable SL', '3: Well-mixed BL', '4: Unstable BL, dSc not o/Cu',
        '5: dSc o/Cu', '6: Cu-capped BL', '7: Shear-dom unstable BL']

    # Type I: Stable boundary layer (with or without cloud)
    # Type II: Boundary layer with stratocumulus over a stable near-surface layer
    # Type III: Well mixed boundary layer
    # Type IV: Unstable boundary layer with a DSC layer not over cumulus
    # Type V: Boundary layer with a DSC layer over cumulus
    # Type VI: Cumulus-capped boundary layer
    # Type VII: Shear-dominated unstable layer

    #################################################################
    ## create figure and axes instances
    #################################################################
    ax = plt.gca()

    plt.bar([0,1], types[1], label = doc[0])
    plt.bar([0,1], types[2], bottom = types[1], label = doc[1]); bars = np.add(types[1], types[2]).tolist()
    for i in range(3,8):
        # print(i)
        if i == 4:
            plt.bar([0,1], types[i], bottom = bars, hatch = '/', label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
        else:
            plt.bar([0,1], types[i], bottom = bars, label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
    plt.xticks([0,1], [label1, label2])
    plt.legend(bbox_to_anchor=(1.1, 0.3, 1., .102), loc=4, ncol=1)
    plt.title('BL type occurrences (normalised)')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLType_oden_metum_casim-100.svg'
    plt.savefig(fileout)
    plt.show()


    #### ------------------------------------------------------------------
    #### ------------------------------------------------------------------
    ####        SEASONAL SPLIT
    #### ------------------------------------------------------------------
    #### ------------------------------------------------------------------


    melt = np.where(data1['time'] < 240.0)
    freeze = np.where(data1['time'] >= 240.0)

    #################################################################
    ## split data into relative contributions from each BL type
    #################################################################

    data1['histogram_nMelt'], data1['histogram_bins'] = np.histogram(data1['bl_type'][melt], bins = np.arange(1,9))
    data2['histogram_nMelt'], data2['histogram_bins'] = np.histogram(data2['bl_type'][melt], bins = np.arange(1,9))

    data1['histogram_nFreeze'], data1['histogram_bins'] = np.histogram(data1['bl_type'][freeze], bins = np.arange(1,9))
    data2['histogram_nFreeze'], data2['histogram_bins'] = np.histogram(data2['bl_type'][freeze], bins = np.arange(1,9))

    ### calculate total number of occurrences for normalised charts
    total1melt = float(np.sum(data1['histogram_nMelt']))
    total2melt = float(np.sum(data2['histogram_nMelt']))
    total1freeze = float(np.sum(data1['histogram_nFreeze']))
    total2freeze = float(np.sum(data2['histogram_nFreeze']))


    ### create arranys of instances for e.g. type I, type II etc. between each simulation
    types = {}
    for i in range(0,7): types[i+1] = [data1['histogram_nMelt'][i]/total1melt, data1['histogram_nFreeze'][i]/total1freeze,
                                        data2['histogram_nMelt'][i]/total2melt, data2['histogram_nFreeze'][i]/total2freeze]

    #### list of BL types from UM documentation
    doc = ['1: Stable BL', '2: Sc over stable SL', '3: Well-mixed BL', '4: Unstable BL, dSc not o/Cu',
        '5: dSc o/Cu', '6: Cu-capped BL', '7: Shear-dom unstable BL']

    # Type I: Stable boundary layer (with or without cloud)
    # Type II: Boundary layer with stratocumulus over a stable near-surface layer
    # Type III: Well mixed boundary layer
    # Type IV: Unstable boundary layer with a DSC layer not over cumulus
    # Type V: Boundary layer with a DSC layer over cumulus
    # Type VI: Cumulus-capped boundary layer
    # Type VII: Shear-dominated unstable layer

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.figure(figsize=(12,7))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.68, left = 0.05,
            hspace = 0.4, wspace = 0.10)
    ax = plt.gca()

    plt.bar(np.arange(0,4), types[1], label = doc[0])
    plt.bar(np.arange(0,4), types[2], bottom = types[1], label = doc[1]); bars = np.add(types[1], types[2]).tolist()
    for i in range(3,8):
        # print(i)
        if i == 4:
            plt.bar(np.arange(0,4), types[i], bottom = bars, hatch = '/', label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
        else:
            plt.bar(np.arange(0,4), types[i], bottom = bars, label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
    plt.xticks(np.arange(0,4), [label1 + '\nMelt', label1 + '\nFreeze', label2 + '\nMelt', label2 + '\nFreeze'])
    plt.legend(bbox_to_anchor=(0.5, 0.3, 1., .102), loc=4, ncol=1)
    plt.title('BL type occurrences (normalised)')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLType_oden_metum_casim-100_splitSeason.svg'
    plt.savefig(fileout)
    plt.show()

def thetaVL_Stats(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm


    print ('******')
    print ('')
    print ('Plotting radiosonde and model theta profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if out_dir4 == 'OUT_25H/':
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data1['pressure'][data1['pressure'] == -9999] = np.nan
    data2['pressure'][data2['pressure'] == -9999] = np.nan
    data3['pressure'][data3['pressure'] <= 0] = np.nan
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### calculate equivalent potential temperature
    #### ---------------------------------------------------------------
    data1['theta'], data1['thetaL'], data1['thetaVL'] = calcThetaVL(data1['temperature'], data1['pressure'], data1['q'],
                data1['qliq'], data1['qice'], data1['time'], data1['height'])
    data2['theta'], data2['thetaL'], data2['thetaVL'] = calcThetaVL(data2['temperature'], data2['pressure'], data2['q'],
                data2['qliq'], data2['qice'], data2['time'], data2['height'])

    # obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
    #     np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3,
    #     obs['sondes']['doy'], obs['sondes']['gpsaltitude'])
    #
    # obs['sondes']['theta'] = np.transpose(obs['sondes']['theta'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data2',data2)
    # np.save('working_dataObs',obs['sondes'])

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    # print ('...')
    # print ('Re-gridding sonde and ifs data...')
    # print ('')
    # data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'theta')
    # print ('')
    # print ('Done!')

def reGrid_Sondes(data1, data2, data3, obs, doy, var):

    from scipy.interpolate import interp1d

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

    ###
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

    #### ---------------------------------------------------------------
    ### build list of variables names wrt input data [OBS, UM, CASIM, IFS]
    #### ---------------------------------------------------------------
    if var == 'temp':
        varlist = ['temperature','temperature','temperature','temperature']
    elif var == 'thetaE':
        # varlist = ['epottemp','thetaE','thetaE','thetaE']     # use sonde file's epottemp
        varlist = ['thetaE','thetaE','thetaE','thetaE']         # use sonde calculated thetaE
    elif var == 'theta':
        # varlist = ['pottemp','theta','theta','theta']     # use sonde file's pottemp
        varlist = ['theta','theta','theta','theta']         # use sonde calculated theta
    elif var == 'q':
        varlist = ['mr','q','q','q']

    #### ---------------------------------------------------------------
    #### save hourly temperature model profiles (using the ii index defined by the time indices)
    #### ---------------------------------------------------------------
    data1[var + '_hrly'] = np.squeeze(data1[varlist[1]][ii,:])
    data2[var + '_hrly'] = np.squeeze(data2[varlist[2]][ii,:])
    data3[var + '_hrly'] = np.squeeze(data3[varlist[3]][ii,:])

    #### ---------------------------------------------------------------
    #### explicitly save 6-hourly temperature model profiles and time binning for ease
    #### ---------------------------------------------------------------
    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_6hrly'] = data1['time_hrly'][::6]
    data2['time_6hrly'] = data2['time_hrly'][::6]
    data3['time_6hrly'] = data3['time_hrly'][::6]
    data1[var + '_6hrly'] = data1[var + '_hrly'][::6]
    data2[var + '_6hrly'] = data2[var + '_hrly'][::6]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6]

    #### ---------------------------------------------------------------
    #### index to only look at altitudes <10km
    #### ---------------------------------------------------------------
    iTim = 0        ### initialised
    iObs = np.where(obs['sondes']['gpsaltitude'][:,iTim] <= 11000)
    iUM = np.where(data1['height'] <= 11000)
    iIFS = np.where(data3['height'][iTim,:] <= 11000)

    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    print ('')
    print ('Defining IFS temperature profile as a function:')
    print ('using ifs.height[0,:] to define temperature profiles...')
    data3[var + '_hrly_UM'] = np.zeros([np.size(data3['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(data3['time_hrly'],0)):
        fnct_IFS = interp1d(np.squeeze(data3['height'][0,iIFS]), np.squeeze(data3[var + '_hrly'][iTim,iIFS]))
        data3[var + '_hrly_UM'][iTim,:] = fnct_IFS(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('IFS(UM Grid) function worked!')
    print (var + ' IFS data now on UM vertical grid')
    print ('*****')

    #### INTERPOLATION TESTING:
    # print data3['temp_hrly_UM'].shape
    # print data3['time_hrly'][::6].shape
    # print data1['temp_hrly'][:,iUM[0][3:]].shape
    # print data1['time_hrly'][::6].shape
    # plt.plot(data3['temp_hrly_UM'][10,:],data1['height'][iUM[0][2:]])
    # plt.plot(np.squeeze(data3['temp_hrly'][10,iIFS]),np.squeeze(data3['height'][10,iIFS]))
    # plt.show()

    print ('')
    print ('Defining Sonde temperature profile as a function:')
    obs['sondes'][var + '_allSondes_UM'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(obs['sondes']['doy'],0)):
        # print 'iTim = ', str(iTim)
        fnct_Obs = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
        obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_Obs(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('Sonde(UM Grid) function worked!')
    print ('All ' + var + ' sonde data now on UM vertical grid.')
    print ('*****')

    #### INTERPOLATION TESTING:
    # print obs['sondes']['temp_hrly_UM'].shape
    # print obs['sondes']['doy'].shape
    # print obs['sondes']['temp_hrly_UM']
    # plt.plot(obs['sondes']['temp_hrly_UM'],data1['height'][iUM[0][2:]])
    # plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,iTim]),np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]))
    # plt.show()

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['sondes']['doy'] >= 225.9, obs['sondes']['doy'] <= 258.0))

    ### save in dict for ease
    obs['sondes']['doy_drift'] = obs['sondes']['doy'][drift]

    #### ---------------------------------------------------------------
    #### make some dictionary assignments for use later
    #### ---------------------------------------------------------------
    data1['universal_height'] = data1['height'][iUM[0][3:]].data
    data1['universal_height_UMindex'] = iUM[0][3:]

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])

    return data1, data2, data3, obs, drift

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
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        misc_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        # misc_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/'
    out_dir2 = '5_u-bl661_RA1M_CASIM/OUT_R0/'
    # out_dir3 = 'MET_DATA/'
    out_dir4 = 'OUT_25H/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R0/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param

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
    print ('Loading observations:')
            # -------------------------------------------------------------
            # Which file does what?
            # -------------------------------------------------------------
            #### ice station: net LW / net SW
                    #### obs['ice_station']/mast_radiation_30min_v2.3.mat
            #### obs['foremast']:
                    #### obs['foremast']/ACAS_AO2018_obs['foremast']_30min_v2_0.nc
            #### 7th deck: temperature, surface temperature, RH, downwelling SW, downwelling LW
                    #### 7thDeck/ACAS_AO2018_WX_30min_v2_0.nc

    obs = {}

    print ('Load temporary ice station data from Jutta...')
    obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_wTemp1p5m.nc','r')

    print ('Load ice station data from Jutta...')
    obs['ice_station'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')
            #### mast_radiation_30min_v2.3.mat
            #### flux30_trhwxrel.mat

    print ('Load radiosonde data from Jutta...')
    obs['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print ('Load foremast data from John...')
    obs['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print ('Load 7th deck weather station data from John...')
    obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print ('...')

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin file read in at ' + time.strftime("%c"))
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
            '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180813_oden_']

    doy = np.arange(226,258)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)          ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):
        filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
        filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
        if out_dir4 == 'OUT_25H/':
            print( '***IFS being compared***')
            ifs_flag = True
            filename_um3 = misc_root_dir + out_dir4 + names[i] + 'ecmwf.nc'
        else:
            print ('***IFS NOT being compared***')
            filename_um3 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            ifs_flag = False
        print (filename_um1)
        print (filename_um2)
        print (filename_um3)
        print ('')

        #### LOAD CUBE
        print( 'Loading first run diagnostics:')
        nc1 = Dataset(filename_um1,'r')
        print ('...')
        print( 'Loading second run diagnostics:')
        nc2 = Dataset(filename_um2,'r')
        print ('...')
        print ('Loading third run diagnostics:')
        nc3 = Dataset(filename_um3,'r')
        print ('...')
        # -------------------------------------------------------------
        # print 'i = ' + str(i)
        print ('')

        #### LOAD IN SPECIFIC DIAGNOSTICS
        # if out_dir == '4_u-bg610_RA2M_CON/OUT_R1/':
        var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux','latent_heat_flux',
            'temp_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice']
        var_list2 = var_list1
        if ifs_flag: var_list3 = ['height','temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
            'sfc_temp_2m','flx_ls_rain','flx_conv_rain','flx_ls_snow','q','pressure','sfc_bl_height']
        if not ifs_flag: var_list3 = var_list1

        if i == 0:
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
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
            nc3.close()
        else:
            if month_flag == -1:
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
                    data1[var_list1[j]] = np.append(data1[var_list1[j]].data,nc1.variables[var_list1[j]][:])
                elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                    data1[var_list1[j]] = np.append(data1[var_list1[j]].data,nc1.variables[var_list1[j]][:],0)
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            for j in range(0,len(var_list2)):
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) == 1:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]].data,nc2.variables[var_list2[j]][:])
                elif np.ndim(nc2.variables[var_list2[j]]) == 2:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]].data,nc2.variables[var_list2[j]][:],0)
            nc2.close()
            ## ------------------
            #### um3 / ifs
            ## ------------------
            for j in range(0,len(var_list3)):
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]].data,nc3.variables[var_list3[j]][:])
                elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]].data,nc3.variables[var_list3[j]][:],0)
            nc3.close()

    #################################################################
    ## save time to dictionary now we're not looping over all diags anymore
    #################################################################
    data1['time'] = time_um1
    data2['time'] = time_um2
    data3['time'] = time_um3

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T'
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'

    label2 = 'undefined_label'
    if out_dir2[:10] == '11_u-bq798': label2 = 'UM_CASIM-100_Meyers'
    if out_dir2[:10] == '10_u-bq791': label2 = 'UM_CASIM-100_Fletcher'
    if out_dir2[:9] == '8_u-bp738': label2 = 'UM_ERAI-GLM'
    if out_dir2[:9] == '7_u-bn068': label2 = 'UM_RA2T'
    if out_dir2[:9] == '6_u-bm410': label2 = 'UM_CASIM-200'
    if out_dir2[:9] == '5_u-bl661': label2 = 'UM_CASIM-100'
    if out_dir2[:9] == '4_u-bg610': label2 = 'UM_RA2M'

    label3 = 'undefined_label'
    if out_dir4 == 'OUT_25H/': label3 = 'ECMWF_IFS'
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
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['sondes'])

    # -------------------------------------------------------------
    # Plot combined column data (5x2 timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_multidate_TS(timem, data, cube, month_flag, missing_files, out_dir)
                ### doesn't matter which cube, just needed for dim_coords

    # -------------------------------------------------------------
    # Plot combined CASIM column data (4x3 timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_multidate_casim_TS(timem, data, cube, month_flag, missing_files, out_dir)
                ### doesn't matter which cube, just needed for dim_coords

    # -------------------------------------------------------------
    # Plot combined timeseries as lineplot
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # figure = plot_line_BLDepth(time_um1, time_um2, data1, data2, cube_um1, cube_um2, month_flag,
    #             missing_files, out_dir1, obs, doy)

    # figure = plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3,
    #     month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # figure = plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # Plot paper figures
    # -------------------------------------------------------------
    # figure = plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_paperRadiation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_Precipitation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_BLType(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesTemperature(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesQ(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    figure = thetaVL_Stats(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesTheta(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_ERAI_GLM(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)


    # -------------------------------------------------------------
    # Plot data (5x2 monthly timeseries)
    # -------------------------------------------------------------
    # figure = plot_multicontour_TS(cube, filename, out_dir)


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
