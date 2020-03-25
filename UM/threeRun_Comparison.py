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

def trackShip(data):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==25,data.values[:,1]==8),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print( '******')
    print( '')
    print( 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    print( 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print( '')

    return trackShip_index

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
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

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

    # time_um1 = data1['time'][:]
    # time_um2 = data2['time'][:]
    # time_um3 = data3['time'][:]

    #################################################################
    ## create figure and axes instances
    #################################################################

    plt.subplot(211)
    ax = plt.gca()
    plt.plot(data1['time'][:], data1['temp_1.5m'].data - 273.15, color = 'steelblue', label = label1)
    plt.plot(data2['time'][:], data2['temp_1.5m'].data - 273.15, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_temp_2m'].data - 273.15, color = 'darkorange', label =  label3)
    else:
        plt.plot(data3['time'][:], data3['temp_1.5m'].data - 273.15, color = 'darkorange')#, label = '2m')
    plt.plot(time_temp,obs['obs_temp'].variables['Tice'][:] - 273.15, color = 'black', label = 'Observations')
    plt.legend()
    plt.title('Temperature [$^{o}C$]')
    plt.ylim([263 - 273,275 - 273])
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
    plt.plot(data1['time'][:], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'][:], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'][:], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.title('Net SW radiation [W/m2]')
    plt.ylim([0,100])
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
    # # data['surface_net_LW_radiation'].data[data['surface_net_LW_radiation'].data == 0] = np.nan
    # plt.plot(data1['time'][:], data1['surface_net_LW_radiation'].data, color = 'steelblue', label = label1)
    # plt.plot(data2['time'][:], data2['surface_net_LW_radiation'].data, color = 'forestgreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'][:], data3['sfc_net_lw'].data, color = 'darkorange')
    # else:
    #     plt.plot(data3['time'][:], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    # plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'Observations')
    # # plt.legend()
    # plt.title('Net LW radiation [W/m2]')
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

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        if out_dir2[0:20] == '6_u-bm410_RA1M_CASIM':
            fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-200_tempoC_SW.png'
        elif out_dir2[0:20] == '5_u-bl661_RA1M_CASIM':
            if ifs_flag == True:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
            else:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-100_tempoC_SW.png'
        elif out_dir2[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_tempoC_SW.png'
        else:
            fileout = '../FIGS/comparisons/' + out_dir2[:18] + '_oden_metum_tempoC_SW.png'
    # print 'Saving as: ' + fileout
    fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
    plt.savefig(fileout)#, dpi=400)
    plt.show()

def plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
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
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'firebrick')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'firebrick', label = label2)
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
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'firebrick')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'firebrick')# * -1.0)
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
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'firebrick')# * -1.0)
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

    fileout = '../FIGS/comparisons/CASIM-NiceTest_C86-F62-M92_DOY243-249.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
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
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'purple')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'darkorange', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'darkorange')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
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
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'purple')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'purple')# * -1.0)
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
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'purple')# * -1.0)
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

    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_ra2t_ifs_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_subSect(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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
    print ('Plotting combined timeseries and PDFs of radiation terms:')
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
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

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

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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

    #### add override for data2 to allow 24h data to be used for testing purposes
    if out_dir2[-4:] == '24h/':
        data2['surface_net_LW_radiation'][data2['surface_net_LW_radiation'] == 0] = np.nan
        data2['surface_net_SW_radiation'][data2['surface_net_SW_radiation'] == 0] = np.nan

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(12,10))

    ax  = fig.add_axes([0.07,0.7,0.55,0.22])   # left, bottom, width, height
    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]
    ax = plt.gca()
    yA = [-65, 85]
    # plt.plot([240.0,240.0],[yA[0],yA[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.legend(bbox_to_anchor=(-0.11, 0.65, 1., .102), loc=4, ncol=2)
    plt.ylim([-60,80])

    ax  = fig.add_axes([0.07,0.4,0.55,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    # plt.plot([240.0,240.0],[yB[0],yB[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    # plt.legend()
    ax.set_xlim([doy[0],doy[-1]])
    plt.ylim([-3,120])

    ax  = fig.add_axes([0.07,0.1,0.55,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yC = [-90, 10]
    # plt.plot([240.0,240.0],[yC[0],yC[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'purple')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('surface_net_LW_radiation [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')
    plt.ylim([-90,5])

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)

    #### only compare obs over model dates available:
    ####        all model data share a timestamp

    subSect = np.where(np.logical_and(time_radice > data1['time_hrly'][0], time_radice <= data1['time_hrly'][-1]))

    sw1 = data1['surface_net_SW_radiation'][data1['hrly_flag']]
    lw1 = data1['surface_net_LW_radiation'][data1['hrly_flag']]
    sw2 = data2['surface_net_SW_radiation'][data2['hrly_flag']]
    lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    if ifs_flag == True:
        sw3 = data3['sfc_net_sw'][data3['hrly_flag']]
        lw3 = data3['sfc_net_lw'][data3['hrly_flag']]
    else:
        sw3 = data3['surface_net_SW_radiation'][data3['hrly_flag']]
        lw3 = data3['surface_net_LW_radiation'][data3['hrly_flag']]

    ax  = fig.add_axes([0.7,0.7,0.25,0.22])   # left, bottom, width, height
    yDmax = 0.1
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1 + lw1
    sns.distplot(crf1, hist=False, color="steelblue", kde_kws={"shade": True})
    crf3 = sw3 + lw3
    sns.distplot(crf3, hist=False, color="darkorange", kde_kws={"shade": True})
    crf2 = sw2 + lw2
    sns.distplot(crf2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netLW[subSect] + netSW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(55,0.07), xytext=(55,0.07), fontsize = 14)
    plt.xlabel('CRF [W/m2]')
    plt.xlim([-50,80])
    plt.ylim([0,yDmax])

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.4,0.25,0.22])   # left, bottom, width, height
    yEmax = 0.14
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(sw1, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(sw3, hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(sw2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netSW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(87,0.07), xytext=(87,0.07), fontsize = 14)
    # plt.legend()
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('$SW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.1,0.25,0.22])   # left, bottom, width, height
    yFmax = 0.11
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    sns.distplot(lw1, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(lw3, hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(lw2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netLW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(0,0.14), xytext=(0,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('$LW_{net,surf}$ [W/m2]')


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/CRF_netSW_netLW_line+PDFS_oden_iceStation_' + label1[3:] + '_' + label3 + '_' + label2[3:] + '.svg'
    plt.savefig(fileout)
    plt.show()

def plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting timeseries of turbulent fluxes:')
    print ('')

    ##################################################
    ##################################################
    #### 	SET AXES PROPERTIES
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 18
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### set diagnostic naming flags for if IFS being used
    ### -------------------------------
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### -------------------------------
    ### for reference in figures
    ### -------------------------------
    zeros = np.zeros(len(data2['time']))

    ### -------------------------------
    ### change matlab time for obs
    ### -------------------------------
    time_iceStation = calcTime_Mat2DOY(np.squeeze(obs['ice_station_fluxes']['mday'][:]))

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12))

    ax  = fig.add_axes([0.07,0.55,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1],
        obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    plt.plot(time_iceStation[np.squeeze(obs['ice_station_fluxes']['taflag'][:]==1)],
        np.squeeze(obs['ice_station_fluxes']['taflux'][obs['ice_station_fluxes']['taflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 7, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'darkorange')
    plt.ylim([-20, 40])
    plt.title('sensible_heat_flux [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.1,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1],
        obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    # index = np.logical_and(obs['ice_station_fluxes']['lrflux']>=-30, obs['ice_station_fluxes']['lrflux']<=70)
    plt.plot(time_iceStation[np.squeeze(obs['ice_station_fluxes']['lrflag'][:]==1)],
        np.squeeze(obs['ice_station_fluxes']['lrflux'][obs['ice_station_fluxes']['lrflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 7, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'forestgreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'darkorange')# * -1.0)
    plt.title('latent_heat_flux [W/m2]')
    plt.ylim([-20, 60])
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)
    ax  = fig.add_axes([0.7,0.55,0.27,0.35])   # left, bottom, width, height
    # zerosC = np.zeros(len(data2['time']))
    yCmax = 0.16
    plt.plot([0,0],[0,yCmax],'--', color='lightgrey')
    ##---
    shf1 = data1['sensible_heat_flux'][data1['hrly_flag']].data
    indextaum = np.logical_and(shf1 >= -50, shf1 <= 50)
    sns.distplot(shf1[indextaum], hist=False, color="steelblue", kde_kws={"shade": True}, label = label1)
    ##---
    shf3 = data3['sfc_down_sens_heat_flx'][data3['hrly_flag']].data * -1.0
    indextaifs = np.logical_and(shf3 >= -50, shf3 <= 50)
    sns.distplot(shf3[indextaifs], hist=False, color="darkorange", kde_kws={"shade": True}, label = label3)
    ##---
    shf2 = data2['sensible_heat_flux'][data2['hrly_flag']].data
    indextacasim = np.logical_and(shf2 >= -50, shf2 <= 50)
    sns.distplot(shf2[indextacasim], hist=False, color="forestgreen", kde_kws={"shade": True}, label = label2)
    ##---
    fmst_taflux = obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1]
    indextafmst = np.logical_and(fmst_taflux>=-50, fmst_taflux<=50)
    sns.distplot(fmst_taflux[indextafmst], hist=False, color="black", label = 'Foremast')#, kde_kws={"shade": True}, label = 'Foremast')
    ##---
    taflux = np.squeeze(obs['ice_station_fluxes']['taflux'][obs['ice_station_fluxes']['taflag'][:] == 1])
    indexta = np.logical_and(taflux>=-50, taflux<=50)
    sns.distplot(taflux[indexta], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3}, label = 'Ice_station')
    plt.title('sensible_heat_flux [W/m2]')
    plt.legend()
    plt.xlim([-20,50])
    plt.ylim([0,yCmax])

    ax  = fig.add_axes([0.7,0.1,0.27,0.35])   # left, bottom, width, height
    yDmax = 0.12
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    ##---
    lhf1 = data1['latent_heat_flux'][data1['hrly_flag']].data
    indexlrum = np.logical_and(lhf1 >= -50, lhf1 <= 50)
    sns.distplot(lhf1[indexlrum], hist=False, color="steelblue", kde_kws={"shade": True})
    ##---
    lhf3 = data3['sfc_down_lat_heat_flx'][data3['hrly_flag']].data * -1.0
    indexlrifs = np.logical_and(lhf3 >= -50, lhf3 <= 50)
    sns.distplot(lhf3[indexlrifs], hist=False, color="darkorange", kde_kws={"shade": True})
    ##---
    lhf2 = data2['latent_heat_flux'][data2['hrly_flag']].data
    indexlrcasim = np.logical_and(lhf2 >= -50, lhf2 <= 50)
    sns.distplot(lhf2[indexlrcasim], hist=False, color="forestgreen", kde_kws={"shade": True})
    ##---
    fmst_lrflux = obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1]
    indexlrfmst = np.logical_and(fmst_lrflux>=-50, fmst_lrflux<=50)
    sns.distplot(fmst_lrflux[indexlrfmst], hist=False, color="black")#, kde_kws={"shade": True})
    ##---
    lrflux = np.squeeze(obs['ice_station_fluxes']['lrflux'][obs['ice_station_fluxes']['lrflag'][:] == 1])
    indexlr = np.logical_and(lrflux>=-50, lrflux<=50)
    sns.distplot(lrflux[indexlr], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3})
    plt.title('latent_heat_flux [W/m2]')
    plt.xlim([-20,50])
    plt.ylim([0,yDmax])

    fileout = '../FIGS/comparisons/SHF_LHF_line+PDFS_oden_foremast+iceStationQC_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_paperRadiation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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
    print ('Plotting combined timeseries and PDFs of radiation terms:')
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
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

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

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12))

    ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot([240.0,240.0],[yB[0],yB[-1]],'--', color='grey')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('$SW_{net,surf}$ [W/m2]')
    # plt.legend()
    ax.set_xlim([doy[0],doy[-1]])
    plt.ylim([-3,120])

    ax  = fig.add_axes([0.07,0.4,0.53,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yC = [-90, 10]
    plt.plot([240.0,240.0],[yC[0],yC[-1]],'--', color='grey')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'forestgreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'darkorange')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'darkorange')
    plt.title('$LW_{net,surf}$ [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.ylim([-90,5])

    ax  = fig.add_axes([0.07,0.1,0.53,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yA = [-65, 85]
    plt.plot([240.0,240.0],[yA[0],yA[-1]],'--', color='grey')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'steelblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'forestgreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'darkorange', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'darkorange', label = label3)
    plt.title('CRF [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.legend(bbox_to_anchor=(-0.11, 0.67, 1., .102), loc=4, ncol=2)
    plt.ylim([-60,80])
    plt.xlabel('Day of year')    

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)

    #### only compare with observations where we have data:
    ####        all model data share a timestamp
    melt = np.where(np.logical_and(data1['time_hrly'] >= time_radice[0], data1['time_hrly'] < 240.0))
    freeze = np.where(data1['time_hrly'] >= 240.0)

    obsmelt = np.where(time_radice < 240.0)
    obsfreeze = np.where(time_radice >= 240.0)

    sw1 = data1['surface_net_SW_radiation'][data1['hrly_flag']]
    lw1 = data1['surface_net_LW_radiation'][data1['hrly_flag']]
    sw2 = data2['surface_net_SW_radiation'][data2['hrly_flag']]
    lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    sw3 = data3['sfc_net_sw'][data3['hrly_flag']]
    lw3 = data3['sfc_net_lw'][data3['hrly_flag']]

    ax  = fig.add_axes([0.64,0.7,0.15,0.22])   # left, bottom, width, height
    yEmax = 0.08
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(sw1[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(sw3[melt], hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(sw2[melt], hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netSW[obsmelt], hist=False, color="black")
    # plt.title('Melt')
    plt.annotate('Melt', xy=(87,0.07), xytext=(87,0.07), fontsize = 14)
    # plt.legend()
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('$SW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.64,0.4,0.15,0.22])   # left, bottom, width, height
    yFmax = 0.16
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    sns.distplot(lw1[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(lw3[melt], hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(lw2[melt], hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW[obsmelt], hist=False, color="black")
    # plt.title('Melt')
    plt.annotate('Melt', xy=(0,0.14), xytext=(0,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('$LW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.64,0.1,0.15,0.22])   # left, bottom, width, height
    yDmax = 0.08
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1[melt] + lw1[melt]
    sns.distplot(crf1, hist=False, color="steelblue", kde_kws={"shade": True})
    crf3 = sw3[melt] + lw3[melt]
    sns.distplot(crf3, hist=False, color="darkorange", kde_kws={"shade": True})
    crf2 = sw2[melt] + lw2[melt]
    sns.distplot(crf2, hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW[obsmelt] + netSW[obsmelt], hist=False, color="black")
    # plt.title('Melt')
    plt.annotate('Melt', xy=(55,0.07), xytext=(55,0.07), fontsize = 14)
    plt.xlabel('CRF [W/m2]')
    plt.xlim([-50,80])
    plt.ylim([0,yDmax])

    ax  = fig.add_axes([0.83,0.7,0.15,0.22])   # left, bottom, width, height
    yEmax = 0.08
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(sw1[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(sw3[freeze], hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(sw2[freeze], hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netSW[obsfreeze], hist=False, color="black")
    plt.annotate('Freeze', xy=(77,0.07), xytext=(77,0.07), fontsize = 14)
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('$SW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.83,0.4,0.15,0.22])   # left, bottom, width, height
    yFmax = 0.16
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    sns.distplot(lw1[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(lw3[freeze], hist=False, color="darkorange", kde_kws={"shade": True})
    sns.distplot(lw2[freeze], hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW[obsfreeze], hist=False, color="black")
    plt.annotate('Freeze', xy=(-8,0.14), xytext=(-8,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('$LW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.83,0.1,0.15,0.22])   # left, bottom, width, height
    yDmax = 0.08
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1[freeze] + lw1[freeze]
    sns.distplot(crf1, hist=False, color="steelblue", kde_kws={"shade": True})
    crf3 = sw3[freeze] + lw3[freeze]
    sns.distplot(crf3, hist=False, color="darkorange", kde_kws={"shade": True})
    crf2 = sw2[freeze] + lw2[freeze]
    sns.distplot(crf2, hist=False, color="forestgreen", kde_kws={"shade": True})
    sns.distplot(netLW[obsfreeze] + netSW[obsfreeze], hist=False, color="black")
    # plt.title('Freeze')
    plt.annotate('Freeze', xy=(45,0.07), xytext=(45,0.07), fontsize = 14)
    plt.xlim([-50,80])
    plt.ylim([0,yDmax])
    plt.xlabel('CRF [W/m2]')


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/netSW_netLW_CRF_line+PDFS-gt230DOY_oden_iceStation_metum_ifs_casim-100_splitSeason.svg'
    plt.savefig(fileout)
    plt.show()

def plot_Precipitation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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
    print ('Plotting combined timeseries and PDFs of precipitation terms:')
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
    plt.figure(figsize=(8,4.5))
    plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

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

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)



    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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

    data1['rainfall_flux'][data1['rainfall_flux'] < 0] = np.nan
    data2['rainfall_flux'][data2['rainfall_flux'] < 0] = np.nan
    # flx_ls_rain = np.nansum(data3['flx_ls_rain'],1)
    # flx_conv_rain = np.nansum(data3['flx_conv_rain'],1)
    # flx_conv_rain[flx_conv_rain < 0] = np.nan
    # flx_ls_snow = np.nansum(data3['flx_ls_snow'],1)

    flx_ls_rain = data3['flx_ls_rain'][:,0]         #### just take surface value
    flx_ls_snow = data3['flx_ls_snow'][:,0]             #### assumes all precip which forms at altitudes
                                                        #### above evaporates/sublimes before it reaches
                                                        #### the surface

    # flx_ls_rain = np.nanmean(data3['flx_ls_rain'],1)         #### take average rate
    # flx_ls_snow = np.nanmean(data3['flx_ls_snow'],1)

    # flx_ls_rain = np.nansum(data3['flx_ls_rain'],1)         #### take total
    # flx_ls_snow = np.nansum(data3['flx_ls_snow'],1)

    #### remove flagged values
    flx_ls_rain[flx_ls_rain < 0] = np.nan
    flx_ls_snow[flx_ls_snow < 0] = np.nan

    ### doy drift flag
    drift = np.squeeze(np.where(np.logical_and(obs['pws']['doy'] >= 226.0, obs['pws']['doy'] <= 258.0)))

    # print (drift.shape)

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    precip1 = data1['rainfall_flux'][data1['hrly_flag']].data*3600 + data1['snowfall_flux'][data1['hrly_flag']].data*3600
    precip2 = data2['rainfall_flux'][data2['hrly_flag']].data*3600 + data2['snowfall_flux'][data2['hrly_flag']].data*3600
    if ifs_flag: precip3 = flx_ls_rain[data3['hrly_flag']]*3600 + flx_ls_snow[data3['hrly_flag']]*3600

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    # fig = plt.figure(figsize=(10,5))

    # ax  = fig.add_axes([0.1,0.12,0.8,0.75])   # left, bottom, width, height
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['pws']['doy'][drift[0]],obs['pws']['prec_int'][drift[0]], 'k', label = 'Obs_PWS')
    plt.plot(data1['time_hrly'][::3], precip1[::3],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(data2['time_hrly'][::3], precip2[::3],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time_hrly'][::3], precip3[::3],
            'd', color = 'darkorange', markeredgecolor = 'saddlebrown', label = label3)
    plt.ylabel('Precipitation flux [mm/hr]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.ylim([0,2])
    plt.xlabel('Day of year')
    plt.legend()


    # ax  = fig.add_axes([0.75,0.15,0.2,0.7])   # left, bottom, width, height


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TotalPrecip_oden-pws_metum_ifs-z0_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

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

    #### make inversion tempvars to allow for easy subsampling
    inv1 = np.squeeze(data1['inversions']['invbase'][data1['hrly_flag'],0])
    inv2 = np.squeeze(data2['inversions']['invbase'][data2['hrly_flag'],0])
    inv3 = np.squeeze(data3['inversions']['invbase'][data3['hrly_flag'],0])

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
    fig = plt.figure(figsize=(13,7))

    #################################################################
    ## create figure and axes instances
    #################################################################
    ax  = fig.add_axes([0.08,0.56,0.45,0.36])   # left, bottom, width, height
    plt.plot(np.squeeze(obs['inversions']['doy_drift']),np.squeeze(obs['inversions']['sfmlheight'][drift]),
        color = 'k', label = 'Obs_Radiosondes')
    plt.plot(data1['time_hrly'][::6], bldepth1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(data2['time_hrly'][::6], bldepth2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(data3['time_hrly'][::6], bldepth3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    plt.legend(bbox_to_anchor=(0.0, 0.67, 1., .102), loc=4, ncol=2)
    plt.title('BL_depth / sfmlheight [m]')
    ax.set_xlim([doy[0],doy[-1]])
    # plt.xlabel('Day of year')
    plt.ylabel('Z [m]')

    ax  = fig.add_axes([0.08,0.1,0.45,0.36])   # left, bottom, width, height
    plt.plot(np.squeeze(obs['inversions']['doy_drift']), np.squeeze(obs['inversions']['invbase'][drift]),
        color = 'k', label = 'Obs: main inversion')
    plt.plot(data1['time_hrly'][::6], inv1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(data2['time_hrly'][::6], inv2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(data1['time_hrly'][::6], inv3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    # plt.legend()
    plt.title('Main inversion height [m]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')
    plt.ylabel('Z [m]')

    #### -----------------------------------------------------------------
    #### scatterplots
    #### -----------------------------------------------------------------
    ax  = fig.add_axes([0.6,0.64,0.15,0.22])   # left, bottom, width, height
    ymax1 = 750
    xmax1 = ymax1
    blmelt1 = bldepth1[melt]
    blmelt2 = bldepth2[melt]
    blmelt3 = bldepth3[melt]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{SML}$ [m]')
    plt.ylabel('Model$_{SML}$ [m]')
    plt.title('Melt')

    ax  = fig.add_axes([0.83,0.64,0.15,0.22])   # left, bottom, width, height
    ymax1 = 1500
    xmax1 = ymax1
    blfreeze1 = bldepth1[freeze]
    blfreeze2 = bldepth2[freeze]
    blfreeze3 = bldepth3[freeze]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{SML}$ [m]')
    plt.ylabel('Model$_{SML}$ [m]')
    plt.title('Freeze')

    ax  = fig.add_axes([0.6,0.18,0.15,0.22])   # left, bottom, width, height
    ymax1 = 2700
    xmax1 = ymax1
    invmelt1 = inv1[melt]
    invmelt2 = inv2[melt]
    invmelt3 = inv3[melt]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{inv}$ [m]')
    plt.ylabel('Model$_{inv}$ [m]')
    plt.title('Melt')

    ax  = fig.add_axes([0.83,0.18,0.15,0.22])   # left, bottom, width, height
    ymax1 = 2700
    xmax1 = ymax1
    invfreeze1 = inv1[freeze]
    invfreeze2 = inv2[freeze]
    invfreeze3 = inv3[freeze]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze1[::6],
        '^', color = 'steelblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze2[::6],
        'v', color = 'forestgreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze3[::6],
        'd', color = 'darkorange', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{inv}$ [m]')
    plt.ylabel('Model$_{inv}$ [m]')
    plt.title('Freeze')


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_oden_metum_ifs_casim-100.svg'
    # fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_metum.svg'
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
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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

def plot_RadiosondesTemperature(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'temp')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')
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
    fig = plt.figure(figsize=(19,10))

    Tmin = -45
    Tmax = 5
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['temperature'][:,drift[0]],
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, T[degC]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['temp_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', T[degC]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['temp_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', T[degC]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['temp_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.xlabel('Day of year')
    plt.title(label2 + ', T[degC]')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), T[degC]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['temp_hrly_UM'][::6])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), T[degC]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', T[degC]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', T[degC]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
        vmin = Tmin, vmax = Tmax)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), T[degC]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), T[K]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), T[K]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), T[K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TemperatureProfiles_REGRID_10km_sondes_metum_ifs_casim-100.png'
    plt.savefig(fileout, dpi = 300)
    plt.show()
    # plt.close()


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
    fig = plt.figure(figsize=(7,10))

    Tmin = -45
    Tmax = 5
    ymax = 8000

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.15,0.78,0.85,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
        vmin = Tmin, vmax = Tmax)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), T[degC]')

    ax  = fig.add_axes([0.15,0.54,0.85,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data3['temp_anomalies'] = dat3
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), T[K]')

    ax  = fig.add_axes([0.15,0.3,0.85,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data1['temp_anomalies'] = dat1
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), T[K]')

    ax  = fig.add_axes([0.15,0.06,0.85,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data2['temp_anomalies'] = dat2
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), T[K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TemperatureProfiles_REGRID_anomOnly_sondes_metum_ifs_casim-100.png'
    plt.savefig(fileout, dpi = 300)
    plt.show()
    # plt.close()

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
    fig = plt.figure(figsize=(11,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)

    plt.subplot(131)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M')
    plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100')
    plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS')
    plt.legend()
    plt.ylim([0,1e4])
    plt.xlim([-1.6,1.0])
    plt.ylabel('Z [m]')
    plt.xlabel('median T anomaly [K]')
    plt.grid('on')
    plt.title('Total drift')

    plt.subplot(132)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M median')
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100 median')
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS median')
    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-1.6,1.0])
    plt.xlabel('median T anomaly [K]')
    plt.title('Melt')

    plt.subplot(133)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M median')
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100 median')
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS median')
    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-1.6,1.0])
    plt.xlabel('median T anomaly [K]')
    plt.title('Freeze')

    fileout = '../FIGS/comparisons/TemperatureMedianProfiles_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesQ(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model thetaE profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

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
    fig = plt.figure(figsize=(19,10))

    qmax = 4.0
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['mr'][:,drift[0]],
        vmin = 0.0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, q [g/kg]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['q_6hrly'])*1e3,
        vmin = 0.0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', q [g/kg]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['q_6hrly'])*1e3,
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', q [g/kg]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['q_6hrly'])*1e3,
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ', q [g/kg]')
    plt.xlabel('Day of year')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:]),
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), q [g/kg]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['q_hrly_UM'][::6])*1e3,
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), q [g/kg]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3,
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', q [g/kg]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3,
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', q [g/kg]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:]),
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID),  q [g/kg]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID),  q [g/kg]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), q [g/kg]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), q [g/kg]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/QProfiles_REGRID_10km_sondes_metum_ifs_casim-100.png'
    plt.savefig(fileout, dpi = 300)
    plt.show()

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    #### anomaly plot on it's own

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
    fig = plt.figure(figsize=(7,10))

    ymax = 8000
    qmax = 4.0

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.15,0.78,0.85,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_driftSondes_UM']),
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID),  q [g/kg]')

    ax  = fig.add_axes([0.15,0.54,0.85,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data3['q_anomalies'] = dat3
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID),  q [g/kg]')

    ax  = fig.add_axes([0.15,0.3,0.85,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data1['q_anomalies'] = dat1
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), q [g/kg]')

    ax  = fig.add_axes([0.15,0.06,0.85,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data2['q_anomalies'] = dat2
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), q [g/kg]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/QProfiles_REGRID_anomOnly_sondes_metum_ifs_casim-100.png'
    plt.savefig(fileout, dpi = 300)
    plt.show()
    # plt.close()

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
    fig = plt.figure(figsize=(11,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)

    plt.subplot(131)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M')
    plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100')
    plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS')
    plt.legend()
    plt.ylim([0,1e4])
    plt.xlim([-0.05,0.45])
    plt.ylabel('Z [m]')
    plt.xlabel('median q anomaly [g/kg]')
    plt.grid('on')
    plt.title('Total drift')

    plt.subplot(132)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M median')
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100 median')
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS median')
    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-0.05,0.45])
    plt.xlabel('median q anomaly [g/kg]')
    plt.title('Melt')

    plt.subplot(133)
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'steelblue', label = 'UM_RA2M median')
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'forestgreen', label = 'UM_CASIM-100 median')
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1),data1['universal_height'],
        color = 'darkorange', label = 'ECMWF_IFS median')
    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-0.05,0.45])
    plt.xlabel('median q anomaly [g/kg]')
    plt.title('Freeze')

    fileout = '../FIGS/comparisons/QMedianProfiles_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model thetaE profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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
    data1['theta'], data1['thetaE'] = calcThetaE(data1['temperature'], data1['pressure'], data1['q'], data1['time'], data1['height'])
    data2['theta'], data2['thetaE'] = calcThetaE(data2['temperature'], data2['pressure'], data2['q'], data2['time'], data2['height'])
    data3['theta'], data3['thetaE'] = calcThetaE(data3['temperature'], data3['pressure'], data3['q'], data3['time'], np.squeeze(data3['height'][0,:]))

    obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
        np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3,
        obs['sondes']['doy'], obs['sondes']['gpsaltitude'])

    obs['sondes']['thetaE'] = np.transpose(obs['sondes']['thetaE'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'thetaE')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

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
    fig = plt.figure(figsize=(19,10))

    Tmin = -5
    Tmax = 45
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['epottemp'][:,drift[0]],
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta_{E}$ [degC]')
    plt.xlabel('Day of year')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['thetaE_hrly_UM'][::6])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['thetaE_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['thetaE_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta_{E}$ [degC]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['thetaE_hrly_UM'][::6]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), $\Theta_{E}$ [K]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['thetaE_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), $\Theta_{E}$ [K]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['thetaE_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), $\Theta_{E}$ [K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/ThetaEProfiles_REGRID_10km_sondes-calculated_metum_ifs_casim-100.png'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesTheta(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model theta profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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
    data1['theta'], data1['thetaE'] = calcThetaE(data1['temperature'], data1['pressure'], data1['q'], data1['time'], data1['height'])
    data2['theta'], data2['thetaE'] = calcThetaE(data2['temperature'], data2['pressure'], data2['q'], data2['time'], data2['height'])
    data3['theta'], data3['thetaE'] = calcThetaE(data3['temperature'], data3['pressure'], data3['q'], data3['time'], np.squeeze(data3['height'][0,:]))

    obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
        np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3,
        obs['sondes']['doy'], obs['sondes']['gpsaltitude'])

    obs['sondes']['theta'] = np.transpose(obs['sondes']['theta'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'theta')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

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
    fig = plt.figure(figsize=(19,10))

    Tmin = -5
    Tmax = 45
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['pottemp'][:,drift[0]],
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta$ [degC]')
    plt.xlabel('Day of year')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['theta_hrly_UM'][::6])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['theta_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['theta_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta$ [degC]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['theta_hrly_UM'][::6]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), $\Theta$ [K]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['theta_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), $\Theta$ [K]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['theta_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), $\Theta$ [K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/ThetaProfiles_REGRID_10km_sondes-calculated_metum_ifs_casim-100.png'
    plt.savefig(fileout)
    plt.show()

def reGrid_Sondes(data1, data2, data3, obs, doy, var):

    from scipy.interpolate import interp1d

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

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

    ### stop double counting of 0000 and 2400 from model data
    temp = np.zeros([len(data1['time'])])
    for i in range(0, len(temp)-1):
        if data1['time'][i] == data1['time'][i+1]:
            continue
        else:
            temp[i] = data1['time'][i]
    ii = np.where(temp != 0.0)      ### picks out where data are non-zero

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
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems
            ####        further on
    data3['height_hrly'] = np.squeeze(data3['height'][ii,:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    print ('')
    print ('Defining IFS temperature profile as a function:')
    print ('using ifs.height[i,:] to define temperature profiles...')
    data3[var + '_hrly_UM'] = np.zeros([np.size(data3['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(data3['time_hrly'],0)):
        # print (iTim)
        iIFSind = np.where(data3['height_hrly'][iTim,:] <= 11000)
        if np.all(data3['height_hrly'][iTim,:] == 0.0):
            data3[var + '_hrly_UM'][iTim,:] = np.nan
        else:
            fnct_IFS = interp1d(np.squeeze(data3['height_hrly'][iTim,iIFSind]), np.squeeze(data3[var + '_hrly'][iTim,iIFSind]))
            data3[var + '_hrly_UM'][iTim,:] = fnct_IFS(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('IFS(UM Grid) function worked!')
    print (var + ' IFS data now on UM vertical grid')
    print ('*****')
    ### assign for easier indexing later
    data3[var + '_6hrly_UM'] = data3[var + '_hrly_UM'][::6,:]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6,:]
    data3['height_6hrly'] = data3['height_hrly'][::6,:]  ### need to explicitly save since height coord changes at each timedump

    #### INTERPOLATION TESTING:
    # print (data3['temp_hrly_UM'].shape)
    # print (data3['time_hrly'][::6].shape)
    # print (data1['temp_hrly'][:,iUM[0][3:]].shape)
    # print (data1['time_hrly'][::6].shape)
    # for i in range(0, np.size(data3['temp_6hrly_UM'],0)):
    #     fig = plt.figure()
    #     plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], label = 'interpd')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), label = 'height indexed')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height'][0,iIFS]), label = 'height0')
    #     plt.title('IFS test ' + str(data3['time_6hrly'][i]))
    #     plt.legend()
    #     plt.savefig('../FIGS/regrid/IFS_test_doy' + str(data3['time_6hrly'][i]) + '.png')
    #     if i == 0:
    #         plt.show()
    #     else:
    #         plt.close()

    print ('')
    print ('Defining Sonde temperature profile as a function for the UM:')
    obs['sondes'][var + '_allSondes_UM'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(obs['sondes']['doy'],0)):
        # print 'iTim = ', str(iTim)
        fnct_Obs = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
        obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_Obs(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('Sonde(UM Grid) function worked!')
    print ('All ' + var + ' sonde data now on UM vertical grid.')
    print ('*****')
    #
    # print ('')
    # print ('Defining Sonde temperature profile as a function for the IFS:')
    # obs['sondes'][var + '_allSondes_IFS'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][0,iIFS])])
    # for iTim in range(0,np.size(obs['sondes']['doy'],0)):
    #     # print 'iTim = ', str(iTim)
    #     iIFS = np.where(data3['height'][iTim,:] <= 11000)
    #     fnct_ObsIFS = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
    #     obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_ObsIFS(data3['height'][iTim,iIFS])
    # print ('...')
    # print ('Sonde(IFS Grid) function worked!')
    # print ('All ' + var + ' sonde data now on IFS_DATA vertical grid.')
    # print ('*****')

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['sondes']['doy'] >= 225.9, obs['sondes']['doy'] <= 258.0))

    ### save in dict for ease
    obs['sondes']['doy_drift'] = obs['sondes']['doy'][drift]
    obs['sondes']['drift'] = drift
    obs['sondes'][var + '_driftSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][drift[0],:]

    #### INTERPOLATION TESTING - IFS + SONDE + UM_RA2M:
    # print (obs['sondes']['doy_drift'].shape)
    # print (obs['sondes']['temp_allSondes_UM'][drift[0],:].shape)
    # if var == 'temp':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,drift[0][i]]) + 273.15,np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['temp_driftSondes_UM'][i,:] + 273.15,data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'darkorange', label = 'ifs-Zindexed')
    #         plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = 'darkorange', label = 'ifs-interpd')
    #         plt.plot(data1['temp_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2m')
    #         plt.plot(data2['temp_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'forestgreen', label = 'um_casim-100')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid/REGRID_Ttest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'q':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['mr'][iObs,drift[0][i]]), np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes'][var + '_driftSondes_UM'][i,:], data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3[var + '_6hrly'][i,iIFS])*1e3,np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'darkorange', label = 'ifs-Zindexed')
    #         plt.plot(data3[var + '_6hrly_UM'][i,:]*1e3,data1['height'][iUM[0][3:]], color = 'darkorange', label = 'ifs-interpd')
    #         plt.plot(data1[var + '_6hrly'][i,iUM[0][3:]]*1e3, data1['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2m')
    #         plt.plot(data2[var + '_6hrly'][i,iUM[0][3:]]*1e3, data2['height'][iUM[0][3:]], color = 'forestgreen', label = 'um_casim-100')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid/REGRID_Qtest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()

    #### ---------------------------------------------------------------
    #### make some dictionary assignments for use later
    #### ---------------------------------------------------------------
    data1['universal_height'] = data1['height'][iUM[0][3:]]
    data1['universal_height_UMindex'] = iUM[0][3:]

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data2',data2)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])
    outfiles = write_reGrid(data1, data2, data3, obs, var)

    return data1, data2, data3, obs, drift

def write_reGrid(data1, data2, data3, obs, var):

    #################################################################
    ## Write regridded data to netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing out data:')
    print ('')

    outfiles = ['REGRID-' + var + '-SONDES.nc','REGRID-' + var + '-UM_RA2M.nc','REGRID-' + var + '-UM_CASIM-100.nc','REGRID-' + var + '-ECMWF_IFS.nc']

    ###################################
    ## Open File
    ###################################
    nc0 = Dataset(outfiles[0], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc0.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc0.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc0.createDimension('time', np.size(obs['sondes']['doy_drift']))
    height = nc0.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    tempvar = np.zeros(np.size(obs['sondes']['doy_drift']))
    tempvar[:] = obs['sondes']['doy_drift'][:]
    times = nc0.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = tempvar[:]

    #### height
    height = nc0.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat0 = nc0.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat0.scale_factor = float(1)
    dat0.add_offset = float(0)
    if var == 'temp':
        dat0.units = 'K'
        dat0.long_name = 'temperature'
        dat0[:,:] = obs['sondes'][var + '_driftSondes_UM'][:,:] + 273.15
    elif var == 'q':
        dat0.units = 'g/kg'
        dat0.long_name = 'water vapour mixing ratio'
        dat0[:,:] = obs['sondes'][var + '_driftSondes_UM'][:,:]

    nc0.title = 'Radiosonde interpolated ' + var + ' data for the AO2018 drift period.'
    nc0.description = var + ' data up to 10 km, interpolated on to the UM vertical grid.'
    nc0.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc0.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc0.comment = 'Revision no. 0: Preliminary data.'
    nc0.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc0.close()


    ###################################
    ## Open File
    ###################################
    nc1 = Dataset(outfiles[1], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc1.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc1.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc1.createDimension('time', np.size(data1['time_6hrly']))
    height = nc1.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc1.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data1['time_6hrly'][:]

    #### height
    height = nc1.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat1 = nc1.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat1.scale_factor = float(1)
    dat1.add_offset = float(0)
    if var == 'temp':
        dat1.units = 'K'
        dat1.long_name = 'temperature'
        dat1[:,:] = data1[var + '_6hrly'][:,data1['universal_height_UMindex']]
    elif var == 'q':
        dat1.units = 'g/kg'
        dat1.long_name = 'water vapour mixing ratio'
        dat1[:,:] = data1[var + '_6hrly'][:,data1['universal_height_UMindex']]*1e3

    nc1.title = 'UM_RA2M ' + var + ' data for the AO2018 drift period.'
    nc1.description = var + ' data up to 10 km, referencing UM vertical grid.'
    nc1.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc1.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc1.comment = 'Revision no. 0: Preliminary data.'
    nc1.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc1.close()

    ###################################
    ## Open File
    ###################################
    nc2 = Dataset(outfiles[2], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc2.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc2.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc2.createDimension('time', np.size(data2['time_6hrly']))
    height = nc2.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc2.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data2['time_6hrly'][:]

    #### height
    height = nc2.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat2 = nc2.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat2.scale_factor = float(1)
    dat2.add_offset = float(0)
    if var == 'temp':
        dat2.units = 'K'
        dat2.long_name = 'temperature'
        dat2[:,:] = data2[var + '_6hrly'][:,data1['universal_height_UMindex']]
    elif var == 'q':
        dat2.units = 'g/kg'
        dat2.long_name = 'water vapour mixing ratio'
        dat2[:,:] = data2[var + '_6hrly'][:,data1['universal_height_UMindex']]*1e3

    nc2.title = 'UM_CASIM-100 ' + var + ' data for the AO2018 drift period.'
    nc2.description = var + ' data up to 10 km, referencing UM vertical grid.'
    nc2.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc2.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc2.comment = 'Revision no. 0: Preliminary data.'
    nc2.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc2.close()


    ###################################
    ## Open File
    ###################################
    nc3 = Dataset(outfiles[3], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc3.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc3.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc3.createDimension('time', np.size(data3['time_6hrly']))
    height = nc3.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc3.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data3['time_6hrly'][:]

    #### height
    height = nc3.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat3 = nc3.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat3.scale_factor = float(1)
    dat3.add_offset = float(0)
    if var == 'temp':
        dat3.units = 'K'
        dat3.long_name = 'temperature'
        dat3[:,:] = data3[var + '_6hrly_UM'][:,:]
    elif var == 'q':
        dat3.units = 'g/kg'
        dat3.long_name = 'water vapour mixing ratio'
        dat3[:,:] = data3[var + '_6hrly_UM'][:,:]*1e3

    nc3.title = 'ECMWF_IFS interpolated ' + var + ' data for the AO2018 drift period.'
    nc3.description = var + ' data up to 10 km, interpolated on to the UM vertical grid.'
    nc3.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc3.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc3.comment = 'Revision no. 0: Preliminary data.'
    nc3.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc3.close()

    return outfiles

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
        platformflag = 'jasmin'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        um_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        misc_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        obs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/'
    if platform == 'LAPTOP':
        platformflag = 'laptop'
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        misc_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    if platform == 'LAPTOP':
        out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/'
        out_dir2 = '5_u-bl661_RA1M_CASIM/OUT_R0/'
        # out_dir3 = 'MET_DATA/'
        out_dir4 = 'OUT_25H/'
    elif platform == 'JASMIN':
        out_dir1 = 'UM_RA2M/'
        out_dir2 = 'UM_CASIM-100/'
        # out_dir3 = 'MET_DATA/'
        out_dir4 = 'ECMWF_IFS/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R0/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param
    ### 12_u-br210_RA1M_CASIM/OUT_24h/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507
    ### 13_u-br409_RA1M_CASIM/OUT_24h/           # 100/cc accum mode aerosol; ARG + Cooper; passive aerosol processing

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

    ### print ('Load ice station radiation data from Jutta...')
    ### obs['ice_station_radiation'] = readMatlabStruct(obs_root_dir + 'ice_station/mast_radiation_30min_v2.3.mat')

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

    doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)        ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)
    # doy = np.arange(243,259)        ## set DOY for CASIM-AeroProf (31st Aug to 14th Sep)
    # doy = np.arange(244,254)        ## set DOY for CASIM-100_AP (1st Sep to 9th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):
        filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
        filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
        if np.logical_or(out_dir4 == 'OUT_25H/', out_dir4 == 'ECMWF_IFS/'):
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
        if ifs_flag: var_list3 = ['height','flx_height','temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
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

            print ('Starting on t=0 UM data:')
            for j in range(0,len(var_list1)):
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) >= 1:
                    data1[var_list1[j]] = nc1.variables[var_list1[j]][:]
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            print ('Starting on t=0 CASIM data:')
            for j in range(0,len(var_list2)):
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
                    data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
            nc2.close()
            ## ------------------
            #### um3
            ## ------------------
            print ('Starting on t=0 IFS data:')
            for j in range(0,len(var_list3)):
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) >= 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
            nc3.close()
            print ('')
        else:
            if month_flag == -1:
                time_um1 = np.append(time_um1, doy[i] + (nc1.variables['forecast_time'][:]/24.0))
                time_um2 = np.append(time_um2, doy[i] + (nc2.variables['forecast_time'][:]/24.0))
                if ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['time'][:]/24.0))
                if not ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['forecast_time'][:]/24.0))
            ## ------------------
            #### UM
            ## ------------------
            print ('Appending UM data:')
            for j in range(0,len(var_list1)):
                # print (var_list1[j])
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:])
                elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:],0)
            # np.save('working_data1',data1)
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            print ('Appending CASIM data:')
            for j in range(0,len(var_list2)):
                # print (var_list2[j])
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
            print ('Appending IFS data:')
            for j in range(0,len(var_list3)):
                # print (var_list3[j])
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:])
                elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:],0)
            nc3.close()
            print ('')

    #################################################################
    ## save time to dictionary now we're not looping over all diags anymore
    #################################################################
    data1['time'] = time_um1
    data2['time'] = time_um2
    data3['time'] = time_um3

    ### stop double counting of 0000 and 2400 from model data
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

    #### add override for data2 to allow 24h data to be used for testing purposes
    if out_dir2[-4:] == '24h/':
        data2['time_hrly'] = data2['time']
        data2['hrly_flag'] = np.arange(len(data2['time_hrly']))

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
    if out_dir1[:10] == '13_u-br409': label1 = 'UM_CASIM-100_AP'
    if out_dir1[:10] == '12_u-br210': label1 = 'UM_CASIM-AeroProf'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T'
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'
    if out_dir1 == 'UM_RA2M/': label1 = 'UM_RA2M'

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
    if out_dir2 == 'UM_CASIM-100/': label2 = 'UM_CASIM-100'

    label3 = 'undefined_label'
    if out_dir4 == 'OUT_25H/': label3 = 'ECMWF_IFS'
    if out_dir4[:10] == '13_u-br409': label3 = 'UM_CASIM-100_AP'
    if out_dir4[:10] == '12_u-br210': label3 = 'UM_CASIM-AeroProf'
    if out_dir4[:10] == '11_u-bq798': label3 = 'UM_CASIM-100_Meyers'
    if out_dir4[:10] == '10_u-bq791': label3 = 'UM_CASIM-100_Fletcher'
    if out_dir4[:9] == '8_u-bp738': label3 = 'UM_ERAI-GLM'
    if out_dir4[:9] == '7_u-bn068': label3 = 'UM_RA2T'
    if out_dir4[:9] == '6_u-bm410': label3 = 'UM_CASIM-200'
    if out_dir4[:9] == '5_u-bl661': label3 = 'UM_CASIM-100'
    if out_dir4[:9] == '4_u-bg610': label3 = 'UM_RA2M'
    if out_dir4 == 'ECMWF_IFS/': label3 = 'ECMWF_IFS'


    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    # np.save('working_dataObs', obs['sondes'])

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
    figure = plot_paperRadiation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_Precipitation(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_BLType(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesTemperature(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesQ(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesTheta(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)
    # figure = plot_line_subSect(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir4, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_dataObs', obs['sondes'])

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
