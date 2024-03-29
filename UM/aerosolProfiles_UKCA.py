###
###
### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
###
###

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
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
from time_functions import calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcAirDensity


def serial_date_to_string(srl_no):

    '''
    Function to convert date (wrt 1981-10-01)
    '''

    new_date = datetime.datetime(1981,10,1,0,0) + datetime.timedelta(srl_no)
    return new_date.strftime("%Y-%m-%d")


def readfile(filename):

    import pandas as pd

    # print('******'
    print('')
    print('Reading .txt file with pandas')
    print('')

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

    print('******')
    print('')
    # print('Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print('Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print('Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4]))
    print('')

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

    print('******')
    print('')
    # print('Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print('Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print('In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print('CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print('')
    print('Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print('Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print('Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print('Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print('')
    print('Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))

    return inIce_index

def trackShip(data):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==25,data.values[:,1]==8),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print('******')
    print('')
    # print('Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print('Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print('Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print('Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print('End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print('trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print('')

    return trackShip_index

def plot_cartmap(ship_data, nc, doy):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.path as mpath

    ###################################
    ## PLOT MAP
    ###################################

    print('******')
    print('')
    print('Plotting cartopy map:')
    print('')

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
    plt.figure(figsize=(6,8))#, dpi=300)
    ax = plt.axes(projection=ccrs.Orthographic(30, 70))    # NP Stereo
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([0, 60, 86.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE
    # ax.set_extent([-180, 180, 70, 90], crs=ccrs.PlateCarree())    ### V LARGE
    # ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())    ### POSTER
    # ax.set_global()

    ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

    #################################################################
    ## add geographic features/guides for reference
    #################################################################
    ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=0)
    ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=0, edgecolor='black')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.set_global()
    ax.gridlines()

    #################################################################
    ## choose date + time to plot
    ################################################################
    plt.pcolormesh(nc.variables['number_concentration_of_soluble_accumulation_mode_aerosol'][0,0,:,:],
             transform = ccrs.PlateCarree())

    #################################################################
    ## plot UM data
    ################################################################
    # iplt.pcolormesh(cube[diag][290:370,150:230])


    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data)

    ### Plot full track as line plot
    plt.plot(ship_data.values[:,6], ship_data.values[:,7], '--',
             color = 'pink', linewidth = 2,
             transform = ccrs.PlateCarree(), label = 'Whole',
             )
    plt.plot(ship_data.values[inIce_index,6], ship_data.values[inIce_index,7],
             color = 'palevioletred', linewidth = 3,
             transform = ccrs.PlateCarree(), label = 'In Ice',
             )
    plt.plot(ship_data.values[inIce_index[0],6], ship_data.values[inIce_index[0],7],
             'k^', markerfacecolor = 'palevioletred', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )
    plt.plot(ship_data.values[inIce_index[-1],6], ship_data.values[inIce_index[-1],7],
             'kv', markerfacecolor = 'palevioletred', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )
    plt.plot(ship_data.values[drift_index,6], ship_data.values[drift_index,7],
             color = 'red', linewidth = 4,
             transform = ccrs.PlateCarree(), label = 'Drift',
             )


    plt.legend()

    print('******')
    print('')
    print('Finished plotting cartopy map! :)')
    print('')

    # plt.savefig('FIGS/HighArctic_vPOSTER.svg', dpi=100)
    plt.show()

def plot_aeroProfiles(nc2, nc3, doy):

    ###################################
    ## PLOT AEROSOL PROFILES
    ###################################

    print('******')
    print('')
    print('Plotting aerosol profiles for whole drift period:')
    print('')

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
    fig = plt.figure(figsize=(8,6))
    plt.subplots_adjust(top = 0.93, bottom = 0.1, right = 0.82, left = 0.08,
            hspace = 0.3, wspace = 0.1)

    plt.subplot(211)
    ax = plt.gca()
    img = plt.pcolormesh(nc3.variables['day_of_year'][:],nc3.variables['level_height'][:],
        np.transpose(np.nanmean(np.nanmean(nc3.variables['number_concentration_of_soluble_coarse_mode_aerosol'][:,:,-2:,:],3),2)),
        vmin = 0, vmax = 0.3)
    plt.ylim([0, 1e4])
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,3.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    ax.set_yticks(axminor, minor = True)
    cbaxes = fig.add_axes([0.85, 0.6, 0.015, 0.3])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    plt.ylabel('N$_{aer, sol, coarse}$ [cm$^{-3}$]', rotation = 270, labelpad = 25)

    plt.subplot(212)
    ax = plt.gca()
    img = plt.pcolormesh(nc2.variables['day_of_year'][:],nc2.variables['level_height'][:],
        np.transpose(np.nanmean(np.nanmean(nc2.variables['number_concentration_of_soluble_accumulation_mode_aerosol'][:,:,-2:,:],3),2)),
        vmin = 0, vmax = 200)
    plt.ylim([0, 1e4])
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.xlabel('Date')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    ax.set_yticks(axminor, minor = True)
    cbaxes = fig.add_axes([0.85, 0.12, 0.015, 0.3])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    plt.ylabel('N$_{aer, sol, accum}$ [cm$^{-3}$]', rotation = 270, labelpad = 25)

    fileout = '../FIGS/UKCA/UKCA_aeroProfiles_Lat-2_LonAll_226-257DOY_Dates.png'
    # plt.savefig(fileout, dpi = 300)
    # plt.show()
    plt.close()

    Z = nc2.variables['level_height'][:]
    naer = np.nanmean(np.nanmean(nc2.variables['number_concentration_of_soluble_accumulation_mode_aerosol'][:,:,-2:,:],3),2)

    Zindex = np.where(Z < 500.)

    print ('UKCA_Naer,sol,accum = ')
    print (np.nanmean(np.nanmean(naer[:,Zindex[0]],1),0))
    print (np.nanstd(np.nanmean(naer[:,Zindex[0]],1),0))

def interpolate_aeroProfiles(nc1, nc2, nc3, doy, ukca_index):

    from scipy.interpolate import interp1d

    ###################################
    ## interpolate UKCA aerosol profiles from Ruth on to my UM height grid and MONC Case study grid
    ###################################

    ### design monc vertical grid:
    monc_Z = np.arange(0., 1000., 50.)
    monc_Z = np.append(monc_Z, np.arange(1000., 2501., 100.))

    print('******')
    print('')
    print('Calculating aerosol profiles on UM nZ=70 grid:')
    print('')

    naer_accum = np.nanmean(np.nanmean(nc2.variables['number_concentration_of_soluble_accumulation_mode_aerosol'][:,:,-2:,:],3),2)
    naer_coarse = np.nanmean(np.nanmean(nc3.variables['number_concentration_of_soluble_coarse_mode_aerosol'][:,:,-2:,:],3),2)

    print('naer_accum shape = ', naer_accum.shape)
    print('naer_coarse shape = ', naer_coarse.shape)
    print('')

    um_height = nc1.variables['height'][:]
    ukca_height = nc2.variables['level_height'][:]

    print('um_height shape = ', um_height.shape)
    print('ukca_height shape = ', ukca_height.shape)
    print('')

    fnct_accum = interp1d(ukca_height, naer_accum[ukca_index,:])
    print('Accumulation mode interpolation function succeeded!')
    print('')
    print('Next: test function on um_height')
    newAccum = fnct_accum(um_height[3:].data)        ### z=3 == 22m, lower altitude bins below 1st UKCA bin
    moncAccum = fnct_accum(monc_Z[1:])
    print('')
    print('Function worked! :)')
    print('')

    fnct_coarse = interp1d(ukca_height, naer_coarse[ukca_index,:])
    print('Coarse mode interpolation function succeeded!')
    print('')
    print('Next: test function on um_height')
    newCoarse = fnct_coarse(um_height[3:].data)        ### z=3 == 22m, lower altitude bins below 1st UKCA bin
    moncCoarse = fnct_coarse(monc_Z[1:])
    print('')
    print('Function worked! :)')
    print('')

    #### ----------------------------------------------------------
    #### ukca timestamps are midpoints of bin (+/- 12h)
    ####        seem to be at midnight when looking at DOY, but forecast_period
    ####        and forecast_period_bnds suggest it's a daily average and the time
    ####        step should be at midday
    #### ----------------------------------------------------------

    # ##################################################
    # ##################################################
    # #### 	SET UP FIGURE PROPERTIES
    # ##################################################
    # ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.figure(figsize=(12,7))
    # plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.13,
    #         hspace = 0.4, wspace = 0.2)
    #
    # plt.subplot(121)
    # plt.plot(naer_accum[ukca_index,:], ukca_height, label = 'UKCA')
    # plt.plot(newAccum, um_height[3:], '--', label = 'UM Interpolated')
    # plt.ylabel('Z [m]')
    # plt.xlabel('N$_{aer, accum}$ [cm$^{-3}$]')
    # plt.ylim([0, 40000])
    # plt.legend()
    #
    # plt.subplot(122)
    # plt.plot(naer_coarse[ukca_index,:], ukca_height, label = 'UKCA')
    # plt.plot(newCoarse, um_height[3:], '--', label = 'UM Interpolated')
    # # plt.ylabel('Z [m]')
    # plt.xlabel('N$_{aer, coarse}$ [cm$^{-3}$]')
    # plt.ylim([0, 40000])
    # # plt.legend()
    #
    # fileout = '../FIGS/UKCA/UKCA_aeroProfiles_UM-interp1d_example.svg'
    # # plt.savefig(fileout)
    # plt.show()

    #### ----------------------------------------------------------
    #### next, print out newAccum and newCoarse on UM binning for
    ####    use in um/rose-app.conf namelist
    #### ----------------------------------------------------------

    print('newAccum.shape = ', newAccum.shape)
    print('newCoarse.shape = ', newCoarse.shape)

    #### build new arrays with extra 2 bins at the surface
    numAccum = np.zeros(len(um_height)-1)   ### 70 height bins
    numCoarse = np.zeros(len(um_height)-1)   ### 70 height bins

    numAccum[0:2] = newAccum[0]*1e6
    numAccum[2:] = newAccum[:]*1e6
    print('numAccum = ', numAccum)

    numCoarse[0:2] = newCoarse[0]*1e6
    numCoarse[2:] = newCoarse[:]*1e6
    print('numCoarse = ', numCoarse)

    print('')
    print('numAccum.shape = ', numAccum.shape)
    print('numCoarse.shape = ', numCoarse.shape)

    #### fill in lowest MONC levels with data at [1]
    moncNumAccum = np.zeros(np.size(monc_Z))
    moncNumAccum[0] = moncAccum[0]*1e6
    moncNumAccum[1:] = moncAccum[:]*1e6

    moncNumCoarse = np.zeros(np.size(monc_Z))
    moncNumCoarse[0] = moncCoarse[0]*1e6
    moncNumCoarse[1:] = moncCoarse[:]*1e6

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
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.96, left = 0.13,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(121)
    plt.plot(naer_accum[ukca_index,:]*1e6, ukca_height, label = 'UKCA')
    plt.plot(numAccum, um_height[1:], '--', label = 'UM Interpolated')
    plt.plot(moncNumAccum, monc_Z, 'k.', label = 'MONC Interpolated')
    plt.plot([1e8, 1e8], [monc_Z[0], monc_Z[-1]], '--', color='grey', label = '100 cm$^{-3}$ approximation')
    plt.ylabel('Z [m]')
    plt.xlabel('N$_{aer, accum}$ [m$^{-3}$]')
    plt.ylim([0, 2500])
    plt.xlim([0, 1.1e8])
    # plt.grid('on')
    plt.legend()

    plt.subplot(122)
    plt.plot(naer_coarse[ukca_index,:]*1e6, ukca_height, label = 'UKCA')
    plt.plot(numCoarse, um_height[1:], '--', label = 'UM Interpolated')
    plt.plot(moncNumCoarse, monc_Z, 'k.', label = 'MONC Interpolated')
    # plt.ylabel('Z [m]')
    plt.xlabel('N$_{aer, coarse}$ [m$^{-3}$]')
    plt.ylim([0, 2500])
    plt.xlim([0, 8e4])
    # plt.grid('on')
    # plt.legend()

    fileout = '../FIGS/UKCA/UKCA_UM_MONC_numSolAccum_numSolCoarse_20180913.svg'
    # plt.savefig(fileout)
    plt.show()

    return numAccum, numCoarse, moncNumAccum, moncNumCoarse

def scaleMass(numAccum, numCoarse):

    #### -------------------------------------------------------------
    #### SCALE AEROSOL MASS (accumulation mode: 1.5*1e-9 for every 1.00*1e8 aerosol particles)
    #### -------------------------------------------------------------

    print('')
    print('****')
    print('Scaling mass by default values used in CASIM suite')
    print('')
    print('accumulation mode: 1.5*1e-9 for every 1.00*1e8 aerosol particles')
    print('')

    ### calculate scaling factor
    factor = 1.5e-9 / 1e8

    massAccum = factor * (numAccum*1e2)
    print('massAccum = ', massAccum)
    print('')

    ### add filler for now
    massCoarse = 0

    return massAccum, massCoarse

def estimateMass(N, rho_air, flag):

    #### -------------------------------------------------------------
    #### SCALE AEROSOL MASS (accumulation mode: 1.5*1e-9 for every 1.00*1e8 aerosol particles)
    #### -------------------------------------------------------------

    print('')
    print('****')
    print('Estimate mass by mean modal radius and assuming spherical particles:')
    print('')

    #### Accumulation mode: 0.1um < d_p < 1um

    ### make dummy variables
    # M = 1.0
    if flag == 1:
        sigma = 1.5         #### == fixed_aerosol_sigma (mphys_constants.F90)
        rho = 1777.0        #### == fixed_aerosol_density (mphys_constants.F90); kg/m3
        Rm = 0.5*1.0e-6     #### == fixed_aerosol_rm (mphys_constants.F90); 500nm
    elif flag == 2:
        sigma = 1.5         #### == fixed_aerosol_sigma (mphys_constants.F90)
        rho = 2000.0        #### == fixed_aerosol_density (mphys_constants.F90); kg/m3
        Rm = 5*1.0e-6       #### == fixed_aerosol_rm (mphys_constants.F90); 5 um
    else:
        print('****Mode option not valid!****')

    print('Calculating aerosol mass mixing ratio assuming: ')
    print('rho_aer = ', rho, ' g/cm3')
    print('Rm = ', Rm*1e6, ' um')
    print('...')

    ### calculation for mean radius given mass and number:
    # MNtoRm = ( 3.0*M*np.exp(-4.5*np.log(sigma)**2) /
    #     (4.0*N*np.pi*rho) )**(1.0/3.0)
                ### just copied from casim/lognormal_funcs.F90

    mass = ( (4.0/3.0)*np.pi*Rm**3 ) * (N*rho) / (np.exp(-4.5*np.log(sigma)**2))
            ### gives mass concentration in kg/m3

    #### need mass concentration in kg/kg for casim input
    M = mass / rho_air

    print('mass = ', M)
    print('')

    return M

def main():

    START_TIME = time.time()
    print('******')
    print('')
    print('Start: ' + time.strftime("%c"))
    print('')

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
        um_root_dir = '/home/gillian/MOCCHA/MOCCHA_GIT/UM/DATA/'
        obs_root_dir = '/home/gillian/MOCCHA/MOCCHA_GIT/ODEN/DATA/'
        ship_filename = '~/MOCCHA/MOCCHA_GIT/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir1 = '4_u-bg610_RA2M_CON/OUT_R1/'
    out_dir2 = 'UKCA/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R0/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param

    print('******')
    print('')
    print('Identifying .nc file: ')
    print('')

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print('******')
    print('')
    print('Load in ship track file:')
    print('')
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    print('Loading observations:')
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

    print('Load temporary ice station data from Jutta...')
    obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_wTemp1p5m.nc','r')

    print('Load ice station data from Jutta...')
    obs['ice_station'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')
            #### mast_radiation_30min_v2.3.mat
            #### flux30_trhwxrel.mat

    print('Load radiosonde data from Jutta...')
    obs['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print('Load foremast data from John...')
    obs['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print('Load 7th deck weather station data from John...')
    obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print('...')

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print('******')
    print('')
    print('Begin cube read in at ' + time.strftime("%c"))
    print(' ')

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

    doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for RA2T  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)          ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Flag for individual file or monthly:
    combine = 1
    ## Choose month:
    names = moccha_names
    missing_files = moccha_missing_files
    month_flag = -1

    #### -------------------------------------------------------------
    #### LOAD UKCA DATA
    #### -------------------------------------------------------------
    filename_um2 = um_root_dir + out_dir2 + 'number_concentration_of_soluble_accumulation_mode_aerosol.nc'
    filename_um3 = um_root_dir + out_dir2 + 'number_concentration_of_soluble_coarse_mode_aerosol.nc'
    #### LOAD DATASET
    print('Loading UKCA aerosol data (accumulation and coarse mode):')
    nc2 = Dataset(filename_um2,'r')
    print('...')
    nc3 = Dataset(filename_um3,'r')
    print('...')

    #### -------------------------------------------------------------
    #### -------------------------------------------------------------
    #### CHOOSE DATE
    #### -------------------------------------------------------------
    #### -------------------------------------------------------------

    date = '20180913'
    doyIndex = calcTime_Date2DOY(date)
    ukca_index = np.where(nc2.variables['day_of_year'][:] == doyIndex)
    # um_index = np.where(data1['time'][:] == doyIndex)

    print('****')
    print('UKCA time = ', doyIndex, ' at index = ', np.squeeze(ukca_index))
    print('Proof: nc2.variables[''day_of_year''][ukca_index] = ', nc2.variables['day_of_year'][ukca_index])
    print('****')
    print('')

    #### -------------------------------------------------------------
    #### LOAD UM DATA
    #### -------------------------------------------------------------
    filename_um1 = um_root_dir + out_dir1 + date + '_oden_metum.nc'
    print(filename_um1)
    print('')

    #### LOAD DATASET
    print('Loading UM diagnostics for reference:')
    nc1 = Dataset(filename_um1,'r')
    print('...')
    # ---------------------

    #### -------------------------------------------------------------
    #### PLOT MAP
    #### -------------------------------------------------------------
    # figure = plot_cartmap(ship_data, nc2, doy)

    #### -------------------------------------------------------------
    #### PLOT AEROSOL PROFILES
    #### -------------------------------------------------------------
    figure = plot_aeroProfiles(nc2, nc3, doy)

    #### -------------------------------------------------------------
    #### CREATE N_AER PROFILES (IN /M)
    #### -------------------------------------------------------------
    # numAccum, numCoarse, moncNumAccum, moncNumCoarse = interpolate_aeroProfiles(nc1, nc2, nc3, doy, np.squeeze(ukca_index))

    # data = {}
    # data['moncNumAccum'] = moncNumAccum
    # data['moncNumCoarse'] = moncNumCoarse
    # np.save('MONC_UKCAInputs-20180913', data)

    ## ### -------------------------------------------------------------
    ## ### SCALE AEROSOL MASS
    ## ###        Accum: 1.5*1e-9 for every 1.00*1e8 aerosol particles
    ## ### -------------------------------------------------------------
    ## massAccum, massCoarse = scaleMass(numAccum, numCoarse)

    #### -------------------------------------------------------------
    #### CALCULATE DAY MEAN AIR DENSITY
    ####        for use in mass conversion calculation
    #### -------------------------------------------------------------
    # # print(nc1.variables['pressure'][0,0])
    # rho_air = calcAirDensity(np.nanmean(nc1.variables['temperature'][:,1:],0),np.nanmean(nc1.variables['pressure'][:,1:]/1e2,0))
    # plt.plot(rho_air[:],nc1.variables['height'][1:]);plt.show()

    #### -------------------------------------------------------------
    #### ESTIMATE AEROSOL MASS
    ####        assume spherical particles
    #### -------------------------------------------------------------

    # print('****')
    # print('Estimate accumulation mode mass:')
    # print('')
    # modeFlag = 1
    # massAccum = estimateMass(numAccum, rho_air, modeFlag)
    # plt.plot(massAccum,nc1.variables['height'][1:]); plt.title('massAccum'); plt.show()
    #
    # print('****')
    # print('Estimate coarse mode mass:')
    # print('')
    # modeFlag = 2
    # massCoarse = estimateMass(numCoarse, rho_air, modeFlag)
    # plt.plot(massCoarse,nc1.variables['height'][1:]);  plt.title('massCoarse');plt.show()

    #### -------------------------------------------------------------
    #### FORMAT OUTPUT FOR UM/ROSE-APP.CONF
    #### -------------------------------------------------------------
    # outfile = 'AUX_DATA/CASIM-AeroProf_input_' + date + '.txt'
    #
    # ### 1) print("testing with i = %i, num is %f" % (i, num[i]))
    # ### 2) for i in range(0,4): print("testing with i = %i, num is %f" % (i, num[i]))
    # ### 3) for i in range(0,len(num)): print("%f, " % (num[i]), end=" ") ### prints dynamically to the same line
    #
    # print('Formatting output for use in run_casim_tracers UM namelist...')
    # print('')
    #
    # ### open text file
    # file = open(outfile,'w+')
    #
    # file.write("accum_sol_num_var=%f, " % (numAccum[0]))#, end=" ")
    # for i in range(1,len(numAccum)-1): file.write("%f, " % (numAccum[i]))#, end=" ")
    # file.write("%f \n" % (numAccum[-1]))
    # file.write('\n')
    # file.write("accum_sol_mass_var=%e, " % (massAccum[0]))#, end=" ")
    # for i in range(1,len(massAccum)-1): file.write("%e, " % (massAccum[i]))#, end=" ")
    # file.write("%e \n" % (massAccum[-1]))
    # file.write('\n')
    # file.write("coarse_sol_num_var=%f, " % (numCoarse[0]))#, end=" ")
    # for i in range(1,len(numCoarse)-1): file.write("%f, " % (numCoarse[i]))#, end=" ")
    # file.write("%f \n" % (numCoarse[-1]))
    # file.write('\n')
    # file.write("coarse_sol_mass_var=%e, " % (massCoarse[0]))#, end=" ")
    # for i in range(1,len(massCoarse)-1): file.write("%e, " % (massCoarse[i]))#, end=" ")
    # file.write("%e " % (massCoarse[-1]))
    # file.write('')
    #
    # ### close text file
    # file.close()
    #
    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print('******')
    print('')
    print('End: ' + time.strftime("%c"))
    print('')

if __name__ == '__main__':

    main()
