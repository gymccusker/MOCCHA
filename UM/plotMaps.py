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

def plot_driftMap(data1, data3, obs, doy):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.path as mpath
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting cartopy map:')
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
    ax.set_extent([0, 60, 86.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE
    # ax.set_extent([-180, 180, 70, 90], crs=ccrs.PlateCarree())    ### V LARGE
    # ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())    ### POSTER
    # ax.set_global()

    ### DON'T USE: PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

    #################################################################
    ## add geographic features/guides for reference
    #################################################################
    ax.add_feature(cartopy.feature.OCEAN, color='white', zorder=0)
    ax.add_feature(cartopy.feature.LAND, color='lightgrey', zorder=0, edgecolor='black')
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.set_global()
    ax.gridlines()

    #################################################################
    ## plot UM data
    ################################################################
    iplt.pcolormesh(cube[diag][290:370,150:230])

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[diag][hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[diag][hour,386:495,211:305])          ### misc
    # qplt.outline(cube[diag][290:370,150:230])



    # gridship = gridShipTrack(cube[diag], xoffset, yoffset)

            #### MID POINT: (433, 258)

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data, date)

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

    print ('******')
    print ('')
    print ('Finished plotting cartopy map! :)')
    print ('')

    # plt.savefig('FIGS/HighArctic_vPOSTER.svg', dpi=100)
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
        platformflag = 'jasmin'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        um_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        ifs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        # misc_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        obs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/'
    if platform == 'LAPTOP':
        platformflag = 'laptop'
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        ifs_root_dir = '/home/gillian/MOCCHA/ECMWF/'
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
    ### 12_u-br210_RA1M_CASIM/OUT_R0/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

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
    ###             just need one date - it's the grid that matters!
    ### -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Identifying iris cubes: ')
    print ('')

    name = '20180901_oden_'
    filename_um = um_root_dir + out_dir1 + name + 'metum.nc'
    filename_ifs = ifs_root_dir + out_dir2 + name + 'ecmwf.nc'
    print (filename_um)
    print (filename_ifs)
    print ('')

    #### LOAD CUBE
    print ('Loading first run diagnostics:')
    cube_um = iris.load(filename_um)
    print ('...')
    print ('Loading second run diagnostics:')
    cube_ifs = iris.load(filename_ifs)
    # -------------------------------------------------------------
    print ('...')

    print (cube_um)
    print ('')
    print (cube_ifs)
    print ('')

    # -------------------------------------------------------------
    # Plot map
    # -------------------------------------------------------------
    # figure = plot_driftMap(data1, data3, obs, doy)


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
