###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import datetime
from iris.time import PartialDateTime

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')

# STASH_CODE = 'm01s04i118'


def readfile(filename):

    import pandas as pd

    # print '******'
    print ('')
    print ('Reading .txt file with pandas')
    print ('')

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data, values

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
    # Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print ('******')
    print ('')
    print ('CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print ('')
    print ('Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print ('Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print ('Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print ('Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))
    print ('')

    return inIce_index

def readGriddedTrack(grid_filename):

    import pandas as pd

    print ('******')
    print ('')
    print ('Reading ' + grid_filename + ' file with pandas')
    print ('')

    data = pd.read_csv(grid_filename, sep = " ")
    values = data.values

    tim = values[:,1]
    ilon = values[:,2]
    ilat = values[:,3]

    return tim, ilat, ilon

def readGlobal(cube, ship_data, date):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy

    print ('******')
    print ('')
    print ('Defining longitude and latitude boundaries:')
    print ('')

    if np.ndim(cube[date][0].data) == 3:      ### 2D data + time
        lats = cube[date][0].dim_coords[1].points
        lons = cube[date][0].dim_coords[2].points
    if np.ndim(cube[date][0].data) == 4:      ### 3D data + time
        lats = cube[date][0].dim_coords[2].points
        lons = cube[date][0].dim_coords[3].points

    ###---------------------------------------------------------------------------------
    ### find northern and southern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    nb_lats = lats + ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    sb_lats = lats - ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    print ('sb_lats.shape = ' + str(sb_lats.shape))

    ###---------------------------------------------------------------------------------
    ### find western and eastern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    wb_lons = lons - ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    eb_lons = lons + ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    #################################################################
    ## find date of interest
    #################################################################
    month = int(date[5])
    day = int(date[6:8])
    endpt = int(day)+2

    trackShip_start = np.where(np.logical_and(np.logical_and(ship_data.values[:,2]==day,ship_data.values[:,1]==month),ship_data.values[:,3]>=12))
    trackShip_end = np.where(np.logical_and(ship_data.values[:,2]==endpt, ship_data.values[:,1]==month))

    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][0])
    day_ind = trackShip_index

    print(day_ind)

    # day_ind = np.array([])
    # day_ind = np.where(np.logical_and(ship_data.values[:,2] == float(date[6:8]),ship_data.values[:,1] == float(date[4:6])))
    print ('Daily ship track for ' + date + ': ' + str(len(day_ind)) + ' pts ')

    #################################################################
    ## save daily lat/lons as temp vars
    #################################################################
    ship_lats = ship_data.values[day_ind,7]
    ship_lons = ship_data.values[day_ind,6]

    print ('ship_lats.shape = ' + str(ship_lats.shape))

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    ### compare hourly lat-lon with GLM grid
    data = {}
    data['ship_lons'] = np.zeros(36)
    data['ship_lats'] = np.zeros(36)
    data['ship_i'] = np.zeros(36); data['ship_i'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_j'] = np.zeros(36); data['ship_j'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_hour'] = np.zeros(36)
    hours = np.arange(0,36)
    jflag = np.zeros(36)        ### flag for if grid boundary is crossed

    for h in hours:
        print ('')
        print ('hour = ' + str(h))
        for j in range(0,len(sb_lats)):     ### for all latitude points
            if np.logical_and(ship_lats[h] >= sb_lats[j], ship_lats[h] < nb_lats[j]):     ### find where ship lat sits on glm lat grid
                for i in range(0,len(wb_lons)):     ### for all longitude points
                    if np.logical_and(ship_lons[h] >= wb_lons[i], ship_lons[h] < eb_lons[i]):     ### find where ship lon sits on glm lon grid
                        print ('lats and lons match at i = ' + str(i) + ', j = ' + str(j))
                        jflag[h] = jflag[h] + 1
                        data['ship_lons'][h] = lons[i]
                        data['ship_hour'][h] = hours[h]
                        data['ship_lats'][h] = lats[j]
                        data['ship_j'][h] = j         # define grid point indices for use later
                        data['ship_i'][h] = i         # define grid point indices for use later

    # print data['ship_lats']
    # print data['ship_j']
    # print data['ship_lons']
    # print data['ship_i']

    # np.save('working_glm_grid', data)

    ### need to constract an hour, lon index, lat index list like used for the lam
    # for h in hours:
    #     if data['ship_i'][h+1] == data['ship_i'][h]:

    ### arguments to be passed back to pullTrack_CloudNet
    tim = hours
    ilon = data['ship_i']
    ilat = data['ship_j']

    print (tim)
    print (ilon)
    print (ilat)

    return tim, ilat, ilon

def trackShip(data, date):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################

    print ('Extracting: ' + date)

    month = int(date[5])
    day = int(date[6:8])
    endpt = int(day)+2

    print ('month = ' + str(month))
    print ('start point = ' + str(day))
    print ('end point = ' + str(endpt))

    print (data)

    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==day,data.values[:,1]==month),data.values[:,3]>=12))
    trackShip_end = np.where(np.logical_and(data.values[:,2]==endpt ,data.values[:,1]==month))

    print (trackShip_start[0])
    print (trackShip_end[0])

    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][0])

    print ('******')
    print ('')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print ('Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print ('trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print ('')

    return trackShip_index

def plot_cartmap(ship_data, cube, date_dir, model):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting cartopy map:')
    print ('')

    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    diag = 0
    dd = date_dir[0]
    if len(cube[dd][diag].dim_coords[-1].points) == 25:
    # if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif len(cube[dd][diag].dim_coords[-1].points) == 56:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif len(cube[dd][diag].dim_coords[-1].points) == 94:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -211
        yoffset = -385
    elif len(cube[dd][diag].dim_coords[-1].points) == 81:          ### 14th and 24th August
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -209
        yoffset = -399
    elif len(cube[dd][diag].dim_coords[-1].points) == 380:         ### needs checked
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -60
        yoffset = -110
    else:
    # elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

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


    #################################################################
    ## create figure and axes instances
    #################################################################
    i = 0
    for date in date_dir:
        if date[:4] == '2018':
            i = i + 1
            plt.figure(figsize=(5,4))
            plt.subplots_adjust(top = 0.85, bottom = 0.15, right = 0.97, left = 0.1,
                    hspace = 0.3, wspace = 0.3)
            # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))


            ### set size
            # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())      ### full track
            # ax.set_extent([20, 50, 88.35, 89.95], crs=ccrs.PlateCarree())     ### SWATH
            ax.set_extent([20, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM

            ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

            #################################################################
            ## add geographic features/guides for reference
            #################################################################
            ax.add_feature(cartopy.feature.OCEAN, zorder=0)
            ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
            # ax.set_global()
            ax.gridlines()

            #################################################################
            ## plot UM data
            #################################################################
            # iplt.pcolormesh(cube[date][0][0,:,:], vmin = 0.8, vmax = 1.0, cmap=mpl_cm.Blues_r)
            # plt.title(cube[date][0].standard_name)
            # plt.colorbar()

            #################################################################
            ## plot UM nest
            #################################################################
            ### draw outline of grid
            if model == 'lam': qplt.outline(cube[date][0][0,:,:])

            #################################################################
            ## plot ship track
            #################################################################
            ### DEFINE DRIFT + IN_ICE PERIODS
            drift_index = iceDrift(ship_data)
            inIce_index = inIce(ship_data)
            trackShip_index = trackShip(ship_data, date)

            ### Plot tracks as line plot
            plt.plot(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7],
                     color = 'darkorange', linewidth = 3,
                     transform = ccrs.PlateCarree(), label = date,
                     )
            plt.plot(ship_data.values[trackShip_index[0],6], ship_data.values[trackShip_index[0],7],
                     'k^', markerfacecolor = 'darkorange', linewidth = 3,
                     transform = ccrs.PlateCarree(),
                     )
            plt.plot(ship_data.values[trackShip_index[-1],6], ship_data.values[trackShip_index[-1],7],
                     'kv', markerfacecolor = 'darkorange', linewidth = 3,
                     transform = ccrs.PlateCarree(),
                     )

            #################################################################
            ## read in and plot gridded ship track
            #################################################################
            date_extension =[0, 1]

            for t in date_extension:
                if model == 'lam':
                    grid_dirname = 'AUX_DATA/'
                    if int(date[6:8]) <= 8: grid_filename = grid_dirname + date[:6] + '0' + str(int(date[6:8])+t) + '_ShipTrack_GRIDDED.csv'
                    if int(date[6:8]) >= 9: grid_filename = grid_dirname + date[:6] + str(int(date[6:8])+t) + '_ShipTrack_GRIDDED.csv'

                    tim, ilat, ilon = readGriddedTrack(grid_filename)

                    ### Plot tracks as line plot
                    if t == 0:
                        loop_index = np.where(tim>=12.0)        ### only from 1200Z on day 0
                    else:
                        loop_index = np.where(tim>=0.0)         ### from 0000Z on day 1
                    times = tim[loop_index]
                    lons = ilon[loop_index]
                    lats = ilat[loop_index]
                    cc = ['black', 'blue']

                    ###
                    ### plot for sanity check
                    ###
                    for i in range(0, len(times)):
                        iplt.scatter(cube[date][0].dim_coords[2][int(lons[i] + xoffset)], cube[date][0].dim_coords[1][int(lats[i] + yoffset)],color=cc[t])
                        print (times[i])

                    ###
                    ### prepare output arrays for writing
                    ###
                    if t == 0:
                        time_forecast = np.copy(times)
                        time_forecast = time_forecast - 12.0
                        lons_forecast = np.copy(lons)
                        lats_forecast = np.copy(lats)
                    else:
                        times_extended = times + 12.0
                        time_forecast = np.append(time_forecast, times_extended)
                        lons_forecast = np.append(lons_forecast, lons)
                        lats_forecast = np.append(lats_forecast, lats)

                    out = writeout36HGrid(time_forecast, lats_forecast, lons_forecast, date, model)

                    # for i in range(0, len(time_forecast)):
                    #     iplt.scatter(cube[date][0].dim_coords[2][int(lons_forecast[i] + xoffset)], cube[date][0].dim_coords[1][int(lats_forecast[i] + yoffset)], color='red')

                elif model == 'glm':

                    qplt.outline(cube[date][0][0,:,::10])

                    tim, ilat, ilon = readGlobal(cube, ship_data, date)
                    out = writeout36HGrid(tim, ilat, ilon, date, model)
                    ###
                    ### plot for sanity check
                    ###
                    for i in range(0, len(tim)):
                        iplt.scatter(cube[date][0].dim_coords[2][int(ilon[i] + xoffset)], cube[date][0].dim_coords[1][int(ilat[i] + yoffset)],color='k')
                        # print (tim[i])

            plt.legend()

            plt.savefig('../../FIGS/Grid_' + model + '_ZoomedTrack_' + date + '.png')
            plt.close()

    print ('******')
    print ('')
    print ('Finished plotting cartopy map! :)')
    print ('')

def writeout36HGrid(tim, lat, lon, date, model):

    import pandas as pd

    # ******
    # write to csv file
    # ******

    print ('******')
    print ('Writing ' + date + ' grid to file:')
    print ('')
    dat = np.zeros([len(tim), 3])
    dat[:,0] = tim
    dat[:,1] = lon
    dat[:,2] = lat
    df = pd.DataFrame(dat)
    if model == 'lam': filename = 'AUX_DATA/' + date + '-36HForecast_ShipTrack_GRIDDED.csv'
    if model == 'glm': filename = 'AUX_DATA/' + date + '-36HForecast-GLM_ShipTrack_GRIDDED.csv'
    df.to_csv(filename,  sep = " ")
    print ('... finished!')
    print ('')
    print ('******')

def plot_basemap(ship_data, cube):

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plot basemap:')
    print ('')

    ##################################################
    ##################################################
    #### 	BASEMAP
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

    ## create figure and axes instances
    fig = plt.figure(figsize=(6,8))

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 2000000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=86,lat_0=86,lon_0=10)
    m.drawcoastlines()
    m.bluemarble()

    # define parallels/meridians
    m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,0,0,0],linewidth=0.8,fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],linewidth=0.8,fontsize=10)
    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='lightgrey')
    # m.fillcontinents(color='white')

    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)

    ### MAP ONTO PROJECTION
    x, y = m(ship_data.values[:,6], ship_data.values[:,7])
    x_inIcePeriod, y_inIcePeriod = m(ship_data.values[inIce_index,6],ship_data.values[inIce_index,7])
    x_driftPeriod, y_driftPeriod = m(ship_data.values[drift_index,6],ship_data.values[drift_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, color = 'yellow', linewidth = 2, label = 'Whole')
    plt.plot(x_inIcePeriod, y_inIcePeriod, color = 'darkorange', linewidth = 3, label = 'In Ice')
    plt.plot(x_driftPeriod, y_driftPeriod, color = 'red', linewidth = 4, label = 'Drift')

    ###########################################
    ### PLOT NEST + SWATH FOR INCREASED FREQ DIAGS VIS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629
        #
        # R.P.: original 1.5km nest -> (0, 86.625) @ 500x500

    ### SWATH
    lats = np.arange(80.9998,89.9998,0.09)
    lons = np.arange(3.0,76.0,0.73)
    print ('******')
    print ('')
    print ('lat/lon vertices of proposed swath: ', lats[0], lats[-1], lons[0], lons[-1])
    print ('')
    x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)

    ### NEST (input)
    grx = float(500)
    gry = float(500)
    centlon = float(0.0)
    centlat = float(86.625)
    latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    print ('******')
    print ('')
    print ('lat/lon vertices of nest (input): ', latn[0], latn[-1], lonn[0], lonn[-1])
    print ('')
    x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    ### NEST (output)
    lono, lato = unrotateGrid(cube)
    print ('******')
    print ('')
    print ('lat/lon vertices of nest (output): ', lato[0], lato[-1], lono[0], lono[-1])
    print ('')
    x1o, x2o, x3o, x4o, y1o, y2o, y3o, y4o = gridSetup(lono, lato, m)

    ### NEST (required input to cover swath)
    # grx = float(5600)
    # gry = float(700)
    # centlon = float(39.5)
    # centlat = float(85.275)
    # latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    # lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    # print '******'
    # print ''
    # print 'lat/lon vertices of nest (input): ', latn[0], latn[-1], lonn[0], lonn[-1]
    # print ''
    # x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    # draw swath
    pols =  Polygon([(x1s,y1s),(x2s,y2s),(x3s,y3s),(x4s,y4s)],\
                  facecolor='none',linestyle='--',edgecolor='g',linewidth=2,label='Swath')
    plt.gca().add_patch(pols)

    # draw nest (input)
    poln =  Polygon([(x1n,y1n),(x2n,y2n),(x3n,y3n),(x4n,y4n)],\
                  facecolor='none',linestyle='-',edgecolor='w',linewidth=2,label='Nest')
    plt.gca().add_patch(poln)

    # draw nest (output)
    polo =  Polygon([(x1o,y1o),(x2o,y2o),(x3o,y3o),(x4o,y4o)],\
                  facecolor='none',linestyle='-',edgecolor='y',linewidth=2,label='Nest (output)')
    plt.gca().add_patch(polo)

    # draw nest (output)
    # polot =  Polygon([(x1ot,y1ot),(x2ot,y2ot),(x3ot,y3ot),(x4ot,y4ot)],\
    #               facecolor='none',linestyle='--',edgecolor='m',linewidth=2,label='Test nest')
    # plt.gca().add_patch(polot)

    ### ADD LEGEND
    plt.legend()

    plt.show()

def gridSetup(lont, latt, m):

    lon, lat = np.meshgrid(lont, latt)

    nx = np.size(lon,0)
    ny = np.size(lat,1)
    x1, y1 = m(lon[nx-1,0],lat[nx-1,0])
    x2, y2 = m(lon[0,0],lat[0,0])
    x3, y3 = m(lon[0,ny-1],lat[0,ny-1])
    x4, y4 = m(lon[nx-1,ny-1],lat[nx-1,ny-1])

    return x1, x2, x3, x4, y1, y2, y3, y4

def unrotateGrid(data):
    ##
    ##      ROTATE GRID BACK TO STANDARD
    ##

    import iris.analysis.cartography as ircrt

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 37.5000
    pole_lon = 177.5000

    rot_lat = data.coord('grid_latitude').points
    rot_lon = data.coord('grid_longitude').points

    rot_pole = data.coord('grid_latitude').coord_system.as_cartopy_crs()

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, frst_lon, frst_lat)

    # Print to check conversion
    print ('******')
    print ('Test of unrotated coordinate grid: ')
    print ('Rotated lon coord = ', rot_lon[0])
    print ('Rotated lat coord = ', rot_lat[0])
    print ('Lon = ', lon[0])
    print ('Lat = ', lat[0])
    print (' ')

    # ******
    # lat/lon vertices of nest (output):  85.960219407715 80.41973098346767 49.567255645848604 -27.55740381723681
    # ******

    return lon, lat

def rotateGrid(data):
    ##
    ##      RELATE GRID TO ROTATED POLE
    ##

    import iris.analysis.cartography as ircrt

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 37.5000
    pole_lon = 177.5000

    lat = np.arange((centlat-(gry*float(0.5)*dlat)),(centlat+(gry*float(0.5)*dlat)),dlat)
    lon = np.arange((centlon-(grx*float(0.5)*dlon)),(centlon+(grx*float(0.5)*dlon)),dlon)

    rot_lon, rot_lat = ircrt.rotate_pole(lon, lat, frst_lon, frst_lat)

    # Print to check conversion
    print ('******')
    print ('Test of rotated coordinate grid: ')
    print ('Lon = ', lon[0])
    print ('Lat = ', lat[0])
    print ('Rotated lon coord = ', rot_lon[0])
    print ('Rotated lat coord = ', rot_lat[0])
    print (' ')

    # ******
    # lat/lon vertices of nest (output):  85.960219407715 80.41973098346767 49.567255645848604 -27.55740381723681
    # ******

    return lon, lat

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    return GlobalStashList

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def assignTimecoord(cube):
    '''
    assign Time as dimension coordinate if not present in cube
    '''
    if not cube.coords('time', dim_coords=True):
        time = cube.coord('time')
        time.bounds = None
        iris.util.promote_aux_coord_to_dim_coord(cube, time)

    return cube

def fixTimecoord(local_cube_list):

    import cf_units
    import iris.coords as icoords
    from datetime import datetime

    period_1 = 0

    time_ref = datetime(1970,1,1,0,0,0,0)
    time_start = datetime(2018,8,11,12,0,0,0)

    time_unit = cf_units.Unit(
        'hours since 1970-01-01', calendar=cf_units.CALENDAR_GREGORIAN)

    time.gmtime((cube.coord('time')[0].points)*3600)    ## TIME IN UTC

    # build ref Time coord
    LBYR = str(local_cube_list[period_1].attributes['LBYR'])
    if local_cube_list[period_1].attributes['LBMON'] >= 10:
        LBMON = str(local_cube_list[period_1].attributes['LBMON'])
    else:
        LBMON = '0' + \
            str(local_cube_list[period_1].attributes['LBMON'])
    if local_cube_list[period_1].attributes['LBDAT'] >= 10:
        LBDAT = str(local_cube_list[period_1].attributes['LBDAT'])
    else:
        LBDAT = '0' + \
            str(local_cube_list[period_1].attributes['LBDAT'])
    if local_cube_list[period_1].attributes['LBHR'] >= 10:
        LBHR = str(local_cube_list[period_1].attributes['LBHR'])
    else:
        LBHR = '0' + \
            str(local_cube_list[period_1].attributes['LBHR'])
    if local_cube_list[period_1].attributes['LBMIN'] > - 10:
        LBMIN = str(local_cube_list[period_1].attributes['LBMIN'])
    else:
        LBMIN = '0' + \
            str(local_cube_list[period_1].attributes['LBMIN'])

    scanTime = LBYR + '/' + LBMON + '/' + \
        LBDAT + ' ' + LBHR + ':' + LBMIN

    refTime = datetime.datetime.strptime(
        scanTime, "%Y/%m/%d %H:%M")

    refTimeCoord = icoords.AuxCoord(
        time_unit.date2num(refTime),
        standard_name='forecast_reference_time',
        units=time_unit)

    # Data Time coord
    LBYR = str(local_cube_list[period_2].attributes['LBYRD'])
    if local_cube_list[period_2].attributes['LBMOND'] >= 10:
        LBMON = str(local_cube_list[period_2].attributes['LBMOND'])
    else:
        LBMON = '0' + \
            str(local_cube_list[period_2].attributes['LBMOND'])
    if local_cube_list[period_2].attributes['LBDATD'] >= 10:
        LBDAT = str(local_cube_list[period_2].attributes['LBDATD'])
    else:
        LBDAT = '0' + \
            str(local_cube_list[period_2].attributes['LBDATD'])
    if local_cube_list[period_2].attributes['LBHRD'] >= 10:
        LBHR = str(local_cube_list[period_2].attributes['LBHRD'])
    else:
        LBHR = '0' + \
            str(local_cube_list[period_2].attributes['LBHRD'])
    if local_cube_list[period_2].attributes['LBMIND'] > - 10:
        LBMIN = str(local_cube_list[period_2].attributes['LBMIND'])
    else:
        LBMIN = '0' + \
            str(local_cube_list[period_2].attributes['LBMIND'])

    scanTime = LBYR + '/' + LBMON + '/' + \
        LBDAT + ' ' + LBHR + ':' + LBMIN

    dataTime = datetime.datetime.strptime(
        scanTime, "%Y/%m/%d %H:%M")

    dataTimeCoord = icoords.AuxCoord(time_unit.date2num(dataTime),
                                     standard_name='time',
                                     units=time_unit)

def testInput(cube):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Investigating input file:')
    print ('')

    print ('Cube attributes:')
    print (cube.attributes)
    print ('')
    print ('Cube aux coords:')
    print (cube.aux_coords)
    print ('')
    print ('Cube dim coords:')
    print (cube.dim_coords)
    print ('')
    print (np.size(cube.shape))

def loadUMStartDump(filename):

    '''
    Loads in data from UM start dump used to drive global model
    Only loads in specific data relevant for calculating cloud + thermodynamic biases
    '''

    varlist = ['air_potential_temperature', 'sea_ice_thickness', 'mass_fraction_of_cloud_liquid_water_in_air',
                'mass_fraction_of_cloud_ice_in_air', 'specific_humidity']

    cube = iris.load(filename, varlist)

    return cube

def combinePP(root_dir, out_dir, date_dir):

    # -------------------------------------------------------------
    # Define output stream filenames to look at:
    #           start at 012 if 3h dumps (a, b)
    #           start at 011 if 1h dumps (c--e)
    # -------------------------------------------------------------
    print ('Checking if .pp files need to be combined:')
    print ('...')

    dummy = []

    names = ['_pa000','_pb000','_pc000','_pd000','_pe000']

    print ('out_dir is: ' + out_dir)
    print ('...')
    print ('')

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Define time and Stash constraints:
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Combine output pp streams
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    for date in date_dir:
        if date[0:4] == '2018':
        # if date[0:8] == '20180903':
            ## -------------------------------------------------------------------------
            ## Set fixed variable constraint (i.e. which variable to load in based on stash code)
            ## -------------------------------------------------------------------------
            # print '******'
            # print ''
            # print 'Use fixed constraint at ' + time.strftime("%c")
            # print ' '
            # var_con = iris.AttributeConstraint(STASH=STASH_CODE)

            ## -------------------------------------------------------------------------
            ## Time: constrain to take data every hour
            ## -------------------------------------------------------------------------
            # time_con = iris.Constraint(forecast_period=lambda xcell: any(np.isclose(xcell.point % 1, [0, 1./60.])))

            print ('******')
            print ('')
            print ('Naming output files: ')
            print ('')
            pp2_filename = str(int(date[:-6])+1) + '_oden_metum.pp'
            nc_filename = str(int(date[:-6])+1) + '_oden_metum.nc'
            print ('Output will be: ' + nc_filename)

            # print '******'
            print ('')
            print ('Identifying .pp files to read in: ')
            print ('')
            for stream in names:
                ### -------------------------------------------------------------------------
                ### define output filenames
                ### -------------------------------------------------------------------------
                # 20180827T1200Z_glm_pc035.pp
                filename1 = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream + '.pp'
                filename2 = root_dir + out_dir + date + '/' + date + '_glm' + stream + '.pp'
                filenames = [filename1, filename2]
                for filename in filenames:
                    print ('Checking: ' + filename)
                    if os.path.exists(filename):
                        pp_filename = filename[:-3] + '_r0.pp'

                        print ('---')
                        print ('Start files exist, continuing:')
                        print ('')

                        # print (filename[-12:-9])

                        ### define range to loop over
                        if filename[-12:-6] == 'glm_pa':
                            looping = range(0,6)
                            print ('Loop over 6 hours')
                        elif np.logical_or(stream[1:3] == 'pa', stream[1:3] == 'pb'):
                            looping = range(0,12)
                            print ('Loop over 3 hours')
                        else:
                            looping = range(0,36)
                            print ('Loop every hour')

                        for i in looping:
                            if np.size(looping) == 12:
                                res = i*3     # how many hourly dumps in file
                            elif filename[-12:-6] == 'glm_pa':
                                res = i*6
                            else:
                                res = i
                            str_i = "%03d" % res # file number
                            # fileout = root_dir + out_dir + date + stream[:-3] + str_i
                            if filename == filename1: fileout = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream[:-3] + str_i + '.pp'
                            if filename == filename2: fileout = root_dir + out_dir + date + '/' + date + '_glm' + stream[:-3] + str_i + '.pp'
                            # fileout = root_dir + out_dir + date + '/umnsaa_pa' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pb' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pc' + str_i
                            # fileout = root_dir + out_dir + date + '/umnsaa_pd' + str_i
                            print (fileout)
                            # # -------------------------------------------------------------
                            # # Load cubes
                            # # -------------------------------------------------------------
                            print ('******')
                            print ('')
                            print ('Begin ' + str_i + ' cube read in at ' + time.strftime("%c"))
                            print (' ')

                            #### LOAD CUBE
                            if 'var_con' in locals():
                                cube = iris.load(fileout, var_con)

                                # -------------------------------------------------------------
                                # Write out data
                                # -------------------------------------------------------------
                                print ('******')
                                print ('')
                                print ('Outputting fixed constraint ' + str_i + ' data:')
                                print ('')
                                iris.save(cube, pp_filename, append=True)

                            elif 'global_con' in locals():
                                cube = iris.load(fileout, global_con, callback)

                                # -------------------------------------------------------------
                                # Write out data
                                # -------------------------------------------------------------
                                print ('******')
                                print ('')
                                print ('Outputting global constraint ' + str_i + ' data at ' + time.strftime("%c"))
                                print (cube)
                                print ('')
                                iris.save(cube, pp_filename, append=True)
                                #### remove file to keep directory tidy
                                print ('Directory clean up: removing ' + fileout)
                                print ('')
                                os.remove(fileout)
                    else:
                        print ('Combined output files already exist, or the directory does not exist')
                        print ('')

    return dummy

def loadPA(root_dir, out_dir, date_dir, model_flag):

    '''
    Load in PA stream, combined pp file (suffix _r0.pp)
    '''
    model = ['_HighArctic_1p5km_','_glm']
    stream = '_pa000_r0'

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop
    #
    # def hourly_data(cell):
    #    # return True or False as to whether the cell in question should be kept
    #    # in this case, should return hourly data
    #    return cell >= 30


    # time_hourly_with_stash = iris.AttributeConstraint(STASH=lambda stash: str(stash) in GlobalStashList) & iris.Constraint(time=hourly_data)
    # time_hourly = iris.Constraint(time=hourly_data)

    cube = {}
    cubea = {}
    if model_flag == 0:
        model_expt = model[model_flag] + expt
    elif model_flag == 1:
        model_expt = model[model_flag]

    print(root_dir)
    print(out_dir)
    print(model_expt)
    print(stream)

    for date in date_dir:
        if date[:4] == '2018':
            print(date)
            filename = root_dir + out_dir + date + '/' + date + model_expt + stream + '.pp'
            print (filename)

            cubea[date] = iris.load(filename, 'm01s01i505')
            # cubea[date] = date
            #
            # for i in range(0, len(cube[date])):
            #     ### only load swath variables (size n=94)
            #      if np.size(cube[date][i].dim_coords[1],0) <= 100.:
            #          # print (cube[date][i].dim_coords[1])
            #          cubea[date] = cube[date][i]


    return cubea

def loadPD(root_dir, out_dir, date_dir):

    '''
    Load in PD stream, combined pp file (suffix _r0.pp)
    '''
    model = ['_HighArctic_1p5km_','_glm']
    stream = '_pd000_r0'

    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    # def hourly_data(cell):
    #    # return True or False as to whether the cell in question should be kept
    #    # in this case, should return hourly data
    #    return cell >= 30


    # time_hourly_with_stash = iris.AttributeConstraint(STASH=lambda stash: str(stash) in GlobalStashList) & iris.Constraint(time=hourly_data)
    # time_hourly = iris.Constraint(time=hourly_data)

    # cube = {}
    cubec = {}
    if model_flag == 0:
        model_expt = model[model_flag] + expt
    elif model_flag == 1:
        model_expt = model[model_flag]

    for date in date_dir:
        filename = root_dir + out_dir + date + '/' + date + model_expt + stream + '.pp'

        cubec[date] = iris.load(filename, 'sea_ice_area_fraction') #global_con, callback)
        # cubea[date] = date
        #
        # for i in range(0, len(cube[date])):
        #     ### only load swath variables (size n=94)
        #      if np.size(cube[date][i].dim_coords[1],0) <= 100.:
        #          # print (cube[date][i].dim_coords[1])
        #          cubea[date] = cube[date][i]


    return cubec

def fixHeight(data, cube):

    print ('******')
    print ('')
    print ('Adjusting height to common vertical grid...')
    print ('')

    # height = cube[1].aux_coords[2].points.data       ### 71 levels

    ### wind fields have Z[0] == 2.5
    ### all other 4D fields have Z[0] >= 5.0

    if np.round(cube.aux_coords[2][0].points) > 3:
    # if np.round(cube.aux_coords[2][0].points) == 5:
        ### making 70 levels into 71 for common grid
        cubedata = np.zeros([71,35])
        cubedata[1:,:] = data
        cubedata[0,:] = np.nan
    elif np.round(cube.aux_coords[2][0].points) == 2:
        ### interpolating to n71 common grid
        ### upper bounds = cube[8].aux_coords[2].bounds[:,1]
        cubedata = np.zeros([71,35])
        for i in range(0,35):
            temp = np.interp(cube.aux_coords[2].bounds[:,1],cube.aux_coords[2].points,data[:,i])
            cubedata[1:,i] = temp
            cubedata[0,i] = np.nan
    else:
        cubedata = data

    return cubedata

def pull36HTrack_CloudNet(cube, grid_filename, con, stream, date, model, ship_data, nc_outfile):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt
    import pandas as pd

    print ('******')
    print ('')
    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    if len(cube[0].dim_coords[-1].points) == 25:
        xoffset = -239
        yoffset = -470
    elif len(cube[0].dim_coords[-1].points) == 56:
        xoffset = -210
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 94:
        xoffset = -211
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 81:          ### 14th and 24th August
        xoffset = -209
        yoffset = -399
    elif len(cube[0].dim_coords[-1].points) == 380:         ### needs checked
        xoffset = -60
        yoffset = -110
    else:
        ### if glm
        xoffset = 0
        yoffset = 0

    print ('Because cube shape = ', str(len(cube[0].dim_coords[-1].points)))
    print ('xoffset = ', xoffset)
    print ('yoffset = ', yoffset)

    #################################################################
    ## load gridded ship track
    #################################################################
    # print '******'
    print ('')
    print ('Pulling gridded track from cube:')
    print ('')

    tim, ilat, ilon = readGriddedTrack(grid_filename)
    # if model == 'lam':
    # elif model == 'glm':
    #     tim, ilat, ilon = readGlobal(cube, ship_data, date)
    # else:
    #     print ('Model option is not valid')

    #################################################################
    ## fix time index
    #################################################################

    if np.size(cube)>1:
        print ('')
        print ('More than one variable constraint. Proceeding...')
        print ('')
        print (np.size(cube))

        #################################################################
        ## CREATE EMPTY CUBE FOR PC COLUMN DIAGNOSTICS
        #################################################################
        ncube = Cube(np.zeros([np.size(cube),70,36]))

        #################################################################
        ## POPULATE CUBE WITH EACH INDIVIDUAL VARIABLE
        #################################################################
        for k in range(0,np.size(cube)):            ### loop over number of variables
            print ('')
            print ('k = ', k) ###', so processing', con[k]   # doesn't work with global_con
            print ('')
            #################################################################
            ## only consider hourly diagnostics
            #################################################################
            if len(np.round(cube[k].coord('forecast_period').points)) > 13:
                ###---------------------------------
                ### CHECK IF OFFSETS NEED TO BE RE-ADJUSTED
                ###---------------------------------
                print ('Double-checking grid:')
                if len(cube[k].dim_coords[-1].points) == 25:
                # if cube[0,0].shape >= 25-1:    # ll = 240, 471
                    xoffset = -239
                    yoffset = -470
                elif len(cube[k].dim_coords[-1].points) == 56:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -210
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 94:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -211
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 81:          ### 14th and 24th August
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -209
                    yoffset = -399
                elif len(cube[k].dim_coords[-1].points) == 380:         ### needs checked
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -60
                    yoffset = -110
                else:
                # elif cube[0,0].shape >= 500-1:
                    xoffset = 0
                    yoffset = 0
                print ('Offsets are: ' + str(xoffset) + ', ' + str(yoffset))

                #################################################################
                ## make hourly time array
                #################################################################
                print ('cubetime = ' + str(np.size(cube[k].coord('forecast_period').points)))

                cubetime = np.round(cube[k].coord('forecast_period').points)      ### full 36H forecast period

                #################################################################
                ## PROBE VARIABLE
                #################################################################
                ### do we want to average exluding zeros?
                stash_flag, stash = excludeZeros(cube[k])

                ### do we need to re-grid?  -- DOESN'T WORK LIKE WRF, GRID NOT SPACED SAME WAY
                # cube[k], wind_stash = checkWind(cube[k])

                #################################################################
                ## CHECK DIMENSIONS
                #################################################################
                if np.logical_and(np.size(cube[k].data,1) > 68, np.size(cube[k].data,1) < 72):
                    print ('Variable is 4D:')
                    print ('')
                    #### create empty arrays to be filled
                    data = np.zeros([len(cube[k].coord('model_level_number').points),len(cubetime)-1])
                    ### make dimension flag
                    dim_flag = 1        ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')
                else:
                    print ('Variable is 3D:')
                    print ('')
                    #### create empty arrays to be filled
                    if stream[1:3] == 'pb':
                        if cube[k].long_name == 'large_scale_ice_water_path':
                            data = np.zeros([len(cubetime)])
                        elif cube[k].long_name == 'large_scale_liquid_water_path':
                            data = np.zeros([len(cubetime)])
                        else:
                            data = np.zeros([len(cubetime)-1])
                    elif stream[1:3] == 'pa':
                        if len(cubetime) == 36:
                            data = np.zeros([len(cubetime)-1])
                        elif len(cubetime) == 35:
                            data = np.zeros([len(cubetime)])
                    elif stream[1:3] == 'pd':
                        data = np.zeros([len(cubetime)-1])
                    dim_flag = 0       ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')

                #################################################################
                ## LOOP OVER TIME INDEX, DECOMPOSE ONTO 24H TIMESERIES
                #################################################################
                for j in range(0,len(cubetime)-1):              ### loop over time
                    if j < len(cubetime[:-1]):
                        itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
                    else:
                        ### end point (23h)
                        itime = np.where(tim >= cubetime[-1])
                    # print ''
                    print ('For ', str(j), 'h, itime = ', itime)
                    if dim_flag == 1: dat = np.zeros([len(cube[k].coord('model_level_number').points),len(itime[0])])
                    if dim_flag == 0: dat = np.zeros([len(itime[0])])
                    for i in range(0, len(itime[0])):                   ### loop over time gridded by ship track
                        if np.size(itime) > 1:
                            # print 'Processing i = ', str(itime[0][i])
                            # print '...'
                            if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                            if dim_flag == 0: temp = cube[k][j,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                        else:
                            # print 'Processing i = ', str(itime[i])
                            # print '...'
                            if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                            if dim_flag == 0: temp = cube[k][j,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                        if dim_flag == 1: dat[:,i] = np.squeeze(temp.data)
                        if dim_flag == 0: dat[i] = np.squeeze(temp.data)
                        if np.size(itime) > 1:
                            if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                            if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                            if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                            # print 'averaging over itime ...'
                            # print ''
                        else:
                            if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                            if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                            # print 'no averaging, itime = 1 ...'
                            print ('')
                    # print data
            # print 'data.shape = ', data.shape

            #################################################################
            ## FIGURES TO TEST OUTPUT
            #################################################################
            ### timeseries of lowest model level
            # plt.figure(figsize=(7,5))
            # plt.plot(cubetime[:-1],data[0:10,:])
            # plt.show()

            ### vertical profile of 1st timestep
            # plt.figure(figsize=(7,5))
            # plt.plot(data[:,0],cube.coord('model_level_number').points)
            # plt.show()

            ### pcolormesh of timeseries
            # plt.figure(figsize=(7,5))
            # plt.pcolormesh(cubetime[:-1], cube.coord('model_level_number').points, data)
            # plt.colorbar()
            # plt.show()

            #################################################################
            ## CREATE CUBE
            #################################################################
            ### ECMWF FIELD NAMES
            # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
            #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
            #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};
                varname = varnames.findfieldName(stash)
                print ('standard_name = ', cube[k].standard_name)
                print ('long name = ', cube[k].long_name)
                print ('varname = ', varname)
                print ('')

                if stream[1:3] == 'pa':
                    a = len(cube[k].aux_coords)
                    for ft in range(0,a):
                        print(cube[k].aux_coords[ft].standard_name)
                        if cube[k].aux_coords[ft].standard_name == 'forecast_period':
                            if np.size(cube[k].aux_coords[ft].points) > 35:          ## accounts for arrays with 25 timesteps (includes t12 and t36)
                                ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                            else:
                                ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                else:
                    if cube[k].long_name == 'large_scale_ice_water_path':
                        ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                    elif cube[k].long_name == 'large_scale_liquid_water_path':
                        ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                    else:
                        ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                print (ntime.shape)
                if dim_flag == 1:         ### 4D VARIABLE
                    if stream[1:3] == 'pd':
                        model_height = DimCoord(cube[k].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
                        comdata = data                    #### leave BL diagnostics on RHO levels
                    else:
                        model_height = DimCoord(cube[1].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
                        comdata = fixHeight(data, cube[k])
                    ncube = Cube(np.transpose(comdata),
                            dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
                            standard_name = cube[k].standard_name,
                            long_name = cube[k].long_name,
                            units = cube[k].units,
                            var_name = varname,
                            attributes = cube[k].attributes,
                            aux_coords_and_dims = None,
                            )
                elif dim_flag == 0:         ### 3D VARIABLE
                    ncube = Cube(np.transpose(data),
                            dim_coords_and_dims=[(ntime, 0)],
                            standard_name = cube[k].standard_name,
                            long_name = cube[k].long_name,
                            units = cube[k].units,
                            var_name = varname,
                            attributes = cube[k].attributes,
                            aux_coords_and_dims = None,
                            )
                # ncube.attributes = cube[k].attributes
                # iris.save(ncube, pp_outfile, append=True)
                if k == 0:
                    print ('Initialising fcube')
                    print ('')
                    fcube = [ncube]
                else:
                    print ('Appending variable to fcube')
                    print ('')
                    fcube.append(ncube)

        print (fcube)

    #################################################################
    ## define output filename
    #################################################################
    # print 'Define pp stream outfile:'
    # pp_outfile = date[:6] + str(int(date[6:8])+1) + '_oden_metum_' + str(stream[2:3]) + '.pp'
    # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
    ### bespoke setup if dir is 20180831T1200Z (for 20180901 data)
    # if date == '20180831T1200Z': nc_outfile = '20180901_oden_metum.nc'
    # print 'Outfile = ', pp_outfile

    ### save cube to netcdf file
    print ('')
    print ('Writing fcube to file:')
    print ('')
    if stream[1:3] == 'pc':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pc stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('')
        if not os.path.exists(nc_outfile):
            if 'fcube' in locals():
                out = writeNetCDF(date, fcube, nc_outfile)
            # if PC outfile already exists, combine other stream data
            # if PC outfile doesn't exist, write new

    if stream[1:3] == 'pd':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pd stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('***file is merged to outfile later***')
        print ('')
        doutfile = nc_outfile[:-3] + '_d.nc'
        if not os.path.exists(doutfile):
            if 'fcube' in locals():
                out = writePD_BL(fcube, doutfile)
            # if PD outfile already exists, combine other stream data
            # if PD outfile doesn't exist, write new

    if stream[1:3] == 'pe':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pd stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('***file is merged to outfile later***')
        print ('')
        eoutfile = nc_outfile[:-3] + '_e.nc'
        if not os.path.exists(eoutfile):
            if 'fcube' in locals():
                out = writeFile_netCDF4(fcube, eoutfile)
            # if PC outfile already exists, combine other stream data
            # if PC outfile doesn't exist, write new

    elif stream[1:3] == 'pb':
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
        print ('***file is merged to outfile later***')
        print ('')
        ## Next, append 1D timeseries (surface) data (pb stream)
        ## Can't use Iris for this as cubes can't be 1D
        ##              -> uses standard netCDF appending function
        boutfile = nc_outfile[:-3] + '_b.nc'
        if not os.path.exists(boutfile):
            if 'fcube' in locals():
                out = writePB_Cloudnet(fcube, boutfile)     ##!!!! NEEDS UPDATING TO ONLY WRITE VARIABLES IN FILE, NOT HARD CODED

    elif stream[1:3] == 'pa':
        print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
        print ('***file is merged to outfile later***')
        print ('')
        ## Next, append 1D timeseries (surface) data (pb stream)
        ## Can't use Iris for this as cubes can't be 1D
        ##              -> uses standard netCDF appending function
        aoutfile = nc_outfile[:-3] + '_a.nc'
        if not os.path.exists(aoutfile):
            if 'fcube' in locals():
                out = writePA_Analysis(fcube, aoutfile)

    return nc_outfile

def writeNetCDF(date, cube, nc_outfile):

    #################################################################
    ## CREATE NETCDF
    #################################################################
    #################################################################
    ## define output filename
    #################################################################
    print ('******')
    print ('Define .nc stream outfile:')
    # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
    print ('Outfile will be = ', nc_outfile)

    #################################################################
    ## load in each stream
    #################################################################
    ### USE IRIS TO SAVE OUT PC CUBE TO NETCDF (CREATING NEW FILE):
    # -------------------------------------------------------------
    # Convert .pp to .nc
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Converting to netCDF:')
    print ('')
    # cube = iris.load(outfile[0], global_con, callback)
    iris.save(cube, nc_outfile)

    return nc_outfile

def writePB_Cloudnet(cube, boutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 1D data to ' + boutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(boutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    print (cube)

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 1200Z UTC initialisation.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write pbXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
        dat.scale_factor = float(1)
        dat.add_offset = float(0)
        dat.units = str(cube[d].units)
        dat.STASH = str(cube[d].attributes['STASH'])
        if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
        if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
        dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writePA_Analysis(cube, aoutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 1D data to ' + aoutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(aoutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    print (cube)
    # print cube[0].dim_coords

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))     ## use net sw to define, 1st in fcube

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 1200Z UTC initialisation.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write paXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        if np.size(cube[d].dim_coords[0],0) > 24:      ### ignore 3-hourly data for now
            print ('Writing ' + cube[d].var_name)
            print ('')
            if not cube[d].var_name in dataset.variables:
                dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                dat.units = str(cube[d].units)
                dat.STASH = str(cube[d].attributes['STASH'])
                if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
                if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
                dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writePD_BL(cube, doutfile):
    #################################################################
    ## Write boundary layer diagnosticsa (PD) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing 3D data to ' + doutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(doutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')
    print ('Cube is: ')
    print (cube)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    #### find first occurrence of 2D variable, then break
    for l in range(0,len(cube)):
        if np.ndim(cube[l]) == 2:
            lind = l
            print ('height dim based on ' )
            print (cube[l])
            break

    forecast_time = dataset.createDimension('forecast_time', np.size(cube[lind].dim_coords[0].points))
    height = dataset.createDimension('height', np.size(cube[lind].dim_coords[1].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 1200Z UTC initialisation.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[lind].dim_coords[0].points      ### forecast time

    #### height
    height = dataset.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = ''
    height.units = 'm'
    height.long_name = 'height'
    height[:] = cube[lind].dim_coords[1].points      ### forecast time

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        print
        if np.ndim(cube[d]) == 2:
            if cube[d].var_name == 'air_pressure': continue
            dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time','height',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            dat.units = str(cube[d].units)
            dat.STASH = str(cube[d].attributes['STASH'])
            if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
            if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
            # if np.size(cube[d].data,1) == 70:         ### need this for UM_RA2T_TKEdissrate
            #     dat[:,:] = cube[d].data
            # elif np.size(cube[d].data,1) == 69:
            #     dat[:,:-1] = cube[d].data
            #     dat[:,-1] = np.nan
            dat[:,:] = cube[d].data
        elif np.ndim(cube[d]) == 1:
            dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            dat.units = str(cube[d].units)
            dat.STASH = str(cube[d].attributes['STASH'])
            if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
            if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
            dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writeFile_netCDF4(cube, eoutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 3D data to ' + eoutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(eoutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    # print cube
    # print cube[0].dim_coords

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))
    height = dataset.createDimension('height', np.size(cube[0].dim_coords[1].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 1200Z UTC initialisation.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time (ignore first 12h)

    #### height
    height = dataset.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = ''
    height.units = 'm'
    height.long_name = 'height'
    height[:] = cube[0].dim_coords[1].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write paXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time','height',), fill_value='-9999')
        dat.scale_factor = float(1)
        dat.add_offset = float(0)
        dat.units = str(cube[d].units)
        dat.STASH = str(cube[d].attributes['STASH'])
        if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
        if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
        dat[:,:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def appendMetaNetCDF(outfile, date, out_dir, model):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Appending metadata to ' + outfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(outfile, 'a', format ='NETCDF4_CLASSIC')
    # infile = net.Dataset("2015%s%s-160000_0.nc" % (month,day), "a")
    # print ''
    # print dataset.file_format
    # print ''

    ###################################
    ## Global Attributes
    ###################################
    dataset.title = 'Met Office Unified Model single-site (Oden) output during MOCCHA'
    revision = 'undefined'
    if out_dir[3:10] == 'u-cc278':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). '
        revision = 'Revision no. 2. '
    elif out_dir[3:10] == 'u-cc324':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. Also set l_subgrid_qcl_mp to .true. to allow for turbulent production of mixed-phase cloud. Extended BL diagnostic list. '
        revision = 'Revision no. 2. '
    elif out_dir[3:10] == 'u-cc568':
        micro = 'Cloud microphysics: Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983). Extended BL diagnostic list. Updated revision of suite u-bg610. '
        revision = 'Revision no. 2. '
    else:
        micro = '<MICROPHYSICS UNDEFINED IN META>'
    wind = 'U and V wind components interpolated on to common vertical grid. '
    forecast_notes = 'Full 36-H forecast, rather than daily files used for main body of study (Young et al., 2021; ACPD). '
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young McCusker <G.Y.McCusker@leeds.ac.uk> using Python (Iris/netCDF4).'
    # dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'Rose suite ID: ' + out_dir[2:10]
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic. '
    if model == 'lam':
        modelnote = 'UM limited area model. '
        dataset.description = 'Hourly data taken from grid box closest to ship location. Where the ship covers more than one grid box within an hour period, data are averaged from all grid boxes crossed. '
    elif model == 'glm':
        modelnote = 'UM global model. '
        dataset.description = 'Hourly data taken from grid box closest to ship location. '
    dataset.comment = revision + modelnote + forecast_notes + micro + wind
    dataset.institution = 'University of Leeds.'
    dataset.initialization_time = date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' ' + date[9:14] + '.'

    ###################################
    ## Additional variables
    ###################################
    #### Model resolution
    if not 'horizontal_resolution' in dataset.variables.keys():
        if model == 'lam':
            res = dataset.createVariable('horizontal_resolution', np.float32, fill_value='-9999')
            res.comment = 'Horizontal grid size of nested region.'
            res.units = 'km'
            res[:] = 1.5
        elif model == 'glm':
            res = dataset.createVariable('horizontal_resolution', np.float32, fill_value='-9999')
            res.comment = 'Horizontal grid size of global model.'
            res.units = 'km'
            res[:] = 17.0

    ###################################
    ## Open pbXXX netCDF file
    ###################################
    boutfile = outfile[:-3] + '_b.nc'

    if os.path.exists(boutfile):
        ncB = Dataset(boutfile, 'r')

        ###################################
        ## Append pbXXX stream diagnostics
        ###################################

        print ('Appending pbXXX diagnostics:')
        print ('---')
        for d in ncB.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                dat = dataset.createVariable(d, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                if getattr(ncB.variables[d],'units', None):
                    dat.units = str(ncB.variables[d].units)
                else:
                    dat.units = 'unknown'
                if getattr(ncB.variables[d],'STASH', None):
                    dat.STASH = str(ncB.variables[d].STASH)
                if getattr(ncB.variables[d],'standard_name', None):
                    dat.standard_name = str(ncB.variables[d].standard_name)
                if getattr(ncB.variables[d],'long_name', None):
                    dat.long_name = str(ncB.variables[d].long_name)
                dat[:] = ncB.variables[d][:]

        ###################################
        ## Close read-only pbXXX file
        ###################################
        ncB.close()

    ###################################
    ## Open paXXX netCDF file
    ###################################
    aoutfile = outfile[:-3] + '_a.nc'

    if os.path.exists(aoutfile):
        ncA = Dataset(aoutfile, 'r')

        ###################################
        ## Append paXXX stream diagnostics
        ###################################
        print ('Appending paXXX diagnostics:')
        print ('---')
        for d in ncA.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                dat = dataset.createVariable(d, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                if getattr(ncA.variables[d],'units', None):
                    dat.units = str(ncA.variables[d].units)
                else:
                    dat.units = 'unknown'
                if getattr(ncA.variables[d],'STASH', None):
                    dat.STASH = str(ncA.variables[d].STASH)
                if getattr(ncA.variables[d],'standard_name', None):
                    dat.standard_name = str(ncA.variables[d].standard_name)
                if getattr(ncA.variables[d],'long_name', None):
                    dat.long_name = str(ncA.variables[d].long_name)
                dat[:] = ncA.variables[d][:]

        ###################################
        ## Close read-only paXXX file
        ###################################
        ncA.close()

    ###################################
    ## Open pdXXX netCDF file
    ###################################
    doutfile = outfile[:-3] + '_d.nc'

    print ('What variables do we have before pdXXX read in?:')
    print (dataset)

    if os.path.exists(doutfile):
        ncD = Dataset(doutfile, 'r')

        print (ncD)

        #### height on rho levels
        height2 = dataset.createDimension('height2', np.size(ncD.variables['height'][:]))
        height2 = dataset.createVariable('height2', np.float64, ('height2',), fill_value='-9999')
        height2.scale_factor = float(1)
        height2.add_offset = float(0)
        height2.comment = 'height coordinate on rho levels'
        height2.units = 'm'
        height2.long_name = 'height'
        height2[:] = ncD.variables['height'][:]      ### forecast time (ignore first 12h)

        ###################################
        ## Append pdXXX stream diagnostics
        ###################################
        print ('Appending pdXXX diagnostics:')
        print ('---')
        for d in ncD.variables:
            print (d)
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                if np.ndim(ncD.variables[d]) == 2:
                    print ('Variable is 2D:')
                    daat = dataset.createVariable(d, np.float64, ('forecast_time', 'height2',), fill_value='-9999')
                    daat.scale_factor = float(1)
                    daat.add_offset = float(0)
                    if getattr(ncD.variables[d],'units', None):
                        daat.units = str(ncD.variables[d].units)
                    else:
                        daat.units = 'unknown'
                    if getattr(ncD.variables[d],'STASH', None):
                        daat.STASH = str(ncD.variables[d].STASH)
                    if getattr(ncD.variables[d],'standard_name', None):
                        daat.standard_name = str(ncD.variables[d].standard_name)
                    if getattr(ncD.variables[d],'long_name', None):
                        daat.long_name = str(ncD.variables[d].long_name)
                    daat[:,:] = ncD.variables[d][:,:]
                elif np.ndim(ncD.variables[d]) == 1:
                    print ('Variable is 1D:')
                    dat = dataset.createVariable(d, np.float64, ('forecast_time', ), fill_value='-9999')
                    dat.scale_factor = float(1)
                    dat.add_offset = float(0)
                    if getattr(ncD.variables[d],'units', None):
                        dat.units = str(ncD.variables[d].units)
                    else:
                        dat.units = 'unknown'
                    if getattr(ncD.variables[d],'STASH', None):
                        dat.STASH = str(ncD.variables[d].STASH)
                    if getattr(ncD.variables[d],'standard_name', None):
                        dat.standard_name = str(ncD.variables[d].standard_name)
                    if getattr(ncD.variables[d],'long_name', None):
                        dat.long_name = str(ncD.variables[d].long_name)
                    dat[:] = ncD.variables[d][:]

        ###################################
        ## Close read-only pdXXX file
        ###################################
        ncD.close()

    ###################################
    ## Open peXXX netCDF file
    ###################################
    eoutfile = outfile[:-3] + '_e.nc'

    if os.path.exists(eoutfile):
        ncE = Dataset(eoutfile, 'r')

        ###################################
        ## Append paXXX stream diagnostics
        ###################################
        print ('Appending peXXX diagnostics:')
        print ('---')
        for d in ncE.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                daat = dataset.createVariable(d, np.float64, ('forecast_time', 'height',), fill_value='-9999')
                daat.scale_factor = float(1)
                daat.add_offset = float(0)
                if getattr(ncE.variables[d],'units', None):
                    daat.units = str(ncE.variables[d].units)
                else:
                    daat.units = 'unknown'
                if getattr(ncE.variables[d],'STASH', None):
                    daat.STASH = str(ncE.variables[d].STASH)
                if getattr(ncE.variables[d],'standard_name', None):
                    daat.standard_name = str(ncE.variables[d].standard_name)
                if getattr(ncE.variables[d],'long_name', None):
                    daat.long_name = str(ncE.variables[d].long_name)
                daat[:,:] = ncE.variables[d][:,:]

        ###################################
        ## Close read-only peXXX file
        ###################################
        ncE.close()

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def excludeZeros(cube):

    print ('')
    print ('Checking stash list:')

    print ('Want to exclude zeros in the following fields:')
    ### list of stash items where we want to exclude zeros
    STASH = ['m01s00i012','m01s00i254','m01s00i075','m01s00i076','m01s00i078',
        'm01s00i079','m01s00i081','m01s00i271','m01s00i273','m01s00i272',
        'm01s00i083','m01s00i084','m01s00i088']
    print (STASH)

    print ('Diag is:')
    str_m = "%02d" % cube.attributes['STASH'][0]
    str_s = "%02d" % cube.attributes['STASH'][1]
    str_i = "%03d" % cube.attributes['STASH'][2]
    stash = str('m' + str_m + 's' + str_s + 'i' + str_i)
    print (stash)

    for i in range(0, len(STASH)):
        if STASH[i] == stash:
            # data[data==0] = np.nan              # set zeros to nans
            flag = 1                           # flagging if in list
            print ('In list, so excluding zeros')
            break
        else:
            flag = 0                           # flagging if not in list
            print ('Not in list, so not excluding zeros')

    # if flag == 1:
    # if flag == 0:
    # print ''

    return flag, stash

def pullTrack(date, root_dir, out_dir, global_con, model, ship_data):

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin cube read in at ' + time.strftime("%c"))
    print (' ')
    # var_con = 'specific_humidity'
    # cube = iris.load_cube(filename1, var_con)
    # global_con = ['atmosphere_downward_eastward_stress','atmosphere_downward_northward_stress']

    grid_dirname = 'AUX_DATA/'
    if model == 'lam': grid_filename = grid_dirname + date + '-36HForecast_ShipTrack_GRIDDED.csv'
    if model == 'glm': grid_filename = grid_dirname + date + '-36HForecast-GLM_ShipTrack_GRIDDED.csv'

    print ('Loading grid:...')
    print (grid_filename)

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    # -------------------------------------------------------------
    # Define output stream filenames to look at:
    #           start at 012 if 3h dumps (a)
    #           start at 009 if 1h dumps in pb
    #           start at 011 if 1h dumps (c--e)
    # -------------------------------------------------------------
    names = ['_pa000_r0','_pb000_r0','_pd000_r0','_pe000_r0','_pc000_r0']
    # names = ['_pa000_r0']         ### only do specific files as a test
    if out_dir[-6:-1] == 'CASIM':
        expt = out_dir[-11:-1]
    elif out_dir[-4:-1] == 'CON':
        expt = out_dir[-9:-1]
    outfiles = [] ### define list to add processed filenames to

    for stream in names:
        ### -------------------------------------------------------------------------
        ### define output filename
        ### -------------------------------------------------------------------------
        if model == 'lam':
            # #### LAM
            filename = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream + '.pp'
            dirout = out_dir[3:10] + '_lam/'
        elif model == 'glm':
            ### GLM
            filename = root_dir + out_dir + date + '/' + date + '_glm' + stream + '.pp'
            dirout = out_dir[3:10] + '_glm/'
        print ('dirout is: ' + dirout)

        # if np.logical_or(np.logical_or(out_dir[0] == '1',out_dir[0] == '2'),out_dir[0] == '3'):
        #     dirout = out_dir[3:10] + '/'
        # else:
        #     dirout = out_dir[2:9] + '/'

        print ('Checking: ' + filename)
        exist_flag = 0 # initialise exist_flag
        if os.path.exists(filename):
            exist_flag = 1
            #### LOAD CUBE
            if 'var_con' in locals():
                print ('Loading single diagnostic:')
                print (var_con)
                cube1 = iris.load_cube(filename, var_con, callback)
                con_flag = 0            # constraint flag
            elif 'global_con' in locals():
                print ('Loading multiple diagnostics:')
                # cube = iris.load_cubes(filename1, global_con)
                cube = iris.load(filename, global_con, callback)
                con_flag = 1            # constraint flag
                print (cube)

            # ------------------------------------------------------------

            ### -------------------------------------------------------------
            ### Use the following to plot quick maps of loaded cubes
            ### -------------------------------------------------------------

            # hour = 0
            # figure = plot_cartmap(ship_data, cube, hour, grid_filename)

            ########################################################################
            ### -------------------------------------------------------------
            ### Pull gridded ship track from cube
            ### -------------------------------------------------------------
            ########################################################################
            # -------------------------------------------------------------
            ### 1. use the following if only want the exact ship position and no variability
            # -------------------------------------------------------------
            ### LOAD CUBE
            nc_outfile = dirout + date[:6] + str(int(date[6:8])).zfill(2) + '-36HForecast_oden_metum.nc'
            aoutfile = nc_outfile[:-3] + '_a.nc'
            boutfile = nc_outfile[:-3] + '_b.nc'
            doutfile = nc_outfile[:-3] + '_d.nc'
            eoutfile = nc_outfile[:-3] + '_e.nc'

            print ('...')
            print (' ')

            if stream[1:3] == 'pa':
                print ('Checking: ' + aoutfile + '...')
                if not os.path.exists(aoutfile):
                    print (aoutfile + ' does not exist, so pulling ship track...')
                    outfile = pull36HTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
            elif stream[1:3] == 'pb':
                print ('Checking: ' + boutfile + '...')
                if not os.path.exists(boutfile):
                    print (boutfile + ' does not exist, so pulling ship track...')
                    outfile = pull36HTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
            elif stream[1:3] == 'pd':
                print ('Checking: ' + doutfile + '...')
                if not os.path.exists(doutfile):
                    print (doutfile + ' does not exist, so pulling ship track...')
                    outfile = pull36HTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
            elif stream[1:3] == 'pe':
                print ('Checking: ' + eoutfile + '...')
                if not os.path.exists(eoutfile):
                    print (eoutfile + ' does not exist, so pulling ship track...')
                    outfile = pull36HTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
            elif stream[1:3] == 'pc':
                print ('Checking: ' + nc_outfile + '...')
                if not os.path.exists(nc_outfile):
                    print (nc_outfile + ' does not exist, so pulling ship track...')
                    outfile = pull36HTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                    print ('******')
                    print ('')
                    print ('stream = ' + stream + ', so appending pa, pb, pd, pe (if present), and metadata')
                    print ('')
                    out = appendMetaNetCDF(nc_outfile, date, out_dir, model)
            else:
                print ('Valid stream not found.')

def loadObservations(obs, platform, obs_root_dir):

    from readMAT import readMatlabStruct
    from time_functions import calcTime_Mat2DOY

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

    if platform == 'LAPTOP':
        print ('Load temporary ice station data from Jutta...')
        obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_V3_30minres.nc','r')
        print ('Load ice station flux data from Jutta...')
        obs['ice_station_fluxes'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')

        print ('Load HATPRO data used by Cloudnet...')
        dirname = '/home/gillian/MOCCHA/MOCCHA_GIT/ODEN/DATA/hatpro/'
        dir = os.listdir(dirname)
        obs['hatpro'] = {}
        for file in dir:
            if file == 'V1': continue
            # print (file)
            IWVtemp = readMatlabStruct(dirname + file)
            # print (IWVtemp.keys())
            if file == '20180814_IWV_30s_V2.mat':       ### if it is the first file
                obs['hatpro']['IWV'] = np.squeeze(IWVtemp['iwv'])
                obs['hatpro']['mday'] = np.squeeze(IWVtemp['mday'])
                obs['hatpro']['LWP'] = np.squeeze(IWVtemp['lwp'])
                obs['hatpro']['rainflag'] = np.squeeze(IWVtemp['rainflag'])
            else:
                obs['hatpro']['IWV'] = np.append(np.squeeze(obs['hatpro']['IWV']),np.squeeze(IWVtemp['iwv']))
                obs['hatpro']['mday'] = np.append(np.squeeze(obs['hatpro']['mday']),np.squeeze(IWVtemp['mday']))
                obs['hatpro']['LWP'] = np.append(np.squeeze(obs['hatpro']['LWP']),np.squeeze(IWVtemp['lwp']))
                obs['hatpro']['rainflag'] = np.append(np.squeeze(obs['hatpro']['rainflag']),np.squeeze(IWVtemp['rainflag']))
        obs['hatpro']['doy'] = calcTime_Mat2DOY(obs['hatpro']['mday'])

        print ('Load albedo estimates from Michael...')
        obs['albedo'] = readMatlabStruct(obs_root_dir + 'MOCCHA_Albedo_estimates_Michael.mat')

        print ('Load MRR rainrate data from Jutta...')
        mrr_dirname = '/home/gillian/MOCCHA/MOCCHA_GIT/ODEN/DATA/MRR/'
        mrr_dir = os.listdir(mrr_dirname)
        obs['mrr'] = {}
        mrr_filename = 'Rainrate_from_jutta.mat'
        obs['mrr'] = readMatlabStruct(mrr_dirname + mrr_filename)
        print (obs['mrr'].keys())
        # i = 0
        # for file in mrr_dir:
        #     doy = np.arange(226,258)
            # mrr = Dataset(mrr_dirname + 'rain/' + file,'r')
            # if file == '20180814_oden_mrr.nc':       ### if it is the first file
            #     obs['mrr']['rain'] = mrr.variables['rain'][:]
            #     obs['mrr']['time'] = mrr.variables['time'][:]/24.0 + doy[0]
            # else:
            #     i = i + 1
            #     obs['mrr']['rain'] = np.append(obs['mrr']['rain'],mrr.variables['rain'][:])
            #     obs['mrr']['time'] = np.append(obs['mrr']['time'],mrr.variables['time'][:]/24.0 + doy[i])

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
    obs['pws'] = readMatlabStruct(obs_root_dir + '7thDeck/ACAS_AO2018_PWD_1min_v1_0.mat')

    print ('...')

    return obs

def loadNCs(nc, data, root_dir, dir, model):

    '''
    Load in netCDF files of 36H forecasts into numpy dictionary (called 'data')
    '''

    filenames = os.listdir(root_dir + dir)

    #### read in netcdf data and initialise postprocessing dictionary
    print ('***')
    print ('Loading in: ' + dir[:2])

    for filename in filenames:
        print (filename + ' + ' + model)
        if model == 'lam':
            nc[dir[:2]][filename[:8]] = Dataset(root_dir + dir + filename)
            data[dir[:2]][filename[:8]] = {}
            # print (nc[dir[:2]][filename[:8]])
        if model == 'glm':
            nc[dir[:2] + '_glm'][filename[:8]] = Dataset(root_dir + dir + filename)
            data[dir[:2] + '_glm'][filename[:8]] = {}
            # print (nc[dir[:2] + '_glm'][filename[:8]])

    return nc, data

def reGrid_Sondes(data, obs, dir, filenames, model_list, var):

    from scipy.interpolate import interp1d
    from readMAT import readMatlabStruct
    from time_functions import calcTime_Mat2DOY

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

    #### ---------------------------------------------------------------
    ### build list of variables names wrt input data [OBS, UM, CASIM, IFS]
    #### ---------------------------------------------------------------
    if var == 'temp':
        varlist = ['temperature','temperature','temperature','temperature', 'temperature', 'temperature']
    elif var == 'thetaE':
        # varlist = ['epottemp','thetaE','thetaE','thetaE']     # use sonde file's epottemp
        varlist = ['thetaE','thetaE','thetaE','thetaE', 'thetaE', 'thetaE']         # use sonde calculated thetaE
    elif var == 'theta':
        # varlist = ['pottemp','theta','theta','theta']     # use sonde file's pottemp
        varlist = ['theta','theta','theta','theta','theta','theta']         # use sonde calculated theta
    elif var == 'q':
        varlist = ['sphum','q','q','q','q','q']
    elif var == 'rh':
        varlist = ['RH','rh','rh','rh','rh','rh']

    #### ---------------------------------------------------------------
    #### save hourly temperature model profiles (using the ii index defined by the time indices)
    #### ---------------------------------------------------------------
    data1[var + '_hrly'] = np.squeeze(data1[varlist[1]][data1['hrly_flag'],:])
    data2[var + '_hrly'] = np.squeeze(data2[varlist[2]][data2['hrly_flag'],:])
    data3[var + '_hrly'] = np.squeeze(data3[varlist[3]][data3['hrly_flag'],:])
    data4[var + '_hrly'] = np.squeeze(data4[varlist[4]][data4['hrly_flag'],:])
    data5[var + '_hrly'] = np.squeeze(data5[varlist[5]][data5['hrly_flag'],:])

    # print(data2[var + '_hrly'].shape)

    #### ---------------------------------------------------------------
    #### explicitly save 6-hourly temperature model profiles and time binning for ease
    #### ---------------------------------------------------------------
    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_6hrly'] = data1['time_hrly'][::6]
    data2['time_6hrly'] = data2['time_hrly'][::6]
    data3['time_6hrly'] = data3['time_hrly'][::6]
    data4['time_6hrly'] = data4['time_hrly'][::6]
    data5['time_6hrly'] = data5['time_hrly'][::6]
    data1[var + '_6hrly'] = data1[var + '_hrly'][::6]
    data2[var + '_6hrly'] = data2[var + '_hrly'][::6]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6]
    data4[var + '_6hrly'] = data4[var + '_hrly'][::6]
    data5[var + '_6hrly'] = data5[var + '_hrly'][::6]

    print(data2[var + '_6hrly'].shape)

    #### ---------------------------------------------------------------
    #### index to only look at altitudes <10km
    #### ---------------------------------------------------------------
    iTim = 0        ### initialised
    iObs = np.where(obs['sondes']['gpsaltitude'][:,iTim] <= 11000)
    iUM = np.where(data1['height'] <= 11000)
    iGLM = np.where(data5['height'] <= 11000)
    if ifs_flag == True:
        iIFS = np.where(data3['height'][iTim,:] <= 11000)
    else:
        iIFS = np.where(data3['height'] <= 11000)

    #### ---------------------------------------------------------------
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    if ifs_flag == True:
        data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems
            ####        further on
        data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    if ifs_flag == True:
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
    else:
        data3['universal_height'] = data3['height'][iUM[0][3:]]
        data3['universal_height_UMindex'] = iUM[0][3:]
        data3['height_6hrly'] = np.zeros([np.size(data3['time_6hrly']), np.size(data3['universal_height_UMindex'])])
        for i in range(0, len(data3['time_6hrly'])): data3['height_6hrly'][i,:] = data3['universal_height'][:]
        data3[var + '_6hrly_UM'] = data3[var + '_hrly'][::6]


    #### INTERPOLATION TESTING:
    # print (data3['temp_hrly_UM'].shape),'radr_refl'
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

    print ('')
    print ('Defining GLM profile as a function:')
    data5[var + '_hrly_UM'] = np.zeros([np.size(data5['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(data5['time_hrly'],0)):
        # print (iTim)
        fnct_GLM = interp1d(np.squeeze(data5['height'][iGLM]), np.squeeze(data5[var + '_hrly'][iTim,iGLM]))
        data5[var + '_hrly_UM'][iTim,:] = fnct_GLM(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('GLM(UM Grid) function worked!')
    print (var + ' GLM data now on UM(lam) vertical grid')
    print ('*****')
    ### assign for easier indexing later
    data5[var + '_6hrly_UM'] = data5[var + '_hrly_UM'][::6,:]

    print (np.size(data5[var + '_hrly_UM'],0))
    print (np.size(data5[var + '_6hrly_UM'],0))

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['sondes']['doy'] >= 225.9, obs['sondes']['doy'] <= 258.0))
    subset = np.where(np.logical_and(obs['sondes']['doy'] >= doy[0], obs['sondes']['doy'] <= doy[-1] + 0.05))
    # drift = np.where(np.logical_and(obs['sondes']['doy'] >= doy[0], obs['sondes']['doy'] <= doy[-1] + 0.05))

    # print (obs['sondes']['doy'][drift[0]])
    # print (obs['sondes']['doy'][drift[-1]])

    ### save in dict for ease
    obs['sondes']['doy_drift'] = obs['sondes']['doy'][drift]
    obs['sondes']['drift'] = drift
    obs['sondes'][var + '_driftSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][drift[0],:]

    # print (obs['sondes'][var + '_driftSondes_UM'].shape)
    # print (data2[var + '_6hrly'].shape)
    # print (data4[var + '_6hrly'].shape)

    ### save subset of drift in dict for ease
    obs['sondes']['doy_subset'] = obs['sondes']['doy'][subset]
    obs['sondes']['subset'] = subset
    obs['sondes'][var + '_subsetSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][subset[0],:]

    # #### INTERPOLATION TESTING - IFS + SONDE + UM_RA2M:
    # print (obs['sondes']['doy_drift'].shape)
    # print (obs['sondes']['temp_allSondes_UM'][drift[0],:].shape)
    # if var == 'temp':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,drift[0][i]]) + 273.15,np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['temp_driftSondes_UM'][i,:] + 273.15,data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = '#FFC107', label = 'ifs-Zindexed')
    #         plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = '#FFC107', label = 'ifs-interpd')
    #         plt.plot(data1['temp_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2['temp_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.plot(data4['temp_6hrly'][i,iUM[0][3:]], data4['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2t')
    #         plt.plot(np.squeeze(data5['temp_6hrly'][i,iGLM]), data5['height'][iGLM], '--', color = 'grey', label = 'um_glm')
    #         plt.plot(data5['temp_6hrly_UM'][i,:], data1['height'][iUM[0][3:]], color = 'grey', label = 'um_glm-interpd')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)) + '\n Model time = ' + str(data1['time_6hrly'][i]))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid_v3/REGRID_Ttest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'q':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['mr'][iObs,drift[0][i]]), np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes'][var + '_driftSondes_UM'][i,:], data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3[var + '_6hrly'][i,iIFS])*1e3,np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = '#FFC107', label = 'ifs-Zindexed')
    #         plt.plot(data3[var + '_6hrly_UM'][i,:]*1e3,data1['height'][iUM[0][3:]], color = '#FFC107', label = 'ifs-interpd')
    #         plt.plot(data1[var + '_6hrly'][i,iUM[0][3:]]*1e3, data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2[var + '_6hrly'][i,iUM[0][3:]]*1e3, data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.plot(data4[var + '_6hrly'][i,iUM[0][3:]]*1e3, data4['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2t')
    #         plt.plot(np.squeeze(data5['q_6hrly'][i,iGLM])*1e3, data5['height'][iGLM], '--', color = 'grey', label = 'um_glm')
    #         plt.plot(data5['q_6hrly_UM'][i,:]*1e3, data1['height'][iUM[0][3:]], color = 'grey', label = 'um_glm-interpd')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)) + '\n Model time = ' + str(data1['time_6hrly'][i]))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid_v3/REGRID_Qtest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'thetaE':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['thetaE'][iObs,drift[0][i]]),np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['thetaE_driftSondes_UM'][i,:],data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['thetaE_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = '#FFC107', label = 'ifs-Zindexed')
    #         plt.plot(data3['thetaE_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = '#FFC107', label = 'ifs-interpd')
    #         plt.plot(data1['thetaE_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2['thetaE_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.title('REGRID test DOY ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.xlabel('$\Theta_{E}$ [K]')
    #         plt.ylabel('Z [m]')
    #         plt.ylim([0,3000])
    #         plt.xlim([260,320])
    #         plt.legend()
    #         plt.savefig('../FIGS/inversionIdent/REGRID_ThetaE_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()

    #### ---------------------------------------------------------------
    #### make some dictionary assignments for use later
    #### ---------------------------------------------------------------
    data1['universal_height'] = data1['height'][iUM[0][3:]]
    data1['universal_height_UMindex'] = iUM[0][3:]

    return data, obs

def radiosondeAnalysis(nc, data, dir, obs, filenames, model_list):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from readMAT import readMatlabStruct
    from time_functions import calcTime_Mat2DOY

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    # obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))
    print (obs['sondes'].keys())

    ### for reference in figures
    print()
    zeros = np.zeros(len(nc[dir[:2]][filenames[0][:8]]['forecast_time']))

    #### set flagged values to nans
    for filename in filenames:
        for model in model_list:
            if model == 'lam':
                data[dir[:2]][filename[:8]]['temperature'] = nc[dir[:2]][filename[:8]]['temperature'][:]
                data[dir[:2]][filename[:8]]['q'] = nc[dir[:2]][filename[:8]]['q'][:]
                data[dir[:2]][filename[:8]]['temperature'][data[dir[:2]][filename[:8]]['temperature'] == -9999] = np.nan
                data[dir[:2]][filename[:8]]['temperature'][data[dir[:2]][filename[:8]]['temperature'] <= 0] = np.nan
                data[dir[:2]][filename[:8]]['q'][data[dir[:2]][filename[:8]]['q'] == -9999] = np.nan
                data[dir[:2]][filename[:8]]['q'][data[dir[:2]][filename[:8]]['q'] <= 0] = np.nan
            elif model == 'glm':
                if dir[:2] + '_glm' in data.keys():
                    data[dir[:2] + '_glm'][filename[:8]]['temperature'] = nc[dir[:2] + '_glm'][filename[:8]]['temperature'][:]
                    data[dir[:2] + '_glm'][filename[:8]]['q'] = nc[dir[:2] + '_glm'][filename[:8]]['q'][:]
                    data[dir[:2] + '_glm'][filename[:8]]['temperature'][data[dir[:2] + '_glm'][filename[:8]]['temperature'] == -9999] = np.nan
                    data[dir[:2] + '_glm'][filename[:8]]['temperature'][data[dir[:2] + '_glm'][filename[:8]]['temperature'] <= 0] = np.nan
                    data[dir[:2] + '_glm'][filename[:8]]['q'][data[dir[:2] + '_glm'][filename[:8]]['q'] == -9999] = np.nan
                    data[dir[:2] + '_glm'][filename[:8]]['q'][data[dir[:2] + '_glm'][filename[:8]]['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    # data, obs = reGrid_Sondes(data, obs, dir, filenames, model_list, 'temp')
    # data, obs = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    # Tmin = -45
    # Tmax = 5
    # ymax = 9000
    # qmax = 4.0
    #
    # data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    #
    # data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])

    # ##################################################
    # ##################################################
    # #### MEDIAN PROFILES
    # ##################################################
    # ##################################################
    # SMALL_SIZE = 12
    # MED_SIZE = 16
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    #
    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(7.5,5.5))
    # plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
    #         hspace = 0.3, wspace = 0.2)
    #
    # ####        all model data share a timestamp
    # melt = np.where(data1['time_hrly'][::6] < 240.0)
    # freeze = np.where(data1['time_hrly'][::6] >= 240.0)
    #
    #
    # ###-------------------------
    # plt.subplot(121)
    # ax1 = plt.gca()
    # # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
    #     np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
    #     np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
    #      color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
    #     np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
    #     np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = '#FFC107', linewidth = 0.5)
    # plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = '#FFC107', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = '#FFC107', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2[:-8], zorder = 1)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    # plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # plt.legend(bbox_to_anchor=(0.9, 1.03, 1., .102), loc=4, ncol=2)
    # plt.legend(bbox_to_anchor=(1.25, 1.03, 1., .102), loc=4, ncol=3)
    # plt.ylabel('Z [km]')
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-2.0,1.0])
    # plt.xlabel('T bias [K]')
    #
    # ###-------------------------
    # plt.subplot(122)
    # ax1 = plt.gca()
    # # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
    #     np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
    #     np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
    #     np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
    #     np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    #     '--', color = '#FFC107', linewidth = 0.5)
    # plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    #     '--', color = '#FFC107', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = '#FFC107', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    # plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # plt.legend()
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    # plt.grid('on')
    #
    # ###-------------------------
    # fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_ifs_casim-100_ra2t_ra2m-25_wUMGlobal.svg'
    # # plt.savefig(fileout)
    # plt.show()

    # ##################################################
    # ##################################################
    # #### MEDIAN PROFILES
    # ##################################################
    # ##################################################
    # SMALL_SIZE = 12
    # MED_SIZE = 16
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    #
    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(7.5,5.5))
    # plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
    #         hspace = 0.3, wspace = 0.2)
    #
    # ####        all model data share a timestamp
    # melt = np.where(data1['time_hrly'][::6] < 240.0)
    # freeze = np.where(data1['time_hrly'][::6] >= 240.0)
    #
    #
    # ###-------------------------
    # plt.subplot(121)
    # ax1 = plt.gca()
    # plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkorange', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'steelblue', label = label2[:-8], zorder = 1)
    # plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # plt.legend(bbox_to_anchor=(0.9, 1.03, 1., .102), loc=4, ncol=2)
    # plt.legend(bbox_to_anchor=(1.25, 1.03, 1., .102), loc=4, ncol=3)
    # plt.ylabel('Z [km]')
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-2.0,1.0])
    # plt.xlabel('T bias [K]')
    #
    # ###-------------------------
    # plt.subplot(122)
    # ax1 = plt.gca()
    # plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkorange', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'steelblue', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # plt.legend()
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    # plt.grid('on')
    #
    # ###-------------------------
    # fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_ifs_casim-100_wUMGlobal.svg'
    # # plt.savefig(fileout)
    # plt.show()

    #
    # ##################################################
    # ##################################################
    # #### PCOLOR TIMESERIES
    # ##################################################
    # ##################################################
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    #
    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(12,11))
    #
    # Tmin = -45
    # Tmax = 5
    # qbiasmin = -1.0
    # qbiasmax = 1.0
    # ymax = 9000
    # yticklist = [0, 3e3, 6e3, 9e3]
    #
    # #########-------------------------------------------------------------------------------------------
    # ####       DEFINE PERIODS
    # ####               all model data share a timestamp
    # p3 = np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] <= 230.0))
    # p4 = np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] <= 240.0))
    # p5 = np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] <= 247.0))
    # p6 = np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] <= 251.0))
    # p7 = np.where(np.logical_and(data1['time_hrly'][::6] >= 251.0, data1['time_hrly'][::6] <= 255.5))
    # p8 = np.where(data1['time_hrly'][::6] >= 255.5)
    # #   p3: np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] < 230.0))
    # #   p4: np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] < 240.0))
    # #   p5: np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] < 247.0))
    # #   p6: np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] < 251.0))
    #
    # #########-------------------------------------------------------------------------------------------
    #
    # ### -------------------------------
    # ### model temp anomalies wrt radiosondes
    # ### ------------------------------
    # ax  = fig.add_axes([0.05,0.81,0.4,0.12])   # left, bottom, width, height
    # sfig1 = plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
    #     vmin = Tmin, vmax = Tmax)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][0]], data1['time_6hrly'][p3[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p3[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p3[0][-1]], 8000, 'w.')
    # plt.annotate('P3', xy=(227,6800), xytext=(227.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][0]], data1['time_6hrly'][p4[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p4[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p4[0][-1]], 8000, 'w.')
    # plt.annotate('P4', xy=(234,6800), xytext=(234.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][0]], data1['time_6hrly'][p5[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p5[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p5[0][-1]], 8000, 'w.')
    # plt.annotate('P5', xy=(242.5,6800), xytext=(242.6,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][0]], data1['time_6hrly'][p6[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p6[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p6[0][-1]], 8000, 'w.')
    # plt.annotate('P6', xy=(248,6800), xytext=(248.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][0]], data1['time_6hrly'][p7[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p7[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p7[0][-1]], 8000, 'w.')
    # plt.annotate('P7', xy=(252.25,6800), xytext=(252.5,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p8
    # plt.plot([data1['time_6hrly'][p8[0][0]], data1['time_6hrly'][p8[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p8[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p8[0][-1]], 8000, 'w.')
    # plt.annotate('P8', xy=(256,6800), xytext=(256.1,6801), fontsize = 12, color = 'w')
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # # plt.colorbar()
    # plt.ylabel('Z [km]')
    # plt.xlabel('Date')
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.title('T [$^{\circ}$C]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel('Radiosondes', rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    # cbaxes1 = fig.add_axes([0.1, 0.98, 0.3, 0.015])
    # cb1 = plt.colorbar(sfig1, cax = cbaxes1, orientation = 'horizontal')
    #
    # ax  = fig.add_axes([0.05,0.56,0.4,0.12])   # left, bottom, width, height
    # # dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # # data3['temp_anomalies'] = dat3
    # sfig2 = plt.pcolor(data3['time_6hrly'], data1['universal_height'], data3['temp_anomalies'],
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # # plt.colorbar()
    # # plt.set_cmap('seismic')
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label3, rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    # plt.title('T [K]')
    # # plt.title(label3 + ' - Radiosondes, T[degC]')
    # cbaxes2 = fig.add_axes([0.1, 0.73, 0.3, 0.015])
    # cb2 = plt.colorbar(sfig2, cax = cbaxes2, orientation = 'horizontal')
    #
    # ax  = fig.add_axes([0.05,0.39,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'],data1['universal_height'], data2['temp_anomalies'],
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # # plt.colorbar()
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label2[:8], rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    #
    # ax  = fig.add_axes([0.05,0.22,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data4['time_6hrly'],data1['universal_height'], data4['temp_anomalies'],
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # # plt.colorbar()
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label4[:-4], rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    #
    # ax  = fig.add_axes([0.05,0.05,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'],data1['universal_height'], data1['temp_anomalies'],
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.xlabel('Date')
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label1, rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    #
    # ###-------------------------------------
    # ### -------------------------------
    # ### model q anomalies wrt radiosondes
    # ### ------------------------------
    # ax  = fig.add_axes([0.55,0.81,0.4,0.12])   # left, bottom, width, height
    # sfig3 = plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_driftSondes_UM']),
    #     vmin = 0, vmax = qmax)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][0]], data1['time_6hrly'][p3[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p3[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p3[0][-1]], 8000, 'w.')
    # plt.annotate('P3', xy=(227,6800), xytext=(227.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][0]], data1['time_6hrly'][p4[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p4[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p4[0][-1]], 8000, 'w.')
    # plt.annotate('P4', xy=(234,6800), xytext=(234.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][0]], data1['time_6hrly'][p5[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p5[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p5[0][-1]], 8000, 'w.')
    # plt.annotate('P5', xy=(242.5,6800), xytext=(242.6,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][0]], data1['time_6hrly'][p6[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p6[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p6[0][-1]], 8000, 'w.')
    # plt.annotate('P6', xy=(248,6800), xytext=(248.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][0]], data1['time_6hrly'][p7[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p7[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p7[0][-1]], 8000, 'w.')
    # plt.annotate('P7', xy=(252.25,6800), xytext=(252.5,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ## p8
    # plt.plot([data1['time_6hrly'][p8[0][0]], data1['time_6hrly'][p8[0][-1]]], [8000, 8000], 'w--')
    # plt.plot(data1['time_6hrly'][p8[0][0]], 8000, 'w.')
    # plt.plot(data1['time_6hrly'][p8[0][-1]], 8000, 'w.')
    # plt.annotate('P8', xy=(256,6800), xytext=(256.1,6801), fontsize = 12, color = 'w')
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # # plt.colorbar()
    # plt.ylabel('Z [km]')
    # plt.xlabel('Date')
    # plt.title('q [g kg$^{-1}$]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel('Radiosondes', rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    # cbaxes3 = fig.add_axes([0.6, 0.98, 0.3, 0.015])
    # cb3 = plt.colorbar(sfig3, cax = cbaxes3, orientation = 'horizontal')
    #
    # ax  = fig.add_axes([0.55,0.56,0.4,0.12])   # left, bottom, width, height
    # sfig4 = plt.pcolor(data3['time_6hrly'], data1['universal_height'], data3['q_anomalies'],
    #     vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label3, rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    # plt.title('q [g kg$^{-1}$]')
    # cbaxes4 = fig.add_axes([0.6, 0.73, 0.3, 0.015])
    # cb4 = plt.colorbar(sfig4, cax = cbaxes4, orientation = 'horizontal')
    #
    # ax  = fig.add_axes([0.55,0.39,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'], data1['universal_height'], data2['q_anomalies'],
    #     vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label2[:8], rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    #
    # ax  = fig.add_axes([0.55,0.22,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data4['time_6hrly'], data1['universal_height'], data4['q_anomalies'],
    #     vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.ylabel('Z [km]')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label4[:-4], rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    #
    # ax  = fig.add_axes([0.55,0.05,0.4,0.12])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'], data1['universal_height'], data1['q_anomalies'],
    #     vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # #### plot periods
    # ## p3
    # plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ## p8
    # # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    # ##
    # plt.ylim([0,ymax])
    # plt.yticks(yticklist)
    # ax.set_yticklabels([0,3,6,9])
    # plt.xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # # plt.colorbar()
    # plt.ylabel('Z [km]')
    # plt.xlabel('Date')
    # ax2 = ax.twinx()
    # ax2.set_ylabel(label1, rotation = 270, labelpad = 15)
    # ax2.set_yticks([])
    # # plt.title(label1 + ' - Radiosondes, T[degC]')
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # # fileout = '../FIGS/comparisons/TimeSeriesProfiles_TandSpHum-BrBG_ifs_casim-100-GA6alb_ra2t_ra2m-25_Dates_fixedRA2T.png'
    # fileout = '../FIGS/PaperSubmission/Figure8.svg'
    # # plt.savefig(fileout, dpi=300)
    # # plt.show()
    # plt.close()

    return nc, data

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### close any leftover figure instances
    plt.close()

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/INITIAL_CONDITIONS_TEST/'
        init_dir = '/gws/nopw/j04/arcticcloud/MOCCHA/UM_STARTFILES/'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '/home/gillian/MOCCHA/MOCCHA_GIT/UM/DATA/INITIAL_CONDITIONS_TEST/'
        ship_filename = '/home/gillian/MOCCHA/MOCCHA_GIT/ODEN/DATA/2018_shipposition_1hour.txt'
        obs_root_dir = '/home/gillian/MOCCHA/MOCCHA_GIT/ODEN/DATA/'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'

    ### CHOSEN RUN
    out_dir = '23_u-cc278_RA1M_CASIM/'
    out_dir2 = '24_u-cc324_RA2T_CON/'
    out_dir3 = '25_u-cc568_RA2M_CON/'
    out_dir_glm = '24_u-cc324_RA2T_CON/'

    ### Define model of interest (for pulling track only)
    model_list = ['lam','glm']
    model_flag = 0 # 0 for LAM, 1 for GLM
    model = model_list[model_flag]

    ## 4_u-bg610_RA2M_CON/              # Wilson and Ballard 1999 uphys
    ## 5_u-bl661_RA1M_CASIM/            # 100/cc accum mode aerosol; ARG + Cooper
    ## 6_u-bm410_RA1M_CASIM/            # 200/cc accum mode aerosol
    ## 7_u-bn068_RA2T_CON/              # RA2T_CON nest + global 4D stash
    ## 8_u-bp738_RA2M_CON/              # ERAI
    ## 10_u-bq791_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Fletcher
    ## 11_u-bq798_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Meyers
    ## 12_u-br210_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507
    ## 13_u-br409_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; passive aerosol processing
    ## 14_u-bu570_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit
    ## 15_u-bu687_RA2M_CON/           # Wilson and Ballard 1999 uphys; new RHcrit
    ## 16_u-bv926_RA2T_CON/              # RA2T_CON nest + global 4D stash + no subgrid mp production
    ## 17_u-bz429_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics
    ## 18_u-ca011_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics; 1A BL Scheme
    ## 19_u-ca012_RA2T_CON/              # RA2T_CON nest + global 4D stash; includes diagnosed turbulent dissipation rate
    ## 20_u-ca362_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; CICE sea ice albedo scheme
    ## 23_u-cc278_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM
    ## 24_u-cc324_RA2T_CON/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included
    ## 25_u-cc568_RA2M_CON/             # Wilson and Ballard 1999 uphys. sea ice albedo and extra BL diags
    ## 26_u-cd847_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. identical suite = u-cd852
    ## 27_u-ce112_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. passive aerosol processing.
    ## 28_u-ce627_RA2T_CON/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included. Mid-level convection switched off in GLM.
    ## 30_u-cg179_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; passive aerosol processing
    ## 31_u-cl349_RA2M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; jules fluxes

    # -------------------------------------------------------------
    # Extract from MASS with:
    # -------------------------------------------------------------
    #### moo select stash_extract_CloudNet.p2 moose:crum/u-bg610/apm.pp 4_u-bg610_RA2M_CON/20180811T1200Z/
    #### moo select stash_extract_CASIM.p5 moose:crum/u-bl661/apm.pp 5_u-bl661_RA1M_CASIM/20180831T1200Z/
    #### moo select stash_extract_CASIM.p5 moose:crum/u-bm410/apm.pp 6_u-bm410_RA1M_CASIM/20180905T1200Z/
    #### moo select stash_extract_all.p4 moose:crum/u-bn068/apm.pp 7_u-bn068_RA2T_CON/20180827T1200Z/
    #### moo select stash_extract_all.p4 moose:crum/u-bp738/apm.pp 8_u-bp738_RA2M_CON/20180904T1200Z/
    #### moo select stash_extract_CASIM_BL.p6 moose:crum/u-br210/apm.pp 12_u-br210_RA1M_CASIM/20180901T1200Z/
    #### moo select stash_extract_CASIM_BL.p6 moose:crum/u-br409/apm.pp 13_u-br409_RA1M_CASIM/20180901T1200Z/

    #### stash CASIM:
    ####    2, 3, 4, 10, 12, 24, 75, 78, 83, 86, 150, 254, 266, 267, 268, 271, 272, 273, 408, 409, 1201,
    ####    2201, 2391, 2392, 3025, 3217, 3234, 3236, 3245, 3247, 3248, 3360, 3361, 3476, 4118, 4203, 4204,
    ####    5216, 9203, 9204, 9205, 16004, 30461

    #### stash CASIM_BL:
    #### 2, 3, 4, 10, 12, 24, 26, 31, 75, 78, 83, 84, 88, 150, 254, 266, 267, 268, 271, 272, 273, 408, 409,
    #### 1201, 2201, 2391, 2392, 3002, 3025, 3208, 3217, 3219, 3220, 3223, 3234, 3236, 3245, 3247, 3248, 3360,
    #### 3361, 3362, 3363, 3460, 3461, 3464, 3465, 3469, 3471, 3473, 3476, 3501, 4118, 4203, 4204, 5216, 9203,
    #### 9204, 9205, 16004, 30461

    ####    RUN SCRIPT IN BACKGROUND (change to executable with chmod +x diags_CloudNet.py)
    #### module load jaspy
    #### nohup python diags_CloudNet.py > nohup_u-bz429_diags_CloudNet.out &

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data, values = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Combine PP files
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # dummy = combinePP(root_dir, out_dir, date_dir)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load combined PP files - pd for sea ice area fraction cartopy plot
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------


    # date_dir = os.listdir(root_dir + out_dir)
    # cube = loadPD(root_dir, out_dir, date_dir)
    # cube = loadPA(root_dir, out_dir, date_dir, model_flag)
    # print (cube)
    # figure = plot_cartmap(ship_data, cube, date_dir, model)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Pull track (LAM AND GLM Functional)
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # date_dir = os.listdir(root_dir + out_dir)
    # for date in date_dir:
    # # date = '20180815T1200Z'
    #     if date[:4] == '2018':
    #         data = pullTrack(date, root_dir, out_dir, global_con, model, ship_data)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load in relevant UM start dump data
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    # ### list start dumps in UM_STARTFILES/
    # umdumps = os.listdir(init_dir)
    #
    # ### test out with first file
    # startdump = loadUMStartDump(umdumps[0])

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load in observations
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    obs = {}
    obs = loadObservations(obs, platform, obs_root_dir)

    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Load pulled track files
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    dir1 = out_dir + 'OUT_R2/'
    dir2 = out_dir2 + 'OUT_R2_LAM/'
    dir3 = out_dir3 + 'OUT_R2/'
    dir_glm = out_dir_glm + 'OUT_R2_GLM/'

    out_dirs = [dir1, dir2, dir3, dir_glm]

    nc = {}   ### load netcdfs into a single dictionary
    data = {}   ### load postprocessed data into a single dictionary
    for dir in out_dirs:
        model = model_list[model_flag]  ## reset model assignment
        # dir = dir1
        if dir[:2] == '24':
            for model in model_list:
                if model == 'lam':
                    nc[dir[:2]] = {}  ### use run number as dictionary index
                    data[dir[:2]] = {}  ### use run number as dictionary index
                    nc, data = loadNCs(nc, data, root_dir, dir, model)
                elif model == 'glm':
                    nc[dir[:2] + '_glm'] = {}  ### use -glm suffix with dictionary index
                    data[dir[:2] + '_glm'] = {}  ### use -glm suffix with dictionary index
                    nc, data = loadNCs(nc, data, root_dir, dir, model)
        else:
            nc[dir[:2]] = {}
            data[dir[:2]] = {}
            nc, data = loadNCs(nc, data, root_dir, dir, model)

    print (nc.keys())
    print (data.keys())



    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    ### Start sonde analysis
    ### -------------------------------------------------------------------------
    ### -------------------------------------------------------------------------
    filenames = os.listdir(root_dir + dir)
    for dir in out_dirs:
        print (data[dir[:2]].keys())
        nc, data = radiosondeAnalysis(nc, data, dir, obs, filenames, model_list)

    np.save('working_data', data)
    # np.save('working_obs', obs['sondes'])


    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')


if __name__ == '__main__':

    main()
