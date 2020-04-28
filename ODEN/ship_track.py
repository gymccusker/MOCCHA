#####################################
##
## Read in and plot ship track, UM nest + swath
##
##
#####################################

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime, date
import time

def readfile(filename):

    import pandas as pd

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data, values

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def readSeaice(seaice_path):

    import h5py as hdf

    seaice = hdf.File(seaice_path, 'r')

    ### LIST HDF KEYS
    seaice_keys = seaice.keys()

    print ''
    print seaice_keys
    print ''

    return seaice

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = np.arange(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print ''
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print 'Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4])
    print ''

    return drift_index

def periodHighlight(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_highlight = np.where(np.logical_and(data.values[:,2]>=28,data.values[:,1]==8))
    Sep_highlight = np.where(np.logical_and(np.logical_and(data.values[:,2]<=2,data.values[:,1]==9),data.values[:,3]<=23))
    highlight = np.arange(Aug_highlight[0][0],Sep_highlight[0][-1])

    print ''
    # print 'Aug drift: ' + str(data.values[Aug_highlight[0][0],0:3]) + ' - ' + str(data.values[Aug_highlight[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_highlight[0][0],0:3]) + ' - ' + str(data.values[Sep_highlight[0][-1],0:3])
    print 'Whole drift: ' + str(data.values[highlight[0],0:4]) + ' - ' + str(data.values[highlight[-1],0:4])
    print ''

    return highlight

def inIce(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    Sep_inIce = np.where(np.logical_and(data.values[:,2]<20,data.values[:,1]==9))
    inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    print ''
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ''

    return inIce_index

def gridSetup(lont, latt, m):

    lon, lat = np.meshgrid(lont, latt)

    nx = np.size(lon,0)
    ny = np.size(lat,1)
    x1, y1 = m(lon[nx-1,0],lat[nx-1,0])
    x2, y2 = m(lon[0,0],lat[0,0])
    x3, y3 = m(lon[0,ny-1],lat[0,ny-1])
    x4, y4 = m(lon[nx-1,ny-1],lat[nx-1,ny-1])

    return x1, x2, x3, x4, y1, y2, y3, y4

def plotmap(data):

    # import cartopy
    # import cartopy.crs as crs
    # import cartopy.feature as cfe
    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

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
    fig = plt.figure(figsize=(7, 10), dpi=300)

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 1750000

    m = Basemap(width=0.7*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=82,lat_0=83,lon_0=0)
    m.drawcoastlines()
    # m.bluemarble()

    # define parallels/meridians
    # m.drawparallels(np.arange(50.,91.,2.),color='k')
    # m.drawmeridians(np.arange(-180.,181.,10.),color='k')

    parallels = np.arange(0.,90,5.)
    # # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,False,False,False], color = 'grey')
    meridians = np.arange(10.,351.,20.)
    m.drawmeridians(meridians,labels=[False,False,False,False], color = 'grey')

    m.drawcoastlines(linewidth=1.)

    m.fillcontinents(color='lightgrey')

    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(data)
    highlight_index = periodHighlight(data)
    inIce_index = inIce(data)

    ### MAP ONTO PROJECTION
    x, y = m(data.values[:,6], data.values[:,7])
    x_inIcePeriod, y_inIcePeriod = m(data.values[inIce_index,6],data.values[inIce_index,7])
    x_driftPeriod, y_driftPeriod = m(data.values[drift_index,6],data.values[drift_index,7])
    x_highlightPeriod, y_highlightPeriod = m(data.values[highlight_index,6],data.values[highlight_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, '--', color = 'pink', linewidth = 2, label = 'Whole')
    plt.plot(x_inIcePeriod, y_inIcePeriod, color = 'palevioletred', linewidth = 4, label = 'In Ice')
    plt.plot(x_driftPeriod, y_driftPeriod, color = 'red', linewidth = 5, label = 'Drift')
    # plt.plot(x_highlightPeriod, y_highlightPeriod, color = 'purple', linewidth = 4, label = 'Case examples')

    ###########################################
    ### PLOT NEST + SWATH FOR INCREASED FREQ DIAGS VIS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629
        # R.P. original (0, 86.625) @ 500x500

    ### SWATH
    # lats = np.arange(80.9998,89.9998,0.09)
    # lons = np.arange(3.0,76.0,0.73)
    # x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)
    #
    # ### NEST (input)
    # grx = float(5600)
    # gry = float(700)
    # centlon = float(39.5)
    # centlat = float(85.275)
    # latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    # lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    # x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    # draw swath
    # pols =  Polygon([(x1s,y1s),(x2s,y2s),(x3s,y3s),(x4s,y4s)],\
    #               facecolor='none',linestyle='--',edgecolor='g',linewidth=2,label='Swath')
    # plt.gca().add_patch(pols)

    # draw nest
    # poln =  Polygon([(x1n,y1n),(x2n,y2n),(x3n,y3n),(x4n,y4n)],\
    #               facecolor='none',linestyle='-',edgecolor='w',linewidth=2,label='Nest')
    # plt.gca().add_patch(poln)

    ### ADD LEGEND
    plt.legend()

    # plt.title('MOCCHA ship track')

    plt.savefig('FIGS/HighArctic_vPRES.svg')
    plt.show()

def plotmap_poster(data):

    # import cartopy
    # import cartopy.crs as crs
    # import cartopy.feature as cfe
    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

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
    fig = plt.figure(figsize=(8.27, 11.69), dpi=300)

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 2500000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=82,lat_0=82,lon_0=-25)
    m.drawcoastlines()
    # m.bluemarble()

    # define parallels/meridians
    # m.drawparallels(np.arange(50.,91.,2.),color='k')
    # m.drawmeridians(np.arange(-180.,181.,10.),color='k')

    parallels = np.arange(0.,90,5.)
    # # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[False,False,False,False], color = 'grey')
    meridians = np.arange(10.,351.,20.)
    m.drawmeridians(meridians,labels=[False,False,False,False], color = 'grey')

    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='lightgrey')
    m.fillcontinents(color='lightgrey')

    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(data)
    highlight_index = periodHighlight(data)
    inIce_index = inIce(data)

    ### MAP ONTO PROJECTION
    x, y = m(data.values[:,6], data.values[:,7])
    x_inIcePeriod, y_inIcePeriod = m(data.values[inIce_index,6],data.values[inIce_index,7])
    x_driftPeriod, y_driftPeriod = m(data.values[drift_index,6],data.values[drift_index,7])
    x_highlightPeriod, y_highlightPeriod = m(data.values[highlight_index,6],data.values[highlight_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, '--', color = 'pink', linewidth = 2, label = 'Whole')
    plt.plot(x_inIcePeriod, y_inIcePeriod, color = 'palevioletred', linewidth = 3, label = 'In Ice')
    plt.plot(x_driftPeriod, y_driftPeriod, color = 'red', linewidth = 3, label = 'Drift')
    # plt.plot(x_highlightPeriod, y_highlightPeriod, color = 'purple', linewidth = 4, label = 'Case examples')

    ###########################################
    ### PLOT NEST + SWATH FOR INCREASED FREQ DIAGS VIS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629
        # R.P. original (0, 86.625) @ 500x500

    ### SWATH
    # lats = np.arange(80.9998,89.9998,0.09)
    # lons = np.arange(3.0,76.0,0.73)
    # x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)
    #
    # ### NEST (input)
    # grx = float(5600)
    # gry = float(700)
    # centlon = float(39.5)
    # centlat = float(85.275)
    # latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    # lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    # x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    # draw swath
    # pols =  Polygon([(x1s,y1s),(x2s,y2s),(x3s,y3s),(x4s,y4s)],\
    #               facecolor='none',linestyle='--',edgecolor='g',linewidth=2,label='Swath')
    # plt.gca().add_patch(pols)

    # draw nest
    # poln =  Polygon([(x1n,y1n),(x2n,y2n),(x3n,y3n),(x4n,y4n)],\
    #               facecolor='none',linestyle='-',edgecolor='w',linewidth=2,label='Nest')
    # plt.gca().add_patch(poln)

    ### ADD LEGEND
    plt.legend()

    # plt.title('MOCCHA ship track')

    plt.savefig('FIGS/HighArctic_vPOSTER_v3.svg',dpi=100)
    plt.show()

def main():

    START_TIME = time.time()
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
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    #### CODE NOT WORKING
    # seaice_rootdir = '/nfs/see-fs-02_users/eargy/MOCCHA/parent/data/seaice/AMSR2/'
    # seaice_file = 'asi-AMSR2-n6250-20180801-v5.hdf'
    # seaice_path = seaice_rootdir + seaice_file

    data, values = readfile(ship_filename)

    columns = assignColumns(data)

    # print ''
    # print data.head
    # print ''

    # map = plotmap(data)
    map = plotmap_poster(data)

    # seaice = readSeaice(seaice_file)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
