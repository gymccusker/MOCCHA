#####################################
##
## Read in and plot ship track
##      -- future release: relate to UM nest
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
    fig = plt.figure(figsize=(6,8))

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 2500000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=86,lat_0=86,lon_0=0)
    m.drawcoastlines()
    m.bluemarble()

    # define parallels/meridians
    m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,0,0,0],linewidth=0.8,fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,20.),labels=[0,0,0,1],linewidth=0.8,fontsize=10)
    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='lightgrey')
    # m.fillcontinents(color='white')

    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(data)
    inIce_index = inIce(data)

    ### MAP ONTO PROJECTION
    x, y = m(data.values[:,6], data.values[:,7])
    x_inIcePeriod, y_inIcePeriod = m(data.values[inIce_index,6],data.values[inIce_index,7])
    x_driftPeriod, y_driftPeriod = m(data.values[drift_index,6],data.values[drift_index,7])

    # Plot tracks as line plot
    plt.plot(x, y, color = 'yellow', linewidth = 2, label = 'Whole')
    plt.plot(x_inIcePeriod, y_inIcePeriod, color = 'darkorange', linewidth = 3, label = 'In Ice')
    plt.plot(x_driftPeriod, y_driftPeriod, color = 'red', linewidth = 4, label = 'Drift')

    ###########################################
    ### PLOT SWATH FOR INCREASED FREQ DIAGS
    ###########################################
        # I.B.:
        # Drift limits are:
        # latitude   88.4502 to 89.6388
        # longitude  4.6830 to 73.7629

    latr = np.arange(88.4502,89.6388,0.01188600000000008)
    lonr = np.arange(4.6830,73.7629,0.690799)

    lon, lat = np.meshgrid(lonr, latr)

    nx = np.size(lon,0)
    ny = np.size(lat,1)
    x1, y1 = m(lon[ny-1,0],lat[ny-1,0])
    x2, y2 = m(lon[0,0],lat[0,0])
    x3, y3 = m(lon[0,nx-1],lat[0,nx-1])
    x4, y4 = m(lon[ny-1,nx-1],lat[ny-1,nx-1])

    # draw nests
    pol =  Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],\
                  facecolor='none',linestyle='--',edgecolor='w',linewidth=2,label='Swath')
    plt.gca().add_patch(pol)

    ### ADD LEGEND
    plt.legend()

    plt.title('MOCCHA ship track')

    plt.show()

def main():

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ####    UM DIRECTORY ON JASMIN
    # root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
    # out_dir = '2_20180801_61DIAGS_TEST/'


    #### FILENAME - DEPENDS ON PLATFORM
    ## LEEDS
    # filename = '/nfs/see-fs-02_users/eargy/MOCCHA/gillian/ship/2018_shipposition_1hour.txt'

    ## JASMIN
    filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'

    #### CODE NOT WORKING
    # seaice_rootdir = '/nfs/see-fs-02_users/eargy/MOCCHA/parent/data/seaice/AMSR2/'
    # seaice_file = 'asi-AMSR2-n6250-20180801-v5.hdf'
    # seaice_path = seaice_rootdir + seaice_file

    data, values = readfile(filename)

    columns = assignColumns(data)

    # print ''
    # print data.head
    # print ''

    map = plotmap(data)

    # seaice = readSeaice(seaice_file)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
