#####################################
##
## Read in and plot ship track
##      -- future release: relate to UM nest
##
#####################################

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime, date
import time

def readfile(filename):

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data, values

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def plotmap(data):

    # import cartopy
    # import cartopy.crs as crs
    # import cartopy.feature as cfe
    from mpl_toolkits.basemap import Basemap

    ###################################
    ## PLOT MAP
    ###################################

    ##################################################
    ##################################################
    #### 	BASEMAP
    ##################################################
    ##################################################

    SMALL_SIZE = 10
    MED_SIZE = 12
    LARGE_SIZE = 14

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=SMALL_SIZE)
    plt.rc('ytick',labelsize=SMALL_SIZE)
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.rc('figure',titlesize=LARGE_SIZE)

    ## create figure and axes instances
    fig = plt.figure(figsize=(5,4))

    #########################################################################################################

    ax  = fig.add_axes([0.1,0.1,0.8,0.8])	# left, bottom, width, height

    ### MAP DIMENSIONS
    dim = 6500000

    m = Basemap(width=0.75*dim,height=dim,
                resolution='l',projection='stere',\
                lat_ts=90,lat_0=90,lon_0=0)
    m.drawcoastlines()
    m.bluemarble()

    # define parallels/meridians
    m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,0,0,0],linewidth=0.8,fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,5.),labels=[0,0,1,0],linewidth=0.8,fontsize=10)
    m.drawcoastlines(linewidth=1.)

    # m.drawmapboundary(fill_color='lightgrey')
    # m.fillcontinents(color='white')

    ######################################################
    ### CARTOPY CODE
    ######################################################
    # # Create a figure
    # fig = plt.figure(figsize=(8,4))
    #
    # # Set the GeoAxes to the projection used by WRF
    # ax = fig.add_axes([0.1,0.1,0.4,0.8], projection=cart_proj)	# left, bottom, width, height
    # # ax = plt.axes(projection=cart_proj)

    # Add coastlines
    # ax.coastlines('50m', linewidth=0.8)

    # Plot ship track as line plot
    x, y = m(data.values[:,6], data.values[:,7])
    plt.plot(x, y, linewidth = 2)

    # Add a color bar
    # cbar = plt.colorbar(ax=ax, shrink=.62)
    # cbar.set_label(qncloud1.name[-5:])

    # Set the map limits.  Not really necessary, but used for demonstration.
    # ax.set_xlim(wrf.cartopy_xlim(qncloud1))
    # ax.set_ylim(wrf.cartopy_ylim(qncloud1))

    # Add the gridlines
    # ax.gridlines(color="black", linestyle="dotted")

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

    filename = '/nfs/see-fs-02_users/eargy/MOCCHA/gillian/ship/2018_shipposition_1hour.txt'

    data, values = readfile(filename)

    columns = assignColumns(data)

    print ''
    print data.head
    print ''

    map = plotmap(data)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
