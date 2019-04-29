###
###
### SCRIPT TO READ IN UM MODEL DATA AS IRIS CUBE
###
###

import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm

def cart_plot(data1, data2):

    import iris.plot as iplt
    import cartopy.crs as crs
    import cartopy.feature as cfe
    import iris.quickplot as qplt

    # Create a figure
    fig = plt.figure(figsize=(8,4))

    ## set axes position
    ax = fig.add_axes([0.1,0.1,0.35,0.8])	# left, bottom, width, height

    ## Allow interactive plotting
    # plt.interactive(True)

    ## Draw the contour with 25 levels.
    # contour = qplt.contourf(data[0,0,:,:], cmap = mpl_cm.Reds)
    # iplt.plot(data[0,:,0,0])
    # iplt.plot(data[0,:,200,200],data.aux_coords[2])
    ## iplt.plot(data[1,0:30,200,200],data.aux_coords[2][0:30])
    for i in range(0,np.size(data1,0)):
        strgi = "%1.f" % (i) # string of timestep
        iplt.plot(data1[i,0:30,200,200],data1.aux_coords[2][0:30],label=strgi)
    plt.legend()
    plt.title(data1.standard_name + ', ' + str(data1.units))

    # iplt.plot(data1[1,0:30,200,200],data1.aux_coords[2][0:30])

    ## set axes position
    ax = fig.add_axes([0.5,0.1,0.35,0.8])	# left, bottom, width, height
    for i in range(0,np.size(data1,0)):
        strgi = "%1.f" % (i) # string of timestep
        iplt.plot(data2[i,0:30,200,200],data2.aux_coords[2][0:30],label=strgi)
    plt.legend()
    plt.title(data2.standard_name + ', ' + str(data2.units))

    plt.show()

def plot_basemap(ship_data, cube):

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plot basemap:'
    print ''

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
    print '******'
    print ''
    print 'lat/lon vertices of proposed swath: ', lats[0], lats[-1], lons[0], lons[-1]
    print ''
    x1s, x2s, x3s, x4s, y1s, y2s, y3s, y4s = gridSetup(lons, lats, m)

    ### NEST (input)
    grx = float(500)
    gry = float(500)
    centlon = float(0.0)
    centlat = float(86.625)
    latn = np.arange((centlat-(gry*float(0.5)*0.0135)),(centlat+(gry*float(0.5)*0.0135)),0.0135)
    lonn = np.arange((centlon-(grx*float(0.5)*0.0135)),(centlon+(grx*float(0.5)*0.0135)),0.0135)
    print '******'
    print ''
    print 'lat/lon vertices of nest (input): ', latn[0], latn[-1], lonn[0], lonn[-1]
    print ''
    x1n, x2n, x3n, x4n, y1n, y2n, y3n, y4n = gridSetup(lonn, latn, m)

    ### NEST (output)
    lono, lato = unrotateGrid(cube)
    print '******'
    print ''
    print 'lat/lon vertices of nest (output): ', lato[0], lato[-1], lono[0], lono[-1]
    print ''
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

def readfile(filename):

    import pandas as pd

    # print '******'
    print ''
    print 'Reading .txt file with pandas'
    print ''

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
    drift_index = np.arange(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print '******'
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

    print '******'
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
    print '******'
    print 'Test of unrotated coordinate grid: '
    print 'Rotated lon coord = ', rot_lon[0]
    print 'Rotated lat coord = ', rot_lat[0]
    print 'Lon = ', lon[0]
    print 'Lat = ', lat[0]
    print ' '

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
    print '******'
    print 'Test of rotated coordinate grid: '
    print 'Lon = ', lon[0]
    print 'Lat = ', lat[0]
    print 'Rotated lon coord = ', rot_lon[0]
    print 'Rotated lat coord = ', rot_lat[0]
    print ' '

    # ******
    # lat/lon vertices of nest (output):  85.960219407715 80.41973098346767 49.567255645848604 -27.55740381723681
    # ******

    return lon, lat

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    print GlobalStashList
    print GlobalStashList[0]

    return GlobalStashList

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def writeNetCDF(cube, outfile):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print '******'
    print ''
    print 'Writing NetCDF file:'
    print ''

    ###################################
    ## Open File
    ###################################
    dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC')

    print ''
    print dataset.file_format
    print ''

    ###################################
    ## Global Attributes
    ###################################
    desc = 'Test netCDF write out'
    micro = 'Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983)'
    dataset.description = desc
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'N/A'
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    dataset.comment = ''
    dataset.institution = 'University of Leeds.'

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    ###################################
    time = dataset.createDimension('time', np.size(cube.coord('time').points))
    # Z = dataset.createDimension('Z', np.size(icenum1,1))
    lat = dataset.createDimension('grid_latitude', np.size(cube.coord('grid_latitude').points))
    lon = dataset.createDimension('grid_longitude', np.size(cube.coord('grid_longitude').points))

    ###################################
    ## Dimensions variables
    ###################################
    #### Time
    time = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999')
    time.comment = cube.coord('time').standard_name
    time.units = str(cube.coord('time').units)
    time[:] = cube.aux_coords[1].points

    #### Latitude
    lat = dataset.createVariable('grid_latitude', np.float32, ('grid_latitude',),fill_value='-9999')
    lat.comment = cube.coord('grid_latitude').standard_name
    lat.units = str(cube.coord('grid_latitude').units)
    lat[:] = cube.coord('grid_latitude').points

    #### Longitude
    lon = dataset.createVariable('grid_longitude', np.float32, ('grid_longitude',),fill_value='-9999')
    lon.comment = cube.coord('grid_latitude').standard_name
    lon.units = str(cube.coord('grid_longitude').units)
    lon[:] = cube.coord('grid_longitude').points

    ###################################
    ###################################
    ## Create 2-d variables
    ###################################
    ###################################
    #### Nisg
    # nisgbar = dataset.createVariable('nisgbar', np.float32, ('time','Z',),fill_value='-9999')
    # nisgbar.long_name = 'domain-averaged total ice number concentration'
    # nisgbar.comment = 'Sum of ice, snow, and graupel particles'
    # nisgbar.units = 'L-1'
    # nisgbar[:] = icenum1[:]

    ###################################
    ###################################
    ## Create 3-d variables
    ###################################
    ###################################
    data = dataset.createVariable(cube.standard_name, np.float32, ('time','grid_latitude','grid_longitude',),fill_value='-9999')
    data.units = str(cube.units)
    # data.comment = cube.metadata
    data.attributes = str(cube.attributes)
    data[:] = cube.data[:]

    ###################################
    ###################################
    ## Create 4-d variables
    ###################################
    ###################################
    # nisg = dataset.createVariable('nisg', np.float32, ('time','Z','X','Y',),fill_value='-9999')
    # nisg.long_name = 'total ice number concentration'
    # nisg.comment = 'Sum of ice, snow, and graupel particles'
    # nisg.units = 'L-1'
    # nisg[:] = ice1[:]

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'JASMIN'
    ### LAPTOP
    ### MONSOON

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'

    ### CHOSEN RUN
    out_dir = '2_20180801_61DIAGS_TEST/'

    ## 1_20160401_61DIAG_TEST/
    ## 2_20180801_61DIAGS_TEST

    #### 	FILE NAMES
    #In [15]: ls
    #umnsaa_cb000  umnsaa_cb009  umnsaa_cb018  umnsaa_cb027     umnsaa_pverb000  umnsaa_pverc024
    #umnsaa_cb001  umnsaa_cb010  umnsaa_cb019  umnsaa_cb028     umnsaa_pverb006  umnsaa_pverd000
    #umnsaa_cb002  umnsaa_cb011  umnsaa_cb020  umnsaa_pa000     umnsaa_pverb012  umnsaa_pverd006
    #umnsaa_cb003  umnsaa_cb012  umnsaa_cb021  umnsaa_pb000     umnsaa_pverb018  umnsaa_pverd012
    #umnsaa_cb004  umnsaa_cb013  umnsaa_cb022  umnsaa_pvera000  umnsaa_pverb024  umnsaa_pverd018
    #umnsaa_cb005  umnsaa_cb014  umnsaa_cb023  umnsaa_pvera006  umnsaa_pverc000  umnsaa_pverd024
    #umnsaa_cb006  umnsaa_cb015  umnsaa_cb024  umnsaa_pvera012  umnsaa_pverc006  umnsa.stash
    #umnsaa_cb007  umnsaa_cb016  umnsaa_cb025  umnsaa_pvera018  umnsaa_pverc012  umnsa.xhist
    #umnsaa_cb008  umnsaa_cb017  umnsaa_cb026  umnsaa_pvera024  umnsaa_pverc018

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print '******'
    print ''
    print 'Make stash list for cube read in at ' + time.strftime("%c")
    print ' '
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load

    print '******'
    print ''
    print 'Identifying .pp files: '
    filename1 = root_dir + out_dir + 'umnsaa_pb000'
    for i in range(0,3):
        res = i*3.0
        str_i = "%03d" % res # file number
        fileout = root_dir + out_dir + 'umnsaa_pc' + str_i
        print fileout
        print ' '

        # # -------------------------------------------------------------
        # # Load cubes
        # # -------------------------------------------------------------
        # print '******'
        # print 'Begin cube read in at ' + time.strftime("%c")
        # print ' '
        #
        # # cube = iris.load(filenames, global_con, callback)
        # # cube = iris.load(filename1, global_con, callback)
        #
        # ## Set variable constraint (i.e. which variable to load in based on stash code)
        # var_con = iris.AttributeConstraint(STASH='m01s16i222')
        # cube = iris.load_cube(fileout, var_con)


    # -------------------------------------------------------------
    # Load cubes
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin cube read in at ' + time.strftime("%c")
    print ' '

    # cube = iris.load(filenames, global_con, callback)
    # cube = iris.load(filename1, global_con, callback)

    ## Set variable constraint (i.e. which variable to load in based on stash code)
    var_con = iris.AttributeConstraint(STASH='m01s16i222')
    cube = iris.load_cube(filename1, var_con)

    print '******'
    print ''
    print 'Cubes read in complete at ' + time.strftime("%c")
    print ' '

    # cube = iris.load(filename1)
    print '******'
    print ''
    print cube # lists all diagnostics in file
    print ''

    # rot_pole = cube1.coord('grid_latitude').coord_system.as_cartopy_crs()

    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data, values = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # map = plot_basemap(ship_data, cube1)

    nc_filename = filename1 + '_r0.nc'
    # out = writeNetCDF(cube, nc_filename)

    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
