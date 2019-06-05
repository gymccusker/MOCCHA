##--------------------------------------------------------------------------
##
##			Script to read in WRF output files, extract necessary data,
##			then save into new NetCDF file (reduces file size for archiving)
##					-- GYoung
##
##--------------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np
import time
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
import matplotlib.pyplot as plt


def readFile(filename):

	###################################
	# Load in data
	###################################

	dat = np.loadtxt(filename)

	###################################
	## Define headers
	###################################

	headers = assignColumns(dat)

	###################################
	## Read in ASCII variables
	###################################

	# ### Create dictionary to contain each variable
	data = {}

	for i in range(0,len(headers)-1):
		data = data.update({headers[i]: dat[:,i]})

	dat.close()

	return data

def assignColumns(data):

    columns = ['idat','itim','vdat','vtim#','lev','P_HIST','U_HIST','V_HIST','T_HIST','Q_HIST','L_HIST','I_HIST','A_HIST','R_HIST','W_HIST','N_HIST','O_HIST']

    return columns

def writeNetCDF(data):

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
    desc = 'ECMWF test netCDF write out'
    micro = ''
    dataset.description = desc
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'N/A'
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic.'
    dataset.comment = ''
    dataset.institution = 'ECMWF/University of Leeds.'

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

	##--------------------------------------------------------------------------
	##--------------------------------------------------------------------------
	##---------------				IN
	##--------------------------------------------------------------------------
	##--------------------------------------------------------------------------


    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'DESKTOP'

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
        root_dir = '/nfs/see-fs-02_users/eargy/MOCCHA/parent/working/data/ecmwf/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'


	###################################
	###################################
	## SET FILENAME TO READ IN
	###################################
	###################################

	filename = root_dir + 'MOCCHA_001_20180812_var'

	data = readFile(filename)



    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''


if __name__ == '__main__':

    main()
