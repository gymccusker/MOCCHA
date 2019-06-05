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
		data.update({headers[i]: dat[:,i]})

	return data

def splitData(data_all, filename):

	###################################
	# Split up data into columns per hourly dump
	###################################

	data = {}
	temp1 = np.arange(data_all['vtim'][0],24.0)
	temp2 = np.arange(0.0,24.0)
	temp3 = np.arange(0.0,data_all['vtim'][-1]+1)
	time_array = np.append(temp1,temp2)
	time_array = np.append(time_array,temp3)
	time_index = np.where(np.logical_and(data_all['vtim'] == data_all['vtim'][0],data_all['vdat'] == data_all['vdat'][0]))
	data = {'model_level': data_all['lev'][time_index], 'time': time_array}
	pres = np.zeros([np.size(time_index),np.size(time_array)])
	temp = np.zeros([np.size(time_index),np.size(time_array)])
	rh = np.zeros([np.size(time_index),np.size(time_array)])
	for i in range(0,len(time_array)):
		time_index = np.where(np.logical_and(data_all['vtim'] == data_all['vtim'][i*137],data_all['vdat'] == data_all['vdat'][i*137]))
		pres[:,i] = data_all['P_HIST'][time_index]
		temp[:,i] = data_all['T_HIST'][time_index]
		rh[:,i] = data_all['R_HIST'][time_index]
		# data.update({'time': data_all['vtim'][time_index][0],
		# 		'pressure': data_all['P_HIST'][time_index],
		# 		'u_wind': data_all['U_HIST'][time_index],
		# 		'v_wind': data_all['V_HIST'][time_index],
		# 		'w_wind': data_all['W_HIST'][time_index],
		# 		'temperature': data_all['T_HIST'][time_index],
		# 		'relative_humidity': data_all['R_HIST'][time_index],
		# 		'cloud_fraction': data_all['A_HIST'][time_index],
		# 		'ice_mixing_ratio': data_all['I_HIST'][time_index],
		# 		'liquid_mixing_ratio': data_all['L_HIST'][time_index],
		# 		'vapour_mixing_ratio': data_all['Q_HIST'][time_index],
		# 		})
		print data['time']

	data.update{'pressure': pres, 'temperature': temp, 'relative_humidity': rh}

	return data

def assignColumns(data):

    columns = ['idat','itim','vdat','vtim','lev','P_HIST','U_HIST','V_HIST','T_HIST','Q_HIST','L_HIST','I_HIST','A_HIST','R_HIST','W_HIST','N_HIST','O_HIST']

    return columns

def writeNetCDF(data, outfile, filename):

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
    dataset.source = ''
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
    time = dataset.createDimension('time', np.size(data['time']))
    mlevel = dataset.createDimension('model_level', np.size(data['model_level']))
	# plevel = dataset.createDimension('pressure_level', np.size(data['P_HIST']))

    ###################################
    ## Dimensions variables
    ###################################
    #### Time
    time = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999')
    time.comment = 'Hour after '
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
        root_dir = '/nfs/see-fs-02_users/eargy/MOCCHA/parent/working/data/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'

	###################################
	###################################
	## SET FILENAME TO READ IN
	###################################
	###################################

	ascii_filename = 'ecmwf/MOCCHA_001_20180812_var'
	filename = root_dir + ascii_filename

	data_all = readFile(filename)
	print data_all
	data = splitData(data_all, filename)

	nc_filename = root_dir + 'ecmwf_netcdf/' + ascii_filename[17:25] + '_oden_ecmwf' + '.nc'
	# out = writeNetCDF(data, nc_filename)

    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''


if __name__ == '__main__':

    main()
