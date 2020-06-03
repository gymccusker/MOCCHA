"""
Functions for loading specific datasets
==============================

"""

import numpy as np
# from __future__ import print_function

def load_radar(proj, day):

    """
    Function to load in radar data
    Set up to take:
     'proj' argument, options:
        1. moccha
         'day' argument, options:
            1. date in format YYYYMMDD (string)
            2. all (string), loads all data

    """

    from netCDF4 import Dataset
    from time_functions import calcTime_Date2DOY

    if proj == 'moccha':
        if day != 'all':

            ### define filename
            data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
            filename = data_dir + 'mmcr/' + day + '_oden_mira.nc'

            ### choose variables to load
            var_list = ['time','range','Zh']

            ### load file
            nc = Dataset(filename,'r')

            ### make empty data dictionary
            data = {}

            ### populate dictionary
            data['time'] = nc.variables['time'][:]
            data['range'] = nc.variables['range'][:]
            data['Zh'] = nc.variables['Zh'][:]

            ### find date in DOY format
            date = calcTime_Date2DOY(day)

            ### update time array to reference base date
            data['time'] = data['time']/24.0 + date

            nc.close()

        else:

            ### define filename
            data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'

            moccha_names = ['20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

            ### choose variables to load
            var_list = ['time','range','Zh']

            ### make empty data dictionary
            data = {}

            for name in moccha_names:

                filename = data_dir + 'mmcr/' + name[:8] + '_oden_mira.nc'

                ### load file
                nc = Dataset(filename,'r')

                ### populate dictionary
                data['time'] = nc.variables['time'][:]
                data['range'] = nc.variables['range'][:]
                data['Zh'] = nc.variables['Zh'][:]

                ### find date in DOY format
                date = calcTime_Date2DOY(name[:8])

                ### update time array to reference base date
                data['time'] = data['time']/24.0 + date

                nc.close()

    else:
        print('*** invalid project name ***')

    return data
