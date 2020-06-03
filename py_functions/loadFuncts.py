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

    if proj == 'moccha':
        if day != 'all':

            ### define filename
            data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
            file_name = data_dir + 'mmcr/' + day + '_Oden_mira.nc'

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

    else:
        print('*** invalid project name ***')

    return data
