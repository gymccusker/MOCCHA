"""
Functions for loading specific datasets
==============================

"""

import numpy as np
from __future__ import print_function

def load_radar(day, proj):

    """
    Function to load in radar data
    Set up to take 'proj' argument, options:
        1. moccha
    """

    if proj == 'moccha':
        data_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        file_name = data_dir + 'mmcr/' day + '_Oden_mira.nc'];

    else:
        print('*** invalid project name ***')

    return
