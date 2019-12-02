"""
Steps to read in Matlab struct files (saved as .mat)
==============================

"""

from scipy.io import loadmat
import numpy as np


def readMAT(filename, struct_name):

    #### --------------------------------------------------------------------
    #### LOAD MATLAB FILE USING SCIPY
    #### --------------------------------------------------------------------
    dat = loadmat(filename)

    #### --------------------------------------------------------------------
    #### USE STRUCT_NAME TO DEFINE INTERMEDIATE STRUCT ARRAY
    #### --------------------------------------------------------------------
    struct = dat[struct_name]

    #### --------------------------------------------------------------------
    #### IDENTIFY DATA AS FIRST ENTRY IN INTERMEDIATE STRUCT
    #### --------------------------------------------------------------------
    data = struct[0,0]

    return data      #### returns dictionary containing matlab struct
