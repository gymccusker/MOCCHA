"""
Steps to read in Matlab struct files (saved as .mat) and associated functions to make them easier to use
==============================

"""

from scipy.io import loadmat
import numpy as np
import scipy.io as sio

def readMatlabStruct(filename):

    #### EXAMPLE OF USE:
    #### data = readMatlabStruct('../../jutta/UserReadyData/radiosondes/SondeData_h10int_V02.mat','RS10intall')
    #### find struct name with:
        #### sio.whosmat(filename)
    #### access data with:
        #### e.g. data['mday']
        ####        data.dtype.names - lists var names in structured array

    ### ----------------------------------
    ### Find struct name from .mat file using sio
    ### ----------------------------------
    dat = sio.whosmat(filename)

    ### ----------------------------------
    ### Extract out struct name
    ### ----------------------------------
    structname = dat[0][0]

    #### --------------------------------------------------------------------
    #### LOAD MATLAB FILE USING SCIPY
    #### --------------------------------------------------------------------
    print 'Reading in .mat file including struct...'
    dat = loadmat(filename)
    print ''

    #### --------------------------------------------------------------------
    #### USE STRUCT_NAME TO DEFINE INTERMEDIATE STRUCT ARRAY
    #### --------------------------------------------------------------------
    print 'Dealing with intermediate data assignments...'
    struct = dat[structname]
    print ''

    #### --------------------------------------------------------------------
    #### IDENTIFY DATA AS FIRST ENTRY IN INTERMEDIATE STRUCT
    #### --------------------------------------------------------------------
    a = struct[0,0]
        #### data.dtype:
            #### returns keys of dictionary (normal python dictionary access
            #### commands don't quite work...). MATLAB structs come back as
            #### numpy structured arrays.

    #### --------------------------------------------------------------------
    #### CHANGE NUMPY STRUCTURED ARRAY TO DICTIONARY FOR EASE OF USE
    ####
    #### --------------------------------------------------------------------
    # aa = a.astype(float)
    b = {name:a[name].astype(float) for name in a.dtype.names}

    print 'Finished! :)'
    print 'Reading out ' + structname + ' struct within .mat file'
    print ''

    return b     #### returns structured numpy array containing matlab struct

# def findMatlabStruct(filename):
#
#     '''
#     Use to find Matlab struct name in .mat file (can never remember the command!)
#     '''
#
#     ### ----------------------------------
#     ### Find struct name from .mat file using sio
#     ### ----------------------------------
#     dat = sio.whosmat(filename)
#
#     ### ----------------------------------
#     ### Extract out struct name
#     ### ----------------------------------
#
#     structname = dat[0][0]
#
#     return structname
