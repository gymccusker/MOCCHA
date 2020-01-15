"""
Steps to read in Matlab struct files (saved as .mat) and associated functions to make them easier to use
==============================

"""

from scipy.io import loadmat
import numpy as np
import scipy.io as sio

def readMatlabStruct(filename, struct_name):

    #### EXAMPLE OF USE:
    #### data = readMatlabStruct('../../jutta/UserReadyData/radiosondes/SondeData_h10int_V02.mat','RS10intall')
    #### find struct name with:
        #### sio.whosmat(filename)

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
    struct = dat[struct_name]
    print ''

    #### --------------------------------------------------------------------
    #### IDENTIFY DATA AS FIRST ENTRY IN INTERMEDIATE STRUCT
    #### --------------------------------------------------------------------
    data = struct[0,0]
        #### data.dtype:
            #### returns keys of dictionary (normal python dictionary access
            #### commands don't quite work...). MATLAB structs come back as
            #### numpy structured arrays.

    #### --------------------------------------------------------------------
    #### CHANGE NUMPY STRUCTURED ARRAY TO DICTIONARY FOR EASE OF USE
    ####            --- come back to this later
    #### --------------------------------------------------------------------
    # dicto = {}
    # for name in dat[struct_name].dtype.names:#len(dat[struct_name].dtype.names)):
    #     dicto[dat[struct_name].dtype.names]= dat[struct_name][]

    print 'Finished! :)'
    print 'Reading out ' + struct_name + ' struct within .mat file'
    print ''

    return data      #### returns structured numpy array containing matlab struct
