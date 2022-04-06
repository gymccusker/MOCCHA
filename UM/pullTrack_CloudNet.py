###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE,
###         PULL SHIP TRACK, AND OUTPUT AS NETCDF FOR CLOUDNET
###

from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import numpy as np
import diags_MOCCHA as diags
import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')

def readfile(filename):

    import pandas as pd

    # print '******'
    print ('')
    print ('Reading .txt file with pandas')
    print ('')

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def unrotateGrid(cube):
    ##
    ##      ROTATE GRID BACK TO STANDARD
    ##

    import iris.analysis.cartography as ircrt
    import pandas as pd

    ### LAM Configuration from suite u-bg610
    dlat = 0.015714
    dlon = 0.016334
    frst_lat = -5.6112
    frst_lon = 353.0345
    pole_lat = 3.375 # UM SUITE: 37.5000
    pole_lon = 210.0 # UM SUITE: 177.5000

    rot_lat = cube.coord('grid_latitude').points
    rot_lon = cube.coord('grid_longitude').points

    rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()

    lon, lat = ircrt.unrotate_pole(rot_lon, rot_lat, pole_lon, pole_lat)

    # Print to check conversion
    print ('******')
    print ('Test of unrotated coordinate grid: ')
    print ('Rotated lon coord = ', rot_lon[250])
    print ('Rotated lat coord = ', rot_lat[250])
    print ('Lon = ', lon[250])
    print ('Lat = ', lat[250])
    print (' ')

    # ******
    # write to csv file
    # ******

    # print '******'
    # print 'Writing unrotated coordinate grid to file:'
    # print ''
    # lonp = pd.DataFrame(lon)
    # latp = pd.DataFrame(lat)
    # dat = {'Latitude': lat, 'Longitude': lon}
    # df = pd.DataFrame(dat,columns=['Latitude','Longitude'])
    # df.to_csv('POSITION_UNROTATED.csv',  sep = " ")
    # print '... finished!'
    # print ''
    # print '******'

    return lon, lat

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print ('Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4]))
    print ('')

    return drift_index

def inIce(data):

    ###################################
    ## DEFINE IN ICE PERIOD
    ###################################
    # Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    Aug_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8),data.values[:,3]>=0))
    # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=13,data.values[:,1]==8),data.values[:,3]>=0))
    # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ('CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print ('')
    print ('Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print ('Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print ('Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print ('Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))
    print ('')

    return inIce_index

def writeoutGrid(tim, lat, lon, date):

    import pandas as pd

    # ******
    # write to csv file
    # ******

    print ('******')
    print ('Writing ' + date + ' grid to file:')
    print ('')
    dat = np.zeros([len(tim), 3])
    dat[:,0] = tim
    dat[:,1] = lon
    dat[:,2] = lat
    df = pd.DataFrame(dat)
    filename = 'AUX_DATA/' + date + '_ShipTrack_GRIDDED.csv'
    df.to_csv(filename,  sep = " ")
    print ('... finished!')
    print ('')
    print ('******')

def gridShipTrack(cube, xoffset, yoffset):

    import iris.plot as iplt

    # cube.dim_coords[1].coord_system

    ###---------------------------------
    ### 12th August 2018
    ###---------------------------------

    date = '20180812'

    ### HOURS 0 TO 19 IN REGION WHERE BCS CAUSE PROBLEMS

    # ### box pick 0-1h
    lon0 = np.array([253,253,252,252,252,251])
    lat0 = np.array([489,490,490,491,492,492])
    tim0 = np.zeros([np.size(lon0)])
    # for i in range(0,np.size(lon0)):
    #     iplt.scatter(cube.dim_coords[2][lon0[i] + xoffset], cube.dim_coords[1][lat0[i] + yoffset],color='yellow')
    #
    # ### box pick 1-2h
    lon1 = np.array([251,251])
    lat1 = np.array([492,493])
    tim1 = np.zeros([np.size(lon1)])
    tim1[:] = 1.0
    tim_128 = np.append(tim0, tim1)
    lat_128 = np.append(lat0, lat1)
    lon_128 = np.append(lon0, lon1)
    # for i in range(0,np.size(lon1)):
    #     iplt.scatter(cube.dim_coords[2][lon1[i] + xoffset], cube.dim_coords[1][lat1[i] + yoffset],color='blue')
    #
    # ### box pick 2-3h
    lon2 = np.array([251,251])
    lat2 = np.array([493,492])
    tim2 = np.zeros([np.size(lon2)])
    tim2[:] = 2.0
    tim_128 = np.append(tim_128, tim2)
    lat_128 = np.append(lat_128, lat2)
    lon_128 = np.append(lon_128, lon2)
    # for i in range(0,np.size(lon2)):
    #     iplt.scatter(cube.dim_coords[2][lon2[i] + xoffset], cube.dim_coords[1][lat2[i] + yoffset],color='green')
    #
    # ### box pick 3-18h
    tim3 = np.arange(3,18)
    lon3 = np.zeros([np.size(tim3)])
    lon3[:] = 251
    lat3 = np.zeros([np.size(tim3)])
    lat3[:] = 492
    tim_128 = np.append(tim_128, tim3)
    lat_128 = np.append(lat_128, lat3)
    lon_128 = np.append(lon_128, lon3)
    # for i in range(0,np.size(lon3)):
    #     iplt.scatter(cube.dim_coords[2][int(lon3[i] + xoffset)], cube.dim_coords[1][int(lat3[i] + yoffset)],color='black')
    #
    # ### box pick 18-19h
    lon18 = np.array([251,251,251])
    lat18 = np.array([492,491,490])
    tim18 = np.zeros([np.size(lon18)])
    tim18[:] = 18.0
    tim_128 = np.append(tim_128, tim18)
    lat_128 = np.append(lat_128, lat18)
    lon_128 = np.append(lon_128, lon18)
    # for i in range(0,np.size(lon18)):
    #     iplt.scatter(cube.dim_coords[2][lon18[i] + xoffset], cube.dim_coords[1][lat18[i] + yoffset],color='red')


    ### plot from here for 12th Aug

    ### box pick 19-20h
    lon19 = np.array([251,251,251,252,252])
    lat19 = np.array([490,489,488,488,487])
    tim19 = np.zeros([np.size(lon19)])
    tim19[:] = 19.0
    tim_128 = np.append(tim_128, tim19)
    lat_128 = np.append(lat_128, lat19)
    lon_128 = np.append(lon_128, lon19)
    # for i in range(0,np.size(lon19)):
    #     iplt.scatter(cube.dim_coords[2][lon19[i] + xoffset], cube.dim_coords[1][lat19[i] + yoffset],color='blue')

    ### box pick 20-21h
    lon20 = np.array([252,252,252,252,252])
    lat20 = np.array([487,486,485,484,483])
    tim20 = np.zeros([np.size(lon20)])
    tim20[:] = 20.0
    tim_128 = np.append(tim_128, tim20)
    lat_128 = np.append(lat_128, lat20)
    lon_128 = np.append(lon_128, lon20)
    # for i in range(0,np.size(lon20)):
    #     iplt.scatter(cube.dim_coords[2][lon20[i] + xoffset], cube.dim_coords[1][lat20[i] + yoffset],color='green')

    ### box pick 21-22h
    lon21 = np.array([252,251,250,249,248])
    lat21 = np.array([483,483,483,483,483])
    tim21 = np.zeros([np.size(lon21)])
    tim21[:] = 21.0
    tim_128 = np.append(tim_128, tim21)
    lat_128 = np.append(lat_128, lat21)
    lon_128 = np.append(lon_128, lon21)
    # for i in range(0,np.size(lon21)):
    #     iplt.scatter(cube.dim_coords[2][lon21[i] + xoffset], cube.dim_coords[1][lat21[i] + yoffset],color='black')

    ### box pick 22-23h
    lon22 = np.array([248,248,248,248,248])
    lat22 = np.array([483,482,481,480,479])
    tim22 = np.zeros([np.size(lon22)])
    tim22[:] = 22.0
    tim_128 = np.append(tim_128, tim22)
    lat_128 = np.append(lat_128, lat22)
    lon_128 = np.append(lon_128, lon22)
    # for i in range(0,np.size(lon22)):
    #     iplt.scatter(cube.dim_coords[2][lon22[i] + xoffset], cube.dim_coords[1][lat22[i] + yoffset],color='red')

    ### box pick 23-00h
    lon23 = np.array([248,249,249,250,250])
    lat23 = np.array([479,479,478,478,477])
    tim23 = np.zeros([np.size(lon23)])
    tim23[:] = 23.0
    tim_128 = np.append(tim_128, tim23)
    lat_128 = np.append(lat_128, lat23)
    lon_128 = np.append(lon_128, lon23)
    # for i in range(0,np.size(lon23)):
    #     iplt.scatter(cube.dim_coords[2][lon23[i] + xoffset], cube.dim_coords[1][lat23[i] + yoffset],color='magenta')

    # ******
    # write out index arrays
    # ******

    # out = writeoutGrid(tim_128, lat_128, lon_128, date)

    ###---------------------------------
    ### 13th August 2018
    ###---------------------------------

    date = '20180813'

    ### box pick 0-1h
    lon0 = np.array([250,251,252,252])
    lat0 = np.array([477,477,477,476])
    tim0 = np.zeros([np.size(lon0)])
    tim0[:] = 0.0
    # for i in range(0,np.size(lon0)):
    #     iplt.scatter(cube.dim_coords[2][lon0[i] + xoffset], cube.dim_coords[1][lat0[i] + yoffset],color='yellow')

    ### box pick 1-2h    (13th aug)
    lon1 = np.array([252,253])
    lat1 = np.array([476,476])
    tim1 = np.zeros([np.size(lon1)])
    tim1[:] = 1.0
    tim_138 = np.append(tim0, tim1)
    lat_138 = np.append(lat0, lat1)
    lon_138 = np.append(lon0, lon1)
    # for i in range(0,np.size(lon1)):
    #     iplt.scatter(cube.dim_coords[2][lon1[i] + xoffset], cube.dim_coords[1][lat1[i] + yoffset],color='black')

    ### box pick 2-3h    (13th aug)
    lon2 = np.array([253,253,254])
    lat2 = np.array([476,475,475])
    tim2 = np.zeros([np.size(lon2)])
    tim2[:] = 2.0
    tim_138 = np.append(tim_138, tim2)
    lat_138 = np.append(lat_138, lat2)
    lon_138 = np.append(lon_138, lon2)
    # for i in range(0,np.size(lon2)):
    #     iplt.scatter(cube.dim_coords[2][lon2[i] + xoffset], cube.dim_coords[1][lat2[i] + yoffset],color='red')

    ### box pick 3-4h    (13th aug)
    lon3 = np.array([254,254,255,255])
    lat3 = np.array([475,474,474,473])
    tim3 = np.zeros([np.size(lon3)])
    tim3[:] = 3.0
    tim_138 = np.append(tim_138, tim3)
    lat_138 = np.append(lat_138, lat3)
    lon_138 = np.append(lon_138, lon3)
    # for i in range(0,np.size(lon3)):
    #     iplt.scatter(cube.dim_coords[2][lon3[i] + xoffset], cube.dim_coords[1][lat3[i] + yoffset],color='blue')

    ### box pick 4-5h    (13th aug)
    lon4 = np.array([255,256])
    lat4 = np.array([473,473])
    tim4 = np.zeros([np.size(lon4)])
    tim4[:] = 4.0
    tim_138 = np.append(tim_138, tim4)
    lat_138 = np.append(lat_138, lat4)
    lon_138 = np.append(lon_138, lon4)
    # for i in range(0,np.size(lon4)):
    #     iplt.scatter(cube.dim_coords[2][lon4[i] + xoffset], cube.dim_coords[1][lat4[i] + yoffset],color='green')

    ### box pick 5-6h
    lon5 = np.array([256,256,255,255])
    lat5 = np.array([473,472,472,471])
    tim5 = np.zeros([np.size(lon5)])
    tim5[:] = 5.0
    tim_138 = np.append(tim_138, tim5)
    lat_138 = np.append(lat_138, lat5)
    lon_138 = np.append(lon_138, lon5)
    # for i in range(0,np.size(lon5)):
    #     iplt.scatter(cube.dim_coords[2][lon5[i] + xoffset], cube.dim_coords[1][lat5[i] + yoffset],color='black')

    ### box pick 6-17h
    tim6 = np.arange(6,17)
    lon6 = np.zeros([np.size(tim6)])
    lon6[:] = 255
    lat6 = np.zeros([np.size(tim6)])
    lat6[:] = 471
    tim_138 = np.append(tim_138, tim6)
    lat_138 = np.append(lat_138, lat6)
    lon_138 = np.append(lon_138, lon6)
    # for i in range(0,np.size(lon6)):
    #     iplt.scatter(cube.dim_coords[2][int(lon6[i]) + xoffset], cube.dim_coords[1][int(lat6[i]) + yoffset],color='red')

    ### box pick 17-18h
    lon17 = np.array([255,255])
    lat17 = np.array([471,472])
    tim17 = np.zeros([np.size(lon17)])
    tim17[:] = 17.0
    tim_138 = np.append(tim_138, tim17)
    lat_138 = np.append(lat_138, lat17)
    lon_138 = np.append(lon_138, lon17)
    # for i in range(0,np.size(lon17)):
    #     iplt.scatter(cube.dim_coords[2][lon17[i] + xoffset], cube.dim_coords[1][lat17[i] + yoffset],color='blue')

    ### box pick 18-23h
    tim18 = np.arange(18,23)
    lon18 = np.zeros([np.size(tim18)])
    lon18[:] = 255
    lat18 = np.zeros([np.size(tim18)])
    lat18[:] = 472
    tim_138 = np.append(tim_138, tim18)
    lat_138 = np.append(lat_138, lat18)
    lon_138 = np.append(lon_138, lon18)
    # for i in range(0,np.size(lon18)):
    #     iplt.scatter(cube.dim_coords[2][int(lon18[i]) + xoffset], cube.dim_coords[1][int(lat18[i]) + yoffset],color='green')

    ### box pick 23-00h
    lon23 = np.array([255,255])
    lat23 = np.array([472,471])
    tim23 = np.zeros([np.size(lon23)])
    tim23[:] = 23.0
    tim_138 = np.append(tim_138, tim23)
    lat_138 = np.append(lat_138, lat23)
    lon_138 = np.append(lon_138, lon23)
    # for i in range(0,np.size(lon23)):
    #     iplt.scatter(cube.dim_coords[2][lon23[i] + xoffset], cube.dim_coords[1][lat23[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_138, lat_138, lon_138, date)

    ###---------------------------------
    ### 14th August 2018
    ###---------------------------------

    date = '20180814'

    ### box pick 0-4h
    tim0 = np.arange(0,4)
    lon0 = np.zeros([np.size(tim0)])
    lon0[:] = 255
    lat0 = np.zeros([np.size(tim0)])
    lat0[:] = 471
    # for i in range(0,np.size(lon0)):
    #     iplt.scatter(cube.dim_coords[2][int(lon0[i]) + xoffset], cube.dim_coords[1][int(lat0[i]) + yoffset],color='yellow')

    ### box pick 4-5h
    lon4 = np.array([255,255,254])
    lat4 = np.array([471,470,470])
    tim4 = np.zeros([np.size(lon4)])
    tim4[:] = 4.0
    tim_148 = np.append(tim0, tim4)
    lat_148 = np.append(lat0, lat4)
    lon_148 = np.append(lon0, lon4)
    # for i in range(0,np.size(lon4)):
    #     iplt.scatter(cube.dim_coords[2][lon4[i] + xoffset], cube.dim_coords[1][lat4[i] + yoffset],color='red')

    ### box pick 5-10h
    tim5 = np.arange(5,11)
    lon5 = np.zeros([np.size(tim5)])
    lon5[:] = 254
    lat5 = np.zeros([np.size(tim5)])
    lat5[:] = 470
    tim_148 = np.append(tim_148, tim5)
    lat_148 = np.append(lat_148, lat5)
    lon_148 = np.append(lon_148, lon5)
    # for i in range(0,np.size(lon5)):
    #     iplt.scatter(cube.dim_coords[2][int(lon5[i]) + xoffset], cube.dim_coords[1][int(lat5[i]) + yoffset],color='blue')

    ### box pick 11-16h
    tim11 = np.arange(11,17)
    lon11 = np.zeros([np.size(tim11)])
    lon11[:] = 254
    lat11 = np.zeros([np.size(tim11)])
    lat11[:] = 469
    tim_148 = np.append(tim_148, tim11)
    lat_148 = np.append(lat_148, lat11)
    lon_148 = np.append(lon_148, lon11)
    # for i in range(0,np.size(lon11)):
    #     iplt.scatter(cube.dim_coords[2][int(lon11[i]) + xoffset], cube.dim_coords[1][int(lat11[i]) + yoffset],color='green')

    ### box pick 16-17h
    lon16 = np.array([254,254])
    lat16 = np.array([469,468])
    tim16 = np.zeros([np.size(lon16)])
    tim16[:] = 16.0
    tim_148 = np.append(tim_148, tim16)
    lat_148 = np.append(lat_148, lat16)
    lon_148 = np.append(lon_148, lon16)
    # for i in range(0,np.size(lon16)):
    #     iplt.scatter(cube.dim_coords[2][lon16[i] + xoffset], cube.dim_coords[1][lat16[i] + yoffset],color='black')

    ### box pick 17-18h
    lon17 = np.array([254,253])
    lat17 = np.array([468,468])
    tim17 = np.zeros([np.size(lon17)])
    tim17[:] = 17.0
    tim_148 = np.append(tim_148, tim17)
    lat_148 = np.append(lat_148, lat17)
    lon_148 = np.append(lon_148, lon17)
    # for i in range(0,np.size(lon17)):
    #     iplt.scatter(cube.dim_coords[2][lon17[i] + xoffset], cube.dim_coords[1][lat17[i] + yoffset],color='black')

    ## box pick 18-21h
    tim18 = np.arange(18,22)
    lon18 = np.zeros([np.size(tim18)])
    lon18[:] = 253
    lat18 = np.zeros([np.size(tim18)])
    lat18[:] = 468
    tim_148 = np.append(tim_148, tim18)
    lat_148 = np.append(lat_148, lat18)
    lon_148 = np.append(lon_148, lon18)
    # for i in range(0,np.size(lon18)):
    #     iplt.scatter(cube.dim_coords[2][int(lon18[i]) + xoffset], cube.dim_coords[1][int(lat18[i]) + yoffset],color='red')

    ## box pick 18-21h
    tim22 = np.arange(22,24)
    lon22 = np.zeros([np.size(tim22)])
    lon22[:] = 253
    lat22 = np.zeros([np.size(tim22)])
    lat22[:] = 467
    tim_148 = np.append(tim_148, tim22)
    lat_148 = np.append(lat_148, lat22)
    lon_148 = np.append(lon_148, lon22)
    # for i in range(0,np.size(lon22)):
    #     iplt.scatter(cube.dim_coords[2][int(lon22[i]) + xoffset], cube.dim_coords[1][int(lat22[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_148, lat_148, lon_148, date)

    ###---------------------------------
    ### 15th August 2018
    ###---------------------------------

    date = '20180815'
    tim_158 = np.array([])
    lat_158 = np.array([])
    lon_158 = np.array([])

    ### box pick 0-1h
    lon = np.array([253])
    lat = np.array([467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([253,253])
    lat = np.array([467,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-4h
    tim = np.arange(2,4)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 4-5h
    lon = np.array([253,253])
    lat = np.array([466,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 5-7h
    tim = np.arange(5,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 465
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 7-8h
    lon = np.array([253,253])
    lat = np.array([465,464])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 8-9h
    lon = np.array([253])
    lat = np.array([464])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([253,253])
    lat = np.array([464,463])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 10-12h
    tim = np.arange(10,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 463
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 12-13h
    lon = np.array([253,253])
    lat = np.array([463,462])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 13-15h
    tim = np.arange(13,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 462
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 15-16h
    lon = np.array([253,252])
    lat = np.array([462,462])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 16-17h
    lon = np.array([252,252])
    lat = np.array([462,461])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 18-21h
    tim = np.arange(18,21)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 461
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 21-22h
    lon = np.array([252,252])
    lat = np.array([461,460])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 18-21h
    tim = np.arange(22,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 460
    tim_158 = np.append(tim_158, tim)
    lat_158 = np.append(lat_158, lat)
    lon_158 = np.append(lon_158, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_158, lat_158, lon_158, date)

    ###---------------------------------
    ### 16th August 2018
    ###---------------------------------

    date = '20180816'
    tim_168 = np.array([])
    lat_168 = np.array([])
    lon_168 = np.array([])

    ### box pick 0-2h
    tim = np.arange(0,2)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 460
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 2-3h
    lon = np.array([252,252])
    lat = np.array([460,459])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 3-4h
    lon = np.array([252,251])
    lat = np.array([459,459])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 4-13h
    tim = np.arange(4,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 251
    lat = np.zeros([np.size(tim)])
    lat[:] = 459
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-15h
    lon = np.array([251,250])
    lat = np.array([459,459])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 4-18h
    tim = np.arange(14,18)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 459
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 18-19h
    lon = np.array([250,250])
    lat = np.array([459,460])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 4-18h
    tim = np.arange(19,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 460
    tim_168 = np.append(tim_168, tim)
    lat_168 = np.append(lat_168, lat)
    lon_168 = np.append(lon_168, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_168, lat_168, lon_168, date)

    ###---------------------------------
    ### 17th August 2018
    ###---------------------------------

    date = '20180817'
    tim_178 = np.array([])
    lat_178 = np.array([])
    lon_178 = np.array([])

    ### box pick 0-1h
    lon = np.array([250,250])
    lat = np.array([460,461])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-5h
    tim = np.arange(1,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 461
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 5-6h
    lon = np.array([250,250])
    lat = np.array([461,462])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
        # iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 6-9h
    tim = np.arange(6,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 462
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 9-10h
    lon = np.array([250,250])
    lat = np.array([462,463])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 10-15h
    tim = np.arange(10,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 463
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([250,251,251])
    lat = np.array([463,463,464])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 16-20h
    tim = np.arange(16,20)
    lon = np.zeros([np.size(tim)])
    lon[:] = 251
    lat = np.zeros([np.size(tim)])
    lat[:] = 464
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 20-21h
    lon = np.array([251,251])
    lat = np.array([464,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 21-22h
    lon = np.array([251,252])
    lat = np.array([465,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 22-0h
    tim = np.arange(22,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 465
    tim_178 = np.append(tim_178, tim)
    lat_178 = np.append(lat_178, lat)
    lon_178 = np.append(lon_178, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_178, lat_178, lon_178, date)

    ###---------------------------------
    ### 18th August 2018
    ###---------------------------------

    date = '20180818'
    tim_188 = np.array([])
    lat_188 = np.array([])
    lon_188 = np.array([])

    ### box pick 0-1h
    lon = np.array([252])
    lat = np.array([465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([252,252])
    lat = np.array([465,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([252,253])
    lat = np.array([466,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-5h
    tim = np.arange(3,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([253,253])
    lat = np.array([466,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-9h
    tim = np.arange(6,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 467
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 9-10h
    lon = np.array([253,254])
    lat = np.array([467,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 10-13h
    tim = np.arange(10,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 254
    lat = np.zeros([np.size(tim)])
    lat[:] = 467
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')


    ### box pick 13-14h
    lon = np.array([254,254])
    lat = np.array([467,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 14-15h
    lon = np.array([254,255])
    lat = np.array([468,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-0h
    tim = np.arange(15,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 255
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_188 = np.append(tim_188, tim)
    lat_188 = np.append(lat_188, lat)
    lon_188 = np.append(lon_188, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_188, lat_188, lon_188, date)

    ###---------------------------------
    ### 19th August 2018
    ###---------------------------------

    date = '20180819'
    tim_198 = np.array([])
    lat_198 = np.array([])
    lon_198 = np.array([])

    ### box pick 0-1h
    lon = np.array([255,255,256])
    lat = np.array([468,469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-7h
    tim = np.arange(1,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 256
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-8h
    lon = np.array([256,255])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 8-11h
    tim = np.arange(8,11)
    lon = np.zeros([np.size(tim)])
    lon[:] = 255
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 11-12h
    lon = np.array([255,255])
    lat = np.array([469,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 12-15h
    tim = np.arange(12,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 255
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([255,254])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 16-21h
    tim = np.arange(16,21)
    lon = np.zeros([np.size(tim)])
    lon[:] = 254
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 21-22h
    lon = np.array([254,254])
    lat = np.array([470,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 16-21h
    tim = np.arange(22,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 254
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_198 = np.append(tim_198, tim)
    lat_198 = np.append(lat_198, lat)
    lon_198 = np.append(lon_198, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_198, lat_198, lon_198, date)

    ###---------------------------------
    ### 20th August 2018
    ###---------------------------------

    date = '20180820'
    tim_208 = np.array([])
    lat_208 = np.array([])
    lon_208 = np.array([])

    ### box pick 0-1h
    lon = np.array([254,253])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-5h
    tim = np.arange(1,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 253
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 5-6h
    lon = np.array([253,253])
    lat = np.array([471,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 6-7h
    lon = np.array([253,252])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 6-12h
    tim = np.arange(6,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 12-13h
    lon = np.array([252,251])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 13-17h
    tim = np.arange(13,17)
    lon = np.zeros([np.size(tim)])
    lon[:] = 251
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 17-18h
    lon = np.array([251,250])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 17-23h
    tim = np.arange(17,23)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 17-18h
    lon = np.array([250,249])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_208 = np.append(tim_208, tim)
    lat_208 = np.append(lat_208, lat)
    lon_208 = np.append(lon_208, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_208, lat_208, lon_208, date)

    ###---------------------------------
    ### 21st August 2018
    ###---------------------------------

    date = '20180821'
    tim_218 = np.array([])
    lat_218 = np.array([])
    lon_218 = np.array([])

    ### box pick 0-1h
    lon = np.array([249,249])
    lat = np.array([472,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([249])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([249,248])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 3-5h
    tim = np.arange(3,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([248,247])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 6-8h
    tim = np.arange(6,8)
    lon = np.zeros([np.size(tim)])
    lon[:] = 247
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 8-9h
    lon = np.array([247,246])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([246,246])
    lat = np.array([471,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 10-11h
    lon = np.array([246])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 11-12h
    lon = np.array([246,245])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 12-13h
    lon = np.array([245])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 13-14h
    lon = np.array([245,244])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 14-16h
    tim = np.arange(14,16)
    lon = np.zeros([np.size(tim)])
    lon[:] = 244
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 16-17h
    lon = np.array([244,243])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 17-18h
    lon = np.array([243])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 18-19h
    lon = np.array([243,242])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 19-21h
    tim = np.arange(19,21)
    lon = np.zeros([np.size(tim)])
    lon[:] = 242
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 21-22h
    lon = np.array([242,241])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 22-23h
    lon = np.array([241])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 22-23h
    lon = np.array([241,240])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_218 = np.append(tim_218, tim)
    lat_218 = np.append(lat_218, lat)
    lon_218 = np.append(lon_218, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_218, lat_218, lon_218, date)

    ###---------------------------------
    ### 22nd August 2018
    ###---------------------------------

    date = '20180822'
    tim_228 = np.array([])
    lat_228 = np.array([])
    lon_228 = np.array([])

    ### box pick 0-1h
    lon = np.array([240])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
        # iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([240,239])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([239])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([239,239])
    lat = np.array([470,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 4-5
    lon = np.array([239,238])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 5-7h
    tim = np.arange(5,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 238
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-8h
    lon = np.array([238,237])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 8-10h
    tim = np.arange(8,10)
    lon = np.zeros([np.size(tim)])
    lon[:] = 237
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 10-11h
    lon = np.array([237,237])
    lat = np.array([469,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 11-12h
    lon = np.array([237])
    lat = np.array([468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 12-13h
    lon = np.array([237,236])
    lat = np.array([468,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 13-18h
    tim = np.arange(13,18)
    lon = np.zeros([np.size(tim)])
    lon[:] = 236
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 18-19h
    lon = np.array([236,236])
    lat = np.array([468,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 18-0h
    tim = np.arange(18,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 236
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_228 = np.append(tim_228, tim)
    lat_228 = np.append(lat_228, lat)
    lon_228 = np.append(lon_228, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_228, lat_228, lon_228, date)

    ###---------------------------------
    ### 23rd August 2018
    ###---------------------------------

    date = '20180823'
    tim_238 = np.array([])
    lat_238 = np.array([])
    lon_238 = np.array([])

    ### box pick 0-1h
    lon = np.array([236])
    lat = np.array([469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
        # iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([236,235])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-9h
    tim = np.arange(2,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 235
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([235,236])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 10-12h
    tim = np.arange(10,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 236
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')
    #
    ### box pick 12-13h
    lon = np.array([236,236])
    lat = np.array([469,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 13-18h
    tim = np.arange(13,18)
    lon = np.zeros([np.size(tim)])
    lon[:] = 236
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 18-19h
    lon = np.array([236,237])
    lat = np.array([468,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 19-21h
    tim = np.arange(19,21)
    lon = np.zeros([np.size(tim)])
    lon[:] = 237
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 21-22h
    lon = np.array([237,237])
    lat = np.array([468,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 22-23h
    lon = np.array([237,238])
    lat = np.array([467,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 23-0h
    lon = np.array([238])
    lat = np.array([467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_238 = np.append(tim_238, tim)
    lat_238 = np.append(lat_238, lat)
    lon_238 = np.append(lon_238, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_238, lat_238, lon_238, date)

    ###---------------------------------
    ### 24th August 2018
    ###---------------------------------

    date = '20180824'
    tim_248 = np.array([])
    lat_248 = np.array([])
    lon_248 = np.array([])

    ### box pick 0-3h
    tim = np.arange(0,3)
    lon = np.zeros([np.size(tim)])
    lon[:] = 238
    lat = np.zeros([np.size(tim)])
    lat[:] = 467
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 3-4h
    lon = np.array([238,238,239])
    lat = np.array([467,466,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 4-7h
    tim = np.arange(4,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 239
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 7-8h
    lon = np.array([239,240])
    lat = np.array([466,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 8-11h
    tim = np.arange(8,11)
    lon = np.zeros([np.size(tim)])
    lon[:] = 240
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 11-12h
    lon = np.array([240,241])
    lat = np.array([466,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 12-13h
    lon = np.array([241,241])
    lat = np.array([466,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 13-15h
    tim = np.arange(13,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 241
    lat = np.zeros([np.size(tim)])
    lat[:] = 465
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 15-16h
    lon = np.array([241,242])
    lat = np.array([465,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 16-19h
    tim = np.arange(16,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 242
    lat = np.zeros([np.size(tim)])
    lat[:] = 465
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 19-20h
    lon = np.array([242,243])
    lat = np.array([465,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 19-23h
    tim = np.arange(19,23)
    lon = np.zeros([np.size(tim)])
    lon[:] = 243
    lat = np.zeros([np.size(tim)])
    lat[:] = 465
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 23-0h
    lon = np.array([243,244])
    lat = np.array([465,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_248 = np.append(tim_248, tim)
    lat_248 = np.append(lat_248, lat)
    lon_248 = np.append(lon_248, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_248, lat_248, lon_248, date)

    ###---------------------------------
    ### 25th August 2018
    ###---------------------------------

    date = '20180825'
    tim_258 = np.array([])
    lat_258 = np.array([])
    lon_258 = np.array([])

    ### box pick 0-1h
    lon = np.array([244])
    lat = np.array([465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([244,244])
    lat = np.array([465,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-4h
    tim = np.arange(2,4)
    lon = np.zeros([np.size(tim)])
    lon[:] = 244
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 4-5h
    lon = np.array([244,245])
    lat = np.array([466,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-7h
    tim = np.arange(5,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 245
    lat = np.zeros([np.size(tim)])
    lat[:] = 466
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 7-8h
    lon = np.array([245,245])
    lat = np.array([466,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 8-10h
    tim = np.arange(8,10)
    lon = np.zeros([np.size(tim)])
    lon[:] = 245
    lat = np.zeros([np.size(tim)])
    lat[:] = 467
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 10-11h
    lon = np.array([245,246])
    lat = np.array([467,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 11-13h
    tim = np.arange(11,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 246
    lat = np.zeros([np.size(tim)])
    lat[:] = 467
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 13-14h
    lon = np.array([246,246])
    lat = np.array([467,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 14-16h
    tim = np.arange(14,16)
    lon = np.zeros([np.size(tim)])
    lon[:] = 246
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 16-17h
    lon = np.array([246,246])
    lat = np.array([468,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 17-19h
    tim = np.arange(17,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 246
    lat = np.zeros([np.size(tim)])
    lat[:] = 469
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 19-20h
    lon = np.array([246,246,247])
    lat = np.array([469,470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-0h
    tim = np.arange(20,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 247
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_258 = np.append(tim_258, tim)
    lat_258 = np.append(lat_258, lat)
    lon_258 = np.append(lon_258, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_258, lat_258, lon_258, date)

    ###---------------------------------
    ### 26th August 2018
    ###---------------------------------

    date = '20180826'
    tim_268 = np.array([])
    lat_268 = np.array([])
    lon_268 = np.array([])

    ### box pick 0-1h
    lon = np.array([247])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([247,247,248])
    lat = np.array([470,471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-5h
    tim = np.arange(2,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 5-6h
    lon = np.array([248,248])
    lat = np.array([471,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 6-10h
    tim = np.arange(6,10)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 10-11h
    lon = np.array([248,249])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 11-13h
    tim = np.arange(11,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 249
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 13-14h
    lon = np.array([249,248])
    lat = np.array([472,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-17h
    tim = np.arange(14,17)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 17-18h
    lon = np.array([248,248])
    lat = np.array([472,473])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 18-0h
    tim = np.arange(18,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 473
    tim_268 = np.append(tim_268, tim)
    lat_268 = np.append(lat_268, lat)
    lon_268 = np.append(lon_268, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_268, lat_268, lon_268, date)

    ###---------------------------------
    ### 27th August 2018
    ###---------------------------------

    date = '20180827'
    tim_278 = np.array([])
    lat_278 = np.array([])
    lon_278 = np.array([])

    ### box pick 0-3h
    tim = np.arange(0,3)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 473
    tim_278 = np.append(tim_278, tim)
    lat_278 = np.append(lat_278, lat)
    lon_278 = np.append(lon_278, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 3-4h
    lon = np.array([248,247])
    lat = np.array([473,473])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_278 = np.append(tim_278, tim)
    lat_278 = np.append(lat_278, lat)
    lon_278 = np.append(lon_278, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 4-18h
    tim = np.arange(4,18)
    lon = np.zeros([np.size(tim)])
    lon[:] = 247
    lat = np.zeros([np.size(tim)])
    lat[:] = 473
    tim_278 = np.append(tim_278, tim)
    lat_278 = np.append(lat_278, lat)
    lon_278 = np.append(lon_278, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 18-19h
    lon = np.array([247,247])
    lat = np.array([473,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_278 = np.append(tim_278, tim)
    lat_278 = np.append(lat_278, lat)
    lon_278 = np.append(lon_278, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 19-0h
    tim = np.arange(19,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 247
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_278 = np.append(tim_278, tim)
    lat_278 = np.append(lat_278, lat)
    lon_278 = np.append(lon_278, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_278, lat_278, lon_278, date)

    ###---------------------------------
    ### 28th August 2018
    ###---------------------------------

    date = '20180828'
    tim_288 = np.array([])
    lat_288 = np.array([])
    lon_288 = np.array([])

    ### box pick 0-1h
    lon = np.array([247])
    lat = np.array([472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([247,247])
    lat = np.array([472,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([247])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([247,248])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 4-7h
    tim = np.arange(4,7)
    lon = np.zeros([np.size(tim)])
    lon[:] = 248
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 7-8h
    lon = np.array([248,249])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 8-10h
    tim = np.arange(8,10)
    lon = np.zeros([np.size(tim)])
    lon[:] = 249
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 10-11h
    lon = np.array([249,249,250])
    lat = np.array([471,470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 11-13h
    tim = np.arange(11,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 250
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 13-14h
    lon = np.array([250,251])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 14-16h
    tim = np.arange(14,16)
    lon = np.zeros([np.size(tim)])
    lon[:] = 251
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 16-17h
    lon = np.array([251,252])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 17-19h
    tim = np.arange(17,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 252
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 19-20h
    lon = np.array([252,253])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 20-21h
    lon = np.array([253,253])
    lat = np.array([470,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 21-22h
    lon = np.array([253,254])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 22-23h
    lon = np.array([254])
    lat = np.array([469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 23-0h
    lon = np.array([254,255])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_288 = np.append(tim_288, tim)
    lat_288 = np.append(lat_288, lat)
    lon_288 = np.append(lon_288, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_288, lat_288, lon_288, date)

    ###---------------------------------
    ### 29th August 2018
    ###---------------------------------

    date = '20180829'
    tim_298 = np.array([])
    lat_298 = np.array([])
    lon_298 = np.array([])

    ### box pick 0-1h
    lon = np.array([255])
    lat = np.array([469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([255,256])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([256,257])
    lat = np.array([469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([257,257])
    lat = np.array([469,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 4-5h
    lon = np.array([257,258])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 5-6h
    lon = np.array([258,259])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 6-7h
    lon = np.array([259])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 7-8h
    lon = np.array([259,260])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 8-9h
    lon = np.array([260,261])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 9-10h
    lon = np.array([261,261])
    lat = np.array([470,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 10-11h
    lon = np.array([261,262])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([262])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 12-13h
    lon = np.array([262,263])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 13-15h
    tim = np.arange(13,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 263
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([263,264])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 16-17h
    lon = np.array([264])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 17-18h
    lon = np.array([264,265])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 18-20h
    tim = np.arange(18,20)
    lon = np.zeros([np.size(tim)])
    lon[:] = 265
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-21h
    lon = np.array([265,266])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 21-23h
    tim = np.arange(21,23)
    lon = np.zeros([np.size(tim)])
    lon[:] = 266
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 20-21h
    lon = np.array([266,267,267])
    lat = np.array([471,471,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_298 = np.append(tim_298, tim)
    lat_298 = np.append(lat_298, lat)
    lon_298 = np.append(lon_298, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_298, lat_298, lon_298, date)

    ###---------------------------------
    ### 30th August 2018
    ###---------------------------------

    date = '20180830'
    tim_308 = np.array([])
    lat_308 = np.array([])
    lon_308 = np.array([])

    ### box pick 0-2h
    tim = np.arange(0,2)
    lon = np.zeros([np.size(tim)])
    lon[:] = 267
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 2-3h
    lon = np.array([267,268])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 3-4h
    lon = np.array([268])
    lat = np.array([470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 4-5h
    lon = np.array([268,269])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([269,269])
    lat = np.array([470,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 6-7h
    lon = np.array([269,270])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 7-8h
    lon = np.array([270,271])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 8-9h
    lon = np.array([271])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 9-10h
    lon = np.array([271,272])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 10-11h
    lon = np.array([272])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 11-12h
    lon = np.array([272,272,273])
    lat = np.array([471,470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 12-14h
    tim = np.arange(12,14)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-15h
    lon = np.array([273,274])
    lat = np.array([470,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 15-16h
    lon = np.array([274,274])
    lat = np.array([470,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 16-17h
    lon = np.array([274])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 17-18h
    lon = np.array([274,275])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 18-19h
    lon = np.array([275])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 19-20h
    lon = np.array([275,276])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 19-22h
    tim = np.arange(19,22)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 22-23h
    lon = np.array([276,277])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 23-0h
    lon = np.array([277])
    lat = np.array([471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_308 = np.append(tim_308, tim)
    lat_308 = np.append(lat_308, lat)
    lon_308 = np.append(lon_308, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_308, lat_308, lon_308, date)

    ###---------------------------------
    ### 31st August 2018
    ###---------------------------------

    date = '20180831'
    tim_318 = np.array([])
    lat_318 = np.array([])
    lon_318 = np.array([])

    ### box pick 0-4h
    tim = np.arange(0,4)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 4-5h
    lon = np.array([277,277])
    lat = np.array([471,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 5-11h
    tim = np.arange(5,11)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([277,277])
    lat = np.array([472,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 12-13h
    lon = np.array([277,276])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 13-15h
    tim = np.arange(13,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([276,275])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 16-19h
    tim = np.arange(16,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 19-20h
    lon = np.array([275,275])
    lat = np.array([471,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 20-22h
    tim = np.arange(20,22)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 22-23h
    lon = np.array([275,275,274])
    lat = np.array([470,469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 23-0h
    lon = np.array([274,274])
    lat = np.array([469,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    # out = writeoutGrid(tim_318, lat_318, lon_318, date)

    ###---------------------------------
    ### 1st September 2018
    ###---------------------------------

    date = '20180901'
    tim_019 = np.array([])
    lat_019 = np.array([])
    lon_019 = np.array([])

    ### box pick 0-2h
    tim = np.arange(0,2)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 468
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 2-3h
    lon = np.array([274,274,273])
    lat = np.array([468,467,467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 3-4h
    lon = np.array([273])
    lat = np.array([467])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 4-5h
    lon = np.array([273,273])
    lat = np.array([467,466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([273])
    lat = np.array([466])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-7h
    lon = np.array([273,273])
    lat = np.array([466,465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-8h
    lon = np.array([273])
    lat = np.array([465])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 8-9h
    lon = np.array([273,273])
    lat = np.array([465,464])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 9-10h
    lon = np.array([273,273])
    lat = np.array([464,463])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 10-11h
    lon = np.array([273,273])
    lat = np.array([463,462])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 11-12h
    lon = np.array([273])
    lat = np.array([462])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 12-13h
    lon = np.array([273,273])
    lat = np.array([462,461])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 13-14h
    lon = np.array([273])
    lat = np.array([461])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 14-15h
    lon = np.array([273,273])
    lat = np.array([461,460])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([273])
    lat = np.array([460])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 16-17h
    lon = np.array([273,273])
    lat = np.array([460,459])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 17-18h
    lon = np.array([273,274])
    lat = np.array([459,459])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 18-20h
    tim = np.arange(18,20)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 459
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-21h
    lon = np.array([274,274])
    lat = np.array([459,458])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 21-22h
    lon = np.array([274,275])
    lat = np.array([458,458])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 22-23h
    lon = np.array([275,275])
    lat = np.array([458,457])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 23-0h
    lon = np.array([275])
    lat = np.array([457])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_019 = np.append(tim_019, tim)
    lat_019 = np.append(lat_019, lat)
    lon_019 = np.append(lon_019, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_019, lat_019, lon_019, date)

    ###---------------------------------
    ### 2nd September 2018
    ###---------------------------------

    date = '20180902'
    tim_029 = np.array([])
    lat_029 = np.array([])
    lon_029 = np.array([])

    ### box pick 0-1h
    lon = np.array([275,275])
    lat = np.array([457,456])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([275])
    lat = np.array([456])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([275,275])
    lat = np.array([456,455])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-5h
    tim = np.arange(3,5)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 455
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([275,276,276])
    lat = np.array([455,455,454])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-7h
    lon = np.array([276])
    lat = np.array([454])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-8h
    lon = np.array([276,276])
    lat = np.array([454,453])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 8-9h
    lon = np.array([276,277])
    lat = np.array([453,453])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 9-10h
    lon = np.array([277,277])
    lat = np.array([453,452])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 10-11h
    lon = np.array([277,277])
    lat = np.array([452,451])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 11-12h
    lon = np.array([277])
    lat = np.array([451])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 12-13h
    lon = np.array([277,277])
    lat = np.array([451,450])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 13-14h
    lon = np.array([277])
    lat = np.array([450])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-15h
    lon = np.array([277,277])
    lat = np.array([450,449])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 15-16h
    lon = np.array([277])
    lat = np.array([449])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 16-17h
    lon = np.array([277,277])
    lat = np.array([449,448])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 17-18h
    lon = np.array([277])
    lat = np.array([448])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 18-19h
    lon = np.array([277,278,278])
    lat = np.array([448,448,447])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 19-20h
    lon = np.array([278])
    lat = np.array([447])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-21h
    lon = np.array([278,278])
    lat = np.array([447,446])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 21-22h
    lon = np.array([278,279])
    lat = np.array([446,446])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 22-23h
    lon = np.array([279,279])
    lat = np.array([446,445])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 23-0h
    lon = np.array([279,279])
    lat = np.array([445,444])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_029 = np.append(tim_029, tim)
    lat_029 = np.append(lat_029, lat)
    lon_029 = np.append(lon_029, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_029, lat_029, lon_029, date)

    ###---------------------------------
    ### 3rd September 2018
    ###---------------------------------

    date = '20180903'
    tim_039 = np.array([])
    lat_039 = np.array([])
    lon_039 = np.array([])

    ### box pick 0-1h
    lon = np.array([279])
    lat = np.array([444])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([279,279])
    lat = np.array([444,443])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([279,279])
    lat = np.array([443,442])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([279])
    lat = np.array([442])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 4-5h
    lon = np.array([279,279])
    lat = np.array([442,441])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 5-6h
    lon = np.array([279])
    lat = np.array([441])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 6-7h
    lon = np.array([279,278,278])
    lat = np.array([441,441,440])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 7-8h
    lon = np.array([278])
    lat = np.array([440])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 8-9h
    lon = np.array([278,278])
    lat = np.array([440,439])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 9-10h
    lon = np.array([278])
    lat = np.array([439])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 10-11h
    lon = np.array([278,278])
    lat = np.array([439,438])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([278])
    lat = np.array([438])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 12-13h
    lon = np.array([278,278])
    lat = np.array([438,437])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 13-14h
    lon = np.array([278])
    lat = np.array([437])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 14-15h
    lon = np.array([278,277])
    lat = np.array([437,437])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 15-16h
    lon = np.array([277,277])
    lat = np.array([437,436])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 16-18h
    tim = np.arange(16,18)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 436
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 18-19h
    lon = np.array([277,276])
    lat = np.array([436,436])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 19-20h
    lon = np.array([276])
    lat = np.array([436])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 20-21h
    lon = np.array([276,276])
    lat = np.array([436,435])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 21-0h
    tim = np.arange(21,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 435
    tim_039 = np.append(tim_039, tim)
    lat_039 = np.append(lat_039, lat)
    lon_039 = np.append(lon_039, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_039, lat_039, lon_039, date)

    ###---------------------------------
    ### 4th September 2018
    ###---------------------------------

    date = '20180904'
    tim_049 = np.array([])
    lat_049 = np.array([])
    lon_049 = np.array([])

    ### box pick 0-1h
    lon = np.array([276,276])
    lat = np.array([435,434])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([276])
    lat = np.array([434])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([276,275])
    lat = np.array([434,434])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([275,275])
    lat = np.array([434,433])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 4-6h
    tim = np.arange(4,6)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 433
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-7h
    lon = np.array([275,275])
    lat = np.array([433,432])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-9h
    tim = np.arange(7,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 432
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([275,275])
    lat = np.array([432,431])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 10-11h
    lon = np.array([275])
    lat = np.array([431])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 11-12h
    lon = np.array([275,275])
    lat = np.array([431,430])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 12-13h
    lon = np.array([275,274])
    lat = np.array([430,430])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 13-14h
    lon = np.array([274])
    lat = np.array([430])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-15h
    lon = np.array([274,274])
    lat = np.array([430,429])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 15-17h
    tim = np.arange(15,17)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 429
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 17-18h
    lon = np.array([274,274])
    lat = np.array([429,428])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 18-21h
    tim = np.arange(18,21)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 428
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 21-22h
    lon = np.array([274,274])
    lat = np.array([428,427])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 18-21h
    tim = np.arange(22,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 427
    tim_049 = np.append(tim_049, tim)
    lat_049 = np.append(lat_049, lat)
    lon_049 = np.append(lon_049, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_049, lat_049, lon_049, date)

    ###---------------------------------
    ### 5th September 2018
    ###---------------------------------

    date = '20180905'
    tim_059 = np.array([])
    lat_059 = np.array([])
    lon_059 = np.array([])

    ### box pick 0-1h
    lon = np.array([274])
    lat = np.array([427])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([274,274,273])
    lat = np.array([427,426,426])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-11h
    tim = np.arange(2,11)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 426
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([273,273])
    lat = np.array([426,425])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 12-19h
    tim = np.arange(12,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 425
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 19-20h
    lon = np.array([273,274])
    lat = np.array([425,425])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-0h
    tim = np.arange(20,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 425
    tim_059 = np.append(tim_059, tim)
    lat_059 = np.append(lat_059, lat)
    lon_059 = np.append(lon_059, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_059, lat_059, lon_059, date)

    ###---------------------------------
    ### 6th September 2018
    ###---------------------------------

    date = '20180906'
    tim_069 = np.array([])
    lat_069 = np.array([])
    lon_069 = np.array([])

    ### box pick 0-1h
    lon = np.array([274,274])
    lat = np.array([425,424])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-14h
    tim = np.arange(1,14)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 424
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 14-15h
    lon = np.array([274,274])
    lat = np.array([424,423])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 15-16h
    lon = np.array([274])
    lat = np.array([423])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 16-17h
    lon = np.array([274,273])
    lat = np.array([423,423])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 17-20h
    tim = np.arange(17,20)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 423
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 20-21h
    lon = np.array([273,273])
    lat = np.array([423,422])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 21-0h
    tim = np.arange(21,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 422
    tim_069 = np.append(tim_069, tim)
    lat_069 = np.append(lat_069, lat)
    lon_069 = np.append(lon_069, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_069, lat_069, lon_069, date)

    ###---------------------------------
    ### 7th September 2018
    ###---------------------------------

    date = '20180907'
    tim_079 = np.array([])
    lat_079 = np.array([])
    lon_079 = np.array([])

    ### box pick 0-1h
    lon = np.array([273,273])
    lat = np.array([422,421])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-3h
    tim = np.arange(1,3)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 421
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 3-4h
    lon = np.array([273,273])
    lat = np.array([421,420])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 4-5h
    lon = np.array([273])
    lat = np.array([420])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 5-6h
    lon = np.array([273,273])
    lat = np.array([420,419])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-8h
    tim = np.arange(6,8)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 419
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 8-9h
    lon = np.array([273,273])
    lat = np.array([419,418])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([273])
    lat = np.array([418])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 10-11h
    lon = np.array([273,273])
    lat = np.array([418,417])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 11-13h
    tim = np.arange(11,13)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 417
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 13-14h
    lon = np.array([273,273])
    lat = np.array([417,416])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 14-15h
    lon = np.array([273])
    lat = np.array([416])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 15-16h
    lon = np.array([273,273])
    lat = np.array([416,415])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 16-17h
    lon = np.array([273,274])
    lat = np.array([415,415])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 17-18h
    lon = np.array([274,274])
    lat = np.array([415,414])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 18-19h
    lon = np.array([274])
    lat = np.array([414])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 19-20h
    lon = np.array([274,274])
    lat = np.array([414,413])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 20-22h
    tim = np.arange(20,22)
    lon = np.zeros([np.size(tim)])
    lon[:] = 274
    lat = np.zeros([np.size(tim)])
    lat[:] = 413
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 22-23h
    lon = np.array([274,274])
    lat = np.array([413,412])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 23-0h
    lon = np.array([274])
    lat = np.array([412])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_079 = np.append(tim_079, tim)
    lat_079 = np.append(lat_079, lat)
    lon_079 = np.append(lon_079, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_079, lat_079, lon_079, date)

    ###---------------------------------
    ### 8th September 2018
    ###---------------------------------

    date = '20180908'
    tim_089 = np.array([])
    lat_089 = np.array([])
    lon_089 = np.array([])

    ### box pick 0-1h
    lon = np.array([274,274])
    lat = np.array([412,411])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([274,275])
    lat = np.array([411,411])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([275])
    lat = np.array([411])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([275,275])
    lat = np.array([411,410])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 4-6h
    tim = np.arange(4,6)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 410
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 6-7h
    lon = np.array([275,275,276])
    lat = np.array([410,409,409])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 7-9h
    tim = np.arange(7,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 409
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 9-10h
    lon = np.array([276,276])
    lat = np.array([409,408])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 10-11h
    lon = np.array([276,277])
    lat = np.array([408,408])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 10-12h
    tim = np.arange(10,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 408
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 12-13h
    lon = np.array([277])
    lat = np.array([408])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 13-14h
    lon = np.array([277,277])
    lat = np.array([408,407])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 14-20h
    tim = np.arange(14,20)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 407
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 20-21h
    lon = np.array([277,277,278])
    lat = np.array([407,406,406])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 21-0h
    tim = np.arange(21,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 278
    lat = np.zeros([np.size(tim)])
    lat[:] = 406
    tim_089 = np.append(tim_089, tim)
    lat_089 = np.append(lat_089, lat)
    lon_089 = np.append(lon_089, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_089, lat_089, lon_089, date)

    ###---------------------------------
    ### 9th September 2018
    ###---------------------------------

    date = '20180909'
    tim_099 = np.array([])
    lat_099 = np.array([])
    lon_099 = np.array([])

    ### box pick 21-0h
    tim = np.arange(0,9)
    lon = np.zeros([np.size(tim)])
    lon[:] = 278
    lat = np.zeros([np.size(tim)])
    lat[:] = 406
    tim_099 = np.append(tim_099, tim)
    lat_099 = np.append(lat_099, lat)
    lon_099 = np.append(lon_099, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 9-10h
    lon = np.array([278,278])
    lat = np.array([406,405])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_099 = np.append(tim_099, tim)
    lat_099 = np.append(lat_099, lat)
    lon_099 = np.append(lon_099, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 21-0h
    tim = np.arange(10,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 278
    lat = np.zeros([np.size(tim)])
    lat[:] = 405
    tim_099 = np.append(tim_099, tim)
    lat_099 = np.append(lat_099, lat)
    lon_099 = np.append(lon_099, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_099, lat_099, lon_099, date)

    ###---------------------------------
    ### 10th September 2018
    ###---------------------------------

    date = '20180910'
    tim_109 = np.array([])
    lat_109 = np.array([])
    lon_109 = np.array([])

    ## box pick 0-12h
    tim = np.arange(0,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 278
    lat = np.zeros([np.size(tim)])
    lat[:] = 405
    tim_109 = np.append(tim_109, tim)
    lat_109 = np.append(lat_109, lat)
    lon_109 = np.append(lon_109, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 12-13h
    lon = np.array([278,278])
    lat = np.array([405,404])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_109 = np.append(tim_109, tim)
    lat_109 = np.append(lat_109, lat)
    lon_109 = np.append(lon_109, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 13-14h
    lon = np.array([278])
    lat = np.array([404])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_109 = np.append(tim_109, tim)
    lat_109 = np.append(lat_109, lat)
    lon_109 = np.append(lon_109, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 14-15h
    lon = np.array([278,277])
    lat = np.array([404,404])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_109 = np.append(tim_109, tim)
    lat_109 = np.append(lat_109, lat)
    lon_109 = np.append(lon_109, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ## box pick 15-0h
    tim = np.arange(15,24)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 404
    tim_109 = np.append(tim_109, tim)
    lat_109 = np.append(lat_109, lat)
    lon_109 = np.append(lon_109, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_109, lat_109, lon_109, date)

    ###---------------------------------
    ### 11th September 2018
    ###---------------------------------

    date = '20180911'
    tim_119 = np.array([])
    lat_119 = np.array([])
    lon_119 = np.array([])

    ## box pick 0-3h
    tim = np.arange(0,3)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 404
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 3-4h
    lon = np.array([277,277])
    lat = np.array([404,403])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ## box pick 4-6h
    tim = np.arange(4,6)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 403
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 6-7h
    lon = np.array([277,276])
    lat = np.array([403,403])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ## box pick 7-12h
    tim = np.arange(7,12)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 403
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 12-13h
    lon = np.array([276,276])
    lat = np.array([403,402])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ## box pick 13-16h
    tim = np.arange(13,16)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 402
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 16-17h
    lon = np.array([276,276])
    lat = np.array([402,401])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 17-18h
    lon = np.array([276])
    lat = np.array([401])
    tim = np.zeros([np.size(lon)])
    tim[:] = 17.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 18-19h
    lon = np.array([276,275])
    lat = np.array([401,401])
    tim = np.zeros([np.size(lon)])
    tim[:] = 18.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 19-20h
    lon = np.array([275])
    lat = np.array([401])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 20-21h
    lon = np.array([275,275])
    lat = np.array([401,400])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 21-22h
    lon = np.array([275])
    lat = np.array([400])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 22-23h
    lon = np.array([275,275])
    lat = np.array([400,399])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 23-0h
    lon = np.array([275])
    lat = np.array([399])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_119 = np.append(tim_119, tim)
    lat_119 = np.append(lat_119, lat)
    lon_119 = np.append(lon_119, lon)
    # for i in range(0,np.size(lon)):
    #     iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_119, lat_119, lon_119, date)

    ###---------------------------------
    ### 12th September 2018
    ###---------------------------------

    date = '20180912'
    tim_129 = np.array([])
    lat_129 = np.array([])
    lon_129 = np.array([])

    ### box pick 0-1h
    lon = np.array([275,275])
    lat = np.array([399,398])
    tim = np.zeros([np.size(lon)])
    tim[:] = 0.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 1-2h
    lon = np.array([275,275])
    lat = np.array([398,397])
    tim = np.zeros([np.size(lon)])
    tim[:] = 1.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 2-3h
    lon = np.array([275])
    lat = np.array([397])
    tim = np.zeros([np.size(lon)])
    tim[:] = 2.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 3-4h
    lon = np.array([275,275])
    lat = np.array([397,396])
    tim = np.zeros([np.size(lon)])
    tim[:] = 3.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 4-5h
    lon = np.array([275])
    lat = np.array([396])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 5-6h
    lon = np.array([275,275])
    lat = np.array([396,395])
    tim = np.zeros([np.size(lon)])
    tim[:] = 5.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 6-7h
    lon = np.array([275])
    lat = np.array([395])
    tim = np.zeros([np.size(lon)])
    tim[:] = 6.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 7-8h
    lon = np.array([275,275])
    lat = np.array([395,394])
    tim = np.zeros([np.size(lon)])
    tim[:] = 7.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 8-9h
    lon = np.array([275,275,274])
    lat = np.array([394,393,393])
    tim = np.zeros([np.size(lon)])
    tim[:] = 8.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 9-10h
    lon = np.array([274])
    lat = np.array([393])
    tim = np.zeros([np.size(lon)])
    tim[:] = 9.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 10-11h
    lon = np.array([274,274])
    lat = np.array([393,392])
    tim = np.zeros([np.size(lon)])
    tim[:] = 10.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([274,274])
    lat = np.array([392,391])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 12-13h
    lon = np.array([274])
    lat = np.array([391])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 13-14h
    lon = np.array([274,274])
    lat = np.array([391,390])
    tim = np.zeros([np.size(lon)])
    tim[:] = 13.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 14-15h
    lon = np.array([274,274])
    lat = np.array([390,389])
    tim = np.zeros([np.size(lon)])
    tim[:] = 14.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 15-16h
    lon = np.array([274])
    lat = np.array([389])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 16-17h
    lon = np.array([274,273,273])
    lat = np.array([389,389,388])
    tim = np.zeros([np.size(lon)])
    tim[:] = 16.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ## box pick 17-19h
    tim = np.arange(17,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 273
    lat = np.zeros([np.size(tim)])
    lat[:] = 388
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 19-20h
    lon = np.array([273,273,272])
    lat = np.array([388,387,387])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 20-21h
    lon = np.array([272])
    lat = np.array([387])
    tim = np.zeros([np.size(lon)])
    tim[:] = 20.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 21-22h
    lon = np.array([272,272])
    lat = np.array([387,388])
    tim = np.zeros([np.size(lon)])
    tim[:] = 21.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='black')

    ### box pick 22-23h
    lon = np.array([272,271])
    lat = np.array([388,388])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 23-0h
    lon = np.array([271])
    lat = np.array([388])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_129 = np.append(tim_129, tim)
    lat_129 = np.append(lat_129, lat)
    lon_129 = np.append(lon_129, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='magenta')

    # out = writeoutGrid(tim_129, lat_129, lon_129, date)

def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==17,data.values[:,1]==9),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==18,data.values[:,1]==9),data.values[:,3]==0))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print ('******')
    print ('')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print ('Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print ('Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print ('trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print ('')

    return trackShip_index

def readGriddedTrack(grid_filename):

    import pandas as pd

    print ('******')
    print ('')
    print ('Reading ' + grid_filename + ' file with pandas')
    print ('')

    data = pd.read_csv(grid_filename, sep = " ")
    values = data.values

    tim = values[:,1]
    ilon = values[:,2]
    ilat = values[:,3]

    return tim, ilat, ilon

def readGlobal(cube, ship_data, date_dir):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy

    print ('******')
    print ('')
    print ('Defining longitude and latitude boundaries:')
    print ('')

    if np.ndim(cube[0].data) == 3:      ### 2D data + time
        lats = cube[0].dim_coords[1].points
        lons = cube[0].dim_coords[2].points
    if np.ndim(cube[0].data) == 4:      ### 3D data + time
        lats = cube[0].dim_coords[2].points
        lons = cube[0].dim_coords[3].points

    ###---------------------------------------------------------------------------------
    ### find northern and southern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    nb_lats = lats + ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    sb_lats = lats - ((lats[1] - lats[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    print ('sb_lats.shape = ' + str(sb_lats.shape))

    ###---------------------------------------------------------------------------------
    ### find western and eastern boundaries of gridpoints
    ###---------------------------------------------------------------------------------
    wb_lons = lons - ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid
    eb_lons = lons + ((lons[1] - lons[0]) / 2.0)    ## use grid diff between 0 and 1 indices since uniform grid

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    #################################################################
    ## find date of interest
    #################################################################
    date = date_dir[0:8]
    if date == '20180831':
        date = '20180901'
    else:
        date = date[:6] + str(int(date[-2:]) + 1).zfill(2)
    day_ind = np.array([])
    day_ind = np.where(np.logical_and(ship_data.values[:,2] == float(date[-2:]),ship_data.values[:,1] == float(date[-4:-2])))
    print ('Daily ship track for ' + date + ': ' + str(len(day_ind[0])) + ' pts ')

    #################################################################
    ## save daily lat/lons as temp vars
    #################################################################
    ship_lats = ship_data.values[day_ind[0],7]
    ship_lons = ship_data.values[day_ind[0],6]

    print ('ship_lats.shape = ' + str(ship_lats.shape))

    #####--------------------------------------------------------------------------------------------------
    #####--------------------------------------------------------------------------------------------------
    ### compare hourly lat-lon with GLM grid
    data = {}
    data['ship_lons'] = np.zeros(24)
    data['ship_lats'] = np.zeros(24)
    data['ship_i'] = np.zeros(24); data['ship_i'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_j'] = np.zeros(24); data['ship_j'][:] = np.nan        ### set default ship_ind to nan so we can easily pick out out-of-grid values
    data['ship_hour'] = np.zeros(24)
    hours = np.arange(0,24)
    jflag = np.zeros(24)        ### flag for if grid boundary is crossed

    for h in hours:
        print ('')
        print ('hour = ' + str(h))
        for j in range(0,len(sb_lats)):     ### for all latitude points
            if np.logical_and(ship_lats[h] >= sb_lats[j], ship_lats[h] < nb_lats[j]):     ### find where ship lat sits on glm lat grid
                for i in range(0,len(wb_lons)):     ### for all longitude points
                    if np.logical_and(ship_lons[h] >= wb_lons[i], ship_lons[h] < eb_lons[i]):     ### find where ship lon sits on glm lon grid
                        print ('lats and lons match at i = ' + str(i) + ', j = ' + str(j))
                        jflag[h] = jflag[h] + 1
                        data['ship_lons'][h] = lons[i]
                        data['ship_hour'][h] = hours[h]
                        data['ship_lats'][h] = lats[j]
                        data['ship_j'][h] = j         # define grid point indices for use later
                        data['ship_i'][h] = i         # define grid point indices for use later

    # print data['ship_lats']
    # print data['ship_j']
    # print data['ship_lons']
    # print data['ship_i']

    # np.save('working_glm_grid', data)

    ### need to constract an hour, lon index, lat index list like used for the lam
    # for h in hours:
    #     if data['ship_i'][h+1] == data['ship_i'][h]:

    ### arguments to be passed back to pullTrack_CloudNet
    tim = hours
    ilon = data['ship_i']
    ilat = data['ship_j']

    print (tim)
    print (ilon)
    print (ilat)

    # #####--------------------------------------------------------------------------------------------------
    # #####--------------------------------------------------------------------------------------------------
    # ##################################################
    # ##################################################
    # #### 	CARTOPY MAP
    # ##################################################
    # ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=SMALL_SIZE)
    # plt.rc('ytick',labelsize=SMALL_SIZE)
    # plt.rc('legend',fontsize=SMALL_SIZE)
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    #
    # #################################################################
    # ## create figure and axes instances
    # #################################################################
    # plt.figure(figsize=(6,8))
    # ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))
    #
    # ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT
    #
    # #################################################################
    # ## add geographic features/guides for reference
    # #################################################################
    # ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    # ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    # ax.gridlines()
    #
    # #################################################################
    # ## plot global grid outline
    # #################################################################
    # ### draw outline of grid
    # # iplt.pcolormesh(cube[0][0,-10:-2,:-70])      ### covers whole drift
    # iplt.pcolormesh(cube[0][0,-7:-2,70:-70])      ### covers 28Aug - 4Sep subset of drift
    #
    # #################################################################
    # ## plot ship track
    # #################################################################
    # ### Plot tracks as line plot
    # plt.plot(ship_data.values[day_ind[0],6], ship_data.values[day_ind[0],7],
    #          color = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'Ship track',
    #          )
    # plt.plot(data['ship_lons'], data['ship_lats'],
    #          'o', color = 'yellow', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'Ship track',
    #          )
    #
    # plt.show()
    #
    # #####--------------------------------------------------------------------------------------------------

    return tim, ilat, ilon

def plot_cartmap(ship_data, cube, hour, grid_filename): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###################################
    ## CHOOSE DIAGNOSTIC
    ###################################
    diag = 9
    print ('')
    print ('Diag is: ', cube[diag].long_name)
    print ('')
    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    if len(cube[diag].dim_coords[-1].points) == 25:
    # if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif len(cube[diag].dim_coords[-1].points) == 56:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif len(cube[diag].dim_coords[-1].points) == 94:
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -211
        yoffset = -385
    elif len(cube[diag].dim_coords[-1].points) == 81:          ### 14th and 24th August
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -209
        yoffset = -399
    elif len(cube[diag].dim_coords[-1].points) == 380:         ### needs checked
    # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -60
        yoffset = -110
    else:
    # elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

    print ('Because cube shape = ', str(len(cube[diag].dim_coords[-1].points)))
    print ('xoffset = ', xoffset)
    print ('yoffset = ', yoffset)

    ### pcXXX
    # 0: total_radar_reflectivity / (unknown) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 1: air_pressure / (Pa)                 (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 2: air_temperature / (K)               (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 3: eastward_wind / (m s-1)             (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 4: large_scale_cloud_area_fraction / (1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 5: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 6: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 7: northward_wind / (m s-1)            (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 8: specific_humidity / (kg kg-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 9: upward_air_velocity / (m s-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting cartopy map:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
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

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.figure(figsize=(6,8))
    # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    # ax.set_extent([30, 60, 89.1, 89.6], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([40, 50, 88.4, 88.6], crs=ccrs.PlateCarree())       ### ZOOM
    ax.set_extent([0, 60, 87.75, 90], crs=ccrs.PlateCarree())     ### SWATH
    # ax.set_extent([-180, 190, 80, 90], crs=ccrs.PlateCarree())    ### WHOLE

    ### DON'T USE PLATECARREE, NORTHPOLARSTEREO (on it's own), LAMBERT

    #################################################################
    ## add geographic features/guides for reference
    #################################################################
    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    # ax.set_global()
    ax.gridlines()

    #################################################################
    ## plot UM data
    ################################################################
    # if np.size(cube[diag].data.shape) == 4:
    #     iplt.pcolormesh(cube[diag][hour,0,:,:])
    # elif np.size(cube[diag].data.shape) == 3:
    #     iplt.pcolormesh(cube[diag][hour,:,:])
    #     # iplt.pcolormesh(cube[hour,471:495,240:264])
    # elif np.size(cube[diag].data.shape) == 2:
    iplt.pcolormesh(cube[diag][hour,:,:])
    plt.title(cube[diag].standard_name)# + ', ' + str(cube[diag].units))
    plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    qplt.outline(cube[diag][hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[diag][hour,386:495,211:305])          ### misc
    # qplt.outline(cube[diag][0,:,:])          ### 14th/24th swath
    # qplt.outline(cube[diag][hour,:,:])

    # gridship = gridShipTrack(cube[diag], xoffset, yoffset)

            #### MID POINT: (433, 258)

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)
    trackShip_index = trackShip(ship_data)

    ### Plot tracks as line plot
    plt.plot(ship_data.values[trackShip_index,6], ship_data.values[trackShip_index,7],
             color = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(), label = 'Ship track',
             )
    plt.plot(ship_data.values[trackShip_index[0],6], ship_data.values[trackShip_index[0],7],
             'k^', markerfacecolor = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )
    plt.plot(ship_data.values[trackShip_index[-1],6], ship_data.values[trackShip_index[-1],7],
             'kv', markerfacecolor = 'darkorange', linewidth = 3,
             transform = ccrs.PlateCarree(),
             )

    #################################################################
    ## read in and plot gridded ship track
    #################################################################
    tim, ilat, ilon = readGriddedTrack(grid_filename)

    ### Plot tracks as line plot
    for i in range(0, len(ilon)-1):
        iplt.scatter(cube[diag].dim_coords[2][int(ilon[i] + xoffset)], cube[diag].dim_coords[1][int(ilat[i] + yoffset)],color='black')


    ### Plot tracks as line plot
    # plt.plot(ship_data.values[:,6], ship_data.values[:,7],
    #          color = 'yellow', linewidth = 2,
    #          transform = ccrs.PlateCarree(), label = 'Whole',
    #          )
    # plt.plot(ship_data.values[inIce_index,6], ship_data.values[inIce_index,7],
    #          color = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(), label = 'In Ice',
    #          )
    # plt.plot(ship_data.values[inIce_index[0],6], ship_data.values[inIce_index[0],7],
    #          'k^', markerfacecolor = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(),
    #          )
    # plt.plot(ship_data.values[inIce_index[-1],6], ship_data.values[inIce_index[-1],7],
    #          'kv', markerfacecolor = 'darkorange', linewidth = 3,
    #          transform = ccrs.PlateCarree(),
    #          )
    # plt.plot(ship_data.values[drift_index,6], ship_data.values[drift_index,7],
    #          color = 'red', linewidth = 4,
    #          transform = ccrs.PlateCarree(), label = 'Drift',
    #          )


    #### test plotting of unrotated grid
    # lon, lat = unrotateGrid(cube)

    # plt.plot(np.nanmin(lon),np.nanmin(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmin(lon),np.nanmax(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmax(lon),np.nanmin(lat),
    #         color='black',transform = ccrs.PlateCarree())
    # plt.plot(np.nanmax(lon),np.nanmax(lat),
    #         color='black',transform = ccrs.PlateCarree())

    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting cartopy map! :)')
    print ('')

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

def excludeZeros(cube):

    print ('')
    print ('Checking stash list:')

    print ('Want to exclude zeros in the following fields:')
    ### list of stash items where we want to exclude zeros
    STASH = ['m01s00i012','m01s00i254','m01s00i075','m01s00i076','m01s00i078',
        'm01s00i079','m01s00i081','m01s00i271','m01s00i273','m01s00i272',
        'm01s00i083','m01s00i084','m01s00i088']
    print (STASH)

    print ('Diag is:')
    str_m = "%02d" % cube.attributes['STASH'][0]
    str_s = "%02d" % cube.attributes['STASH'][1]
    str_i = "%03d" % cube.attributes['STASH'][2]
    stash = str('m' + str_m + 's' + str_s + 'i' + str_i)
    print (stash)

    for i in range(0, len(STASH)):
        if STASH[i] == stash:
            # data[data==0] = np.nan              # set zeros to nans
            flag = 1                           # flagging if in list
            print ('In list, so excluding zeros')
            break
        else:
            flag = 0                           # flagging if not in list
            print ('Not in list, so not excluding zeros')

    # if flag == 1:
    # if flag == 0:
    # print ''

    return flag, stash

def checkWind(cube):

    # 0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

    print ('')
    print ('Checking stash list:')

    print ('Want to change gridding of the following fields:')
    ### list of stash items where we want to change gridding (u/v wind)
    STASH = ['m01s00i002','m01s00i003']
    print (STASH)

    print ('Diag is:')
    str_m = "%02d" % cube.attributes['STASH'][0]
    str_s = "%02d" % cube.attributes['STASH'][1]
    str_i = "%03d" % cube.attributes['STASH'][2]
    stash = str('m' + str_m + 's' + str_s + 'i' + str_i)
    print (stash)

    for i in range(0, len(STASH)):
        if stash == STASH[i]:
            flag = 1                           # flagging if in list
        else:
            flag = 0                           # flagging if not in list

    tempvar = np.zeros(cube.shape)
    if flag == 1:
        print ('In list, so changing vertical grid')
        tempvar[:,:-1,:,:] = 0.5*(cube.data[:,:-1,:,:] + cube.data[:,1:,:,:])
        print ('Cube = ' + str(cube.data[0,0:9,10,10]))
        print ('Tempvar = ' + str(tempvar[0,0:9,10,10]))
    if flag == 0:
        print ('Not in list, so not changing vertical grid')

    tempvar[:,-1,:,:] = np.nan

    cube.data = tempvar

    return cube, stash

def fixHeight(data, cube):

    print ('******')
    print ('')
    print ('Adjusting height to common vertical grid...')
    print ('')

    # height = cube[1].aux_coords[2].points.data       ### 71 levels

    ### wind fields have Z[0] == 2.5
    ### all other 4D fields have Z[0] >= 5.0

    if np.round(cube.aux_coords[2][0].points) > 3:
    # if np.round(cube.aux_coords[2][0].points) == 5:
        ### making 70 levels into 71 for common grid
        cubedata = np.zeros([71,24])
        cubedata[1:,:] = data
        cubedata[0,:] = np.nan
    elif np.round(cube.aux_coords[2][0].points) == 2:
        ### interpolating to n71 common grid
        ### upper bounds = cube[8].aux_coords[2].bounds[:,1]
        cubedata = np.zeros([71,24])
        for i in range(0,24):
            temp = np.interp(cube.aux_coords[2].bounds[:,1],cube.aux_coords[2].points,data[:,i])
            cubedata[1:,i] = temp
            cubedata[0,i] = np.nan
    else:
        cubedata = data

    return cubedata

def pullSwath_CloudNet(cube, grid_filename, con, stream, date, model, ship_data, nc_outfile):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt
    import pandas as pd

    print ('******')
    print ('')
    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    if len(cube[0].dim_coords[-1].points) == 25:
        xoffset = -239
        yoffset = -470
    elif len(cube[0].dim_coords[-1].points) == 56:
        xoffset = -210
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 94:
        xoffset = -211
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 81:          ### 14th and 24th August
        xoffset = -209
        yoffset = -399
    elif len(cube[0].dim_coords[-1].points) == 380:         ### needs checked
        xoffset = -60
        yoffset = -110
    else:
        ### if glm
        xoffset = 0
        yoffset = 0

    print ('Because cube shape = ', str(len(cube[0].dim_coords[-1].points)))
    print ('xoffset = ', xoffset)
    print ('yoffset = ', yoffset)

    #################################################################
    ## load gridded ship track
    #################################################################
    # print '******'
    print ('')
    print ('Pulling gridded swath from cube:')
    print ('')

    # if model == 'lam':
    #     tim, ilat, ilon = readGriddedTrack(grid_filename)
    # elif model == 'glm':
    #     tim, ilat, ilon = readGlobal(cube, ship_data, date)
    # else:
    #     print ('Model option is not valid')

    #################################################################
    ## fix time index
    #################################################################

    if np.size(cube)>1:
        print ('')
        print ('More than one variable constraint. Proceeding...')
        print ('')
        print (np.size(cube))

        #################################################################
        ## CREATE EMPTY CUBE FOR PC COLUMN DIAGNOSTICS
        #################################################################
        ### swath definition based on date, change of shape at 2 Sept
        if date[5] == 8:
            ncube = Cube(np.zeros([np.size(cube),70,24,26,50]))
            lat0 = 70
            lat1 = -1
            lon0 = 20
            lon1 = 70
        elif np.logical_and(date[5] == 9, int(date[6:8]) <= 2):
            ncube = Cube(np.zeros([np.size(cube),70,24,26,50]))
            lat0 = 70
            lat1 = -1
            lon0 = 20
            lon1 = 70
        else:
            ncube = Cube(np.zeros([np.size(cube),70,24,72,20]))
            lat0 = 0
            lat1 = 72
            lon0 = 50
            lon1 = 70

        print (date)
        print (ncube.shape)
        print (np.size(ncube),2)

        #################################################################
        ## POPULATE NP ARRAY WITH DATA
        #################################################################
        for k in range(0,np.size(cube)):            ### loop over number of variables
            print ('')
            print ('k = ', k) ###', so processing', con[k]   # doesn't work with global_con
            print ('')
            ################################################################
            # only consider hourly diagnostics
            ################################################################
            if len(np.round(cube[k].coord('forecast_period').points)) > 10:
                ###---------------------------------
                ### CHECK IF OFFSETS NEED TO BE RE-ADJUSTED
                ###---------------------------------
                print ('Double-checking grid:')
                if len(cube[k].dim_coords[-1].points) == 25:
                # if cube[0,0].shape >= 25-1:    # ll = 240, 471
                    xoffset = -239
                    yoffset = -470
                elif len(cube[k].dim_coords[-1].points) == 56:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -210
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 94:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -211
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 81:          ### 14th and 24th August
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -209
                    yoffset = -399
                elif len(cube[k].dim_coords[-1].points) == 380:         ### needs checked
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -60
                    yoffset = -110
                else:
                # elif cube[0,0].shape >= 500-1:
                    xoffset = 0
                    yoffset = 0
                print ('Offsets are: ' + str(xoffset) + ', ' + str(yoffset))

                #################################################################
                ## make hourly time array
                #################################################################
                print ('cubetime = ' + str(np.size(cube[k].coord('forecast_period').points)))
                # if np.size(cube[k].coord('forecast_period').points) == 24:
                cubetime = np.round(cube[k].coord('forecast_period').points - 12.0)      ### 24h forecast period (ignore first 12h)
                # else:
                #     cubetime = np.round(cube[k].coord('forecast_period').points[:-1] - 12.0)      ### 25h forecast period (ignore first 12h, exclude 36h)

                        #### will need to add code for dealing with forecast_period
                        ####    i.e. 1.00016 hrs past 0000 UTC
                        ####     will need to concatenate from next day's file?

                print ('')
                print ('Cube times relative to forecast start:', cubetime[:])
                print ('')

                #################################################################
                ## PROBE VARIABLE
                #################################################################
                ### do we want to average exluding zeros?
                stash_flag, stash = excludeZeros(cube[k])

                ### do we need to re-grid?  -- DOESN'T WORK LIKE WRF, GRID NOT SPACED SAME WAY
                # cube[k], wind_stash = checkWind(cube[k])

                #################################################################
                ## CHECK DIMENSIONS
                #################################################################
                if np.logical_and(np.size(cube[k].data,1) > 68, np.size(cube[k].data,1) < 72):
                    print ('Variable is 4D:')
                    print ('')
                    #### create empty arrays to be filled
                    data = np.zeros([len(cube[k].coord('model_level_number').points),len(cubetime)-1, np.size(ncube,3), np.size(ncube,4)])
                    ### make dimension flag
                    dim_flag = 1        ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')
                else:
                    print ('Variable is 3D:')
                    print ('')
                    #### create empty arrays to be filled
                    if stream[1:3] == 'pb':
                        if cube[k].long_name == 'large_scale_ice_water_path':
                            data = np.zeros([len(cubetime), np.size(ncube,3), np.size(ncube,4)])
                        elif cube[k].long_name == 'large_scale_liquid_water_path':
                            data = np.zeros([len(cubetime), np.size(ncube,3), np.size(ncube,4)])
                        else:
                            data = np.zeros([len(cubetime)-1, np.size(ncube,3), np.size(ncube,4)])
                    elif stream[1:3] == 'pa':
                        if len(cubetime) == 25:
                            data = np.zeros([len(cubetime)-1, np.size(ncube,3), np.size(ncube,4)])
                        elif len(cubetime) == 24:
                            data = np.zeros([len(cubetime), np.size(ncube,3), np.size(ncube,4)])
                    elif stream[1:3] == 'pd':
                        data = np.zeros([len(cubetime)-1, np.size(ncube,3), np.size(ncube,4)])
                    dim_flag = 0       ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')

                #################################################################
                ## LOOP OVER TIME INDEX, DECOMPOSE ONTO 24H TIMESERIES
                #################################################################

                varname = varnames.findfieldName(stash)
                print ('standard_name = ', cube[k].standard_name)
                print ('long name = ', cube[k].long_name)
                print ('varname = ', varname)
                print ('')

                # for j in range(0,len(cubetime)-1):              ### loop over time
                #     if j < len(cubetime[:-1]):
                #         itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
                #     else:
                #         ### end point (23h)
                #         itime = np.where(tim >= cubetime[-1])
                #     # print ''
                #     print ('For ', str(j), 'h, itime = ', itime)

                if dim_flag == 1: dat = np.zeros([len(cube[k].coord('model_level_number').points), len(cubetime), np.size(ncube,3), np.size(ncube,4)])
                if dim_flag == 0: dat = np.zeros([len(cubetime), np.size(ncube,3), np.size(ncube,4)])
                for i in range(0, len(cubetime)):                   ### loop over time gridded by ship track
                    if dim_flag == 1:
                        temp = cube[k][i,:,lat0:lat1,lon0:lon1]
                        dat[:,i,:,:] = np.squeeze(temp.data)
                    if dim_flag == 0:
                        temp = cube[k][i,lat0:lat1,lon0:lon1]
                        dat[i,:,:] = np.squeeze(temp.data)
                    # if np.size(itime) > 1:
                    #     if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                    #     if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                    #     if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                    #     # print 'averaging over itime ...'
                    #     # print ''
                    # else:
                    #     if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                    #     if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                    #     # print 'no averaging, itime = 1 ...'
                    #     print ('')
                # print data
            print ('dat.shape = ', dat.shape)
    #
    #         #################################################################
    #         ## CREATE CUBE
    #         #################################################################
    #         ### ECMWF FIELD NAMES
    #         # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
    #         #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
    #         #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};

    #
    #             if stream[1:3] == 'pa':
    #                 a = len(cube[k].aux_coords)
    #                 for ft in range(0,a):
    #                     print(cube[k].aux_coords[ft].standard_name)
    #                     if cube[k].aux_coords[ft].standard_name == 'forecast_period':
    #                         if np.size(cube[k].aux_coords[ft].points) > 24:          ## accounts for arrays with 25 timesteps (includes t12 and t36)
    #                             ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
    #                         else:
    #                             ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
    #             else:
    #                 if cube[k].long_name == 'large_scale_ice_water_path':
    #                     ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
    #                 elif cube[k].long_name == 'large_scale_liquid_water_path':
    #                     ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
    #                 else:
    #                     ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
    #             print (ntime.shape)
    #             if dim_flag == 1:         ### 4D VARIABLE
    #                 if stream[1:3] == 'pd':
    #                     model_height = DimCoord(cube[k].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
    #                     comdata = data                    #### leave BL diagnostics on RHO levels
    #                 else:
    #                     model_height = DimCoord(cube[1].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
    #                     comdata = fixHeight(data, cube[k])
    #                 ncube = Cube(np.transpose(comdata),
    #                         dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
    #                         standard_name = cube[k].standard_name,
    #                         long_name = cube[k].long_name,
    #                         units = cube[k].units,
    #                         var_name = varname,
    #                         attributes = cube[k].attributes,
    #                         aux_coords_and_dims = None,
    #                         )
    #             elif dim_flag == 0:         ### 3D VARIABLE
    #                 ncube = Cube(np.transpose(data),
    #                         dim_coords_and_dims=[(ntime, 0)],
    #                         standard_name = cube[k].standard_name,
    #                         long_name = cube[k].long_name,
    #                         units = cube[k].units,
    #                         var_name = varname,
    #                         attributes = cube[k].attributes,
    #                         aux_coords_and_dims = None,
    #                         )
    #             # ncube.attributes = cube[k].attributes
    #             # iris.save(ncube, pp_outfile, append=True)
    #             if k == 0:
    #                 print ('Initialising fcube')
    #                 print ('')
    #                 fcube = [ncube]
    #             else:
    #                 print ('Appending variable to fcube')
    #                 print ('')
    #                 fcube.append(ncube)
    #
    #     # print fcube
    #
    # #################################################################
    # ## define output filename
    # #################################################################
    # # print 'Define pp stream outfile:'
    # # pp_outfile = date[:6] + str(int(date[6:8])+1) + '_oden_metum_' + str(stream[2:3]) + '.pp'
    # # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
    # ### bespoke setup if dir is 20180831T1200Z (for 20180901 data)
    # # if date == '20180831T1200Z': nc_outfile = '20180901_oden_metum.nc'
    # # print 'Outfile = ', pp_outfile
    #
    # ### save cube to netcdf file
    # print ('')
    # print ('Writing fcube to file:')
    # print ('')
    # if stream[1:3] == 'pc':
    #     ## Combine track-pulled pp output files to one netCDF
    #     ## First, make netCDF with pc stream (using Iris cubes)
    #     print ('fcube = ')
    #     print (fcube)
    #     print ('')
    #     print ('******')
    #     print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
    #     print ('')
    #     if not os.path.exists(nc_outfile):
    #         if 'fcube' in locals():
    #             out = writeNetCDF(date, fcube, nc_outfile)
    #         # if PC outfile already exists, combine other stream data
    #         # if PC outfile doesn't exist, write new
    #
    # if stream[1:3] == 'pd':
    #     ## Combine track-pulled pp output files to one netCDF
    #     ## First, make netCDF with pd stream (using Iris cubes)
    #     print ('fcube = ')
    #     print (fcube)
    #     print ('')
    #     print ('******')
    #     print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
    #     print ('***file is merged to outfile later***')
    #     print ('')
    #     doutfile = nc_outfile[:-3] + '_d.nc'
    #     if not os.path.exists(doutfile):
    #         if 'fcube' in locals():
    #             out = writePD_BL(fcube, doutfile)
    #         # if PD outfile already exists, combine other stream data
    #         # if PD outfile doesn't exist, write new
    #
    # if stream[1:3] == 'pe':
    #     ## Combine track-pulled pp output files to one netCDF
    #     ## First, make netCDF with pd stream (using Iris cubes)
    #     print ('fcube = ')
    #     print (fcube)
    #     print ('')
    #     print ('******')
    #     print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
    #     print ('***file is merged to outfile later***')
    #     print ('')
    #     eoutfile = nc_outfile[:-3] + '_e.nc'
    #     if not os.path.exists(eoutfile):
    #         if 'fcube' in locals():
    #             out = writeFile_netCDF4(fcube, eoutfile)
    #         # if PC outfile already exists, combine other stream data
    #         # if PC outfile doesn't exist, write new
    #
    # elif stream[1:3] == 'pb':
    #     print ('fcube = ')
    #     print (fcube)
    #     print ('')
    #     print ('******')
    #     print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
    #     print ('***file is merged to outfile later***')
    #     print ('')
    #     ## Next, append 1D timeseries (surface) data (pb stream)
    #     ## Can't use Iris for this as cubes can't be 1D
    #     ##              -> uses standard netCDF appending function
    #     boutfile = nc_outfile[:-3] + '_b.nc'
    #     if not os.path.exists(boutfile):
    #         if 'fcube' in locals():
    #             out = writePB_Cloudnet(fcube, boutfile)     ##!!!! NEEDS UPDATING TO ONLY WRITE VARIABLES IN FILE, NOT HARD CODED
    #
    # elif stream[1:3] == 'pa':
    #     print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
    #     print ('***file is merged to outfile later***')
    #     print ('')
    #     ## Next, append 1D timeseries (surface) data (pb stream)
    #     ## Can't use Iris for this as cubes can't be 1D
    #     ##              -> uses standard netCDF appending function
    #     aoutfile = nc_outfile[:-3] + '_a.nc'
    #     if not os.path.exists(aoutfile):
    #         if 'fcube' in locals():
    #             out = writePA_Analysis(fcube, aoutfile)
    #
    # return nc_outfile

def pullTrack_CloudNet(cube, grid_filename, con, stream, date, model, ship_data, nc_outfile):

    from iris.coords import DimCoord
    from iris.cube import Cube
    import iris.plot as iplt
    import pandas as pd

    print ('******')
    print ('')
    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    print ('What grid are we looking at?')
    if len(cube[0].dim_coords[-1].points) == 25:
        xoffset = -239
        yoffset = -470
    elif len(cube[0].dim_coords[-1].points) == 56:
        xoffset = -210
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 94:
        xoffset = -211
        yoffset = -385
    elif len(cube[0].dim_coords[-1].points) == 81:          ### 14th and 24th August
        xoffset = -209
        yoffset = -399
    elif len(cube[0].dim_coords[-1].points) == 380:         ### needs checked
        xoffset = -60
        yoffset = -110
    else:
        ### if glm
        xoffset = 0
        yoffset = 0

    print ('Because cube shape = ', str(len(cube[0].dim_coords[-1].points)))
    print ('xoffset = ', xoffset)
    print ('yoffset = ', yoffset)

    #################################################################
    ## load gridded ship track
    #################################################################
    # print '******'
    print ('')
    print ('Pulling gridded track from cube:')
    print ('')

    if model == 'lam':
        tim, ilat, ilon = readGriddedTrack(grid_filename)
    elif model == 'glm':
        tim, ilat, ilon = readGlobal(cube, ship_data, date)
    else:
        print ('Model option is not valid')

    #################################################################
    ## fix time index
    #################################################################

    if np.size(cube)>1:
        print ('')
        print ('More than one variable constraint. Proceeding...')
        print ('')
        print (np.size(cube))

        #################################################################
        ## CREATE EMPTY CUBE FOR PC COLUMN DIAGNOSTICS
        #################################################################
        ncube = Cube(np.zeros([np.size(cube),70,24]))

        #################################################################
        ## POPULATE NP ARRAY WITH DATA
        #################################################################
        ### populate 0th dimension with time field
        # data[:,0] = cubetime[:,:-1]

        # #################################################################
        # ## if our diagnostics are 3-hourly, ignore
        # #################################################################
        # if len(np.round(cube[k].coord('forecast_period').points)) <= 10:
        #     print cube[k].standard_name
        #     print 'Diagnostic is 3-hourly, break from loop'
        #     break
        # # elif len(np.round(cube[k].coord('forecast_period').points)) > 10:
        # #     if xoffset == 0:
        # #         print cube[k].standard_name
        # #         print 'Diagnostic is 1-hourly, BUT this is a STASH typo since diagnostic covers whole nest.'
        # #         ok = False
        # else:
        #     print cube[k].standard_name
        #     # if int(xoffset) != 0:
        #     print 'Diagnostic is 1-hourly, pull ship track...'
        #     ok = True
        #
        # if not ok: continue

        for k in range(0,np.size(cube)):            ### loop over number of variables
            print ('')
            print ('k = ', k) ###', so processing', con[k]   # doesn't work with global_con
            print ('')
            #################################################################
            ## only consider hourly diagnostics
            #################################################################
            if len(np.round(cube[k].coord('forecast_period').points)) > 10:
                ###---------------------------------
                ### CHECK IF OFFSETS NEED TO BE RE-ADJUSTED
                ###---------------------------------
                print ('Double-checking grid:')
                if len(cube[k].dim_coords[-1].points) == 25:
                # if cube[0,0].shape >= 25-1:    # ll = 240, 471
                    xoffset = -239
                    yoffset = -470
                elif len(cube[k].dim_coords[-1].points) == 56:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -210
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 94:
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -211
                    yoffset = -385
                elif len(cube[k].dim_coords[-1].points) == 81:          ### 14th and 24th August
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -209
                    yoffset = -399
                elif len(cube[k].dim_coords[-1].points) == 380:         ### needs checked
                # elif cube[0,0].shape >= 93-1:    # ll = 211, 386
                    xoffset = -60
                    yoffset = -110
                else:
                # elif cube[0,0].shape >= 500-1:
                    xoffset = 0
                    yoffset = 0
                print ('Offsets are: ' + str(xoffset) + ', ' + str(yoffset))

                #################################################################
                ## make hourly time array
                #################################################################
                print ('cubetime = ' + str(np.size(cube[k].coord('forecast_period').points)))
                # if np.size(cube[k].coord('forecast_period').points) == 24:
                cubetime = np.round(cube[k].coord('forecast_period').points - 12.0)      ### 24h forecast period (ignore first 12h)
                # else:
                #     cubetime = np.round(cube[k].coord('forecast_period').points[:-1] - 12.0)      ### 25h forecast period (ignore first 12h, exclude 36h)

                        #### will need to add code for dealing with forecast_period
                        ####    i.e. 1.00016 hrs past 0000 UTC
                        ####     will need to concatenate from next day's file?

                # print ''
                # print 'Cube times relative to forecast start:', cubetime[:-1]
                # print ''

                #################################################################
                ## PROBE VARIABLE
                #################################################################
                ### do we want to average exluding zeros?
                stash_flag, stash = excludeZeros(cube[k])

                ### do we need to re-grid?  -- DOESN'T WORK LIKE WRF, GRID NOT SPACED SAME WAY
                # cube[k], wind_stash = checkWind(cube[k])

                #################################################################
                ## CHECK DIMENSIONS
                #################################################################
                if np.logical_and(np.size(cube[k].data,1) > 68, np.size(cube[k].data,1) < 72):
                    print ('Variable is 4D:')
                    print ('')
                    #### create empty arrays to be filled
                    data = np.zeros([len(cube[k].coord('model_level_number').points),len(cubetime)-1])
                    ### make dimension flag
                    dim_flag = 1        ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')
                else:
                    print ('Variable is 3D:')
                    print ('')
                    #### create empty arrays to be filled
                    if stream[1:3] == 'pb':
                        if cube[k].long_name == 'large_scale_ice_water_path':
                            data = np.zeros([len(cubetime)])
                        elif cube[k].long_name == 'large_scale_liquid_water_path':
                            data = np.zeros([len(cubetime)])
                        else:
                            data = np.zeros([len(cubetime)-1])
                    elif stream[1:3] == 'pa':
                        if len(cubetime) == 25:
                            data = np.zeros([len(cubetime)-1])
                        elif len(cubetime) == 24:
                            data = np.zeros([len(cubetime)])
                    elif stream[1:3] == 'pd':
                        data = np.zeros([len(cubetime)-1])
                    dim_flag = 0       ### for next loops
                    print ('data.shape = ', str(data.shape))
                    print ('')

                #################################################################
                ## LOOP OVER TIME INDEX, DECOMPOSE ONTO 24H TIMESERIES
                #################################################################
                for j in range(0,len(cubetime)-1):              ### loop over time
                    if j < len(cubetime[:-1]):
                        itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
                    else:
                        ### end point (23h)
                        itime = np.where(tim >= cubetime[-1])
                    # print ''
                    print ('For ', str(j), 'h, itime = ', itime)
                    if dim_flag == 1: dat = np.zeros([len(cube[k].coord('model_level_number').points),len(itime[0])])
                    if dim_flag == 0: dat = np.zeros([len(itime[0])])
                    for i in range(0, len(itime[0])):                   ### loop over time gridded by ship track
                        if np.size(itime) > 1:
                            # print 'Processing i = ', str(itime[0][i])
                            # print '...'
                            if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                            if dim_flag == 0: temp = cube[k][j,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                        else:
                            # print 'Processing i = ', str(itime[i])
                            # print '...'
                            if dim_flag == 1: temp = cube[k][j,:,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                            if dim_flag == 0: temp = cube[k][j,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                        if dim_flag == 1: dat[:,i] = np.squeeze(temp.data)
                        if dim_flag == 0: dat[i] = np.squeeze(temp.data)
                        if np.size(itime) > 1:
                            if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                            if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                            if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                            # print 'averaging over itime ...'
                            # print ''
                        else:
                            if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                            if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                            # print 'no averaging, itime = 1 ...'
                            print ('')
                    # print data
            # print 'data.shape = ', data.shape

            #################################################################
            ## FIGURES TO TEST OUTPUT
            #################################################################
            ### timeseries of lowest model level
            # plt.figure(figsize=(7,5))
            # plt.plot(cubetime[:-1],data[0:10,:])
            # plt.show()

            ### vertical profile of 1st timestep
            # plt.figure(figsize=(7,5))
            # plt.plot(data[:,0],cube.coord('model_level_number').points)
            # plt.show()

            ### pcolormesh of timeseries
            # plt.figure(figsize=(7,5))
            # plt.pcolormesh(cubetime[:-1], cube.coord('model_level_number').points, data)
            # plt.colorbar()
            # plt.show()

            #################################################################
            ## CREATE CUBE
            #################################################################
            ### ECMWF FIELD NAMES
            # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
            #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
            #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};
                varname = varnames.findfieldName(stash)
                print ('standard_name = ', cube[k].standard_name)
                print ('long name = ', cube[k].long_name)
                print ('varname = ', varname)
                print ('')

                if stream[1:3] == 'pa':
                    a = len(cube[k].aux_coords)
                    for ft in range(0,a):
                        print(cube[k].aux_coords[ft].standard_name)
                        if cube[k].aux_coords[ft].standard_name == 'forecast_period':
                            if np.size(cube[k].aux_coords[ft].points) > 24:          ## accounts for arrays with 25 timesteps (includes t12 and t36)
                                ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                            else:
                                ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                else:
                    if cube[k].long_name == 'large_scale_ice_water_path':
                        ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                    elif cube[k].long_name == 'large_scale_liquid_water_path':
                        ntime = DimCoord(cubetime[:], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                    else:
                        ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
                print (ntime.shape)
                if dim_flag == 1:         ### 4D VARIABLE
                    if stream[1:3] == 'pd':
                        model_height = DimCoord(cube[k].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
                        comdata = data                    #### leave BL diagnostics on RHO levels
                    else:
                        model_height = DimCoord(cube[1].aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
                        comdata = fixHeight(data, cube[k])
                    ncube = Cube(np.transpose(comdata),
                            dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
                            standard_name = cube[k].standard_name,
                            long_name = cube[k].long_name,
                            units = cube[k].units,
                            var_name = varname,
                            attributes = cube[k].attributes,
                            aux_coords_and_dims = None,
                            )
                elif dim_flag == 0:         ### 3D VARIABLE
                    ncube = Cube(np.transpose(data),
                            dim_coords_and_dims=[(ntime, 0)],
                            standard_name = cube[k].standard_name,
                            long_name = cube[k].long_name,
                            units = cube[k].units,
                            var_name = varname,
                            attributes = cube[k].attributes,
                            aux_coords_and_dims = None,
                            )
                # ncube.attributes = cube[k].attributes
                # iris.save(ncube, pp_outfile, append=True)
                if k == 0:
                    print ('Initialising fcube')
                    print ('')
                    fcube = [ncube]
                else:
                    print ('Appending variable to fcube')
                    print ('')
                    fcube.append(ncube)

        # print fcube

    else:
        print ('')
        print ('Only one variable constraint. Proceeding...')
        print ('')

        cubetime = np.round(cube.coord('forecast_period').points - 12.0)      ### forecast period (ignore first 12h)
        print ('')
        print ('Cube times relative to forecast start (excluding first 12H):', cubetime[:-1])
        print ('')

        #################################################################
        ## CREATE EMPTY CUBE
        #################################################################
        ncube = Cube(np.zeros([len(cube.coord('model_level_number').points),len(cubetime)-1]))

        #################################################################
        ## PROBE VARIABLE
        #################################################################
        ### do we want to average exluding zeros?
        stash_flag, stash = excludeZeros(cube)

        #################################################################
        ## FIND ARRAY SIZE AND CREATE EMPTY NP ARRAY
        #################################################################
        if np.logical_and(np.size(cube.data,1) >= 69, np.size(cube.data,1) < 71):
            print ('Variable is 4D:')
            print ('')
            #### create empty arrays to be filled
            data = np.zeros([len(cube.coord('model_level_number').points),len(cubetime)-1])
            dim_flag = 1        ### for next loops
            print ('data.shape = ', str(data.shape))
            print ('')
        else:
            print ('Variable is 3D:')
            print ('')
            #### create empty arrays to be filled
            data = np.zeros([len(cubetime)-1])
            dim_flag = 0       ### for next loops
            print ('data.shape = ', str(data.shape))
            print ('')

        #################################################################
        ## POPULATE NP ARRAY WITH DATA
        #################################################################
        ### populate 0th dimension with time field
        # data[:,0] = cubetime[:,:-1]

        for j in range(0,len(cubetime)-1):
            if j < len(cubetime[:-1]):
                itime = np.where(np.logical_and(tim >= cubetime[j], tim < cubetime[j+1]))
            else:
                ### end point (23h)
                itime = np.where(tim >= cubetime[-1])
            print ('For ', str(j), 'h, itime = ', itime)
            if dim_flag == 1: dat = np.zeros([len(cube.coord('model_level_number').points),len(itime[0])])
            if dim_flag == 0: dat = np.zeros([len(itime[0])])
            for i in range(0, len(itime[0])):
                if np.size(itime) > 1:
                    # print 'Processing i = ', str(itime[0][i])
                    if dim_flag == 1: temp = cube[j,:,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                    if dim_flag == 0: temp = cube[j,int(ilat[itime[0][i]] + yoffset),int(ilon[itime[0][i]] + xoffset)]
                else:
                    # print 'Processing i = ', str(itime[i])
                    if dim_flag == 1: temp = cube[j,:,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                    if dim_flag == 0: temp = cube[j,int(ilat[itime[i]] + yoffset),int(ilon[itime[i]] + xoffset)]
                if dim_flag == 1: dat[:,i] = temp.data
                if dim_flag == 0: dat[i] = temp.data
                if np.size(itime) > 1:
                    if stash_flag == 1: dat[dat==0] = np.nan              # set zeros to nans
                    if dim_flag == 1: data[:,j] = np.nanmean(dat,1)     # mean over time indices
                    if dim_flag == 0: data[j] = np.nanmean(dat)     # mean over time indices
                    # print 'averaging over itime...'
                else:
                    if dim_flag == 1: data[:,j] = np.squeeze(dat)                   # if only one index per hour
                    if dim_flag == 0: data[j] = np.squeeze(dat)                   # if only one index per hour
                    # print 'no averaging, itime = 1...'
        # print data
        # print 'data.shape = ', data.shape

        #################################################################
        ## FIGURES TO TEST OUTPUT
        #################################################################
        ### timeseries of lowest model level
        # plt.figure(figsize=(7,5))
        # plt.plot(cubetime[:-1],data[0:10,:])
        # plt.show()

        ### vertical profile of 1st timestep
        # plt.figure(figsize=(7,5))
        # plt.plot(data[:,0],cube.coord('model_level_number').points)
        # plt.show()

        ### pcolormesh of timeseries
        # plt.figure(figsize=(7,5))
        # plt.pcolormesh(cubetime[:-1], cube.coord('model_level_number').points, data)
        # plt.colorbar()
        # plt.show()

        #################################################################
        ## CREATE CUBE
        #################################################################
        ### ECMWF FIELD NAMES
        # field_names = {'forecast_time','pressure','height','temperature','q','rh','ql','qi','uwind','vwind','cloud_fraction',
        #             'wwind','gas_atten','specific_gas_atten','specific_dry_gas_atten','specific_saturated_gas_atten','K2',
        #             'specific_liquid_atten','sfc_pressure','sfc_height_amsl'};

        varname = varnames.findfieldName(stash)
        print ('standard_name = ', cube.standard_name)
        print ('long name = ', cube.long_name)
        print ('varname = ', varname)
        print ('')

        ntime = DimCoord(cubetime[:-1], var_name = 'forecast_time', standard_name = 'time', units = 'h')
        if dim_flag == 1:             ### 4D VARIABLE
            model_height = DimCoord(cube.aux_coords[2].points, var_name = 'height', standard_name = 'height', units='m')
            comdata = fixHeight(data, cube)
            ncube = Cube(np.transpose(data),
                    dim_coords_and_dims=[(ntime, 0),(model_height, 1)],
                    standard_name = cube.standard_name,
                    long_name = cube.long_name,
                    units = cube.units,
                    var_name = varname,
                    )
        elif dim_flag == 0:             ### 3D VARIABLE
            ncube = Cube(np.transpose(data),
                    dim_coords_and_dims=[(ntime, 0)],
                    standard_name = cube.standard_name,
                    long_name = cube.long_name,
                    units = cube.units,
                    var_name = varname,
                    )
        ncube.attributes = cube.attributes
        ### for consistency with multi-diag option
        fcube = ncube

    #################################################################
    ## define output filename
    #################################################################
    # print 'Define pp stream outfile:'
    # pp_outfile = date[:6] + str(int(date[6:8])+1) + '_oden_metum_' + str(stream[2:3]) + '.pp'
    # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
    ### bespoke setup if dir is 20180831T1200Z (for 20180901 data)
    # if date == '20180831T1200Z': nc_outfile = '20180901_oden_metum.nc'
    # print 'Outfile = ', pp_outfile

    ### save cube to netcdf file
    print ('')
    print ('Writing fcube to file:')
    print ('')
    if stream[1:3] == 'pc':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pc stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('')
        if not os.path.exists(nc_outfile):
            if 'fcube' in locals():
                out = writeNetCDF(date, fcube, nc_outfile)
            # if PC outfile already exists, combine other stream data
            # if PC outfile doesn't exist, write new

    if stream[1:3] == 'pd':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pd stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('***file is merged to outfile later***')
        print ('')
        doutfile = nc_outfile[:-3] + '_d.nc'
        if not os.path.exists(doutfile):
            if 'fcube' in locals():
                out = writePD_BL(fcube, doutfile)
            # if PD outfile already exists, combine other stream data
            # if PD outfile doesn't exist, write new

    if stream[1:3] == 'pe':
        ## Combine track-pulled pp output files to one netCDF
        ## First, make netCDF with pd stream (using Iris cubes)
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so making netCDF file with iris')
        print ('***file is merged to outfile later***')
        print ('')
        eoutfile = nc_outfile[:-3] + '_e.nc'
        if not os.path.exists(eoutfile):
            if 'fcube' in locals():
                out = writeFile_netCDF4(fcube, eoutfile)
            # if PC outfile already exists, combine other stream data
            # if PC outfile doesn't exist, write new

    elif stream[1:3] == 'pb':
        print ('fcube = ')
        print (fcube)
        print ('')
        print ('******')
        print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
        print ('***file is merged to outfile later***')
        print ('')
        ## Next, append 1D timeseries (surface) data (pb stream)
        ## Can't use Iris for this as cubes can't be 1D
        ##              -> uses standard netCDF appending function
        boutfile = nc_outfile[:-3] + '_b.nc'
        if not os.path.exists(boutfile):
            if 'fcube' in locals():
                out = writePB_Cloudnet(fcube, boutfile)     ##!!!! NEEDS UPDATING TO ONLY WRITE VARIABLES IN FILE, NOT HARD CODED

    elif stream[1:3] == 'pa':
        print ('Stream = ' + stream[1:] + ', so writing to new netCDF file with netCDF4.Dataset')
        print ('***file is merged to outfile later***')
        print ('')
        ## Next, append 1D timeseries (surface) data (pb stream)
        ## Can't use Iris for this as cubes can't be 1D
        ##              -> uses standard netCDF appending function
        aoutfile = nc_outfile[:-3] + '_a.nc'
        if not os.path.exists(aoutfile):
            if 'fcube' in locals():
                out = writePA_Analysis(fcube, aoutfile)

    return nc_outfile

def writeNetCDF(date, cube, nc_outfile):

    #################################################################
    ## CREATE NETCDF
    #################################################################
    #################################################################
    ## define output filename
    #################################################################
    print ('******')
    print ('Define .nc stream outfile:')
    # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
    print ('Outfile will be = ', nc_outfile)

    #################################################################
    ## load in each stream
    #################################################################
    ### USE IRIS TO SAVE OUT PC CUBE TO NETCDF (CREATING NEW FILE):
    # -------------------------------------------------------------
    # Convert .pp to .nc
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Converting to netCDF:')
    print ('')
    # cube = iris.load(outfile[0], global_con, callback)
    iris.save(cube, nc_outfile)

    return nc_outfile

def writePB_Cloudnet(cube, boutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 1D data to ' + boutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(boutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    print (cube)

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 0000 UTC.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write pbXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
        dat.scale_factor = float(1)
        dat.add_offset = float(0)
        dat.units = str(cube[d].units)
        dat.STASH = str(cube[d].attributes['STASH'])
        if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
        if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
        dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writePA_Analysis(cube, aoutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 1D data to ' + aoutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(aoutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    print (cube)
    # print cube[0].dim_coords

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))     ## use net sw to define, 1st in fcube

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 0000 UTC.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write paXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        if np.size(cube[d].dim_coords[0],0) == 24:      ### ignore 3-hourly data for now
            print ('Writing ' + cube[d].var_name)
            print ('')
            if not cube[d].var_name in dataset.variables:
                dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                dat.units = str(cube[d].units)
                dat.STASH = str(cube[d].attributes['STASH'])
                if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
                if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
                dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writePD_BL(cube, doutfile):
    #################################################################
    ## Write boundary layer diagnosticsa (PD) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing 3D data to ' + doutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(doutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')
    print ('Cube is: ')
    print (cube)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    #### find first occurrence of 2D variable, then break
    for l in range(0,len(cube)):
        if np.ndim(cube[l]) == 2:
            lind = l
            print ('height dim based on ' )
            print (cube[l])
            break

    forecast_time = dataset.createDimension('forecast_time', np.size(cube[lind].dim_coords[0].points))
    height = dataset.createDimension('height', np.size(cube[lind].dim_coords[1].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 0000 UTC.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[lind].dim_coords[0].points      ### forecast time (ignore first 12h)

    #### height
    height = dataset.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = ''
    height.units = 'm'
    height.long_name = 'height'
    height[:] = cube[lind].dim_coords[1].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        print
        if np.ndim(cube[d]) == 2:
            if cube[d].var_name == 'air_pressure': continue
            dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time','height',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            dat.units = str(cube[d].units)
            dat.STASH = str(cube[d].attributes['STASH'])
            if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
            if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
            # if np.size(cube[d].data,1) == 70:         ### need this for UM_RA2T_TKEdissrate
            #     dat[:,:] = cube[d].data
            # elif np.size(cube[d].data,1) == 69:
            #     dat[:,:-1] = cube[d].data
            #     dat[:,-1] = np.nan
            dat[:,:] = cube[d].data
        elif np.ndim(cube[d]) == 1:
            dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time',), fill_value='-9999')
            dat.scale_factor = float(1)
            dat.add_offset = float(0)
            dat.units = str(cube[d].units)
            dat.STASH = str(cube[d].attributes['STASH'])
            if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
            if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
            dat[:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def writeFile_netCDF4(cube, eoutfile):
    #################################################################
    ## Write 1D timeseries Cloudnet data (PB) to newly created netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    # print 'Appending 1D data to ' + outfile
    print ('Writing 3D data to ' + eoutfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(eoutfile, 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (dataset.file_format)
    print ('')

    # print cube
    # print cube[0].dim_coords

    ###################################
    ## Switch off automatic filling
    ###################################
    dataset.set_fill_off()

    ###################################
    ## Data dimensions
    # ###################################
    # forecast_period = dataset.createDimension('forecast_period', 24)
    forecast_time = dataset.createDimension('forecast_time', np.size(cube[0].dim_coords[0].points))
    height = dataset.createDimension('height', np.size(cube[0].dim_coords[1].points))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    timem = dataset.createVariable('forecast_time', np.float64, ('forecast_time',), fill_value='-9999')
    timem.scale_factor = float(1)
    timem.add_offset = float(0)
    timem.comment = 'Hours since 0000 UTC.'
    timem.units = 'hours'
    timem.long_name = 'forecast_time'
    timem[:] = cube[0].dim_coords[0].points      ### forecast time (ignore first 12h)

    #### height
    height = dataset.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = ''
    height.units = 'm'
    height.long_name = 'height'
    height[:] = cube[0].dim_coords[1].points      ### forecast time (ignore first 12h)

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    ###################################
    ## Write paXXX stream diagnostics
    ###################################
    for d in range(0,len(cube)):
        print ('Writing ' + cube[d].var_name)
        print ('')
        dat = dataset.createVariable(cube[d].var_name, np.float64, ('forecast_time','height',), fill_value='-9999')
        dat.scale_factor = float(1)
        dat.add_offset = float(0)
        dat.units = str(cube[d].units)
        dat.STASH = str(cube[d].attributes['STASH'])
        if not cube[d].standard_name == None: dat.standard_name = str(cube[d].standard_name)
        if not cube[d].long_name == None: dat.long_name = str(cube[d].long_name)
        dat[:,:] = cube[d].data

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def appendMetaNetCDF(outfile, date, out_dir, model, swath):

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Appending metadata to ' + outfile)
    print ('')

    ###################################
    ## Open File
    ###################################
    dataset = Dataset(outfile, 'a', format ='NETCDF4_CLASSIC')
    # infile = net.Dataset("2015%s%s-160000_0.nc" % (month,day), "a")
    # print ''
    # print dataset.file_format
    # print ''

    ###################################
    ## Global Attributes
    ###################################
    dataset.title = 'Met Office Unified Model single-site (Oden) output during MOCCHA'
    revision = 'undefined'
    if out_dir[2:9] == 'u-bg610':
        micro = 'Cloud microphysics: Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983). '
        revision = 'Revision no. 1. '
    elif out_dir[2:9] == 'u-bl661':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. No aerosol processing. '
        revision = 'Revision no. 0. '
    elif out_dir[2:9] == 'u-bm410':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 2.00e8 /m3, mass = 3.00e-9 kg/kg. No aerosol processing. '
        revision = 'Revision no. 0. '
    elif out_dir[2:9] == 'u-bn068':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. Also set l_subgrid_qcl_mp to .true. to allow for turbulent production of mixed-phase cloud. '
        revision = 'Revision no. 3. '
    elif out_dir[2:9] == 'u-bp738':
        micro = 'Global model initialised with ERA-Interim reanalyses, LAM run with RA2M_CON configuration (as u-bg610, default run). Cloud microphysics: Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983). '
        revision = 'Revision no. 0. '
    elif out_dir[2:9] == 'u-bq791':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Fletcher (1962) for consistency with Wilson and Ballard (1999) microphysics]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. No aerosol processing. '
        revision = 'Revision no. 0. '
    elif out_dir[2:9] == 'u-bq798':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Meyers et al., (1992)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. No aerosol processing. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-br210':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation and coarse mode soluble aerosol number concentration input taken from UKCA daily average profiles provided by Ruth Price (University of Leeds) <eersp@leeds.ac.uk>. No aerosol processing. '
        revision = 'Revision no. 2. '
    elif out_dir[3:10] == 'u-br409':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. Passive aerosol processing. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-bu570':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. Sea ice albedo available between 28 Aug and 4 Sep. '
        revision = 'Revision no. 1. '
    elif out_dir[2:9] == 'u-bu687':
        micro = 'Cloud microphysics: Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983). Updated RHcrit profile for vn11.4. '
        revision = 'Revision no. 0. '
    elif out_dir[2:9] == 'u-bv926':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. l_subgrid_qcl_mp to .false. to not allow for turbulent production of mixed-phase cloud. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-bz429':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. 15 min temporal resolution of key 3D diagnostics. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-ca011':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. 15 min temporal resolution of key 3D diagnostics. Prognostic 1A BL scheme. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-ca012':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. Also set l_subgrid_qcl_mp to .true. to allow for turbulent production of mixed-phase cloud. Includes diagnosed turbulent dissipation rate diagnostic. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-ca362':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. CICE sea ice scheme. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-cc278':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-cc324':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. Also set l_subgrid_qcl_mp to .true. to allow for turbulent production of mixed-phase cloud. Extended BL diagnostic list. '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-cc568':
        micro = 'Cloud microphysics: Smith (1990) but includes a cloud/precipitation microphysical scheme with prognostic ice (Wilson and Ballard, 1999), based on Rutledge and Hobbs (1983). Extended BL diagnostic list. Updated revision of suite u-bg610. '
        revision = 'Revision no. 1. '
    elif out_dir[3:10] == 'u-cd847':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation and coarse mode soluble aerosol number concentration input taken from UKCA daily average profiles provided by Ruth Price (University of Leeds) <eersp@leeds.ac.uk>. No aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-ce112':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation and coarse mode soluble aerosol number concentration input taken from UKCA daily average profiles provided by Ruth Price (University of Leeds) <eersp@leeds.ac.uk>. Passive aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-ce627':
        micro = 'Cloud microphysics: Both the global model and LAM use the PC2 (Wilson et al., 2008) cloud scheme (i_cld_vn = 2); specifically, the LAM uses the RA2T_CON configuration. Also set l_subgrid_qcl_mp to .true. to allow for turbulent production of mixed-phase cloud. Extended BL diagnostic list. Mid-level convection switched off in GLM (LAM untouched, as u-cc324). '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-cg179':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. Passive aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). '
        revision = 'Revision no. 0. '
    elif out_dir[3:10] == 'u-cl349':
        micro = 'CASIM microphysics + cloud scheme (i_cld_vn = 1). Double-moment [droplet activation = Abdul-Razzak and Ghan (2000); ice nucleation = Cooper (1986)]. 3 modes of soluble aerosol, no insoluble aerosol. Accumulation mode soluble aerosol: num = 1.00e8 /m3, mass = 1.50e-9 kg/kg. Aitken and coarse modes = 0. No aerosol processing. Updated RHcrit profile for vn11.4. Uses sea ice options from the global model (alpham = 0.72 [from 0.5], dtice = 2.0 [from 5.0]). Surface fluxes from JULES. '
        revision = 'Revision no. 0. '
    else:
        micro = '<MICROPHYSICS UNDEFINED IN META>'
    wind = 'U and V wind components interpolated on to common vertical grid. '
    dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young McCusker <G.Y.McCusker@leeds.ac.uk> using Python (Iris/netCDF4).'
    # dataset.source = 'UK Met Office Unified Model, version 11.1. Microphysics = ' + micro
    dataset.references = 'Rose suite ID: ' + out_dir[2:10]
    dataset.project = 'MOCCHA: Microbiology-Ocean-Cloud Coupling in the High Arctic. '
    if model == 'lam':
        modelnote = 'UM limited area model. '
        if swath == False:
            dataset.description = 'Hourly data taken from grid box closest to ship location. Where the ship covers more than one grid box within an hour period, data are averaged from all grid boxes crossed. '
        if swath == True:
            dataset.description = 'Data are from swath area surrounding ship track. Swath setup changes at 2 Sep. '
    elif model == 'glm':
        modelnote = 'UM global model. '
        dataset.description = 'Hourly data taken from grid box closest to ship location. '
    dataset.comment = revision + modelnote + micro + wind
    dataset.institution = 'University of Leeds.'
    dataset.initialization_time = date[0:4] + '-' + date[4:6] + '-' + date[6:8] + ' ' + date[9:14] + '.'

    ###################################
    ## Additional variables
    ###################################
    #### Model resolution
    if not 'horizontal_resolution' in dataset.variables.keys():
        if model == 'lam':
            res = dataset.createVariable('horizontal_resolution', np.float32, fill_value='-9999')
            res.comment = 'Horizontal grid size of nested region.'
            res.units = 'km'
            res[:] = 1.5
        elif model == 'glm':
            res = dataset.createVariable('horizontal_resolution', np.float32, fill_value='-9999')
            res.comment = 'Horizontal grid size of global model.'
            res.units = 'km'
            res[:] = 17.0

    ###################################
    ## Open pbXXX netCDF file
    ###################################
    boutfile = outfile[:-3] + '_b.nc'

    if os.path.exists(boutfile):
        ncB = Dataset(boutfile, 'r')

        ###################################
        ## Append pbXXX stream diagnostics
        ###################################

        print ('Appending pbXXX diagnostics:')
        print ('---')
        for d in ncB.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                dat = dataset.createVariable(d, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                if getattr(ncB.variables[d],'units', None):
                    dat.units = str(ncB.variables[d].units)
                else:
                    dat.units = 'unknown'
                if getattr(ncB.variables[d],'STASH', None):
                    dat.STASH = str(ncB.variables[d].STASH)
                if getattr(ncB.variables[d],'standard_name', None):
                    dat.standard_name = str(ncB.variables[d].standard_name)
                if getattr(ncB.variables[d],'long_name', None):
                    dat.long_name = str(ncB.variables[d].long_name)
                dat[:] = ncB.variables[d][:]

        ###################################
        ## Close read-only pbXXX file
        ###################################
        ncB.close()

    ###################################
    ## Open paXXX netCDF file
    ###################################
    aoutfile = outfile[:-3] + '_a.nc'

    if os.path.exists(aoutfile):
        ncA = Dataset(aoutfile, 'r')

        ###################################
        ## Append paXXX stream diagnostics
        ###################################
        print ('Appending paXXX diagnostics:')
        print ('---')
        for d in ncA.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                dat = dataset.createVariable(d, np.float64, ('forecast_time',), fill_value='-9999')
                dat.scale_factor = float(1)
                dat.add_offset = float(0)
                if getattr(ncA.variables[d],'units', None):
                    dat.units = str(ncA.variables[d].units)
                else:
                    dat.units = 'unknown'
                if getattr(ncA.variables[d],'STASH', None):
                    dat.STASH = str(ncA.variables[d].STASH)
                if getattr(ncA.variables[d],'standard_name', None):
                    dat.standard_name = str(ncA.variables[d].standard_name)
                if getattr(ncA.variables[d],'long_name', None):
                    dat.long_name = str(ncA.variables[d].long_name)
                dat[:] = ncA.variables[d][:]

        ###################################
        ## Close read-only paXXX file
        ###################################
        ncA.close()

    ###################################
    ## Open pdXXX netCDF file
    ###################################
    doutfile = outfile[:-3] + '_d.nc'

    print ('What variables do we have before pdXXX read in?:')
    print (dataset)

    if os.path.exists(doutfile):
        ncD = Dataset(doutfile, 'r')

        print (ncD)

        #### height on rho levels
        height2 = dataset.createDimension('height2', np.size(ncD.variables['height'][:]))
        height2 = dataset.createVariable('height2', np.float64, ('height2',), fill_value='-9999')
        height2.scale_factor = float(1)
        height2.add_offset = float(0)
        height2.comment = 'height coordinate on rho levels'
        height2.units = 'm'
        height2.long_name = 'height'
        height2[:] = ncD.variables['height'][:]      ### forecast time (ignore first 12h)

        ###################################
        ## Append pdXXX stream diagnostics
        ###################################
        print ('Appending pdXXX diagnostics:')
        print ('---')
        for d in ncD.variables:
            print (d)
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                if np.ndim(ncD.variables[d]) == 2:
                    print ('Variable is 2D:')
                    daat = dataset.createVariable(d, np.float64, ('forecast_time', 'height2',), fill_value='-9999')
                    daat.scale_factor = float(1)
                    daat.add_offset = float(0)
                    if getattr(ncD.variables[d],'units', None):
                        daat.units = str(ncD.variables[d].units)
                    else:
                        daat.units = 'unknown'
                    if getattr(ncD.variables[d],'STASH', None):
                        daat.STASH = str(ncD.variables[d].STASH)
                    if getattr(ncD.variables[d],'standard_name', None):
                        daat.standard_name = str(ncD.variables[d].standard_name)
                    if getattr(ncD.variables[d],'long_name', None):
                        daat.long_name = str(ncD.variables[d].long_name)
                    daat[:,:] = ncD.variables[d][:,:]
                elif np.ndim(ncD.variables[d]) == 1:
                    print ('Variable is 1D:')
                    dat = dataset.createVariable(d, np.float64, ('forecast_time', ), fill_value='-9999')
                    dat.scale_factor = float(1)
                    dat.add_offset = float(0)
                    if getattr(ncD.variables[d],'units', None):
                        dat.units = str(ncD.variables[d].units)
                    else:
                        dat.units = 'unknown'
                    if getattr(ncD.variables[d],'STASH', None):
                        dat.STASH = str(ncD.variables[d].STASH)
                    if getattr(ncD.variables[d],'standard_name', None):
                        dat.standard_name = str(ncD.variables[d].standard_name)
                    if getattr(ncD.variables[d],'long_name', None):
                        dat.long_name = str(ncD.variables[d].long_name)
                    dat[:] = ncD.variables[d][:]

        ###################################
        ## Close read-only pdXXX file
        ###################################
        ncD.close()

    ###################################
    ## Open peXXX netCDF file
    ###################################
    eoutfile = outfile[:-3] + '_e.nc'

    if os.path.exists(eoutfile):
        ncE = Dataset(eoutfile, 'r')

        ###################################
        ## Append paXXX stream diagnostics
        ###################################
        print ('Appending peXXX diagnostics:')
        print ('---')
        for d in ncE.variables:
            if d == 'forecast_time': continue
            if not d in dataset.variables:
                print ('Writing ' + d)
                print ('')
                daat = dataset.createVariable(d, np.float64, ('forecast_time', 'height',), fill_value='-9999')
                daat.scale_factor = float(1)
                daat.add_offset = float(0)
                if getattr(ncE.variables[d],'units', None):
                    daat.units = str(ncE.variables[d].units)
                else:
                    daat.units = 'unknown'
                if getattr(ncE.variables[d],'STASH', None):
                    daat.STASH = str(ncE.variables[d].STASH)
                if getattr(ncE.variables[d],'standard_name', None):
                    daat.standard_name = str(ncE.variables[d].standard_name)
                if getattr(ncE.variables[d],'long_name', None):
                    daat.long_name = str(ncE.variables[d].long_name)
                daat[:,:] = ncE.variables[d][:,:]

        ###################################
        ## Close read-only peXXX file
        ###################################
        ncE.close()

    ###################################
    ## Write out file
    ###################################
    dataset.close()

    return dataset

def callback(cube, field, filename):
    '''
    rename cube diagnostics per list of wanted stash diags
    '''

    iStash = cube.attributes['STASH'].__str__()
    if diags.findfieldName(iStash):
        if cube.name() != diags.findfieldName(iStash):
            cube.rename(diags.findfieldName(iStash))

def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags.returnWantedStash()

    # print GlobalStashList
    # print GlobalStashList[0]

    return GlobalStashList

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'JASMIN'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'LAPTOP':
        root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        ship_filename = '/home/gillian/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir = '25_u-cc568_RA2M_CON/'
    date_dir = os.listdir(root_dir + out_dir)

    ## 4_u-bg610_RA2M_CON/              # Wilson and Ballard 1999 uphys
    ## 5_u-bl661_RA1M_CASIM/            # 100/cc accum mode aerosol; ARG + Cooper
    ## 6_u-bm410_RA1M_CASIM/            # 200/cc accum mode aerosol
    ## 7_u-bn068_RA2T_CON/              # RA2T_CON nest + global 4D stash
    ## 8_u-bp738_RA2M_CON/              # ERAI
    ## 10_u-bq791_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Fletcher
    ## 11_u-bq798_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Meyers
    ## 12_u-br210_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507
    ## 13_u-br409_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; passive aerosol processing
    ## 14_u-bu570_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit
    ## 15_u-bu687_RA2M_CON/           # Wilson and Ballard 1999 uphys; new RHcrit
    ## 16_u-bv926_RA2T_CON/              # RA2T_CON nest + global 4D stash + no subgrid mp production
    ## 17_u-bz429_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics
    ## 18_u-ca011_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics; 1A BL Scheme
    ## 19_u-ca012_RA2T_CON/              # RA2T_CON nest + global 4D stash; includes diagnosed turbulent dissipation rate
    ## 20_u-ca362_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; CICE sea ice scheme
    ## 23_u-cc278_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM
    ## 24_u-cc324_RA2T_CON/              # RA2T_CON nest + global 4D stash; all BL diagnostics + sea ice albedo
    ## 25_u-cc568_RA2M_CON/             # Wilson and Ballard 1999 uphys. sea ice albedo and extra BL diags
    ## 26_u-cd847_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. identical suite = u-cd852
    ## 27_u-ce112_RA1M_CASIM/           # UKCA daily averaged aerosol profiles, GA6 albedo options. passive aerosol processing.
    ## 28_u-ce627_RA2T_CON/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included. Mid-level convection switched off in GLM.
    ## 30_u-cg179_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; passive aerosol processing
    ## 31_u-cl349_RA2M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM; jules fluxes

    #### run with nohup:
    ####    nohup python pullTrack_CloudNet.py > nohup_u-bn068_pullTrack_CloudNet.out &

    #### run on lotus (with batch_pullTrack.bsub):
    ####    bsub < batch_pullTrack.bsub

    #### check lotus run status with:
    ####    bjobs

    #### batch_pullTrack:
    #######         #!/bin/bash
    #######         #BSUB -q short-serial
    #######         #BSUB -J u-br210_pullTrack
    #######         #BSUB -o %J.out
    #######         #BSUB -e %J.err
    #######         #BSUB -W 06:00
    #######
    #######         python2.7 pullTrack_CloudNet.py


    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    print ('******')
    print ('')
    print ('Make stash list for cube read in at ' + time.strftime("%c"))
    print (' ')
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)
            ### defines which stash variables to load - should be within a loop

    for date in date_dir:
        ### just do 2018 dates
        # if date[0:4] == '2018':
        ### just do first date:
        if date == date_dir[0]:
        ### just do specific date
        # if date[0:8] == '20180913':
            # # -------------------------------------------------------------
            # # Load cube
            # # -------------------------------------------------------------
            print ('******')
            print ('')
            print ('Begin cube read in at ' + time.strftime("%c"))
            print (' ')
            # var_con = 'specific_humidity'
            # cube = iris.load_cube(filename1, var_con)
            # global_con = ['atmosphere_downward_eastward_stress','atmosphere_downward_northward_stress']

            grid_dirname = 'AUX_DATA/'
            if int(date[6:8]) <= 8: grid_filename = grid_dirname + date[:6] + '0' + str(int(date[6:8])+1) + '_ShipTrack_GRIDDED.csv'
            if int(date[6:8]) >= 9: grid_filename = grid_dirname + date[:6] + str(int(date[6:8])+1) + '_ShipTrack_GRIDDED.csv'

            ### bespoke setup if dir is 20180831T1200Z (for 20180901 data)
            if date == '20180831T1200Z': grid_filename = grid_dirname + '20180901_ShipTrack_GRIDDED.csv'

            ### -------------------------------------------------------------------------
            ### define input filename
            ### -------------------------------------------------------------------------
            # -------------------------------------------------------------
            # Define output stream filenames to look at:
            #           start at 012 if 3h dumps (a)
            #           start at 009 if 1h dumps in pb
            #           start at 011 if 1h dumps (c--e)
            # -------------------------------------------------------------
            # names = ['_pa009','_pb009','_pd011','_pe011','_pc011']
            # names = ['_pa012','_pb012','_pd012','_pe012','_pc012']
            names = ['_pc011']         ### only do specific files as a test
            if out_dir[-6:-1] == 'CASIM':
                expt = out_dir[-11:-1]
            elif out_dir[-4:-1] == 'CON':
                expt = out_dir[-9:-1]
            outfiles = [] ### define list to add processed filenames to

            for stream in names:
                ### -------------------------------------------------------------------------
                ### define output filename
                ### -------------------------------------------------------------------------
                if np.logical_or(np.logical_or(out_dir == '7_u-bn068_RA2T_CON/', out_dir == '24_u-cc324_RA2T_CON/'), out_dir == '28_u-ce627_RA2T_CON/'):    ## choose lam or global for 7_u-bn068/24_u-cc324/28_u-ce627
                    # #### LAM
                    filename = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream + '_r0.pp'
                    model = 'lam'
                    dirout = out_dir[-17:-10] + '_lam/'
                    ### GLM
                    if stream == '_pb009': stream = '_pb012'  ## hard fix for glm, pb stream starts at 012
                    if stream == '_pa009': stream = '_pa012'  ## hard fix for glm, pa stream *sometimes* starts at 012
                    # filename = root_dir + out_dir + date + '/' + date + '_glm' + stream + '_r0.pp'
                    # model = 'glm'
                    # dirout = out_dir[-17:-10] + '_glm/'
                else:
                    filename = root_dir + out_dir + date + '/' + date + '_HighArctic_1p5km_' + expt + stream + '_r0.pp'
                    model = 'lam'
                    if np.logical_or(np.logical_or(out_dir[0] == '1',out_dir[0] == '2'),out_dir[0] == '3'):
                        dirout = out_dir[3:10] + '/'
                    else:
                        dirout = out_dir[2:9] + '/'

                print ('Checking: ' + filename)
                exist_flag = 0 # initialise exist_flag
                if os.path.exists(filename):
                    exist_flag = 1
                    #### LOAD CUBE
                    if 'var_con' in locals():
                        print ('Loading single diagnostic:')
                        print (var_con)
                        cube1 = iris.load_cube(filename, var_con, callback)
                        con_flag = 0            # constraint flag
                    elif 'global_con' in locals():
                        print ('Loading multiple diagnostics:')
                        # cube = iris.load_cubes(filename1, global_con)
                        cube = iris.load(filename, global_con, callback)
                        con_flag = 1            # constraint flag
                        print (cube)

                    # ------------------------------------------------------------

                    ### -------------------------------------------------------------
                    ### Use the following to plot quick maps of loaded cubes
                    ### -------------------------------------------------------------

                    # hour = 0
                    # figure = plot_cartmap(ship_data, cube, hour, grid_filename)

                    ########################################################################
                    ### -------------------------------------------------------------
                    ### Pull gridded ship track from cube
                    ### -------------------------------------------------------------
                    ########################################################################
                    # -------------------------------------------------------------
                    ### 1. use the following if only want the exact ship position and no variability
                    # -------------------------------------------------------------
                    ### LOAD CUBE
                    nc_outfile = dirout + date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum.nc'
                    if date == '20180831T1200Z': nc_outfile = dirout + '20180901_oden_metum.nc'
                    aoutfile = nc_outfile[:-3] + '_a.nc'
                    boutfile = nc_outfile[:-3] + '_b.nc'
                    doutfile = nc_outfile[:-3] + '_d.nc'
                    eoutfile = nc_outfile[:-3] + '_e.nc'

                    ### 1-to-1 or swath perspective?
                    swath = True

                    if swath == False:
                        if stream[1:3] == 'pa':
                            if not os.path.exists(aoutfile):
                                print (aoutfile + ' does not exist, so pulling ship track...')
                                outfile = pullTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pb':
                            if not os.path.exists(boutfile):
                                print (boutfile + ' does not exist, so pulling ship track...')
                                outfile = pullTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pd':
                            if not os.path.exists(doutfile):
                                print (doutfile + ' does not exist, so pulling ship track...')
                                outfile = pullTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pe':
                            if not os.path.exists(eoutfile):
                                print (eoutfile + ' does not exist, so pulling ship track...')
                                outfile = pullTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pc':
                            if not os.path.exists(nc_outfile):
                                print (nc_outfile + ' does not exist, so pulling ship track...')
                                outfile = pullTrack_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        else:
                            print ('Valid stream not found.')
                    elif swath == True:
                        if stream[1:3] == 'pa':
                            if not os.path.exists(aoutfile):
                                print (aoutfile + ' does not exist, so pulling swath...')
                                outfile = pullSwath_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pb':
                            if not os.path.exists(boutfile):
                                print (boutfile + ' does not exist, so pulling swath...')
                                outfile = pullSwath_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pd':
                            if not os.path.exists(doutfile):
                                print (doutfile + ' does not exist, so pulling swath...')
                                outfile = pullSwath_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pe':
                            if not os.path.exists(eoutfile):
                                print (eoutfile + ' does not exist, so pulling swath...')
                                outfile = pullSwath_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        elif stream[1:3] == 'pc':
                            if not os.path.exists(nc_outfile):
                                print (nc_outfile + ' does not exist, so pulling swath...')
                                outfile = pullSwath_CloudNet(cube, grid_filename, global_con, stream, date, model, ship_data, nc_outfile)
                        else:
                            print ('Valid stream not found.')

                    # -------------------------------------------------------------
                    ### 2. use the following if only want the variability over a certain grid size
                    # -------------------------------------------------------------
                    ### LOAD CUBE
                    # nc_outfile = date[:6] + str(int(date[6:8])+1).zfill(2) + '_oden_metum_VAR.nc'
                    # if date == '20180831T1200Z': nc_outfile = '20180901_oden_metum_VAR.nc'
                    # aoutfile = nc_outfile[:-3] + '_a.nc'
                    # boutfile = nc_outfile[:-3] + '_b.nc'
                    # eoutfile = nc_outfile[:-3] + '_e.nc'
                    #
                    # if stream == '_pa012':
                    #     if not os.path.exists(aoutfile):
                    #         print aoutfile + ' does not exist, so pulling ship track...'
                    #         outfile = pullTrack_CloudNet_VAR(cube, grid_filename, global_con, stream, date)
                    # elif stream == '_pb009':
                    #     if not os.path.exists(boutfile):
                    #         print boutfile + ' does not exist, so pulling ship track...'
                    #         outfile = pullTrack_CloudNet_VAR(cube, grid_filename, global_con, stream, date)
                    # elif stream == '_pe011':
                    #     if not os.path.exists(eoutfile):
                    #         print eoutfile + ' does not exist, so pulling ship track...'
                    #         outfile = pullTrack_CloudNet_VAR(cube, grid_filename, global_con, stream, date)
                    # elif stream == '_pc011':
                    #     if not os.path.exists(nc_outfile):
                    #         print nc_outfile + ' does not exist, so pulling ship track...'
                    #         outfile = pullTrack_CloudNet_VAR(cube, grid_filename, global_con, stream, date)
                    # else:
                    #     print 'Valid stream not found.'

                    ########################################################################

                else:
                    print ('')
                    print ('****File does not exist****')
                    print ('')

                # if stream[1:3] == 'pc':
                #     if exist_flag == 1:
                #         ##-------------------------------------------------------------
                #         ## For each date, append metadata to netCDF
                #         ## -------------------------------------------------------------
                #         print ('******')
                #         print ('')
                #         print ('stream = ' + stream + ', so appending pa, pb, pd, pe (if present), and metadata')
                #         print ('')
                #         # outfile = '20180902_oden_metum.nc'
                #         out = appendMetaNetCDF(nc_outfile, date, out_dir, model, swath)
                #             ### final_outfile = root_dir + out_dir + 'OUT/' + nc_outfile
                #             ## os.rename(nc_outfile, final_outfile)

    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')

    #### DIAGNOSTICS TO CHOOSE FROM:

    ### paXXX
    # <iris 'Cube' of air_pressure_at_sea_level / (Pa) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of air_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of dew_point_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of relative_humidity / (%) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of specific_humidity / (1) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_air_pressure / (Pa) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_net_downward_longwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_net_downward_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of surface_temperature / (K) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of toa_incoming_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of toa_outgoing_shortwave_flux / (W m-2) (time: 8; grid_latitude: 500; grid_longitude: 500)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 8; grid_latitude: 501; grid_longitude: 500)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 8; grid_latitude: 501; grid_longitude: 500)>]

    #### 12 AUG ONLY - NO FULL NEST DIAGNOSTICS
    # <iris 'Cube' of surface_downwelling_longwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_downwelling_shortwave_flux_in_air / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_longwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_net_downward_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_incoming_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of toa_outgoing_shortwave_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>]


    ### pbXXX
    # 0: large_scale_ice_water_path / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 1: large_scale_liquid_water_path / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 2: eastward_wind_at_10m / (m s-1)      (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 3: northward_wind_at_10m / (m s-1)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 4: air_temperature_at_1.5m / (K)       (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 5: specific_humidity_at_1.5m / (1)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 6: relative_humidity_at_1.5m / (%)     (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 7: dew_point_temperature_at_1.5m / (K) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 8: turbulent mixing height after boundary layer / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 9: height_of_decoupled_layer_base / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 10: height_of_stratocumulus_cloud_base / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 11: combined_boundary_layer_type / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 12: cloud_area_fraction_assuming_random_overlap / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 13: cloud_area_fraction_assuming_maximum_random_overlap / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 14: wet_bulb_freezing_level_altitude / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 15: total_column_q / (unknown)          (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 16: air_pressure_at_sea_level / (Pa)    (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 17: atmosphere_boundary_layer_thickness / (m) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 18: high_type_cloud_area_fraction / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 19: low_type_cloud_area_fraction / (1)  (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 20: medium_type_cloud_area_fraction / (1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 21: stratiform_rainfall_flux / (kg m-2 s-1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 22: stratiform_snowfall_flux / (kg m-2 s-1) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 23: surface_air_pressure / (Pa)         (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 24: surface_temperature / (K)           (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 25: surface_upward_latent_heat_flux / (W m-2) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 26: surface_upward_sensible_heat_flux / (W m-2) (time: 3; grid_latitude: 25; grid_longitude: 25)
    # 27: water_evaporation_amount / (unknown) (time: 3; grid_latitude: 25; grid_longitude: 25)

    ### pcXXX
    # 0: total_radar_reflectivity / (unknown) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 1: air_pressure / (Pa)                 (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 2: air_temperature / (K)               (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 3: eastward_wind / (m s-1)             (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 4: large_scale_cloud_area_fraction / (1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 5: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 6: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 7: northward_wind / (m s-1)            (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 8: specific_humidity / (kg kg-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 9: upward_air_velocity / (m s-1)       (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)


    ### pdXXX
    # 0: entrainment_rate_for_surface_mixed_layer / (unknown) (grid_latitude: 25; grid_longitude: 25)
    # 1: entrainment_rate_for_boundary_layer / (unknown) (grid_latitude: 25; grid_longitude: 25)
    # 2: obukhov_length / (unknown)          (grid_latitude: 25; grid_longitude: 25)
    # 3: atmosphere_downward_eastward_stress / (Pa) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 4: atmosphere_downward_northward_stress / (Pa) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 5: turbulent_kinetic_energy / (unknown) (model_level_number: 69; grid_latitude: 25; grid_longitude: 25)
    # 6: air_pressure / (Pa)                 (model_level_number: 70; grid_latitude: 25; grid_longitude: 25)
    # 7: surface_downward_eastward_stress / (Pa) (grid_latitude: 25; grid_longitude: 25)
    # 8: surface_downward_northward_stress / (Pa) (grid_latitude: 25; grid_longitude: 25)
    # 9: surface_upward_water_flux / (kg m-2 s-1) (grid_latitude: 25; grid_longitude: 25)

if __name__ == '__main__':

    main()
