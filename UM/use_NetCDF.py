###
###
### SCRIPT TO READ IN UM MODEL DATA IN NETCDF FORMAT AS IRIS CUBE
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
import os

def readfile(filename):

    import pandas as pd

    # print '******'
    print ''
    print 'Reading .txt file with pandas'
    print ''

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print 'Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4])
    print ''

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
    # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print '******'
    print ''
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print ''
    print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')'
    print 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6]))
    print 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7]))
    print ''

    return inIce_index

def writeoutGrid(tim, lat, lon, date):

    import pandas as pd

    # ******
    # write to csv file
    # ******

    print '******'
    print 'Writing ' + date + ' grid to file:'
    print ''
    dat = np.zeros([len(tim), 3])
    dat[:,0] = tim
    dat[:,1] = lon
    dat[:,2] = lat
    df = pd.DataFrame(dat)
    filename = 'AUX_DATA/' + date + '_ShipTrack_GRIDDED.csv'
    df.to_csv(filename,  sep = " ")
    print '... finished!'
    print ''
    print '******'

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
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='yellow')

    ### box pick 4-5h
    lon = np.array([277,277])
    lat = np.array([471,472])
    tim = np.zeros([np.size(lon)])
    tim[:] = 4.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='red')

    ### box pick 5-11h
    tim = np.arange(5,11)
    lon = np.zeros([np.size(tim)])
    lon[:] = 277
    lat = np.zeros([np.size(tim)])
    lat[:] = 472
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='blue')

    ### box pick 11-12h
    lon = np.array([277,277])
    lat = np.array([472,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 11.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='green')

    ### box pick 12-13h
    lon = np.array([277,276])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 12.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 13-15h
    tim = np.arange(13,15)
    lon = np.zeros([np.size(tim)])
    lon[:] = 276
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 15-16h
    lon = np.array([276,275])
    lat = np.array([471,471])
    tim = np.zeros([np.size(lon)])
    tim[:] = 15.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 16-19h
    tim = np.arange(16,19)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 471
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='green')

    ### box pick 19-20h
    lon = np.array([275,275])
    lat = np.array([471,470])
    tim = np.zeros([np.size(lon)])
    tim[:] = 19.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='black')

    ### box pick 20-22h
    tim = np.arange(20,22)
    lon = np.zeros([np.size(tim)])
    lon[:] = 275
    lat = np.zeros([np.size(tim)])
    lat[:] = 470
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][int(lon[i]) + xoffset], cube.dim_coords[1][int(lat[i]) + yoffset],color='red')

    ### box pick 22-23h
    lon = np.array([275,275,274])
    lat = np.array([470,469,469])
    tim = np.zeros([np.size(lon)])
    tim[:] = 22.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='blue')

    ### box pick 23-0h
    lon = np.array([274,274])
    lat = np.array([469,468])
    tim = np.zeros([np.size(lon)])
    tim[:] = 23.0
    tim_318 = np.append(tim_318, tim)
    lat_318 = np.append(lat_318, lat)
    lon_318 = np.append(lon_318, lon)
    for i in range(0,np.size(lon)):
        iplt.scatter(cube.dim_coords[2][lon[i] + xoffset], cube.dim_coords[1][lat[i] + yoffset],color='magenta')

    out = writeoutGrid(tim_318, lat_318, lon_318, date)

def trackShip(data):

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==12,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==13,data.values[:,1]==8),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print '******'
    print ''
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')'
    print 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')'
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4])
    print ''

    return trackShip_index

def findLatLon(ship_data, cube, hour):

    print ''
    print 'Finding lat/lon of ship track'
    print '...'

    #################################################################
    ## find ship track coordinates
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    drift_index = iceDrift(ship_data)
    inIce_index = inIce(ship_data)

    # -------------------------------------------------------------
    # Define unrotated coordinate model grid
    # -------------------------------------------------------------
    #### the following uses iris to unrotate the coordinate grid.
    ####    this only works with square domains (i.e. paXXX files)
    ####    only needs to be performed once -- saved grid as .csv file
    lon, lat = unrotateGrid(cube)

    print 'findLatLon testing:'
    print 'Ship (lon,lat): ' + str(ship_data.values[drift_index,7][0]) + ', ' + str(ship_data.values[drift_index,6][0])
    print 'Model unrotated [max, median], (lat,lon): ' + str(np.max(lat)) + ', ' + str(np.median(lon))


    ship_index = np.where(np.logical_and(np.greater_equal(lat[:],ship_data.values[drift_index,7][0]), np.less_equal(lat[:],ship_data.values[drift_index,7][1])))
    print 'Ship index test'
    print ship_index
    # print lat[ship_index[0]


    print 'test complete!'

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

    return lat, lon

def readGriddedTrack(grid_filename):

    import pandas as pd

    # print '******'
    print ''
    print 'Reading ' + grid_filename + ' file with pandas'
    print ''

    data = pd.read_csv(grid_filename, sep = " ")
    values = data.values

    tim = values[:,1]
    ilat = values[:,2]
    ilon = values[:,3]

    return tim, ilat, ilon

def plot_cartmap(ship_data, cube, hour, grid_filename): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
        # from matplotlib.patches import Polygon

    ###---------------------------------
    ### DEFINE OFFSETS DEPENDENT ON NEST ROI
    ###---------------------------------
    if cube[0,0].shape >= 25-1:    # ll = 240, 471
        xoffset = -239
        yoffset = -470
    elif cube[0,0].shape >= 93-1:    # ll = 211, 386
        xoffset = -210
        yoffset = -385
    elif cube[0,0].shape >= 500-1:
        xoffset = 0
        yoffset = 0

    print 'xoffset = ', xoffset
    print 'yoffset = ', yoffset

    ###################################
    ## PLOT MAP
    ###################################

    print '******'
    print ''
    print 'Plotting cartopy map:'
    print ''

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
    plt.figure(figsize=(10,9))
    # ax = plt.axes(projection=ccrs.Orthographic(0, 90))    # NP Stereo
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=30))

    ### set size
    ax.set_extent([20, 40, 89.6, 89.9], crs=ccrs.PlateCarree())       ### ZOOM
    # ax.set_extent([0, 60, 87.75, 90], crs=ccrs.PlateCarree())     ### SWATH
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
    if np.size(cube.shape) == 4:
        iplt.pcolormesh(cube[hour,0,:,:])
    elif np.size(cube.shape) == 3:
        iplt.pcolormesh(cube[hour,:,:])
        # iplt.pcolormesh(cube[hour,471:495,240:264])
    elif np.size(cube.shape) == 2:
        iplt.pcolormesh(cube[:,:])
    plt.title(cube.standard_name + ', ' + str(cube.units))
    plt.colorbar()

    #################################################################
    ## plot UM nest
    #################################################################
    ### draw outline of grid
    # qplt.outline(cube[hour,380:500,230:285])          ### original swath
    # qplt.outline(cube[hour,386:479,211:305])          ### redesigned swath (>13th)
    # qplt.outline(cube[hour,471:495,240:264])          ### 12-13th Aug swath
    # qplt.outline(cube[hour,450:495,220:305])          ### misc
    # qplt.outline(cube[hour,:,:])

    # gridship = gridShipTrack(cube, xoffset, yoffset)

            #### MID POINT: (433, 258)

    #################################################################
    ## plot ship track
    #################################################################
    ### DEFINE DRIFT + IN_ICE PERIODS
    # drift_index = iceDrift(ship_data)
    # inIce_index = inIce(ship_data)
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
    ## plot gridded ship track
    #################################################################
    ###
    tim, ilat, ilon = readGriddedTrack(grid_filename)

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

    print '******'
    print ''
    print 'Finished plotting cartopy map! :)'
    print ''

    # plt.savefig('FIGS/12-13Aug_Outline_wShipTrackMAPPED.svg')
    plt.show()

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
    print '******'
    print 'Test of unrotated coordinate grid: '
    print 'Rotated lon coord = ', rot_lon[250]
    print 'Rotated lat coord = ', rot_lat[250]
    print 'Lon = ', lon[250]
    print 'Lat = ', lat[250]
    print ' '

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

def main():

    START_TIME = time.time()
    print '******'
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'JASMIN'

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/UM/'
        ship_filename = '~/GWS/MOCCHA/ODEN/2018_shipposition_1hour.txt'
        position_filename = 'POSITION_UNROTATED.csv'
    if platform == 'LAPTOP':
        root_dir = '~/MOCCHA/UM/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    out_dir = '3_12AUG_SWATH_2FCSTS/'

    ## 1_20160401_61DIAG_TEST/
    ## 2_20180801_61DIAGS_TEST/2_30_86.625/
    ## 3_12AUG_SWATH_2FCSTS/

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print '******'
    print ''
    print 'Load in ship track file:'
    print ''
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    grid_dirname = 'AUX_DATA/'
    date = '20180812'
    grid_filename = grid_dirname + date + '_ShipTrack_GRIDDED.csv'

    print '******'
    print ''
    print 'Identifying .nc file: '
    print ''

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    filename1 = root_dir + out_dir + 'umnsaa_pa012_r0.nc'
    print filename1
    print ''

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print '******'
    print ''
    print 'Begin cube read in at ' + time.strftime("%c")
    print ' '
    var = 'surface_net_downward_shortwave_flux'
    cube = iris.load_cube(filename1, var)
    # data = Dataset(filename1,'r')

    print cube
    print ''

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
    # [<iris 'Cube' of cloud_area_fraction_assuming_maximum_random_overlap / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of cloud_area_fraction_assuming_random_overlap / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of m01s03i241 / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of Turbulent mixing height after boundary layer / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of wet_bulb_freezing_level_altitude / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of air_pressure_at_sea_level / (Pa) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of air_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of atmosphere_boundary_layer_thickness / (m) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of dew_point_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of high_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of low_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of medium_type_cloud_area_fraction / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of relative_humidity / (%) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of specific_humidity / (1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of stratiform_rainfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of stratiform_snowfall_flux / (kg m-2 s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_air_pressure / (Pa) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_temperature / (K) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_upward_latent_heat_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of surface_upward_sensible_heat_flux / (W m-2) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 24; grid_latitude: 25; grid_longitude: 25)>]


    ### pcXXX
    # <iris 'Cube' of cloud_volume_fraction_in_atmosphere_layer / (1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of m01s04i118 / (1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of air_pressure / (Pa) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of air_temperature / (K) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of specific_humidity / (kg kg-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of upward_air_velocity / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of x_wind / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>,
    # <iris 'Cube' of y_wind / (m s-1) (time: 25; model_level_number: 70; grid_latitude: 121; grid_longitude: 56)>]

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

    # FORECAST_PERIOD = cube1.aux_coords[1][:]

    # -------------------------------------------------------------
    # Define unrotated coordinate grid
    # -------------------------------------------------------------
    #### the following uses iris to unrotate the coordinate grid.
    ####    this only works with square domains (i.e. paXXX files)
    ####    only needs to be performed once -- saved grid as .csv file
    # lon, lat = unrotateGrid(cube)

    # hour = 0
    # test = findLatLon(ship_data, cube, hour)

    ############## DOESN'T WORK
    #### read in saved unrotated coordinate grid
    # position_data = readfile(position_filename)
    # lon = position_data.values[:,2]     ### unrotated longitude
    # lat = position_data.values[:,1]     ### unrotated latitude
    # lon, lat = np.meshgrid(lon,lat)     ### mesh for use with model diags

    # -------------------------------------------------------------
    # Plot data (map)
    # -------------------------------------------------------------
    ### select hour to plot
    hour = 0
    map = plot_cartmap(ship_data, cube, hour, grid_filename)#, lon, lat)



    END_TIME = time.time()
    print '******'
    print ''
    print 'End: ' + time.strftime("%c")
    print ''


if __name__ == '__main__':

    main()
