"""
Extracting single site data from the UKV
========================================

"""

import time
import datetime
import numpy as np
import iris
import cartopy.crs as ccrs
import os
import komodo.util.path as kmdP
import iris.util
import sys
sys.path.append('/home/h06/rdgman/ProfilesPy/diags')
import diags_Cardington as diags_ukv


ACC_STACH_CODE = ['m01s04i201', 'm01s04i202']


class columnClass():

    def __init__(self):
        self.site = None
        self.k = None
        self.rotated_lon = None
        self.rotated_lat = None
        self.col_rotated_lon = None
        self.col_rotated_lat = None
        self.col_index_lon = None
        self.col_index_lat = None
        self.lon = None
        self.lat = None
        #self.site_id = False


def findLatLonIndex(lat, lon, Lat_array, Long_array, k=None):
    '''
    Find the correct gate within a radar scan given a lat and long.
    This routine is using: scipy.spatial.cKDTree.

    Arg:

        + lat,lat:
            wanted latitude and longitude
        + cube:
            iris cube containing the radar scan
        + k: default: none
            number of neighbour
    return:
        + array of index from the cube.data array.
    '''

    import scipy

    # get lat/long from radar scan
    Rad_lat = Lat_array
    Rad_lon = Long_array

    # Shoe-horn data into format required by cKDTree
    source_grid_data = np.dstack([Rad_lat.ravel(), Rad_lon.ravel()])[0]
    destination_grid_data = list(np.dstack([lat, lon])[0])

    mytree = scipy.spatial.cKDTree(source_grid_data)

    dist, indices = mytree.query(destination_grid_data)

    idy = indices / int(Rad_lat.shape[1])
    idx = indices % Rad_lat.shape[1]

    windowLength = np.sqrt(k)
    a = windowLength - 2

    W = np.linspace(-a, a, num=windowLength)

    jdy = []
    jdx = []
    for i in W:
        for j in W:
            jdy.append(int(idy[0] + i))
            jdx.append(int(idx[0] + j))

    return (jdy[:], jdx[:])


def getRotatedUKSite(site_name, lon, lat, file_ext, rot_pole):
    '''
    '''

    # Find the coordinates of the station
    # Transform real lat, lon point into rotated coords
    #rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()
    ll = ccrs.Geodetic()
    target_xy = rot_pole.transform_point(lon, lat, ll)  # lower left corner
    rot_lon = target_xy[0] + 360.
    rot_lat = target_xy[1]

    # Print to check sites
    print '******'
    print site_name
    print 'Rotated lon  = ', lon
    print 'Rotated lat  = ', lat
    print 'Rotated lon coord = ', rot_lon
    print 'Rotated lat coord = ', rot_lat
    print ' '

    return (rot_lon, rot_lat, file_ext)


def GetColumnLocation(filename, site_name, lon, lat, file_ext, K=None):
    """
    the the x,y of all the column for each site location
    """

    Columns = {}

    var_con = iris.AttributeConstraint(STASH='m01s16i222')
    cube = iris.load_cube(filename, var_con)

    rot_pole = cube.coord('grid_latitude').coord_system.as_cartopy_crs()

    # get lat/long from radar scan
    Rad_lat = cube.coord('grid_latitude').points
    Rad_lon = cube.coord('grid_longitude').points

    m = len(Rad_lat)
    n = len(Rad_lon)

    Lat_array = np.ones((n, m)) * Rad_lat
    Long_array = np.ones((m, n)) * Rad_lon

    Lat_array = np.transpose(Lat_array)

    rot_lon, rot_lat, site = getRotatedUKSite(site_name, lon, lat, file_ext, rot_pole)

    k1 = findLatLonIndex(rot_lon, rot_lat, Long_array, Lat_array, k=K)

    Columns[site] = columnClass()
    Columns[site].site = site
    Columns[site].k = K
    Columns[site].rotated_lon = rot_lon
    Columns[site].rotated_lat = rot_lat
#    Columns[site].col_rotated_lon = Rad_lon[k1[1]][0]
#    Columns[site].col_rotated_lat = Rad_lat[k1[0]][0]
#    Columns[site].col_index_lon = k1[1][0]
#    Columns[site].col_index_lat = k1[0][0]
    Columns[site].col_rotated_lon = Rad_lon[k1[1]]
    Columns[site].col_rotated_lat = Rad_lat[k1[0]]
    Columns[site].col_index_lon = k1[1]
    Columns[site].col_index_lat = k1[0]
#    Columns[site].site_id = site_i

    return Columns


def makeGlobalStashList():
    '''
    make a list of all the stash code we want to load
    '''

    GlobalStashList = diags_ukv.returnWantedStash()

    return GlobalStashList


def fix_timeccord(cube):

    if not cube.coords('time', dim_coords=True):
        time = cube.coord('time')
        time.bounds = None
        iris.util.promote_aux_coord_to_dim_coord(cube, time)

    return cube


def set_attribute(cube, site_name, lon, lat):
    '''
    set attribute in the cube
    '''
    Site_full_name = site_name
    Site_longitude = lon
    Site_latitude = lat

    tmp_attribute = Site_full_name + '; Coordinates of RefSite: ' +\
        str(Site_latitude) + ' N, ' +\
        str(Site_longitude) + ' E'

    cube.attributes['Location'] = tmp_attribute

    return cube


def callback(cube, field, filename):
    iStash = cube.attributes['STASH'].__str__()
    if diags_ukv.findfieldName(iStash):
        if cube.name() != diags_ukv.findfieldName(iStash):
            cube.rename(diags_ukv.findfieldName(iStash))

    if iStash in ACC_STACH_CODE:
        cube.attributes['LBYR'] = field.lbyr
        cube.attributes['LBMON'] = field.lbmon
        cube.attributes['LBDAT'] = field.lbdat
        cube.attributes['LBHR'] = field.lbhr
        cube.attributes['LBMIN'] = field.lbmin

        cube.attributes['LBYRD'] = field.lbyrd
        cube.attributes['LBMOND'] = field.lbmond
        cube.attributes['LBDATD'] = field.lbdatd
        cube.attributes['LBHRD'] = field.lbhrd
        cube.attributes['LBMIND'] = field.lbmind

        cube.attributes['LBFT'] = field.lbft

    # to deal with accumulated field
    if isinstance(cube.coord('time').bounds, (np.ndarray, np.generic)):
        D = cube.coord('time').bounds
        if max(D.shape) == 2:
            cube.coords('time')[0].points = cube.coord('time').bounds[0, 1]
            cube.coords('forecast_period')[0].points = cube.coord(
                'forecast_period').bounds[0, 1]


def dealWithAcc(local_cube_list):

    import cf_units
    import iris.coords as icoords
    acc_CubeList = iris.cube.CubeList()
    Result = iris.cube.CubeList()
    for i in np.array(range(len(local_cube_list) - 1)):
        period_1 = i
        period_2 = i + 1
        Cvacc = iris.analysis.maths.subtract(local_cube_list[period_2],
                                             local_cube_list[period_1])
        time_unit = cf_units.Unit(
            'hours since 1970-01-01', calendar=cf_units.CALENDAR_GREGORIAN)

        # build ref Time coord
        LBYR = str(local_cube_list[period_1].attributes['LBYR'])
        if local_cube_list[period_1].attributes['LBMON'] >= 10:
            LBMON = str(local_cube_list[period_1].attributes['LBMON'])
        else:
            LBMON = '0' + \
                str(local_cube_list[period_1].attributes['LBMON'])
        if local_cube_list[period_1].attributes['LBDAT'] >= 10:
            LBDAT = str(local_cube_list[period_1].attributes['LBDAT'])
        else:
            LBDAT = '0' + \
                str(local_cube_list[period_1].attributes['LBDAT'])
        if local_cube_list[period_1].attributes['LBHR'] >= 10:
            LBHR = str(local_cube_list[period_1].attributes['LBHR'])
        else:
            LBHR = '0' + \
                str(local_cube_list[period_1].attributes['LBHR'])
        if local_cube_list[period_1].attributes['LBMIN'] > - 10:
            LBMIN = str(local_cube_list[period_1].attributes['LBMIN'])
        else:
            LBMIN = '0' + \
                str(local_cube_list[period_1].attributes['LBMIN'])

        scanTime = LBYR + '/' + LBMON + '/' + \
            LBDAT + ' ' + LBHR + ':' + LBMIN

        refTime = datetime.datetime.strptime(
            scanTime, "%Y/%m/%d %H:%M")

        refTimeCoord = icoords.AuxCoord(
            time_unit.date2num(refTime),
            standard_name='forecast_reference_time',
            units=time_unit)

        # Data ref Time coord
        LBYR = str(local_cube_list[period_2].attributes['LBYRD'])
        if local_cube_list[period_2].attributes['LBMOND'] >= 10:
            LBMON = str(local_cube_list[period_2].attributes['LBMOND'])
        else:
            LBMON = '0' + \
                str(local_cube_list[period_2].attributes['LBMOND'])
        if local_cube_list[period_2].attributes['LBDATD'] >= 10:
            LBDAT = str(local_cube_list[period_2].attributes['LBDATD'])
        else:
            LBDAT = '0' + \
                str(local_cube_list[period_2].attributes['LBDATD'])
        if local_cube_list[period_2].attributes['LBHRD'] >= 10:
            LBHR = str(local_cube_list[period_2].attributes['LBHRD'])
        else:
            LBHR = '0' + \
                str(local_cube_list[period_2].attributes['LBHRD'])
        if local_cube_list[period_2].attributes['LBMIND'] > - 10:
            LBMIN = str(local_cube_list[period_2].attributes['LBMIND'])
        else:
            LBMIN = '0' + \
                str(local_cube_list[period_2].attributes['LBMIND'])

        scanTime = LBYR + '/' + LBMON + '/' + \
            LBDAT + ' ' + LBHR + ':' + LBMIN

        dataTime = datetime.datetime.strptime(
            scanTime, "%Y/%m/%d %H:%M")

        dataTimeCoord = icoords.AuxCoord(time_unit.date2num(dataTime),
                                         standard_name='time',
                                         units=time_unit)

        # Create a set to contain the axis names for each data dimension.
        dim_names = []
        for dim in range(len(Cvacc.shape)):
            dim_coords = Cvacc.coords(contains_dimension=dim, dim_coords=True)
            dim_names.append(dim_coords[0].name())

        for iCoord in Cvacc.aux_coords:
            dim_names.append(iCoord.name())

        if 'time' not in dim_names:
            Cvacc.add_aux_coord(dataTimeCoord)
        if 'forecast_reference_time' not in dim_names:
            Cvacc.add_aux_coord(refTimeCoord)
        acc_CubeList.append(Cvacc)

    A = acc_CubeList.merge_cube()
    A.attributes['STASH'] = local_cube_list[0].attributes['STASH']
    A.attributes['um_version'] = local_cube_list[0].attributes['um_version']
    A.attributes['source'] = local_cube_list[0].attributes['source']
    STASH = A.attributes['STASH'].__str__()
    A.rename(diags_ukv.findfieldName(STASH))
    Result.append(A)

    return Result


def main():


    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    RUN = str(sys.argv[1])

    site_name = 'Cardington'
    lon = -0.41596544
    lat = 52.1015311
    file_ext = 'Cardington'

    # Date and Time of Run
    today = datetime.datetime.now()
    yesterday = datetime.datetime.now() - datetime.timedelta(days=1)
    if RUN == '21':
        YYYYMMDD = yesterday.strftime("%Y%m%d")
    else:
        YYYYMMDD = today.strftime("%Y%m%d")

    print 'Date and time of run: ' + YYYYMMDD + ' ' + RUN + 'Z'
    print ' '

    # File Name
    #pathname = '/hpc/scratch/d01/rdgman/tempqv/'
    pathname = '/project/rdgdata/rdgman/CurrentWeather/'
    filename = pathname + 'qwumqva_' + YYYYMMDD + '_qv' + RUN + '_pb000.pp'

    # Make directory for site
    makedir = 'mkdir /project/rdgdata/rdgman/PROFILE/netcdf/'+file_ext
    os.system(makedir)

    cube_dic = {}
    # -------------------------------------------------------------------------
    # If file exist
    # -------------------------------------------------------------------------
    if not kmdP.isFileExist(filename):
        tmpError = '\nThe File:\n       --> ' + filename + '\nDoes not exist!'
        raise NameError(tmpError)

    # -------------------------------------------------------------------------
    # Get location of a 3x3 box of model column
    # -------------------------------------------------------------------------
    Columns = GetColumnLocation(filename, site_name, lon, lat, file_ext, K=9)

    # -------------------------------------------------------------------------
    # make global stash list and constraint
    # -------------------------------------------------------------------------
    GlobalStashList = makeGlobalStashList()
    global_con = iris.AttributeConstraint(
        STASH=lambda stash: str(stash) in GlobalStashList)

    # -------------------------------------------------------------------------
    # Create time constrain to take data every hour
    # -------------------------------------------------------------------------
    #time_con = iris.Constraint(
    #    forecast_period=lambda xcell:  xcell in range(40))

    time_con = iris.Constraint(forecast_period=lambda xcell: any(np.isclose(xcell.point % 1, [0, 1./60.])))


    # -------------------------------------------------------------------------
    # Loop through forecast times
    # each time we need to load a new file as each file only contain 3 hours
    # -------------------------------------------------------------------------
    for tt in range(1, 19):
        CYCLE_TIME = time.time()
        # File Names
        TIME = str(tt * 2).zfill(3)
        filename_pb = pathname + 'qwumqva_' + \
            YYYYMMDD + '_qv' + RUN + '_pb' + TIME + '.pp'
        filename_pc = pathname + 'qwumqva_' + \
            YYYYMMDD + '_qv' + RUN + '_pc' + TIME + '.pp'
        filenames = [pathname + 'qwumqva_' + \
                     YYYYMMDD + '_qv' + RUN + '_p*' + TIME + '.pp']

        print 'Reading in files: '
        print filenames
        print ' '

        print 'Cubes at ' + TIME + ':'

        # ---------------------------------------------------------------------
        # Load all stash
        # ---------------------------------------------------------------------
        if kmdP.isFileExist(filename_pb):

            # -------------------------------------------------------------
            # FIX: 1
            # Some field are unknown, this will stop the concatenation
            # So fix name from stash list
            # Solve in the callback
            # -------------------------------------------------------------
            cube = iris.load(filenames, global_con, callback)

            # -------------------------------------------------------------
            # FIX: 2
            # fix time coordinate
            # it happen for all the snow field
            # -------------------------------------------------------------
            for icube in cube:
                STASH = icube.attributes['STASH'].__str__()
                if STASH not in ACC_STACH_CODE:
                    icube = fix_timeccord(icube)

        else:
            tmpError = '\nThe File:\n       --> ' + \
                filename_pb + '\nDoes not exist!'
            raise NameError(tmpError)

        # ---------------------------------------------------------------------
        # for each stash, extract a profile for each station at a time and save
        # ---------------------------------------------------------------------
        for iStash in GlobalStashList:
            STASH_TIME = time.time()
            print 'Stash: ', iStash, 'cycle: ', tt

            local_stash_constrain = iris.AttributeConstraint(STASH=iStash)
            local_cube_list = cube.extract(local_stash_constrain & time_con)
            if iStash in ACC_STACH_CODE:
                local_cube_list = dealWithAcc(local_cube_list)

            local_cube = None
            if len(local_cube_list) == 1:
                local_cube = local_cube_list[0]
            else:
                cube_list_final_tmp = iris.cube.CubeList()
                for icube in local_cube_list:
                    if len(icube.cell_methods) == 0:
                        cube_list_final_tmp.append(icube)
                if cube_list_final_tmp:
                    local_cube = cube_list_final_tmp.merge_cube()

            if local_cube:
                for key in Columns:
                    local_site = Columns[key]
                    cube_list_final = iris.cube.CubeList()

                    for a, b in zip(local_site.col_index_lat,
                                    local_site.col_index_lon):
                        if local_cube.ndim == 2:
                            cube_list_final.append(local_cube[a, b])
                        if local_cube.ndim == 3:
                            cube_list_final.append(local_cube[:, a, b])
                        if local_cube.ndim == 4:
                            cube_list_final.append(local_cube[:, :, a, b])

                    final_cube = cube_list_final.merge_cube()
                    final_cube = set_attribute(final_cube,
                                               site_name, lon, lat)
                    key = iStash + '_' + local_site.site
                    if key not in cube_dic.keys():
                        cube_dic[key] = iris.cube.CubeList()

                    cube_dic[key].append(final_cube)
            print '*** Stash time: ', time.time() - STASH_TIME #strftime("%d %b %Y %H:%M:%S")
        # end of one cycle
        print '*** Cycle time: ', time.time() - CYCLE_TIME

    print '*** FINAL time: ', time.time() - START_TIME

    print '**** output cubes... '
    for key in cube_dic:
        print key
        if len(cube_dic[key]) > 1:
            A = cube_dic[key].concatenate_cube()
        else:
            A = cube_dic[key][0]

        try:
            alt_coord = A.coord('altitude')
            alt_dims = A.coord_dims(alt_coord)
            alt_factory = A.aux_factory('altitude')
            A.remove_aux_factory(alt_factory)
            A.add_aux_coord(alt_coord, alt_dims)
            print 'Altitude added'
        except:
            print 'No altitude'

        ofile = '/project/rdgdata/rdgman/PROFILE/netcdf/'+file_ext+'/MOUKV_FC' + \
            YYYYMMDD + RUN + 'Z_' + key + '.nc'

        print 'Saving file: ' + ofile
        iris.save(A, ofile)
        print A

        print 'Copying to Reading'
        scp = 'scp '+ofile+' sws17m@oak.reading.ac.uk:/net/stor-nex-pool1.rdg.ac.uk/pool1/gold/met/radar/data/model3d/met-office-ukv'
        os.system(scp)

    print 'done'


if __name__ == '__main__':

    main()
