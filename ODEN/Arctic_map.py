from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm


# filename = '/data/scihub-users/giyoung/MAC/Seaice/ObsSeaice_NASATeam_NSIDC-0051_v1-1_JuntoAug2017_SeasonalMean.nc'
#
# nc = NetCDFFile(filename, 'r')
#
# lat = nc.variables['lat'][:]
# lon = nc.variables['lon'][:]
# data = np.squeeze(nc.variables['sic_mean'][:])

# Create a wider than normal figure to support subplots
fig = plt.figure(figsize=(7, 8.5), dpi=300)

# Also manually adjust the spacings which are used when creating subplots
plt.gcf().subplots_adjust(hspace=0.05, wspace=0.1, top=0.96, bottom=0.1,
                          left=0.05, right=0.95)


dim = 3500000

# for i in range(0,8):
	# plt.subplot(2, 4, i+1)
m = Basemap(width=0.75*dim,height=dim,
            resolution='l',projection='stere',\
            lat_ts=75,lat_0=75,lon_0=0)
m.drawcoastlines()
m.bluemarble()

# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,10.),color='k')
m.drawmeridians(np.arange(-180.,181.,20.),color='k')

# make lat/lon grid
# lons, lats = np.meshgrid(lon,lat)
# x, y = m(lons, lats)

# plot data
# csf = m.contourf(x,y,mslp[i,:,:]/float(100),contourf_levels,cmap=mpl_cm.viridis)
# cs = m.contour(x,y,mslp[i,:,:]/float(100),contour_levels,colors='k',linewidth=0.1)
# plt.clabel(cs,fontsize=8, inline=200,fmt='%1.f')


# m.fillcontinents(color='lightgrey')

x, y = m(lon, lat)
csf = m.pcolor(x,y,data,cmap=mpl_cm.Blues_r)

# add colorbar.
# cbaxes = fig.add_axes([0.2,0.07, 0.6, 0.03])  # This is the position for the colorbar
# cb = plt.colorbar(csf, cax = cbaxes, orientation='horizontal')
# cb.ax.xaxis.set_ticks_position('bottom')
# cb.ax.xaxis.set_label_position('bottom')
# cb.ax.axes.set_xlabel('MSLP [hPa]')

plt.savefig('FIGS/EuropeanArctic_vPOSTER.svg',dpi=100)
plt.show()
