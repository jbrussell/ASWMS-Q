# %% codecell
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# workingdir = './ORCA_M5.5_Zcorr_100km_snr3/'
workingdir = '../ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/'
azi_evlist = 'azi_evlist.txt'

# %% codecell
evs = pd.read_csv(workingdir+azi_evlist, delimiter= '\s+', index_col=False, header=None)
evs.columns = ["id", "evla", "evlo", "evdp", "mag"]
evs

m = Basemap(projection='aeqd',lat_0=lat_0,lon_0=lon_0,resolution='i')
# fill background.
m.drawmapboundary(fill_color='w')
# draw coasts and fill continents.
m.drawcoastlines(linewidth=0.5)

xev,yev=m(evs.evlo.values,evs.evla.values)
m.plot(xev, yev, 'o', color="r",markeredgecolor="black")

# r = 90
# x,y=m(lat_0,lon_0)
# x2,y2 = m(lat_0,lon_0+r) 
# circle1 = plt.Circle((x, y), y2-y, color='red',fill=False)
# ax.add_patch(circle1)

# m.fillcontinents(color='gray',lake_color='gray')
# 20 degree graticule.
# m.drawparallels(np.arange(-80,81,20))
# m.drawmeridians(np.arange(-180,180,20))
# draw a black dot at the center.
xpt, ypt = m(lon_0, lat_0)
m.plot([xpt],[ypt],'*',markersize=20,color='blue',markeredgecolor="black")
# draw the title.
plt.savefig(workingdir+'events_azi.pdf')
plt.show()
# %% codecell

# %% codecell
