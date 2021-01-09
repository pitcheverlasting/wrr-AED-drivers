__author__ = 'lpeng'

"""
Construct the stations database with relatively complete records for spectral analysis
"""

from pylab import *
import pandas as pd
import pickle
import scipy.io
import scipy.stats
import Plotting

##=============================path================================
datadir = '../Data'
workspace = '../workspace/201609'
figdir = '/home/water5/lpeng/Figure/pan_spectral/201609'
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']  # time 1,2,3 year, month, day # Note that original data is multiplied by 10 from the normal unit. This process data devide 10 so that it becomes the original unit
varlongs = ['time', 'daily averaged pressure [hpa]', 'daily averaged air temperature [degree]', 'daily minimum air temperature [degree]', 'daily maximum air temperature [degree]', 'daily averaged vapor pressure [hpa]', 'daily averaged relative humidity [%]', 'daily averaged total cloud cover', 'daily averaged low cloud cover', 'daily averaged wind speed [m/s]', 'daily averaged surface temperature(at 0cm level) [degree]', 'daily averaged sunshine hour [hr]', 'daily averaged rainfall [mm]', 'daily averaged pan evaporation [mm]', 'daily averaged vapor pressure deficit [hpa]', 'daily saturated vapor pressure using daily averaged air temperature [hpa]', 'daily averaged relative humidity calculated [%]']

variables = ['tavg', 'wind', 'ea', 'vpd']
vars_pan_var = ['pan', 'tavg', 'sun', 'wind', 'ea', 'vpd']
vars_penpan = ['tavg', 'tmax', 'tmin', 'p', 'ea', 'wind', 'sun', 'lat', 'elev']  # 'tc'
basinlongs=['Songhuajiang', 'Liaohe', 'Northwestern', 'Haihe', 'Yellow', 'Yangtze', 'Huaihe', 'Southeastern', 'Southwestern', 'Pearl']
geoinfo = load('%s/station_geoinfo' %workspace)
station_number = load('%s/basin_station_number' %workspace)

## Time
styr = 1961
edyr = 2001
stdy = datetime.datetime(styr, 1, 1)
eddy = datetime.datetime(edyr, 12, 31)
dates = pd.date_range(stdy, eddy, freq='D')
tstep = len(dates)
dyears = dates.year
dmonths = dates.month
doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])  # julian day

# Step1. Select good station: no missing years and missing months 80% coverage
# Step1.1 Select good met station
def Quality_Check_station():
	"quality check stations that calculate Penpan models"
	flag_basins = []
	for ibasin in xrange(0, 10):
		flag = zeros((station_number[ibasin], 2))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			# for a station to count it must have more than 80% of coverage for each year and for all variables for Epan calculation
			for v in vars_penpan[:-2]:
				ts = data[v][0, istation][0:tstep]

				# Case1: If the whole series is 0, then throw away
				if np.nansum(ts) == 0.0:
					flag[istation, 0] = 1
					flag[istation, 1] = 1
					break
				# # Case2: For every month, if the missing value is greater than the threshold 20 % then mask down, stricter than year
				npoints_mon = 31.0*(1-80.0/100.0)
				for year in xrange(styr, edyr+1):
					for month in xrange(1, 13):
						count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()

						if count > int(npoints_mon):
							flag[istation, 0] = 1
							flag[istation, 1] = 2
							break

		flag_basins.append(flag)

	return flag_basins

# station_flag = Quality_Check_station()
# pickle.dump(station_flag, open('station_sunhours_80_flag', 'wb'))
# exit()

# Step1.2 Select good pan station
def Quality_Check_station_pan(percentage):
	"calculate the total number of missing records for each station"
	flag_basins = []
	for ibasin in xrange(0, 10):
		flag = zeros((station_number[ibasin], 2))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			# for a station to count it must have more than 80% of coverage for each year and for all variables for Epan calculation
			ts = data['pan'][0, istation][0:tstep]

			# Case1: If the whole series is 0, then throw away
			if np.nansum(ts) == 0.0:
				flag[istation, 0] = 1
				flag[istation, 1] = 1
				break
			# # Case2: For every month, if the missing value is greater than the threshold 20 % then mask down, stricter than year
			npoints_mon = 31.0*(1-percentage/100.0)
			for year in xrange(styr, edyr+1):
				for month in xrange(1, 13):
					count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()

					if count > int(npoints_mon):
						flag[istation, 0] = 1
						flag[istation, 1] = 2

		flag_basins.append(flag)

	pickle.dump(flag_basins, open('station_pan_%2d_flag','wb') % percentage)

	return

# Quality_Check_station_pan(60)
# exit()

# Step1.3 Select good pan station and good met station with t, sun, u, e for coherence spectral analysis
def Quality_Check_station_coherence(percentage):
	"calculate the total number of missing records for each station"
	flag_basins = []
	for ibasin in xrange(0, 10):
		flag = zeros((station_number[ibasin], 2))
		for istation in xrange(0, station_number[ibasin]):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			for v in vars_pan_var:
				ts = data[v][0, istation][0:tstep]

				# Case1: If the whole series is 0, then throw away
				if np.nansum(ts) == 0.0:
					flag[istation, 0] = 1
					flag[istation, 1] = 1
					break
				# # Case2: For every month, if the missing value is greater than the threshold 20 % then mask down, stricter than year
				npoints_mon = 31.0*(1-percentage/100.0)
				for year in xrange(styr, edyr+1):
					for month in xrange(1, 13):
						count = np.isnan(ts[(dyears==year)&(dmonths==month)]).sum()

						if count > int(npoints_mon):
							flag[istation, 0] = 1
							flag[istation, 1] = 2

		flag_basins.append(flag)

	pickle.dump(flag_basins, open('station_coh_%2d_flag'% percentage,'wb') )

	return

# Quality_Check_station_coherence(80)
# exit()

station_list_names = ['_all_stations', '_good_stations', '_for_stations', '_pan_stations'] #, '_coh_stations']
# station_lists = [station_number, good_stations, station_qc, station_pan_qc] #, station_coh]
fignames = ['all', 'good', 'forcing', 'pan']
cols = ['grey', 'springgreen', 'aquamarine', 'lightcoral']

# Step1.4 Plot good stations
def Get_Selected_Station_Locations(j):
	lons = []
	lats = []
	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in station_lists[j][ibasin]:
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			lons.append(geoinfo[index, 2])
			lats.append(geoinfo[index, 1])
		del data
	lons = vstack(lons)
	lats = vstack(lats)
	lons.dump('%s/lons%s' %(workspace, station_list_names[j]))
	lats.dump('%s/lats%s' %(workspace, station_list_names[j]))

	return

def Plot_Selected_Station_Locations(j):
	lons = load('%s/lons%s' %(workspace, station_list_names[j]))
	lats = load('%s/lats%s' %(workspace, station_list_names[j]))
	Plotting.Mapshow(cols[j], lons, lats, 45, 0.5, None, None, None, 'Qualified pan evaporation met stations (%s)' %len(lons), 'Met Stations', '', figdir, '%s_stations_scattermap.png' %fignames[j])

	return

# Get_Selected_Station_Locations(1)
# Plot_Selected_Station_Locations(1)
# exit()
