__author__ = 'lpeng'


from pylab import *
import pickle, itertools
import scipy.io
import scipy.stats, scipy.signal
import pandas as pd
import FFT, Plotting
from PET_Library import Data
import statsmodels.api as sm
import gc
from scipy.interpolate import interp1d
gc.collect()

from matplotlib import rc
rc('font', family='serif') #')'Times New Roman'

# #=============================path================================
datadir = '../Data'
workspace = '../workspace/201609'
# figdir = '/home/water5/lpeng/Figure/pan_spectral/201609'
figdir = '/home/water5/lpeng/Figure/pan_spectral/201705'
# figdir = '/Users/pitch/Google Drive/Figure/pan_spectral/201705'

##=============================variable============================
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']
variables = ['tavg', 'wind', 'ea', 'vpd']  # rh
vars_pan_var = ['pan', 'tavg', 'sun', 'wind', 'ea', 'vpd']
vars_penpan = ['tavg', 'tmax', 'tmin', 'p', 'ea', 'wind', 'sun', 'lat', 'elev']  # 'tc'
basinlongs=['Songhuajiang', 'Liaohe', 'Northwestern', 'Haihe', 'Yellow', 'Yangtze', 'Huaihe', 'Southeastern', 'Southwestern', 'Pearl']
geoinfo = load('%s/station_geoinfo' %workspace)
station_number = load('%s/basin_station_number' %workspace)

styr = 1961
edyr = 2001
stdy = datetime.datetime(styr, 12, 1)
eddy = datetime.datetime(edyr, 12, 31)
dates = pd.date_range(stdy, eddy, freq='D')
tstep = len(dates)
dyears = dates.year
dmonths = dates.month
doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])  # julian day

## quality check to get the station list, for Observed pan evaporation and the coherence needed variables, t, r, u, e
station_flag = pickle.load(open('station_sunhours_80_flag','rb'))
station_qc = [np.where(station_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
station_pan_flag = pickle.load(open('station_pan_80_flag','rb'))
station_pan_qc = [np.where(station_pan_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
good_stations = [intersect1d(station_qc[i], station_pan_qc[i]) for i in xrange(0, 10)]

def ann_missing_days():
	"seasonal cycle climatology"
	# vars_pan_var = ['pan', 'tavg', 'sun', 'wind', 'ea', 'vpd']
	keynames = ['Epan', 'Tair', 'sun hours',  'wind', '$e_a$', 'VPD']
	numbers = ['a', 'b', 'c', 'd', 'e', 'f']

	def annual_all_records(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.month)]
		return array(ann)

	def missing_data(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.isnan(g).sum() for _, g, in df.groupby(df.index.year)]
		return array(ann)

	ist = 0
	met_station = []
	arid = []

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))

		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			dataplot = [data[v][0, istation][0:tstep].flatten() for v in vars_pan_var]
			met_station.append(array([missing_data(pd.Series(dataplot[i], index=dates)) for i in xrange(6)]))
			ist = ist + 1

	met_station = array(met_station)

	fig, ax = plt.subplots(6, 1, figsize=(8, 10), sharex=True)
	Year_ticks = arange(0, 41, 10); Year_labels = arange(1961, 2002, 10)
	Stn_ticks = arange(0, 220, 50); Stn_labels = arange(0, 220, 50)
	vmaxs = (365, 12, 12, 12, 12, 12)
	pos = (0.08, 0.235, 0.39, 0.55, 0.705, 0.865)
	for i in xrange(6):
		im = ax[i].imshow(met_station[:, i, :].T, vmin=0, vmax=vmaxs[i], aspect='equal', cmap='YlOrRd', interpolation='none')
		ax[i].set_yticks(Year_ticks); ax[i].set_yticklabels(Year_labels)
		ax[i].set_xticks(Stn_ticks); ax[i].set_xticklabels(Stn_labels)
		ax[i].tick_params(axis='both', top='off', labelsize=11)
		ax[i].text(5, 10, '('+numbers[i]+') ' + keynames[i], ha='left', fontsize=12)
		cax = fig.add_axes([0.93, pos[5-i], 0.015, 0.094])
		cb = fig.colorbar(im, cax)  # adjust the size
		cb.ax.tick_params(labelsize=10)   # change the colorbar fontsize
		cb.ax.set_title('Days', fontsize=10)
		cax.yaxis.set_ticks_position('left')

	ax[5].set_xlabel('Station', fontsize=12)
	fig.tight_layout()

	plt.show();exit()
	# savefig('%s/figR1_ann_missing_day_228_station.png' %figdir)

	return

# ann_missing_days()
# exit()


# Import R and packages
# from rpy2.robjects import r
# import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate()
# r('source("/home/water5/lpeng/script/UTILS/pettit_test.r")')
# r('library("zyp")')

def Detect_Step_Changes(data,time_year,r):

	#Compute annual mean and variance
	mean = []
	std = []
	years = []
	for year in xrange(time_year[0],time_year[-1]+1):
		idx = np.where((time_year == year) & (np.isnan(data) == 0))
		if idx[0].size > 0:
			mean.append(np.mean(data[idx]))
			std.append(np.std(data[idx]))
		else:
			mean.append(float('NaN'))
			std.append(float('NaN'))
		years.append(year)
	mean = np.array(mean)
	std = np.array(std)
	years = np.array(years)

	#Detect step changes
	if np.where(np.isnan(data) == 0)[0].size > 1:
		p_mean = Pettitt_Test(years,mean,r)
		p_std = Pettitt_Test(years,std,r)

	#Note if station has a step change in mean or variance (p < 0.1)
	p_threshold = 0.01
	if (p_mean <= p_threshold) | (p_std <= p_threshold):
		flag_step_change = True
	else:
		flag_step_change = False

	return flag_step_change

def Pettitt_Test(time,x,r):

	idx = np.where(np.isnan(x) == 0)
	x = x[idx]
	time = time[idx]
	#Load data to R
	r.assign("time",time)
	r.assign("x",x)
	#Fit and remove kendall theil robust line
	r('slope <- zyp.sen(x~time)')
	b = r("slope$coefficients[1]")[0]
	m = r("slope$coefficients[2]")[0]
	x = x - (b + m*time)
	#Run pettitt test
	coeffs = r('pettitt.test(x,time)')
	p = coeffs[0][0]
	t = coeffs[1][0]

	return p


def ann_time_series():
	"seasonal cycle climatology"
	# vars_pan_var = ['pan', 'tavg', 'sun', 'wind', 'ea', 'vpd']
	keynames = ['Epan', 'Tair', 'sun hours',  'wind', '$e_a$', 'VPD']
	numbers = ['a', 'b', 'c', 'd', 'e', 'f']

	def annual_all_records(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.month)]
		return array(ann)

	def ann_mean(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.year)]
		return array(ann)

	ist = 0
	met_station = []
	arid = []

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))

		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			dataplot = [data[v][0, istation][0:tstep].flatten() for v in vars_pan_var]
			met_station.append(array([ann_mean(pd.Series(dataplot[i], index=dates)) for i in xrange(6)]))
			ist = ist + 1

	met_station = array(met_station)

	fig, ax = plt.subplots(6, 1, figsize=(8, 10), sharex=True)
	Year_ticks = arange(0, 41, 10); Year_labels = arange(1961, 2002, 10)
	Stn_ticks = arange(0, 220, 50); Stn_labels = arange(0, 220, 50)
	pos = (0.08, 0.235, 0.39, 0.55, 0.705, 0.865)
	for i in xrange(6):
		# plt.plot(met_station[10, i, :]);plt.show();exit()
		# print Detect_Step_Changes(met_station[20, i, :], arange(1961, 2002, 1) , r)

		im = ax[i].imshow(met_station[:, i, :].T, aspect='equal', cmap='YlGnBu', interpolation='none')
		ax[i].set_yticks(Year_ticks); ax[i].set_yticklabels(Year_labels)
		ax[i].set_xticks(Stn_ticks); ax[i].set_xticklabels(Stn_labels)
		ax[i].tick_params(axis='both', top='off', labelsize=11)
		ax[i].text(5, 10, '('+numbers[i]+') ' + keynames[i], ha='left', fontsize=16)
		cax = fig.add_axes([0.93, pos[5-i], 0.015, 0.1])
		cb = fig.colorbar(im, cax)  # adjust the size
		cb.ax.tick_params(labelsize=10)   # change the colorbar fontsize
		cax.yaxis.set_ticks_position('left')

	ax[5].set_xlabel('Station', fontsize=12)
	fig.tight_layout()

	# plt.show();exit()
	savefig('%s/figR3_ann_mean_var_228_station.png' %figdir)

	return

ann_time_series()
exit()




def ffill_plot():
	"seasonal cycle climatology"
	# vars_pan_var = ['pan', 'tavg', 'sun', 'wind', 'ea', 'vpd']

	def annual_all_records(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.month)]
		return array(ann)

	def missing_data(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.isnan(g).sum() for _, g, in df.groupby(df.index.year)]
		return array(ann)

	def Gapfill(daily):
		"return the array, not pandas object"
		ts = pd.Series(daily, index=dates).fillna(method='pad').values
		return ts

	fig, ax = plt.subplots(1, 1, figsize=(16, 3))

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		istation = 40
		st_idx = (stdy - datetime.datetime(1961, 1, 1)).days  # find out the start date
		ed_idx = (eddy - datetime.datetime(1961, 1, 1)).days + 1
			# [ 1  4  8 11 12 13 16 17 18 19 20 22 23 24 25 26 27 28 29 30 32 34 35 36 37 38 40 41 42 44 46]
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
		dataplot = data[vars_pan_var[2]][0, istation][st_idx:ed_idx].flatten()
		gapfill = Gapfill(dataplot)
		plt.title('Daily sun hours')
		plt.ylabel('hrs', fontsize=16)
		plt.plot_date(x=dates, y=gapfill, fmt='-o', c='red', label='Gap-filled')
		plt.plot_date(x=dates, y=dataplot, fmt='-o', c='blue', label='Original')
		plt.legend(loc=3, frameon=False)
		plt.show()
		# savefig('%s/figR2_ffill_example.png' %figdir)
		exit()
	return

# ffill_plot()

# Get the example lat, lon
# print np.where(good_stations[0]==40)[0]
# exit()
# lons = load('%s/lons_good_stations' %workspace)
# lats = load('%s/lats_good_stations' %workspace)
# print lons[26], lats[26]
# [ 125.13] [ 43.54]
# exit()
