__author__ = 'lpeng'

from pylab import *
import pickle, itertools
import scipy.io
import scipy.stats
from scipy.optimize import minimize
from scipy import signal
import statsmodels.tools.eval_measures as evaluate
import pandas as pd
import IO, FFT, Plotting
from PET_Library import Data
import gc
from scipy.interpolate import interp1d
gc.collect()

# #=============================path================================
datadir = '../Data'
workspace = '../workspace/201609'
# figdir = '/home/water5/lpeng/Figure/pan_spectral/201609'
figdir = '/home/water5/lpeng/Figure/pan_spectral/201705'

##=============================variable============================
vars = ['time', 'p', 'tavg', 'tmin', 'tmax', 'ea', 'rh', 'tc', 'lc', 'wind', 'ts', 'sun', 'rain', 'pan', 'vpd',  'estavg', 'rh_test']
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
doys = dates.dayofyear
# doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])  # julian day

## quality check using data availability as criteria
# the station list is for PET models
station_flag = pickle.load(open('station_sunhours_80_flag','rb'))
station_qc = [np.where(station_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
station_pan_flag = pickle.load(open('station_pan_80_flag','rb'))
station_pan_qc = [np.where(station_pan_flag[ibasin][:, 0]==0)[0] for ibasin in xrange(0, 10)]
good_stations = [intersect1d(station_qc[i], station_pan_qc[i]) for i in xrange(0, 10)]


###======================Datetime Toolkits=============================
def Gapfill(daily):
	"return the array, not pandas object"
	ts = pd.Series(daily, index=dates).fillna(method='pad').values
	return ts

def daily2annual(daily):
	ts = pd.Series(daily, index=dates).resample('A', how='mean').values
	return ts

def daily2monthly(daily):
	ts = pd.Series(daily, index=dates).resample('M', how='mean').values
	return ts

def daily2monthly_df(array):
	arr = pd.DataFrame(array.T, index=dates).resample('M', how='mean').values
	return arr.T

def daily2weekly_df(array):
	arr = pd.DataFrame(array.T, index=dates).resample('W', how='mean').values
	return arr.T

def daily2annual_df(array):
	arr = pd.DataFrame(array.T, index=dates).resample('A', how='mean').values
	return arr.T

###====================================================================

def msc_groupby_DI(arid):
	# aridrange = [0, 1, 1.5, 2.5, 5, 250]
	aridrange = [0, 2, 4, 8, 250]

	index_DI = []
	for igroup in xrange(0, len(aridrange)-1):
		low = aridrange[igroup]
		high = aridrange[igroup+1]
		index_DI.append(np.where((arid>=low) & (arid<high))[0])
	return index_DI

def msc_groupby_DI_more(arid):
	# aridrange = [0, 1, 1.5, 2.5, 5, 250]
	aridrange = [0, 1, 1.5, 2, 4, 8, 20, 40, 80, 160, 250]

	index_DI = []
	for igroup in xrange(0, len(aridrange)-1):
		low = aridrange[igroup]
		high = aridrange[igroup+1]
		index_DI.append(np.where((arid>=low) & (arid<high))[0])
	return index_DI

##########################################################################
# for spectral coherece analysis
##########################################################################
nf = 513
nbasin = 10
nvar = 8
sampling_frequency = 1/(24.0 * 3600.0) # unit: per day


def Coherence_Frequency():

	data = scipy.io.loadmat('%s/1_AP.mat' %(datadir))
	input = data[variables[0]][0, 0][0:tstep].flatten()
	pan = data['pan'][0, 0][0:tstep].flatten()
	freq = FFT.Coherence(input, pan, sampling_frequency, 'linear')[0]

	return freq

def cel2Kelvin(input):
	input['tmax'] = input['tmax']+273.16
	input['tmin'] = input['tmin']+273.16
	input['tavg'] = input['tavg']+273.16
	return input

def Kelvin2cel(input):
	input['tmax'] = input['tmax']-273.16
	input['tmin'] = input['tmin']-273.16
	input['tavg'] = input['tavg']-273.16
	return input

from matplotlib import rc
rc('font', family='serif') #')'Times New Roman'

##########################################################################
# Calculate penpan modelled PE
##########################################################################
def chunk_st_ed(npts):
	if npts%2 == 0: st = npts/2; ed= npts/2-1
	else: st = npts/2; ed = npts/2
	return st, ed

def Permute_each_month_station(input):

	input_rand = []
	# Set for a month
	for mon in xrange(1, 13):
		# pool out all the data within a month into a Dataframe
		idx_mon = np.where((dyears >= styr) & (dyears <= edyr) & (dmonths == mon))[0]  # input.loc[input.index.month==mon]
		df = input.iloc[idx_mon]
		# retrieve all the index in the original order
		days = df.index

		def shuffle(df):
			"each shuffle will produce a new dataframe with original order"
			# resample all the data by permutation
			df = df.sample(frac=1)   # other method: np.random.permutation()
			# store the original time step, it looks like not very useful
			# df['dates'] = df.index
			# reset the time index with original order
			df.index = days
			# df = df.reset_index(drop=True) this is to remove the original index
			return df

		input_rand.append(shuffle(df))
		# [shuffle(df) for i in xrange(10)]
		del df

	input_rand = pd.concat(input_rand).sort_index()

	# scale the mean to the original monthly mean
	# For temperature convert to kelvin
	input = cel2Kelvin(input)
	input_rand = cel2Kelvin(input_rand)

	month_mean_obs = input.resample('MS', how='mean')
	month_mean_exp = input_rand.resample('MS', how='mean')

	ratio = (month_mean_obs/month_mean_exp).reindex(dates, method='ffill')
	input_rand = input_rand * ratio

	# input = Kelvin2cel(input)
	input_rand = Kelvin2cel(input_rand)

	return input_rand.to_dict(orient='series')

def Sample_each_week_station(input):

	input_rand = []

	# select each day of year
	for d in xrange(1, 367):
		# For each single day, sample an ensemble with a 7-day window across multiple years
		wdays = np.arange(d-3, d+4)
		constant = ones((7)) * 366
		mask1 = constant * (wdays>366)
		mask2 = constant * (wdays<1)
		wdays_crop = wdays - mask1 + mask2

		idx_doy = np.where((doys == wdays[3]))[0]
		days = input.iloc[idx_doy].index  # retrieve all the index in the original order
		N = len(days)  # number of day of year across all years, especially for 366
		idx_doy_wind = [np.where((doys == wdays_crop[i]))[0] for i in xrange(0, 7)]
		idx_doy_wind = list(itertools.chain.from_iterable(idx_doy_wind))
		df = input.iloc[idx_doy_wind]

		def shuffle(df):
			"each shuffle will produce a new dataframe with original order"
			# resample all the data by permutation
			df = df.sample(n=N)
			# reset the time index with original order
			df.index = days
			return df
		input_rand.append(shuffle(df))
		del df

	input_rand = pd.concat(input_rand).sort_index()

	# scale the mean to the original monthly mean
	# For temperature convert to kelvin
	input = cel2Kelvin(input)
	input_rand = cel2Kelvin(input_rand)

	month_mean_obs = input.resample('MS', how='mean')
	month_mean_exp = input_rand.resample('MS', how='mean')

	ratio = (month_mean_obs/month_mean_exp).reindex(dates, method='ffill')
	input_rand = input_rand * ratio

	# input = Kelvin2cel(input)
	input_rand = Kelvin2cel(input_rand)

	# plt.plot(input_rand['wind'])
	# plt.plot(input['wind'])

	return input_rand.to_dict(orient='series')

def Sample_station_window(input, wind):

	input_rand = []

	for d in xrange(1, 367):
		# For each single day, sample an ensemble with a 7-day window across multiple years
		wdays = np.arange(d-wind/2, d+wind/2+1)
		constant = ones((wind)) * 366
		mask1 = constant * (wdays>366)
		mask2 = constant * (wdays<1)
		wdays_crop = wdays - mask1 + mask2

		idx_doy = np.where((doys == wdays[wind/2]))[0]
		days = input.iloc[idx_doy].index  # retrieve all the index in the original order
		N = len(days)  # number of day of year across all years, especially for 366
		idx_doy_wind = [np.where((doys == wdays_crop[i]))[0] for i in xrange(0, wind)]
		idx_doy_wind = list(itertools.chain.from_iterable(idx_doy_wind))
		df = input.iloc[idx_doy_wind]

		def shuffle(df):
			"each shuffle will produce a new dataframe with original order"
			# resample all the data by permutation
			df = df.sample(n=N)
			# reset the time index with original order
			df.index = days
			return df
		input_rand.append(shuffle(df))
		del df

	input_rand = pd.concat(input_rand).sort_index()

	# scale the mean to the original monthly mean
	# For temperature convert to kelvin
	input = cel2Kelvin(input)
	input_rand = cel2Kelvin(input_rand)

	# month_mean_obs = input.resample('MS', how='mean')
	# month_mean_exp = input_rand.resample('MS', how='mean')
	# ratio = (month_mean_obs/month_mean_exp).reindex(dates, method='ffill')

	mean_obs = input.resample('%sD' %wind, how='mean')
	mean_exp = input_rand.resample('%sD' %wind, how='mean')

	ratio = (mean_obs/mean_exp).reindex(dates, method='ffill')
	input_rand = input_rand * ratio

	# input = Kelvin2cel(input)
	input_rand = Kelvin2cel(input_rand)

	return input_rand.to_dict(orient='series')

def Remove_interannual_variability(input):

	# For temperature convert to kelvin
	input = cel2Kelvin(input)

	ann_mean_obs = input.resample('AS', how='mean')
	clim_obs = ann_mean_obs.mean()

	ratio = (clim_obs/ann_mean_obs).reindex(dates, method='ffill')
	input_niav = input * ratio

	# input = Kelvin2cel(input)
	input_niav = Kelvin2cel(input_niav)

	return input_niav.to_dict(orient='series')

def Remove_shortterm_variability(input, vars, npts):

	# For temperature convert to kelvin
	output = input.copy()
	output = cel2Kelvin(output)

	def moving_average(y, npts):
		return np.convolve(y, np.ones(npts)/npts, mode='same')

	def running_mean(y, npts):
		return pd.rolling_mean(y, npts, center=True) # [npts-1:]

	for var in vars:
		# output[var] = moving_average(output[var], npts)
		output[var] = running_mean(output[var], npts)

	input_smooth = Kelvin2cel(output)

	return input_smooth.to_dict(orient='series')

def Remove_variable_shortterm_variability(input, var, npts):

	# For temperature convert to kelvin
	output = input.copy()
	output = cel2Kelvin(output)

	def moving_average(y, npts):
		return pd.rolling_mean(y, npts)[npts-1:]

	output[var] = moving_average(output[var], npts)


	input_smooth = Kelvin2cel(output)

	return input_smooth.to_dict(orient='series')


def Calculate_Ep_daily():
	PENPAN = []
	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))

		for istation in good_stations[ibasin]:
			print ibasin, istation
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			# Read all the necessary input into a dataframe
			input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			input = pd.DataFrame.from_dict(input)
			input.index = dates

			# Remove inter-annual variability
			# INPUT = Remove_interannual_variability(input)

			npts = 15
			# INPUT = Remove_shortterm_variability(input, vars_penpan[:-2], npts)

			# Radiation
			# INPUT = Remove_variable_shortterm_variability(input, 'sun', npts)
			# ea
			# INPUT = Remove_variable_shortterm_variability(input, 'ea', npts)
			# wind
			# INPUT = Remove_variable_shortterm_variability(input, 'wind', npts)

			# tair
			INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg'], npts)

			INPUT['doy'] = doys
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]
			### Calculate Epan
			res = Data(INPUT, 'sunhours')
			PENPAN.append(res.penpan)

	pe_model = array(PENPAN)
	# pe_model.dump('%s/pe_mod_penpan_removeiav_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort15_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort7_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort15_sun_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort15_ea_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort15_wind_good_stations' %(workspace))
	# pe_model.dump('%s/pe_mod_penpan_removeshort15_tair_good_stations' %(workspace))

	return

# Calculate_Ep_daily()
# exit()


def Calculate_Ep_daily_smoothwindow():
	"only take tair as example"
	for npts in (7, 15, 30):
		PENPAN = []
		for ibasin in xrange(0, 10):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))

			for istation in good_stations[ibasin]:
				print ibasin, istation
				index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
				# Read all the necessary input into a dataframe
				input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
				input = pd.DataFrame.from_dict(input)
				input.index = dates

				# tair
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg'], npts)
				# wind
				# INPUT = Remove_shortterm_variability(input, ['wind'], npts)
				# humidity
				INPUT = Remove_shortterm_variability(input, ['ea'], npts)
				# sun
				# INPUT = Remove_shortterm_variability(input, [], npts)

				INPUT['doy'] = doys
				INPUT['lat'] = geoinfo[index, 1]
				INPUT['elev'] = geoinfo[index, 3]
				### Calculate Epan
				res = Data(INPUT, 'sunhours') #, npts)

				PENPAN.append(res.penpan)

		pe_model = array(PENPAN)
		# pe_model.dump('%s/pe_mod_penpan_removeshort%s_tair_good_stations' %(workspace, npts))
		# pe_model.dump('%s/pe_mod_penpan_removeshort%s_rnet_good_stations' %(workspace, npts))
		# pe_model.dump('%s/pe_mod_penpan_removeshort%s_wind_good_stations' %(workspace, npts))
		# pe_model.dump('%s/pe_mod_penpan_removeshort%s_ea_good_stations' %(workspace, npts))

	return

# Calculate_Ep_daily_smoothwindow()
# exit()

def Calculate_Ep_daily_smoothvariable():
	"for one period, four variables"
	npts = 7
	PENPAN = []
	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))

		for istation in good_stations[ibasin]:
			print ibasin, istation
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			# Read all the necessary input into a dataframe
			input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			input = pd.DataFrame.from_dict(input)
			input.index = dates

			# tair
			# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg'], npts)
			# wind
			# INPUT = Remove_shortterm_variability(input, ['wind'], npts)
			# humidity
			# INPUT = Remove_shortterm_variability(input, ['ea'], npts)
			# sun
			INPUT = Remove_shortterm_variability(input, [], npts)

			INPUT['doy'] = doys
			INPUT['lat'] = geoinfo[index, 1]
			INPUT['elev'] = geoinfo[index, 3]
			### Calculate Epan
			# res = Data(INPUT, 'sunhours')

			res = Data(INPUT, 'sunhours', npts)  # for rnet

			PENPAN.append(res.penpan)

	pe_model = array(PENPAN)
	pe_model.dump('%s/pe_mod_penpan_removeshort%s_rnet_good_stations' %(workspace, npts))

	return

# Calculate_Ep_daily_smoothvariable()
# exit()


def Calculate_Ep_daily_smooth_test():

	for npts in (7, 15, 31, 61):
		PENPAN = []
		for ibasin in xrange(0, 10):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))

			for istation in good_stations[ibasin]:
				print ibasin, istation
				index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
				# Read all the necessary input into a dataframe
				input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
				input = pd.DataFrame.from_dict(input)
				input.index = dates

				# Remove inter-annual variability
				# INPUT = Remove_interannual_variability(input)
				# vars_penpan = ['tavg', 'tmax', 'tmin', 'p', 'ea', 'wind', 'sun', 'lat', 'elev']  # 'tc'

				# INPUT = Remove_shortterm_variability(input, vars_penpan[:-2], npts)

				# tair
				INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg'], npts)

				# tair+wind
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg', 'wind'], npts)

				# tair+wind+humidity
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg', 'wind', 'ea'], npts)

				# tair+wind+humidity+pressure
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg', 'wind', 'ea', 'p'], npts)

				# tair+wind+humidity+pressure+sun
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg', 'wind', 'ea', 'p'], npts)

				# tair+sun
				# INPUT = Remove_shortterm_variability(input, ['tmax', 'tmin', 'tavg', 'wind', 'ea'], npts)

				INPUT['doy'] = doys
				INPUT['lat'] = geoinfo[index, 1]
				INPUT['elev'] = geoinfo[index, 3]
				### Calculate Epan
				res = Data(INPUT, 'sunhours') #, npts)

				PENPAN.append(res.penpan)

		pe_model = array(PENPAN)
	return

# Calculate_Ep_daily_smooth_test()
# exit()

def Calculate_Ep_daily_ensemble():
	PENPAN = []
	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			print ibasin, istation
			index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
			# Read all the necessary input into a dataframe
			input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
			input = pd.DataFrame.from_dict(input)
			input.index = dates

			# Run the permutation program with multi-ensemble
			penpan_ens = []
			for i in xrange(10):  # Set # of ensemble to 10
				# Monthy permutation
				# INPUT = Permute_each_month_station(input)
				# Weekly permutation
				INPUT = Sample_each_week_station(input)

				INPUT['doy'] = doys
				INPUT['lat'] = geoinfo[index, 1]
				INPUT['elev'] = geoinfo[index, 3]
				### Calculate Epan
				res = Data(INPUT, 'sunhours')

				penpan_ens.append(res.penpan)

			# Collect all the ensembles
			PENPAN.append(penpan_ens)
			del penpan_ens

	pe_model = array(PENPAN)
	# pe_model.dump('%s/pe_mod_penpan_monthpermute_ens10_good_stations' %(workspace))
	pe_model.dump('%s/pe_mod_penpan_weeksample_ens10_good_stations' %(workspace))

	return

# Calculate_Ep_daily_ensemble()
# exit()

def Calculate_Ep_daily_samplewindow():

	for window in (7, 31, 91):
		PENPAN = []

		for ibasin in xrange(0, 10):
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			for istation in good_stations[ibasin]:
				print ibasin, istation
				index = np.where(geoinfo[:, 0]==data['station_name'][0, istation])[0]
				# Read all the necessary input into a dataframe
				input = {vars_penpan[i]: Gapfill(data[v][0, istation][0:tstep].flatten()) for i, v in enumerate(vars_penpan[:-2])}
				input = pd.DataFrame.from_dict(input)
				input.index = dates

				# Run the permutation program with multi-ensemble
				penpan_ens = []
				for i in xrange(10):  # Set # of ensemble to 10
					# Monthy permutation
					INPUT = Sample_station_window(input, window)

					INPUT['doy'] = doys
					INPUT['lat'] = geoinfo[index, 1]
					INPUT['elev'] = geoinfo[index, 3]
					### Calculate Epan
					res = Data(INPUT, 'sunhours')

					penpan_ens.append(res.penpan)

				# Collect all the ensembles
				PENPAN.append(penpan_ens)
				del penpan_ens

		pe_model = array(PENPAN)
		pe_model.dump('%s/pe_mod_penpan_sample%sd_ens10_good_stations' %(workspace, window))

	return

# Calculate_Ep_daily_samplewindow()
# exit()

def Coherence_obs_permute():

	"Compare the observed pan with the modelled pan"
	obs = load('%s/pe_mod_penpan_good_stations' %(workspace))
	# permute = load('%s/pe_mod_penpan_monthpermute_ens10_good_stations' %(workspace))
	permute = load('%s/pe_mod_penpan_weeksample_ens10_good_stations' %(workspace))
	cohere = []
	for ist in xrange(0, 228):
		print ist
		cohere.append(array([FFT.Coherence(obs[ist, :], permute[ist, i, :], sampling_frequency, 'linear')[1] for i in xrange(10)])) #.reshape(1, 5, nf))
	cohere = array(cohere)
	# cohere.dump('%s/coherence_penpan_monthpermute_ens10_good_stations' %(workspace))
	cohere.dump('%s/coherence_penpan_weeksample_ens10_good_stations' %(workspace))

	return

# Coherence_obs_permute()
# exit()

def Coherence_obs_remove_variability_variable():

	"Compare the effects of removing shortterm variability in each variable"
	obs = load('%s/pe_mod_penpan_good_stations' %(workspace))
	# exp = load('%s/pe_mod_penpan_removeiav_good_stations' %(workspace))

	npt = 7
	st, ed = chunk_st_ed(npt)
	vars = ['tair', 'ea', 'rnet', 'wind']
	for var in vars[1:]:
		exp = load('%s/pe_mod_penpan_removeshort%s_%s_good_stations' %(workspace, npt, var))

		cohere = []
		for ist in xrange(0, 228):
			print ist
			# cohere.append(FFT.Coherence(obs[ist, npt:-npt], exp[ist, npt:-npt], sampling_frequency, 'linear')[1])  # for removeshort
			cohere.append(FFT.Coherence(obs[ist, st:-ed], exp[ist, st:-ed], sampling_frequency, 'linear')[1])  # for removeshort

		cohere = array(cohere)
		cohere.dump('%s/coherence_penpan_removeshort%sd_%s_good_stations' %(workspace, npt, var))

	return

# Coherence_obs_remove_variability_variable()
# exit()

def Coherence_obs_remove_variability_window():

	"Test how the window of moving averaging affect the results"
	obs = load('%s/pe_mod_penpan_good_stations' %(workspace))

	for npt in (7, 15, 30):
		st, ed = chunk_st_ed(npt)
		# exp = load('%s/pe_mod_penpan_removeiav_good_stations' %(workspace))
		# exp = load('%s/pe_mod_penpan_removeshort%s_tair_good_stations' %(workspace, npt))
		# exp = load('%s/pe_mod_penpan_removeshort%s_wind_good_stations' %(workspace, npt))
		# exp = load('%s/pe_mod_penpan_removeshort%s_rnet_good_stations' %(workspace, npt))
		# exp = load('%s/pe_mod_penpan_removeshort%s_ea_good_stations' %(workspace, npt))

		cohere = []
		for ist in xrange(0, 228):
			print ist
			cohere.append(FFT.Coherence(obs[ist, st:-ed], exp[ist, st:-ed], sampling_frequency, 'linear')[1])  # for removeshort

		cohere = array(cohere)
		# cohere.dump('%s/coherence_penpan_removeshort%sd_tair_good_stations' %(workspace, npt))
		# cohere.dump('%s/coherence_penpan_removeshort%sd_wind_good_stations' %(workspace, npt))
		# cohere.dump('%s/coherence_penpan_removeshort%sd_rnet_good_stations' %(workspace, npt))
		# cohere.dump('%s/coherence_penpan_removeshort%sd_ea_good_stations' %(workspace, npt))


	# for npt in (7, 31, 91):
	# 	exp = load('%s/pe_mod_penpan_sample%sd_ens10_good_stations' %(workspace, npt))
	# 	cohere = []
	# 	for ist in xrange(0, 228):
	# 		print ist
	# 		cohere.append(array([FFT.Coherence(obs[ist, :], exp[ist, i, :], sampling_frequency, 'linear')[1] for i in xrange(10)]))
	# 	cohere = array(cohere)
		# cohere.dump('%s/coherence_penpan_sample%sd_good_stations' %(workspace, npt))

	return

# Coherence_obs_remove_variability_window()
# exit()


def Plot_PSD_obs_remove_variability_window():

	"Compare the observed pan with the modelled pan"
	freq = Coherence_Frequency()

	fig, ax = plt.subplots(figsize=(8, 4))

	data = []
	for npt in (7, 15, 30):
		# exp = load('%s/pe_mod_penpan_removeiav_good_stations' %(workspace))
		exp = load('%s/pe_mod_penpan_removeshort%s_tair_good_stations' %(workspace, npt))

		psd = []
		for ist in xrange(0, 228):
			print ist
			psd.append(FFT.Power_Spectrum(exp[ist, npt-1:], sampling_frequency, 'linear')[1])
			# psd.append(FFT.Power_Spectrum(obs[ist, npt:-npt], sampling_frequency, 'linear')[1])

		psd = array(psd)
		data.append(mean(psd, axis=0))
	data = array(data)
	Plotting.CoherenceWindowPlot(ax, data, sampling_frequency, freq)

	plt.show()

	return

# Plot_PSD_obs_remove_variability_window()
# exit()


def Plot_crossspectrum_obs_remove_variability_window():

	"Compare the observed pan with the modelled pan"
	freq = Coherence_Frequency()

	obs = load('%s/pe_mod_penpan_good_stations' %(workspace))
	fig, ax = plt.subplots(figsize=(8, 4))

	data = []
	for npt in (7, 15, 30):
		# exp = load('%s/pe_mod_penpan_removeiav_good_stations' %(workspace))
		exp = load('%s/pe_mod_penpan_removeshort%s_tair_good_stations' %(workspace, npt))

		cross = []
		for ist in xrange(0, 228):
			print ist
			cross.append(FFT.CrossPowerSpectrum(obs[ist, npt-1:], exp[ist, npt-1:], sampling_frequency, 'linear')[1])
			# cross.append(FFT.CrossPowerSpectrum(obs[ist, npt:-npt], exp[ist, npt:-npt], sampling_frequency, 'linear')[1])

		cross = array(cross)
		data.append(mean(cross, axis=0))
	data = array(data)
	Plotting.CoherenceWindowPlot(ax, data, sampling_frequency, freq)

	plt.show()

	return
# Plot_crossspectrum_obs_remove_variability_window()
# exit()

def Plot_Coherence_Average():

	fig = plt.figure(figsize=(9, 5))
	freq = Coherence_Frequency()
	# cohere = load('%s/coherence_penpan_monthpermute_ens10_good_stations' %(workspace))

	# cohere = load('%s/coherence_penpan_removeiav_good_stations' %(workspace))
	# cohere = load('%s/coherence_penpan_removeshort15_good_stations' %(workspace))
	# cohere = load('%s/coherence_penpan_removeshort7_good_stations' %(workspace))
	# cohere = load('%s/coherence_penpan_removeshort15_sun_good_stations' %(workspace))
	# cohere = load('%s/coherence_penpan_removeshort15_ea_good_stations' %(workspace))
	# cohere = load('%s/coherence_penpan_removeshort15_wind_good_stations' %(workspace))
	cohere = load('%s/coherence_penpan_removeshort15_tair_good_stations' %(workspace))

	# cohere = load('%s/coherence_penpan_weeksample_ens10_good_stations' %(workspace))

	# for all stations average
	ax = fig.add_subplot(1, 1, 1)

	# plot all ensemble: found they can be averaged
	# Plotting.CoherenceEnsemblePlot(ax, mean(cohere, axis=0), sampling_frequency, freq, 'Average')

	# plot ensemble mean and then all aridity
	# cohere_avg = mean(cohere, axis=1)

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)
	Plotting.CoherenceAridityPlot(ax, cohere, index_DI, sampling_frequency, freq, '', '')  # for ensemble: cohere_avg
	# ax.legend(loc=2, fontsize=15)
	plt.show()
	# savefig('%s/coh_penpan_removeiav_good_stations.tif' %(figdir), dpi=200)
	# savefig('%s/coh_penpan_removeshort15_good_stations.tif' %(figdir), dpi=200)
	# savefig('%s/coh_penpan_removeshort7_good_stations.tif' %(figdir), dpi=200)
	# savefig('%s/coh_penpan_weeksample_ens10_good_stations.tif' %(figdir), dpi=200)

	return

# Plot_Coherence_Average()
# exit()

def Plot_Coherence_multiple():

	freq = Coherence_Frequency()

	rnet = load('%s/coherence_penpan_removeshort30d_rnet_good_stations' %(workspace))
	wind = load('%s/coherence_penpan_removeshort30d_wind_good_stations' %(workspace))
	tair = load('%s/coherence_penpan_removeshort30d_tair_good_stations' %(workspace))
	ea = load('%s/coherence_penpan_removeshort30d_ea_good_stations' %(workspace))
	all = load('%s/coherence_penpan_removeshort15_good_stations' %(workspace))
	sample = load('%s/coherence_penpan_weeksample_ens10_good_stations' %(workspace))
	# plot ensemble mean and then all aridity
	sample_avg = mean(sample, axis=1)

	vars = [tair, rnet, ea, wind, all, sample_avg]
	labels = ('a', 'b', 'c', 'd', 'e', 'f')
	varnames = ['Tair', 'Solar', r'$e_a$', 'Wind', 'All', 'Randomized']

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	fig, axes = plt.subplots(3, 2, figsize=(11, 10))
	for i, var in enumerate(vars):
		Plotting.CoherenceAridityPlot(axes[i/2, i%2], var, index_DI, sampling_frequency, freq, labels[i], varnames[i])  # for ensemble: cohere_avg
		if i == 5:
			axes[i/2, i%2].legend(loc=2, fontsize=13)
	fig.tight_layout()
	plt.show()
	# savefig('%s/coh_penpan_removeshort15_multi_good_stations.tif' %(figdir), dpi=200)

	return

# Plot_Coherence_multiple()
# exit()


def Plot_Coherence_multiwindow():

	freq = Coherence_Frequency()
	fig, ax = plt.subplots(2, 2, figsize=(9, 6), sharey=True, sharex=True)
	varnames = ['tair', 'rnet', 'wind', 'ea']
	nums = ['(a) ', '(b) ', '(c) ', '(d) ']
	labels = [r'$T_a$', r'$R_n$', r'$u_2$', r'$e_a$']

	for i in range(0,4):

		vars = [load('%s/coherence_penpan_removeshort%sd_%s_good_stations' %(workspace, npts, varnames[i])) for npts in (7, 15, 30)]
		coh = array([mean(var, axis=0) for var in vars])

		# vars = [load('%s/coherence_penpan_sample%sd_good_stations' %(workspace, npt)) for npt in (7, 31, 91)]
		# coh = array([mean(mean(var, axis=1), axis=0) for var in vars])

		Plotting.CoherenceWindowPlot(ax[i/2, i%2], 1-coh, nums[i]+labels[i], sampling_frequency, freq)

	ax[1,1].legend(loc='upper right', frameon=False, fontsize=14)
	fig.tight_layout()
	plt.show()
	# savefig('%s/figS1_coh_penpan_removeshort7-30_4var_good_stations.tif' %(figdir), dpi=300)

	return

# Plot_Coherence_multiwindow()
# exit()


def Plot_Coherence_samplewindow():

	freq = Coherence_Frequency()
	names = ['(a) Sampling window', '(b) Dryness']
	fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5))

	# The first figure
	vars = [load('%s/coherence_penpan_sample%sd_good_stations' %(workspace, npt)) for npt in (91, 31, 7)]
	coh = array([mean(mean(var, axis=1), axis=0) for var in vars])
	i = 0
	Plotting.CoherenceWindow2Plot(axes[i], coh, sampling_frequency, freq, names[i])
	axes[i].legend(loc=2, fontsize=14)

	sample = load('%s/coherence_penpan_sample7d_good_stations' %(workspace))
	# plot ensemble mean and then all aridity
	sample_avg = mean(sample, axis=1)

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	# for i, var in enumerate(vars):
	i = 1
	Plotting.CoherenceAridityPlot(axes[i], sample_avg, index_DI, sampling_frequency, freq, names[i], '')  # for ensemble: cohere_avg
	axes[i].legend(loc=2, fontsize=12)

	fig.tight_layout()
	plt.show()
	# savefig('%s/coh_penpan_sample_window_dryness_good_stations.tif' %(figdir), dpi=300)

	return

# Plot_Coherence_samplewindow()
# exit()



def Plot_Coherence_loss_samplewindow():

	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	names = ['(a) Sampling window', '(b) Dryness']
	fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5))

	# The first figure
	vars = [load('%s/coherence_penpan_sample%sd_good_stations' %(workspace, npt)) for npt in (91, 31, 7)]
	coh = array([mean(mean(var, axis=1), axis=0) for var in vars])
	coh_point = array([interp1d(freq_ts[:], 1-coh)(day)[()] for day in [7, 30, 90, 120, 180, 365]])
	i = 0
	Plotting.CoherenceWindowPointPlot(axes[i], coh_point, names[i])
	axes[i].legend(loc=3, frameon=False, fontsize=14)

	# The second figure
	sample = load('%s/coherence_penpan_sample7d_good_stations' %(workspace))
	sample_avg = mean(sample, axis=1)  # plot ensemble mean and then all aridity
	sample_point = array([interp1d(freq_ts[:], 1-sample_avg)(day)[()] for day in [7, 30, 90, 120, 180, 365]])

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	i = 1
	Plotting.CoherenceWindowAridityPointPlot(axes[i], sample_point, index_DI, names[i])  # for ensemble: cohere_avg
	axes[i].legend(loc=3, frameon=False, fontsize=12)

	fig.tight_layout()
	plt.show()
	# savefig('%s/fig9_coh_penpan_sample_window_dryness_point_good_stations.tif' %(figdir), dpi=300)

	return

# Plot_Coherence_loss_samplewindow()
# exit()



def Plot_Coherence_loss_samplewindow_extra():

	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	names = ['(a) Sampling window', '(b) Dryness (7d window)', '(c) Dryness (30d window)', '(d) Dryness (90d window)']
	fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.5))

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)


	# # The first figure
	vars = [load('%s/coherence_penpan_sample%sd_good_stations' %(workspace, npt)) for npt in (91, 31, 7)]
	coh = array([mean(mean(var, axis=1), axis=0) for var in vars])
	coh_point = array([interp1d(freq_ts[:], 1-coh)(day)[()] for day in [7, 30, 90, 120, 180, 365]])
	i = 0
	Plotting.CoherenceWindowPointPlot(axes[0,0], coh_point, names[i])
	axes[0,0].legend(loc=3, frameon=False, fontsize=14)

	# The second figure
	sample = load('%s/coherence_penpan_sample7d_good_stations' %(workspace))
	sample_avg = mean(sample, axis=1)  # plot ensemble mean and then all aridity
	sample_point = array([interp1d(freq_ts[:], 1-sample_avg)(day)[()] for day in [7, 30, 90, 120, 180, 365]])
	i = 1
	Plotting.CoherenceWindowAridityPointPlot(axes[0,1], sample_point, index_DI, names[i])  # for ensemble: cohere_avg

	# The 3rd figure
	sample = load('%s/coherence_penpan_sample31d_good_stations' %(workspace))
	sample_avg = mean(sample, axis=1)  # plot ensemble mean and then all aridity
	sample_30 = array([interp1d(freq_ts[:], 1-sample_avg)(day)[()] for day in [7, 30, 90, 120, 180, 365]])
	i = 2
	Plotting.CoherenceWindowAridityPointPlot(axes[1,0], sample_30, index_DI, names[i])

	# The 4th figure
	sample = load('%s/coherence_penpan_sample91d_good_stations' %(workspace))
	sample_avg = mean(sample, axis=1)  # plot ensemble mean and then all aridity
	sample_90 = array([interp1d(freq_ts[:], 1-sample_avg)(day)[()] for day in [7, 30, 90, 120, 180, 365]])

	i = 3
	Plotting.CoherenceWindowAridityPointPlot(axes[1,1], sample_90, index_DI, names[i])  # for ensemble: cohere_avg
	axes[1,1].legend(loc=3, frameon=False, fontsize=12)

	fig.tight_layout()
	plt.show()
	# savefig('%s/fig10_coh_penpan_sample_window_dryness_point_good_stations_all.tif' %(figdir), dpi=300)

	return

# Plot_Coherence_loss_samplewindow_extra()
# exit()

def Plot_Coherence_grid_test():

	freq = Coherence_Frequency()

	rnet = load('%s/coherence_penpan_removeshort30d_rnet_good_stations' %(workspace))
	wind = load('%s/coherence_penpan_removeshort30d_wind_good_stations' %(workspace))
	tair = load('%s/coherence_penpan_removeshort30d_tair_good_stations' %(workspace))
	ea = load('%s/coherence_penpan_removeshort30d_ea_good_stations' %(workspace))

	vars = [tair, rnet, wind, ea]

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	matrix = np.zeros((4, 2, 3))
	freq_ts = sampling_frequency/freq
	for iv, var in enumerate(vars):
		res = array([[interp1d(freq_ts[:], var[istation, :])(day)[()] for day in [7, 15, 30]] for istation in xrange(0, var.shape[0])])
		wet = vstack((res[index_DI[0], :], res[index_DI[1], :]))
		dry = vstack((res[index_DI[2], :], res[index_DI[3], :]))
		matrix[iv, 0, :] = mean(wet, axis=0)
		matrix[iv, 1, :] = mean(dry, axis=0)

	# rearrange the matrix into table
	imshow_data = np.zeros((4, 6))
	for id in range(0,2):
		for it in range(0,3):
			for iv in range(0,4):
				row = id*2 + iv/2
				col = it*2 + iv%2
				print row, col, matrix[iv, id, it]
				imshow_data[row, col] = 1 - matrix[iv, id, it]
	plt.imshow(imshow_data, cmap='hot_r', interpolation='nearest')
	plt.show()
	# savefig('%s/coh_penpan_removeshort15_multi_good_stations.tif' %(figdir), dpi=200)

	return

# Plot_Coherence_grid_test()
# exit()


def Plot_Coherence_grid_average(fig, ax, npt, title):
	"try my best to extract the information"

	rnet = load('%s/coherence_penpan_removeshort%sd_rnet_good_stations' %(workspace, npt))
	wind = load('%s/coherence_penpan_removeshort%sd_wind_good_stations' %(workspace, npt))
	tair = load('%s/coherence_penpan_removeshort%sd_tair_good_stations' %(workspace, npt))
	ea = load('%s/coherence_penpan_removeshort%sd_ea_good_stations' %(workspace, npt))

	vars = [tair, rnet, wind, ea]
	varnames = [r'$T_a$', r'$R_n$', r'$u_2$', r'$e_a$']

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	# make a 3D matrix for summary table
	matrix = np.zeros((4, 2, 2))
	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	for iv, var in enumerate(vars):
		index_week = (freq_ts>=2) & (freq_ts<=7)
		index_month = (freq_ts>7) & (freq_ts<=30)
		res = vstack((mean(var[:, index_week], axis=1), mean(var[:, index_month], axis=1)))
		wet = hstack((res[:, index_DI[0]], res[:, index_DI[1]]))
		dry = hstack((res[:, index_DI[2]], res[:, index_DI[3]]))
		matrix[iv, 0, :] = mean(wet, axis=1)
		matrix[iv, 1, :] = mean(dry, axis=1)

	# rearrange the matrix into table
	imshow_data = np.zeros((4, 4)); imshow_label = np.empty((4, 4), dtype=int)
	for id in range(0,2):
		for it in range(0,2):
			for iv in range(0,4):
				row = id*2 + iv/2
				col = it*2 + iv%2
				print row, col
				imshow_data[row, col] = 1 - matrix[iv, id, it]  # This is the influence 1-MSC
				imshow_label[row, col] = iv

	im = ax.imshow(imshow_data, vmax=0.45, vmin=0.0, cmap='YlOrRd', interpolation='nearest')

	# colorbar
	if npt == 30:
		# set up the axis
		cax = fig.add_axes([0.91, 0.12, 0.02, 0.78])
		cb = fig.colorbar(im, cax)  # adjust the size
		# cb = ax.colorbar(im) #, fraction=0.046, pad=0.04)  # magic number!!!!!
		cb.ax.tick_params(labelsize=14)   # change the colorbar fontsize

	# Text portion
	ind_array = np.arange(0, 4, 1)
	x, y = meshgrid(ind_array, ind_array)
	for xloc, yloc in zip(x.flatten(), y.flatten()):
		ax.text(xloc, yloc, varnames[imshow_label[yloc, xloc]], va='center', ha='center', fontsize=20)

	# Two separate lines
	ax.plot([1.5, 1.5], [-0.5, 3.5], c='black', linewidth=2)
	ax.plot([-0.5, 3.5],[1.5, 1.5], c='black', linewidth=2)
	# x y label
	fig.subplots_adjust(bottom=0.12)
	ax.set_xticks((0.5, 2.5))
	ax.set_yticks((0.5, 2.5))
	ax.set_xticklabels(['Weekly cycle\n(2-7d)', 'Monthly cycle\n(7-30d)'], fontsize=16)
	ax.set_yticklabels(["Wet" "\n" r"($\phi$<4)", "Dry" "\n" r"($\phi$>4)"], fontsize=16) # treat this as special string
	ax.set_title(title, fontsize=16)
	# savefig('%s/coh_grid_average_scale_climate_4var_removeshort30.tif' %(figdir), dpi=300)

	return

def Plot_Coherence_grid_average_multiple():
	fig, ax = plt.subplots(1, 2, figsize=(12, 5))
	Plot_Coherence_grid_average(fig, ax[0], 7, '(a) Window = 7d')
	Plot_Coherence_grid_average(fig, ax[1], 30, '(b) Window = 30d')

	plt.show()
	# savefig('%s/coh_grid_average_scale_climate_4var_removeshort7-30.tif' %(figdir), dpi=300)
	return
# Plot_Coherence_grid_average_multiple()
# exit()



def Plot_Coherence_grid_average_update(fig, ax, npt, title):
	"try my best to extract the information"

	rnet = load('%s/coherence_penpan_removeshort%sd_rnet_good_stations' %(workspace, npt))
	wind = load('%s/coherence_penpan_removeshort%sd_wind_good_stations' %(workspace, npt))
	tair = load('%s/coherence_penpan_removeshort%sd_tair_good_stations' %(workspace, npt))
	ea = load('%s/coherence_penpan_removeshort%sd_ea_good_stations' %(workspace, npt))

	vars = [tair, rnet, wind, ea]
	varnames = [r'$T_a$', r'$R_n$', r'$u_2$', r'$e_a$']

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	# make a 3D matrix for summary table
	matrix = np.zeros((4, 2, 2))
	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	for iv, var in enumerate(vars):
		index_week = (freq_ts>=2) & (freq_ts<=7)
		index_month = (freq_ts>7) & (freq_ts<=30)
		res = vstack((mean(var[:, index_week], axis=1), mean(var[:, index_month], axis=1)))
		wet = hstack((res[:, index_DI[0]], res[:, index_DI[1]]))
		dry = hstack((res[:, index_DI[2]], res[:, index_DI[3]]))
		matrix[iv, 0, :] = mean(wet, axis=1)
		matrix[iv, 1, :] = mean(dry, axis=1)

	# rearrange the matrix into table
	imshow_data = np.zeros((4, 4)); imshow_label = np.empty((4, 4), dtype=int)
	for id in range(0,2):
		for it in range(0,2):
			for iv in range(0,4):
				row = id*2 + iv/2
				col = it*2 + iv%2
				print row, col
				imshow_data[row, col] = 1 - matrix[iv, id, it]  # This is the influence 1-MSC
				imshow_label[row, col] = iv

	im = ax.imshow(imshow_data, vmax=0.45, vmin=0.0, cmap='YlOrRd', interpolation='nearest')

	# colorbar
	# set up the axis
	# cax = fig.add_axes([0.91, 0.12, 0.02, 0.78])
	# cb = fig.colorbar(im, cax)  # adjust the size
	cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)  # magic number!!!!!
	cb.ax.tick_params(labelsize=14)   # change the colorbar fontsize

	# Text portion
	ind_array = np.arange(0, 4, 1)
	x, y = meshgrid(ind_array, ind_array)
	for xloc, yloc in zip(x.flatten(), y.flatten()):
		ax.text(xloc, yloc, varnames[imshow_label[yloc, xloc]], va='center', ha='center', fontsize=20)

	# Two separate lines
	ax.plot([1.5, 1.5], [-0.5, 3.5], c='black', linewidth=2)
	ax.plot([-0.5, 3.5],[1.5, 1.5], c='black', linewidth=2)
	# x y label
	fig.subplots_adjust(bottom=0.12)
	ax.set_xticks((0.5, 2.5))
	ax.set_yticks((0.5, 2.5))
	ax.set_xticklabels(['Weekly cycle\n(2-7d)', 'Monthly cycle\n(7-30d)'], fontsize=14)
	if npt==7:
		ax.set_yticklabels(["Wet" "\n" r"($\phi$<4)", "Dry" "\n" r"($\phi$>4)"], fontsize=16) # treat this as special string
	else:
		ax.set_yticklabels(["", ""], fontsize=16) # treat this as special string

	ax.set_title(title, fontsize=16)
	# savefig('%s/coh_grid_average_scale_climate_4var_removeshort30.tif' %(figdir), dpi=300)

	return


def Coherence_grid_average(npt):
	"try my best to extract the information"

	rnet = load('%s/coherence_penpan_removeshort%sd_rnet_good_stations' %(workspace, npt))
	wind = load('%s/coherence_penpan_removeshort%sd_wind_good_stations' %(workspace, npt))
	tair = load('%s/coherence_penpan_removeshort%sd_tair_good_stations' %(workspace, npt))
	ea = load('%s/coherence_penpan_removeshort%sd_ea_good_stations' %(workspace, npt))

	vars = [tair, rnet, wind, ea]

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)

	# make a 3D matrix for summary table
	matrix = np.zeros((4, 2, 2))
	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	for iv, var in enumerate(vars):
		index_week = (freq_ts>=2) & (freq_ts<=7)
		index_month = (freq_ts>7) & (freq_ts<=30)
		res = vstack((mean(var[:, index_week], axis=1), mean(var[:, index_month], axis=1)))
		wet = hstack((res[:, index_DI[0]], res[:, index_DI[1]]))
		dry = hstack((res[:, index_DI[2]], res[:, index_DI[3]]))
		matrix[iv, 0, :] = mean(wet, axis=1)
		matrix[iv, 1, :] = mean(dry, axis=1)

	# rearrange the matrix into table
	imshow_data = np.zeros((4, 4)); imshow_label = np.empty((4, 4), dtype=int)
	for id in range(0,2):
		for it in range(0,2):
			for iv in range(0,4):
				row = id*2 + iv/2
				col = it*2 + iv%2
				print row, col
				imshow_data[row, col] = 1 - matrix[iv, id, it]  # This is the influence 1-MSC
				imshow_label[row, col] = iv


	return imshow_data, imshow_label


from matplotlib.colors import Normalize
class MidpointNormalize(Normalize):

	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))

def Plot_difference(fig, ax, label, data1, data2, title):
	varnames = [r'$T_a$', r'$R_n$', r'$u_2$', r'$e_a$']
	# my_cmap = cm.YlOrRd
	# my_cmap.set_under('white')

	norm = MidpointNormalize(vmax=0.4, midpoint=0, vmin=-0.1)

	im = ax.imshow(data2-data1, cmap='bwr', norm=norm, interpolation='nearest')

	# colorbar
	# set up the axis
	cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)  # magic number!!!!!
	cb.ax.tick_params(labelsize=14)   # change the colorbar fontsize

	# Text portion
	ind_array = np.arange(0, 4, 1)
	x, y = meshgrid(ind_array, ind_array)
	for xloc, yloc in zip(x.flatten(), y.flatten()):
		ax.text(xloc, yloc, varnames[label[yloc, xloc]], va='center', ha='center', fontsize=20)

	# Two separate lines
	ax.plot([1.5, 1.5], [-0.5, 3.5], c='black', linewidth=2)
	ax.plot([-0.5, 3.5],[1.5, 1.5], c='black', linewidth=2)
	# x y label
	fig.subplots_adjust(bottom=0.12)
	ax.set_xticks((0.5, 2.5))
	ax.set_yticks((0.5, 2.5))
	ax.set_xticklabels(['Weekly cycle\n(2-7d)', 'Monthly cycle\n(7-30d)'], fontsize=14)
	ax.set_yticklabels(["", ""], fontsize=16) # treat this as special string
	ax.set_title(title, fontsize=16)

	return

def Plot_Coherence_grid_average_multiple_update():
	fig, ax = plt.subplots(1, 3, figsize=(12, 4))
	Plot_Coherence_grid_average_update(fig, ax[0], 7, '(a) Window = 7d')
	Plot_Coherence_grid_average_update(fig, ax[1], 30, '(b) Window = 30d')
	data7, labels = Coherence_grid_average(7)
	data30, labels = Coherence_grid_average(30)
	Plot_difference(fig, ax[2], labels, data7, data30, '(c) 30d - 7d')
	fig.tight_layout()
	plt.show()
	# savefig('%s/Fig9_coh_grid_average_scale_climate_4var_diff_removeshort7-30.tif' %(figdir), dpi=300)
	return
# Plot_Coherence_grid_average_multiple_update()
# exit()


def Print_Coherence_Average():

	cohere = []
	freq = Coherence_Frequency()
	freq_ts = sampling_frequency/freq
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_5model_good_station_%s' %(workspace, basinlongs[ibasin]))
		cohere.append(cohere_basin)

	# for all stations average
	res = [interp1d(freq_ts[:], mean(vstack(cohere), axis=0))(day)[()] for day in [250]]
	print res
	return

# Print_Coherence_Average()
# exit()




