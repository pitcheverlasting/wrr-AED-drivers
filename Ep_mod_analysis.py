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
figdir = '/Users/pitch/Google Drive/Figure/pan_spectral/201705'
# figdir = '/home/water5/lpeng/Figure/pan_spectral/201705'

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
doys = vstack([dates[i].timetuple().tm_yday for i in xrange(0, tstep)])  # julian day

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

from matplotlib import rc
rc('font', family='serif') #')'Times New Roman'


def Coherence_obs_5mod():

	"Compare the observed pan with the modelled pan"
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	penpan = load('%s/pe_mod_%s_good_stations' %(workspace, 'penpan'))
	penman = load('%s/pe_mod_%s_good_stations' %(workspace, 'penman'))
	priest = load('%s/pe_mod_%s_good_stations' %(workspace, 'priestley-taylor'))
	hamon = load('%s/pe_mod_%s_good_stations' %(workspace, 'hamon'))
	turc = load('%s/pe_mod_%s_good_stations' %(workspace, 'turc'))
	pedata = [penpan, penman, priest, hamon, turc]

	ist = 0
	for ibasin in xrange(0, 10):
		print ibasin
		cohere_obs_basin = []
		for istation in good_stations[ibasin]:
			cohere_obs_basin.append(vstack([FFT.Coherence(pan_obs_gapfill[ist, :], pe[ist, :], sampling_frequency, 'linear')[1] for pe in pedata]).reshape(1, 5, nf))
			ist = ist + 1
		cohere_obs_basin = vstack(cohere_obs_basin)
		cohere_obs_basin.dump('%s/coherence_obs_5model_good_station_%s' %(workspace, basinlongs[ibasin]))

	return

# Coherence_obs_5mod()
# exit()

def Coherence_mod_5var():

	"Compare the modelled pan with all climate drivers"
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	penpan = load('%s/pe_mod_%s_good_stations' %(workspace, 'penpan'))
	penman = load('%s/pe_mod_%s_good_stations' %(workspace, 'penman'))
	priest = load('%s/pe_mod_%s_good_stations' %(workspace, 'priestley-taylor'))
	hamon = load('%s/pe_mod_%s_good_stations' %(workspace, 'hamon'))
	turc = load('%s/pe_mod_%s_good_stations' %(workspace, 'turc'))
	pedata = [pan_obs_gapfill, penpan, penman, priest, hamon, turc]

	ist = 0
	for ibasin in xrange(0, 10):
		cohere_obs_basin = []
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			input.insert(1, rn[ist, :].flatten())
			cohere_obs_basin.append(vstack([vstack([FFT.Coherence(v, pe[ist, :], sampling_frequency, 'linear')[1] for v in input]).reshape(1, 5, nf) for pe in pedata]).reshape(1, 6, 5, nf))
			ist = ist + 1
		cohere_obs_basin = vstack(cohere_obs_basin)
		cohere_obs_basin.dump('%s/coherence_obs_5model_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		del data

	return

# Coherence_mod_5var()
# exit()


def Plot_Coherence_Point():
	"Newest"
	msc = []
	fig = plt.figure(figsize=(9, 4))
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_5model_good_station_%s' %(workspace, basinlongs[ibasin]))
		# print cohere_basin.shape;exit()
		freq_ts = sampling_frequency/freq
		msc.append(array([[interp1d(freq_ts[:], cohere_basin[istation, :, :])(day)[()] for day in [7, 30, 90, 120, 180, 365]] for istation in xrange(0, cohere_basin.shape[0])]))
	msc = vstack(msc)
	ax = fig.add_subplot(1, 1, 1)
	Plotting.CoherenceModelPointPlot(ax, msc)

	# fig.tight_layout()
	plt.show()
	# savefig('%s/fig7_coh_obs_all5model_point_good_station.tif' % figdir, dpi=300)

	return

# Plot_Coherence_Point()
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


def Coherence_Groupby_Aridity():

	msc = []
	arid = []
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_5model_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		msc.append(cohere_basin)
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	msc = vstack(msc)

	index_DI = msc_groupby_DI(arid)
	freq = Coherence_Frequency()

	drivername = ['Tair','Rn','Wind','Humidity','VPD']
	pe_modelnames = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	labels = ('a', 'b', 'c', 'd', 'e', 'f')
	# varorder = (4, 1, 0, 3, 2)
	for ivar in xrange(0, 5):
		fig = plt.figure(figsize=(14, 12))
		for imod in xrange(0, 6):
			ax = fig.add_subplot(3, 2, imod+1)
			Plotting.CoherenceAridityPlot(ax, msc[:, imod, ivar, :], index_DI, sampling_frequency, freq, labels[imod], pe_modelnames[imod])
			# if imod == 4:
			# 	ax.legend(bbox_to_anchor=(1.35, 0.5), loc='center left', fontsize=20)

		fig.suptitle('%s' % drivername[ivar], fontsize=24)
		plt.subplots_adjust(top=0.92)
		# fig.tight_layout()
		plt.show()
		# savefig('%s/coh_mod_%s_aridity_good_station.tif' %(figdir, drivername[ivar]), dpi=200)

	return

# Coherence_Groupby_Aridity()
# exit()


def Coherence_Groupby_Model_timescales_dry_wet():
	"Even newer for Fig.8: remove wind and put two dryness together"
	msc = []
	arid = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_5model_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		freq_ts = sampling_frequency/freq
		msc.append(array([[interp1d(freq_ts[:], cohere_basin[istation, :, :, :])(day)[()] for day in [7, 30, 90, 120, 180, 365]] for istation in xrange(0, cohere_basin.shape[0])]))
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	msc = vstack(msc)

	index_DI = msc_groupby_DI(arid)
	drivername = [r'$T_a$',r'$R_n$',r'$u_2$',r'$e_a$','VPD']

	index_model = [0, 1, 2, 3, 4, 5]  # include all
	labels = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j')
	arid_class = [r'Wet ($\phi$<2)', r'Dry ($\phi$>4)']

	fig, ax = plt.subplots(5, 2, figsize=(9, 10))

	# plot_loc = range(0, 8, 2)
	for i, ivar in enumerate((4,1,0,3,2)):
		print i, ivar
		Plotting.CoherencePointSingleAridityPlot(ax[i, 0], msc[:, :, :, ivar], index_DI[0], index_model, labels[i], drivername[ivar])
	for i, ivar in enumerate((4,1,0,3,2)):
		Plotting.CoherencePointSingleAridityPlot(ax[i, 1], msc[:, :, :, ivar], hstack((index_DI[2], index_DI[3])), index_model, labels[i+5], drivername[ivar])
	ax[0,0].set_title(arid_class[0])
	ax[0,1].set_title(arid_class[1])
	ax[4,0].legend(loc='upper left', ncol=2, frameon=False, fontsize=11.5)

	fig.tight_layout()
	plt.show()
	# savefig('%s/fig7_coh_point_period_dry_wet_obs_all5model_5var.tif' %(figdir), dpi=300)

	return

# Coherence_Groupby_Model_timescales_dry_wet()
# exit()

def Model_obs_taylor_diagram():

	"Compare the observed pan with the modelled pan"
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	penpan = load('%s/pe_mod_%s_good_stations' %(workspace, 'penpan'))
	penman = load('%s/pe_mod_%s_good_stations' %(workspace, 'penman'))
	priest = load('%s/pe_mod_%s_good_stations' %(workspace, 'priestley-taylor'))
	hamon = load('%s/pe_mod_%s_good_stations' %(workspace, 'hamon'))
	turc = load('%s/pe_mod_%s_good_stations' %(workspace, 'turc'))

	pe_daily = [pan_obs_gapfill, penpan, penman, priest, hamon, turc]

	# daily case:
	STEP = 'daily'
	pe_data = pe_daily
	# weekly case
	# STEP = 'weekly'
	# pe_data = [daily2weekly_df(pe) for pe in pe_daily]
	# monthly case
	# STEP = 'monthly'
	# pe_data = [daily2monthly_df(pe) for pe in pe_daily]

	import skill_metrics as sm
	statnames = ['sdev', 'crmsd', 'ccoef']
	modelnames = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	intervalsCOR = np.concatenate((np.arange(0,1.0,0.2), [0.8, 0.9, 0.95, 0.99, 1]))

	# Calculate statistics for Taylor diagram
	#  The first array element corresponds to the reference series
	#  for the while the second is that for the predicted series.

	def Taylor_stat(model, obs):
		return sm.taylor_statistics(model, obs)

	def array_one_station(stat, stn):
		res = [Taylor_stat(pe[stn, :], pe_data[0][stn, :])[stat][1] for pe in pe_data[1:]]
		res.insert(0, Taylor_stat(pe_data[1][stn, :], pe_data[0][stn, :])[stat][0])
		return array(res)

	def array_all_station(stat):
		res = [Taylor_stat(pe.ravel(), pe_data[0].ravel())[stat][1] for pe in pe_data[1:]]
		res.insert(0, Taylor_stat(pe_data[1].ravel(), pe_data[0].ravel())[stat][0])
		return array(res)

	# stats = [array_one_station(stat, 5) for stat in statnames]
	# stats = [array_all_station(stat) for stat in statnames]

	def index_aridity():
		arid = []
		for ibasin in xrange(0, 10):
			arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
		arid = array(list(itertools.chain(*arid)))
		index_DI = msc_groupby_DI(arid)

		return	index_DI[0], index_DI[1], hstack((index_DI[2], index_DI[3]))

	def array_station_aridity(stat, index):
		res = [Taylor_stat(pe[index, :].ravel(), pe_data[0][index, :].ravel())[stat][1] for pe in pe_data[1:]]
		res.insert(0, Taylor_stat(pe_data[1][index, :].ravel(), pe_data[0][index, :].ravel())[stat][0])
		return array(res)

	indices = index_aridity()
	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	colors = [my_colors[0], my_colors[1], my_colors[3]]
	regions = ['humid (100 stations)', 'sub-humid (78 stations)', 'arid (50 stations)']

	plt.figure(figsize=(15, 4.5))
	for i, index in enumerate(indices):
		stats = [array_station_aridity(stat, index) for stat in statnames]
		plt.subplot(1,3,i+1)
		# sm.taylor_diagram(stats[0], stats[1], stats[2], markerLabel=modelnames, styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan', tickRMS = np.arange(0,6,1), tickCOR = intervalsCOR) #, tickSTD = np.arange(0,8.1,2), )
		if (i == 0) & (STEP != 'daily'):
			sm.taylor_diagram(stats[0], stats[1], stats[2],
                      styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan',
                      markerLabel=modelnames, markerLabelColor='k', markerColor=colors[i],
                      tickRMS=np.arange(0,5,1), tickRMSangle = 120.0,
                      styleRMS=':', widthRMS=2.0,
                      styleSTD='--', tickSTD=np.arange(0,6,1),
                      tickCOR=intervalsCOR, widthCOR=1.0, styleCOR='-.') # colCOR='m',
		else:
			sm.taylor_diagram(stats[0], stats[1], stats[2],
                      styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan',
                      markerLabel=modelnames, markerLabelColor='k', markerColor=colors[i],
                      tickRMS=np.arange(0,5,1), tickRMSangle = 120.0,
                      styleRMS=':', widthRMS=2.0,
                      styleSTD='--',
                      tickCOR=intervalsCOR, widthCOR=1.0, styleCOR='-.') # colCOR='m',
		plt.title(regions[i], y=0.90)
	plt.suptitle(STEP, fontsize=20, y=1.001)
	plt.tight_layout()
	# savefig('%s/taylor_obs_model_%s_good_station.tif' %(figdir, STEP), dpi=200)

	plt.show()

	return

# Model_obs_taylor_diagram()
# exit()

def Model_obs_taylor_diagram_multiscale():

	"Compare the observed pan with the modelled pan"
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	penpan = load('%s/pe_mod_%s_good_stations' %(workspace, 'penpan'))
	penman = load('%s/pe_mod_%s_good_stations' %(workspace, 'penman'))
	priest = load('%s/pe_mod_%s_good_stations' %(workspace, 'priestley-taylor'))
	hamon = load('%s/pe_mod_%s_good_stations' %(workspace, 'hamon'))
	turc = load('%s/pe_mod_%s_good_stations' %(workspace, 'turc'))
	# daily case:
	pe_daily = [pan_obs_gapfill, penpan, penman, priest, hamon, turc]
	# weekly case
	# pe_weekly = [daily2weekly_df(pe) for pe in pe_daily]
	# monthly case
	pe_monthly = [daily2monthly_df(pe) for pe in pe_daily]
	# annual case
	pe_annual = [daily2annual_df(pe) for pe in pe_daily]
	# print pe_annual[0].shape;exit()

	import skill_metrics as sm
	statnames = ['sdev', 'crmsd', 'ccoef']
	modelnames = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	# scales = ['daily', 'weekly', 'monthly', 'annual']
	# datascales = [pe_daily, pe_weekly, pe_monthly, pe_annual]
	# scales = ['daily', 'monthly', 'annual']
	# datascales = [pe_daily, pe_monthly, pe_annual]
	scales = ['Daily', 'Monthly']
	datascales = [pe_daily, pe_monthly]

	# intervalsCOR = np.concatenate((np.arange(0,1.0,0.2), [0.8, 0.9, 0.95, 0.99, 1]))

	# Calculate statistics for Taylor diagram
	#  The first array element corresponds to the reference series
	#  for the while the second is that for the predicted series.

	def Taylor_stat(model, obs):
		return sm.taylor_statistics(model, obs)

	def array_all_station(stat, pe_data):

		# res = [Taylor_stat(signal.detrend(pe, type='linear').ravel(), signal.detrend(pe_data[0], type='linear').ravel())[stat][1] for pe in pe_data[1:]] # "Add detrend component"
		res = [Taylor_stat(pe.ravel(), pe_data[0].ravel())[stat][1] for pe in pe_data[1:]]
		res.insert(0, Taylor_stat(pe_data[1].ravel(), pe_data[0].ravel())[stat][0])
		return array(res)

	# stats = [array_one_station(stat, 5) for stat in statnames]
	# stats = [array_all_station(stat) for stat in statnames]

	plt.figure(figsize=(11, 5))
	for i, data in enumerate(datascales):
		stats = [array_all_station(stat, data) for stat in statnames]
		plt.subplot(1,2,i+1)
		sm.taylor_diagram(stats[0], stats[1], stats[2],
					  markerLabel=modelnames, markerLabelColor='k', markerColor='k',
					  styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan',
					  tickRMS=np.arange(0,5,1), tickRMSangle = 120.0,
					  styleRMS=':', widthRMS=2.0,
					  styleSTD='--', widthSTD=1.5, tickSTD=np.arange(0,8,2),
					  widthCOR=1.5, styleCOR='-.')  #tickCOR=intervalsCOR,
		plt.title(scales[i], y=0.90, fontsize=20)
	plt.tight_layout()
	# savefig('%s/taylor_obs_model_all_multiscale_good_station.tif' %(figdir), dpi=200)
	# savefig('%s/taylor_obs_model_all_multiscale_detrend_good_station.tif' %(figdir), dpi=200)

	plt.show()

	return

# Model_obs_taylor_diagram_multiscale()
# exit()

def Model_obs_taylor_diagram_multiscale_new():

	"Compare the observed pan with the modelled pan"
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	penpan = load('%s/pe_mod_%s_good_stations' %(workspace, 'penpan'))
	penman = load('%s/pe_mod_%s_good_stations' %(workspace, 'penman'))
	priest = load('%s/pe_mod_%s_good_stations' %(workspace, 'priestley-taylor'))
	hamon = load('%s/pe_mod_%s_good_stations' %(workspace, 'hamon'))
	turc = load('%s/pe_mod_%s_good_stations' %(workspace, 'turc'))
	# daily case:
	pe_daily = [pan_obs_gapfill, penpan, penman, priest, hamon, turc]
	# weekly case
	pe_weekly = [daily2weekly_df(pe) for pe in pe_daily]
	# monthly case
	pe_monthly = [daily2monthly_df(pe) for pe in pe_daily]
	# annual case
	pe_annual = [daily2annual_df(pe) for pe in pe_daily]

	import skill_metrics as sm
	statnames = ['sdev', 'crmsd', 'ccoef']
	modelnames = ['Pan', 'Penpan', 'Penman', 'PT', 'Hamon', 'Turc']
	scales = ['Daily', 'Weekly', 'Monthly', 'Annual']
	datascales = [pe_daily, pe_weekly, pe_monthly, pe_annual]

	# Calculate statistics for Taylor diagram
	#  The first array element corresponds to the reference series
	#  for the while the second is that for the predicted series.

	def Taylor_stat(model, obs):
		return sm.taylor_statistics(model, obs)

	def array_all_station(stat, pe_data):
		res = [Taylor_stat(signal.detrend(pe, type='linear').ravel(), signal.detrend(pe_data[0], type='linear').ravel())[stat][1] for pe in pe_data[1:]] # "Add detrend component"
		res.insert(0, Taylor_stat(signal.detrend(pe_data[1], type='linear').ravel(), signal.detrend(pe_data[0], type='linear').ravel())[stat][0])

		# res = [Taylor_stat(pe.ravel(), pe_data[0].ravel())[stat][1] for pe in pe_data[1:]]
		# res.insert(0, Taylor_stat(pe_data[1].ravel(), pe_data[0].ravel())[stat][0])

		return array(res)

	plt.figure(figsize=(10, 9))
	for i, data in enumerate(datascales):
		stats = [array_all_station(stat, data) for stat in statnames]
		plt.subplot(2,2,i+1)
		# if i == 1:
		# 	sm.taylor_diagram(stats[0], stats[1], stats[2],
		# 			  markerLabel=modelnames, markerLabelColor='k', markerColor='k', markerLegend='on',
		# 			  styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan',
		# 			  styleRMS=':', widthRMS=2.0,
		# 			  styleSTD='--', widthSTD=1.5,
		# 			  widthCOR=1.5, styleCOR='-.',
		# 			  titleCOR='off')
		# else:
		sm.taylor_diagram(stats[0], stats[1], stats[2],
					  markerLabel=modelnames, markerLabelColor='k', markerColor='k',
					  styleOBS='-', colOBS='r', markerobs='o', titleOBS='Pan',
					  styleRMS=':', widthRMS=2.0,
					  styleSTD='--', widthSTD=1.5,
					  widthCOR=1.5, styleCOR='-.')
		# plt.title(scales[i], y=0.90, fontsize=18)
		if i < 3:
			plt.text(4.6, 6, scales[i], fontsize=18)
		else:
			plt.text(0.32, 0.4, scales[i], fontsize=18)

	plt.tight_layout()
	# savefig('%s/fig6_taylor_obs_model_four_multiscale_detrend_good_station.tif' %(figdir), dpi=300)
	plt.show()

	return

# Model_obs_taylor_diagram_multiscale_new()
# exit()



