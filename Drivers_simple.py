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

styr = 1962
edyr = 2001
stdy = datetime.datetime(styr, 1, 1)
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
###====================================================================

def Get_station_aridity():

	for ibasin in xrange(0, 10):
		aridity_basin = []
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
			prec = Gapfill(data['rain'][0, istation][0:tstep].flatten())
			aridity_basin.append(nanmean(pan_obs_gapfill)/nanmean(prec))
		aridity_basin = array(aridity_basin)
		aridity_basin.dump('%s/aridity_station_%s' %(workspace, basinlongs[ibasin]))

	return

# Get_station_aridity()

def Get_station_mean():

		# 	pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
		# 	prec = Gapfill(data['rain'][0, istation][0:tstep].flatten())
		# 	aridity_basin.append(nanmean(pan_obs_gapfill)/nanmean(prec))
		# aridity_basin = array(aridity_basin)
	var_avg = []
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))
	ist = 0

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			input.insert(1, rn[ist, :].flatten())
			input.insert(0, pan_obs_gapfill[ist, :].flatten())
			var_avg.append(nanmean(vstack(input), axis=1))
			ist = ist + 1

	avg = vstack(var_avg)
	avg.dump('%s/6var_avg_good_station' %(workspace))

	return

# Get_station_mean()
# exit()

def Get_basin_aridity():

	Aridity = []
	pan_all = []
	prec_all = []
	for ibasin in xrange(0, 10):
		pan_basin = []
		prec_basin = []
		for istation in good_stations[ibasin]:
			print ibasin, istation
			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
			pan_obs_gapfill = Gapfill(data['pan'][0, istation][0:tstep].flatten())
			prec = Gapfill(data['rain'][0, istation][0:tstep].flatten())
			pan_basin.append(nanmean(pan_obs_gapfill))
			prec_basin.append(nanmean(prec))
		pan_all.append(pan_basin)
		prec_all.append(prec_basin)

		Aridity.append(nanmean(vstack(pan_basin))/nanmean(vstack(prec_basin)))
	Aridity.append(nanmean(vstack(flatten(pan_all)))/nanmean(vstack(flatten(prec_all))))
	return Aridity
# arid = vstack(Get_basin_aridity())
# arid.dump('%s/aridity_basin' %workspace)

def msc_groupby_DI(arid):
	# aridrange = [0, 1, 1.5, 2.5, 5, 250]
	aridrange = [0, 2, 4, 8, 250]

	index_DI = []
	for igroup in xrange(0, len(aridrange)-1):
		low = aridrange[igroup]
		high = aridrange[igroup+1]
		index_DI.append(np.where((arid>=low) & (arid<high))[0])
	return index_DI


def Plot_Station_Data():

	lons = load('%s/lons_good_stations' %workspace)
	lats = load('%s/lats_good_stations' %workspace)

	arid = []
	for ibasin in xrange(0, 10):
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)
	# print [len(index_DI[i]) for i in xrange(0, 4)]
	# print mean(arid), median(arid)
	# print [median(arid[index_DI[i]]) for i in xrange(0, 4)]
	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	arid_class = ['humid', 'sub-humid', 'semi-arid', 'arid']

	[Plotting.Mapshow_Aridity(lons[index_DI[i]], lats[index_DI[i]], 50, 0.8, 0, 1000, my_colors[i], arid_class[i]) for i in xrange(0, 4)]
	# plt.show()
	fig = plt.gcf()
	fig.set_size_inches(6, 5)
	plt.show()
	# savefig('%s/map_aridity_good_station.tif' %figdir, dpi=300)

	return

# Plot_Station_Data()
# exit()

def Plot_Station_Mean():

	lons = load('%s/lons_good_stations' %workspace)
	lats = load('%s/lats_good_stations' %workspace)

	avg = load('%s/6var_avg_good_station' %(workspace))
	vars = ['Pan', 'Tair','Rn','Wind','Humidity','VPD']
	units = ['(mm)', r'($^oC$)', '(MJ/m2/d)', '(m/s)', '(hpa)', '(hpa)']

	for i in xrange(0, 6):
		Plotting.Mapshow(avg[:, i], lons, lats, 30, 0.7, None, None, 'jet', vars[i], None, units[i], None, None)
		savefig('%s/map_%s_avg_good_station.tif' %(figdir, vars[i]), dpi=200)
		plt.clf()
	return

# Plot_Station_Mean()
# exit()

##########################################################################
# Coherece analysis
##########################################################################
nf = 513
# nf = 1025
nbasin = 10
nvar = 8
sampling_frequency = 1/(24.0 * 3600.0) # unit: per day


def Coherence_Frequency():

	data = scipy.io.loadmat('%s/1_AP.mat' %(datadir))
	input = data[variables[0]][0, 0][0:tstep].flatten()
	pan = data['pan'][0, 0][0:tstep].flatten()
	freq = FFT.Coherence(input, pan, sampling_frequency, 'linear')[0]

	return freq

# print len(Coherence_Frequency()), sampling_frequency/Coherence_Frequency();exit()


def PSD_5var():

	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))
	ist = 0
	for ibasin in xrange(0, 10):
		psd_basin = []
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			input.insert(1, rn[ist, :].flatten())
			psd_basin.append(vstack([FFT.Power_Spectrum(v, sampling_frequency, 'linear')[1] for v in input]).reshape(1,5, nf))
			ist = ist + 1

		psd_basin = vstack(psd_basin)
		psd_basin.dump('%s/psd_5var_good_station_%s' %(workspace, basinlongs[ibasin]))

	return

def Crosspectrum_5var():
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))
	ist = 0
	for ibasin in xrange(0, 10):
		cross_basin = []
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			input.insert(1, rn[ist, :].flatten())
			panobs = pan_obs_gapfill[ist, :].flatten()
			cross_basin.append(vstack([FFT.CrossPowerSpectrum(panobs, v, sampling_frequency, 'linear')[1] for v in input]).reshape(1, 5, nf))
			ist = ist + 1

		cross_basin = vstack(cross_basin)
		cross_basin.dump('%s/cross_5var_good_station_%s' %(workspace, basinlongs[ibasin]))

	return
# PSD_5var()
# Crosspectrum_5var()
# exit()

def Plot_Cross_spectrum_Average():

	fig = plt.figure(figsize=(11, 5))
	cohere = []
	freq = Coherence_Frequency()
	arid = load('%s/aridity_basin' %workspace)
	for ibasin in xrange(0, 10):
		# Plot coherence
		# cohere_basin = load('%s/coherence_obs_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		# Plot cospectra
		cohere_basin = load('%s/cross_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		cohere.append(cohere_basin)

	# for all stations average
	ax = fig.add_subplot(1, 1, 1)
	col = ['Crimson', 'Orange', 'DeepSkyblue', 'SpringGreen', 'Purple']
	drivername = ['Tair','Rn','Wind','Humidity','VPD']
	[ax.semilogx(sampling_frequency/freq, 10*np.log10(mean(vstack(cohere), axis=0)[i, :]), color=col[i], linewidth=2.5, alpha=0.6, label=drivername[i]) for i in (4, 1, 0, 3, 2)]  # xrange(nvar-1, -1, -1)
	## add yaxis labels and change font size
	ax.set_xlabel(r"$\tau$ (day)", fontsize=18)
	ax.xaxis.set_label_coords(1.05, -0.01)
	ax.set_ylabel("Cross-spectrum", fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	plt.grid(True)
	plt.show()
	return
# Plot_Cross_spectrum_Average()
# exit()


def ann_cycle_obs_5var_aridity():
	"seasonal cycle climatology"
	keys = ['pan', 'rn', 'vpd', 'ea', 'tavg', 'wind'] # 'rh',
	keynames = ['Epan', r'$R_n$', 'VPD',  r'$e_a$', r'$T_a$', r'$u_2$']  # 'RH',
	cols = ['black', 'Orange', 'Purple', 'springgreen', 'Crimson', 'DeepSkyblue'] # 'olive',
	numbers = ['a', 'b', 'c', 'e', 'f']
	arid_class = ['average', r'humid (0<$\phi$<2)', r'sub-humid (2<$\phi$<4)', r'semi-arid (4<$\phi$<8)', r'arid ($\phi$>8)']

	def annual_all_records(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.month)]
		return array(ann)
	def zscore(ts, mean, std):
		# return (ts - np.nanmean(ts))/np.nanstd(ts)
		return (ts - mean)/std

	ist = 0
	met_station = []
	arid = []
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))

		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			dataplot = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in keys[2:6]]
			dataplot.insert(0, pan_obs_gapfill[ist, :].flatten())
			dataplot.insert(1, rn[ist, :].flatten())
			met_station.append(array([annual_all_records(pd.Series(dataplot[i], index=dates)) for i in xrange(6)]))
			ist = ist + 1

	met_station = array(met_station)
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)
	met_avg = [mean(met_station[index, :, :], axis=0) for index in index_DI]
	met_avg.insert(0, mean(met_station, axis=0))
	met_mean = mean(met_avg[0], axis=1)
	met_std = std(met_avg[0], axis=1)

	fig = plt.figure(figsize=(12, 7.5))

	lons = load('%s/lons_good_stations' %workspace)
	lats = load('%s/lats_good_stations' %workspace)

	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	arid_labels = ['humid', 'sub-humid', 'semi-arid', 'arid']
	ax = fig.add_subplot(2, 3, 4)
	[Plotting.Mapshow_Aridity_subplot(ax, lons[index_DI[i]], lats[index_DI[i]], 30, 0.8, 0, 1000, my_colors[i], arid_labels[i]) for i in xrange(0, 4)]

	plot_loc = (1, 2, 3, 5, 6)
	for i, mets in enumerate(met_avg):
		ax = fig.add_subplot(2, 3, plot_loc[i])
		lines = [ax.plot(zscore(met, met_mean[im], met_std[im]), '--', dashes=(12, 5), color=cols[im], label=keynames[im], linewidth=3, alpha=0.9, rasterized=True) for im, met in enumerate(mets)]  #  '-o', ms=10, mew=0,

		ax.set_xlim(0, 11)
		ax.set_xticks(range(12))
		ax.set_xticklabels(range(1,13))
		ax.set_ylim(-3, 4)
		ymin, ymax = ax.get_ylim()
		ax.text(10.8, ymin+0.03*(ymax-ymin), arid_class[i], ha='right', fontsize=16)
		ax.text(0.8, ymax-0.13*(ymax-ymin), '(%s)' %numbers[i], fontsize=16)
		ax.tick_params(axis='y', labelsize=14)
		ax.tick_params(axis='x', labelsize=14)
		ax.grid('on')

	lines = list(itertools.chain(*lines))
	plt.figlegend(lines, keynames, loc='lower center',  ncol=6, fontsize=16)   # doesn't work bbox_to_anchor=(0, -0.1, 1, 1),, bbox_transform=plt.gcf().transFigure
	# plt.subplots_adjust(bottom=0.2)
	plt.tight_layout(rect=[0, 0.07, 1, 1])
	# savefig('%s/fig1_ann_cycle_5var_aridity_good_station.tif' %figdir, dpi=300)
	plt.show()

	return

# ann_cycle_obs_5var_aridity()
# exit()


def ann_cycle_obs_5var_aridity_update():
	"seasonal cycle climatology"
	keys = ['pan', 'rn', 'vpd', 'ea', 'tavg', 'wind'] # 'rh',
	keynames = ['Epan', r'$R_n$', 'VPD',  r'$e_a$', r'$T_a$', r'$u_2$']  # 'RH',
	cols = ['black', 'Orange', 'Purple', 'springgreen', 'Crimson', 'DeepSkyblue'] # 'olive',
	numbers = ['d', 'b', 'c', 'e', 'f']
	arid_class = ['average', r'humid (0<$\phi$<2)', r'sub-humid (2<$\phi$<4)', r'semi-arid (4<$\phi$<8)', r'arid ($\phi$>8)']

	def annual_all_records(df):
		"This is a new function to calculate mean of each hours across all sites, not independently"
		ann = [np.nanmean(g) for _, g, in df.groupby(df.index.month)]
		return array(ann)
	def zscore(ts):
		return (ts - np.nanmean(ts))/np.nanstd(ts)

	ist = 0
	met_station = []
	arid = []
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))

	for ibasin in xrange(0, 10):
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))

		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			dataplot = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in keys[2:6]]
			dataplot.insert(0, pan_obs_gapfill[ist, :].flatten())
			dataplot.insert(1, rn[ist, :].flatten())
			# met_station.append(array([zscore(annual_all_records(pd.Series(dataplot[i], index=dates))) for i in xrange(6)]))
			met_station.append(array([annual_all_records((pd.Series(zscore(dataplot[i]), index=dates))) for i in xrange(6)]))
			ist = ist + 1

	met_station = array(met_station)
	arid = array(list(itertools.chain(*arid)))
	index_DI = msc_groupby_DI(arid)
	met_avg = [mean(met_station[index, :, :], axis=0) for index in index_DI]
	met_avg.insert(0, mean(met_station, axis=0))

	fig = plt.figure(figsize=(12, 7.5))

	lons = load('%s/lons_good_stations' %workspace)
	lats = load('%s/lats_good_stations' %workspace)

	colors = r_[linspace(0, 1, 4), linspace(0, 1, 4)]
	mymap = plt.get_cmap("coolwarm")
	my_colors = mymap(colors)
	arid_labels = ['humid', 'sub-humid', 'semi-arid', 'arid']
	ax = fig.add_subplot(2, 3, 1)
	[Plotting.Mapshow_Aridity_subplot(ax, lons[index_DI[i]], lats[index_DI[i]], 30, 0.8, 0, 1000, my_colors[i], arid_labels[i]) for i in xrange(0, 4)]

	plot_loc = (4, 2, 3, 5, 6)
	for i, mets in enumerate(met_avg):
		ax = fig.add_subplot(2, 3, plot_loc[i])
		lines = [ax.plot(met, '--', dashes=(12, 5), color=cols[im], label=keynames[im], linewidth=3, alpha=0.9, rasterized=True) for im, met in enumerate(mets)]  #  '-o', ms=10, mew=0,

		ax.set_xlim(0, 11)
		ax.set_xticks(range(12))
		ax.set_xticklabels(range(1,13))
		ax.set_ylim(-2, 2)
		ymin, ymax = ax.get_ylim()
		ax.text(10.8, ymin+0.03*(ymax-ymin), arid_class[i], ha='right', fontsize=16)
		ax.text(0.8, ymax-0.13*(ymax-ymin), '(%s)' %numbers[i], fontsize=16)
		ax.tick_params(axis='y', labelsize=14)
		ax.tick_params(axis='x', labelsize=14)
		ax.grid('on')

	lines = list(itertools.chain(*lines))
	plt.figlegend(lines, keynames, loc='lower center',  ncol=6, fontsize=16)   # doesn't work bbox_to_anchor=(0, -0.1, 1, 1),, bbox_transform=plt.gcf().transFigure
	# plt.subplots_adjust(bottom=0.2)
	plt.tight_layout(rect=[0, 0.07, 1, 1])
	# savefig('%s/fig1_ann_cycle_5var_aridity_good_station.tif' %figdir, dpi=300)
	plt.show()

	return

ann_cycle_obs_5var_aridity_update()
exit()


def Coherence_obs_5var():

	"Prepare the data for plotting"
	pan_obs_gapfill = load('%s/pan_obs_gapfill_good_stations' %workspace)
	rn = load('%s/Rn_0.2_0.5_good_stations' %(workspace))
	ist = 0
	for ibasin in xrange(0, 10):
		cohere_obs_basin = []
		data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
		for istation in good_stations[ibasin]:
			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
			input = [Gapfill(data[v][0, istation][0:tstep].flatten()).flatten() for v in variables]
			input.insert(1, rn[ist, :].flatten())
			panobs = pan_obs_gapfill[ist, :].flatten()
			# Compute the coherence
			cohere_obs_basin.append(vstack([FFT.Coherence(v, panobs, sampling_frequency, 'linear')[1] for v in input]).reshape(1, 5, nf))
			ist = ist + 1

		# store basin average
		cohere_obs_basin = vstack(cohere_obs_basin)
		cohere_obs_basin.dump('%s/coherence_obs_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		# cohere_obs_basin.dump('%s/coherence_obs_5var_test2048_%s' %(workspace, basinlongs[ibasin]))

	return
# Coherence_obs_5var()
# exit()

def Plot_Coherence_Average():

	fig = plt.figure(figsize=(10, 4))
	cohere = []
	freq = Coherence_Frequency()
	for ibasin in xrange(0, 10):
		# Plot coherence
		cohere_basin = load('%s/coherence_obs_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		# cohere_basin = load('%s/coherence_obs_5var_test2048_%s' %(workspace, basinlongs[ibasin]))

		cohere.append(cohere_basin)

	# for all stations average
	ax = fig.add_subplot(1, 1, 1)
	Plotting.CoherenceVarPlot(ax, mean(vstack(cohere), axis=0), sampling_frequency, freq)
	plt.show()
	# savefig('%s/fig2_coh_obs_5var_all_good_station.tif' %figdir, dpi=300)
	return
# Plot_Coherence_Average()
# exit()


############################################
# combine all the old figures into new figures
############################################
def Epan_VPD_Rn_Partial_Coherence_Groupby_Aridity():
	"Figure 3, put the Epan, Rn-VPD and partial together into the same figure"

	pcoh = []
	arid = []
	Epan_cohere = []
	VPD_cohere = []

	for ibasin in xrange(0, 10):
		epan_cohere = load('%s/coherence_obs_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		vpd_cohere = load('%s/coherence_vpd_other_vars_good_station_%s' %(workspace, basinlongs[ibasin]))
		pcoh_basin = load('%s/pcoh_pan_vpd_rn_good_station_%s' %(workspace, basinlongs[ibasin]))
		Epan_cohere.append(epan_cohere)
		VPD_cohere.append(vpd_cohere)
		pcoh.append(pcoh_basin)
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	Epan_cohere = vstack(Epan_cohere)
	VPD_cohere = vstack(VPD_cohere)
	pcoh = vstack(pcoh)

	index_DI = msc_groupby_DI(arid)


	dataplot = [Epan_cohere[:, 4, :], Epan_cohere[:, 1, :], pcoh[:, 0, :], pcoh[:, 1, :], VPD_cohere[:, 0, :]]
	freq = Coherence_Frequency()

	fig = plt.figure(figsize=(10, 8.5))
	drivername = [r'$C_{Epan,VPD}$', r'$C_{Epan,R_n}$', r'$C_{Epan,VPD/R_n}$', r'$C_{Epan,R_n/VPD}$', r'$C_{VPD,R_n}$']
	labels = ('a', 'b', 'c', 'd', 'e')
	positions = (1, 2, 3, 4, 5)


	for i, var in enumerate(drivername):
		ax = fig.add_subplot(3, 2, positions[i])
		Plotting.CoherenceAridityPlot(ax, dataplot[i], index_DI, sampling_frequency, freq, labels[i], var)
		if i == 4:
			ax.legend(bbox_to_anchor=(1.45, 0.5), loc='center left', fontsize=16)


	fig.tight_layout()
	plt.show()
	# savefig('%s/fig3_epan_vpd_rn_coh_pcoh_aridity_good_station.tif' %figdir, dpi=300)

	return

# Epan_VPD_Rn_Partial_Coherence_Groupby_Aridity()
# exit()


def Coherence_noVPDRn_Groupby_Aridity():
	"Figure 4, put the VPD, Rn, Tair, ea together "
	arid = []
	cohere_1 = []; 	cohere_2 = []; cohere_3 = []

	for ibasin in xrange(0, 10):
		cohere_basin_1 = load('%s/coherence_vpd_other_vars_good_station_%s' %(workspace, basinlongs[ibasin]))
		cohere_basin_2 = load('%s/coherence_rn_other_vars_good_station_%s' %(workspace, basinlongs[ibasin]))
		cohere_basin_3 = load('%s/coherence_tavg_other_vars_good_station_%s' %(workspace, basinlongs[ibasin]))

		cohere_1.append(cohere_basin_1)
		cohere_2.append(cohere_basin_2)
		cohere_3.append(cohere_basin_3)

		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	cohere_1 = vstack(cohere_1)
	cohere_2 = vstack(cohere_2)
	cohere_3 = vstack(cohere_3)


	index_DI = msc_groupby_DI(arid)

	fig = plt.figure(figsize=(10, 8.5))
	drivername = [r'$C_{VPD,Tair}$', r'$C_{VPD,e_a}$', r'$C_{R_n,Tair}$', r'$C_{R_n,e_a}$', r'$C_{Tair,e_a}$']
	labels = ('a', 'b', 'c', 'd', 'e')
	freq = Coherence_Frequency()

	dataplot = [cohere_1[:, 1, :], cohere_1[:, 2, :], cohere_2[:, 0, :], cohere_2[:, 1, :], cohere_3[:, 0, :]]

	for i, var in enumerate(drivername):
		ax = fig.add_subplot(3, 2, i+1)

		Plotting.CoherenceAridityPlot(ax, dataplot[i], index_DI, sampling_frequency, freq, labels[i], var)
		if i == 4:
			ax.legend(bbox_to_anchor=(1.45, 0.5), loc='center left', fontsize=16)

	fig.tight_layout()
	plt.show()
	# savefig('%s/fig4_coh_inter_var_aridity_good_station.tif' %figdir, dpi=300)

	return

# Coherence_noVPDRn_Groupby_Aridity()
# exit()


def Wind_Coherence_Groupby_Aridity():

	cohere = []
	arid = []
	for ibasin in xrange(0, 10):
		cohere_basin = load('%s/coherence_obs_5var_good_station_%s' %(workspace, basinlongs[ibasin]))
		cohere.append(cohere_basin)
		arid.append(load('%s/aridity_station_%s' %(workspace, basinlongs[ibasin])))
	arid = array(list(itertools.chain(*arid)))
	cohere = vstack(cohere)

	index_DI = msc_groupby_DI(arid)
	freq = Coherence_Frequency()

	fig = plt.figure(figsize=(6, 3.5))
	ax = fig.add_subplot(1, 1, 1)
	Plotting.CoherenceAridityPlot(ax, cohere[:, 2, :], index_DI, sampling_frequency, freq, '', r'$C_{Epan,u_2}$')

	fig.tight_layout()
	plt.show()
	# savefig('%s/fig5_coh_Epan_u2_aridity_good_station.tif' %figdir, dpi=300)

	return

# Wind_Coherence_Groupby_Aridity()
# exit()

