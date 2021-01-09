__author__ = 'lpeng'

from pylab import *
from scipy import signal
import scipy
import math

def Power_Spectrum(ts, fs, dt):

	"""Compute power spectral density (using psd welch)"""
	"ts: time series of measurement values"
	"fs: sampling frequency"
	"dt: detrend = None/'linear'/'constant'"
	"scal = 'density', normalized with units of V**2/hz"
	"scal = 'spectrum', units of V**2"

	ts[np.isnan(ts)==True] = nanmean(ts)		# Gapfill the missing value
	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length

	############################### nfft: Number of DFT points ###############################
	# The default nfft is the greater of 256 or the next power of 2 greater than the length of the segments.
	# The relation between segment length and the Number of DFT points: segment_length <= nfft
	# If nfft is greater than the segment length, the data is zero-padded. If nfft is less than the segment length, the segment is wrapped using datawrap to make the length equal to nfft
	##########################################################################################

	###############################   welch's method for PSD   ###############################
	# MATLAB function: pwelch(ts, segment_length, Number of overlapped window noverlap, Number of DFT points, fs,'onesided')
	# scipy method: pwelch(ts, fs, nperseg, noverlap, nfft, detrend, scaling)
	##########################################################################################

	f, Pxx_den = signal.welch(ts, fs, nperseg=window, detrend=dt)

	return [f, Pxx_den]


def CrossPowerSpectrum(ts1, ts2, fs, dt):

	"""Compute Cross spectral density (using csd)"""

	# Gapfill the missing value
	mk1 = np.isnan(ts1)
	mk2 = np.isnan(ts2)

	ts1[mk1==True] = nanmean(ts1)
	ts1[mk2==True] = nanmean(ts1)
	ts2[mk1==True] = nanmean(ts2)
	ts2[mk2==True] = nanmean(ts2)

	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts1)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length

	f, Pxy_den = signal.csd(ts1, ts2, fs, nperseg=window, detrend=dt)

	return [f, Pxy_den]

def Coherence_calculate(ts1, ts2, fs, dt):
	"""Compute magnitude squared coherence (calculate)"""

	f, Pxy_den = CrossPowerSpectrum(ts1, ts2, fs, dt)
	Pxx_den = Power_Spectrum(ts1, fs, dt)[1]
	Pyy_den = Power_Spectrum(ts2, fs, dt)[1]
	Cxy = abs(Pxy_den)**2/(Pxx_den*Pyy_den)

	return [f, Cxy, Pxy_den, Pxx_den, Pyy_den]

def Coherence(ts1, ts2, fs, dt):

	"""Compute magnitude squared coherence (using coherence)"""

	# Gapfill the missing value
	mk1 = np.isnan(ts1)
	mk2 = np.isnan(ts2)
	ts1[mk1==True] = nanmean(ts1)
	ts1[mk2==True] = nanmean(ts1)
	ts2[mk1==True] = nanmean(ts2)
	ts2[mk2==True] = nanmean(ts2)

	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts1)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	window = int(2 ** nwindow_fl)   			# segment_length
	# test
	# window = 2048

	f, Cxy = signal.coherence(ts1, ts2, fs, nperseg=window, detrend=dt)

	return [f, Cxy]

def Partial_Coherence(ts1, ts2, ts3, fs, dt):
	"x, and y is the two signals of interest, z is the variable that you wants to get rid of"
	f, Pxy_den = CrossPowerSpectrum(ts1, ts2, fs, dt)
	f, Pxz_den = CrossPowerSpectrum(ts1, ts3, fs, dt)
	f, Pyz_den = CrossPowerSpectrum(ts2, ts3, fs, dt)
	f, Pxx_den = Power_Spectrum(ts1, fs, dt)
	f, Pyy_den = Power_Spectrum(ts2, fs, dt)
	f, Pzz_den = Power_Spectrum(ts3, fs, dt)
	# f, Cxy = Coherence(ts1, ts2, fs, dt)
	# f, Cxz = Coherence(ts1, ts3, fs, dt)
	# f, Cyz = Coherence(ts2, ts3, fs, dt)
	# PGxy_z = abs(Pxy_den - Pxz_den*Pyz_den/Pzz_den)**2 / (Pxx_den - abs(Pxz_den)**2/Pzz_den) / (Pyy_den - abs(Pyz_den)**2/Pzz_den)
	PGxy_z = (abs(Pxy_den) - abs(Pxz_den)*abs(Pyz_den)/abs(Pzz_den))**2 / ((abs(Pxx_den) - abs(Pxz_den)**2/abs(Pzz_den)) * (abs(Pyy_den) - abs(Pyz_den)**2/abs(Pzz_den)))

	return [f, PGxy_z]

def Confidence_interval(ts1, ts2, fs, dt):
	"Wang and Tang, 2004"
	nperseg = 1024.0
	f, Cxy = Coherence(ts1, ts2, fs, dt)
	error = (1 - Cxy) * sqrt(4*Cxy/nperseg)

	return [(1 - 2*error)*Cxy, (1 + 2*error)*Cxy]

def Confidence_Limits(ts, fs, alpha):

	"ts: time series of measurement values"
	"fs: sampling frequency"
	"dt: detrend = None/'linear'/'constant'"
	"scal = 'density', normalized with units of V**2/hz"
	"scal = 'spectrum', units of V**2"

	ts[np.isnan(ts)==True] = nanmean(ts)		# Gapfill the missing value
	nwindow = 9          						# Number of windows to smooth data
	length = math.floor(len(ts)/nwindow)  		# Length calculated by deviding the window
	nwindow_fl = math.floor(log2(length))		# Number of windows with floor window length
	nperseg = int(2 ** nwindow_fl)   		    # window
	# res_bandwidth = fs/nperseg					# Reslution bandwidth
	# ndisjoint = res_bandwidth * len(ts)	* 360 * 2400		# Number of disjoint segments used in spectra estimates
	ndisjoint = len(ts) / float(nperseg)
	E = 1 - (1 - alpha)	** (1/(ndisjoint-1))	# confidence interval when alpha = 0.95 confidence level

	return E


##======================================================
# no use
## method2: average the time series first and then do the spectrum
# def Coherence_average2Basin():
#
# 	fig = plt.figure(figsize=(12, 30))
# 	input = []; pan = []; cohere = []
#
# 	for ibasin in xrange(0, 10):
# 		input_basin = zeros((station_number[ibasin], nvar, tstep))
# 		pan_basin = zeros((station_number[ibasin], tstep))
# 		for istation in xrange(0, station_number[ibasin]):
# 			data = scipy.io.loadmat('%s/%s_AP.mat' %(datadir, ibasin+1))
# 			## Quality check
# 			input_mask = [np.isnan(nanmean(data[v][0, istation])) for v in variables]
# 			flag = []
# 			flag = [1 for i in input_mask if i == True] # which means one variable is missing
# 			if len(flag) > 0:
# 				print "Ignore %s station %s!" % (basinlongs[ibasin], istation)
# 				continue
#
# 			# the PowerSpectrum method take the matrix as different segment, so shoule be a 1d array
# 			input_basin[istation, :, :] = vstack([data[v][0, istation].flatten() for v in variables])
# 			pan_basin[istation, :] = data['pan'][0, istation].flatten()
#
# 		input.append(mean(input_basin, axis=0).reshape(1, nvar, tstep))
# 		pan.append(mean(pan_basin, axis=0))
# 		del input_basin, pan_basin
# 		# Compute the coherence
# 		freq = FFT.Coherence(input[ibasin][0, 0, :], pan[ibasin], sampling_frequency, 'linear')[0]
# 		coh = vstack([FFT.Coherence(input[ibasin][0, v, :], pan[ibasin], sampling_frequency, 'linear')[1] for v in xrange(0, nvar)])  # nvar, nf
# 		# store basin average
# 		cohere.append(coh.reshape(1, nvar, nf))
#
# 		## DRAW FIGURE---------------
# 	 	ax = fig.add_subplot(4, 3, ibasin+1)
# 		Plotting.CoherenceBasinPlot(ax, coh, sampling_frequency, freq, basinlongs[ibasin])
#
# 	# for national average
# 	# method 1: calculate national average first then coherence
# 	input_all = mean(vstack(input), axis=0) # nvar, tstep
# 	pan_all = mean(vstack(pan), axis=0).flatten() # 1, tstep
# 	coh1 = vstack([FFT.Coherence(input_all[v, :].flatten(), pan_all, sampling_frequency, 'linear')[1] for v in xrange(0, nvar)])
# 	ax = fig.add_subplot(4, 3, 11)
# 	Plotting.CoherenceBasinPlot(ax, coh1, sampling_frequency, freq, 'Average1')
# 	# method 2: average the basin coherence to get the national one
# 	coh2 = mean(vstack(cohere), axis=0)
# 	ax = fig.add_subplot(4, 3, 12)
# 	Plotting.CoherenceBasinPlot(ax, coh2, sampling_frequency, freq, 'Average2')
# 	# ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left')
# 	# fig.tight_layout()
# 	plt.show()
#
# 	return