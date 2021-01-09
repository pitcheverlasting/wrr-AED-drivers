__author__ = 'lpeng'
"""Penpan library
refering Mcmahon, 2013; Rotstayn, 2006; Li, 2013; Yang, 2012; Roderick, 2007"""

from pylab import *

class Data:

	# Initialization
	def __init__(self, INPUT, solar): #, npt):

		# Define the incoming grid data variables
		self.Tair = INPUT['tavg']
		self.Tmin = INPUT['tmax']
		self.Tmax = INPUT['tmin']
		self.Pres = self.Convert_Unit_Pres(INPUT['p'])
		self.e = self.Convert_Unit_Pres(INPUT['ea'])    # vapor pressure: hpa
		self.Wind = INPUT['wind']  # Wind: m/s
		self.sunhour = INPUT['sun']  # sunhour: hour
		# self.CF = self.Convert_Unit_CF(INPUT['tc'])
		self.lat = INPUT['lat']
		self.elev = INPUT['elev']
		self.doy = INPUT['doy']
		self.Albedo = 0.23

		# Calculate some variables
		self.Calculate_Tmean()
		# self.Calculate_Saturated_Vapor_Pressure()
		self.Calculate_Mean_Saturated_Vapor_Pressure()
		self.Calculate_Slope_Saturation_Vapor_Pressure_Curve()
		self.Calculate_Gamma()
		self.Calculate_VPD()
		self.Calculate_Rso()
		self.Calculate_Rs(solar)
		self.Calculate_Rs_pan()
		self.Calculate_LWnet()
		self.Calculate_Rnet_pan() # npt temporary with npt
		self.Penpan()
		# self.Penman()
		# self.Priestley_Taylor()
		# self.Hamon()
		# self.Turc()

	def Calculate_Tmean(self):

		self.Tmean = (self.Tmax + self.Tmin) / 2.0

		return


	def Calculate_Saturated_Vapor_Pressure(self):

		self.estar = 0.6108 * np.exp((17.27 * self.Tair) / (237.3 + self.Tair))

		return

	def Calculate_Mean_Saturated_Vapor_Pressure(self):

		estar_max = 0.6108 * np.exp((17.27 * self.Tmax) / (237.3 + self.Tmax))
		estar_min = 0.6108 * np.exp((17.27 * self.Tmin) / (237.3 + self.Tmin))
		self.estar = (estar_max + estar_min) / 2

		return

	def Calculate_Slope_Saturation_Vapor_Pressure_Curve(self):

		# self.DELTA = 4098 * 0.6108 * np.exp((17.27 * self.Tmean) / (237.3 + self.Tmean)) / (237.3 + self.Tmean) ** 2
		self.DELTA = 4098 * self.estar / (237.3 + self.Tmean) ** 2

		return

	def Calculate_Rs(self, solar):
		"Compute the SW irradiance from the input data"
		"R_s : downward solar irradiance at the surface"

		if solar == "sunhours":
			"Sunshine hour data have been used for calculating incoming solar radiation"
			## using Angstrom-Prescott equation
			# a_s = 0.25
			# b_s = 0.5
			# yellow river
			a_s = 0.195
			b_s = 0.5125
			self.SWin = (a_s + b_s * (self.sunhour/self.daylength)) * self.Ra


		elif solar == "cloud":  # this is from Linacre, 1993, equation 19
			"Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
			self.SWin = (0.85 - 0.047 * self.CF) * self.Ra

		return

	def Calculate_Rso(self):

		# dr is inverse relative distance Earth-Sun
		dr = 1 + 0.033 * np.cos(2 * np.pi * self.doy/365.0)
		# delta is solar declination
		delta = 0.409 * np.sin((2*np.pi * self.doy/365.0) - 1.39)
		# from decimal degree to radians
		phi = (np.pi/180) * self.lat
		# omega: sunset hour angle
		omega = np.arccos(-np.tan(phi) * np.tan(delta))
		# daylength: the maximum daylight hours
		self.daylength = 24 / np.pi * omega
		if self.elev:
			z2 = self.elev
		else:
			z2 = 30  # set the station elevation above the sea level as 30m
		Gsc = 0.082   # solar constant = 0.082 [MJ m^-2 min^-1]
		# Ra: extraterrestrial radiation for daily period for different location, different daytime
		self.Ra = 24 * 60 / np.pi * Gsc * dr * (omega * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(omega))
		self.Rso = (0.75 + 2e-5 * z2) * self.Ra

		return


	def Calculate_Rs_pan(self):
		"Compute the total SW irradiance of the pan"
		"Prad : pan evaporation factor, which accounts for the extra direct irradiance intercepted by the walls of the pan when the sun is not directly overhead"
		"f_dir : fraction of fs that is direct"

		# Linacre, 1994
		# P_rad = 1.32 + 4 * 10**(-4) * abs(self.lat) + 8 * 10 ** (-5) * self.lat ** 2
		P_rad = 1.70 + 3 * 10 ** (-4) * self.lat ** 2
		f_dir = -0.11 + 1.31 * self.SWin / self.Ra  		#  extraterrestrial radiation
		self.Rsp = (f_dir * P_rad + 2.0 * (1 - f_dir) + 2.0 * self.Albedo) * self.SWin

		return

	def Calculate_LWnet(self):

		stefan_b = 4.903e-9  # [MJ K-4 m-2 day-1]
		epsilon_s = 0.98
		self.LWnet = stefan_b * ((self.Tmax+273.16)**4 + (self.Tmin+273.16)**4) / 2.0 * (0.34-0.14 * np.sqrt(self.e)) * (1.35 * self.SWin/self.Rso - 0.35)
		# self.LWnet = stefan_b * (self.Tair+273.16)**4 * (0.34-0.14 * np.sqrt(self.e)) * (1.35 * self.SWin/self.Rso - 0.35)
		# self.PDLWnet_PDTair = - 4 * stefan_b * (self.Tair+273.16)**3 * (0.34-0.14 * np.sqrt(self.e)) * (1.35 * self.SWin/self.Rso - 0.35)

		return

	def Calculate_Rnet_pan(self):  #, npt):

		Ap = 0.14  # Class A pan albedo (Linacre, 1992; Rotstayn, 2006; Roderick, 2007; Yang and Yang, 2011)
		self.Rn_pan = (1 - Ap) * self.Rsp - self.LWnet

		# # temporary
		# import pandas as pd
		# def running_mean(y, npts):
		# 	return pd.rolling_mean(y, npts, center=True)
		# Rn_pan = (1 - Ap) * self.Rsp - self.LWnet
		# self.Rn_pan = running_mean(Rn_pan, npt)
		# # temporary

		return

	def Calculate_Gamma(self):

		cp = 1.013   # Specific heat of moist air at constant pressure [kJ kg-1 C-1]
		self.lv = 2.501 - 0.002361 * self.Tair  # Latent heat vapor (MJ/kg)
		self.gamma = ((cp * self.Pres) / (0.622 * self.lv)) * 10 ** (-3)
		return self.gamma

	## Convert unit

	# W/m2 to MJ/m2/day
	def Convert_Unit_Rad(self, input):
		watt2jule = 10e5/86400.0
		data = input / float(watt2jule)
		return data

	# pressure: hpa to kpa
	def Convert_Unit_Pres(self, input):
		data = input / 10.0
		return data

	# 10m wind to 2m wind
	def Convert_Unit_Wind(self, input):

		zh = 10  # 10 m wind field
		data = (4.87 / (np.log(67.8 * zh - 5.42))) * input

		return data

	def Convert_Unit_CF(self, input):
		data = input / 10.0
		return data

	## Calculation for each components

	def Calculate_VPD(self):

		self.vpd = (self.estar - self.e) * (self.estar>=self.e) + 0 * (self.estar<self.e)

		return self.vpd

	# Reference-surface Models
	def Penpan(self):

		"ET_type = D20 Pan Evaporation"

		# These parameters are from Yang 2012
		coeff_hq = 5  # for D20 pan
		# f_pan_u = 2.6 * (1 + 0.536 * self.Wind)
		# Beijing experiment for vapor transfer function
		# f_pan_u = 5.4 * (1 + 0.73 * self.Wind)/self.lv
		# f_pan_u = 1.201 + 1.621 * self.Wind
		f_pan_u = 1.313 + 1.381 * self.Wind
		# f_pan_u = 2.626 + 1.381 * self.Wind
		# f_pan_u = 0.35*(1+9.8e-3 * self.Wind)
		# f_pan_u = 1.39e-8*(1+1.35 * self.Wind)

		PET_R = self.DELTA/(self.DELTA + coeff_hq * self.gamma) * self.Rn_pan / self.lv
		PET_A = coeff_hq * self.gamma/(self.DELTA + coeff_hq * self.gamma) * f_pan_u * self.vpd
		self.penpan = PET_R + PET_A

		return

	def Penpan_2parts(self):

		"ET_type = D20 Pan Evaporation"

		# These parameters are from Yang 2012
		coeff_hq = 5  # for D20 pan
		f_pan_u = 1.313 + 1.381 * self.Wind

		PET_R = self.DELTA/(self.DELTA + coeff_hq * self.gamma) * self.Rn_pan / self.lv
		PET_A = coeff_hq * self.gamma/(self.DELTA + coeff_hq * self.gamma) * f_pan_u * self.vpd

		return PET_R, PET_A

	def Penpan_u2(self, a, b):

		"ET_type = D20 Pan Evaporation"

		# These parameters are from Yang 2012
		coeff_hq = 5  # for D20 pan
		# f_pan_u = 2.6 * (1 + 0.536 * self.Wind)
		# Beijing experiment for vapor transfer function
		# f_pan_u = 1.313 + 1.381 * self.Wind
		f_pan_u = a + b * self.Wind

		PET_R = self.DELTA/(self.DELTA + coeff_hq * self.gamma) * self.Rn_pan / self.lv
		PET_A = coeff_hq * self.gamma/(self.DELTA + coeff_hq * self.gamma) * f_pan_u * self.vpd / self.lv
		penpan = PET_R + PET_A

		return penpan

	def Penman(self):

		# Hydrology Book
		PET_R = (self.DELTA / (self.DELTA + self.gamma)) * self.Rn_pan / self.lv
		PET_A = (self.gamma / (self.DELTA + self.gamma)) * ((6.43 * (1 + 0.536 * self.Wind) * self.vpd) / self.lv)
		self.penman = PET_R + PET_A

		return

	def Priestley_Taylor(self):

		alpha = 1.26
		self.PT = alpha * self.DELTA / (self.DELTA + self.gamma) * self.Rn_pan / self.lv

		return

	def Hamon(self):

		# Tair is daily mean air temperature
		mask = np.zeros(len(self.Tair))
		mask[self.Tair > 0] = 1
		calibration = 2
		estar = 6.108 * np.exp(17.27 * self.Tair / (self.Tair + 237.3)) # you must use daily average temperature
		self.hamon = 0.1651 * 216.7 * estar / (self.Tair + 273.3) * (self.daylength / 12.0) * calibration * mask

		return

	def Turc(self):

		mask = np.zeros(len(self.Tair))
		mask[self.Tair > 0] = 1
		RH = self.e/self.estar*100
		turc = (0.013 * (23.88 * self.Rsp + 50) * (self.Tair) / (self.Tair + 15) * (RH>=50) + 0.013 * (23.88 * self.Rsp + 50) * (self.Tair) / (self.Tair + 15) * (1 + (50-RH)/70) * (RH<50)) * mask

		# this mask doesn't work for the Tair=-15 case, because this will make the denum zero, so the result is nan at first place, it is still nan after multiplying a number
		turc[np.isnan(turc)==True] = 0
		self.turc = turc

		return
