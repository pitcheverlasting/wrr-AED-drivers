__author__ = 'lpeng'
import netCDF4 as netcdf
import numpy as np
import datetime, sys, os
import scipy.io
## start GrADS
##------------
# sys.path.insert(0, "/home/water5/lpeng/Programming/python/mylib/")
# from grads import *
# from gradsen import *
# gradsen()
# ga = GrADS(Bin='/home/latent2/mpan/local/opengrads/Linux/x86_64/grads', Verb=False, Echo=False, Port=False, Window=False, Opts="-c 'q config'")

def Create_NETCDF_File(dims, file, vars, varname, data, tinitial, tstep, nt):

	nlat = dims['nlat']
	nlon = dims['nlon']
	res = dims['res']
	minlon = dims['minlon']
	minlat = dims['minlat']
	undef = dims['undef']
	t = np.arange(0, nt)  #dims['nyr'])

	# Prepare the netcdf file
	# Create file
	f = netcdf.Dataset(file, 'w', format='NETCDF4')

	# Define dimensions
	f.createDimension('lon', nlon)
	f.createDimension('lat', nlat)
	f.createDimension('t', len(t))

	# Longitude
	f.createVariable('lon', 'd', ('lon',))
	f.variables['lon'][:] = np.linspace(minlon, minlon+res*(nlon-1), nlon)
	f.variables['lon'].units = 'degrees_east'
	f.variables['lon'].long_name = 'Longitude'
	f.variables['lon'].res = res

	# Latitude
	f.createVariable('lat', 'd', ('lat',))
	f.variables['lat'][:] = np.linspace(minlat, minlat+res*(nlat-1), nlat)
	f.variables['lat'].units = 'degrees_north'
	f.variables['lat'].long_name = 'Latitude'
	f.variables['lat'].res = res

	# Time
	f.createVariable('t', 'd', ('t', ))
	f.variables['t'][:] = t
	f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (tstep,tinitial.year,tinitial.month,tinitial.day,tinitial.hour)
	f.variables['t'].long_name = 'Time'

	# Data
	if type(vars) is str:
		datafield = f.createVariable(vars, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
		f.variables[vars].long_name = varname
		datafield[:] = data
	else:
		for v, var in enumerate(vars):
			datafield = f.createVariable(var, 'f', ('t', 'lat', 'lon',), fill_value=undef, zlib=True)
			f.variables[var].long_name = varname[v]
			datafield[:] = data[v, :, :]

	f.sync()
	f.close()

	return f


def Grads2netcdf(ctlfile, outdir, filename):

	os.system('~/anaconda/bin/cdo -r -f nc import_binary %s %s/%s.nc' % (ctlfile, outdir, filename))

	return

def mat2python():

	station_geoinfo = scipy.io.loadmat('station_geoinfo.mat')
	geoinfo = np.empty((485, 5))
	geoinfo[:,0] = station_geoinfo['stationsbybasin'][:,3] # station ID
	geoinfo[:,1] = station_geoinfo['stationsbybasin'][:,1] # station lat
	geoinfo[:,2] = station_geoinfo['stationsbybasin'][:,0] # station lon
	geoinfo[:,3] = station_geoinfo['stationsbybasin'][:,2] # station elv
	geoinfo[:,4] = station_geoinfo['stationsbybasin'][:,4] # station basin ID

	return

def daily2monthly(dates, daily_time_series):
	dmonths = dates.month
	res = np.diff(dmonths)
	print res
# 	for