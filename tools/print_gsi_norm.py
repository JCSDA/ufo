from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma

d = Dataset('sondes_obs_2018041500_m.nc4', 'r')
t = d['air_temperature@GsiHofx'][:]
u = d['eastward_wind@GsiHofx'][:]
v = d['northward_wind@GsiHofx'][:]

t_ = ma.compressed(t)
u_ = ma.compressed(u)
v_ = ma.compressed(v)

print "GSI norm (T only): ",  np.sqrt(np.mean(t_**2))
print "GSI norm (T, u): ",    np.sqrt(np.mean(np.concatenate([t_, u_])**2))
print "GSI norm (T, u, v): ", np.sqrt(np.mean(np.concatenate([t_, u_, v_])**2))
