import netCDF4 as nc
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) != 2):
  print("One argument needed: filename of the fil with obs error correlation diagnostics")
  sys.exit()
filename = sys.argv[1]
ds = nc.Dataset(filename)
coords = ds.variables['correlationCoordinate'][:]
corrs = ds.variables['correlations'][:]
rand = ds.variables['randomVector'][:]
randMult = ds.variables['randomVectorMultipliedByCov'][:]
# plot correlation matrix
plt.imshow(corrs,cmap='hot')
plt.title('correlation matrix for selected record')
plt.colorbar()
plt.savefig("corr_record.png")
plt.clf()
# plot several different correlations as function of coordinate
for icorr in range(0, len(coords), 10):
  plt.plot(coords, corrs[icorr])
plt.title('correlations between different locations in one record')
plt.savefig("corr_location.png")
plt.clf()
# plot random vector, and random vector multiplied by correlation
plt.plot(coords, rand, label='Random vector')
plt.plot(coords, randMult, label='Random vector multiplied by covariance')
plt.legend()
plt.savefig("corr_random.png")
