import geospatialtools.gdal_tools as gdal_tools
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import propr_tools
import time
np.random.seed(1)
#Read in the properties (top layer)
dtype={'names': ('code','ud','ld','min','mid','max'),
       'formats': ('S100','f4','f4','f4','f4','f4')}
data = np.loadtxt('properties/20130606_dlr_reference_and_uncertainty_propr.txt',
     skiprows=1,dtype=dtype,usecols=(0,1,2,5,6,7),delimiter=',')
mask = data['ud'] == 0.0
tmp = {}
for var in data.dtype.names:
 tmp[var] = data[var][mask]
data = tmp
data['id'] = np.arange(len(data['code']))
nc = len(data['code'])
nd = 1000
metadata = gdal_tools.retrieve_metadata('probabilities/%s.tif' % data['code'][0])
probabilities = np.zeros((nc,metadata['ny'],metadata['nx']))
for ic in xrange(nc):
 file = 'probabilities/%s.tif' % data['code'][ic]
 probabilities[ic,:,:] = gdal_tools.read_raster(file)[:]
#Sort the array
argsort = np.flipud(np.argsort(probabilities,axis=0)).astype(np.int32)
probs = np.flipud(np.sort(probabilities,axis=0)).astype(np.float32)
#Draw n draws per soil class
draws = np.zeros((nc,nd)).astype(np.float32)
for ic in data['id']:
 #c = 0.1
 loc = data['min'][ic]
 scale = data['max'][ic]-data['min'][ic]
 mode = (data['mid'][ic] - data['min'][ic])/(data['max'][ic] - data['min'][ic])
 tmp = stats.triang.rvs(mode,loc=loc,scale=scale,size=nd)
 draws[ic,:] = tmp
#Create the draws array
'''output = np.zeros((nd,metadata['ny'],metadata['nx']))
for i in xrange(output.shape[1]):
 print i
 for j in xrange(output.shape[2]):
  samples = np.ceil(nd*probs[:,i,j])
  mask = samples > 0
  if np.sum(samples) > nd:
   dif = np.sum(samples) - nd
   for k in xrange(samples.size):
    samples[k] = samples[k] - 1
    dif = dif - 1
    if dif == 0:break
  args = argsort[:,i,j]
  tmp = 0
  for k in xrange(np.sum(mask)):
   output[tmp:tmp+samples[k],i,j] = draws[args[k],0:samples[k]]
   tmp = tmp + samples[k]'''
tic = time.time()
output = propr_tools.assign_draws(probs,argsort,draws)
toc = time.time()
print toc - tic
output = np.sort(output,axis=0)

#Construct the distributions
plt.imshow(np.mean(output,axis=0))
plt.colorbar()
plt.show()
