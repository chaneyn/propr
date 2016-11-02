import propr_tools_fortran
import json
import netCDF4 as nc
import numpy as np
import scipy.stats as stats
import time
import warnings
warnings.filterwarnings("ignore")

def initialize_netcdf_output(instance):

 ud = instance.fp_properties.variables['upper_depth'][:]
 ld = instance.fp_properties.variables['lower_depth'][:]
 fp = nc.Dataset(instance.metadata['output_file'],'w')
 fp.createDimension('z',instance.nl)
 fp.createDimension('x',instance.nx)
 fp.createDimension('y',instance.ny)
 x = fp.createVariable('x','f8',('x',))
 x.Axis = 'X'
 x[:] = instance.x
 y = fp.createVariable('y','f8',('y',))
 y.Axis = 'Y'
 y[:] = instance.y
 z = fp.createVariable('z','f8',('z',))
 z.Axis = 'Z'
 z[:] = (ud[0:instance.nl] + ld[0:instance.nl])/2
 z = fp.createVariable('upper_depth','f8',('z',))
 z[:] = ud[0:instance.nl]
 z = fp.createVariable('lower_depth','f8',('z',))
 z[:] = ld[0:instance.nl]

 #Create the variables for the maps
 for var in instance.metadata['vars']:
  grp = fp.createGroup(var)
  grp.units = instance.variables[var]['units']
  grp.description = instance.variables[var]['description']
  #Create new variables
  shape = instance.variables[var]['shape']
  if len(shape) == 2:
   grp.createVariable('mean','f4',('z','y','x'))
   grp.createVariable('max','f4',('z','y','x'))
   grp.createVariable('min','f4',('z','y','x'))
   grp.createVariable('alpha','f4',('z','y','x'))
   grp.createVariable('beta','f4',('z','y','x'))
  else:
   grp.createVariable('mean','f4',('y','x'))
   grp.createVariable('max','f4',('y','x'))
   grp.createVariable('min','f4',('y','x'))
   grp.createVariable('alpha','f4',('y','x'))
   grp.createVariable('beta','f4',('y','x'))

 return fp

def calculate_properties_on_each_block(instance,fp,ix,iy):

 #Define the arrays
 mr = instance.metadata['maxrank']
 instance.prob = instance.fp_probabilities.variables[instance.metadata['vm']['prob']][0:mr,iy,ix]
 instance.prob = instance.metadata["vm"]["prob_scalar_multiplier"]*instance.prob
 instance.rank = instance.fp_probabilities.variables[instance.metadata['vm']['rank']][0:mr,iy,ix]
 instance.block_ids = np.ma.getdata(np.unique(instance.rank))

 #Curate the probabilities (always needs to sum to 100)
 instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 
 #Create the input data for a variable for each layer
 #vars = instance.fp_properties.groups.keys()
 vars = instance.metadata['vars']
 for var in vars:

  t0 = time.time()
  print "Calculating %s maps" % var
  #Retrieve the variables
  grp = fp[var]
  #vpc5 = fp.variables['%s_pc5' % var]
  #vpc95 = fp.variables['%s_pc95' % var]
  #vmean = fp.variables['%s_mean' % var]
 
  if len(instance.variables[var]['shape']) == 2:

   #Iterate through each layer
   for il in xrange(instance.nl):

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['min'] = instance.fp_properties.groups[var].variables['min'][:,il]
    data['max'] = instance.fp_properties.groups[var].variables['max'][:,il]
    data['mean'] = instance.fp_properties.groups[var].variables['mean'][:,il]
    data['alpha'] = instance.fp_properties.groups[var].variables['alpha'][:,il]
    data['beta'] = instance.fp_properties.groups[var].variables['beta'][:,il]

    #Calculate the properties
    output = instance.calculate_properties(data)

    #Output the properties
    for stat in output:
     fp[var][stat][il,iy,ix] = output[stat]
    #vpc5[il,iy,ix] = output['pc5']
    #vpc95[il,iy,ix] = output['pc95']
    #vmean[il,iy,ix] = output['mean']

  else:

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['min'] = instance.fp_properties.groups[var].variables['min'][:]
    data['max'] = instance.fp_properties.groups[var].variables['max'][:]
    data['mean'] = instance.fp_properties.groups[var].variables['mean'][:]
    data['alpha'] = instance.fp_properties.groups[var].variables['alpha'][:]
    data['beta'] = instance.fp_properties.groups[var].variables['beta'][:]

    #Calculate the properties
    output = instance.calculate_properties(data)

    #Output the properties
    for stat in output:
     fp[var][stat][iy,ix] = output[stat]
    #vpc5[iy,ix] = output['pc5']
    #vpc95[iy,ix] = output['pc95']
    #vmean[iy,ix] = output['mean']
  
  print time.time() - t0

 return

def run(file):

 #Initialize the propr class
 instance = initialize(file)

 #Initialize the output file
 print "Initializing output file"
 fp = initialize_netcdf_output(instance)

 #Iterate per block
 n = instance.metadata['nblocks']+1
 xbounds = np.ceil(np.linspace(instance.ix[0],instance.ix[-1],n)).astype(np.int32)
 ybounds = np.ceil(np.linspace(instance.iy[0],instance.iy[-1],n)).astype(np.int32)

 iblock = 0
 for i in xrange(len(xbounds)-1):
  ix = np.arange(xbounds[i],xbounds[i+1]+1)
  for j in xrange(len(ybounds)-1):
   iy = np.arange(xbounds[j],ybounds[j+1]+1)
   iblock = iblock + 1
   print "Working on block %d/%d" % (iblock,(n-1)*(n-1))
   #Calculate the properties for this block
   calculate_properties_on_each_block(instance,fp,ix,iy)

 #Finalize netcdf file
 print "Finalizing output file"
 fp.close()
 
 return

class initialize:


 #Define functions
 def __init__(self,file):

  #Read in the metadata
  self.metadata = json.load(open(file))

  #Open access to input data
  self.fp_properties = nc.Dataset(self.metadata['properties_file'])
  self.fp_probabilities = nc.Dataset(self.metadata['input_file'])

  #Read some metadata
  self.x = self.fp_probabilities.variables[self.metadata['vm']['x']][:]
  self.y = self.fp_probabilities.variables[self.metadata['vm']['y']][:]
  self.nx = len(self.x)
  self.ny = len(self.y)
  self.ix = np.arange(self.nx)
  self.iy = np.arange(self.ny)
  self.nl = len(self.fp_properties.variables['upper_depth'][:])
  if self.nl > self.metadata['maxnl']:self.nl = self.metadata['maxnl']

  #Construct the shapes of all the input variables
  self.variables = {}
  for var in self.metadata['vars']:
   self.variables[var] = {'shape':self.fp_properties[var]['mean'].shape,
                          'units':self.fp_properties[var].units,
                          'description':self.fp_properties[var].description}

  return

 #Draw data per soil class
 def draw_from_distribution(self,data):
 
  #Define parameters
  undef = -9999
  ids = self.block_ids
  ids = ids[ids >= 0]
  nc = len(ids)#len(data['id'])
  nd = self.metadata['nd']
  mids = data['id']
  #Construct an array of undefined draws
  draws = np.zeros((nc,nd)).astype(np.float32)
  draws[:] = undef
  #Construct mapping array
  mapping = np.zeros(np.max(data['id'])+1).astype(np.int32)
  mapping[ids] = np.arange(nc)
  #Compute the index in the data
  idx = np.in1d(mids,ids)
  #Determine the classes that have information
  m = idx & (data['mean'] != undef)
  #Draw from beta distribution
  tmp = np.random.beta(data['alpha'][m],data['beta'][m],(nd,np.sum(m)))
  #Rescale using min and max
  tmp = tmp*(data['max'][m] - data['min'][m]) + data['min'][m]
  #Place the samples in the draws array
  m = np.in1d(ids,mids[m])
  draws[m,:] = tmp.T
  
  #m = data['mean']
  '''ic = 0
  print time.time() - tic
  for id in ids:
   #Find index of id
   if id in mids:
    idx = mids.index(id)
    if data['mean'][idx] == -9999:
     draws[ic,:] = -9999
     mapping[id] = ic
    else:
     alpha = data['alpha'][idx]
     beta = data['beta'][idx]
     min = data['min'][idx]
     max = data['max'][idx]
     #draw from beta distribution
     tmp = np.random.beta(alpha,beta,(nd,))
     #rescale using min and max
     tmp = tmp*(max - min) + min
     draws[ic,:] = tmp
     mapping[id] = ic
   else:
    draws[ic,:] = -9999
    mapping[id] = ic
   ic += 1
  print time.time() - tic'''

  return (draws,mapping)

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Create the array of draws
  t0 = time.time()
  (draws,mapping) = self.draw_from_distribution(data)
  print 'sampling',time.time() - t0

  #Create the array of draws
  t0 = time.time()
  ncmax = self.metadata['ncmax']
  array = propr_tools_fortran.assign_draws(self.prob,self.rank,draws,mapping,ncmax)
  print 'placing data',time.time() - t0

  #Compute the beta parameters
  t0 = time.time()
  min = np.min(array,axis=0)
  max = np.max(array,axis=0)
  mean = np.mean(array,axis=0)
  nmean = np.mean((array - min[np.newaxis,:,:])/(max[np.newaxis,:,:]-min[np.newaxis,:,:]),axis=(0,))
  nvar = np.var((array - min[np.newaxis,:,:])/(max[np.newaxis,:,:]-min[np.newaxis,:,:]),axis=(0,))
  alpha = ((1-nmean)/nvar - (1/nmean))*nmean**2
  beta = alpha*(1/nmean - 1)
  print 'creating beta parameters',time.time() - t0
  
  #Sort the data
  #array = np.sort(array,axis=0)

  #Compute the 5th pct, 95th pct, and mean
  #ipc5 = np.percentile(np.arange(self.metadata['nd']),5,interpolation='nearest')
  #ipc95 = np.percentile(np.arange(self.metadata['nd']),95,interpolation='nearest')
  #output = {}
  #output['mean'] = np.mean(array,axis=0)
  #output['pc5'] = array[ipc5,:,:]
  #output['pc95'] = array[ipc95,:,:]

  #Assemble output
  output = {}
  output['mean'] = mean
  output['min'] = min
  output['max'] = max
  output['alpha'] = alpha
  output['beta'] = beta

  return output

