import propr_tools_fortran
import json
import netCDF4 as nc
import numpy as np
import scipy.stats as stats
import time
import warnings
warnings.filterwarnings("ignore")
undef = -9999.0

def initialize_netcdf_output(instance,nx,ny):

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
   grp.createVariable('mean','f4',('z','y','x'),chunksizes=(1,ny,nx))
   grp.createVariable('max','f4',('z','y','x'),chunksizes=(1,ny,nx))
   grp.createVariable('min','f4',('z','y','x'),chunksizes=(1,ny,nx))
   grp.createVariable('alpha','f4',('z','y','x'),chunksizes=(1,ny,nx))
   grp.createVariable('beta','f4',('z','y','x'),chunksizes=(1,ny,nx))
  else:
   grp.createVariable('mean','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('max','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('min','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('alpha','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('beta','f4',('y','x'),chunksizes=(ny,nx))
  #Initialize to undef
  stats = ['mean','min','max','alpha','beta']
  #for stat in stats:
  # grp[stat][:] = undef

 return fp

def calculate_properties_on_each_block(instance,fp,ix,iy,iblock):

 #Define the arrays
 mr = instance.metadata['maxrank']
 instance.prob = instance.fp_probabilities.variables[instance.metadata['vm']['prob']][0:mr,iy,ix]
 instance.prob = instance.metadata["vm"]["prob_scalar_multiplier"]*instance.prob
 instance.rank = instance.fp_probabilities.variables[instance.metadata['vm']['rank']][0:mr,iy,ix]
 instance.block_ids = np.ma.getdata(np.unique(instance.rank))
 #instance.array = np.zeros((instance.metadata['nd'],instance.prob.shape[1],instance.prob.shape[2]),dtype='f4', order='F')
 
 #If there are no unique ids then return
 bids = instance.block_ids
 bids = bids[bids >= 0]
 if len(np.unique(bids)) == 0:
  for var in instance.metadata['vars']:
   if len(instance.variables[var]['shape']) == 2:
    for il in xrange(instance.nl):
     for stat in ['mean','min','max','alpha','beta']:
       fp[var][stat][il,iy,ix] = undef
   else:
     for stat in ['mean','min','max','alpha','beta']:
       fp[var][stat][iy,ix] = undef
  return

 #Curate the probabilities (always needs to sum to 100)
 instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 
 #Create the input data for a variable for each layer
 vars = instance.metadata['vars']
 for var in vars:

  t0 = time.time()
  #print instance.metadata['process_id'],"Calculating %s maps" % var
  #Retrieve the variables
  grp = fp[var]
 
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

  t1 = time.time()
  print instance.metadata['process_id'],"Finished %s maps in %f seconds" % (var,(t1-t0))
  
 return

def run(file):

 #Initialize the propr class
 instance = initialize(file)

 #Get the block info
 n = instance.metadata['nblocks']+1
 xbounds = np.ceil(np.linspace(instance.ix[0],instance.ix[-1],n)).astype(np.int32)
 ybounds = np.ceil(np.linspace(instance.iy[0],instance.iy[-1],n)).astype(np.int32)

 #Initialize the output file
 print instance.metadata['process_id'],"Initializing output file"
 nx = xbounds[1] - xbounds[0]
 ny = ybounds[1] - ybounds[0]
 fp = initialize_netcdf_output(instance,nx,ny)

 iblock = 0
 for i in xrange(len(xbounds)-1):
  ix = np.arange(xbounds[i],xbounds[i+1]+1)
  for j in xrange(len(ybounds)-1):
   iy = np.arange(xbounds[j],ybounds[j+1]+1)
   iblock = iblock + 1
   print instance.metadata['process_id'],"Working on block %d/%d" % (iblock,(n-1)*(n-1))
   #Calculate the properties for this block
   calculate_properties_on_each_block(instance,fp,ix,iy,iblock)

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
  #mapping = np.zeros(np.max(data['id'])+1).astype(np.int32)
  mapping = np.zeros(np.max(ids)+1).astype(np.int32)
  mapping[ids] = np.arange(nc)
  #Compute the index in the data
  idx = np.in1d(mids,ids)
  #Determine the classes that have information
  m = idx & (data['mean'] != undef)
  #Draw from beta distribution
  #tmp = np.random.beta(data['alpha'][m],data['beta'][m],(nd,np.sum(m))) THIS CREATES BOUNDARIES
  tmp = []
  for i in xrange(np.sum(m)):
   #Define the random seed
   rnd = np.random.RandomState(self.metadata['seed'])
   tmp.append(rnd.beta(data['alpha'][m][i],data['beta'][m][i],(nd,)))
  tmp = np.array(tmp).T
  #Rescale using min and max
  tmp = tmp*(data['max'][m] - data['min'][m]) + data['min'][m]
  #Place the samples in the draws array
  m = np.in1d(ids,mids[m])
  draws[m,:] = tmp.T
  
  return (draws,mapping)

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Create the array of draws
  t0 = time.time()
  (draws,mapping) = self.draw_from_distribution(data)
  #print 'sampling',time.time() - t0

  #Create the array of draws
  t0 = time.time()
  ncmax = self.metadata['ncmax']
  #array = self.array
  array = propr_tools_fortran.assign_draws(self.prob,self.rank,draws,mapping,ncmax)
  #print 'placing data',time.time() - t0

  #Compute the beta parameters
  t0 = time.time()
  min = np.min(array,axis=0)
  max = np.max(array,axis=0)
  mean = np.mean(array,axis=0)
  #narray = (array - min[np.newaxis,:,:])/(max[np.newaxis,:,:]-min[np.newaxis,:,:])
  narray = (array - min)/(max-min)
  nmean = np.mean(narray,axis=(0,))
  nvar = np.var(narray,axis=(0,))
  alpha = ((1-nmean)/nvar - (1/nmean))*nmean**2
  beta = alpha*(1/nmean - 1)
  #print 'creating beta parameters',time.time() - t0
  
  #Assemble output
  output = {}
  output['mean'] = mean
  output['min'] = min
  output['max'] = max
  output['alpha'] = alpha
  output['beta'] = beta

  return output

