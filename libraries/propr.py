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
 fp.createDimension('b',instance.nb)
 fp.createDimension('eb',instance.nb+1)
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
  if len(shape) == 3:
   grp.createVariable('hist','f4',('y','x','z','b'),chunksizes=(ny,nx,1,1),zlib=True,complevel=1,shuffle=True,least_significant_digit=3)
   grp.createVariable('ebins','f4',('eb'))
  else:
   grp.createVariable('mean','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('max','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('min','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('alpha','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('beta','f4',('y','x'),chunksizes=(ny,nx))
   grp.createVariable('var','f4',('y','x'),chunksizes=(ny,nx))
  #Initialize to undef
  #stats = ['mean','min','max','alpha','beta','var']
  #for stat in stats:
  # grp[stat][:] = undef

 return fp

def initialize_netcdf_output_point(instance):

 ud = instance.fp_properties.variables['upper_depth'][:]
 ld = instance.fp_properties.variables['lower_depth'][:]
 fp = nc.Dataset(instance.metadata['output_file'],'w')
 fp.createDimension('z',instance.nl)
 fp.createDimension('p',instance.np)
 fp.createDimension('b',instance.nb)
 fp.createDimension('eb',instance.nb+1)
 lat = fp.createVariable('lat','f4',('p',))
 lat[:] = instance.fp_probabilities['lat'][:]
 lon = fp.createVariable('lon','f4',('p',))
 lon[:] = instance.fp_probabilities['lon'][:]
 p = fp.createVariable('p','i4',('p',))
 p[:] = instance.fp_probabilities['site'][:]
 z = fp.createVariable('z','f8',('z',))
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
  if len(shape) == 3:
   grp.createVariable('hist','f4',('p','z','b'))
   grp.createVariable('ebins','f4',('eb'))
  else:
   grp.createVariable('mean','f4',('p',))
   grp.createVariable('max','f4',('p',))
   grp.createVariable('min','f4',('p',))
   grp.createVariable('alpha','f4',('p',))
   grp.createVariable('beta','f4',('p',))
   grp.createVariable('var','f4',('p',))
  #Initialize to undef
  #stats = ['mean','min','max','alpha','beta','var']
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
 
 #If there are no unique ids then return
 bids = instance.block_ids
 bids = bids[bids >= 0]
 if len(np.unique(bids)) == 0:
  for var in instance.metadata['vars']:
   if len(instance.variables[var]['shape']) == 3:
    for il in xrange(instance.nl):
     fp[var]['hist'][iy,ix,il,:] = undef
   else:
     fp[var]['hist'][iy,ix,:] = undef
  return

 #Curate the probabilities (always needs to sum to 1)
 instance.prob[instance.prob == -9999] = 0
 instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 
 #Create the input data for a variable for each layer
 vars = instance.metadata['vars']
 for var in vars:

  t0 = time.time()
  #print instance.metadata['process_id'],"Calculating %s maps" % var
  #Retrieve the variables
  grp = fp[var]
  print var
 
  if len(instance.variables[var]['shape']) == 3:

   #Iterate through each layer
   for il in xrange(instance.nl):

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['hist'] = instance.fp_properties.groups[var].variables['hist'][:,il,:]
    data['bins'] = instance.fp_properties.groups[var].variables['bins'][:]
    '''data['min'] = instance.fp_properties.groups[var].variables['min'][:,il]
    data['max'] = instance.fp_properties.groups[var].variables['max'][:,il]
    data['mean'] = instance.fp_properties.groups[var].variables['mean'][:,il]
    data['alpha'] = instance.fp_properties.groups[var].variables['alpha'][:,il]
    data['beta'] = instance.fp_properties.groups[var].variables['beta'][:,il]
    data['var'] = instance.fp_properties.groups[var].variables['var'][:,il]'''

    #Calculate the properties
    output = instance.calculate_properties(data)

    #Reduce the hist
    m = (np.diff(data['bins'])[np.newaxis,np.newaxis,:]*output['hist']) < 0.01
    output['hist'][m] = 0.0

    #Normalize
    m = np.sum(output['hist'],axis=2) > 0
    output['hist'][m,:] = output['hist'][m,:]/np.sum(output['hist'],axis=2)[m,np.newaxis]

    #output['hist'] = output['hist']/np.diff(data['bins'])[np.newaxis,np.newaxis,:]
    #Reduce to three significant digits
    tmp = np.round(1000*(output['hist']))/1000
    
    #Set zeros to -9999
    m = np.sum(tmp,axis=2) == 0
    tmp[m,:] = -9999

    #Output the properties
    fp[var]['hist'][iy,ix,il,:] = tmp[:]#output['hist']
    fp[var]['ebins'][:] = data['bins']

  else:

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['min'] = instance.fp_properties.groups[var].variables['min'][:]
    data['max'] = instance.fp_properties.groups[var].variables['max'][:]
    data['mean'] = instance.fp_properties.groups[var].variables['mean'][:]
    data['alpha'] = instance.fp_properties.groups[var].variables['alpha'][:]
    data['beta'] = instance.fp_properties.groups[var].variables['beta'][:]
    data['var'] = instance.fp_properties.groups[var].variables['var'][:]

    #Calculate the properties
    output = instance.calculate_properties(data)

    #Output the properties
    for stat in output:
     fp[var][stat][iy,ix] = output[stat]

  t1 = time.time()
  print instance.metadata['process_id'],"Finished %s maps in %f seconds" % (var,(t1-t0))
  
 return

def calculate_properties_point(instance,fp):

 #Define the arrays
 mr = instance.metadata['maxrank']
 instance.prob = instance.fp_probabilities.variables[instance.metadata['vm']['prob']][0:mr,:]
 instance.prob = instance.metadata["vm"]["prob_scalar_multiplier"]*instance.prob
 instance.rank = instance.fp_probabilities.variables[instance.metadata['vm']['rank']][0:mr,:]
 instance.block_ids = np.ma.getdata(np.unique(instance.rank))
 
 #Curate the probabilities (always needs to sum to 100)
 instance.prob[instance.prob == -9999] = 0
 instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 #instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 
 #Create the input data for a variable for each layer
 vars = instance.metadata['vars']
 for var in vars:

  t0 = time.time()
  #Retrieve the variables
  grp = fp[var]
  if len(instance.variables[var]['shape']) == 3:

   #Iterate through each layer
   for il in xrange(instance.nl):

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['hist'] = instance.fp_properties.groups[var].variables['hist'][:,il,:]
    data['bins'] = instance.fp_properties.groups[var].variables['bins'][:]
 
    #Calculate the properties
    output = instance.calculate_properties(data)

    #Set zeros to -9999
    m = np.sum(output['hist'],axis=1) == 0
    output['hist'][m,:] = -9999

    #Output the properties
    fp[var]['hist'][:,il] = output['hist']
    fp[var]['ebins'][:] = data['bins']

  else:

    #Define the soil property data
    data = {}
    data['id'] = instance.fp_properties.variables['id'][:]
    data['min'] = instance.fp_properties.groups[var].variables['min'][:]
    data['max'] = instance.fp_properties.groups[var].variables['max'][:]
    data['mean'] = instance.fp_properties.groups[var].variables['mean'][:]
    data['alpha'] = instance.fp_properties.groups[var].variables['alpha'][:]
    data['beta'] = instance.fp_properties.groups[var].variables['beta'][:]
    data['var'] = instance.fp_properties.groups[var].variables['var'][:]

    #Calculate the properties
    output = instance.calculate_properties(data)

    #Output the properties
    for stat in output:
     fp[var][stat][:] = output[stat]

  t1 = time.time()
  print instance.metadata['process_id'],"Finished %s points in %f seconds" % (var,(t1-t0))
  
 return


def run(file):

 #Read the metadata
 metadata = json.load(open(file))
 
 if metadata['type'] == 'point':
 
  #Initialize the propr class
  instance = initialize_point(metadata)

  #Initialize the output file
  print "Initializing output file"
  fp = initialize_netcdf_output_point(instance)

  #Calculating the properties for the points
  print "Calculating properties"
  calculate_properties_point(instance,fp)

 elif metadata['type'] == 'grid':

  #Initialize the propr class
  instance = initialize(metadata)

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

class initialize_point:

 #Define functions
 def __init__(self,metadata):

  #Read in the metadata
  self.metadata = metadata

  #Open access to input data
  self.fp_properties = nc.Dataset(self.metadata['properties_file'])
  self.fp_probabilities = nc.Dataset(self.metadata['input_file'])

  #Read some metadata
  #self.np = len(self.fp_probabilities.variables[self.metadata['vm']['p']][:])
  self.np = self.fp_probabilities.variables[self.metadata['vm']['rank']].shape[1]
  self.nl = len(self.fp_properties.variables['upper_depth'][:])
  if self.nl > self.metadata['maxnl']:self.nl = self.metadata['maxnl']
  self.nb = len(self.fp_properties.dimensions['bins'])

  #Construct the shapes of all the input variables
  self.variables = {}
  for var in self.metadata['vars']:
   self.variables[var] = {'shape':self.fp_properties[var]['hist'].shape,
                          'units':self.fp_properties[var].units,
                          'description':self.fp_properties[var].description}

  return

 #Draw data per soil class
 def extract_parameters(self,data):
 
  #Define parameters
  undef = -9999
  ids = self.block_ids
  ids = ids[ids >= 0]
  nc = len(ids)#len(data['id'])
  #nd = self.metadata['nd']
  mids = data['id']
  #Construct mapping array
  mapping = np.zeros(np.max(ids)+1).astype(np.int32)
  mapping[ids] = np.arange(nc)
  #Compute the index in the data
  idx = np.in1d(mids,ids)
  #Construct parameters
  m0 = idx & (np.mean(data['hist'],axis=1) != -9999.0)
  m1 = np.in1d(ids,mids[m0])
  #hist
  hist = np.zeros((nc,100))
  hist[:] = -9999.0
  hist[m1] = data['hist'][m0,:]

  #place all the parameters together
  params = {'hist':hist,}
  
  return (mapping,params)

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Extract the parameters
  (mapping,params) = self.extract_parameters(data)

  #Compute the histogram per cell
  hist = np.zeros((self.prob.shape[1],100))
  hist[:] = -9999.0
  for p in xrange(self.rank.shape[1]):
   #Extract the probabilities
   probs = self.prob[:,p]
   #Extract the indices for the classes at this cell
   idx = mapping[self.rank[:,p]]
   m = np.mean(params['hist'][idx,:],axis=1) != -9999.0
   tmp = probs[m,np.newaxis]*params['hist'][idx,:][m]
   #Compute the weighted histograms
   hist[p,:] = np.sum(tmp,axis=0)
 
  '''
  #Calculate the beta parameters
  max_in = params['max']
  ncmax = self.metadata['ncmax']
  min_in = params['min']
  mean_in = params['mean']
  var_in = params['var']
  t0 = time.time() 
  #(min,max,mean,var,alpha,beta) = propr_tools_fortran.compute_parameters(max_in,min_in,mean_in,var_in,
  #                      self.prob,self.rank,mapping,ncmax)
  (min,max,mean,var,alpha,beta) = propr_tools_fortran.compute_parameters_point(max_in,min_in,mean_in,var_in,
                        self.prob,self.rank,mapping,ncmax)
  #print 'creating beta parameters',time.time() - t0
  '''
  
  #Assemble output
  #Assemble output
  output = {}
  output['hist'] = hist
  '''output = {}
  output['mean'] = mean
  output['min'] = min
  output['max'] = max
  output['alpha'] = alpha
  output['beta'] = beta
  output['var'] = var'''

  return output

class initialize:


 #Define functions
 def __init__(self,metadata):

  #Read in the metadata
  self.metadata = metadata#json.load(open(file))

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
  self.nb = len(self.fp_properties.dimensions['bins'])

  #Construct the shapes of all the input variables
  self.variables = {}
  for var in self.metadata['vars']:
   self.variables[var] = {'shape':self.fp_properties[var]['hist'].shape,
                          'units':self.fp_properties[var].units,
                          'description':self.fp_properties[var].description}

  return

 #Draw data per soil class
 def extract_parameters(self,data):
 
  #Define parameters
  undef = -9999
  ids = self.block_ids
  ids = ids[ids >= 0]
  nc = len(ids)#len(data['id'])
  nd = self.metadata['nd']
  mids = data['id']
  #Construct mapping array
  mapping = np.zeros(np.max(ids)+1).astype(np.int32)
  mapping[ids] = np.arange(nc)
  #Compute the index in the data
  idx = np.in1d(mids,ids)
  #Construct parameters
  m0 = idx & (np.mean(data['hist'],axis=1) != -9999)
  m1 = np.in1d(ids,mids[m0])
  #hist
  hist = np.zeros((nc,100))
  hist[:] = undef
  hist[m1] = data['hist'][m0,:]
  '''#max
  max = np.zeros(nc)
  max[:] = undef
  max[m1] = data['max'][m0]
  #min
  min = np.zeros(nc)
  min[:] = undef
  min[m1] = data['min'][m0]
  #mean
  mean = np.zeros(nc)
  mean[:] = undef
  mean[m1] = data['mean'][m0]
  #var
  var = np.zeros(nc)
  var[:] = undef
  var[m1] = data['var'][m0]'''
  #place all the parameters together
  #params = {'min':min,'max':max,'mean':mean,'var':var}
  params = {'hist':hist,}
  
  return (mapping,params)

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Extract the parameters
  (mapping,params) = self.extract_parameters(data)

  #Compute the weighted histograms
  ncmax = self.metadata['ncmax']
  hist = propr_tools_fortran.compute_weighted_histogram(params['hist'],self.prob,self.rank,mapping,ncmax)
 
  '''
  #Compute the histogram per cell
  for x in xrange(self.rank.shape[1]):
   print x
   for y in xrange(self.rank.shape[2]):
    #Extract the probabilities
    probs = self.prob[:,x,y]
    #Extract the indices for the classes at this cell
    idx = mapping[self.rank[:,x,y]]
    #Compute the weighted histograms
    hist[x,y,:] = np.sum(probs[:,np.newaxis]*params['hist'][idx,:],axis=0)
  '''

  #Assemble output
  output = {}
  output['hist'] = hist

  return output
