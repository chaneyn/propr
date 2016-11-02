import propr_tools_fortran
import json
import netCDF4 as nc
import numpy as np
import scipy.stats as stats
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
 #vars = instance.fp_properties.groups.keys()
 vars = instance.metadata['vars']
 for var in vars:

  #Create new variables
  fp.createVariable('%s_pc5' % var,'f4',('z','y','x'))
  fp.createVariable('%s_pc95' % var,'f4',('z','y','x'))
  fp.createVariable('%s_mean' % var,'f4',('z','y','x'))

 return fp

def calculate_properties_on_each_block(instance,fp,ix,iy):

 #Define the arrays
 mr = instance.metadata['maxrank']
 instance.prob = instance.fp_probabilities.variables[instance.metadata['vm']['prob']][0:mr,iy,ix]
 instance.prob = instance.metadata["vm"]["prob_scalar_multiplier"]*instance.prob
 instance.rank = instance.fp_probabilities.variables[instance.metadata['vm']['rank']][0:mr,iy,ix]

 #Curate the probabilities (always needs to sum to 100)
 instance.prob = instance.prob/np.sum(instance.prob,axis=0)
 
 #Create the input data for a variable for each layer
 #vars = instance.fp_properties.groups.keys()
 vars = instance.metadata['vars']
 for var in vars:

  print "Calculating %s maps" % var
  #Retrieve the variables
  vpc5 = fp.variables['%s_pc5' % var]
  vpc95 = fp.variables['%s_pc95' % var]
  vmean = fp.variables['%s_mean' % var]

  #Iterate through each layer
  for il in xrange(instance.nl):

   #Define the soil property data
   data = {}
   data['id'] = instance.fp_properties.variables['id'][:]
   data['min'] = instance.fp_properties.groups[var].variables['min'][:,il]
   data['max'] = instance.fp_properties.groups[var].variables['max'][:,il]
   data['mode'] = instance.fp_properties.groups[var].variables['mode'][:,il]

   #Calculate the properties
   output = instance.calculate_properties(data)

   #Output the properties
   vpc5[il,iy,ix] = output['pc5']
   vpc95[il,iy,ix] = output['pc95']
   vmean[il,iy,ix] = output['mean']

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

  return

 #Draw data per soil class
 def draw_from_distribution(self,data):
 
  ids = np.ma.getdata(np.unique(self.rank))
  ids = ids[ids >= 0]
  nc = len(ids)#len(data['id'])
  nd = self.metadata['nd']
  draws = np.zeros((nc,nd)).astype(np.float32)
  mapping = np.zeros(np.max(data['id'])+1).astype(np.int32)
  ic = 0
  for id in ids:
   rel_min = 100*np.abs((data['mode'][id] - data['min'][id])/data['mode'][id])
   rel_max = 100*np.abs((data['mode'][id] - data['max'][id])/data['mode'][id])
   if data['mode'][id] == -9999:
    draws[ic,:] = -9999
    mapping[id] = ic
   elif data['max'][id] == data['min'][id]:
    draws[ic,:] = data['mode'][id]
    mapping[id] = ic
   elif((rel_min < 0.01) & (rel_max < 0.01)):
    draws[ic,:] = data['mode'][id]
    mapping[id] = ic
   else:
    loc = data['min'][id]
    scale = data['max'][id]-data['min'][id]
    mode = (data['mode'][id] - data['min'][id])/(data['max'][id] - data['min'][id])
    np.random.seed(self.metadata['seed'])
    tmp = stats.triang.rvs(mode,loc=loc,scale=scale,size=nd)
    draws[ic,:] = tmp
    mapping[id] = ic
   ic += 1

  return (draws,mapping)

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Create the array of draws
  (draws,mapping) = self.draw_from_distribution(data)

  #Create the array of draws
  ncmax = self.metadata['ncmax']
  array = propr_tools_fortran.assign_draws(self.prob,self.rank,draws,mapping,ncmax)

  #Sort the data
  array = np.sort(array,axis=0)

  #Compute the 5th pct, 95th pct, and mean
  ipc5 = np.percentile(np.arange(self.metadata['nd']),5,interpolation='nearest')
  ipc95 = np.percentile(np.arange(self.metadata['nd']),95,interpolation='nearest')
  output = {}
  output['mean'] = np.mean(array,axis=0)
  output['pc5'] = array[ipc5,:,:]
  output['pc95'] = array[ipc95,:,:]

  return output

