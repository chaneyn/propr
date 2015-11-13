import propr_tools_fortran
import json
import netCDF4 as nc
import numpy as np
import scipy.stats as stats

def initialize_netcdf_output(instance):

 ud = instance.fp_properties.variables['upper_depth'][:]
 ld = instance.fp_properties.variables['lower_depth'][:]
 fp = nc.Dataset(instance.metadata['output_file'],'w')
 fp.createDimension('z',instance.nl)
 fp.createDimension('x',instance.prob.shape[2])
 fp.createDimension('y',instance.prob.shape[1])
 x = fp.createVariable('x','f8',('x',))
 x.Axis = 'X'
 x[:] = instance.x
 y = fp.createVariable('y','f8',('y',))
 y.Axis = 'Y'
 y[:] = instance.y
 z = fp.createVariable('z','f8',('z',))
 z.Axis = 'Z'
 z[:] = (ud + ld)/2
 z = fp.createVariable('upper_depth','f8',('z',))
 z[:] = ud
 z = fp.createVariable('lower_depth','f8',('z',))
 z[:] = ld

 return fp

def run(file):

 #Initialize the propr class
 instance = initialize(file)

 #Initialize the output file
 print "Initializing output file"
 fp = initialize_netcdf_output(instance)

 #Create the input data for a variable for each layer
 vars = instance.fp_properties.groups.keys()
 for var in vars:

  print "Calculating %s maps" % var
  #Create new variables
  vpc5 = fp.createVariable('%s_pc5' % var,'f4',('z','y','x'))
  vpc95 = fp.createVariable('%s_pc95' % var,'f4',('z','y','x'))
  vmean = fp.createVariable('%s_mean' % var,'f4',('z','y','x'))

  #Iterate through each layer
  for il in xrange(instance.nl):
   #Define the soil property data
   data = {}
   data['id'] = instance.fp_properties.variables['soil_classes'][:]
   data['min'] = instance.fp_properties.groups[var].variables['min'][:,il]
   data['max'] = instance.fp_properties.groups[var].variables['max'][:,il]
   data['mode'] = instance.fp_properties.groups[var].variables['mode'][:,il]

   #Calculate the properties
   output = instance.calculate_properties(data)

   #Output the properties
   vpc5[il,:,:] = output['pc5']
   vpc95[il,:,:] = output['pc95']
   vmean[il,:,:] = output['mean']

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

  #Read the probabilities
  self.prob = self.fp_probabilities.variables['prob'][:]
  self.rank = self.fp_probabilities.variables['rank'][:]
  self.x = self.fp_probabilities.variables['x'][:]
  self.y = self.fp_probabilities.variables['y'][:]
  self.nl = len(self.fp_properties.variables['upper_depth'][:])

  return

 #Draw data per soil class
 def draw_from_distribution(self,data):
 
  nc = len(data['id'])
  nd = self.metadata['nd']
  draws = np.zeros((nc,nd)).astype(np.float32)
  for ic in data['id']:
   loc = data['min'][ic]
   scale = data['max'][ic]-data['min'][ic]
   mode = (data['mode'][ic] - data['min'][ic])/(data['max'][ic] - data['min'][ic])
   tmp = stats.triang.rvs(mode,loc=loc,scale=scale,size=nd)
   draws[ic,:] = tmp

  return draws

 #Calculate the mapped properties
 def calculate_properties(self,data):

  #Create the array of draws
  draws = self.draw_from_distribution(data)

  #Create the array of draws
  array = propr_tools_fortran.assign_draws(self.prob,self.rank,draws)

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

