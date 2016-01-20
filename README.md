#Installation
1. Install Conda (Python 2.7)
	* Download the installer from http://conda.pydata.org/miniconda.html
	* Set the default python environment to this one
2. Install all the required libraries (e.g., conda install numpy)
	* numpy
	* netCDF4
	* scipy
3. Install propr
	* python setup.py install

#Metadata
The metadata.json file contains all the metadata necessary to run the software. Avoid changing the python code at all costs. Put your data into the formats that this approach requires. I am including a brief explanation for each parameter in the metadata.json file below:

1. input_file - This contains two 3d arrays. Dimension 1 is the id and dimensions 2 and 3 are the latitude and longitude dimensions respectively. One array is the probability array and the other is the soil class id array. 

2. output_file - This is what the algorithm outputs. It inclueads the pc5, pc95, mean for all the variables that you are mapping

3. properties_file - This has all the properties information that is associated to each soil class

4. nd - Number of samples per soil class

5. vm - This sets parameters and maps the names of the input variables to what the algorithm expects
	* rank - soil class ids ranked according to probabilities. 
	* prob - probability of the corresponding soil class
	* x - longitude
	* y - latitude

6. seed - Numpy random seed

7. nblocks - Split the domain into nblocks by nblocks. This limits that amount of memory that needs to be used. 

8. maxnl - Number of layers to map. Starting from the surface. 

9. vars - List of variables that you wish to map. These variables must be available in the properties file.

10. ncmax - Maximum number of soil classes to include.

11. maxrank - Maximum rank to include

#Tutorial

1. Download the data (wget https://www.dropbox.com/s/fyx6xtdoghqn7lx/propr_tutorial.tar.gz?dl=0)

2. Untar (tar -xvzf propr_tutorial.tar.gz)

3. Run the algorithm (python driver.py)

4. Plot the results (cd data & python test.py)



