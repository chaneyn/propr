#First compile the source code
import os
cwd = os.getcwd()
os.chdir('src')
os.system('python compile.py')
os.chdir(cwd)

#Write the libraries
from distutils.core import setup
setup(name='propr',
      version='1.0',
      package_dir={'propr': 'libraries'},
      packages=['propr'],
      package_data={'propr': ['*.py','*.so']}
      )
