import os
os.system('rm -rf *.dSYM')
#subroutines
subroutines = 'assign_draws'

#Create library
cmd = 'f2py -c only: %s : -m propr_tools_fortran propr_tools.f90 -lgomp --fcompiler=gnu95 --f90flags="-w -fopenmp -O3"' % subroutines
print cmd
os.system(cmd)

#Move to the previos directory
os.system('mv propr_tools_fortran.so ../libraries/.')
