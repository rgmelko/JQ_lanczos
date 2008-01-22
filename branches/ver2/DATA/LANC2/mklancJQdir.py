#!/usr/bin/python
# This python script creates a directory hierarchy with the correct input files


#from Numeric import *
import os
import numpy

execut_pre = 'lancJQ_'
Pi=4.0*numpy.arctan(1.0)

Lcl = '26A'
Sz = 4
Nindex = 200

index = numpy.arange(0,Nindex,1)
	
L_dir_name = 'L' + str(Lcl)
executable = execut_pre + str(Lcl)

# open the condor submit file 
sub_name = L_dir_name + os.sep +'shell.submit' 
mk_dir_name = 'mkdir -p ' +  L_dir_name + os.sep
try:
	cfile=open(sub_name,'w') 
except:
	os.system(mk_dir_name)
	cfile=open(sub_name,'w') 


for cindex in index:		
	# Creates the required input directory names and path
	job_dir_name = 'theta'+str(cindex) +'of'+str(Nindex) + os.sep + 'Sz' + str(Sz) 
	
	# write the approprite lines to the condor submit file
	changedir = 'cd ' + job_dir_name
	cfile.write(changedir+'\n')
	
	cfile.write(executable+' > data.out\n')
	
	dir_name = L_dir_name + os.sep + job_dir_name
	mk_dir_name = 'mkdir -p ' +  dir_name
	os.system(mk_dir_name)
	
	# Create the input file param.dat
	file_name = dir_name + os.sep + 'param.dat'
	pfile = open(file_name,'w')
	cJ=numpy.cos(2.0*Pi*float(cindex)/float(Nindex))
	cQ=numpy.sin(2.0*Pi*float(cindex)/float(Nindex))
	pfile.write('%f\n%f\n%d\n8\n1\n' % (cJ,cQ,Sz)) 
	pfile.close
	
	cfile.write('cd ../../\n')
	
cfile.close        
        
        
