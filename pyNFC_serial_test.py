#!/usr/bin/env python
"""
 testing serial implementation of pyNFC
"""


import numpy as np
import math
import time


import sys
sys.path.insert(1,'.')
from vtkHelper import saveVelocityAndPressureVTK_binary as writeVTK
#import pyNFC
import FluidChannel as fc
from pyNFC_serial import pyNFC_serial_problem as Prob


# Information about the LBM run that produced the data - I should get this from params.lbm
# Read data from params.lbm
input_file_name = 'params.lbm'
input_data = open(input_file_name,'r')
lattice_type = str(input_data.readline().rstrip()) # perhaps I can actually use this as intended now
Num_ts = int(input_data.readline())
ts_rep_freq = int(input_data.readline())
Warmup_ts = int(input_data.readline())
plot_freq = int(input_data.readline())
Cs = float(input_data.readline())
rho_lbm = float(input_data.readline())
u_lbm = float(input_data.readline())
omega = float(input_data.readline())
Nx = int(input_data.readline())
Ny = int(input_data.readline())
Nz = int(input_data.readline())
Restart_flag = int(input_data.readline())
Lx_p = float(input_data.readline())
Ly_p = float(input_data.readline())
Lz_p = float(input_data.readline())
t_conv_fact = float(input_data.readline())
l_conv_fact = float(input_data.readline())
p_conv_fact = float(input_data.readline())
input_data.close()

# initialize the pyNFC_serial_problem object
myProb = Prob(Nx,Ny,Nz,rho_lbm,u_lbm,omega,Cs,lattice_type)
print "Problem Intialized!"


# each process initialize their partition:
#myPart = pyNFC.NFC_LBM_partition(rank,size,comm,Nx,Ny,Nz,rho_lbm,u_lbm,omega,Cs,lattice_type)
#myPart.write_node_sorting() # should this be done in the constructor?

# do some time stepping
#numTs = 10
#plot_freq = 5;
#if rank == 0:
#    time1 = time.time()
#
#for ts in range(Num_ts):
#
#    print "executing time step %d" % ts
#    isEven = (ts%2 == 0)
##    myPart.take_LBM_timestep(isEven)
#
#    if ((ts > 0) and (ts % plot_freq == 0)):
##        myPart.write_data(isEven)
#    
#
#time2 = time.time()
#ex_time = time2 - time1
#print "approximate execution time is %g seconds" % (ex_time)
#numLP = Nx*Ny*Nz
#LPUs = numLP*Num_ts/(ex_time)
#print "approximate LPU/sec = %g " % LPUs
