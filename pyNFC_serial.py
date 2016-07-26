#pyNFC_serial.py
"""
  implementation file for Python classes for serial pyNFC
"""
import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTKpt
import pyLattices as pl


class pyNFC_serial_problem(object):
    def __init__(self,Nx,Ny,Nz,rho_lbm,u_bc,omega,Cs,lattice_type = 'D3Q15'):
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz;
        self.rho_lbm = rho_lbm; self.u_bc = u_bc;
        self.omega = omega;
        self.Cs = Cs;
        self.lattice_type = lattice_type
        
        if lattice_type == 'D3Q15':
            self.lattice = pl.D3Q15Lattice(self.Nx,self.Ny,self.Nz)