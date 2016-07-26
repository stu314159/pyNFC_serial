#pyLattices.py
"""
  implementation file for pyLattices object.
  note that this is different from the pyLattice object
  developed for pyNFC.  This will be a serial version
  that, perhaps, can later be parallelized
"""

import numpy as np

class Lattice(object):
    
    def __init__(self,Nx,Ny,Nz):
        """
           basic class constructor
        """
        
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz
        self.ex = []; self.ey = []; self.ez = []
        self.bbSpd = []; self.w = [];
        
        
        
class D3Q15Lattice(Lattice):
    """
      D3Q15 Lattice
    """
    def __init__(self,Nx,Ny,Nz):
        """
          D3Q15 Lattice
        """
        super(D3Q15Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1]; self.ex = np.array(self.ex,dtype=np.float32);
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1]; self.ey = np.array(self.ey,dtype=np.float32);
        self.ez = [0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1]; self.ez = np.array(self.ez,dtype=np.float32);

        self.bbSpd = [0,2,1,4,3,6,5,14,13,12,11,10,9,8,7]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
            1./72.,1./72.,1./72.,1./72.,
	    1./72.,1./72.,1./72.,1./72.]
        self.w = np.array(self.w,dtype=np.float32)
        
        
class D3Q19Lattice(Lattice):
    def __init__(self,Nx,Ny,Nz):
        super(D3Q19Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0]; self.ex = np.array(self.ex,dtype=np.float32)
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1]; self.ey = np.array(self.ey,dtype=np.float32)
        self.ez = [0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1]; self.ez = np.array(self.ez,dtype=np.float32)
        self.bbSpd = [0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
	       1./72.,1./72.,1./72.,1./72.,
	       1./72.,1./72.,1./72.,1./72.]
        self.w = np.array(self.w,dtype=np.float32)