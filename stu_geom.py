import FluidChannel as fc
import numpy as np

ocStu = fc.FluidChannel(Lx_p = 2., Ly_p = 2., Lz_p = 8.,
                        N_divs = 19,
                        obst = fc.EmptyChannel(2))
ocStu.write_mat_file('openChannel')



