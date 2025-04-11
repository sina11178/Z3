'''

WE WILL CALL TO RUN RESPECTIVE SIMULATION HERE

'''

from f_phi_autocorr_SIM import f_phi_autocorrelator_SIMULATION_EXCEL, f_phi_autocorrelator_SIMULATION_NPY
from size_comparison_autocorr_SIM import size_autocorrelator_SIMULATION_EXCEL, size_autocorrelator_SIMULATION_NPY
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import sys


array_parameters = [[.1, (.7*cm.pi)/6], [.1, (.8*cm.pi)/6], [.1, (.9*cm.pi)/6],
                    [.1, (cm.pi)/6], [.1, (1.1*cm.pi)/6],
                    [.1, (1.2*cm.pi)/6], [.1, (1.3*cm.pi)/6]]
L = 4
array=[.1, 1.3*(cm.pi)/6]
L_array = [4, 5, 6, 7, 8]

#size_autocorrelator_SIMULATION_NPY(0.05, 100000, 10, L_array, array, 1)
#f_phi_autocorrelator_SIMULATION_NPY(0, 30, .5, L, [array_parameters[1]], 1)

# FOR RUNNING ON CLUSTER

size_autocorrelator_SIMULATION_NPY(0, 30, .5, L_array, array_parameters[int(sys.argv[1])], 1)