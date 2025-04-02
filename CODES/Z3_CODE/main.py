'''

WE WILL CALL TO RUN RESPECTIVE SIMULATION HERE

'''

from f_phi_autocorr_SIM import f_phi_autocorrelator_SIMULATION_EXCEL, f_phi_autocorrelator_SIMULATION_NPY
from size_comparison_autocorr_SIM import size_autocorrelator_SIMULATION_EXCEL, size_autocorrelator_SIMULATION_NPY
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd


array_parameters = [[.1, 1.3*(cm.pi)/6]]

f_phi_autocorrelator_SIMULATION_NPY(0.05, 100000, 10, 4, array_parameters, 1)
