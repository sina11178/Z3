'''
PYTHON CODE FOR Z3 PROJECT

REFERENCE FOR HAMILTONIAN: https://arxiv.org/pdf/1407.6376

Z3 symmetry --> N = 3
'''
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 1
phi = 1
f = 1
L = 1  # SYSTEM SIZE
gamma = 1  # DETERMINES HERMICITY

w = np.exp((2* cm.pi * 1.j)/3)  # For defining Sigma operator

sigma = [
    [1, 0, 0],
    [0, w, 0],
    [0, 0, w**2]
]

clock = [
    [0, 0, 1],
    [1, 0, 0],
    [0, 1, 0]
]

I = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]

# CONVERT MATRICES INTO NP MATRICES

np_SIG = np.matrix(sigma)
np_CLK = np.matrix(clock)
np_I = np.matric(I)