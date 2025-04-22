import numpy as np
from model_NearlyZero_Z3_EM import *
from scipy import sparse
from scipy.sparse import identity
import scipy
from scipy.sparse.linalg import eigs
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import itertools
import sys
import os



L = int(sys.argv[1])
J = float(sys.argv[2])
f = float(sys.argv[3])
phi = float(sys.argv[4])*(np.pi/6)
Nt = int(sys.argv[5])
T = float(sys.argv[6])



model1 = model_NearlyZero_Z3_EM(L)
'''
#testing the model
H = model1.H_Pott(L, f, J, phi)
x_H = model1.x_H_Pott(L, f, J, phi)
U, V = np.linalg.eig(H.todense())
xU, xV = np.linalg.eig(x_H.todense())
print(max(abs(np.array(sorted(U))-np.array(sorted(xU)))))
'''

time_steps = model1.Saeed_sequence(Nt+1, 15, 2.6)
np.save("Time_steps.npy", time_steps)


A_t = []
if os.path.exists("A_L_"+str(L)+".npy") == True:
    A_t = np.load("A_L_"+str(L)+".npy")
    if A_t.shape[0] >= Nt+1:
        print("this job was completed")
        sys.exit()
    if A_t.shape[0] < Nt+1:
        time_steps = time_steps[A_t.shape[0]:]

A_t = list(A_t)
H = model1.H_Pott(L, f, J, phi)
print("Hermitian test of H:", np.matrix.max(abs(H.todense()-np.conjugate(H.todense().T))))
U, V = np.linalg.eigh(H.todense())
for nt in time_steps:
    UT = np.dot(V,np.dot(np.diag(np.exp(-1.j*(nt*T)*U)),np.conjugate(V.T)))
    A_t.append(np.trace(np.dot(np.dot(np.dot(np.conjugate(UT.T), np.conjugate(model1.Sigma_z_1(L))),UT),model1.Sigma_z_1(L)))*(1.0/(3.**L)))
    np.save("A_L_"+str(L)+".npy",A_t)
    print("nt=", nt, " is done!")

    
print("######################################")
print("Calculation detail:"+"\nL = ", L, ", J= ", J, ", f= ", f, ", phi= ", phi, "T=", T)
print("######################################")
