import numpy as np
import cmath
from random import randint
import random
import scipy.misc
import scipy.sparse as sparse
#from sympy.physics.quantum import TensorProduct
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm, expm_multiply
from numpy import linalg
from scipy import sparse
from scipy.sparse import identity
import time
import cmath
import itertools

class model_NearlyZero_Z3_EM:
    '''Simple class to model Z3 time crystal'''
    #def __init__(self, L,  Omega, h, g, Gamma,  deg=2):
    def __init__(self, L ):
        self.L = L + 0
        self.J = 1.0
        self.deg = 3
##################################################
## HEZ Hamiltonian in Tao_z basis
##################################################

    def Find_ni(self, Num, i):
        return ((Num%(self.deg**(i+1)))-(Num%(self.deg**i)))/self.deg**i
    
    def ni_list(self, myL, mystate):
        res_list = []
        for l in range(myL):
            res_list.append(self.Find_ni(mystate, l))
        return res_list

    # stolen from "Physical Review B 90.16 (2014): 165106" equation (4)
    def H_Pott(self, myL, myf, myJ, myphi):
        w=np.e**(1.j*2*np.pi/3.0)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate) 
            for n in range(myL-1):
                myMTX[mystate,mystate] = myMTX[mystate,mystate] - myJ*(np.e**(1.j*myphi))*(w**(ni_all[n] +2*ni_all[n+1]))  - myJ*(np.e**(-1.j*myphi))*(w**(2*ni_all[n] +ni_all[(n)+1]))
                myMTX[mystate, self.tau_x_dag(mystate, n)] = myMTX[mystate, self.tau_x_dag(mystate, n)] -myf
                myMTX[mystate, self.tau_x(mystate, n)] = myMTX[mystate, self.tau_x(mystate, n)]-myf
            n=myL-1
            myMTX[mystate, self.tau_x_dag(mystate, n)] = myMTX[mystate, self.tau_x_dag(mystate, n)] -myf
            myMTX[mystate, self.tau_x(mystate, n)] = myMTX[mystate, self.tau_x(mystate, n)]-myf
        return myMTX

    def Sigma_z_1(self, myL):
        sigma = np.diag([1,np.exp(1.j*2*np.pi/3.0), np.exp(1.j*4*np.pi/3.0)])
        #sigma = np.array([[0,0,1.],[1.,0,0],[0,1.,0]])
        #print "trace in function:", np.trace(sigma)
        return np.kron(np.eye(self.deg**(myL-1)), sigma)


    def tau_x(self, state, myl):
        ni = self.Find_ni(state, myl)
        return state - (ni*(self.deg**myl)) + (((ni+1)%3)*(self.deg**myl))


    def tau_x_dag(self, state, myl):
        ni = self.Find_ni(state, myl)
        return state - (ni*(self.deg**myl)) + (((ni+2)%3)*(self.deg**myl))

    
    def x_H_Pott(self, myL, myf, myJ, myphi):
        w=np.e**(1.j*2*np.pi/3)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate)
            for n in range(myL-1):
                myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf*(w**(ni_all[n]) + w**(2*ni_all[n]))
                myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n+1), n)] = myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n+1), n)] - myJ*(np.e**(1.j*myphi))
                myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n), n+1)] = myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n), n+1)] - myJ*(np.e**(-1.j*myphi))
            n = myL-1
            myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf*(w**(ni_all[n]) + w**(2*ni_all[n]))
        return myMTX



    def Saeed_sequence(self, myn, mys, mym):
        sequ=[]
        i=0
        for x in range(myn):
            if x < mys:       
                i = (i + x)
            if x >= mys:      
                i = (i + x)*mym
            sequ.append(int(i))
        return sequ


