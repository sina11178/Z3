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
from scipy.sparse.linalg import norm
from scipy.sparse import identity
import time
import cmath
import pandas
import math
import itertools

class model_NearlyZero_Z3_EM:
    '''Simple class to model Z3 parafermions'''
    #def __init__(self, L,  Omega, h, g, Gamma,  deg=2):
    def __init__(self, L ):
        self.L = L + 0
        self.J = 1.0
        self.deg = 3

    def Find_ni(self, Num, i):
        return int(((Num%(self.deg**(i+1)))-(Num%(self.deg**i)))/self.deg**i)
    
    def ni_list(self, myL, mystate):
        res_list = []
        for l in range(myL):
            res_list.append(self.Find_ni(mystate, l))
        return res_list

    ############################################################################################################################################
    # Baxter model for non-periodic case Borrowed form Fendly paper
    #in Z basis
    def H_Baxter(self, myL, myf, myJ, myphi, myhc): #H_pott(myL::Int, myf::Float 1D array, myJ::Float, myh::Float 1D array, myphi::Float, myhc::Float, myg::Float)
        w=np.e**(1.j*2*np.pi/3.0)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate) 
            for n in range(myL-1):
                myMTX[mystate,mystate] = (myMTX[mystate,mystate] - (myJ*(np.e**(1.j*myphi))*(w**(2*ni_all[n] + ni_all[n+1]))))
                myMTX[mystate, self.tau_x(mystate, n)] = myMTX[mystate, self.tau_x(mystate, n)] - myf[n]
            n=myL-1
            myMTX[mystate, self.tau_x(mystate, n)] = myMTX[mystate, self.tau_x(mystate, n)] - myf[n]
        return (myMTX + myhc*np.conjugate(myMTX.T))

    # Baxter model for non-periodic case Borrowed form Fendly paper
    #in x basis
    def x_H_Baxter(self, myL, myf, myJ, myphi, myhc): ##@@@###@@@###@@@###@@@
        w=np.e**(1.j*2*np.pi/3)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate)
            for n in range(myL-1):
                myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf[n]*(w**(ni_all[n])) 
                myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n), n+1)] =  (myMTX[mystate, self.tau_x(self.tau_x_dag(mystate, n), n+1)] - myJ*(np.e**(1.j*myphi)))
            n = myL-1
            myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf[n]*(w**(ni_all[n]))
        return (myMTX + myhc*np.conjugate(myMTX.T))


    ############################################################################################################################################
    # Modified Baxter model for non-periodic case 
    # In Z basis
    def H_modf_Baxter(self, myL, myf, myJ, myphi, myhc): #H_pott(myL::Int, myf::Float 1D array, myJ::Float, myh::Float 1D array, myphi::Float, myhc::Float, myg::Float)
        w=np.e**(1.j*2*np.pi/3.0)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        #print(self.tau_x)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate) 
            for n in range(myL-1):
                myMTX[mystate,mystate] = myMTX[mystate,mystate] - myJ*(np.e**(-1.j*myphi))*(w**(ni_all[n] + 2*ni_all[n+1]))
                myMTX[mystate, self.tau_x(self.tau_x(mystate, n), n)] = myMTX[mystate, self.tau_x(mystate, n)] - myf[n]
            n=myL-1
            myMTX[mystate, self.tau_x(self.tau_x(mystate, n), n)] = myMTX[mystate, self.tau_x(mystate, n)] - myf[n]
        return (myMTX + myhc*np.conjugate(myMTX.T))
    
    # Modified Baxter model for non-periodic case 
    # in X basis
    def x_H_modf_Baxter(self, myL, myf, myJ, myphi, myhc): ##@@@###@@@###@@@###@@@
        w=np.e**(1.j*2*np.pi/3)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate)
            for n in range(myL-1):
                #print("st: ", self.tau_x_dag(self.tau_x(mystate, n), n+1), self.tau_x(self.tau_x_dag(mystate, n), n+1))
                myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf[n]*(w**(ni_all[n]))
                myMTX[mystate, self.tau_x_dag(self.tau_x(mystate, n), n+1)] =  (myMTX[mystate, self.tau_x_dag(self.tau_x(mystate, n), n+1)] - myJ*(np.e**(-1.j*myphi)))
            n = myL-1
            #print("myf:", myf[n]*(w**(2*ni_all[n])), myf[n]*(w**(ni_all[n])), ni_all[n])
            myMTX[mystate, mystate] =  myMTX[mystate,mystate] - myf[n]*(w**(ni_all[n]))
        return (myMTX + myhc*np.conjugate(myMTX.T))


    ############################################################################################################################################
    # Parafermionic models
    # Baxter Parafermionic
    def Baxter_Parafermions(self, myL, myJ, myf, myphi, myhc):
        w=np.e**(1.j*2*np.pi/3)
        Res_List = []
        for l in range(2*myL-1):
            Res_List.append(l%2*((w**2)*-1*myJ*(np.e**(1.j*myphi)))+((l%2)-1)*w*-1*myf[int(l/2)])
        myH = np.diag(Res_List, 1)
        return  myH + myhc*np.conjugate(myH.T)
        
    # Modified Baxter Parafermions
    def ModBaxter_Parafermions(self, myL, myJ, myf, myphi, myhc):
        w=np.e**(1.j*2*np.pi/3)
        J_list = []
        f_list = []
        for l in range(2*myL-1):
            J_list.append(l%2*((w**2)*-1*myJ)*(np.e**(-1.j*myphi)))
            f_list.append(((l%2)-1)*w*-1*myf[int(l/2)])
        myH = np.diag(f_list, 1) + np.diag(J_list, -1)
        return myH + myhc*np.conjugate(myH.T)


    ###########################################################################################################################################
    # Fermionic parity
    def FermiParity_sectors(self, myL):
        myH = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        Parity0 = []
        Parity1 = []
        Parity2 = []
        for mystate in range(self.deg**myL):
            ni_all = self.ni_list(myL, mystate)
            parity = sum(np.array(ni_all))
            if parity%3 == 0:
                Parity0.append(mystate)
            if parity%3 == 1:
                Parity1.append(mystate)
            if parity%3 == 2:
                Parity2.append(mystate)
        #print(Parity0, Parity1, Parity2)
        return np.array(Parity0), np.array(Parity1), np.array(Parity2)

    ############################################################################################################################################
    # Baxter model for non-periodic case and different sectros
    #in x basis
    def x_H_Baxter_sectors(self, mysector, myL, myf, myJ, myphi, myhc): ##@@@###@@@###@@@###@@@
        w=np.e**(1.j*2*np.pi/3)
        myMTX = sparse.lil_matrix((mysector.shape[0], mysector.shape[0]),dtype=complex)
        for i in range(mysector.shape[0]):
            mystate = mysector[i]
            ni_all = self.ni_list(myL, mystate)
            for n in range(myL-1):
                myMTX[i,i] =  myMTX[i,i] - myf[n]*(w**(ni_all[n])) 
                myMTX[i, list(mysector).index(self.tau_x(self.tau_x_dag(mystate, n), n+1))] =  (myMTX[i, list(mysector).index(self.tau_x(self.tau_x_dag(mystate, n), n+1))] - myJ*(np.e**(1.j*myphi)))
            n = myL-1
            myMTX[i,i] =  myMTX[i,i] - myf[n]*(w**(ni_all[n]))
        return (myMTX + myhc*np.conjugate(myMTX.T))











    

    ############################################################################################################################################
    def normed_LeftRighteigens(self, myVR, myVL):
        myVL_norm = np.zeros(myVL.shape,dtype=complex)
        myVR_norm = np.zeros(myVR.shape,dtype=complex)
        diagonals = np.diagonal(np.dot(np.conjugate(myVL).T, myVR))
        for i in range(myVL.shape[0]):
            myVL_norm[:,i] = (myVL[:,i]/np.conjugate(np.sqrt(diagonals[i])))[:]
            myVR_norm[:,i] = (myVR[:,i]/np.sqrt(diagonals[i]))[:]
        return myVR_norm, myVL_norm

    
    def Sigma_z_1(self, myL):
        sigma = np.diag([1,np.exp(1.j*2*np.pi/3.0), np.exp(1.j*4*np.pi/3.0)])
        return np.kron(np.eye(self.deg**(myL-1)), sigma)


    def tau_x(self, state, myl):
        ni = self.Find_ni(state, myl)
        return state - (ni*(self.deg**myl)) + (((ni+1)%3)*(self.deg**myl))


    def tau_x_dag(self, state, myl):
        ni = self.Find_ni(state, myl)
        return state - (ni*(self.deg**myl)) + (((ni+2)%3)*(self.deg**myl))

    

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

    def Commute(self, A,B):
        return (np.dot(A,B)-np.dot(B,A))


    def Krylov_bn(self, myL, myH, myO, myN, mybn, myphin, myetan, myOn, myOrtho_nn):
        HSp_dim = myH.shape[0]
        b = [0.0]; On = []
        phi = [1.0]; eta = [1.0]
        On.append(myO)
        LbO0 = (myH)*myO-myO*myH
        b1 =  np.sqrt((np.conjugate(LbO0.T)*LbO0).diagonal().sum())/(3**(myL/2.0))
        b.append(b1)
        On.append(LbO0/b1)
        Ortho = []
        if mybn.shape[0]>0:
            b = list(mybn); On = list(myOn); Ortho = list(myOrtho_nn) 
            phi = list(myphin); eta = list(myetan)
        ni = np.array(b).shape[0]
        print("ni in recursive: ", ni)
        for n in range(ni,myN+1,1):
            if n>10:
                On = On[-7:]
            if n>7:
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-3]).A)/(3**myL))
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-5]).A)/(3**myL))
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-2]).A)/(3**myL))
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-6]).A)/(3**myL))
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-1]).A)/(3**myL))
                Ortho.append(np.trace((np.conjugate(On[-4].T)*On[-7]).A)/(3**myL))
            #An = (myH*On[n-1]-On[n-1]*myH)-b[n-1]*On[n-2]
            An = ((myH)*On[-1]-On[-1]*myH)-b[-1]*On[-2]
            bn = np.sqrt((np.conjugate(An.T)*An).diagonal().sum())/(3**(myL/2.0))
            On.append(An/bn)
            b.append(bn)
            #if n%10==0:
            print("n= ", n, ", bn: ", np.round(b[-1],5), ", Othogonal test: ", np.round(Ortho[-6:],10))
            if (n+1)%2==0 and n>2:
                eta.append(-1.0*eta[-1]*(b[n-1]/b[n]))
            if (n)%2==0:
                phi.append(-1.0*phi[-1]*(b[n-1]/b[n]))
            np.save("bn_L_"+str(myL)+".npy",b)
            np.save("phi_L_"+str(myL)+".npy",phi)
            np.save("eta_L_"+str(myL)+".npy",eta)
            np.save("On_L_"+str(myL)+".npy", np.array(On))
            np.save("Ortho_nn_L_"+str(myL)+".npy", Ortho)
        return b[1:], phi, eta



    def Gram_Schmidt_bn(self, myL, myH, myO, myN):
        HSp_dim = myH.shape[0]
        On_GS = []
        On_GS.append(np.array((myO.A).T).flatten())
        LbO0 = myH*myO-myO*myH
        On_GS.append(((LbO0.A).T).flatten())
        Ortho = []
        ni = np.array(On_GS).shape[0]
        print("ni in GS: ", ni)
        for n in range(ni,myN+1,1):
            #print("n is: ", n)
            On_GS.append(((np.dot((myH.A),np.reshape(On_GS[-1], (3**myL, 3**myL)).T)-np.dot(np.reshape(On_GS[-1], (3**myL, 3**myL)).T,myH.A)).T).flatten())
            #np.save("GS_On_L_"+str(myL)+".npy", np.array(On_GS))
            #np.save("Ortho_nn_L_"+str(myL)+".npy",Ortho)
        print("curly L loop finished")
        Q, R = np.linalg.qr(np.array(On_GS).T)
        print("Gram-Shmit completed")
        myCurly_L = np.zeros((myN+1, myN+1), dtype=complex)
        #print("Q shape: ", np.array(Q).shape)
        #print("myO: \n", np.round(myO.A,3))
        #print(np.linalg.norm(Q.T[0]), np.linalg.norm(np.array((myO.A).T).flatten()))
        #print("Q: \n", np.round(np.reshape(Q.T[0], (3**myL, 3**myL)),3))
        for n in range(0,myN+1,1):
            for m in range(0,myN+1,1):
                myCurly_L[n,m] = np.trace(np.dot(np.conjugate(np.reshape(Q.T[n], (3**myL, 3**myL)).T), np.dot((myH.A),np.reshape(Q.T[m], (3**myL, 3**myL)))-np.dot(np.reshape(Q.T[m], (3**myL, 3**myL)),myH.A)))
        #print("Curly_L: \n", np.round(myCurly_L,3))
        return abs(myCurly_L)
            
    

    def h0_htild0(self, mybn):
        max_n = np.array(mybn).shape[0]
        print("max_n: ", max_n)
        myhn = []
        myhtildn = []
        for nn in range(max_n-1):
            myhn.append(0.5*(mybn[nn]+mybn[nn+1]))
            myhtildn.append(((-1)**nn)*0.5*(mybn[nn]-mybn[nn+1]))
        return myhn, myhtildn 


    def AutoCorr_bn(self, myL, myH, myO, myN, myt_steps, myt):
        O_t = myO[:]
        A_inf_GS = []
        for t_step in np.arange(1, 4, 1.0):
            if t_step-1 in myt_steps:
                A_inf_GS.append(np.trace((np.conjugate(myO.T)*O_t).A)/(3**myL))
            print("t = ", t_step)
            print("O_t: \n", np.round(O_t, 3))
            mydt = np.round(1.0/10,3)
            for i in range(10):
                O_dt = sparse.lil_matrix((3**myL, 3**myL),dtype=complex)
                LnO = []
                LnO = O_t[:]
                #exp_test = 0
                for n in range(0,myN+1,1):
                    #exp_test = exp_test +  ((1.j*mydt)**n)/math.factorial(n)
                    O_dt = (O_dt + LnO*(((1.j*mydt)**n)/math.factorial(n)))[:]
                    LnO = ((myH)*LnO-LnO*myH)[:]
                    #print("LnO: \n", np.round(LnO.A, 3))
                    #LnO.append(LnO_res)
                #print("exp_test : ", exp_test, ", exp(i*0.01) = ", np.exp(1.j*mydt))
                O_t = O_t + O_dt
            #np.save("O_t_"+str(myL)+".npy", O_t)
        return A_inf_GS

    def FermiParity(self, myL):
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        for j in range(3**myL):
            newstat = j+0.0
            for l in range(myL):
                newstat = self.tau_x(self.tau_x(newstat, l), l)
            myMTX[j, newstat] = 1
        return myMTX


    def Level_stat(self, myEs):
        myRn_list = []
        #print("myEs: ", myEs)
        for n in range(1,int(np.array(myEs).shape[0])-1, 1):
            E1 = myEs[n-1]
            E2 = myEs[n]
            E3 = myEs[n+1]
            #print("E1: ", E1, ", E2: ", E2,", E3: ", E3, ", position: ", int(np.array(myEs).shape[0]/3))
            myRn_list.append(min([E3-E2, E2-E1])/max([E3-E2, E2-E1]))
        return (myRn_list)
    
    
    def PasitySec_LevelStat(self, myU, myV, myL, myReturn):
        sector0_E = []
        sector1_E = []
        sector2_E = []
        for FP, E in  zip(np.round((np.log(np.array(self.FermiParity_sectors(myV, myL)).flatten())/(2*1.j*np.pi/3)).real,5), myU):
            #print("FP , E = ", FP, np.round(E,20))
            if FP != 0 and FP != 1 and FP!=2:
                if np.round(E,10) not in np.round(sector1_E,10):
                    sector1_E.append(E)
                if np.round(E,10) not in np.round(sector2_E,10):
                    sector2_E.append(E)
            if FP == 0 and E not in sector0_E:
                sector0_E.append(E)
            if FP == 1 and E not in sector1_E:
                sector1_E.append(E)
            if FP == -1 and E not in sector2_E:
                sector2_E.append(E)
        #print("sector0_E: ", sorted(np.round(abs(np.array(sector0_E)),3)))#min(abs(np.array(sector0_E))))
        #print("sector1_E: ", sorted(np.round(abs(np.array(sector1_E)),3)))#min(abs(np.array(sector1_E))))
        #print("sector2_E: ", sorted(np.round(abs(np.array(sector2_E)),3)))#min(abs(np.array(sector2_E))))
        if myReturn == 'LStat':
            return self.Level_stat(sorted(np.array(sector0_E))), self.Level_stat(sorted(np.array(sector1_E))), self.Level_stat(sorted(np.array(sector2_E)))
        if myReturn == 'Sspace' and myL > 2:
            return sorted(np.array(sector0_E))[int(3**(myL-2))+1:2*int(3**(myL-2))+1] - sorted(np.array(sector0_E))[int(3**(myL-2)):2*int(3**(myL-2))], sorted(np.array(sector1_E))[int(3**(myL-2))+1:2*int(3**(myL-2))+1] - sorted(np.array(sector1_E))[int(3**(myL-2)):2*int(3**(myL-2))], sorted(np.array(sector2_E))[int(3**(myL-2))+1:2*int(3**(myL-2))+1] - sorted(np.array(sector2_E))[int(3**(myL-2)):2*int(3**(myL-2))]
        if myReturn == 'Sspace' and myL <= 2:
            return sorted(np.array(sector0_E))[1:] - sorted(np.array(sector0_E))[:-1], sorted(np.array(sector1_E))[1:] - sorted(np.array(sector1_E))[:-1], sorted(np.array(sector2_E))[1:] - sorted(np.array(sector2_E))[:-1]



    def alpha_j(self, myL, myl):
        print('myl: ', myl)
        w=np.e**(1.j*2*np.pi/3.0)
        myMTX = sparse.lil_matrix((self.deg**myL, self.deg**myL),dtype=complex)
        if myl%2!=0:
            for mystate in range(self.deg**myL):
                new_state = mystate
                site = int((myl-1)/2)
                for site in range(int((myl-1)/2)+1):
                    #print('site: ', site)
                    new_state = self.tau_x(new_state, site)
                #print('last site: ', site+1)
                myMTX[mystate, new_state] = w**self.Find_ni(new_state, site)
        if myl%2==0:
            for mystate in range(self.deg**myL):
                new_state = mystate
                if myl>0:
                    for site in range(int((myl)/2)):
                        #print('In alpha site: ', site, ", myl = ", myl)
                        new_state = self.tau_x(new_state, site)
                #print('last site: ',  int((myl)/2))
                myMTX[mystate, new_state] = w**self.Find_ni(new_state, int(myl/2))
        return myMTX

    
    def Psi_Left(self, myL, myf, myJ, myphi):
        w=np.e**(1.j*2*np.pi/3.0)
        myPsi_L = self.alpha_j(myL, 0)
        for site in range(2, 2*myL, 2):
            print("site: ", site)#, ", Psi_L: \n", np.round((self.alpha_j(myL, site)).A,3))
            myPsi_L = myPsi_L + self.alpha_j(myL, site)*(((myf*w)/(myJ*np.exp(1.j*myphi)))**(site/2))
        print(np.round(myPsi_L.A,2))
        return myPsi_L


    def find_closest_numbers(self, array1, SSH_n, size):
        myArg = abs(array1).argmin()
        print("myArg: ", myArg)
        if myArg>=SSH_n-size+1:
            return np.arange(SSH_n,myArg-size, -1)
        if myArg<=size-1:
            return np.arange(0,myArg+size, 1)
        if myArg>size-1 and myArg<SSH_n-size+1:
            return np.arange(myArg-size,myArg+size, 1)


    def complex_sort(self, A):
        B = np.vstack((np.round(A,12).real, np.round(A,12).imag, A)).T
        #print("B:\b", np.round(B,3))
        sorted_B = np.array(sorted(B,key=lambda x: (x[0],x[1])))
        #print("sorted_B:\b", np.round(sorted_B, 3))
        sorted_A = [_[2] for _ in sorted_B]
        return sorted_A
