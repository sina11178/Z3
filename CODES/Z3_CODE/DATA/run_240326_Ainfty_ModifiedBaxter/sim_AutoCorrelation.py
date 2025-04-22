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
#import seaborn as sns
from scipy.fft import fft, fftfreq



L = int(sys.argv[1])
J = float(sys.argv[2])
f = float(sys.argv[3])
phi = float(sys.argv[4])*(np.pi/6)
Nt = int(sys.argv[5]) #number of times selection to be calculated
Finding = str(sys.argv[6]) #'yes' --> Krylov yes, 'no' --> Krylov no, 'GS, Gram-Schmit, 'Diagonal' --> diagonalize
Wf = float(sys.argv[7])
hc = float(sys.argv[8])
SSH_n = int(sys.argv[9])
Model = str(sys.argv[10])
T = 1.0
print("L=", L, ", Wf = ", Wf)
#####################
#L = int(sys.argv[1])
#J = float(sys.argv[2])
#f = float(sys.argv[3])
#phi = float(sys.argv[4])*(np.pi/6)
#Nt = int(sys.argv[5]) #number of times selection to be calculated
#T = float(sys.argv[6]) #Period
#Krylov = str(sys.argv[7])



model1 = model_NearlyZero_Z3_EM(L)

Dis_f = [_+f for _ in Wf*(np.random.random(L)*2.-1.)]


### testing the Baxter model ###
H = model1.H_Baxter(L, Dis_f, J, phi, hc)
x_H = model1.x_H_Baxter(L, Dis_f, J, phi, hc)
Ub, Vb = np.linalg.eig(H.todense())
xUb, xVb = np.linalg.eig(x_H.todense())
#print("Hb:\n ", np.round(H.A, 3))
#print("x_Hb:\n ", np.round(x_H.A, 3))
print("Testing the [Spin] Baxter model: ", max(abs(np.array(sorted(Ub))-np.array(sorted(xUb)))))


### testing the modified Baxter model ###
Hmb = model1.H_modf_Baxter(L, Dis_f, J, phi, hc)
x_Hmb = model1.x_H_modf_Baxter(L, Dis_f, J,  phi, hc)
Umb, Vmb = np.linalg.eig(Hmb.todense())
xUmb, xVmb = np.linalg.eig(x_Hmb.todense())
#print("Hmb:\n ", np.round(Hmb.A, 3))
#print("x_Hmb:\n ", np.round(x_Hmb.A, 2)) 
#print(np.round(model1.complex_sort(Umb), 3))
#print(np.round(model1.complex_sort(xUmb), 3))
print("Testing the [Spin] modified Baxter model: ", max(abs(np.array(model1.complex_sort(Umb))-np.array(model1.complex_sort(xUmb)))))
#print("Comparing [Spin] modified Baxter model and Baxter model: ", max(abs(np.array(sorted(Ub))-np.array(sorted(Umb)))))

### Comparing modified Baxter with Baxter(should be hc=1) 
Ub, Vb = np.linalg.eigh(H.todense())
Umb, Vmb = np.linalg.eigh(Hmb.todense())
print("Comparing the [Spin] modified Baxter model with Baxter model: ", max(abs(np.array(sorted(Umb))-np.array(sorted(Ub)))))


'''
## comparing Baxter and modified Baxter in Parafermion picture ###
# here they are hermitian: hc = 1
Hbp =  model1.Baxter_Parafermions(L, J, Dis_f, phi, 1)
Ubp, Vbp = np.linalg.eigh(Hbp)
#print("Hbp:\n", np.round(Hbp,3))
Hmbp = model1.ModBaxter_Parafermions(L, J, Dis_f, phi, 1)
Umbp, Vmbp = np.linalg.eigh(Hmbp)
#print("Hmbp:\n", np.round(Hmbp,3))
print("Comparing [Parafermion] modified Baxter model and Baxter model: ", max(abs(np.array(sorted(Ubp))-np.array(sorted(Umbp)))))
'''



"""
### Checking Left and right eigenvalues ###
H = model1.H_Baxter(L, Dis_f, J, phi, hc).A
print("Hdagger*H - H*Hdagger = ", max(abs(np.dot(H, np.conjugate(H).T) - np.dot(np.conjugate(H).T, H)).flatten()))
UR, VR = np.linalg.eig(H)
UL, VL = np.linalg.eig(np.conjugate(H).T)
print("UR:\n", np.round(UR,3))
print("UL:\n", np.round(UL,3))
print("VL.VR:\n", np.diagonal(np.round(np.dot(np.conjugate(VL).T, VR),10)))
print("VLn.UR.VRn-H: ", max(abs(np.dot(VR, np.dot(np.diag(UR), np.conjugate(VL).T))-H).flatten()))
VRn, VLn = model1.normed_LeftRighteigens(VR, VL)
print("VLn.VRn: ", max(abs(np.dot(np.conjugate(VLn).T, VRn)-np.eye(3**L)).flatten()))
print("VLn.VRn:\n", (np.round(np.dot(np.conjugate(VLn).T, VRn),10)))
"""


time_steps = [int(np.e**x) for x in range(Nt)]
np.save("Time_steps.npy", time_steps)


A_t = []
A_t_nom = []
A_t_den = []
if os.path.exists("A_L_"+str(L)+".npy") == True:
    A_t = np.load("A_L_"+str(L)+".npy")
    if hc != 1 and Finding == 'Ainfty_ED':
        A_t_nom = list(np.load("A_nom_L_"+str(L)+".npy"))
        A_t_den = list(np.load("A_den_L_"+str(L)+".npy"))
    if A_t.shape[0] >= Nt+1:
        print("this job was completed")
        sys.exit()
    if A_t.shape[0] < Nt+1:
        time_steps = time_steps[A_t.shape[0]:]

#n_LnO_ = [0]
#O_t_ =  sparse.lil_matrix((3**L, 3**L),dtype=complex)
#LnO_ = []; LnO_.append(sparse.lil_matrix((3**L, 3**L),dtype=complex))
#if os.path.exists("n_LnO_L_"+str(L)+".npy") == True and Finding == "bn":
#    n_LnO_ = list(np.load("n_LnO_L_"+str(L)+".npy"))
#    O_t_ = np.load("O_t_"+str(L)+".npy", allow_pickle=True)
#    LnO_ = list(np.load("LnO_"+str(L)+".npy", allow_pickle=True))
if Model == "Baxter":
    print("Baxter model was created")
    H = model1.H_Baxter(L, Dis_f, J, phi, hc)
if Model == "Baxter_Parafermions":
    print("Baxter Parafermions model was created")
    H =  model1.Baxter_Parafermions(L, J, Dis_f, phi, hc)
if Model == "modf_Baxter":
    print("Modified Baxter was created")
    H = model1.H_modf_Baxter(L, Dis_f, J, phi, hc)
if Model == "ModBaxter_Parafermions":
    print("Modified Baxter Parafermions was created")
    Hmbp = model1.ModBaxter_Parafermions(L, J, Dis_f, phi, hc)

print("Hermitian test of H:", np.matrix.max(abs(H.todense()-np.conjugate(H.todense().T))))

if Model == "Baxter" and Finding == 'Levelstat':
    #plt.rcParams["figure.figsize"] = (10,3)
    #fig, ax = plt.subplots()
    #ax.locator_params(nbins=3, axis='y')
    P0, P1, P2 = model1.FermiParity_sectors(L)
    p = 0
    Ep = []
    Rn = []
    for P in [P0, P1, P2]:
        Hp = model1.x_H_Baxter_sectors(P, L, Dis_f, J, phi, hc)
        Up, Vp = np.linalg.eigh(Hp.todense())
        y = [p for _ in Up]
        #plt.plot(Up,y, '.')
        #np.save("EigenParity_"+str(p)+"_Baxter_L"+str(L)+"_f"+str(f)+".npy", Up)
        p = p + 1
        Ep.append(Up)
        Rn.append(model1.Level_stat(Up[3**(L-2)-1:-(3**(L-2))+1]))
    H = model1.H_Baxter(L, Dis_f, J, phi, hc)
    np.save("Rn_list.npy", Rn)
    Ub, Vb = np.linalg.eigh(H.todense())
    print("Comparing full Baxter model and sectors of Baxter model: ", max(abs(np.array(sorted(Ub))-np.array(sorted(np.array(Ep).flatten())))))
    #plt.savefig("BaxterEnergy"+str(np.round(float(sys.argv[4])))+"_f"+str(f)+".png", dpi=250)
    print("######################################")
    print("Calculation detail:"+"\nL = ", L, ", J= ", J, ", f= ", f, ", phi= ", float(sys.argv[4]),"*(np.pi/6), T=", T, ", hc = ", hc, " and SSH_n = ", SSH_n)
    print("######################################")
    sys.exit("Level statistic was calculated")


    
if hc ==1:
    A_t = list(A_t)
    if Finding == 'Ainfty_ED':
        U, V = np.linalg.eigh(H.todense())
        for nt in time_steps:
            UT = np.dot(V,np.dot(np.diag(np.exp(-1.j*(nt*T)*U)),np.conjugate(V.T)))
            A_t.append(np.trace(np.dot(np.dot(np.dot(np.conjugate(UT.T), np.conjugate(model1.Sigma_z_1(L))),UT),model1.Sigma_z_1(L)))*(1.0/(3.**L)))
            np.save("A_L_"+str(L)+".npy",A_t)
            print("nt=", nt, " is done!,  A_t=", np.round(abs(np.array(A_t)[-1]),3))
        print("AutoCorr function A was created using ED")

        
    if Finding == 'KrylovHmlt_bn':
        #U, V = np.linalg.eigh(H.todense())
        print("Krylov Hamiltonian loop started")
        bn = np.array([]); phin = []; etan = []; O = []; Ortho_nn = []
        if os.path.exists("./bn_L_"+str(L)+".npy") == True:
            bn = np.load("./bn_L_"+str(L)+".npy")
            phin = np.load("./phi_L_"+str(L)+".npy")
            etan = np.load("./eta_L_"+str(L)+".npy")
            O = np.load("./On_L_"+str(L)+".npy", allow_pickle=True)
            Ortho_nn = np.load("Ortho_nn_L_"+str(L)+".npy")
        bns, phis, etas = model1.Krylov_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n, np.array(bn), phin, etan, O, Ortho_nn)
        print(">>>>  bn, phin, etan was calculated!")
        #print("bns: \n", np.round(bns,3))
        bn_Curly_L = np.diag(bns, 1)+np.diag(bns, -1)
        Curly_L = model1.Gram_Schmidt_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n)
        print("bn_Curly_L:\n", np.round(np.real(bn_Curly_L[:10,:10]),2))
        print("GS_Curly_L:\n", np.real(np.round(Curly_L[:10,:10],2)))
        print("GS and recrusive HK: ", max(abs(bn_Curly_L - (Curly_L)).flatten()))
        ### Find hn and htildn
        hn, htildn = model1.h0_htild0(bns)
        np.save("hn_L_"+str(L)+"_"+str(J)+"_"+str(f)+"_"+str(phi)+".npy", hn)
        np.save("htildn_L_"+str(L)+"_"+str(J)+"_"+str(f)+"_"+str(phi)+".npy", htildn)
        print(">>>>  hn and htildn was calculated!")
        ### Finding Kylov Hamiltonina eigenvalues
        print("Finding Kylov Hamiltonina eigenvalues")
        HK = np.diag(bns, -1)+np.diag(bns, 1)
        EK, VK = np.linalg.eig(HK)
        np.save("EK_L"+str(L)+"_"+str(J)+"_"+str(f)+"_"+str(phi)+".npy", EK)
        print(">>>> Kylov Hamiltonina eigenvalues was found!")

        
    if Finding == 'Ainfty_bn':
        print("Krylov Hamiltonian procedure started")
        A_t = model1.AutoCorr_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n, time_steps, 100)
        print("A_in using bn method: ", np.round(A_t, 3))
        np.save("A_bn_L_"+str(L)+".npy", (A_t))
        print("AutoCorr function A was created using Kryolv process")

        
    if Finding == 'GS':
        #fig, ax = plt.subplots()
        #ax.set_title("($O_n|\mathcal{L}|O_m$) matrix elements for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ", hc = "+ str(hc))
        Curly_L = model1.Gram_Schmidt_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n)
        print("Curly_L: \n", np.real(np.round(Curly_L[:10,:10],2)))
        np.save("Curly_L_L"+str(L)+"_hc_"+str(hc)+".npy", Curly_L)
        #sns.heatmap(Curly_L, cmap="mako")
        #plt.savefig("Curly_L_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        Uk, Vk = np.linalg.eigh(Curly_L)
        #fig, ax = plt.subplots()
        #plt.ylim(0.0, 0.1)
        Uclose_0 = sorted(abs(Uk))[:1]
        print(np.round(Uk,3))
        #for s in Uclose_0:
        #    s_index = list(abs(Uk)).index(s)
        #    plt.plot([abs(_)**2 for _ in Vk[:,s_index]], label = "E = "+str(np.round(s,5)))
        #plt.legend()
        #plt.savefig("Vk_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        ##plt.show()
        
print("sigma_z*sigma_z: ", np.trace(np.dot(np.conjugate(model1.Sigma_z_1(L)), model1.Sigma_z_1(L))))


###########################################################################################################
###########################################################################################################
###########################################################################################################


if hc !=1:
    A_t = list(A_t)
    if Finding == 'Ainfty_ED':
        U, V = np.linalg.eig(H.todense())
        Hpr = scipy.sparse.csr_matrix(H - 1.j*max(np.array(U).imag)*np.eye(3**L))
        U, V = np.linalg.eig(Hpr.todense())
        for nt in time_steps:
            UT = np.dot(V,np.dot(np.diag(np.exp(-1.j*(nt*T)*U)),scipy.linalg.inv(V)))
            UTdag = np.dot(np.conjugate(scipy.linalg.inv(V).T), np.dot(np.diag(np.exp(1.j*(nt*T)*np.conjugate(U))),np.conjugate(V.T)))
            A_nom = np.trace(np.dot(np.dot(np.dot(UTdag, np.conjugate(model1.Sigma_z_1(L))),UT),model1.Sigma_z_1(L)))
            A_den = (np.trace(np.dot(UTdag, UT)))
            A_t.append(A_nom/A_den)
            A_t_nom.append(A_nom)
            A_t_den.append(A_den)
            np.save("A_L_"+str(L)+".npy",A_t)
            np.save("A_nom_L_"+str(L)+".npy", A_t_nom)
            np.save("A_den_L_"+str(L)+".npy", A_t_den)
            #print("A_nom: ", A_nom, ", Aden: ", A_den)
            print("nt=", nt, " is done!,  A_t=", np.round(abs(A_t[-1]),3))
        print("AutoCorr function A was created using ED")

        
    if Finding == 'Ainfty_bn':
        print("Krylov Hamiltonian procesure started")
        A_t = model1.AutoCorr_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n, time_steps, 10)
        print("A_in using bn method: ", np.round(A_t, 3))
        np.save("A_bn_L_"+str(L)+".npy", (A_t))
        print("AutoCorr function A was created using Kryolv process")

        
    if Finding == 'GS' or Finding == 'KrylovHmlt_bn':
        #fig, ax = plt.subplots()
        #ax.set_title("($O_n|\mathcal{L}|O_m$) matrix elements for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ",  $\Gamma$ = "+ str(hc))
        Curly_L = model1.Gram_Schmidt_bn(L, H,  sparse.csr_matrix(model1.Sigma_z_1(L)), SSH_n)
        print("Curly_L: \n", np.real(np.round(Curly_L[:13,:13], 2)))
        np.save("Curly_L_L"+str(L)+"_hc_"+str(hc)+".npy", Curly_L)
        #sns.heatmap(Curly_L, cmap="mako")
        #plt.savefig("Curly_L_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        #fig, ax = plt.subplots()
        #ax.set_title("right eigenvector for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ",  $\Gamma$ = "+ str(hc))
        #plt.ylim(0.0, 0.05)
        Uk, Vk = np.linalg.eig(Curly_L)
        #print("Right Uk: ", np.abs(model1.complex_sort(Uk)))
        Uclose_0 = np.argmin(np.abs(Uk))
        #print("Uclose_0: ", Uclose_0, Uk[Uclose_0])
        #plt.plot([abs(_)**2 for _ in Vk[:,Uclose_0]], label = "right E = "+str(Uk[Uclose_0]))
        ##plt.plot([abs(_)**2 for _ in (Vk[:, Uclose_0[0]]+Vk[:, Uclose_0[1]])/np.sqrt(2)], label = "superpos")
        #plt.legend()
        #plt.savefig("Vk_right_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        ######   plottinh right eigenvalues
        #fig, ax = plt.subplots()
        #ax.set_title("Right eigenvalues for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ",  $\Gamma$ = "+ str(hc)) 
        #plt.plot([-10, 10], [0,0], 'k', linewidth=0.5)
        #plt.plot([0, 0], [-10,10], 'k', linewidth=0.5)
        #plt.ylim(min(Uk.imag)-0.5, max(Uk.imag)+0.5)
        #plt.xlim(min(Uk.real)-0.5, max(Uk.real)+0.5)
        #plt.plot(Uk.real, Uk.imag, '.')
        #plt.savefig("Uk_right_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        ######   plt.savefig("Vk_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        #fig, ax = plt.subplots()
        #ax.set_title("Left eigenvector for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ", $\Gamma$ = "+ str(hc))
        #plt.ylim(0.0, 0.05)
        Uk, Vk = np.linalg.eig(np.conjugate(Curly_L.T))
        #print("left Uk: ", np.abs(model1.complex_sort(Uk)))
        Uclose_0 = np.argmin(np.abs(Uk))
        #print("Uclose_0: ", Uclose_0, Uk[Uclose_0])
        #plt.plot([abs(_)**2 for _ in Vk[:,Uclose_0]], label = "left E = "+str(Uk[Uclose_0]))
        ##plt.plot([abs(_)**2 for _ in (Vk[:, Uclose_0[0]]+Vk[:, Uclose_0[1]])/np.sqrt(2)], label = "superpose")
        #plt.legend()
        #plt.savefig("Vk_left_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        #### FFT
        #fig, ax = plt.subplots()
        #ax.set_title("FFT of Left eigenvector for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ", $\Gamma$ = "+ str(hc))
        #plt.ylim(0.0, 0.2)
        y = [abs(_)**2 for _ in Vk[8:30,Uclose_0]]
        yf = fft(y)
        np.save("yFFT.npy", yf)
        np.save("xFFT.npy", np.linspace(0.0, 1, 22, endpoint=False))
        #plt.plot([0.5,0.5], [0, 0.5], '--k')
        #plt.plot(np.linspace(0.0, 1, 22, endpoint=False)[4:-4], yf[4:-4])
        #plt.savefig("FFT_left_L"+str(L)+"_hc_"+str(hc)+"_f"+str(f)+".png", dpi=250)
        ######   plottinh left eigenvalues
        #fig, ax = plt.subplots()
        #ax.set_title("Left eigenvalues for \nL = "+str(L)+", J= "+ str(J)+ ", f= "+ str(f)+ ", phi= "+ str(float(sys.argv[4]))+ "*$\pi/6$, T="+ str(T)+ ",  $\Gamma$ = "+ str(hc)) 
        #plt.plot([-10, 10], [0,0], 'k', linewidth=0.5)
        #plt.plot([0, 0], [-10,10], 'k', linewidth=0.5)
        #plt.ylim(min(Uk.imag)-0.5, max(Uk.imag)+0.5)
        #plt.xlim(min(Uk.real)-0.5, max(Uk.real)+0.5)
        #plt.plot(Uk.real, Uk.imag, '.')
        #plt.savefig("Uk_left_L"+str(L)+"_hc_"+str(hc)+".png", dpi=250)
        ##plt.show()

        
    if Finding == "Psi_left":
        w = np.e**(1.j*2*np.pi/3.0)
        Hpfr = model1.modf_Baxter_Pfermion(L, J, f, phi, hc)
        Psi_L = model1.Psi_Left(L, f, J, phi)
        #print("Sz_1:\n", np.round(model1.Sigma_z_1(L),2))
        #print("Psi_L:\n", np.round(Psi_L.A,2))
        print("Sz_1 - Psi_L = ", max(abs(model1.Sigma_z_1(L)-Psi_L.A).flatten()))
        print("[H, Psi_Left] =", max(abs(np.array((H*Psi_L-Psi_L*H).A).flatten())))
        CommuteTheory = (f*(w-1)*model1.alpha_j(L, 2*L-1)*(((f*w)/(J*np.exp(1.j*phi)))**(L-1)))
        print("[H2, Psi_Left] =", max(abs(np.array((model1.alpha_j(L, (2*L)-1)*(Hpfr*Psi_L-Psi_L*Hpfr) - (Hpfr*Psi_L-Psi_L*Hpfr)*model1.alpha_j(L, (2*L)-1) ).A).flatten())))

        
print("######################################")
print("Calculation detail:"+"\nL = ", L, ", J= ", J, ", f= ", f, ", phi= ", float(sys.argv[4]),"*(np.pi/6), T=", T, ", hc = ", hc, " and SSH_n = ", SSH_n)
print("######################################")
