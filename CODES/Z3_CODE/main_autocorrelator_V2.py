import matplotlib.pyplot as plt
#import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
from scipy.sparse.linalg import eigsh, eigs
import numpy as np
import cmath as cm
import time
#import scipy

def Ham(J, phi, Len, f, gamma):
    w = np.exp(2j * np.pi / 3)
    D = 3**Len
    powers = 3**np.arange(Len)  # precompute 3^j for each site

    rows, cols, data = [], [], []

    for s in range(D):
        # Extract all site digits at once
        digits = (s // powers) % 3  # shape (L,)

        # --- Diagonal: -J e^{i phi} sum_j omega^{-digits[j]} * omega^{digits[j+1]} ---
        diag = 0.0 + 0.0j
        for j in range(Len - 1):
            diag += -J * np.exp(1.j * phi) * (w**(-digits[j])) * (w**(digits[j+1]))

        rows.append(s); cols.append(s); data.append(diag)

        # --- Off-diagonal: -f tau†_j for each site ---
        # tau† sends digit d -> (d+1)%3, changing state s
        for j in range(Len):
            old_d = digits[j]
            new_d = (old_d + 1) % 3
            # New basis state: replace digit j
            s2 = s + (new_d - old_d) * powers[j]
            rows.append(s2); cols.append(s); data.append(-f)

    H = sci.coo_matrix((data, (rows, cols)), shape=(D, D), dtype=np.complex128).tocsr()
    return H + gamma * H.getH()  # Add Hermitian conjugate to make it non-Hermitian if gamma != 0


def make_sigma_site0(length):
    w = np.exp(2j * np.pi / 3)
    D = 3**length
    powers = 3**np.arange(length)
    
    data = np.zeros(D, dtype=np.complex128)
    
    for s in range(D):
        digit_0 = (s // powers[0]) % 3   # = s % 3, just the site-0 digit
        data[s] = w**digit_0
    
    return sci.diags(data, format='csr')



# NOTE: NON-HERMITIAN HAMILTONIAN
def autocorr_nonherm_all_times(A, U, U_inv, times, len):
    results = []
    sigma = make_sigma_site0(len).toarray()
    sigma_conj = np.conj(sigma)
    A_conj = np.conj(A)
    U_inv_conj = np.conj(U_inv)
    U_conj = np.conj(U)
    for t in times:
        exp_conj_A = np.diag(np.exp(1.j * t * A_conj))
        exp_A = np.diag(np.exp(-1.j * t * A))
        #numerator = np.trace(U_inv_dag @ exp_conj_A @ U_dag @ sigma_conj @ U @ exp_A @ U_inv @ sigma)
        #denominator = np.trace(U_inv_dag @ exp_conj_A @ U_dag @ U @ exp_A @ U_inv)
        result = np.einsum('ai,ka,bk,mb,mn,nl,lj,ji->', U_inv_conj, exp_conj_A, U_conj, sigma_conj, U, exp_A, U_inv, sigma, optimize='optimal')
        denominator = np.einsum('ai,ka,bk,bn,nl,li->', U_inv_conj, exp_conj_A, U_conj, U, exp_A, U_inv, optimize='optimal')

        results.append(result / denominator)
    return results


# NOTE: HERMITIAN HAMILTONIAN
def autocorr_herm_all_times(A, U, times, len):
    results = []
    sigma = make_sigma_site0(len).toarray()
    sigma_conj = np.conj(sigma)
    U_conj = np.conj(U)
    for t in times:
        exp_conj_A = np.diag(np.exp(1.j * t * A))
        exp_A = np.diag(np.exp(-1.j * t * A))
        result = np.einsum('ia,ak,bk,mb,mn,nl,lj,ji->', U, exp_conj_A, U_conj, sigma_conj, U, exp_A, U_conj, sigma, optimize='optimal')
        denominator = 3**len

        results.append(result / denominator)
    return results





pd = {}
pd["size"] = 6
pd["J"] = 1
pd["phi"] = np.pi/6
pd["f"] = 0.1
pd["gamma"] = 1
pd["H"] = Ham(pd["J"], pd["phi"], pd["size"], pd["f"], pd["gamma"])
pd["times_log"] = np.arange(-10, 25, 0.5)
pd["times"] = np.exp(pd["times_log"])



def main(pd):
    H = pd["H"]


    # NOTE: non-hermitian
    if pd["gamma"] != 1:
        H_dense = H.toarray()
        A, U = np.linalg.eig(H_dense)      # non-Hermitian: eig not eigh
        U_inv = np.linalg.inv(U)
        #Hpr = (H - 1.j*max(np.array(U).imag)*np.eye(3**pd["size"]))
        A = A - np.max(A.imag)*1.j

        results = autocorr_nonherm_all_times(A, U, U_inv, pd["times"], pd["size"])
        npy_name = "autocorr_gamma_" + str(pd["gamma"]) + "_L" + str(pd["size"]) +"_f_"+str(round(pd["f"],3))+"_phi_"+str(round(pd["phi"],3))+"_NPY.npy" #NOTE: ONLY FOR ARRAY SIZE OF 1
        all_results = [np.exp(pd["times_log"]), results]
        np.save(npy_name, np.array(all_results))
        print('SUCCESS')
    

    #NOTE: hermitian
    if pd["gamma"] == 1:
        H_dense = H.toarray()
        A, U = np.linalg.eigh(H_dense)      # non-Hermitian: eig not eigh
        U_inv = np.linalg.inv(U)
        #Hpr = (H - 1.j*max(np.array(U).imag)*np.eye(3**pd["size"]))
        #A = A - np.max(A.imag)*1.j

        results = autocorr_nonherm_all_times(A, U, U_inv, pd["times"], pd["size"])
        npy_name = "autocorr_gamma_" + str(pd["gamma"]) + "_L" + str(pd["size"]) +"_f_"+str(round(pd["f"],3))+"_phi_"+str(round(pd["phi"],3))+"_NPY.npy" #NOTE: ONLY FOR ARRAY SIZE OF 1
        all_results = [np.exp(pd["times_log"]), results]
        np.save(npy_name, np.array(all_results))
        print('SUCCESS')

#start_time = time.perf_counter()
main(pd)