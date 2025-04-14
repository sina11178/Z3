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
np_I = np.matrix(I)

'''

HAMILTONIAN CONSTRUCTION

'''

def tensor_product(array, length):  # Calculates tensor product of all matrices needed

    final_matrix = np.zeros((3**length, 3**length))
    for i in range(length):
        if i == 0:
            final_matrix = array[i]
            #print(final_matrix)
        else:
            final_matrix = np.kron(final_matrix, array[i])
    return final_matrix


def adding_matrixes_in_list(array, length):  #  Used in all terms
    matrix = np.zeros((3**length, 3**length))
    
    for i in range(len(array)):
        #print(len(array))
        matrix = matrix + array[i]

    return matrix

def term_one(length, Jay, PHI):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)
    list_2 = []  # Has all terms in a list

    SIG_dag = np.conjugate(np_SIG.T)
    # LIST OF DIFFERENT MATRICES ASSOCIATED WITH SITES CREATED HERE
    for i in range(length - 1):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I

        list[i] = SIG_dag
        # if i == (length - 1):
        #     list[0] = np_x
        # else:
        list[i+1] = np_SIG

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_XX = -Jay * np.exp(PHI* 1.j) * adding_matrixes_in_list(list_2, length)

    return H_XX

def term_two(F, length):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)

    list_2 = []  # Has all terms in a list

    for i in range(length):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I
        list[i] = np_CLK
        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
    #print(list_2)
    CLOCK_SUM = -F * adding_matrixes_in_list(list_2, length)

    return CLOCK_SUM


'''

RUNNING SIMULATIONS

'''


#  AUTOCORRELATOR

def time_evol_op_at_time_t(hamiltonian, time):
    A, U = np.linalg.eigh(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
    U_T = np.conjugate(U.T)  # Conjugate transpose of U
    #A_matrix = exp_diag_matrix(A, time)
    A_matrix = np.diag(np.exp(A* time*(-1.j)))
    evol_to_t = U @ A_matrix @ U_T
    return evol_to_t

def autocorrelator(time, size, hamiltonian):  # ONLY FOR SIG1_X
    list = []
    for a in range(size):
        list.append(np_I)


    list[0] = np_SIG # first site association with sig_x set
    #list[size-1] = np_SIG
    time_evol_op = time_evol_op_at_time_t(hamiltonian, time)  # Finds time evol operator at time t
    non_time_evold_op_sig = tensor_product(list, size) #  SIGX_1(0)
    time_evolved_sig = (np.conjugate(time_evol_op.T)) @ non_time_evold_op_sig @ (time_evol_op)  #  SIGX_1(time)

    # print(time_evol_op.getH() * time_evol_op)
    prefactor = 1/(3**size)
    product_SIG = np.conjugate(time_evolved_sig.T) @ non_time_evold_op_sig
    autocorr_SIG = prefactor * np.trace(product_SIG)
    print(autocorr_SIG)

    return autocorr_SIG

def list_inputs(init_time, fin_time, interval):
    list_non_exp = np.arange(init_time, fin_time, interval)
    list_exp = np.exp(list_non_exp)
    #num_of_runs = int((fin_time-init_time)/interval)
    #for i in range(num_of_runs + 1):
    #    list.append(init_time + i*interval)

    return [list_non_exp, list_exp] # SECOND INDEX HAS LIST OF ACTUAL TIMES; FIRST HAS LOG TIMES

def output_list(input_list, size, hamiltonian):# TOOK OUT .REAL
    out_list = [0] * len(input_list)
    for i in range(len(input_list)):
        out_list[i] = autocorrelator(input_list[i], size, hamiltonian).real  #  Stores respective outputs in list
        #out_list[i] = autocorrelator(input_list[i], size, hamiltonian) --> DOES NOT FORCE REAL...
    #print(out_list)
    return np.array(out_list)

