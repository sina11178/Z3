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

def time_evol_op_at_time_t(hamiltonian, time, gamma):
    if gamma == 1:
        A, U = np.linalg.eigh(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
        U_T = np.conjugate(U.T)  # Conjugate transpose of U
        #A_matrix = exp_diag_matrix(A, time)
        A_matrix = np.diag(np.exp(A* time*(-1.j)))
        evol_to_t = U @ A_matrix @ U_T
        return evol_to_t
    else:
        A, U = np.linalg.eig(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
        U_inv = np.linalg.inv(U)  # Conjugate transpose of U
        #A_matrix = exp_diag_matrix(A, time)
        A_matrix = np.diag(np.exp((A - (1.j)*np.imag(A).max())* time*(-1.j)))
        evol_to_t = U @ A_matrix @ U_inv
        return evol_to_t

def autocorrelator(time, size, hamiltonian, gamma):  # ONLY FOR SIG1_X
    list = []
    for a in range(size):
        list.append(np_I)
    list[0] = np_SIG # first site association with sig_x set
    #list[size-1] = np_SIG
    
    # print(time_evol_op.getH() * time_evol_op)
    if gamma == 1:
        time_evol_op = time_evol_op_at_time_t(hamiltonian, time, gamma)  # Finds time evol operator at time t
        non_time_evold_op_sig = tensor_product(list, size) #  SIGX_1(0)
        time_evol_op_dag = (np.conjugate(time_evol_op.T))
        time_evolved_sig = time_evol_op_dag @ non_time_evold_op_sig @ (time_evol_op)  #  SIGX_1(time)
        #print(time_evolved_sig)
        prefactor = 1/(3**size)
        product_SIG = np.conjugate(time_evolved_sig.T) @ non_time_evold_op_sig
        autocorr_SIG = prefactor * np.trace(product_SIG)
        print(autocorr_SIG)

        return autocorr_SIG
    else:
        time_evol_op = time_evol_op_at_time_t(hamiltonian, time, gamma)  # Finds time evol operator at time t
        #norm = np.linalg.norm(time_evol_op)
        #norm_time_evol_op = (1/norm)*time_evol_op
        non_time_evold_op_sig = tensor_product(list, size) #  SIGX_1(0)
        #norm_time_evol_op_dag = (np.conjugate(norm_time_evol_op.T))
        #time_evolved_sig = norm_time_evol_op_dag @ non_time_evold_op_sig @ (norm_time_evol_op)  #  SIGX_1(time)
        time_evol_op_dag = (np.conjugate(time_evol_op.T))
        time_evolved_sig = time_evol_op_dag @ non_time_evold_op_sig @ (time_evol_op)
        #print(time_evolved_sig)
        prefactor = 1/(np.trace(time_evol_op_dag @ time_evol_op))
        product_SIG = np.conjugate(time_evolved_sig.T) @ non_time_evold_op_sig
        #product_SIG = (time_evolved_sig) @ non_time_evold_op_sig

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

def output_list(input_list, size, hamiltonian, gamma):# TOOK OUT .REAL
    out_list = [0] * len(input_list)
    for i in range(len(input_list)):
        out_list[i] = autocorrelator(input_list[i], size, hamiltonian, gamma).real  #  Stores respective outputs in list
        #out_list[i] = autocorrelator(input_list[i], size, hamiltonian) --> DOES NOT FORCE REAL...
    #print(out_list)
    return np.array(out_list)


'''

RECURSION METHOD FOR HERMITIAN HAMILTONIAN

'''


def matrix_2_vector(matrix):
    matrix_transpose_1 = matrix.T
    matrix_transpose = matrix_transpose_1.tolist()
    vector = matrix_transpose[0]
    for i in range(1, len(matrix_transpose)):
        vector = vector + matrix_transpose[i]
    vector_return = np.array(vector)
    #print(vector_return)
    return (vector_return)

def Liovillan(Hamiltonian, size): #  Liovillan only applied to a matrix converted to vector
    I = np.eye(3**size)
    S_1 = np.kron(I, Hamiltonian)
    S_2 = np.kron((Hamiltonian.T), I)
    #S_1 = np.kron(Hamiltonian, I)
    #S_2 = np.kron(I, Hamiltonian)
    #print(S_1)

    return (S_1+(-1*S_2))

def A_vectors(liovillan, b_array, O_basis_array, n_index):  # Used for A_2 and on
    #first_term = L @ O_basis_array[n_index]
    first_term = liovillan @ O_basis_array[n_index]
    second_term = (b_array[n_index - 1]) * O_basis_array[(n_index - 1)]

    return (first_term +(-1*second_term))


def norm (L_1, L_2, size):  # L_1 is ket, L_2 is the bra
    ket = L_1
    bra = np.conjugate(L_2)

    norm_squared = np.dot(bra, ket)
    norm = cm.sqrt((1/(3**size))*norm_squared)
    print(norm)
    #return np.linalg.norm(norm)
    return norm.real

def recursion_method(n, H, size):
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_SIG
    m = np.sqrt((1/(3**size)))

    strong_mode = tensor_product(list, size)
    L = (Liovillan(H, size)).tolist()
    b = []  # ORDER: 1, 2, ...
    O_0 = ((strong_mode.T).flatten()).tolist() # Returns as list w/ vector inside of it;  ORDER: 0, 1, 2, ..., n

    #print(O_1)
    krylov_basis = [np.array(O_0[0])]  # YOU CAN ALSO JUST NORMALIZE THIS AND LEAVE OUT THE PREFACTOR IN THE NORM FUNCTION
    A_basis = []  # ORDER: 1, 2, ...
    for i in range(0, n):
        if (i == 0):  # BASE CASE
            A_1 = np.array(L) @ krylov_basis[0]
            #A_1 = krylov_basis[0] @ np.array(L)
            A_basis.append(A_1)
            #b_1 = norm(A_1, A_1, size)
            b_1 = m * np.linalg.norm(A_1)
            b.append(b_1)
            #print(type(A_1))
            O_1 = (1/(b_1))*A_1
            krylov_basis.append(O_1)

        else:
            k = A_vectors(L, b, krylov_basis, i)
            A_basis.append(k)
            #b.append(norm(k, k, size))
            b.append(m * np.linalg.norm(k))
            krylov_basis.append((1/b[i])*A_basis[i])

    #print(b)


    matrix_final = np.eye(len(krylov_basis))
    for j in range(len(krylov_basis)):
        k = L @ krylov_basis[j]
        for i in range(len(krylov_basis)):
            matrix_final[i, j] = np.conjugate(krylov_basis[i]) @ k
    return matrix_final
    #return b


'''

GRAHM SCHMIDT - FOR NON HERMITIAN HAMILTONIANS

'''



def grahm_schmidt(n, H, size):
    L = np.array((Liovillan(H, size)).tolist())
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_SIG
    strong_mode = tensor_product(list, size)
    O_0 = ((strong_mode.T).flatten()).tolist()
    gs_basis_non_ortho = [np.array(O_0[0])]

    # CREATION OF VECTORS I WANT TO ORTHOGONALIZE
    U = []  # Creation of list that will hold orthogonal vectors
    for j in range(1, n):
        gs_basis_non_ortho.append(L @ gs_basis_non_ortho[j-1])

    for i in range(len(gs_basis_non_ortho)):
        v = gs_basis_non_ortho[i]
        for k in range(i):
            proj = ((np.dot(np.conjugate(v), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
            v -= proj
        U.append(v / np.linalg.norm(v))

    # COLLECT B VALUES
    matrix_final = np.eye(len(U))
    for j in range(len(U)):
        k = L @ U[j]
        for i in range(len(U)):
            matrix_final[i, j] = np.conjugate(U[i]) @ k

    return matrix_final

    # b = []  # B VALUES
    # for p in range(0, len(U) - 1):
    #     k = L @ U[p+1]
    #     b.append(np.conjugate(U[p]) @ k)
    # print(b)
    # return b