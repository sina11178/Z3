'''
PYTHON CODE FOR Z3 PROJECT

REFERENCE FOR HAMILTONIAN: https://arxiv.org/pdf/1407.6376

Z3 symmetry --> N = 3
'''
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
#import scipy

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 1
phi = 1
f = 1
L = 1  # SYSTEM SIZE
gamma = 1  # DETERMINES HERMICITY

w = np.exp((2 * np.pi * 1.j)/3)  # For defining Sigma operator

sigma = [
    [1, 0, 0],
    [0, w, 0],
    [0, 0, w**2]
]

clock = [  # TAU SQUARED
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

#np_SIG = np.array(sigma)
#np_CLK = np.array(clock)
#np_I = np.array(I)

np_SIG = sci.csr_matrix(sigma, dtype=np.complex128)  # Sigma operator
np_CLK = sci.csr_matrix(clock, dtype=np.complex128)  # Clock operator
np_I = sci.csr_matrix(I, dtype=np.complex128)  # Identity operator


'''

HAMILTONIAN CONSTRUCTION

'''

'''
def tensor_product(array, length):  # Calculates tensor product of all matrices needed

    #final_matrix = np.zeros((3, 3), dtype=complex)
    for i in range(length):
        if i == 0:
            final_matrix = array[i]
            #print(final_matrix)
        else:
            final_matrix = np.kron(final_matrix, array[i])
    return final_matrix
'''

def tensor_product(array, length):  # Calculates tensor product of all matrices needed

    #final_matrix = np.zeros((3, 3), dtype=complex)
    for i in range(length):
        if i == 0:
            final_matrix = array[i]
            #print(final_matrix)
        else:
            # final_matrix = np.kron(final_matrix, array[i])
            final_matrix = sci.kron(final_matrix, array[i])
    return final_matrix

#print(-1*tensor_product([np.conjugate(np_SIG.T), np_SIG], 2))    

'''
def adding_matrixes_in_list(array, length):  #  Used in all terms
    matrix = np.zeros((3**length, 3**length), dtype=complex)
    
    for i in range(len(array)):
        #print(len(array))
        matrix = matrix + array[i]

    return matrix
'''

def adding_matrixes_in_list(array, length):  #  Used in all terms
    # matrix = np.zeros((3**length, 3**length), dtype=np.complex128)
    matrix = sci.csr_matrix((3**length, 3**length), dtype=np.complex128)

    for i in range(len(array)):
        #print(len(array))
        matrix = matrix + array[i]

    return matrix



def term_one(length, Jay, PHI):

    #list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)
    list_2 = []  # Has all terms in a list

    #SIG_dag = np.conjugate(np_SIG.T)
    SIG_dag = np_SIG.conjugate().T
    # LIST OF DIFFERENT MATRICES ASSOCIATED WITH SITES CREATED HERE
    for i in range(length - 1):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators

        list[i] = SIG_dag
        # if i == (length - 1):
        #     list[0] = np_x
        # else:
        list[i+1] = np_SIG

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_XX = -Jay * np.exp(PHI * 1.j) * adding_matrixes_in_list(list_2, length)

    return H_XX

def term_two(F, length):

    list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)

    list_2 = []  # Has all terms in a list

    for i in range(length):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators
        list[i] = np_CLK  # Clock operator applied to Sigma operator
        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
    #print(list_2)
    CLOCK_SUM = -F * adding_matrixes_in_list(list_2, length)  # NOTE: SAEED DOES NOT HAVE A w PREFACTOR FOR THIS TERM

    #print('Clock Sum:')
    #print(CLOCK_SUM)
    return CLOCK_SUM


'''

RUNNING SIMULATIONS

'''


#  AUTOCORRELATOR

def time_evol_op_at_time_t(hamiltonian, time, gamma, size):
    if gamma == 1:
        A, U = np.linalg.eigh(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
        U_T = np.conjugate(U.T)  # Conjugate transpose of U
        #A_matrix = exp_diag_matrix(A, time)
        A_matrix = np.diag(np.exp(A* time*(-1.j)))
        evol_to_t = U @ A_matrix @ U_T
        return evol_to_t
    else:
        # A, U = np.linalg.eig(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
        # U_inv = np.linalg.inv(U)  # Conjugate transpose of U
        U, V = np.linalg.eig(hamiltonian)
        #print(hamiltonian.shape)
        #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
        Hpr = (hamiltonian - 1.j*max(np.array(U).imag)*np.eye(3**size))

        A, U = np.linalg.eig(Hpr)
        U_inv = np.linalg.inv(U)
        #A_matrix = exp_diag_matrix(A, time)
        #A_matrix = np.diag(np.exp((A - (1.j)*np.imag(A).max())* time*(-1.j)))
        A_matrix = np.diag(np.exp(A*(-1.j)*time))
        evol_to_t = U @ A_matrix @ U_inv
        print('SUBTRACTION: \n')
        print(np.max(np.abs(np.conjugate(evol_to_t.T) - (np.conjugate(U_inv.T) @ np.conjugate(A_matrix.T) @ np.conjugate(U.T)))))
        return evol_to_t

def autocorrelator(time, size, hamiltonian, gamma):  # ONLY FOR SIG1_X
    list = []
    for a in range(size):
        list.append(np_I)
    list[0] = np_SIG # first site association with sig_x set
    #list[size-1] = np_SIG
    
    # print(time_evol_op.getH() * time_evol_op)
    if gamma == 1:
        time_evol_op = time_evol_op_at_time_t(hamiltonian, time, gamma, size)  # Finds time evol operator at time t
        non_time_evold_op_sig = tensor_product(list, size) #  SIGX_1(0)
        time_evol_op_dag = (np.conjugate(time_evol_op.T))
        time_evolved_sig = time_evol_op_dag @ non_time_evold_op_sig @ (time_evol_op)  #  SIGX_1(time)
        #print(time_evolved_sig)
        prefactor = 1/(3**size)
        product_SIG = np.conjugate(time_evolved_sig.T) @ non_time_evold_op_sig
        autocorr_SIG = prefactor * np.trace(product_SIG)
        #print(autocorr_SIG)

        return autocorr_SIG
    else:
        time_evol_op = time_evol_op_at_time_t(hamiltonian, time, gamma, size)  # Finds time evol operator at time t
        #norm = np.linalg.norm(time_evol_op)
        #norm_time_evol_op = (1/norm)*time_evol_op
        non_time_evold_op_sig = tensor_product(list, size) #  SIGX_1(0)
        #norm_time_evol_op_dag = (np.conjugate(norm_time_evol_op.T))
        #time_evolved_sig = norm_time_evol_op_dag @ non_time_evold_op_sig @ (norm_time_evol_op)  #  SIGX_1(time)
        time_evol_op_dag = (np.conjugate(time_evol_op.T))  # ERROR THIS LINE
        time_evolved_sig = time_evol_op_dag @ np.conjugate(non_time_evold_op_sig) @ (time_evol_op)  # TAKING DAGGER OF SIGMA ONLY (RATHER THAN ENTIRE THING)
        #print(time_evolved_sig)
        prefactor = 1/(np.trace(time_evol_op_dag @ time_evol_op))
        #product_SIG = np.conjugate(time_evolved_sig.T) @ non_time_evold_op_sig
        product_SIG = (time_evolved_sig) @ non_time_evold_op_sig
        #product_SIG = (time_evolved_sig) @ non_time_evold_op_sig

        autocorr_SIG = prefactor * np.trace(product_SIG)
        #print(autocorr_SIG)

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


def Liovillan(Hamiltonian, size): #  Liovillan only applied to a matrix converted to vector
    I = sci.eye(3**size, format='csr')
    Hamiltonian = sci.csr_matrix(Hamiltonian)
    S_1 = sci.kron(I, Hamiltonian)
    S_2 = sci.kron((Hamiltonian.T), I)
    #S_1 = np.kron(Hamiltonian, I)
    #S_2 = np.kron(I, Hamiltonian)
    #print(S_1)

    return (S_1+(-1*S_2))

def A_vectors(liovillan, b_array, O_basis_array, n_index, ham, size):  # Used for A_2 and on
    #first_term = (((ham @ np.reshape(O_basis_array[n_index], (3**size, 3**size)).T) - 
    #               (np.reshape(O_basis_array[n_index], (3**size, 3**size)).T @ ham)).T).flatten()
    first_term = liovillan @ O_basis_array[n_index]
    second_term = (b_array[n_index - 1]) * O_basis_array[(n_index - 1)]

    return (first_term +(-1*second_term))


def norm (L_1, L_2, size):  # L_1 is ket, L_2 is the bra
    ket = L_1
    bra = np.conjugate(L_2)

    norm_squared = np.dot(bra, ket)
    norm = cm.sqrt((1/(3**size))*norm_squared)
    #print(norm)
    #return np.linalg.norm(norm)
    return norm.real

def recursion_method(n, H, size):
    list = []
    list_squared = []
    for i in range(size):
        list.append(np_I)
        list_squared.append(np_I)
    list[0] = np_SIG
    #list_squared = list
    list_squared[0] = np_SIG.conjugate().T  # FOR ALPHA^2 OPERATOR
    m = np.sqrt((1/(3**size)))

    seeding = tensor_product(list, size).toarray()
    seeding_squared = tensor_product(list_squared, size).toarray()  # FOR ALPHA^2 OPERATOR
    Lio = Liovillan(H, size)
    b = []  # ORDER: 1, 2, ...
    b_squared = []  # FOR ALPHA^2 OPERATOR
    O_0 = ((seeding.T).flatten()) # Returns as list w/ vector inside of it;  ORDER: 0, 1, 2, ..., n
    O_0_squared = ((seeding_squared.T).flatten())  # FOR ALPHA^2 OPERATOR

    #print(O_1)
    krylov_basis = [O_0]  # YOU CAN ALSO JUST NORMALIZE THIS AND LEAVE OUT THE PREFACTOR IN THE NORM FUNCTION
    krylov_basis_squared = [O_0_squared]  # FOR ALPHA^2 OPERATOR
    A_basis = []  # ORDER: 1, 2, ...
    A_basis_squared = []  # FOR ALPHA^2 OPERATOR
    for i in range(0, n):
        if (i == 0):  # BASE CASE
            A_1 = Lio @ krylov_basis[0]
            #A_1 = (((H @ np.reshape(krylov_basis[0], (3**size, 3**size)).T) - 
             #      (np.reshape(krylov_basis[0], (3**size, 3**size)).T @ H)).T).flatten()
            #A_1 = krylov_basis[0] @ np.array(L)
            A_basis.append(A_1)
            #b_1 = norm(A_1, A_1, size)
            b_1 = m * np.linalg.norm(A_1)
            b.append(b_1)
            #print(type(A_1))
            O_1 = (1/(b_1))*A_1
            krylov_basis.append(O_1)

            # NOTE: FOR ALPHA^2 OPERATOR
            A_1_squared = Lio @ krylov_basis_squared[0]
            #A_1_squared = (((H @ np.reshape(krylov_basis_squared[0], (3**size, 3**size)).T) -
            #               (np.reshape(krylov_basis_squared[0], (3**size, 3**size)).T @ H)).T).flatten()
            #A_1_squared = krylov_basis_squared[0] @ np.array(L)
            A_basis_squared.append(A_1_squared)
            b_1_squared = m * np.linalg.norm(A_1_squared)  # FOR ALPHA^2 OPERATOR
            b_squared.append(b_1_squared)
            O_1_squared = (1/b_1_squared) * A_1_squared  # FOR ALPHA^2 OPERATOR
            krylov_basis_squared.append(O_1_squared)  # FOR ALPHA^2 OP

        else:
            k = A_vectors(Lio, b, krylov_basis, i, H, size)
            A_basis.append(k)
            #b.append(norm(k, k, size))
            b.append(m * np.linalg.norm(k))
            krylov_basis.append((1/b[i])*A_basis[i])

            # NOTE: FOR ALPHA^2 OPERATOR

            k_squared = A_vectors(Lio, b_squared, krylov_basis_squared, i, H, size)
            A_basis_squared.append(k_squared)
            #b_squared.append(norm(k_squared, k_squared, size))  # FOR ALPHA^2 OPERATOR
            b_squared.append(m * np.linalg.norm(k_squared))  # FOR ALPHA^2
            krylov_basis_squared.append((1/b_squared[i]) * A_basis_squared[i])  # FOR ALPHA^2 OPERATOR

    #print(b)

    l = len(krylov_basis)
    matrix_final = np.zeros((l,l), dtype=complex)
    matrix_final_squared = np.zeros((l,l), dtype=complex)  # FOR ALPHA^2 OPERATOR
    for j in range(len(krylov_basis)):
        k = Lio @ krylov_basis[j]
        k_squared = Lio @ krylov_basis_squared[j]  # FOR ALPHA^2 OPERATOR
        for i in range(len(krylov_basis)):
            matrix_final[i, j] = np.conjugate(krylov_basis[i].T) @ k
            matrix_final_squared[i, j] = np.conjugate(krylov_basis_squared[i].T) @ k_squared
    #return matrix_final
    return [matrix_final, matrix_final_squared]  # FOR ALPHA^2 OPERATOR
    #return b


'''

GRAHM SCHMIDT - FOR NON HERMITIAN HAMILTONIANS

'''



def grahm_schmidt(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    list_squared = []  # FOR ALPHA^2 OPERATOR
    for i in range(size):
        list.append(np_I)
        list_squared.append(np_I)  # FOR ALPHA^2 OPERATOR
    list[0] = np_SIG
    #list[0] = np.conjugate(np_SIG).T # FOR ALPHA^2 OPERATOR
    #strong_mode = alph_3  # tensor_product(list, size)
    #list[0] = np.conjugate(np_SIG.T)  # FOR ALPHA^2 OPERATOR
    list_squared[0] = np_SIG.conjugate().T  # FOR ALPHA^2 OPERATOR
    seeding = tensor_product(list, size).toarray()
    seeding_squared = tensor_product(list_squared, size).toarray()  # FOR ALPHA^2 OPERATOR
    #strong_mode = tensor_product(list, size) + (f*w/(J*np.exp(1.j*phi))) * np.conjugate(alph_3).T
    O_0 = ((seeding.T).flatten())
    gs_basis_non_ortho = [O_0]
    #print(list == list_squared)

    O_0_squared = ((seeding_squared.T).flatten())  # FOR ALPHA^2 OPERATOR
    gs_basis_non_ortho_squared = [O_0_squared]  # FOR ALPHA^2 OPERATOR

    for P in range(1, n+1):
        K = Lio @ gs_basis_non_ortho[P-1]
        K_squared = Lio @ gs_basis_non_ortho_squared[P-1]
        #f = (((H @ np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T @ H)).T).flatten()
        gs_basis_non_ortho.append(K)
        gs_basis_non_ortho_squared.append(K_squared)

    gs_basis_non_ortho_matrix = np.column_stack(gs_basis_non_ortho)
    gs_basis_non_ortho_matrix_squared = np.column_stack(gs_basis_non_ortho_squared)  # FOR ALPHA^2 OPERATOR
    Q, R = np.linalg.qr(gs_basis_non_ortho_matrix)
    Q_squared, R_squared = np.linalg.qr(gs_basis_non_ortho_matrix_squared)  # FOR ALPHA^2 OPERATOR
    #print(Q == Q_squared)

    l = len(gs_basis_non_ortho)
    matrix_final = np.zeros((l,l), dtype=complex)  # FOR ALPHA OPERATOR
    matrix_final_squared = np.zeros((l,l), dtype=complex)  # FOR ALPHA^2 OPERATOR
    for j in range(l):
        c = Lio @ Q[:, j]
        c_squared = Lio @ Q_squared[:, j]
        for i in range(l):
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ c
            matrix_final_squared[i, j] = np.conjugate(Q_squared[:, i]) @ c_squared
    return [matrix_final, matrix_final_squared]
    # return matrix_final




#print(tensor_product([np.array([[0,0], [0,0]]), np.array([[0,0], [0,0]])], 2))
#print(type(tensor_product([np.array([[0,0], [0,0]]), np.array([[0,0], [0,0]])], 2)))
