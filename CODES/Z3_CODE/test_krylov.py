import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
from scipy.sparse.linalg import eigs
#import scipy

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 0
phi = 0  # PHASE FACTOR
f = -1
L = 2 # SYSTEM SIZE
#gamma = 1  # DETERMINES HERMICITY

w = np.exp((2 * np.pi * 1.j)/3)  # For defining Sigma operator

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

np_SIG = np.array(sigma, dtype=np.complex128)  # Sigma operator
np_CLK = np.array(clock, dtype=np.complex128)  # Clock operator
np_I = np.array(I, dtype=np.complex128)  # Identity operator

'''

HAMILTONIAN CONSTRUCTION

'''

# [I, tau, I, I] --> I (x) tau (x) I (x) I

def tensor_product(array, length):  # Calculates tensor product of all matrices needed

    #final_matrix = np.zeros((3, 3), dtype=complex)
    for i in range(length):
        if i == 0:
            final_matrix = array[i]
            #print(final_matrix)
        else:
            final_matrix = np.kron(final_matrix, array[i])
    return final_matrix


#print(-1*tensor_product([np.conjugate(np_SIG.T), np_SIG], 2))    


def adding_matrixes_in_list(array, length):  #  Used in all terms
    matrix = np.zeros((3**length, 3**length), dtype=np.complex128)
    
    for i in range(len(array)):
        #print(len(array))
        matrix = matrix + array[i]

    return matrix


# NOTE: THE FOLLOWING IS THE HAMILTONIAN IN SIGMA BASIS


def term_one_SIGMA(length, Jay, PHI):

    #list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)
    list_2 = []  # Has all terms in a list

    SIG_dag = np.conjugate(np_SIG.T)
    # LIST OF DIFFERENT MATRICES ASSOCIATED WITH SITES CREATED HERE
    for i in range(length - 1):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators

        list[i] = np_SIG
        # if i == (length - 1):
        #     list[0] = np_x
        # else:
        list[i+1] = SIG_dag

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_XX = -Jay * np.exp(PHI * 1.j) * adding_matrixes_in_list(list_2, length)

    return H_XX

def term_two_SIGMA(F, length):

    list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)

    list_2 = []  # Has all terms in a list

    for i in range(length):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators
        list[i] = np_CLK
        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
    #print(list_2)
    CLOCK_SUM = -F * adding_matrixes_in_list(list_2, length)  # NOTE: SAEED DOES NOT HAVE A w PREFACTOR FOR THIS TERM

    print('Clock Sum:')
    print(CLOCK_SUM)
    return CLOCK_SUM

def Liovillan(Hamiltonian, size): #  Liovillan only applied to a matrix converted to vector
    I = sci.eye(3**size, format='csr')
    Hamiltonian = sci.csr_matrix(Hamiltonian)
    S_1 = sci.kron(I, Hamiltonian)
    S_2 = sci.kron((Hamiltonian.T), I)
    #S_1 = np.kron(Hamiltonian, I)
    #S_2 = np.kron(I, Hamiltonian)
    #print(S_1)

    return (S_1+(-1*S_2))


def grahm_schmidt_a(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    #list[0] = np_SIG
    list[0] = np.conjugate(np_SIG.T)  # FOR ALPHA^2 OPERATOR

    strong_mode = tensor_product(list, size)
    O_0 = ((strong_mode.T).flatten())
    gs_basis_non_ortho = [O_0]

    for j in range(1, n+1):
        p = Lio @ gs_basis_non_ortho[j-1]
        #f = (((H @ np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T @ H)).T).flatten()
        gs_basis_non_ortho.append(p)
    for k in range(len(gs_basis_non_ortho)):
        gs_basis_non_ortho[k] = gs_basis_non_ortho[k].tolist()  # Convert to list for easier manipulation
    gs_basis_non_ortho = np.array(gs_basis_non_ortho, dtype=complex).T  # Convert to numpy array for calculations

    n = gs_basis_non_ortho.shape[1]  # Number of vectors in the basis
    Q = np.zeros_like(gs_basis_non_ortho, dtype=np.complex128)  # Initialize orthogonal basis
    R = np.zeros((n, n), dtype=np.complex128)

    for i in range(n):
        v = gs_basis_non_ortho[:, i]
        #print('MAX V ABS:')
        #print(np.max(np.abs(v)))
        #print(v.shape)
        for j in range(i):
            R[j, i] = np.vdot(Q[:, j], v)  # Inner product
            v -= R[j, i] * Q[:, j]
        R[i, i] = np.linalg.norm(v)
        #if R[i, i] < 1e-12:
        #    raise ValueError("Vectors are linearly dependent or nearly so")
    #print('MIN R DIAG: (SIG_Dagger)')
    #print(np.min(np.abs(np.diag(R))))
        Q[:, i] = v / R[i, i]

    matrix_final = np.zeros((n,n), dtype=complex)
    for j in range(n):
        k = Lio @ Q[:, j]
        for i in range(n):
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ k

    return matrix_final

def grahm_schmidt_a_squared(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    #list[0] = np_SIG
    #list[0] = np.conjugate(np_SIG.T)  # FOR ALPHA^2 OPERATOR
    list[0] = np_SIG @ np_SIG

    strong_mode = tensor_product(list, size)
    O_0 = ((strong_mode.T).flatten())
    gs_basis_non_ortho = [O_0]

    for j in range(1, n+1):
        p = Lio @ gs_basis_non_ortho[j-1]
        #f = (((H @ np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[j-1], (3**size, 3**size)).T @ H)).T).flatten()
        gs_basis_non_ortho.append(p)
    for k in range(len(gs_basis_non_ortho)):
        gs_basis_non_ortho[k] = gs_basis_non_ortho[k].tolist()  # Convert to list for easier manipulation
    gs_basis_non_ortho = np.array(gs_basis_non_ortho, dtype=complex).T  # Convert to numpy array for calculations

    n = gs_basis_non_ortho.shape[1]  # Number of vectors in the basis
    Q = np.zeros_like(gs_basis_non_ortho, dtype=np.complex128)  # Initialize orthogonal basis
    R = np.zeros((n, n), dtype=np.complex128)

    for i in range(n):
        v = gs_basis_non_ortho[:, i]
        #print('MAX V ABS:')
        #print(np.max(np.abs(v)))
        for j in range(i):
            R[j, i] = np.vdot(Q[:, j], v)  # Inner product
            v -= R[j, i] * Q[:, j]
        R[i, i] = np.linalg.norm(v)
        #if R[i, i] < 1e-12:
        #    raise ValueError("Vectors are linearly dependent or nearly so")
        Q[:, i] = v / R[i, i]
    #print('MIN R DIAG: (Sig_Squared)')
    #print(np.min(np.abs(np.diag(R))))
    matrix_final = np.zeros((n,n), dtype=complex)
    for j in range(n):
        k = Lio @ Q[:, j]
        for i in range(n):
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ k

    return matrix_final


H_SIGMA = term_one_SIGMA(length=L, Jay=J, PHI=phi) + term_two_SIGMA(F=f, length=L)

matrix_final_a = grahm_schmidt_a(24, H_SIGMA, L)
matrix_final_a_squared = grahm_schmidt_a_squared(24, H_SIGMA, L)
difference_1 = (matrix_final_a - matrix_final_a_squared)
print(np.max(np.abs(difference_1)))

#print(np.abs(np.max(difference_1)))

#print(np.max(np.imag(matrix_final_a)))

# A_a, U_a = np.linalg.eig(matrix_final_a)
# A_a_squared, U_a_squared = np.linalg.eig(matrix_final_a_squared)

# eig_a = U_a[:, np.argmin(np.abs(A_a))]
# eig_a_squared = U_a_squared[:, np.argmin(np.abs(A_a_squared))]
# eig_a = (1/np.linalg.norm(eig_a)) * eig_a
# eig_a_squared = (1/np.linalg.norm(eig_a_squared)) * eig_a_squared

# matrix_a = eig_a.reshape(8, 8).T
# matrix_a_squared = eig_a_squared.reshape(8, 8).T
# difference = ((matrix_a @ matrix_a) - matrix_a_squared)

# print(np.max(np.abs(difference)))

