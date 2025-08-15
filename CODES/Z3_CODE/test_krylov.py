import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix

#import scipy

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 1
phi = 0 * np.pi/6  # PHASE FACTOR
f = .1
L = 5 # SYSTEM SIZE
#gamma = 1  # DETERMINES HERMICITY

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

# np_SIG = np.array(sigma, dtype=np.complex128)  # Sigma operator
# np_CLK = np.array(clock, dtype=np.complex128)  # Clock operator
# np_I = np.array(I, dtype=np.complex128)  # Identity operator

np_SIG = sci.csr_matrix(sigma, dtype=np.complex128)  # Sigma operator
np_CLK = sci.csr_matrix(clock, dtype=np.complex128)  # Clock operator
np_I = sci.csr_matrix(I, dtype=np.complex128)  # Identity operator

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
            # final_matrix = np.kron(final_matrix, array[i])
            final_matrix = sci.kron(final_matrix, array[i])
    return final_matrix


#print(-1*tensor_product([np.conjugate(np_SIG.T), np_SIG], 2))    


def adding_matrixes_in_list(array, length):  #  Used in all terms
    # matrix = np.zeros((3**length, 3**length), dtype=np.complex128)
    matrix = sci.csr_matrix((3**length, 3**length), dtype=np.complex128)
    
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
    SIG_dag = np_SIG.conjugate().T
    #SIG_dag = np.conjugate(np_SIG.T)
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

def term_two_SIGMA(F, length):

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

def Liovillan(Hamiltonian, size): #  Liovillan only applied to a matrix converted to vector
    ID = sci.eye(3**size, format='csr')
    Hamiltonian = sci.csr_matrix(Hamiltonian)
    S_1 = sci.kron(ID, Hamiltonian)
    S_2 = sci.kron((Hamiltonian.T), ID)
    #S_1 = np.kron(Hamiltonian, I)
    #S_2 = np.kron(I, Hamiltonian)
    #print(S_1)

    return (S_1+(-1*S_2))


def GS_NEW(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_SIG
    #list[0] = np.conj(np_SIG.T)
    #list[0] = np_SIG @ np_SIG  # FOR ALPHA^2 OPERATOR
    strong_mode = tensor_product(list, size)
    O_0 = ((strong_mode.T).flatten())
    gs_basis_non_ortho_zero_charge = [O_0]
    gs_basis_non_ortho_one_charge = []
    gs_basis_non_ortho_two_charge = []
    for i in range(1, n+1):
        if i%3 == 0:
            gs_basis_non_ortho_zero_charge.append(Lio @ gs_basis_non_ortho_two_charge[len(gs_basis_non_ortho_two_charge) - 1])
        if i%3 == 1:
            gs_basis_non_ortho_one_charge.append(Lio @ gs_basis_non_ortho_zero_charge[len(gs_basis_non_ortho_zero_charge) - 1])
        if i%3 == 2:
            gs_basis_non_ortho_two_charge.append(Lio @ gs_basis_non_ortho_one_charge[len(gs_basis_non_ortho_one_charge) - 1])
    
    U = []
    one_len = len(gs_basis_non_ortho_one_charge)
    two_len = len(gs_basis_non_ortho_two_charge)
    for j in range(len(gs_basis_non_ortho_zero_charge)):
        v_zero = gs_basis_non_ortho_zero_charge[j]
        if j < one_len:
            v_one = gs_basis_non_ortho_one_charge[j]
        if j < two_len:
            v_two = gs_basis_non_ortho_two_charge[j]
        for k in range(j):
            proj = ((np.dot(np.conjugate(v_zero), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
            v_zero -= proj
            if j < one_len:
                proj = ((np.dot(np.conjugate(v_one), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
                v_one -= proj
            if j < two_len:
                proj = ((np.dot(np.conjugate(v_two), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
                v_two -= proj
        U.append(v_zero/np.linalg.norm(v_zero))
        if j < one_len:
            U.append(v_one/np.linalg.norm(v_one))
        if j < two_len:
            U.append(v_two/np.linalg.norm(v_two))

    p = len(U)
    matrix_final = np.zeros((p,p), dtype=complex)
    for j in range(p):
        k = Lio @ U[j]
        for i in range(p):
            matrix_final[i, j] = np.conjugate(U[i]) @ k

    return matrix_final

def GS_NEW_squared(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    #list[0] = np_SIG
    list[0] = np_SIG @ np_SIG  # FOR ALPHA^2 OPERATOR
    strong_mode = tensor_product(list, size)
    O_0 = ((strong_mode.T).flatten())
    gs_basis_non_ortho_zero_charge = [O_0]
    gs_basis_non_ortho_one_charge = []
    gs_basis_non_ortho_two_charge = []
    for i in range(1, n+1):
        if i%3 == 0:
            gs_basis_non_ortho_zero_charge.append(Lio @ gs_basis_non_ortho_two_charge[len(gs_basis_non_ortho_two_charge) - 1])
        if i%3 == 1:
            gs_basis_non_ortho_one_charge.append(Lio @ gs_basis_non_ortho_zero_charge[len(gs_basis_non_ortho_zero_charge) - 1])
        if i%3 == 2:
            gs_basis_non_ortho_two_charge.append(Lio @ gs_basis_non_ortho_one_charge[len(gs_basis_non_ortho_one_charge) - 1])
    
    U = []
    one_len = len(gs_basis_non_ortho_one_charge)
    two_len = len(gs_basis_non_ortho_two_charge)
    for j in range(len(gs_basis_non_ortho_zero_charge)):
        v_zero = gs_basis_non_ortho_zero_charge[j]
        if j < one_len:
            v_one = gs_basis_non_ortho_one_charge[j]
        if j < two_len:
            v_two = gs_basis_non_ortho_two_charge[j]
        for k in range(j):
            proj = ((np.dot(np.conjugate(v_zero), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
            v_zero -= proj
            if j < one_len:
                proj = ((np.dot(np.conjugate(v_one), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
                v_one -= proj
            if j < two_len:
                proj = ((np.dot(np.conjugate(v_two), U[k]))/(np.dot(np.conjugate(U[k]), U[k]))) * U[k]
                v_two -= proj
        U.append(v_zero/np.linalg.norm(v_zero))
        if j < one_len:
            U.append(v_one/np.linalg.norm(v_one))
        if j < two_len:
            U.append(v_two/np.linalg.norm(v_two))

    p = len(U)
    matrix_final = np.zeros((p,p), dtype=complex)
    for j in range(p):
        k = Lio @ U[j]
        for i in range(p):
            matrix_final[i, j] = np.conjugate(U[i]) @ k

    return matrix_final

def GS_NEWER(n, H, size):
    Lio = Liovillan(H, size)

    # Build initial operator with a single-site non-identity
    op_list = [np_I for _ in range(size)]
    op_list[0] = np_SIG
    strong_mode = tensor_product(op_list, size)
    O_0 = strong_mode.T.flatten()

    # Initialize Krylov vectors in different charge sectors
    gs_basis = {
        0: [O_0],
        1: [],
        2: []
    }

    # Bootstrap first two steps manually
    gs_basis[1].append(Lio @ O_0)
    gs_basis[2].append(Lio @ gs_basis[1][0])

    # Generate Krylov vectors cyclically
    for i in range(3, n+1):
        sector = i % 3
        prev_sector = (i - 1) % 3
        gs_basis[sector].append(Lio @ gs_basis[prev_sector][-1])

    # Flatten and orthonormalize all vectors in round-robin order
    U = []
    max_len = max(len(gs_basis[0]), len(gs_basis[1]), len(gs_basis[2]))
    for j in range(max_len):
        for sector in [0, 1, 2]:
            if j < len(gs_basis[sector]):
                v = gs_basis[sector][j]
                for u in U:
                    proj = (np.vdot(u, v) / np.vdot(u, u)) * u
                    v -= proj
                norm = np.linalg.norm(v)
                if norm > 1e-12:
                    U.append(v / norm)

    # Build effective Liouvillian matrix in Krylov basis
    p = len(U)
    matrix_final = np.zeros((p, p), dtype=complex)
    for j in range(p):
        k = Lio @ U[j]
        for i in range(p):
            matrix_final[i, j] = np.vdot(U[i], k)

    return matrix_final

def GS_NEWER_squared(n, H, size):
    Lio = Liovillan(H, size)

    # Build initial operator with a single-site non-identity
    op_list = [np_I for _ in range(size)]
    op_list[0] = np_SIG @ np_SIG
    strong_mode = tensor_product(op_list, size)
    O_0 = strong_mode.T.flatten()

    # Initialize Krylov vectors in different charge sectors
    gs_basis = {
        0: [O_0],
        1: [],
        2: []
    }

    # Bootstrap first two steps manually
    gs_basis[1].append(Lio @ O_0)
    gs_basis[2].append(Lio @ gs_basis[1][0])

    # Generate Krylov vectors cyclically
    for i in range(3, n+1):
        sector = i % 3
        prev_sector = (i - 1) % 3
        gs_basis[sector].append(Lio @ gs_basis[prev_sector][-1])

    # Flatten and orthonormalize all vectors in round-robin order
    U = []
    max_len = max(len(gs_basis[0]), len(gs_basis[1]), len(gs_basis[2]))
    for j in range(max_len):
        for sector in [0, 1, 2]:
            if j < len(gs_basis[sector]):
                v = gs_basis[sector][j]
                for u in U:
                    proj = (np.vdot(u, v) / np.vdot(u, u)) * u
                    v -= proj
                norm = np.linalg.norm(v)
                if norm > 1e-12:
                    U.append(v / norm)

    # Build effective Liouvillian matrix in Krylov basis
    p = len(U)
    matrix_final = np.zeros((p, p), dtype=complex)
    for j in range(p):
        k = Lio @ U[j]
        for i in range(p):
            matrix_final[i, j] = np.vdot(U[i], k)

    return matrix_final


n = []
for i in range(L):
    n.append(np_I)
n[0] = np_CLK @ np_CLK
n[1] = np_SIG

alph_3 = tensor_product(n, L)

n[0] = np_CLK @ np_CLK
n[1] = np_CLK @ np_CLK
n[2] = np_SIG
alpha_5 = tensor_product(n, L)



def grahm_schmidt(n, H, size):
    #Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_SIG

    #seeding = tensor_product(list, size)
    seeding = tensor_product(list, size).toarray()
    #O_0 = (seeding.T.reshape((-1, 1))).toarray()
    O_0 = ((seeding.T).flatten())
    gs_basis_non_ortho = [O_0]

    for P in range(1, n+1):
        #K = Lio @ gs_basis_non_ortho[P-1]
        K = (((H @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()
        #K = (((H.conjugate().T @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()
        gs_basis_non_ortho.append(K)
    #print('EXPERIMENT/TEST')
    #print(type(gs_basis_non_ortho[0]))
    #print(gs_basis_non_ortho[0].shape)
    gs_basis_non_ortho_matrix = np.column_stack(gs_basis_non_ortho)

    # CHECKS RANK OF MATRIX W/ NON-ORTHOGONALIZED BASIS VECTORS
    rank = rank_test(gs_basis_non_ortho_matrix)
    print('Matrix Rank:')
    print(rank)

    Q, R = np.linalg.qr(gs_basis_non_ortho_matrix)

    l = len(gs_basis_non_ortho)
    matrix_final = np.zeros((l,l), dtype=complex)
    for j in range(l):
        #c = Lio @ Q[:, j]
        c = (((H @ np.reshape(Q[:, j], (3**size, 3**size)).T) -
            (np.reshape(Q[:,j], (3**size, 3**size)).T @ H)).T).flatten()
        for i in range(l):
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ c

    #print(np.abs(matrix_final[16, :]))
    #print(np.abs(matrix_final[:, 16]))
    truncation = 0
    for i in range(l):
        if np.all(np.abs(matrix_final[i, :]) < 1e-1):
            print('Row where all elements are zero:')
            print(str(i) + ' (zero indexing)')
            if truncation == 0:
                truncation = i
            #print(truncation)
        if np.all(np.abs(matrix_final[:, i]) < 1e-1):
            print('column where all elements are zero:')
            print(str(i) + ' (zero indexing)')
            if truncation == 0:
                truncation = i
    
    range_list = np.arange(0, n+1, 1)
    plt.figure()
    plt.plot(range_list, np.abs(matrix_final[truncation, :]), marker = ".", linestyle = '-', color = 'r', label = 'ROW ' + str(truncation))
    plt.plot(range_list, np.abs(matrix_final[:, truncation]), marker = ".", linestyle = '-', color = 'g', label = 'COLUMN ' + str(truncation))
            #print(truncation)
    plt.legend()
    plt.title('NOTE: THESE SHOW ABS VALUES: L = ' + str(size))
    plt.show()
    return matrix_final



def grahm_schmidt_squared(n, H, size):
    Lio = Liovillan(H, size)
    list = []
    for i in range(size):
        list.append(np_I)
    #list[0] = np_SIG
    #list[0] = np.conjugate(np_SIG).T # FOR ALPHA^2 OPERATOR
    list[0] = np_SIG.conjugate().T

    #strong_mode = tensor_product(list, size)
    #strong_mode = tensor_product(list, size).toarray()
    #O_0 = (strong_mode.T.reshape((-1, 1))).T
    O_0 = ((tensor_product(list, size).toarray().T).flatten())
    gs_basis_non_ortho = [O_0]

    for P in range(1, n+1):
        #K = Lio @ gs_basis_non_ortho[P-1]
        K = (((H @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()
        #K = (((H.conjugate().T @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()

        gs_basis_non_ortho.append(K)

    gs_basis_non_ortho_matrix = np.column_stack(gs_basis_non_ortho)
    Q, R = np.linalg.qr(gs_basis_non_ortho_matrix)

    l = len(gs_basis_non_ortho)
    matrix_final = np.zeros((l,l), dtype=complex)
    for j in range(l):
        #c = Lio @ Q[:, j]
        c = (((H @ np.reshape(Q[:, j], (3**size, 3**size)).T) -
            (np.reshape(Q[:,j], (3**size, 3**size)).T @ H)).T).flatten()
        for i in range(l):
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ c
    return matrix_final

def rank_test(stacked_matrix):
    return np.linalg.matrix_rank(stacked_matrix)


def Liovillan_new(Hamiltonian, size): #  Liovillan only applied to a matrix converted to vector
    ID = sci.eye(3**size, format='csr')
    Hamiltonian = sci.csr_matrix(Hamiltonian)
    S_1 = sci.kron(ID, Hamiltonian.conjugate().T)
    S_2 = sci.kron((Hamiltonian.T), ID)
    #S_1 = np.kron(Hamiltonian, I)
    #S_2 = np.kron(I, Hamiltonian)
    #print(S_1)

    return (S_1+(-1*S_2))

def Experiment(n, H, size):
    Lio = Liovillan_new(H, size)
    list = []
    list_2 = []
    for i in range(size):
        list.append(np_I)
        #list_2.append(np_I)
    #list[0] = np_SIG
    #list[0] = np.conjugate(np_SIG).T # FOR ALPHA^2 OPERATOR
    list[0] = np_SIG.conjugate().T
    #list_2[0] = np_SIG

    #strong_mode = tensor_product(list, size)
    #strong_mode = tensor_product(list, size).toarray()
    #O_0 = (strong_mode.T.reshape((-1, 1))).T
    O_0 = ((tensor_product(list, size).toarray().T).flatten())
    #O_0_2 = ((tensor_product(list_2, size).toarray().T).flatten())
    gs_basis_non_ortho = [O_0]
    #gs_basis_non_ortho_2 = [O_0_2]

    for P in range(1, n+1):
        K = Lio @ gs_basis_non_ortho[P-1]
        #K = (((H @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()
        #K = (((H.conjugate().T @ np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T) - (np.reshape(gs_basis_non_ortho[P-1], (3**size, 3**size)).T @ H)).T).flatten()

        gs_basis_non_ortho.append(K)

    gs_basis_non_ortho_matrix = np.column_stack(gs_basis_non_ortho)
    Q, R = np.linalg.qr(gs_basis_non_ortho_matrix)

    l = len(gs_basis_non_ortho)
    matrix_final = np.zeros((l,l), dtype=complex)
    for j in range(l):
        c = Lio @ Q[:, j]
        #c = (((H @ np.reshape(Q[:, j], (3**size, 3**size)).T) -
        #    (np.reshape(Q[:,j], (3**size, 3**size)).T @ H)).T).flatten()
        for i in range(l):
            #print(np.conjugate(Q[:, i]) @ Q[:, j])
            matrix_final[i, j] = np.conjugate(Q[:, i]) @ c
    return matrix_final


GAMM = 0
H_SIGMA = term_one_SIGMA(length=L, Jay=J, PHI=phi) + term_two_SIGMA(F=f, length=L)
H_SIGMA = H_SIGMA + GAMM * (H_SIGMA.conjugate().T)  # Make Hermitian
#print('Difference between SIG^2 and SIG^t:')
#print(np.max(np.abs((np_SIG @ np_SIG) - np.conjugate(np_SIG).T)))
#print(np.max(((0.5*(np_SIG @ np_SIG + np.conjugate(np_SIG).T)) - np.conjugate(np_SIG).T)))

p = 200 # NUMBER OF ITERATIONS
#matrix_final_a = GS_NEWER(p, H_SIGMA, L)
matrix_final_a = Experiment(p, H_SIGMA, L)  # SEED OPERATOR IS ALPHA_1
#matrix_final_a = grahm_schmidt(p, H_SIGMA, L)  # SEED OPERATOR IS ALPHA_1
matrix_final_a_squared = Experiment(p, H_SIGMA, L)  # SEED OPERATOR IS (ALPHA_1)^2
#matrix_final_a_squared = grahm_schmidt_squared(p, H_SIGMA, L)  # SEED Operator IS (ALPHA_1)^2

#print(type(matrix_final_a[1,1]))
if GAMM == 1:
    A, U = np.linalg.eigh(matrix_final_a)
    A_sqr, U_sqr = np.linalg.eigh(matrix_final_a_squared)
else:
    A, U = np.linalg.eig(matrix_final_a)  # Get eigenvalues and eigenvectors
    A_sqr, U_sqr = np.linalg.eig(matrix_final_a_squared)

EIG_SQR = U_sqr[:, np.argmin(np.abs(A_sqr))]  # Smallest eigenvalue eigenvector of squared operator
EIG_SINGLE = U[:, np.argmin(np.abs(A))]  # Smallest eigenvalue eigenvector of single operator

plt.figure()
plt.plot(np.real(A_sqr), np.imag(A_sqr), marker = ".", linestyle = 'None', color = 'r', label = 'Eigenvalues of (alpha_1)^2')
plt.legend()
plt.show()
#print('\nABS OF DIFF OF MIN EIGENVALUE EIGENVECTOR')
#print(np.max(np.abs(EIG_SQR) - np.abs(EIG_SINGLE)))


difference_1 = (np.abs(matrix_final_a) - np.abs(matrix_final_a_squared))
#print('\nDifference between Krylov Matrix Outputs:')
#print(np.max(np.abs(difference_1)))
#rint(np.max(matrix_final_a.real - matrix_final_a_squared.real))
F_2 = np.abs(matrix_final_a)
F_1 = np.abs(matrix_final_a_squared)
plt.subplot(1, 2, 1)
plt.imshow(F_1, cmap='seismic', aspect='auto')
plt.colorbar(label='Value')
#plt.clim(0, 3.2)
plt.title('alpha_1_squared (abs value of matrix)')
plt.xlabel('Column Index')
plt.ylabel('Row Index')

plt.subplot(1, 2, 2)
plt.imshow(F_2, cmap='seismic', aspect='auto')
plt.colorbar(label='Value')
#plt.clim(0, 1.8)
plt.title('alpha_1 (abs value of matrix)')
plt.xlabel('Column Index')
plt.ylabel('Row Index')

#plt.clim(-.1, .1)
plt.tight_layout()
plt.savefig('Krylov_Matrices.pdf')
plt.show()

#A, U = np.linalg.eig((matrix_final_a))  # Get eigenvalues and eigenvectors
#A_sqr, U_sqr = np.linalg.eig(matrix_final_a_squared)
#print(np.imag(A))
in_list = np.arange(0, p+1, 1)
EIG_STATES = []
# EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
sorted = np.sort(np.abs(A))
sort_sqr = np.sort(np.abs(A_sqr))

for i in range(len(A)):
    #EIG_STATE = (U[:, np.abs(A).argsort()[i]])
    EIG_STATE = U[:, np.abs(A).argsort()[i]]
    EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
    #EIG_STATE = np.abs(EIG_STATE)
    EIG_STATE = np.abs(EIG_STATE)
    #print(np.linalg.norm(EIG_STATE))
    EIG_STATES.append(EIG_STATE)

plt.plot(in_list, EIG_STATES[0], marker = ".", linestyle = '-', color = 'r', label = 'a_1 E = ' + str(np.min(np.abs(A))))
plt.plot(in_list, np.abs(EIG_SQR), marker = ".", linestyle = '-', color = 'g', label = '(a_1)^dagger E = ' + str(np.min(np.abs(A_sqr))))
plt.legend()
plt.ylim(0, 1)
plt.savefig('Krylov_Eigenvector_0.pdf')


# ncols = 4  # columns per row
# nrows = int(np.ceil(len(EIG_STATES) / ncols))
# fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows))
# axes = axes.flatten()
# plt.suptitle('Hello')
# for i in range(int(len(EIG_STATES))):
#         ax = axes[i]
#         ax.plot(in_list, EIG_STATES[i], marker = ".", linestyle = '-', label = 'E = ' + str(np.round(sorted[i], 5)) + " (" + str(i) + ")")
#         #title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Eigenvector graph (L = " + str(size) + ")"

#         ax.set_xlabel("n")
#         ax.set_ylabel("Eigenvector Components")

#         ax.set_ylim(0, 1)
#         ax.legend()
# plt.show()


