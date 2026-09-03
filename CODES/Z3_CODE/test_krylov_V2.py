import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
from scipy.sparse.linalg import eigs
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh, eigs, LinearOperator
from scipy.sparse import issparse
#import scipy

#import scipy

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 1
phi = 0 * np.pi/6  # PHASE FACTOR
f = .1
L = 4 # SYSTEM SIZE
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

def liouv_twisted(Hamiltonian, size):
    Hamiltonian  = sci.csr_matrix(Hamiltonian)
    ID = sci.eye(3**size, format='csr')
    S_1 = sci.kron(ID, Hamiltonian.conjugate().T)
    S_2 = sci.kron(Hamiltonian.T, ID)
    return (S_1+(-1*S_2))

def liouvillian_linop_twisted(Hamiltonian, size):
    """
    Return a LinearOperator that performs the commutator action on flattened operators:
        x -> (H @ M - M @ H).ravel(order='F')
    without building the full d^2 x d^2 superoperator.
    """
    #H = sci.csr_matrix(Hamiltonian)
    H = Hamiltonian.toarray()
    d = H.shape[0]
    def mv(x):
        # ensure x is a 1-D numpy array
        x = np.asarray(x, dtype=np.complex128).ravel()
        # reshape to (d,d) using the same order you use when flattening elsewhere
        M = x.reshape((d, d), order='F')   # use 'F' if you ravel(..., order='F')
        # ensure M is dense ndarray (needed to do M @ H where H is sparse)
        if issparse(M):
            M = M.toarray()
        # compute commutator; H is sparse so H.dot(M) returns dense ndarray
        left = H.conjugate().T @ M     # shape (d,d) dense
        right = M @ H       # shape (d,d) dense
        C = left - right
        return C.ravel(order='F')
    return LinearOperator((d*d, d*d), matvec=mv, dtype=np.complex128)


def grahm_schmidt(n, H, size):
    Lio = liouvillian_linop_twisted(H, size) 
    d = 3**size
    list = []
    list_squared = []  # FOR ALPHA^2 OPERATOR
    for i in range(size):
        list.append(np_I)
        list_squared.append(np_I)  # FOR ALPHA^2 OPERATOR
    list[0] = np_SIG
    list_squared[0] = (np_SIG.conjugate().T)  # FOR ALPHA^2 OPERATOR
    seeding = tensor_product(list, size).toarray().ravel(order = 'F')  # Converts to array for matrix operations
    seeding_squared = tensor_product(list_squared, size).toarray().ravel(order = 'F')  # FOR ALPHA^2 OPERATOR

    Q = [seeding/np.linalg.norm(seeding)]  # ORTHONORMAL BASIS
    Q_squared = [seeding_squared/np.linalg.norm(seeding_squared)]  # FOR ALPHA^2 OPERATOR
    v = seeding
    v_squared = seeding_squared
    for k in range(n):
        w = Lio.matvec(v)   # L sparse
        w_squared = Lio.matvec(v_squared)
        for q in Q:
            alpha =  np.vdot(q, w)   # inner product
            w = w - alpha * q           # sparse subtraction
        for q_squared in Q_squared:
            alpha_sq = np.vdot(q_squared, w_squared)
            w_squared = w_squared - alpha_sq * q_squared

        w = w / np.linalg.norm(w)
        w_squared = w_squared/np.linalg.norm(w_squared)
        Q.append(w)
        Q_squared.append(w_squared)
        v = w
        v_squared = w_squared
    #Q = [seeding/np.linalg.norm(seeding)] + Q
    #Q_squared = [seeding_squared/np.linalg.norm(seeding_squared)] + Q_squared
    l = len(Q)
    matrix_final = np.zeros((l,l), dtype=complex)  # FOR ALPHA OPERATOR
    matrix_final_squared = np.zeros((l,l), dtype=complex)  # FOR ALPHA^2 OPERATOR
    for j in range(l):
        c = Lio.matvec(Q[j])
        c_squared = Lio.matvec(Q_squared[j])
        for i in range(l):
            matrix_final[i, j] = np.vdot(Q[i], c)
            matrix_final_squared[i, j] = np.vdot(Q_squared[i], c_squared)
    return [matrix_final, matrix_final_squared]


def main_1():
    fig, ax = plt.subplots(2, 2)
    for i, l in enumerate([3, 5]):
        H = term_one_SIGMA(l, Jay=1, PHI = 0*np.pi/6) + term_two_SIGMA(0.1, l)
        Gamm = 0
        H = H + Gamm * H.conjugate().T
        n = 26
        matrix_final, matrix_final_squared = grahm_schmidt(n, H, l)
        ax[1][i].imshow(np.real(matrix_final_squared), cmap='viridis', interpolation='nearest')

        if Gamm != 0:
            A, U = np.linalg.eig(matrix_final_squared)
        else:
            A, U = np.linalg.eigh(matrix_final_squared)

        eig_vector = U[:, np.argmin(np.abs(A))]
        ax[0][0].plot(np.abs(eig_vector), label=f'L={l}, E = {A[np.argmin(np.abs(A))]:.4f}')
    ax[0][0].legend()

    plt.tight_layout()
    plt.show()

main_1()