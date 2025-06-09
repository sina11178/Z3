import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import scipy.sparse as sci
from scipy.linalg import qr
#import scipy

# PARAMETERS TO DEFINE IN HAMILTONIAN
J = 1
phi = np.pi/6  # PHASE FACTOR
f = .1
L = 5  # SYSTEM SIZE
gamma = 1  # DETERMINES HERMICITY

w = np.exp((2 * np.pi * 1.j)/3)  # For defining Sigma operator

sigma = [
    [1, 0, 0],
    [0, w, 0],
    [0, 0, w**2]
]

clock = [
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]
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

    return CLOCK_SUM


H_SIGMA = term_one_SIGMA(length=L, Jay=J, PHI=phi) + term_two_SIGMA(F=f, length=L)



# NOTE: THE FOLLOWING IS THE HAMILTONIAN IN TAU BASIS


def term_one_TAU(length, Jay, PHI):

    list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)
    list_2 = []  # Has all terms in a list

    TAU_dag = np.conjugate(np_CLK.T)
    # LIST OF DIFFERENT MATRICES ASSOCIATED WITH SITES CREATED HERE
    for i in range(length - 1):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators
        list[i] = np_CLK
        # if i == (length - 1):
        #     list[0] = np_x
        # else:
        list[i+1] = TAU_dag

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_XX = -Jay * np.exp(PHI * 1.j) * adding_matrixes_in_list(list_2, length)

    return H_XX

def term_two_TAU(F, length):

    list = []
    #for k in range(length): # Initializes list to idententity operators
    #    list.append(np_I)

    list_2 = []  # Has all terms in a list

    for i in range(length):
        #for j in range(length): # Initializes list to idententity operators
        #    list[j] = np_I
        list = [np_I] * length  # Initializes list to identity operators
        list[i] = np_SIG
        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
    #print(list_2)
    SIG_SUM = -F * adding_matrixes_in_list(list_2, length)  # NOTE: SAEED DOES NOT HAVE A w PREFACTOR FOR THIS TERM

    return SIG_SUM

H_TAU = term_one_TAU(length=L, Jay=J, PHI=phi) + term_two_TAU(F=f, length=L)

#print(H_TAU.dtype)

#H_SIGMA = np.round(H_SIGMA, decimals=14)
#H_TAU = np.round(H_TAU, decimals=14)


A_SIGMA, U_SIGMA = np.linalg.eig(H_SIGMA)  # NOTE: EIGENVALUES AND EIGENVECTORS OF HAMILTONIAN IN SIGMA BASIS
A_TAU, U_TAU = np.linalg.eig(H_TAU)      # NOTE: EIGENVALUES AND EIGENVECTORS OF HAMILTONIAN IN TAU BASIS

#A_SIGMA_ABS_SORTED = np.argsort(np.abs(A_SIGMA))  # SORTS EIGENVALUES IN ABSOLUTE VALUE
#A_TAU_ABS_SORTED = np.argsort(np.abs(A_TAU))      # SORTS EIGENVALUES IN ABSOLUTE VALUE

#print(np.abs(A_SIGMA))
#print(np.abs(A_TAU))
#print(A_SIGMA)
#print(A_TAU)
#idx_SIGMA = np.lexsort((A_SIGMA.imag, A_SIGMA.real))  # SORTS EIGENVALUES IN SIGMA BASIS
#idx_TAU = np.lexsort((A_TAU.imag, A_TAU.real))      # SORTS EIGENVALUES IN TAU BASIS


A_SIGMA_SORT = np.sort_complex(A_SIGMA)  # SORTS EIGENVALUES IN ABSOLUTE VALUE
A_TAU_SORT = np.sort_complex(A_TAU)      # SORTS EIGENVALUES IN ABSOLUTE VALUE
#print(A_SIGMA_SORT)
#print(A_TAU_SORT)
difference = A_SIGMA_SORT - A_TAU_SORT  # DIFF BETWEEN EIGENVALUES IN SIGMA AND TAU BASIS

#print(difference)

# def eigenvalues_match_unordered(a, b, tol=1e-13):
#     a = np.array(a)
#     b = np.array(b)
#     if len(a) != len(b):
#         return False

#     used = np.full(len(b), False)

#     for x in a:
#         diffs = np.abs(b - x)
#         mask = (diffs < tol) & (~used)
#         if not np.any(mask):
#             return False
#         used[np.argmax(mask)] = True  # Mark one match as used

#     return True

def eigenvalues_match(arr1, arr2, tol=1e-15):
    arr1 = list(arr1)
    arr2 = list(arr2)

    if len(arr1) != len(arr2):
        return False

    used = [False] * len(arr2)
    for val1 in arr1:
        matched = False
        for i, val2 in enumerate(arr2):
            if not used[i] and np.isclose(val1, val2, atol=tol):
                used[i] = True
                matched = True
                break
        if not matched:
            return False
    return True


#diff = np.all(difference.real < 1e-15) and np.all(difference.imag < 1e-15)  # CHECKS IF DIFFERENCE IS SMALLER THAN 1e-10

# for i in difference.real:
#        if abs(i) > 1e-14:
#            print(i)

# for j in difference.imag:
#        if abs(j) > 1e-14:
#            print(j) 


#A_TAU_SORT[2] = 2 


diff = eigenvalues_match(A_SIGMA_SORT, A_TAU_SORT)

print(diff)