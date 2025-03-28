'''
Reference: https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.124.206803
L = 5
J_z = 0
g = 0.2
gamma = 0.8

'''
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd

g = .7
J_z = .7
L = 6
gam = 1
Hamiltonian_components = []

parameter_array = [[g, J_z, gam], 
[1, .4, .4], [1, .5, .5], [1, .6, .6], [1, .7, .7]]  # PARAMETERS [g, J_z, gam]

sig_x = [  # PAULI X
    [0, 1],
    [1, 0]
]

sig_y = [
    [0, -1.j],
    [1.j, 0]
]

sig_z = [  # PAULI Z
    [1, 0],
    [0, -1]
]

I = [
    [1, 0],
    [0, 1]
]

# NUMPY MATRICES
np_x = np.matrix(sig_x)
np_y = np.matrix(sig_y)
np_z = np.matrix(sig_z)
np_I = np.matrix(I)


# Hamiltonian construction

def tensor_product(array, length):  # Calculates tensor product of all matrices needed

    final_matrix = np.zeros((2**length, 2**length))
    for i in range(length):
        if i == 0:
            final_matrix = array[i]
            #print(final_matrix)
        else:
            final_matrix = np.kron(final_matrix, array[i])
    return final_matrix

def adding_matrixes_in_list(array, length):  #  Used in all terms
    matrix = np.zeros((2**length, 2**length))
    
    for i in range(len(array)):
        #print(len(array))
        matrix = matrix + array[i]

    return matrix


def term_one(gamma, length):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)
    list_2 = []  # Has all terms in a list

    # LIST OF DIFFERENT MATRICES ASSOCIATED WITH SITES CREATED HERE
    for i in range(length - 1):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I

        list[i] = np_x
        # if i == (length - 1):
        #     list[0] = np_x
        # else:
        list[i+1] = np_x

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_XX = ((1+gamma)/2) * adding_matrixes_in_list(list_2, length)

    return H_XX

def term_two(gamma, length):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)
    list_2 = []  # Has all terms in a list

    for i in range(length - 1):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I

        list[i] = np_y
        # if i == (length - 1):
        #     list[0] = np_y
        # else:
        list[i+1] = np_y

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list

    H_YY = ((1-gamma)/2) * adding_matrixes_in_list(list_2, length)

    return H_YY

def term_three(Jay_Z, length):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)
    list_2 = []  # Has all terms in a list

    for i in range(length - 1):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I

        list[i] = np_z
        # if i == (length - 1):
        #     list[0] = np_z
        # else:
        list[i+1] = np_z

        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
        #print(list_2[i])

    H_ZZ = Jay_Z * adding_matrixes_in_list(list_2, length)

    return H_ZZ

def term_four(G, length):

    list = []
    for k in range(length): # Initializes list to idententity operators
        list.append(np_I)

    list_2 = []  # Has all terms in a list

    for i in range(length):
        for j in range(length): # Initializes list to idententity operators
            list[j] = np_I
        list[i] = np_z
        list_2.append(tensor_product(list, length)) #  Performs tensor product and puts them in a list
    #print(list_2)
    H_Z = G * adding_matrixes_in_list(list_2, length)

    return H_Z

H_XX1 = term_one(gam, L)
H_YY1 = term_two(gam, L)
H_ZZ1 = term_three(J_z, L)
H_Z1 = term_four(g, L)

H_1 = H_XX1 + H_YY1 + H_ZZ1 + H_Z1  #  HAMILTONIAN ONE  L
'''
H_XX2 = term_one(parameter_array[1][0], L)
H_YY2 = term_two(parameter_array[1][0], L)
H_ZZ2 = term_three(parameter_array[1][1], L)
H_Z2 = term_four(parameter_array[1][2], L)

H_2 = H_XX2 + H_YY2 + H_ZZ2 + H_Z2  #  HAMILTONIAN ONE  L

H_XX3 = term_one(parameter_array[2][0], L)
H_YY3 = term_two(parameter_array[2][0], L)
H_ZZ3 = term_three(parameter_array[2][1], L)
H_Z3 = term_four(parameter_array[2][2], L)

H_3 = H_XX3 + H_YY3 + H_ZZ3 + H_Z3  #  HAMILTONIAN ONE  L

H_XX4 = term_one(parameter_array[3][0], L)
H_YY4 = term_two(parameter_array[3][0], L)
H_ZZ4 = term_three(parameter_array[3][1], L)
H_Z4 = term_four(parameter_array[3][2], L)

H_4 = H_XX4 + H_YY4 + H_ZZ4 + H_Z4  #  HAMILTONIAN ONE  L

H_XX5 = term_one(parameter_array[4][0], L)
H_YY5 = term_two(parameter_array[4][0], L)
H_ZZ5 = term_three(parameter_array[4][1], L)
H_Z5 = term_four(parameter_array[4][2], L)

H_5 = H_XX5 + H_YY5 + H_ZZ5 + H_Z5  #  HAMILTONIAN ONE  L
'''
#H_XX6 = term_one(parameter_array[5][0], L)
#H_YY6 = term_two(parameter_array[5][0], L)
#H_ZZ6 = term_three(parameter_array[5][1], L)
#H_Z6 = term_four(parameter_array[5][2], L)

#H_6 = H_XX6 + H_YY6 + H_ZZ6 + H_Z6  #  HAMILTONIAN ONE  L

#H_XX7 = term_one(parameter_array[6][0], L)
#H_YY7 = term_two(parameter_array[6][0], L)
#H_ZZ7 = term_three(parameter_array[6][1], L)
#H_Z7 = term_four(parameter_array[6][2], L)

#H_7 = H_XX7 + H_YY7 + H_ZZ7 + H_Z7  #  HAMILTONIAN ONE  L

#H_XX8 = term_one(parameter_array[7][0], L)
#H_YY8 = term_two(parameter_array[7][0], L)
#H_ZZ8 = term_three(parameter_array[7][1], L)
#H_Z8 = term_four(parameter_array[7][2], L)

#H_8 = H_XX8 + H_YY8 + H_ZZ8 + H_Z8  #  HAMILTONIAN ONE  L

# H_XX2 = term_one(gam, L+2)
# H_YY2 = term_two(gam, L+2)
# H_ZZ2 = term_three(J_z, L+2)
# H_Z2 = term_four(g, L+2)

# H_2 = H_XX2 + H_YY2 + H_ZZ2 + H_Z2  # HAMILTONIAN L+1

# H_XX3 = term_one(gam, L+4)
# H_YY3 = term_two(gam, L+4)
# H_ZZ3 = term_three(J_z, L+4)
# H_Z3 = term_four(g, L+4)

# H_3 = H_XX3 + H_YY3 + H_ZZ3 + H_Z3  # HAMILTONIAN L+2

#H_XX4 = term_one(gam, L+3)
#H_YY4 = term_two(gam, L+3)
#H_ZZ4 = term_three(J_z, L+3)
#H_Z4 = term_four(g, L+3)

#H_4 = H_XX4 + H_YY4 + H_ZZ4 + H_Z4  # HAMILTONIAN L+3

#H_XX5 = term_one(gam, L+4)
#H_YY5 = term_two(gam, L+4)
# H_ZZ5 = term_three(J_z, L+4)
# H_Z5 = term_four(g, L+4)

# H_5 = H_XX5 + H_YY5 + H_ZZ5 + H_Z5  # HAMILTONIAN L+3

#H_tot = [H_1, H_2, H_3, H_4, H_5]  # List of all hamiltonians I want to simulate
#L_tot = [L, L, L, L, L]  # Site numbers I want to simulate


# I now must define the time evolution operator

def time_evol_op_at_time_t(hamiltonian, time):
    A, U = np.linalg.eigh(hamiltonian)  # A is vector of eigenvalues; U is matrix of eigenvectors
    U_T = np.conjugate(U.T)  # Conjugate transpose of U
    #A_matrix = exp_diag_matrix(A, time)
    A_matrix = np.diag(np.exp(A* time*(-1.j)))
    evol_to_t = U @ A_matrix @ U_T
    return evol_to_t

def exp_diag_matrix(eigen_array, t):
    matrix = np.identity(len(eigen_array))
    matrix = np.matrix.tolist(matrix) # Convert back to normal matrix
    for k in range(len(eigen_array)):
        matrix[k][k] = np.exp((eigen_array[k])*(-1.j)*t)  # I NEED TO FIGURE OUT HOW TO ADD INDEP VARIABLE

    matrix = np.matrix(matrix)
    return matrix

def conjugate_transpose_time_op(operator):
    return operator.getH()

# Now, we will create the autocorrelator

def autocorrelator(time, size, hamiltonian):  # ONLY FOR SIG1_X
    list = []
    for a in range(size):
        list.append(np_I)


    list[0] = np_x # first site association with sig_x set
    time_evol_op = time_evol_op_at_time_t(hamiltonian, time)  # Finds time evol operator at time t
    non_time_evold_op_sig_x = tensor_product(list, size) #  SIGX_1(0)
    time_evolved_sig_x = (np.conjugate(time_evol_op.T)) @ non_time_evold_op_sig_x @ (time_evol_op)  #  SIGX_1(time)

    # print(time_evol_op.getH() * time_evol_op)
    prefactor = 1/(2**size)
    product_x = time_evolved_sig_x @ non_time_evold_op_sig_x
    autocorr_x = prefactor * np.trace(product_x)

    return autocorr_x

def autocorrelator_xyz_compare(time, size, hamiltonian):  # COMPARES SIG1_X, SIG1_Y, SIG1_Z)
    list = []
    list_y = []
    list_z = []
    autocorr = []
    for a in range(size):
        list.append(np_I)
        list_y.append(np_I)
        list_z.append(np_I)

    list[0] = np_x # first site association with sig_x set
    list_y[0] = np_y
    list_z[0] = np_z
    time_evol_op = time_evol_op_at_time_t(hamiltonian, time)  # Finds time evol operator at time t

    non_time_evold_op_sig_x = tensor_product(list, size) #  SIGX_1(0)
    non_time_evold_op_sig_y = tensor_product(list_y, size) #  SIGY_1(0)
    non_time_evold_op_sig_z = tensor_product(list_z, size) #  SIGZ_1(0)
    time_evolved_sig_x = (np.conjugate(time_evol_op.T)) @ non_time_evold_op_sig_x @ (time_evol_op)  #  SIGX_1(time)
    time_evolved_sig_y = (np.conjugate(time_evol_op.T)) @ non_time_evold_op_sig_y @ (time_evol_op)  #  SIGY_1(time)
    time_evolved_sig_z = (np.conjugate(time_evol_op.T)) @ non_time_evold_op_sig_z @ (time_evol_op)  #  SIGZ_1(time)
    # print(time_evol_op.getH() * time_evol_op)
    prefactor = 1/(2**size)
    product_x = time_evolved_sig_x @ non_time_evold_op_sig_x
    product_y = time_evolved_sig_y @ non_time_evold_op_sig_y
    product_z = time_evolved_sig_z @ non_time_evold_op_sig_z

    autocorr_x = prefactor * np.trace(product_x)
    autocorr_y = prefactor * np.trace(product_y)
    autocorr_z = prefactor * np.trace(product_z)
    
    autocorr.append(autocorr_x)
    autocorr.append(autocorr_y)
    autocorr.append(autocorr_z)

    return autocorr
#print(autocorrelator(1, L, H_1).real)  # Add .real (complex number shown is negligible, but
                                       # to make graphs we must ignore negligible complex nums)

def list_inputs(init_time, fin_time, interval):
    list = []
    num_of_runs = int((fin_time-init_time)/interval)
    for i in range(num_of_runs + 1):
        list.append(init_time + i*interval)


    return list

def output_list(input_list, size, hamiltonian):# TOOK OUT .REAL
    out_list = [0] * len(input_list)
    for i in range(len(input_list)):
        out_list[i] = autocorrelator(input_list[i], size, hamiltonian).real  #  Stores respective outputs in list
        #out_list[i] = autocorrelator(input_list[i], size, hamiltonian) --> DOES NOT FORCE REAL...
    return np.array(out_list)

def output_list_xyz_comp(input_list, size, hamiltonian):
    out_list_x = [0] * len(input_list)
    out_list_y = [0] * len(input_list)
    out_list_z = [0] * len(input_list)
    out_list = []
    for i in range(len(input_list)):
        auto_val_array = autocorrelator_xyz_compare(input_list[i], size, hamiltonian)
        out_list_x[i] = auto_val_array[0].real  #  Stores respective outputs in list
        out_list_y[i] = auto_val_array[1].real
        out_list_z[i] = auto_val_array[2].real
    out_list.append(np.array(out_list_x))
    out_list.append(np.array(out_list_y))
    out_list.append(np.array(out_list_z))
    return out_list

def simulation(init_time, fin_time, interval, size_array, hamiltonian):  #  COMPARES DIFFERENT SIZES

    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    return_list = []
    for k in range(len(hamiltonian)):
        return_list.append(output_list(in_list, size_array[k], hamiltonian[k]))

    in_list = np.array(in_list)  # LOG values of time
    
    plt.plot(in_list, return_list[0], marker = ".", linestyle = '-', color = 'b', label = '.2, .2, 1')
    plt.plot(in_list, return_list[1], marker = ".", linestyle = '-', color = 'g', label = '.4, .4, 1')
    plt.plot(in_list, return_list[2], marker = ".", linestyle = '-', color = 'r', label = '.5, .5, 1')
    plt.plot(in_list, return_list[3], marker = ".", linestyle = '-', color = 'm', label = '.6, .6, 1')
    plt.plot(in_list, return_list[4], marker = ".", linestyle = '-', color = 'y', label = '.7, .7, 1')

    plt.legend()
    plt.xscale('log')

    plt.xlabel('Time(s) - Log')
    plt.ylabel('Autocorrelator')
    plt.title('Test Autocorrelator graph')

    # plt.grid(True)
    plt.show()

def simulation_EXCEL_PARAM_COMPARISON(init_time, fin_time, interval, size_array, hamiltonian, array_params): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    return_list = [in_list]
    for k in range(len(hamiltonian)):
        return_list.append(np.ndarray.tolist(output_list(in_list, size_array[0], hamiltonian[k])))  # LIST W/ OUTPUTS OF DIFF PARAMS

    df = pd.DataFrame(return_list).T  #  TRANSPOSES ROWS TO MAKE COLUMNS WHICH CORRESPOND TO INPUTS AND OUTPUTS
    column_names = ['INPUT']
    for k in range(len(array_params)):
        input_params = array_params[k]
        input_name = ','.join(str(x) for x in input_params)
        name = "OUTPUT: " + input_name
        column_names.append(name)

    df.columns = column_names
    df.to_excel("C:/Users/frsin/OneDrive/Desktop/Grad school stuff/Research_Notes/Kolodrubetz Group/Data/CSV stuff/EXCEL/Mitra Paper/MY_OUTPUTS_MITRA_PAPER.xlsx") 

    print('SUCCESS')

def simulation_xyz_comp(init_time, fin_time, interval, size, hamiltonian): # Compares SIG_X1, SIG_Y1, SIG_Z1
    in_list = list_inputs(init_time, fin_time, interval)
    outputs = output_list_xyz_comp(in_list, size, hamiltonian)

    in_list = np.array(in_list)  # LOG values of time
    
    plt.plot(in_list, outputs[0], marker = ".", linestyle = '-', color = 'b', label = 'SIG_X')
    plt.plot(in_list, outputs[1], marker = ".", linestyle = '-', color = 'g', label = 'SIG_Y')
    plt.plot(in_list, outputs[2], marker = ".", linestyle = '-', color = 'r', label = 'SIG_Z')
    #plt.plot(in_list, return_list[3], marker = ".", linestyle = '-', color = 'm', label = 'L = 6')
    #plt.plot(in_list, return_list[4], marker = ".", linestyle = '-', color = 'y', label = 'L = 7')

    plt.legend()
    plt.xscale('log')

    plt.xlabel('Time(s) - Log')
    plt.ylabel('Autocorrelator')
    Title = "Test Autocorrelator graph: L = " + str(size)
    plt.title(Title)

    # plt.grid(True)
    plt.show()

def sim_autocorrelator_comparison_varying_params_GRAPHS(init_time, fin_time, interval, size, hamiltonian_array, array_params): # Makes graphs of varied params

    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    return_list = []    
    for k in range(len(hamiltonian_array)):
        return_list.append(output_list(in_list, size, hamiltonian_array[k]))

    in_list = np.array(in_list)
    #fig, axes = plt.subplot(len(hamiltonian_array), 1, figsize = (10, 4))

    for i in range(len(hamiltonian_array)):
        title = str(array_params[i][0]) + ", " + str(array_params[i][1]) + ", " + str(array_params[i][2])
        plt.subplot(int(len(hamiltonian_array)), 1, i+1)
        plt.plot(in_list, return_list[i])
        plt.title(title)
        plt.legend("L = " + str(size))
        plt.xscale('log')

    plt.tight_layout()
    plt.show()

#simulation_xyz_comp(0, 1000, .1, L, H_1)
#simulation(0, 500, .1, L_tot, H_tot)
#simulation_EXCEL_PARAM_COMPARISON(0, 500, .1, L_tot, H_tot, parameter_array)
#sim_autocorrelator_comparison_varying_params_GRAPHS(0, 1000, .1, L, [H_1], parameter_array)


#  KRYLOV METHOD

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
    I = np.eye(2**size)
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
    norm = cm.sqrt((1/(2**size))*norm_squared)

    #return np.linalg.norm(norm)
    return norm.real

def recursion_method(n, H, size):
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_x

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
            b_1 = norm(A_1, A_1, size)
            b.append(b_1)
            #print(type(A_1))
            O_1 = (1/(b_1))*A_1
            krylov_basis.append(O_1)

        else:
            A_basis.append(A_vectors(L, b, krylov_basis, i))
            b.append(norm(A_vectors(L, b, krylov_basis, i), A_vectors(L, b, krylov_basis, i), size))
            krylov_basis.append((1/b[i])*A_basis[i])

    print(np.array(b))
    return np.array(b)

#(recursion_method(50, H_1, L))

# def commutator(A, B):
#     first = A @ B
#     second = B @ A

#     return (first + (-1*second))


# def recursion_method_new(n, H, size):
#     list = []
#     for i in range(size):
#         list.append(np_I)
#     list[0] = np_x
#     strong_mode = tensor_product(list, size)
#     b = []
#     A_basis = []
#     O_0 = ((strong_mode.T).flatten()).tolist()
#     krylov_basis = [np.array(O_0[0])]

#     for j in range(0, n):
#         if (j == 0):
#             O_0_matrix_T = krylov_basis[0].reshape(strong_mode.shape)
#             O_0_matrix = O_0_matrix_T.T
#             A_1_matrix_T = (commutator(H, O_0_matrix)).T
#             A_1 = A_1_matrix_T.flatten().tolist()
#             A = np.array(A_1[0])
#             #print(A)
#             A_basis.append(A)
#             b_1 = norm(A, A, size).real
#             b.append(b_1)
#             O_1 = (1/b_1)*A
#             krylov_basis.append(O_1)
#             print(b)

#         else:
#             O_matrix_T = krylov_basis[j].reshape(strong_mode.shape)
#             O_matrix = O_matrix_T.T
#             A_matrix_T = (commutator(H, O_matrix)).T
#             A_basis.append(np.array(A_matrix_T.flatten().tolist()[0]) + (-1*b[j-1]*krylov_basis[j-1]))
#             b.append(norm(A_basis[j], A_basis[j], size))
#             print(b)
#             krylov_basis.append((1/b[j])*A_basis[j])

#         print(b)
#     return b




# in_list = []
# for i in range(1, 151):
#     in_list.append(i)

# plt.plot(np.array(in_list), recursion_method(150, H_1, L), marker = ".", linestyle = '-', color = 'b', label = '.2, .2, 1')
# plt.xscale('log')
# plt.ylim(0, 10)
# plt.yticks(np.arange(0, 11, 2))

# plt.show()



# GRAHM-SCHMIDT PROCEDURE

def grahm_schmidt(n, H, size):
    L = np.array((Liovillan(H, size)).tolist())
    list = []
    for i in range(size):
        list.append(np_I)
    list[0] = np_x
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
    b = []  # B VALUES
    for p in range(0, len(U) - 1):
        k = L @ U[p+1]
        b.append(U[p] @ k)

    return np.array(b)

# GRAPH TESTING
# in_list = []
# for i in range(1, 51):
#     in_list.append(i)

# plt.plot(np.array(in_list), recursion_method(50, H_1, L), marker = ".", linestyle = '-', color = 'r', label = '.2, .2, 1')
# plt.plot(np.array(in_list), grahm_schmidt(51, H_1, L), marker = ".", linestyle = '-', color = 'b', label = '.2, .2, 1')
# plt.xscale('log')
# plt.ylim(0, 10)
# plt.yticks(np.arange(0, 11, 2))

# plt.show()

def to_excel_recursion(n, H, size):
    in_list = []
    for i in range(1, n+1):
        in_list.append(i)
    
    b_values = recursion_method(n, H, size)
    return_list = [in_list, b_values]
    df = pd.DataFrame(return_list).T
    column_names = ["ITERATION NUMBER", "LACZOS COEFFICIENTS"]
    df.columns = column_names
    df.to_excel("C:/Users/frsin/OneDrive/Desktop/Grad school stuff/Research_Notes/Kolodrubetz Group/Data/CSV stuff/EXCEL/Mitra Paper/LACZOS/LACZOS_COEFF_FOR_HAM.xlsx") 
    print("SUCCESS")

    plt.plot(np.array(in_list), b_values, marker = ".", linestyle = '-', color = 'r', label = '.2, .2, 1')
    #plt.xscale('log')
    plt.ylim(0, 10)
    plt.yticks(np.arange(0, 11, 2))
    plt.show()


to_excel_recursion(50, H_1, L)


# b_rec = recursion_method(50, H_1, L).tolist()
# b_gs = grahm_schmidt(51, H_1, L).tolist()
# arr_rec = np.array(b_rec[:40])
# arr_gs = np.array(b_gs[:40])

# print(np.max(np.abs(arr_gs-arr_rec)))

