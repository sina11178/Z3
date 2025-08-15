from Z3_code import term_one, term_two, list_inputs, output_list
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime


def decay_rates(init_time, fin_time, interval, size, params, gamma_array):
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    in_list_exp = np.ndarray.tolist(in_list[1])
    return_list = [in_list_exp]
    hamiltonians = []
    A_list = []
    U_list = []
    U_inv_list = []
    for j in range(len(gamma_array)):
        H_1 = term_one(size, Jay=1, PHI=params[1]) + term_two(F=params[0], length=size)
        if gamma_array[j] == 0:
            H = H_1
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size))

            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        elif gamma_array[j] == 1:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + H_dag
            A, U = np.linalg.eigh(H.toarray())  # A is vector of eigenvalues; U is matrix of eigenvectors
            U_inv = np.conjugate(U).T  # Conjugate transpose of U

        else:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + (gamma_array[j] * H_dag)
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size))
            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        hamiltonians.append(H)
        A_list.append(A)
        U_list.append(U)
        U_inv_list.append(U_inv)

    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist(output_list(A_list[k], U_list[k], U_inv_list[k], in_list_exp, size, hamiltonians[k], gamma_array[k])))  # LIST W/ OUTPUTS OF DIFF PARAMS

    decay_rates = []

    for i in range(1, len(return_list)):
        decay_rates.append(-1*np.polyfit(return_list[0], np.log(np.abs(return_list[i])), 1)[0])

    plt.figure()
    plt.plot(gamma_array, decay_rates, marker='o')
    plt.xlabel('Gamma')
    plt.ylabel('Decay Rate')
    plt.title('Decay Rates vs Gamma')
    plt.grid()
    plt.show()