import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from Z3_code import term_one, term_two, recursion_method, grahm_schmidt
import os 




def LACZOS_COEF_vary_gam_SIM(n, size, gamma, params):
    #in_list = []
    #in_list_2 = []
    
    H_1 = term_one(size, Jay=1, PHI=params[1]) + term_two(F=params[0], length=size)
    H_dag = np.conjugate(H_1.T)
    H = H_1 + (gamma * H_dag)

    if gamma == 1:

        list_matrix_final = recursion_method(n, H, size)
        matrix_final = list_matrix_final[0]
        matrix_final_squared = list_matrix_final[1]
    else:

        list_matrix_final = grahm_schmidt(n, H, size)
        matrix_final = list_matrix_final[0]
        matrix_final_squared = list_matrix_final[1]
        truncated_matrix_size = (3*size) - 2
        matrix_final_trunc = matrix_final[:truncated_matrix_size, :truncated_matrix_size]

    os.makedirs('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_'+ str(size) + '_gamm_' + str(gamma), exist_ok=True)

    npy_name = "krylov_matrix/f_" + str(params[0]) + "_phi_" + str(params[1]) + "_L_"+ str(size) + "_gamm_" + str(gamma) + "/alpha_1_NPY.npy"
    np.save(npy_name, matrix_final)
    npy_name_2 = "krylov_matrix/f_" + str(params[0]) + "_phi_" + str(params[1]) + "_L_"+ str(size) + "_gamm_" + str(gamma) + "_alpha_1_squared_NPY.npy"
    np.save(npy_name_2, matrix_final_squared)

  
    if gamma != 0:
        F = np.abs(matrix_final)
    else:
        F = np.abs(matrix_final_trunc)

    F_squared = np.abs(matrix_final_squared)
    plt.figure()
    plt.imshow(F, cmap='seismic', aspect='auto')
    #plt.clim(-np.max(np.abs(F)), np.max(np.abs(F)))
    #plt.clim(-0.1, 0.1)
    plt.colorbar(label='Value')  # ← this adds the side bar
    plt.title('alpha_1 (abs value of matrix)')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/krylov_matrix_alpha_1.pdf')
    #plt.show()

    plt.figure()
    plt.imshow(F_squared, cmap='seismic', aspect='auto')
    plt.colorbar(label='Value')
    #plt.clim(0, 1.8)
    plt.title('alpha_1 squared (abs value of matrix)')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/krylov_matrix_alpha_1_squared.pdf')
    #plt.show()


    # NOTE: GRAPH OF EIGENVALUES


    if gamma != 0:
        # b_values = recursion_method(n, H, size)
        if gamma == 1:
            A, U = np.linalg.eigh(matrix_final)
            A_squared, U_squared = np.linalg.eigh(matrix_final_squared)
        else:
            A, U = np.linalg.eig(matrix_final)
            A_squared, U_squared = np.linalg.eig(matrix_final_squared)

        EIG_STATES = []
        EIG_STATES_SQUARED = []
        in_list = np.arange(0, n+1, 1)
        # EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        sorted = np.sort(np.abs(A))
        sorted_squared = np.sort(np.abs(A_squared))
        # EIG_STATE_2 = np.abs(U[:, np.abs(A).argsort()[1]])  # Second smallest eigenvalue
        # EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
        # EIG_STATE_2 = (1/np.linalg.norm(EIG_STATE_2)) * EIG_STATE_2
        for i in range(len(A)):
            EIG_STATE = np.abs(U[:, np.abs(A).argsort()[i]])
            EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
            EIG_STATES.append(EIG_STATE)

            EIG_STATE_SQUARED = np.abs(U_squared[:, np.abs(A_squared).argsort()[i]])
            EIG_STATE_SQUARED = (1/np.linalg.norm(EIG_STATE_SQUARED)) * EIG_STATE_SQUARED
            EIG_STATES_SQUARED.append(EIG_STATE_SQUARED)

    else:
        A, U = np.linalg.eig(matrix_final_trunc)  # NOTE: WE ONLY LOOK AT TRUNCATED MATRIX FOR FULLY NON-HERMITIAN
        A_squared, U_squared = np.linalg.eig(matrix_final_squared)  
        #print(np.imag(A))
        in_list = np.arange(0, n+1, 1)
        in_list_trunc = np.arange(0, truncated_matrix_size, 1)
        EIG_STATES = []
        EIG_STATES_SQUARED = []
        # EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        sorted = np.sort(np.abs(A))
        sorted_squared = np.sort(np.abs(A_squared))
        print(np.min(np.abs(A)))
        # EIG_STATE_2 = np.abs(U[:, np.abs(A).argsort()[1]])
        # EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
        # EIG_STATE_2 = (1/np.linalg.norm(EIG_STATE_2)) * EIG_STATE_2
        # Abs_A = np.abs(A)
        for i in range(len(A)):
            EIG_STATE = (U[:, np.abs(A).argsort()[i]])
            EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
            EIG_STATE = np.abs(EIG_STATE)
            #print(np.linalg.norm(EIG_STATE))
            EIG_STATES.append(EIG_STATE)

            EIG_STATE_SQUARED = (U_squared[:, np.abs(A_squared).argsort()[i]])
            EIG_STATE_SQUARED = (1/np.linalg.norm(EIG_STATE_SQUARED)) * EIG_STATE_SQUARED
            EIG_STATE_SQUARED = np.abs(EIG_STATE_SQUARED)
            EIG_STATES_SQUARED.append(EIG_STATE_SQUARED)


    #plt.subplot(2, 1, 1)
    if gamma != 0:
        plt.figure()
        plt.plot(in_list, EIG_STATES[0], marker = ".", linestyle = '-', color = 'r', label = 'a_1 E = ' + str(np.min(np.abs(A))))
        #title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Zero Eigenvector graphs (L = " + str(size) + ")"
        #plt.suptitle(title)
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.2)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1.pdf')
        plt.close()
    else:
        plt.figure()
        plt.plot(in_list_trunc, EIG_STATES[0], marker = ".", linestyle = '-', color = 'r', label = 'a_1 E = ' + str(np.min(np.abs(A))))
        #title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Zero Eigenvector graphs (L = " + str(size) + ")"
        #plt.suptitle(title)
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1 (NOTE: TRUNCATED)')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.2)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1.pdf')
        plt.close()       


    #plt.subplot(2, 1, 2)
    
    plt.figure()
    plt.plot(in_list, EIG_STATES_SQUARED[0], marker = ".", linestyle = '-', color = 'g', label = 'a_1^2 E = ' + str(np.min(np.abs(A_squared))))
    plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1 squared')
    plt.xlabel("n")
    plt.ylabel("Eigenvector Components")
    plt.ylim(0, 0.2)
    plt.legend()
    plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1_squared.pdf')

    # ncols = 4  # columns per row
    # nrows = int(np.ceil(len(EIG_STATES) / ncols))
    # nrows_squared = int(np.ceil(len(EIG_STATES_SQUARED) / ncols))  
    # fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows))
    # fig_squared, axes_squared = plt.subplots(nrows=nrows_squared, ncols=ncols, figsize=(4*ncols, 3*nrows_squared))
    # axes = axes.flatten()
    # axes_squared = axes_squared.flatten()
    # plt.suptitle(title)
    # for i in range(int(len(EIG_STATES))):
    # # #     #plt.figure()
    #      ax = axes[i]
    #      #ax_squared = axes_squared[i]
    #      ax.plot(in_list, EIG_STATES[i], marker = ".", linestyle = '-', label = 'E = ' + str(np.round(sorted[i], 3)) + " (" + str(i) + ")")
    #      #ax_squared.plot(in_list, EIG_STATES_SQUARED[i], marker = ".", linestyle = '-', label = 'E = ' + str(np.round(sorted_squared[i], 3)) + " (" + str(i) + ")")
    #      title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Eigenvector graph (L = " + str(size) + ")"

    #      ax.set_xlabel("n")
    #      #ax_squared.set_xlabel("n")
    #      ax.set_ylabel("Eigenvector Components")
    #      #ax_squared.set_ylabel("Eigenvector Components")


    #      ax.set_ylim(0, 0.2)
    #      #ax_squared.set_ylim(0, 0.2)
    #      #ax.set_yticks(np.arange(0, 0.21, 0.025))
    #      ax.legend()
    #      #ax_squared.legend()
    # plt.savefig("eigenval_graphs_alpha_1.pdf")
