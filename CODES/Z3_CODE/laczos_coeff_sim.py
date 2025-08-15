import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from Z3_code import term_one, term_two, recursion_method, grahm_schmidt
import os 




def LACZOS_COEF_vary_gam_SIM(n, size, gamma, params):

    H_1 = term_one(size, Jay=1, PHI=params[1]) + term_two(F=params[0], length=size)
    H_dag = np.conjugate(H_1.T)
    H = H_1 + (gamma * H_dag)

    if gamma == 1:

        list_matrix_final = recursion_method(n, H, size)
        matrix_final = list_matrix_final[0]
        matrix_final_squared = list_matrix_final[1]
        list_matrix_final_2 = grahm_schmidt(n, H, size)

        # print(np.max(np.abs((matrix_final)-(list_matrix_final_2[0]))))
        # arr_1 = []
        # arr_2 = []
        # for i in range(len(matrix_final)-1):
        #     arr_1.append(np.abs(matrix_final[i, i+1]))
        #     arr_2.append(np.abs((list_matrix_final_2[0][i, i+1])))
        # in_list = np.arange(0, n, 1)
        # plt.figure()
        # plt.plot(in_list, arr_1, marker = ".", linestyle = '-', color = 'r', label = 'a_1')
        # plt.plot(in_list, arr_2, marker = ".", linestyle = '-', color = 'b', label = 'a_1 (Graham Schmidt)')
        # #plt.ylim(7, 8.5, .05)
        # plt.legend()
        # plt.show()

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

    plt.imshow(F, cmap='seismic', aspect='auto')
    #plt.clim(-np.max(np.abs(F)), np.max(np.abs(F)))
    #plt.clim(-0.1, 0.1)
    plt.colorbar(label='Value')  # ← this adds the side bar
    plt.title('alpha_1 (abs value of matrix)')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/krylov_matrix_alpha_1.pdf')
    plt.show()


    plt.imshow(F_squared, cmap='seismic', aspect='auto')
    plt.colorbar(label='Value')
    #plt.clim(0, 1.8)
    plt.title('alpha_1 squared (abs value of matrix)')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/krylov_matrix_alpha_1_squared.pdf')
    plt.show()


    # NOTE: GRAPH OF EIGENVALUES


    if gamma != 0:
        # b_values = recursion_method(n, H, size)
        if gamma == 1:
            A, U = np.linalg.eigh(matrix_final)
            A_2, U_2 = np.linalg.eigh(list_matrix_final_2[0])
            A_squared, U_squared = np.linalg.eigh(matrix_final_squared)
            #print(np.min(np.abs(A)))
        else:
            A, U = np.linalg.eig(matrix_final)
            A_squared, U_squared = np.linalg.eig(matrix_final_squared)


        in_list = np.arange(0, n+1, 1)

        EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
        #EIG_STATE_2 = np.abs(U_2[:, np.argmin(np.abs(A_2))])
        #EIG_STATE_2 = (1/np.linalg.norm(EIG_STATE_2)) * EIG_STATE_2
        EIG_STATE_SQUARED = np.abs(U_squared[:, np.argmin(np.abs(A_squared))])
        EIG_STATE_SQUARED = (1/np.linalg.norm(EIG_STATE_SQUARED)) * EIG_STATE_SQUARED
        
        plt.figure()
        plt.plot(in_list, EIG_STATE, marker = ".", linestyle = '-', color = 'r', label = 'a_1 E = ' + str(np.min(np.abs(A))))
        #plt.plot(in_list, EIG_STATE_2, marker = ".", linestyle = '-', color = 'b', label = 'a_1 (Graham Schmidt) E = ' + str(np.min(np.abs(A_2))))
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.5)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1.pdf')
        plt.close()

        plt.figure()
        plt.plot(in_list, EIG_STATE_SQUARED, marker = ".", linestyle = '-', color = 'g', label = 'a_1^2 E = ' + str(np.min(np.abs(A_squared))))
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1 squared')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.5)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1_squared.pdf')

    else:
        A, U = np.linalg.eig(matrix_final_trunc)  # NOTE: WE ONLY LOOK AT TRUNCATED MATRIX FOR FULLY NON-HERMITIAN
        A_squared, U_squared = np.linalg.eig(matrix_final_squared)  
        #print(np.imag(A))
        in_list = np.arange(0, n+1, 1)
        in_list_trunc = np.arange(0, truncated_matrix_size, 1)


        EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
        EIG_STATE_SQUARED = np.abs(U_squared[:, np.argmin(np.abs(A_squared))])
        EIG_STATE_SQUARED = (1/np.linalg.norm(EIG_STATE_SQUARED)) * EIG_STATE_SQUARED
        
        plt.figure()
        plt.plot(in_list_trunc, EIG_STATE, marker = ".", linestyle = '-', color = 'r', label = 'a_1 E = ' + str(np.min(np.abs(A))))
        #title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Zero Eigenvector graphs (L = " + str(size) + ")"
        #plt.suptitle(title)
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1 (NOTE: TRUNCATED)')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.2)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1.pdf')
        plt.close()    

        plt.figure()
        plt.plot(in_list, EIG_STATE_SQUARED, marker = ".", linestyle = '-', color = 'g', label = 'a_1^2 E = ' + str(np.min(np.abs(A_squared))))
        plt.title('Eigenvector Components for Zero Eigenvalue: alpha_1 squared')
        plt.xlabel("n")
        plt.ylabel("Eigenvector Components")
        plt.ylim(0, 0.2)
        plt.legend()
        plt.savefig('krylov_matrix/f_' + str(params[0]) + '_phi_' + str(params[1]) + '_L_' + str(size) + '_gamm_' + str(gamma) + '/zero_eigenvector_alpha_1_squared.pdf')
    
