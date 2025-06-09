import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from Z3_code import term_one, term_two, recursion_method, grahm_schmidt




def LACZOS_COEF_vary_gam_SIM(n, size, gamma, params):
    #in_list = []
    #in_list_2 = []
    
    H_1 = term_one(size, Jay=1, PHI=params[1]) + term_two(F=params[0], length=size)
    H_dag = np.conjugate(H_1.T)
    H = H_1 + (gamma * H_dag)

    if gamma == 1:
        # b_values = recursion_method(n, H, size)
        #for i in range(1, n+1):
        #     in_list.append(i)
        matrix_final = recursion_method(n, H, size)
    else:
        # b_values = grahm_schmidt(n, H, size)
        #for i in range(1, n):
        #     in_list.append(i)
        matrix_final = grahm_schmidt(n, H, size)

    
    #return_list = [gamma, in_list, b_values]  # FIRST COMPONENT IS GAMMA VALUE

    #now = datetime.datetime.now()
    #npy_name = "vary_gamma_laczos_coeff" + now.strftime("%Y-%m-%d") + "_NPY.npy"
    #np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs
    npy_name = "krylov_matrix_gamma" + str(gamma) + "_NPY.npy"
    np.save(npy_name, matrix_final)

    plt.subplot()
    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " krylov hamiltonian (L = " + str(size) + ")"
    #plt.plot(np.array(in_list), np.array(b_values), marker = ".", linestyle = '-', color = 'r', label = str(gamma))
    #plt.xscale('log')
    #plt.plot(np.array(in_list), np.array(b_values), marker = ".", linestyle = '-', color = 'r', label = 'Recursive')
    #plt.plot(np.array(in_list_2), np.array(b_values_2), marker = ".", linestyle = '-', color = 'g', label = 'GS')
    F = matrix_final.real
    plt.imshow(F, cmap='seismic', aspect='auto')
    #plt.clim(-np.max(np.abs(F)), np.max(np.abs(F)))
    #plt.clim(-0.1, 0.1)
    plt.colorbar(label='Value')  # ← this adds the side bar
    plt.title(title)
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    plt.tight_layout()
    file_name = "gamma_" + str(gamma) + "_params_" + str(params[0])+"_"+str(params[1])+"_matrix.png"
    plt.savefig(file_name, dpi=300, bbox_inches='tight')


    # NOTE: GRAPH OF EIGENVALUES


    if gamma == 1:
        # b_values = recursion_method(n, H, size)
        A, U = np.linalg.eigh(matrix_final)
        EIG_STATES = []
        in_list = np.arange(0, n+1, 1)
        # EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        sorted = np.sort(np.abs(A))
        # EIG_STATE_2 = np.abs(U[:, np.abs(A).argsort()[1]])  # Second smallest eigenvalue
        # EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
        # EIG_STATE_2 = (1/np.linalg.norm(EIG_STATE_2)) * EIG_STATE_2
        for i in range(len(A)):
            EIG_STATE = np.abs(U[:, np.abs(A).argsort()[i]])
            EIG_STATE = (1/np.linalg.norm(EIG_STATE)) * EIG_STATE
            EIG_STATES.append(EIG_STATE)

    else:
        A, U = np.linalg.eig(matrix_final)
        #print(np.imag(A))
        in_list = np.arange(0, n+1, 1)
        EIG_STATES = []
        # EIG_STATE = np.abs(U[:, np.argmin(np.abs(A))])
        sorted = np.sort(np.abs(A))
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


    #plt.plot(in_list, EIG_STATE, marker = ".", linestyle = '-', color = 'r', label = 'E = ' + str(np.min(np.abs(A))))
    #plt.plot(in_list, EIG_STATE_2, marker = ".", linestyle = '-', color = 'g', label = 'E = ' + str(sorted[69]))
    ncols = 4  # columns per row
    nrows = int(np.ceil(len(EIG_STATES) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows))
    axes = axes.flatten()
    plt.suptitle(title)
    for i in range(int(len(EIG_STATES))):
    # #     #plt.figure()
         ax = axes[i]
         ax.plot(in_list, EIG_STATES[i], marker = ".", linestyle = '-', label = 'E = ' + str(np.round(sorted[i], 3)) + " (" + str(i) + ")")
         title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Eigenvector graph (L = " + str(size) + ")"

         ax.set_xlabel("n")
         ax.set_ylabel("Eigenvector Components")

         ax.set_ylim(0, 0.1)
         ax.legend()
    plt.savefig("eigenval_graphs_krylov_matrix.png", dpi=300, bbox_inches='tight')

    #print(EIG_STATES[0][3])
    # #plt.plot(in_list, EIG_STATES[62], marker = ".", linestyle = '-', label = 'E = ' + str(np.round(sorted[21], 3)) + " (" + str(11) + ")")

'''
    plt.subplot()
    plt.plot(in_list, EIG_STATES[0], marker = ".", linestyle = '-', color = 'r', label = 'E = ' + str(np.min(np.abs(A))))
    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma) + " Eigenvector graph (L = " + str(size) + ")"
    plt.suptitle(title)
    
    plt.ylim(0, 0.1)
    #plt.ylim(0, 0.01)
    plt.xlabel("n")
    plt.ylabel("Eigenvector Components")
    plt.tight_layout()
    plt.savefig("eigenval_graphs_krylov_matrix.png", dpi=300, bbox_inches='tight')
'''





    #plt.figure()
    #plt.plot(in_list, np.sort(Abs_A), marker = ".")
    #plt.tight_layout()
    #plt.show()