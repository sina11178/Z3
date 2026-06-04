'''

SIMULATION CODE FOR COMPARING AUTOCORRELATORS OF HAMILTONIANS WITH DIFFERENT SIZES (SAME PARAMETERS)

'''

from Z3_code import term_one, term_two, list_inputs, output_list, output_list_corr
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime


# Array Params --> (f, phi)
def size_autocorrelator_SIMULATION_EXCEL(init_time, fin_time, interval, size_array, params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    in_list_exp = np.ndarray.tolist(in_list[1])
    return_list = [in_list_exp]
    hamiltonians = []
    for j in range(len(size_array)):
        H_1 = term_one(size_array[j], Jay=1, PHI=params[1]) + term_two(F=params[0], length=size_array[j])
        H_dag = np.conjugate(H_1.T)
        H = H_1 + (gamma * H_dag)
        hamiltonians.append(H)

    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist(output_list(in_list_exp, size_array[k], hamiltonians[k], gamma)))  # LIST W/ OUTPUTS OF DIFF PARAMS

    df = pd.DataFrame(return_list).T  #  TRANSPOSES ROWS TO MAKE COLUMNS WHICH CORRESPOND TO INPUTS AND OUTPUTS
    column_names = ['INPUT']
    for k in range(len(size_array)):
        input_size = size_array[k]
        #input_name = ','.join(str(x) for x in input_params)
        name = "OUTPUT(size): (" + str(input_size) + ")"
        column_names.append(name)

    df.columns = column_names

    #  Name EXCEL based on date
    now = datetime.datetime.now()
    excel_name = "./CODES/Z3_CODE/DATA/size_autocorr_comparison_" + now.strftime("%Y-%m-%d") + ".xlsx"
    df.to_excel(excel_name) 

    print('SUCCESS')


def size_autocorrelator_SIMULATION_NPY(init_time, fin_time, interval, size_array, params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    in_list_exp = np.ndarray.tolist(in_list[1])
    return_list = [in_list_exp]
    hamiltonians = []
    A_list = []
    U_list = []
    U_inv_list = []
    for j in range(len(size_array)):
        H_1 = term_one(size_array[j], Jay=1, PHI=params[1]) + term_two(F=params[0], length=size_array[j])
        if gamma == 0:
            H = H_1
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size_array[j]))

            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        elif gamma == 1:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + H_dag
            A, U = np.linalg.eigh(H.toarray())  # A is vector of eigenvalues; U is matrix of eigenvectors
            U_inv = np.conjugate(U).T  # Conjugate transpose of U

        else:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + (gamma * H_dag)
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size_array[j]))
            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        hamiltonians.append(H)
        A_list.append(A)
        U_list.append(U)
        U_inv_list.append(U_inv)

    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist(output_list(A_list[k], U_list[k], U_inv_list[k], in_list_exp, size_array[k], hamiltonians[k], gamma)))  # LIST W/ OUTPUTS OF DIFF PARAMS

    
    #  Name EXCEL based on date
    #now = datetime.datetime.now()
    if len(size_array) == 1:
        npy_name = "autocorr_gamma_" + str(gamma) + "_L" + str(size_array[0]) +"_f_"+str(round(params[0],3))+"_phi_"+str(round(params[1],3))+"_NPY.npy" #NOTE: ONLY FOR ARRAY SIZE OF 1
    else:
        npy_name = "size_autocorr_comparison_gamma_" + str(gamma) +"_NPY.npy"
    np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs

    print('SUCCESS')

    #  FOR GRAPHING
#    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma)
#    for i in range(len(hamiltonians)):
#        #plt.subplot(1, int(len(hamiltonians)), i+1)
#        plt.plot(in_list[0], np.abs(return_list[i+1]), label= "L = " + str(size_array[i]))
#        plt.title(title)
#        plt.xlabel("Time: ln(t)")
#        plt.ylabel("Autocorr.")
#        plt.legend()
        #plt.xscale('log')

#    plt.tight_layout()
#    plt.savefig("vary_params.png", dpi=300, bbox_inches='tight')
    #plt.show()



def size_correlator_SIM_NPY(init_time, fin_time, interval, size_array, params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    in_list_exp = np.ndarray.tolist(in_list[1])
    return_list = [in_list_exp]
    hamiltonians = []
    A_list = []
    U_list = []
    U_inv_list = []
    for j in range(len(size_array)):
        H_1 = term_one(size_array[j], Jay=1, PHI=params[1]) + term_two(F=params[0], length=size_array[j])
        if gamma == 0:
            H = H_1
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size_array[j]))

            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        elif gamma == 1:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + H_dag
            A, U = np.linalg.eigh(H.toarray())  # A is vector of eigenvalues; U is matrix of eigenvectors
            U_inv = np.conjugate(U).T  # Conjugate transpose of U

        else:
            H_dag = np.conjugate(H_1.T)
            H = H_1 + (gamma * H_dag)
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size_array[j]))
            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        hamiltonians.append(H)
        A_list.append(A)
        U_list.append(U)
        U_inv_list.append(U_inv)

    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist(output_list_corr(A_list[k], U_list[k], U_inv_list[k], in_list_exp, size_array[k], hamiltonians[k], gamma)))  # LIST W/ OUTPUTS OF DIFF PARAMS

    
    #  Name EXCEL based on date
    now = datetime.datetime.now()
    npy_name = "size_autocorr_comparison_gamma_" + str(gamma) + "_NPY.npy"
    np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs

    print('SUCCESS')

    #  FOR GRAPHING
    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma)
    for i in range(len(hamiltonians)):
        #plt.subplot(1, int(len(hamiltonians)), i+1)
        plt.plot(in_list[0], np.abs(return_list[i+1]), label= "L = " + str(size_array[i]))
        plt.title(title)
        plt.xlabel("Time: ln(t)")
        plt.ylabel("Autocorr.")
        plt.legend()
        #plt.xscale('log')

    plt.tight_layout()
    plt.savefig("vary_params.png", dpi=300, bbox_inches='tight')
    #plt.show()
