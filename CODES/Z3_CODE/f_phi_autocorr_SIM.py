'''

SIMULATION CODE FOR COMPARING AUTOCORRELATORS OF HAMILTONIANS WITH DIFFERENT PARAMETERS

'''

from Z3_code import term_one, term_two, list_inputs, output_list
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime


# Array Params --> (f, phi)
def f_phi_autocorrelator_SIMULATION_EXCEL(init_time, fin_time, interval, size, array_params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    return_list = [np.ndarray.tolist(in_list[1])]
    hamiltonians = []
    for j in range(len(array_params)):
        H_1 = term_one(size, Jay=1, PHI=array_params[j][1]) + term_two(F=array_params[j][0], length=size)
        H_dag = np.conjugate(H_1.T)
        H = H_1 + (gamma * H_dag)
        hamiltonians.append(H)

    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist(output_list(np.ndarray.tolist(in_list[1]), size, hamiltonians[k], gamma)))  # LIST W/ OUTPUTS OF DIFF PARAMS

    df = pd.DataFrame(return_list).T  #  TRANSPOSES ROWS TO MAKE COLUMNS WHICH CORRESPOND TO INPUTS AND OUTPUTS
    column_names = ['INPUT']
    for k in range(len(array_params)):
        input_params = array_params[k]
        input_name = ','.join(str(x) for x in input_params)
        name = "OUTPUT(f, phi): (" + str(input_name) + ")"
        column_names.append(name)

    df.columns = column_names

    #  Name EXCEL based on date
    now = datetime.datetime.now()
    excel_name = "./CODES/Z3_CODE/DATA/f_phi_autocorr_comparison_" + now.strftime("%Y-%m-%d") + ".xlsx"
    df.to_excel(excel_name) 

    print('SUCCESS')

def time_sequence(myn, mys, mym):
    sequ=[]
    i=0
    for x in range(myn):
        if x < mys:
            i = (i + x)
        if x >= mys:
            i = (i + x)*mym
        sequ.append(int(i))
    return sequ





def f_phi_autocorrelator_SIMULATION_NPY(init_time, fin_time, interval, size, array_params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
    #in_list = list_inputs(init_time, fin_time, interval)  # GIVES LIST OF NUMPY ARRAYS
    #in_list = [np.log(np.array(time_sequence(42, 15, 2.6))), np.array(time_sequence(42, 15, 2.6))]
    in_list = list_inputs(init_time, fin_time, interval)

    in_list_exp = np.ndarray.tolist(in_list[1])
    return_list = [in_list_exp]
    hamiltonians = []
    A_list = []
    U_list = []
    U_inv_list = []
    for j in range(len(array_params)):
        '''
        H_1 = term_one(size, Jay=1, PHI=array_params[j][1]) + term_two(F=array_params[j][0], length=size)
        H_dag = np.conjugate(H_1.T)
        H = H_1 + (gamma * H_dag)
        hamiltonians.append(H)
        np.save('Sina_ModifiedBaxter.npy', H)
        '''
        H_1 = term_one(size, Jay=1, PHI=array_params[j][1]) + term_two(F=array_params[j][0], length=size)
        if gamma == 0:
            H = H_1
            ham = H.toarray()
            U, V = np.linalg.eig(ham)
            #print(hamiltonian.shape)
            #print((1.j*max(np.array(U).imag)*np.eye(3**L)).shape)
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size))

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
            Hpr = (ham - 1.j*max(np.array(U).imag)*np.eye(3**size))
            A, U = np.linalg.eig(Hpr)
            U_inv = np.linalg.inv(U)

        hamiltonians.append(H)
        A_list.append(A)
        U_list.append(U)
        U_inv_list.append(U_inv)


    for k in range(len(hamiltonians)):
        return_list.append(np.ndarray.tolist((output_list(A_list[k], U_list[k], U_inv_list[k], in_list_exp, size, hamiltonians[k], gamma))))  # LIST W/ OUTPUTS OF DIFF PARAMS

    
    #  Name EXCEL based on date
    now = datetime.datetime.now()
    npy_name = "f_phi_autocorr_comparison_" + now.strftime("%Y-%m-%d") + "_NPY.npy"  # SAVE IN SCRATCH
    np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs

    print('SUCCESS')

    #  FOR GRAPHING
    for i in range(len(hamiltonians)):  # GIVES GRAPH WITH LOG TIME VALUES
        title = "L= " + str(size) + "; gam = " + str(gamma)
        #plt.subplot(1, int(len(hamiltonians)), i+1)
        plt.plot(in_list[1], return_list[i+1], label=str(array_params[i][0]) + ", " + str(array_params[i][1])+ "  " + "(f, PHI)")
        plt.title(title)
        plt.xlabel("Time: ln(t)")
        plt.ylabel("Autocorr.")
        plt.legend()
        plt.xscale('log')

    plt.tight_layout()
    #plt.savefig("vary_params.png", dpi=300, bbox_inches='tight')
    plt.show()
    #plt.close()