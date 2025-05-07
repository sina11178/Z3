'''

SIMULATION CODE FOR COMPARING AUTOCORRELATORS OF HAMILTONIANS WITH DIFFERENT SIZES (SAME PARAMETERS)

'''

from Z3_code import term_one, term_two, list_inputs, output_list
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


def size_autocorrelator_SIMULATION_NPY(init_time, fin_time, interval, size_array, params, gamma): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
   # in_list = list_inputs(init_time, fin_time, interval)
    #out_list = [0] * len(in_list)
    k = np.array(time_sequence(42, 15, 2.6))
    in_list = [np.log(k), k]
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

    
    #  Name EXCEL based on date
    now = datetime.datetime.now()
    npy_name = "size_autocorr_comparison_gamma_" + str(gamma) + "_NPY.npy"
    np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs

    print('SUCCESS')

    #  FOR GRAPHING
    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma)
    for i in range(len(hamiltonians)):
        #plt.subplot(1, int(len(hamiltonians)), i+1)
        plt.plot(in_list[0], return_list[i+1], label= "L = " + str(size_array[i]))
        plt.title(title)
        plt.xlabel("Time: ln(t)")
        plt.ylabel("Autocorr.")
        plt.legend()
        #plt.xscale('log')

    plt.tight_layout()
    plt.savefig("vary_params.png", dpi=300, bbox_inches='tight')
    #plt.show()
