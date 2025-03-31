from Z3_code import term_one, term_two, list_inputs, output_list
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd



def f_phi_autocorrelator_SIMULATION(init_time, fin_time, interval, size_array, hamiltonian, array_params): # COMPARES SYSTEMS W DIFFERENT PARAMS AT SAME SIZE
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
