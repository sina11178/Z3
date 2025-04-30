import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from Z3_code import term_one, term_two, recursion_method, grahm_schmidt




def LACZOS_COEF_vary_gam_SIM(n, size, gamma, params):
    in_list = []
    in_list_2 = []
    
    H_1 = term_one(size, Jay=1, PHI=params[1]) + term_two(F=params[0], length=size)
    H_dag = np.conjugate(H_1.T)
    H = H_1 + (gamma * H_dag)

    if gamma == 1:
        b_values = recursion_method(n, H, size)
        for i in range(1, n+1):
            in_list.append(i)
    else:
        b_values = grahm_schmidt(n, H, size)
        for i in range(1, n):
            in_list.append(i)
    
    return_list = [gamma, in_list, b_values]  # FIRST COMPONENT IS GAMMA VALUE

    now = datetime.datetime.now()
    npy_name = "vary_gamma_laczos_coeff" + now.strftime("%Y-%m-%d") + "_NPY.npy"
    np.save(npy_name, np.array(return_list))  #  NOTE: First row of array is input; the rest are outputs

    title = str(params[0]) + ", " + str(params[1]) + "  " + "(f, PHI)" + " gam = " + str(gamma)
    plt.plot(np.array(in_list), np.array(b_values), marker = ".", linestyle = '-', color = 'r', label = str(gamma))
    #plt.xscale('log')
    #plt.plot(np.array(in_list), np.array(b_values), marker = ".", linestyle = '-', color = 'r', label = 'Recursive')
    #plt.plot(np.array(in_list_2), np.array(b_values_2), marker = ".", linestyle = '-', color = 'g', label = 'GS')


    plt.title(title)
    plt.xlabel("n")
    plt.ylabel("b_n")
    #plt.ylim(-0.05, 0)
    plt.legend()
    plt.tight_layout()
    plt.savefig("vary_gamma_params.png", dpi=300, bbox_inches='tight')