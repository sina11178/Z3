import numpy as np
import os
import pandas as pd
import shutil



param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [0, .2, .4, .6, .8, 1]
L_array = [4, 6, 8, 9, 10]


for i in range(60):
     if os.path.isfile('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str( gamma_array[int(i/5)%6]) + '_size_'+str(L_array[i % 5])+'_NPY.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[int(i/30)][0]
        phi = round(param_array[int(i/30)][1], 3)
        
        destination = '../Autocorrelation_data/f_'+str(f)+'_phi_'+str(phi)+'/'
        f_n = './b'+str(i)+'/size_autocorr_comparison_gamma_' + str( gamma_array[int(i/5)%6]) + '_size_'+str(L_array[i % 5])+'_NPY.npy'
        shutil.copy(f_n, destination)


print("successful")
