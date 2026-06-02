import numpy as np
import os
import pandas as pd
import shutil



param_array=[[0.1, 0.524]]
#gamma_array = [0, .2, .4, .6, .8, 1]
#L_array = [4, 5, 6, 7, 8]
gamma_array = [0, .1, .5, .9, 1]

for i in range(5):
    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i])+'_f_'+str(param_array[0][0])+'_phi_'+str(round(param_array[0][1], 3))+'.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[0][0]
        phi = round(param_array[0][1], 3)
        
        destination = '../krylov_matrices_V2/normal/f_'+str(f)+'_phi_'+str(phi)+'/seed/'
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i])+'_f_'+str(f)+'_phi_'+str(phi)+'.npy'
        os.makedirs(destination, exist_ok=True)
        shutil.copy(f_n, destination)


    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i])+'_f_'+str(param_array[0][0])+'_phi_'+str(round(param_array[0][1], 3))+'_squared.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[0][0]
        phi = round(param_array[0][1], 3)

        destination = '../krylov_matrices_V2/normal/f_'+str(f)+'_phi_'+str(phi)+'/seed_squared/'
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i])+'_f_'+str(f)+'_phi_'+str(phi)+'_squared.npy'
        os.makedirs(destination, exist_ok=True)
        shutil.copy(f_n, destination)


print("successful")
