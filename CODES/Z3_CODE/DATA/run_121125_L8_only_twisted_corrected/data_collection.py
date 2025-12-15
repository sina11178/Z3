import numpy as np
import os
import pandas as pd
import shutil



param_array=[[0.1, 0], [0.1, np.pi/6]]
#gamma_array = [0, .2, .4, .6, .8, 1]
#L_array = [4, 5, 6, 8]
gamma_array = [0, .4, .8, .9]

for i in range(8):
    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i%4])+'_twisted.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[int(i/4)][0]
        phi = round(param_array[int(i/4)][1], 3)
        
        destination = '../krylov_matrices/twisted/f_'+str(f)+'_phi_'+str(phi)+'/seed/'
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i%4])+'_twisted.npy'
        #f_n = './b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 4]) + '_gamma_'+str(gamma_array[int(int(i)/4)%6])+'_twisted.npy'
        shutil.copy(f_n, destination)


    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i%4])+'_twisted_squared.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[int(i/4)][0]
        phi = round(param_array[int(i/4)][1], 3)

        destination = '../krylov_matrices/twisted/f_'+str(f)+'_phi_'+str(phi)+'/seed_squared/'
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(8) + '_gamma_'+str(gamma_array[i%4])+'_twisted_squared.npy'

        shutil.copy(f_n, destination)


print("successful")
