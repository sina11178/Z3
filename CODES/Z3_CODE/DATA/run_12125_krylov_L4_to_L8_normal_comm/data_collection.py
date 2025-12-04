import numpy as np
import os
import pandas as pd
import shutil



param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [0, .2, .4, .6, .8, 1]
L_array = [4, 5, 6, 7, 8]


for i in range(60):
    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str( param_array[int(int(i)/30)][1])+'.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[int(i/30)][0]
        phi = round(param_array[int(i/30)][1], 3)
        
        destination = '../krylov_matrices/normal/krylov_matrices_120425_vary_gamma/f_'+str(f)+'_phi_'+str(phi)+'/seed/'
<<<<<<< HEAD
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str( param_array[int(int(i)/30)][1])+'.npy'
=======
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str(param_array[int(int(i)/30)][1])+'.npy'
>>>>>>> 149ed3def87fcc79f8f1613d8cb25716845f20d0
        shutil.copy(f_n, destination)


    if os.path.isfile('./b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str(param_array[int(int(i)/30)][1])+'_squared.npy') == True:

        #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        #data_list.append(p[1, :])
        f = param_array[int(i/30)][0]
        phi = round(param_array[int(i/30)][1], 3)

        destination = '../krylov_matrices/normal/krylov_matrices_120425_vary_gamma/f_'+str(f)+'_phi_'+str(phi)+'/seed_squared/'
<<<<<<< HEAD
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str( param_array[int(int(i)/30)][1])+'_squared.npy'
=======
        f_n = './b'+str(i)+'/krylov_matrix_L' + str(L_array[int(i) % 5]) + '_gamma_'+str(gamma_array[int(int(i)/5)%6])+'_f_'+str( param_array[int(int(i)/30)][0])+'_phi_'+str(param_array[int(int(i)/30)][1])+'_squared.npy'
>>>>>>> 149ed3def87fcc79f8f1613d8cb25716845f20d0
        shutil.copy(f_n, destination)


print("successful")
