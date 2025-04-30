import matplotlib.pylab as plt
import numpy as np 
import os
import pandas as pd

array_param=[.1, 0]
gamma_list = [0, .2, .4, .6, .8, .9, 1, 1.2]
size_list = [4, 5, 6, 7]
data_list = []
only_hermitian_L8 = []
time_list = np.arange(0, 30, 1)
time_list_exp = np.exp(time_list)



for i in np.arange(32):
     if os.path.isfile('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy') == True:
            
        p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
        data_list.append(p[1, :])



for k in range(len(gamma_list)):
    fig, ax = plt.subplots()
    plt.suptitle("Paramater (f=.1, phi=0): gamma = " + str(gamma_list[k]))
    for j in range(4):
        plt.plot(time_list, data_list[len(size_list)*k + j], label= "L = " + str(size_list[j]))
        
    plt.legend()
    plt.savefig("parameter_gamma_"+str(gamma_list[k])+"_graph.png")


for j in range(32):
    data_list[j] = data_list[j].tolist()
data_list = np.concatenate(time_list_exp, np.array(data_list))
df = pd.DataFrame(data_list).T
df.to_excel('autocorr_data.xlsx', index=False, header=False)