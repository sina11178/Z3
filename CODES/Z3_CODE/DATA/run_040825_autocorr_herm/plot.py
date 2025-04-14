import matplotlib.pylab as plt
import numpy as np 
import os

size_list = [4, 5, 6, 7, 8]

for i in np.arange(16):
    if os.path.isfile('./b'+str(i)+'/size_autocorr_comparison_2025-04-13_NPY.npy') == True:
        dat = np.load('./b'+str(i)+'/size_autocorr_comparison_2025-04-13_NPY.npy')
        fig, ax = plt.subplots()
        plt.suptitle("Paramater" + str(i))
        time_list = dat[0, :]
        for j in range(1, dat.shape[0]):
            plt.plot(np.log(time_list), dat[j, :], label= "L = " + str(size_list[j-1]))
        plt.legend()
        #plt.savefig("parameter"+str(i)+".pdf") # --> FOR PUBLICATION, USE PDF
        plt.savefig("parameter"+str(i)+".png")

