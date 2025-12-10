import numpy as np
import matplotlib.pyplot as plt


parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 6, 8]
gamma = [0, 0.4, 0.8, 1]

input = np.arange(0, 30, .5)

fig, axes = plt.subplots(
    nrows=2, ncols=4, 
    figsize=(18, 8),
    sharex=False, sharey=False
)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
       ax = axes[i, j]
       for L in size:
            data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')
            ax.plot(np.exp(input), data[1], label=f'L = {L}')
            ax.set_xscale('log')

       if i == 0 and j == 0:
            ax.legend()
       if i == 1:
            ax.set_xlabel('t - log scale')
       if j == 0:
            ax.set_ylabel(f'$A_∞(t)$')
       if i == 0:
            ax.set_title(f'$Γ={gam}$')


fig.text(0.03, 0.75, rotation =90, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
         va='center', ha='left', fontsize=14)

fig.text(0.03, 0.27, rotation =90, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
         va='center', ha='left', fontsize=14)
fig.text(0.05, .9, '(a)', fontsize=18)
fig.text(0.05, 0.45, '(b)', fontsize=18)

plt.savefig("my_figure.pdf", format="pdf")
