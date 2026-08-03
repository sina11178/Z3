import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 5, 6, 7, 8]
gamma = [0, 0.1, 0.5, 0.9, 1]
input = np.arange(0, 30, .5)


colors = ['#1f77b4', '#ff7f0e', '#2ca02c', "#5107d0", "#be008c"]  # colorblind-friendly
labels = ['a.)', 'b.)', 'c.)', 'd.)', 'e.)']
row_labels = [
    r'$f = 0.1, \phi = 0$',
    r'$f = 0.1, \phi = \pi/6$'
]

fig, axes = plt.subplots(
    nrows=2, ncols=5,
    figsize=(18, 8),
    sharex=False, sharey=False
)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
        ax = axes[i, j]
        for k, L in enumerate(size):
            data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')
            ax.plot(np.exp(input), data[1], label=f'$L = {L}$', color=colors[k], linewidth=1.5)
            ax.set_xscale('log')
        #ax.set_ylim(-0.5, 1.1)
        #ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
        if i == 0 and j == 0:
            ax.legend(fontsize=9,
        loc='upper right',
        frameon=True,
        handlelength=1.2,
        handletextpad=0.4,
        borderpad=0.3,
        labelspacing=0.2)
        if i == 1 and j == 0:
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


plt.savefig("my_figure_V2.pdf", format="pdf", bbox_inches='tight')
plt.savefig("my_figure_V2.png", dpi=600, bbox_inches='tight')
plt.show()