import matplotlib.ticker as ticker
import numpy as np
import matplotlib.pyplot as plt
from torch import lt


parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 5, 6, 7, 8]
gamma = [0, 0.1, 0.5, 0.9, 1]

input = np.arange(0, 30, .5)

#plt.figure(figsize=(5, 5))
'''
for i, param in enumerate(parameters):
     for j, gam in enumerate(gamma):
          for L in size:
               data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')[1]
               plt.plot(np.exp(input), data, label=f'L = {L}')
          plt.xscale('log')
          plt.ylim(-0.4, 1.1)
          #plt.xticks([])
          #plt.yticks([])
          if i == 0 and j == 0:
               plt.legend()
          #if j == 0 and i == 1:
          #     plt.ylabel(f'$A_∞(t)$')
          if j != 0:
               plt.yticks([])
          if i != 1:
               plt.xticks([])
          plt.savefig(f"f_{param[0]}_phi_{param[1]}_gamma_{gam}_size_{L}.pdf", format="pdf")
          plt.savefig(f"f_{param[0]}_phi_{param[1]}_gamma_{gam}_size_{L}.png", dpi = 600, bbox_inches='tight')
          plt.clf()
'''

fig, axes = plt.subplots(
    nrows=2, ncols=5, 
    figsize=(18, 7),
    sharex='col', sharey='row', constrained_layout=True
)

for i, param in enumerate(parameters):
     for j, gam in enumerate(gamma):
          ax = axes[i, j]
          for k, L in enumerate(size):
               data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')
               ax.plot(np.exp(input), data[1], label=f'L = {L}')
               ax.tick_params(axis='both', labelsize=13)
               plt.locator_params(axis='y', nbins=4)
               #plt.locator_params(axis='x', nbins=4)
               plt.setp(ax.spines.values(), linewidth=1.3)
               ax.set_xscale('log')
               ax.xaxis.set_major_locator(ticker.LogLocator(base=10, numticks=4))
               #ax.yaxis.set_major_locator(ticker.Locator(base=10, numticks=4))

               ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
               ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
               ax.tick_params(axis="x",direction="in", which= 'both')
               ax.tick_params(axis="y",direction="in", which= 'both')
               ax.tick_params(axis="x", which= 'major', length = 5.0, width = 1.4)
               ax.tick_params(axis="x", which= 'minor', length = 3.0, width = 1.4)
               ax.tick_params(axis="y", which= 'major', length = 5.0, width = 1.4)
               ax.tick_params(axis="y", which= 'minor', length = 3.0, width = 1.4)
               plt.xticks(fontsize=16)
               plt.yticks(fontsize=16)
            #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
            #ax.set_yticks([-.4, 0, .4, .8, 1])
          if i == 0 and j == 0:
               ax.legend(fontsize=13, framealpha=0.3)
       #if j != 0:
       #     ax.set_yticks([])
       #if i != 1:
       #     ax.set_xticks([])
       #if i == 1:
       #     ax.set_xlabel('t - log scale')
       #if j == 0:
       #     ax.set_ylabel(f'$A_∞(t)$')
       #if i == 0:
       #     ax.set_title(f'$Γ={gam}$')
       #ax.spines['top'].set_visible(False)
       #ax.spines['right'].set_visible(False)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# fig.text(0.03, 0.75, rotation =90, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
#          va='center', ha='left', fontsize=14)

# fig.text(0.03, 0.27, rotation =90, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
#          va='center', ha='left', fontsize=14)
# fig.text(0.05, .9, '(a)', fontsize=18)
# fig.text(0.05, 0.45, '(b)', fontsize=18)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
plt.savefig("my_figure.pdf", format="pdf")
plt.savefig("my_figure.png", dpi = 600, bbox_inches='tight')