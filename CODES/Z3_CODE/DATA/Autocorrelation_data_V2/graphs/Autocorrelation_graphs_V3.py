import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec

parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 5, 6, 7, 8]
#gamma = [0, 0.1, 0.5, 0.9, 1]
gamma = [1, 0.9, 0.5, 0.1, 0]
input = np.arange(0, 30, .5)


colors = ["#EE5A24", "#009432", "#0652DD", "#9980FA", "#833471"]  # colorblind-friendly

colors = ["#12CBC4", "#009432", "#0652DD", "#FFC312", "#EA2027"]  # colorblind-friendly
labels = ['a.)', 'b.)', 'c.)', 'd.)', 'e.)']
row_labels = [
    r'$f = 0.1, \phi = 0$',
    r'$f = 0.1, \phi = \pi/6$'
]

fig, axes = plt.subplots(
    nrows=2, ncols=5,
    figsize=(14, 5), # Originally (10, 5)
    sharex=False, sharey=False
)


fig.subplots_adjust(
    wspace=0.1,   # horizontal spacing
    hspace=0.15    # vertical spacing
)

gap = 0.1
gs = gridspec.GridSpec(2, 5, figure=fig, wspace=gap)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
        ax = axes[i, j]
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
        plt.setp(ax.spines.values(), linewidth=1.6)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.tick_params(axis="x",direction="in", which= 'both')
        ax.tick_params(axis="y",direction="in", which= 'both')
        ax.tick_params(axis="x", which= 'major', length = 5.0, width = 1.2)
        ax.tick_params(axis="x", which= 'minor', length = 3.0, width = 1.2)
        ax.tick_params(axis="y", which= 'major', length = 5.0, width = 1.2)
        ax.tick_params(axis="y", which= 'minor', length = 3.0, width = 1.2)
        #ax.tick_params(axis='y', labelrotation=90)
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        ax.set_xlim(1, np.exp(29.5))
        if i==0:
            ax.set_ylim(-0.5, 1.05)
            ax.set_yticks([-0.5, 0, 0.5, 1])
        if i==1:
            ax.set_ylim(-0.75, 1.05)
            ax.set_yticks([-0.5, 0, 0.5, 1])
        ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        for k, L in enumerate(size):
            data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')
            ax.plot(np.exp(input), data[1], label=f'$L = {L}$', color=colors[k], linewidth=2)
            ax.set_xscale('log')
            ax.xaxis.set_major_locator(ticker.LogLocator(base=10))
        #ax.set_ylim(-0.5, 1.1)
        #ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)


        
        if i == 0 and j == 0:
            ax.legend(fontsize=11,
        loc='upper right',
        frameon=False, bbox_to_anchor=(0.8, 1),
        handlelength=1.2,
        handletextpad=0.4,
        borderpad=0.3,
        labelspacing=0.2)
        if i == 1 and j == 0:
            #ax.set_xlabel('t')
            None
        if i == 1 and j == 0:
            #ax.set_ylabel(f'$A_∞(t)$')
            None

        if j != 0:
            ax.set_yticklabels([])
        if i != 1:
            ax.set_xticklabels([])
        ax.set_xticks([1e3, 1e7, 1e11])
        if i == 1:
            ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
            #ax.set_xticklabels([1, None, 1e7, None])
        #if i == 0:
        #    ax.set_title(f'$Γ={gam}$')

#fig.text(0.03, 0.75, rotation =90, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
#         va='center', ha='left', fontsize=14)

#fig.text(0.03, 0.27, rotation =90, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
#         va='center', ha='left', fontsize=14)
#fig.text(0.05, .9, '(a)', fontsize=18)
#fig.text(0.05, 0.45, '(b)', fontsize=18)

plt.savefig("my_figure_V3.pdf", format="pdf", transparent=True)
plt.savefig("my_figure_V3.png", dpi=600, transparent=True)
#plt.show()