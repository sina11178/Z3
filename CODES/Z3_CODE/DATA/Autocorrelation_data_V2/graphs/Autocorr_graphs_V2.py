import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 6, 8]
gamma = [0, 0.1, 0.5, 0.9, 1]
input = np.arange(0, 30, .5)

# Publication style settings
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 13,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # colorblind-friendly
labels = ['a.)', 'b.)', 'c.)', 'd.)', 'e.)']
row_labels = [
    r'$f = 0.1, \phi = 0$',
    r'$f = 0.1, \phi = \pi/6$'
]

fig, axes = plt.subplots(
    nrows=2, ncols=5,
    figsize=(18, 7),
    sharex='col', sharey='row',
    constrained_layout=True
)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
        ax = axes[i, j]
        for k, L in enumerate(size):
            data = np.load(f'../f_{param[0]}_phi_{param[1]}/size_autocorr_comparison_gamma_{gam}_size_{L}_NPY.npy')
            ax.plot(np.exp(input), data[1], label=f'$L = {L}$', color=colors[k], linewidth=1.5)

        ax.set_xscale('log')
        ax.set_ylim(-0.5, 1.1)
        ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

        # Column titles (top row only)
        if i == 0:
            ax.set_title(f'$\\Gamma = {gam}$', fontsize=15, pad=6)
            # Subplot labels (a, b, c...)
            ax.text(0.05, 0.95, labels[j], transform=ax.transAxes,
                    fontsize=13, fontweight='bold', va='top', ha='left')

        # Row labels (rightmost column only)
        if j == 4:
            ax.text(1.05, 0.5, row_labels[i], transform=ax.transAxes,
                    fontsize=13, va='center', ha='left', rotation=90)

        # x-axis label (bottom row only)
        if i == 1:
            ax.set_xlabel('$t$', fontsize=14)

        # y-axis label (leftmost column only)
        if j == 0:
            ax.set_ylabel('$C(t)$', fontsize=14)

        # Vertical dashed separator after gamma=0.5 (col index 2)
        if j == 2:
            ax.axvline(x=ax.get_xlim()[1], color='black', linewidth=1,
                       linestyle='--', alpha=0.4, clip_on=False)

        # Legend — top-left panel only, transparent background
        if i == 0 and j == 0:
            legend = ax.legend(
                loc='lower left',
                framealpha=0.3,        # transparent
                edgecolor='gray',
                fontsize=13,
                handlelength=1.5,
            )

plt.savefig("my_figure_V2.pdf", format="pdf", bbox_inches='tight')
plt.savefig("my_figure_V2.png", dpi=600, bbox_inches='tight')
plt.show()