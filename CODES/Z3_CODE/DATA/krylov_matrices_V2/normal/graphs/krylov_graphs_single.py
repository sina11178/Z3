import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec



'''
base = plt.cm.RdBu
cmap = mcolors.LinearSegmentedColormap.from_list(
    "custom",
    base(np.linspace(0.5, 1, 256))   # start at white
)
'''

color_kryl = mcolors.LinearSegmentedColormap.from_list(
'custom',["#FFFFFF","#0D74DB","#F01B1B"]
)


parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 6, 8]
#gamma = [0, 0.1, 0.5, 0.9, 1]
gamma = [1, 0]

n_max = 25
n_data_shown = n_max
#n_max = 31
plt.figure()
BASE = Path(__file__).resolve().parent
path = (
    BASE
    / ".."
    / f"f_{0.1}_phi_{0}"
    / "seed_squared"
    / f"krylov_matrix_L{8}_gamma_{1}_f_{0.1}_phi_{0}_squared.npy"
).resolve()
data = np.load(path)
print(data.shape)
data_abs = np.abs(data)
bn_max = data_abs.max()
plt.imshow(np.abs(data[:25, :25]), cmap=color_kryl, aspect='equal', vmin = 0, vmax = bn_max/3)
plt.colorbar(label='|Matrix|')  # adds color scale
plt.title("Heatmap of Matrix Magnitude")
plt.xlabel("Column")
plt.ylabel("Row")
plt.show()

colors_line = ["#2C0DDA", "#09C96C", "#B81B07",]

fig, axes = plt.subplots(
    nrows=2, ncols=2, 
    figsize=(14, 5),
    sharex=False, sharey=False
)

fig.subplots_adjust(
    wspace=-0.8,   # horizontal spacing
    hspace=0.3    # vertical spacing
)


for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
       ax = axes[i, j]
       #ax_matrix = axes[2*i+1, j]
       for k, L in enumerate(size):
            BASE = Path(__file__).resolve().parent
            path = (
                BASE
                / ".."
                / f"f_{param[0]}_phi_{param[1]}"
                / "seed_squared"
                / f"krylov_matrix_L{L}_gamma_{gam}_f_{param[0]}_phi_{param[1]}_squared.npy"
            ).resolve()
            data = np.load(path)
            #data = np.load(f'../f_{param[0]}_phi_{param[1]}/seed_squared/krylov_matrix_L{L}_gamma_{gam}_f_{param[0]}_phi_{param[1]}_squared.npy')
            data = data[:n_max, :n_max]
            A, U = np.linalg.eig(data)
            index = np.argmin(np.abs(A))
            eigval_min = A[index]
            EIG_STATE_SQUARED = np.abs(U[:, index])**2
            if L == 4:
                ax.plot(np.arange(n_max)[:n_data_shown], EIG_STATE_SQUARED[:n_data_shown], c=colors_line[k], label=f'L = {L}', linewidth= 4, alpha=1)
                #ax.get_legend().remove()
            if L == 6:
                ax.plot(np.arange(n_max)[:n_data_shown], EIG_STATE_SQUARED[:n_data_shown],  c=colors_line[k], label=f'L = {L}', linewidth= 9.5, alpha=0.7)
                #ax.get_legend().remove()
            if L == 8:
                ax.plot(np.arange(n_max)[:n_data_shown], EIG_STATE_SQUARED[:n_data_shown], c=colors_line[k], label=f'L = {L}', linestyle='-', linewidth= 3, alpha=1)
            #ax.legend()
            if  i == 0:
                ax.set_ylim(0,0.003)# EIG_STATE_SQUARED[:20].max()*.9)
                #pass
            if  i == 1:
                ax.set_ylim(0, 0.003)# EIG_STATE_SQUARED[:20].max()*.9)
                #pass
            ax.set_box_aspect(1)
            # ax.set_yscale('log')
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
       ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
       #ax.set_yscale('log')
       #ax.set_ylim(1e-32, 1e15)
       if i == 0:
           ax.set_xticks([])
       if j==1:
           ax.set_yticks([])
       if j == 0:
           #ax.set_yticks([1e-30, 1e1])
           pass
       if i == 0 and j == 1:
           ax.legend(fontsize=8)
       


#fig.text(0.94, 0.69, rotation =270, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
#         va='center', ha='left', fontsize=27)

#fig.text(0.94, 0.29, rotation =270, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
#         va='center', ha='left', fontsize=27)
#fig.text(0.05, .9, '(a)', fontsize=22)
#fig.text(0.05, 0.69, '(b)', fontsize=22)
#fig.text(0.05, 0.49, '(c)', fontsize=22)
#fig.text(0.05, 0.29, '(d)', fontsize=22)

plt.savefig("4_graphs_ZM.pdf", format="pdf", transparent=True)
plt.savefig("4_graphs_ZM.png", transparent=True)
plt.savefig("4_graphs_ZM.svg", format="svg", transparent=True)
