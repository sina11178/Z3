import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as mcolors

'''
base = plt.cm.RdBu
cmap = mcolors.LinearSegmentedColormap.from_list(
    "custom",
    base(np.linspace(0.5, 1, 256))   # start at white
)
'''

cmap = mcolors.LinearSegmentedColormap.from_list(
"cmap",
[
"#FFFFFF", 
"#186DE5", 
"#186DE5", 
"#186DE5", 
"#0D56E9", 
"#08348B", 
"#042368",
"#000C18" 
]
)


parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 6, 8]
gamma = [0, 0.1, 0.5, 0.9, 1]

n_max = 31
plt.figure()
BASE = Path(__file__).resolve().parent
path = (
    BASE
    / ".."
    / f"f_{0.1}_phi_{0}"
    / "seed_squared"
    / f"krylov_matrix_L{8}_gamma_{0}_f_{0.1}_phi_{0}_squared.npy"
).resolve()
data = np.load(path)
print(data.shape)

plt.imshow(np.abs(data[:6, :6]), cmap='seismic', aspect='auto')
plt.colorbar(label='|Matrix|')  # adds color scale
plt.title("Heatmap of Matrix Magnitude")
plt.xlabel("Column")
plt.ylabel("Row")
plt.show()



fig, axes = plt.subplots(
    nrows=4, ncols=5, 
    figsize=(28, 24),
    sharex=False, sharey=False
)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
       ax = axes[2*i, j]
       ax_matrix = axes[2*i+1, j]
       for L in size:
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
            EIG_STATE_SQUARED = np.abs(U[:, index])
            ax.plot(np.arange(n_max), EIG_STATE_SQUARED, label=f'L = {L}, $|λ_{L}|$={np.abs(eigval_min):.3e}')
            ax.legend()
            ax.set_ylim(0, .3)
            ax.set_box_aspect(1)
        
       ax_matrix.imshow(np.abs(data), cmap='bwr', aspect='equal')
       im = ax_matrix.images[0]  # grab the imshow object from the last column
       cbar = fig.colorbar(im, ax=ax_matrix, fraction=0.046, pad=0.04)
       cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #    if i == 0 and j == 0:
    #         ax.legend()
       #if i == 1:
       ax.set_xlabel('n', labelpad=15, fontsize=19) # ITERATION INDEX
       if j == 0:
            ax.set_ylabel(r'$|ψ^{(n)}|$', labelpad=15, fontsize=22)
       if i == 0:
            ax.set_title(f'$Γ={gam}$', pad=15, fontsize=22)
       #if j == 4:
       #     im = ax_matrix.images[0]  # grab the imshow object from the last column
       #     fig.colorbar(im, ax=ax_matrix, fraction=0.046, pad=0.04)


fig.text(0.94, 0.69, rotation =270, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
         va='center', ha='left', fontsize=27)

fig.text(0.94, 0.29, rotation =270, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
         va='center', ha='left', fontsize=27)
fig.text(0.05, .9, '(a)', fontsize=22)
fig.text(0.05, 0.69, '(b)', fontsize=22)
fig.text(0.05, 0.49, '(c)', fontsize=22)
fig.text(0.05, 0.29, '(d)', fontsize=22)

plt.savefig("my_figure.pdf", format="pdf")
