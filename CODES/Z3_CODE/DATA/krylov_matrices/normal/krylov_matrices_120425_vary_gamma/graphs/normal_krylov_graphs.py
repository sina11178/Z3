import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


parameters = [[0.1, 0], [0.1, 0.524]]
size = [4, 6, 8]
gamma = [0, 0.4, 0.8, 1]

n_max = 50
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
    nrows=2, ncols=4, 
    figsize=(18, 8),
    sharex=False, sharey=False
)

for i, param in enumerate(parameters):
    for j, gam in enumerate(gamma):
       ax = axes[i, j]
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
            ax.set_ylim(0, 0.1)

    #    if i == 0 and j == 0:
    #         ax.legend()
       if i == 1:
            ax.set_xlabel('iteration index - n')
       if j == 0:
            ax.set_ylabel(f'$|ψ_{min}^n|$')
       if i == 0:
            ax.set_title(f'$Γ={gam}$')


fig.text(0.03, 0.75, rotation =90, s=f"(f = {parameters[0][0]}, $φ$ = {parameters[0][1]})", 
         va='center', ha='left', fontsize=14)

fig.text(0.03, 0.27, rotation =90, s=f"(f = {parameters[1][0]}, $φ$ = $π/6$)", 
         va='center', ha='left', fontsize=14)
fig.text(0.05, .9, '(a)', fontsize=18)
fig.text(0.05, 0.45, '(b)', fontsize=18)

plt.savefig("my_figure.pdf", format="pdf")
