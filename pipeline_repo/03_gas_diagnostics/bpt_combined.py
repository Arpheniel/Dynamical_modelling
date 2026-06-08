"""Совместная BPT-диаграмма PGC + LEDA из bpt_points_*.npz."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = ".."  # запуск из handoff/chemistry; npz лежат в рабочих папках
P = np.load(r"../../pgc_final_C/bpt_points_PGC.npz", allow_pickle=True)
L = np.load(r"../../leda_final_C/bpt_points_LEDA.npz", allow_pickle=True)

xx = np.linspace(-1.8, 0.0, 200)
kau = 0.61/(xx-0.05)+1.30          # Kauffmann+03 (SF)
xk  = np.linspace(-1.8, 0.40, 200)
kew = 0.61/(xk-0.47)+1.19          # Kewley+01 (max starburst)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharex=True, sharey=True)
for ax, D, nm, med in [
    (axes[0], P, "PGC 35706",   "LINER-доминирует"),
    (axes[1], L, "LEDA 2220522","SF-доминирует"),
]:
    x = D['x']; y = D['y']; cls = D['cls']
    col = {'SF':'tab:blue','COMP':'tab:green','LINER':'tab:red','SEY':'tab:orange'}
    for c in ['SF','COMP','LINER','SEY']:
        m = cls == c
        if m.any():
            ax.scatter(x[m], y[m], s=10, c=col[c], label=f"{c} ({m.sum()})", alpha=0.6)
    ax.plot(xx, kau, 'k--', lw=1.2, label='Kauffmann+03')
    ax.plot(xk, kew, 'k-',  lw=1.2, label='Kewley+01')
    ax.set_title(f"{nm}: {med}", fontsize=12)
    ax.set_xlabel(r"$\log\,[\mathrm{NII}]6584/\mathrm{H}\alpha$")
    ax.set_xlim(-1.6, 0.6); ax.set_ylim(-1.3, 1.3)
    ax.legend(fontsize=8, loc='lower left')
axes[0].set_ylabel(r"$\log\,[\mathrm{OIII}]5007/\mathrm{H}\beta$")
fig.suptitle("BPT: ионизация газа (бины S/N>3)", fontsize=13)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig("bpt_combined.png", dpi=120)
print("saved bpt_combined.png")
