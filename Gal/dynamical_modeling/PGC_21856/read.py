import numpy as np, agama as _agama

N = 12000  # number of orbits
M_total = 2e9  # in units of M_sun
Filename = "M1e+07_O80_Rh150_Vh180_i55_a40_N12000_R5.00_GH_DensityCylindricalTopHat.npz"
R_scale = 0.436

archive = np.load(Filename, allow_pickle=True, encoding='latin1')

LOSVD = archive['LOSVD']
Lambda_z = archive['MOD_Lambda_z']
Vel = archive["Velocity"]
circ2 = archive["circ2"]
Upsilon = archive["Upsilon"]
weights = archive["weights"]
Rmean = archive["Rmean"]
cosincl = archive["cosincl"]
ic = archive["ic"]
inttime = archive["inttime"]

Units = _agama.getUnits()

orbits = []
mass_c_p = []
mass_c_n = []
mass_sigma = []
M_p = 0
M_n = 0
M_s = 0

for i in range(0, N):
    if weights[4][i] >= 0.00001:
        orbits.append([circ2[i], weights[1][i] * M_total, Rmean[i] * R_scale, cosincl[i], Lambda_z[i], Vel[i]])

for r in range(0, 30):
    for i in range(0, len(orbits)):
        if orbits[i][3] > 0.5 and orbits[i][2] > (r - 0.5) and orbits[i][2] < (r + 0.5):
            M_p = M_p + orbits[i][1]
        if orbits[i][3] < -0.5 and orbits[i][2] > (r - 0.5) and orbits[i][2] < (r + 0.5):
            M_n = M_n + orbits[i][1]
        if orbits[i][2] > (r - 0.5) and orbits[i][2] < (r + 0.5):
            M_s = M_s + orbits[i][1]

    print("#################", M_s, M_p, M_n)
    mass_sigma.append(M_s)
    mass_c_p.append(M_p)
    mass_c_n.append(M_n)
    M_p = 0
    M_n = 0
    M_s = 0

import matplotlib.pyplot as plt

matrix = [[], [], [], [], []]

for i in range(0, len(orbits)):
    matrix[0].append(orbits[i][0])
    matrix[1].append(orbits[i][1])
    matrix[2].append(orbits[i][2])
    matrix[3].append(orbits[i][3])
    matrix[4].append(orbits[i][4])

fig, ax = plt.subplots()

ax.plot(mass_c_p, label="mass of co-rotating component")

ax.plot(mass_c_n, label="mass of counter-rotating component")

ax.plot(mass_sigma, label="mass of all stars")

fig.legend()
ax.set_xlabel("R_mean,kpc")
ax.set_ylabel("weight,M_sun")
ax.set_title("M_Total=%.0f M_sun" % M_total)
fig.savefig("mass_dest")
fig.show

Vel_c_p = []
Vel_c_n = []
Vel_sigma = []
V_p = []
V_n = []
V_s = []

for r in range(0, 30):
    r = r / 2
    for i in range(0, len(orbits)):
        if orbits[i][3] > 0.5 and orbits[i][2] > (r - 0.25) and orbits[i][2] < (r + 0.25):
            V_p.append(orbits[i][5])
        if orbits[i][3] < -0.5 and orbits[i][2] > (r - 0.25) and orbits[i][2] < (r + 0.25):
            V_n.append(orbits[i][5])
        if orbits[i][2] > (r - 0.25) and orbits[i][2] < (r + 0.25):
            V_s.append(orbits[i][5])
    print("#################", np.mean(V_s), np.mean(V_p), np.mean(V_n))
    Vel_sigma.append(np.mean(V_s))
    Vel_c_p.append(np.mean(V_p))
    Vel_c_n.append(np.mean(V_n))
    V_p = []
    V_n = []
    V_s = []

fig, ax = plt.subplots()

ax.plot(Vel_c_p, label="mean velocity of co-rotating component")

ax.plot(Vel_c_n, label="mean velocity of counter-rotating component")

ax.plot(Vel_sigma, label="mean velocity of all stars")

fig.legend()
ax.set_xlabel("R_mean,kpc")
ax.set_ylabel("Velcity,km/s")
ax.set_title("Velocity")
fig.savefig("vel")

fig, ax = plt.subplots()

lambda_z_y = []
lambda_z_x = []
Rad = np.linspace(0,10,20)
lambda_z_pl= np.linspace(-1,1,20)
for r in Rad:
    lambda_z_y = []
    for lambda_z in lambda_z_pl:
        lambda_z_pr = 0
        for i in range(0, len(orbits)):
            if orbits[i][2] > (r - 0.25) and orbits[i][2] < (r + 0.25) and orbits[i][4] > (lambda_z - 0.05) and orbits[i][4] < (lambda_z + 0.05):
                lambda_z_pr = lambda_z_pr + orbits[i][1]
        lambda_z_y.append(lambda_z_pr)
    lambda_z_x.append(lambda_z_y)
ax.contourf(lambda_z_x)
fig.savefig("lambda_z")

fig, axs = plt.subplots(2, 2, layout='constrained')
pc = axs[0, 0].scatter(matrix[1], matrix[2], c=matrix[3], s=5, cmap='plasma')
fig.colorbar(pc, ax=axs[0, 0], extend='both', label='cos(incl)')
# axs[0, 0].set_title('c=cos(incl)')
axs[0, 0].set_xlabel('weights,M_sun')
axs[0, 0].set_ylabel('R_mean,kpc')

pc = axs[1, 0].scatter(matrix[1], matrix[2], c=matrix[0], s=5, cmap='plasma')
fig.colorbar(pc, ax=axs[1, 0], extend='both', label='circ2')
# axs[1, 0].set_title('c=circ2')
axs[1, 0].set_xlabel('weights,M_sun')
axs[1, 0].set_ylabel('R_mean,kpc')

pc = axs[0, 1].scatter(matrix[2], matrix[3], c=matrix[1], s=5, cmap='plasma')
fig.colorbar(pc, ax=axs[0, 1], extend='both', label='weights,M_sun')
# axs[0, 1].set_title('c=weights,M_sun')
axs[0, 1].set_xlabel('R_mean,kpc')
axs[0, 1].set_ylabel('cos(incl)')

pc = axs[1, 1].scatter(matrix[1], matrix[0], c=matrix[3], s=5, cmap='plasma')
fig.colorbar(pc, ax=axs[1, 1], extend='both', label='cos(incl)')
# axs[1, 1].set_title('c=cos(incl)')
axs[1, 1].set_xlabel('weights,M_sun')
axs[1, 1].set_ylabel('circ2')

fig.savefig("plots")

# data = _numpy.transpose(_numpy.array([circ2,weights[5],Rmean,cosincl]))
# np.savetxt("/home/denis/Documents/agama/Gal/data.txt",sorted(orbits, key = lambda orbits: orbits[1]), fmt="%8.3f", header="circ2 weights Rmean Lz/l")
