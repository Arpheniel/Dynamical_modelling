"""
MC-значимость per-component средних возраста/[M/H].

N реализаций: наблюдаемые карты возмущаются гауссом по их ошибкам, BVLS
пере-решается на готовом (кэшированном) weight_matrix (с уже зашитой
регуляризацией λ), считаются массово-взвешенные средние по компонентам.
Печатаются mean±σ по компонентам и РАЗНИЦА коротация−контрвращение по
реализациям (учёт корреляции — diff считается внутри каждой реализации).

Запуск из папки галактики: python mc_components.py
"""
import numpy as np
import scipy.optimize
import chemo_unified as cu

np.random.seed(42)          # фиксированный seed (воспроизводимость)
N = 50

W = np.load(cu._rpath("weight_matrix"))["weight_matrix"]
nb = cu.num_bin
age_obs, age_err = cu.bin_age[:nb], cu.bin_age_err
met_obs, met_err = cu.bin_met[:nb], cu.bin_met_err

w = np.asarray(cu.weights[cu.model_index], float)
mass = np.array([float(np.sum(w[cu.sorted_orbs[ir][ilz]])) if cu.sorted_orbs[ir][ilz] else 0.0
                 for ir, ilz in cu.dyn_comps_data])
idx = cu.dyn_comps_data[:, 1]
comps = {"коротация": idx >= 14, "сфероид": (idx >= 7) & (idx < 14), "контрвращ": idx < 7}


def solve(obs_pert, lo, hi):
    rhs = np.zeros(W.shape[0]); rhs[:nb] = obs_pert
    return scipy.optimize.lsq_linear(W, rhs, bounds=(lo, hi), method="bvls",
                                     tol=1e-30, max_iter=200, verbose=0).x


def comp_mean(phi, sel):
    s = sel & (mass > 0); Wm = mass[s].sum()
    return float(np.sum(phi[s] * mass[s]) / Wm)


age_s = {c: [] for c in comps}
met_s = {c: [] for c in comps}
d_age, d_met = [], []          # коротация − контрвращ по реализации
ea = np.where(age_err > 0, age_err, 0.0)
em = np.where(met_err > 0, met_err, 0.0)
for k in range(N):
    pa = solve(age_obs + np.random.normal(0, ea), cu.age_min, cu.age_max)
    pm = solve(met_obs + np.random.normal(0, em), cu.met_min, cu.met_max)
    for c, sel in comps.items():
        age_s[c].append(comp_mean(pa, sel)); met_s[c].append(comp_mean(pm, sel))
    d_age.append(comp_mean(pa, comps["коротация"]) - comp_mean(pa, comps["контрвращ"]))
    d_met.append(comp_mean(pm, comps["коротация"]) - comp_mean(pm, comps["контрвращ"]))

print(f"=== {cu.GALAXY}: per-component mean±σ_MC (N={N}) ===")
print(f"{'компонент':>12} {'age':>14} {'[M/H]':>14}")
for c in comps:
    print(f"{c:>12} {np.mean(age_s[c]):6.2f}±{np.std(age_s[c]):.2f}   "
          f"{np.mean(met_s[c]):6.2f}±{np.std(met_s[c]):.2f}")
da, dm = np.array(d_age), np.array(d_met)
print(f"\nРазница коротация − контрвращ (по реализациям):")
print(f"  Δage   = {da.mean():+.2f} ± {da.std():.2f}  ->  {abs(da.mean())/da.std():.1f}σ")
print(f"  Δ[M/H] = {dm.mean():+.2f} ± {dm.std():.2f}  ->  {abs(dm.mean())/dm.std():.1f}σ")
