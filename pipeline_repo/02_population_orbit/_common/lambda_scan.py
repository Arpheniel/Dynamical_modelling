"""
Скан силы регуляризации lambda_reg для BVLS-декомпозиции (финал PGC).

Идея: data-блок weight_matrix (строки 0..num_bin-1) НЕ зависит от lambda_reg,
рег-блок (строки num_bin..) линеен по lambda_reg. Поэтому берём data-блок из
готового кэша (results_cache_*/weight_matrix.npz), рег-матрицу строим из
смежности ячеек (reg_matrix_1 при lambda=1), и для каждого lambda решаем BVLS
заново — без чтения орбит.

Для каждого lambda печатаем:
  - chi2_red (ошибко-взвешенный, по data-бинам) для возраста и [M/H];
  - долю ячеек, упёртых в границы сетки (railing) — индикатор вырождения.
Выбираем lambda на «колене»: railing уже сбит, chi2_red ещё приемлем.
"""
import numpy as np
import scipy.optimize
import chemo_unified as cu

W0 = np.load(cu._rpath("weight_matrix"))["weight_matrix"]   # (num_bin*2-3, n_cells)
num_bin = cu.num_bin
ncell = cu.dyn_comp_num
data = W0[:num_bin]                                          # λ-независимый data-блок

# Рег-матрица при lambda=1 (масштабируется λ): строки bin_i, столбцы ячеек
nreg = num_bin - 3
L1 = np.zeros((nreg, ncell))
for bi in range(nreg):
    for ci in range(ncell):
        L1[bi, ci] = cu.reg_matrix_1(ci, bi, cu.dyn_comps_data_map, cu.dyn_comps_data, 1.0)

age_obs, age_err = cu.bin_age[:num_bin], cu.bin_age_err
met_obs, met_err = cu.bin_met[:num_bin], cu.bin_met_err
va = age_err > 0
vm = met_err > 0


def solve(L_scaled, rhs_obs, lo, hi):
    W = np.vstack([data, L_scaled])
    rhs = np.concatenate([rhs_obs, np.zeros(nreg)])
    sol = scipy.optimize.lsq_linear(W, rhs, bounds=(lo, hi), method="bvls",
                                    tol=1e-30, max_iter=200, verbose=0)
    return sol.x


def chi2_red(phi, obs, err, valid):
    pred = data @ phi
    r = (pred[valid] - obs[valid]) / err[valid]
    return float(np.sum(r * r) / valid.sum())


def rail(phi, lo, hi):
    return float(np.mean(np.isclose(phi, lo) | np.isclose(phi, hi)) * 100)


lambdas = [0.0, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
print(f"{'lambda':>8} | {'chi2r_age':>9} {'rail_age%':>9} | {'chi2r_met':>9} {'rail_met%':>9}")
print("-" * 56)
for lam in lambdas:
    Ls = lam * L1
    pa = solve(Ls, age_obs, cu.age_min, cu.age_max)
    pm = solve(Ls, met_obs, cu.met_min, cu.met_max)
    print(f"{lam:>8.4f} | {chi2_red(pa, age_obs, age_err, va):>9.1f} "
          f"{rail(pa, cu.age_min, cu.age_max):>9.0f} | "
          f"{chi2_red(pm, met_obs, met_err, vm):>9.2f} "
          f"{rail(pm, cu.met_min, cu.met_max):>9.0f}")
