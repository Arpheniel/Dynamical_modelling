"""
Независимая перепроверка головных чисел и привязки компонент.
Запускать из рабочей папки (pgc_final_C / leda_final_C).
"""
import sys, numpy as np
import chemo_unified as cu

GAL = cu.GALAXY
cache = f"results_cache_{GAL}_m{cu.model_index}"
sol = np.load(f"{cache}/lsq_solution.npz")
model_age = sol["model_age"]; model_met = sol["model_met"]

# масса каждой dyn-компоненты (как в compute_mass_fractions)
comp_mass = np.array([
    np.sum(cu.weights[cu.model_index][cu.sorted_orbs[ir][ilz]]) if cu.sorted_orbs[ir][ilz] else 0.0
    for ir, ilz in cu.dyn_comps_data])
idx_lz = cu.dyn_comps_data[:, 1]
total = comp_mass.sum()

print(f"\n===== {GAL}: перепроверка =====")
print("Доли масс (спутники регионов LZ_REGIONS):")
regions = [("counterrotating (λz<0)", 0, 7), ("spherical", 7, 14), ("corotating (λz>0)", 14, 21)]
for nm, a, b in regions:
    m = (idx_lz >= a) & (idx_lz < b)
    frac = comp_mass[m].sum()/total
    wage = np.average(model_age[m], weights=comp_mass[m])
    wmet = np.average(model_met[m], weights=comp_mass[m])
    print(f"  {nm:22s} доля={frac*100:5.1f}%  <age>={wage:5.2f}  <[M/H]>={wmet:+.3f}")

# --- привязка: какой λz-регион = бо́льшая часть света (chem[6]=MAP_VEL1) ---
corotV   = np.array(cu.DYN_COMP_LOSVD_MAP_INDICES(14, 21)[0], float)
counterV = np.array(cu.DYN_COMP_LOSVD_MAP_INDICES(0, 7)[0], float)
obsV = np.array(cu.chem[6].data, float)  # MAP_VEL1 = бо́льшая часть света
def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    x = a[m]-a[m].mean(); y = b[m]-b[m].mean()
    d = np.sqrt((x*x).sum()*(y*y).sum())
    return float((x*y).sum()/d) if d else np.nan
print("\nПривязка (chem[6]=MAP_VEL1 = бо́льшая часть света):")
print(f"  corr(MAP_VEL1, corot_model.T)   = {corr(obsV, corotV.T):+.3f}  (ждём +, corot=осн.свет)")
print(f"  corr(MAP_VEL1, counter_model.T) = {corr(obsV, counterV.T):+.3f}  (ждём −)")
print(f"  corr(MAP_VEL1, corot_model asis)= {corr(obsV, corotV):+.3f}  (ждём ~0: без .T диполь ⟂)")
