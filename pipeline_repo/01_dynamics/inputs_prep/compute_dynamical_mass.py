"""
compute_dynamical_mass.py — оценка динамической массы внутри сферы
заданного радиуса, для любой галактики из нашего пайплайна.

Динамическая масса:
    M_dyn(<R) = M_stars(<R) + M_DM(<R) + M_BH

Все три компонента уже доступны как agama.Density-объекты:
  - densityStars: создаётся makeDensityMGE(mge_file, distance, arcsec2kpc, beta).
                  Значения в MGE-файле = СВЕТИМОСТЬ в L_sun/pc^2, так что
                  для перевода в массу — умножить на Upsilon (M/L).
  - densityHalo:  создаётся makeDensityLogHalo(rhalo, vhalo) или makeDensityNFWHalo.
                  Сразу в массовых единицах.
  - M_BH: константа.

Метод agama.Density.enclosedMass(R) возвращает массу внутри СФЕРЫ радиуса R.
Для axisymmetric distributions это близко к 3D-массе в цилиндре проекции
на R_IFU при условии не слишком высокого inclination и компактности галактики.

Использование:
    python compute_dynamical_mass.py --galaxy LEDA
    python compute_dynamical_mass.py --galaxy PGC --r 6 --r 8 --r 10
    python compute_dynamical_mass.py --galaxy LEDA --r-kpc 2 --r-kpc 5
"""

import sys
import argparse
import numpy as np

try:
    import agama
except ImportError:
    print("[FATAL] agama не доступен. Запускайте в окружении пайплайна.")
    sys.exit(1)


# ===== ПАРАМЕТРЫ ГАЛАКТИК =====
# Все значения взяты из run_forstand_grid_*.py файлов.
GALAXY_PARAMS = {
    "LEDA": {
        "name":      "LEDA 2220522",
        "mge_file":  "mge_LEDA_2220522_z_legacy.txt",  # ← скорректировать если другое имя
        "distance":  108643.0,    # kpc
        "rhalo":     131.0,       # kpc (RHALO из best-fit)
        "vhalo":     324.0,       # km/s (VHALO из best-fit)
        "halotype":  "LOG",       # LOG или NFW
        "incl_deg":  27.0,        # INCL = beta из run-script
        "Upsilon":   3.6,         # из таблицы 1 статьи
        "Mbh":       1e6,         # масса чёрной дыры (если есть; иначе 0 безопасно)
        "ifu_arcsec": 8.0,        # радиус MaNGA IFU в arcsec
    },
    "PGC": {
        "name":      "PGC 35706",
        "mge_file":  "mge_PGC35706_z_legacy.txt",
        "distance":  117490.0,    # kpc
        "rhalo":     107.0,       # kpc
        "vhalo":     154.0,       # km/s
        "halotype":  "LOG",
        "incl_deg":  42.0,
        "Upsilon":   12.4,        # из таблицы 1 статьи
        "Mbh":       1e6,
        "ifu_arcsec": 11.0,       # радиус MaNGA IFU
    },
}


def arcsec2kpc(distance_kpc):
    """Угл. сек. -> кпк, при заданной дистанции в кпк."""
    return distance_kpc * np.pi / (180 * 3600)


def compute_masses(params, radii_arcsec=None, radii_kpc=None):
    """Вычисление M_stars, M_DM, M_BH, M_total в каждом из заданных радиусов.

    params: dict с параметрами галактики (см. GALAXY_PARAMS)
    radii_arcsec: список/массив радиусов в arcsec (опционально)
    radii_kpc: список/массив радиусов в kpc (опционально)
    Если ни тот, ни другой не задан — используется ifu_arcsec из params.

    Возвращает list of dicts: для каждого радиуса {R_arcsec, R_kpc,
        M_stars, M_DM, M_BH, M_total}.
    """
    # Подготовка радиусов
    a2k = arcsec2kpc(params["distance"])  # arcsec → kpc; также Agama "length unit"

    # ВАЖНО: agama.setUnits должен быть тот же, что в run-script:
    # mass=1 M_sun, length=arcsec2kpc kpc, velocity=1 km/s.
    # Иначе плотности halo посчитаются с неправильной нормировкой:
    # в makeDensityLogHalo формула rho0 ∝ (v/r)² / G, а G зависит от units.
    # И enclosedMass(R) принимает R в length-units, не в kpc.
    agama.setUnits(mass=1, length=a2k, velocity=1)
    print(f"[Agama units] mass=1 M_sun, length={a2k:.4f} kpc, velocity=1 km/s")
    print(f"  → 1 length-unit = 1 arcsec = {a2k:.4f} kpc")

    if radii_arcsec is None and radii_kpc is None:
        radii_arcsec = [params["ifu_arcsec"]]

    # Радиусы в arcsec (= length-units Agama) и в kpc.
    radii_list = []  # (R_arcsec, R_kpc, R_agama_units)
    if radii_arcsec is not None:
        for r_a in radii_arcsec:
            # 1 arcsec = a2k kpc = 1 length-unit (потому что length=arcsec2kpc)
            radii_list.append((r_a, r_a * a2k, r_a))
    if radii_kpc is not None:
        for r_k in radii_kpc:
            r_a = r_k / a2k
            radii_list.append((r_a, r_k, r_a))

    # Загрузка MGE
    try:
        mge = np.loadtxt(params["mge_file"])
    except (FileNotFoundError, OSError) as e:
        print(f"[ERROR] Не найден MGE-файл '{params['mge_file']}': {e}")
        return None

    # densityStars в условных (luminosity) единицах
    beta = np.radians(params["incl_deg"])
    densityStars = agama.schwarzlib.makeDensityMGE(
        mge, params["distance"], a2k, beta)

    # densityHalo в массовых единицах
    if params["halotype"].upper() == "LOG":
        densityHalo = agama.schwarzlib.makeDensityLogHalo(
            params["rhalo"], params["vhalo"])
    elif params["halotype"].upper() == "NFW":
        densityHalo = agama.schwarzlib.makeDensityNFWHalo(
            params["rhalo"], params["vhalo"])
    else:
        raise ValueError(f"Неизвестный halotype: {params['halotype']}")

    # Полные массы для контекста
    L_total = densityStars.totalMass()  # в условных (luminosity) единицах
    M_stars_total = L_total * params["Upsilon"]
    M_halo_total = densityHalo.totalMass()

    print(f"\n=== {params['name']} ===")
    print(f"Параметры: D={params['distance']/1e3:.0f} Мпк, "
          f"RHALO={params['rhalo']:.0f} kpc, VHALO={params['vhalo']:.0f} km/s, "
          f"Υ={params['Upsilon']:.2f}, halotype={params['halotype']}")
    print(f"L_total (звёздная светимость): {L_total:.3e} L_sun-условных")
    print(f"M_stars_total = Υ × L_total:    {M_stars_total:.3e} M_sun "
          f"(= {M_stars_total/1e10:.2f} × 10^10)")
    print(f"M_halo_total:                   {M_halo_total:.3e} M_sun "
          f"(= {M_halo_total/1e10:.2f} × 10^10)")
    print(f"M_BH:                           {params['Mbh']:.3e} M_sun")
    print()
    print(f"{'R [arcsec]':>10s} {'R [kpc]':>10s} {'M_stars':>15s} "
          f"{'M_DM':>15s} {'M_BH':>15s} {'M_dyn':>15s}    [×10^10 M_sun]")
    print("-" * 95)

    results = []
    for r_a, r_k, r_agama in radii_list:
        # enclosedMass(R) — масса внутри СФЕРЫ радиуса R, R в Agama length-units.
        # Звёзды: умножаем на Υ (M/L отношение).
        M_stars_in_R = densityStars.enclosedMass(r_agama) * params["Upsilon"]
        # DM: по тому же setUnits, в массовых единицах M_sun.
        M_DM_in_R    = densityHalo.enclosedMass(r_agama)
        M_BH_in_R    = params["Mbh"] if r_agama > 1e-3 else 0.0  # BH at origin
        M_dyn_in_R   = M_stars_in_R + M_DM_in_R + M_BH_in_R

        result = {
            "R_arcsec": r_a, "R_kpc": r_k,
            "M_stars": M_stars_in_R, "M_DM": M_DM_in_R,
            "M_BH": M_BH_in_R, "M_total": M_dyn_in_R,
        }
        results.append(result)

        print(f"{r_a:10.2f} {r_k:10.3f} "
              f"{M_stars_in_R/1e10:15.4f} "
              f"{M_DM_in_R/1e10:15.4f} "
              f"{M_BH_in_R/1e10:15.4f} "
              f"{M_dyn_in_R/1e10:15.4f}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Расчёт динамической массы галактики внутри сферы радиуса R.")
    parser.add_argument("--galaxy", required=True, choices=list(GALAXY_PARAMS.keys()),
                        help="Какую галактику считать (LEDA/PGC).")
    parser.add_argument("--r", type=float, action="append", default=None,
                        help="Радиус(ы) в arcsec. Можно указать несколько раз.")
    parser.add_argument("--r-kpc", type=float, action="append", default=None,
                        help="Радиус(ы) в kpc. Можно указать несколько раз.")
    parser.add_argument("--upsilon", type=float, default=None,
                        help="Переопределить Υ (по умолчанию из таблицы).")
    parser.add_argument("--mbh", type=float, default=None,
                        help="Переопределить массу BH в M_sun.")
    args = parser.parse_args()

    params = dict(GALAXY_PARAMS[args.galaxy])
    if args.upsilon is not None:
        params["Upsilon"] = args.upsilon
    if args.mbh is not None:
        params["Mbh"] = args.mbh

    compute_masses(params,
                   radii_arcsec=args.r,
                   radii_kpc=args.__dict__["r_kpc"])


if __name__ == "__main__":
    main()
