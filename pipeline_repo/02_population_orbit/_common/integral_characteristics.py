"""
integral_characteristics.py — таблица интегральных характеристик
Schwarzschild-модели (Cappellari+2013 Tab.1; Boardman+2016 Tab.2;
Santucci+2022 Tab.3): M(<Re), M(<5Re), v_circ(Re), M*/L, f_DM(<Re),
а также форма (p, q, T) для каждой галактики.

  - Лучшая ячейка (r_halo, v_halo, Upsilon) берётся из grid_chi2_GH.npz
    (argmin chi^2), либо переопределяется ключами --rhalo/--vhalo/--upsilon.
  - Масштабы (distance, incl, arcsec2kpc, Mbh) читаются из grid_chi2_GH.npz.
  - Re и q_obs считаются из MGE-файла (без schwarzlib).
  - v_circ(Re) — из M_dyn(<Re) и Re по G (astropy), без schwarzlib.

╔══════════════════════════════════════════════════════════════════════╗
║  ВНИМАНИЕ. Массы (M_stars, M_DM) считаются через schwarzlib. Локальная ║
║  agama.schwarzlib — БЕЗ пользовательских правок, поэтому числа масс,   ║
║  полученные не в pipeline-окружении, НЕ АВТОРИТЕТНЫ. Скрипт пытается   ║
║  импортировать модифицированный schwarzlib_really_fixed.py; если его   ║
║  нет — падает на agama.schwarzlib и громко предупреждает. Итоговую     ║
║  таблицу следует пересчитать в pipeline-окружении на финальной модели. ║
╚══════════════════════════════════════════════════════════════════════╝

Запуск из папки галактики:
    python integral_characteristics.py --galaxy PGC
    python integral_characteristics.py --galaxy LEDA --rhalo 131 --vhalo 324 --upsilon 3.6
"""

import os
import glob
import argparse
import numpy as np

# G в единицах kpc * (км/с)^2 / M_sun (astropy CODATA) — для v_circ, без agama-юнитов
G_KPC_KMS2_MSUN = 4.300917270e-6


GALAXY_PARAMS = {
    "PGC": {
        "name":       "PGC 35706",
        "mge_glob":   "mge_PGC35706*z*legacy.txt",
        "halotype":   "LOG",
        "ifu_arcsec": 11.0,
    },
    "LEDA": {
        "name":       "LEDA 2220522",
        "mge_glob":   "mge_LEDA*z*legacy.txt",
        "halotype":   "LOG",
        "ifu_arcsec": 8.0,
    },
}


def import_schwarzlib():
    """Предпочесть модифицированный schwarzlib_really_fixed.py; иначе agama.schwarzlib."""
    try:
        import schwarzlib_really_fixed as sl
        return sl, "schwarzlib_really_fixed.py (локальный модифицированный файл)", True
    except Exception:
        import agama
        if hasattr(agama, "schwarzlib"):
            return (agama.schwarzlib,
                    "agama.schwarzlib (ВНИМАНИЕ: вероятно БЕЗ правок — числа не авторитетны)",
                    False)
        raise ImportError("schwarzlib недоступен ни как файл, ни как agama.schwarzlib")


def load_mge(mge_glob):
    cands = glob.glob(mge_glob)
    if not cands:
        raise FileNotFoundError(f"MGE по шаблону '{mge_glob}' не найден в {os.getcwd()}")
    mge = np.loadtxt(cands[0])
    return mge, cands[0]


def effective_radius_arcsec(mge):
    """Круговой эффективный радиус (половина проектированного света) из MGE.
    mge: колонки [I (Lsun/pc^2), sigma (arcsec), q]. Возвращает Re в arcsec и
    светимостно-взвешенное q_obs."""
    I, sig, q = mge[:, 0], mge[:, 1], mge[:, 2]
    Lk = I * sig**2 * q                      # ∝ полный свет гауссианы (нормировка сокращается)
    Ltot = np.sum(Lk)
    q_obs = np.sum(Lk * q) / Ltot            # светимостно-взвешенное q

    def Lcum(R):
        return np.sum(Lk * (1.0 - np.exp(-R**2 / (2.0 * sig**2))))

    target = 0.5 * Ltot
    lo, hi = 1e-3, 1e3
    for _ in range(200):                     # бисекция
        mid = 0.5 * (lo + hi)
        if Lcum(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), q_obs


def intrinsic_q(q_obs, incl_deg):
    """Внутреннее сжатие из наблюдаемого q и наклона (oblate axisymmetric)."""
    cosi = np.cos(np.radians(incl_deg))
    val = q_obs**2 - cosi**2
    if val <= 0:
        return np.nan
    return np.sqrt(val) / np.sin(np.radians(incl_deg))


def main():
    ap = argparse.ArgumentParser(description="Таблица интегральных характеристик.")
    ap.add_argument("--galaxy", default=None, choices=list(GALAXY_PARAMS) + [None])
    ap.add_argument("--grid", default="grid_chi2_GH.npz")
    ap.add_argument("--rhalo", type=float, default=None, help="Переопределить r_halo (arcsec).")
    ap.add_argument("--vhalo", type=float, default=None, help="Переопределить v_halo (км/с).")
    ap.add_argument("--upsilon", type=float, default=None, help="Переопределить Upsilon.")
    ap.add_argument("--halotype", default=None, choices=["LOG", "NFW", None])
    args = ap.parse_args()

    galaxy = args.galaxy
    if galaxy is None:
        base = os.path.basename(os.getcwd()).upper()
        galaxy = "LEDA" if "LEDA" in base else "PGC"
    P = GALAXY_PARAMS[galaxy]
    halotype = (args.halotype or P["halotype"]).upper()

    import agama
    sl, sl_src, sl_fixed = import_schwarzlib()

    print("=" * 70)
    print(f"Интегральные характеристики: {P['name']}")
    print(f"schwarzlib: {sl_src}")
    if not sl_fixed:
        print("!!! Числа масс НЕ авторитетны локально — пересчитать в pipeline-окружении.")
    print("=" * 70)

    # ── grid: best-fit + масштабы ────────────────────────────────────────────
    d = np.load(args.grid, allow_pickle=True)
    rhalo_arr = np.asarray(d["rhalo"], float)
    vhalo_arr = np.asarray(d["vhalo"], float)
    chi2 = np.asarray(d["chi2"], float)
    ups_arr = np.asarray(d["upsilon"], float)
    i_b, j_b = np.unravel_index(np.nanargmin(chi2), chi2.shape)

    rhalo = args.rhalo if args.rhalo is not None else float(rhalo_arr[i_b])
    vhalo = args.vhalo if args.vhalo is not None else float(vhalo_arr[j_b])
    upsilon = args.upsilon if args.upsilon is not None else float(ups_arr[i_b, j_b])

    distance = float(d["distance"]) if "distance" in d.files else None
    incl = float(d["incl"]) if "incl" in d.files else None
    Mbh = float(d["Mbh"]) if "Mbh" in d.files else 0.0
    a2k = float(d["arcsec2kpc"]) if "arcsec2kpc" in d.files else \
        (distance * np.pi / 648000.0 if distance else None)

    print(f"best-fit грида: r_halo={rhalo:g} arcsec, v_halo={vhalo:g} км/с, "
          f"Upsilon={upsilon:.3f}  (chi2_min={np.nanmin(chi2):.0f})")
    print(f"масштабы: distance={distance:g} kpc, incl={incl:g} deg, "
          f"arcsec2kpc={a2k:.4f}, Mbh={Mbh:g}")

    # ── MGE: Re, q ───────────────────────────────────────────────────────────
    mge, mge_file = load_mge(P["mge_glob"])
    Re_arcsec, q_obs = effective_radius_arcsec(mge)
    Re_kpc = Re_arcsec * a2k
    q_int = intrinsic_q(q_obs, incl)
    print(f"MGE: {mge_file}  ->  Re={Re_arcsec:.2f} arcsec = {Re_kpc:.3f} kpc, "
          f"q_obs={q_obs:.3f}, q_intr={q_int:.3f}")

    # ── плотности через schwarzlib ───────────────────────────────────────────
    agama.setUnits(mass=1, length=a2k, velocity=1)   # length-unit = 1 arcsec
    beta = np.radians(incl)
    densityStars = sl.makeDensityMGE(mge, distance, a2k, beta)
    if halotype == "LOG":
        densityHalo = sl.makeDensityLogHalo(rhalo, vhalo)
    else:
        densityHalo = sl.makeDensityNFWHalo(rhalo, vhalo)

    def masses_within(R_arcsec):
        M_s = densityStars.enclosedMass(R_arcsec) * upsilon
        M_d = densityHalo.enclosedMass(R_arcsec)
        M_b = Mbh if R_arcsec > 1e-3 else 0.0
        return M_s, M_d, M_b, M_s + M_d + M_b

    radii = [("Re", Re_arcsec), ("5 Re", 5 * Re_arcsec), ("R_IFU", P["ifu_arcsec"])]

    print()
    hdr = f"{'апертура':>8s} {'R[arcsec]':>10s} {'R[kpc]':>8s} " \
          f"{'M*[1e10]':>10s} {'M_DM[1e10]':>11s} {'M_dyn[1e10]':>12s} {'f_DM':>6s}"
    print(hdr); print("-" * len(hdr))
    rows = {}
    for tag, R_a in radii:
        M_s, M_d, M_b, M_t = masses_within(R_a)
        fdm = M_d / M_t if M_t > 0 else np.nan
        rows[tag] = dict(R_a=R_a, R_k=R_a * a2k, Ms=M_s, Md=M_d, Mt=M_t, fdm=fdm)
        print(f"{tag:>8s} {R_a:10.2f} {R_a*a2k:8.3f} "
              f"{M_s/1e10:10.4f} {M_d/1e10:11.4f} {M_t/1e10:12.4f} {fdm:6.3f}")

    # ── производные ──────────────────────────────────────────────────────────
    M_dyn_Re = rows["Re"]["Mt"]
    vcirc_Re = np.sqrt(G_KPC_KMS2_MSUN * M_dyn_Re / Re_kpc)   # M[Msun], Re[kpc] -> км/с
    p_shape, T_shape = 1.0, 0.0                               # axisymmetric oblate

    print()
    print(f"v_circ(Re)   = {vcirc_Re:7.1f} км/с")
    print(f"M*/L         = {upsilon:7.3f} M_sun/L_sun")
    print(f"f_DM(<Re)    = {rows['Re']['fdm']:7.3f}")
    print(f"форма:  p = {p_shape:.2f},  q = {q_int:.3f},  T = {T_shape:.2f}  (axisymmetric)")

    # ── сохранить таблицу и LaTeX-строку ─────────────────────────────────────
    out_txt = f"{galaxy}_integral_characteristics.txt"
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"# {P['name']} — интегральные характеристики Schwarzschild-модели\n")
        f.write(f"# schwarzlib: {sl_src}\n")
        if not sl_fixed:
            f.write("# ВНИМАНИЕ: массы не авторитетны — пересчитать в pipeline-окружении\n")
        f.write(f"# best-fit: r_halo={rhalo:g} arcsec, v_halo={vhalo:g} км/с, "
                f"Upsilon={upsilon:.3f}, halotype={halotype}\n")
        f.write(f"Re_arcsec {Re_arcsec:.3f}\nRe_kpc {Re_kpc:.4f}\n")
        for tag in ("Re", "5 Re", "R_IFU"):
            r = rows[tag]
            f.write(f"{tag.replace(' ', '')}_Mstar_Msun {r['Ms']:.4e}\n")
            f.write(f"{tag.replace(' ', '')}_MDM_Msun {r['Md']:.4e}\n")
            f.write(f"{tag.replace(' ', '')}_Mdyn_Msun {r['Mt']:.4e}\n")
            f.write(f"{tag.replace(' ', '')}_fDM {r['fdm']:.4f}\n")
        f.write(f"vcirc_Re_kms {vcirc_Re:.2f}\nML {upsilon:.4f}\n")
        f.write(f"q_obs {q_obs:.4f}\nq_intr {q_int:.4f}\np 1.0\nT 0.0\n")

    # одна строка для таблицы статьи (LaTeX)
    latex = (f"{P['name']} & {Re_kpc:.2f} & {M_dyn_Re/1e10:.2f} & "
             f"{rows['5 Re']['Mt']/1e10:.2f} & {vcirc_Re:.0f} & {upsilon:.2f} & "
             f"{rows['Re']['fdm']:.2f} & {p_shape:.2f} & {q_int:.2f} & {T_shape:.2f} \\\\")
    out_tex = f"{galaxy}_integral_characteristics_row.tex"
    with open(out_tex, "w", encoding="utf-8") as f:
        f.write("% столбцы: Галактика & Re[kpc] & M(<Re)[1e10] & M(<5Re)[1e10] & "
                "vcirc(Re)[км/с] & M*/L & fDM(<Re) & p & q & T\n")
        f.write(latex + "\n")

    print(f"\nсохранено: {out_txt}")
    print(f"сохранено: {out_tex}")
    print("LaTeX-строка:", latex)


if __name__ == "__main__":
    main()
