"""
ПРЯМАЯ проверка по данным MaNGA DAP (без Schwarzschild-модели, без λz, без .T):
вращается ли газ Hα в ту же сторону, что и СВЕТИМОСТНО-ВЗВЕШЕННАЯ масса звёзд.

STELLAR_VEL  = средняя (по свету) звёздная скорость, одна гауссиана DAP.
EMLINE_GVEL[Ha] = скорость газа.
Сравниваем знак: corr(stars, gas) > 0 -> со-вращение; < 0 -> противовращение.
Плюс сверяем с NBursts VEL1(=газ по LOSVD1)/VEL2 из results-fits (тот же грид).
"""
import sys
import numpy as np
from astropy.io import fits

GAL  = sys.argv[1] if len(sys.argv) > 1 else "PGC"
CFG = {
    "PGC":  dict(maps="pgc_35706_CD/pgc_35706_CD/manga-8992-3704-MAPS-VOR10-MILESHC-MASTARSSP.fits",
                 res ="pgc_final_C/results_8992-3704_vorb020_md19_ad-1_nmom4.fits"),
    "LEDA": dict(maps="leda_2220522_CD/leda_2220522_CD/manga-8254-1902-MAPS-VOR10-MILESHC-MASTARSSP.fits",
                 res ="leda_final_C/results_8254-1902_vorb020_md19_ad-1_nmom4.fits"),
}[GAL]

h = fits.open(CFG["maps"])
stV   = np.array(h["STELLAR_VEL"].data, float)
stMsk = np.array(h["STELLAR_VEL_MASK"].data, int)
HA = 23  # 0-based channel C24 = Ha-6564
gV    = np.array(h["EMLINE_GVEL"].data[HA], float)
gMsk  = np.array(h["EMLINE_GVEL_MASK"].data[HA], int)
gFlux = np.array(h["EMLINE_GFLUX"].data[HA], float)

# NBursts (results fits), тот же 42x42 грид
r = fits.open(CFG["res"])
nb1 = np.array(r[6].data, float)    # MAP_VEL1 (газ = LOSVD1)
nb2 = np.array(r[14].data, float)   # MAP_VEL2

# валидные пиксели DAP: mask==0
stV_ok = stMsk == 0
gV_ok  = (gMsk == 0) & (gFlux > 0)

def demean(a, m):
    a = a.copy(); a[~m] = np.nan
    a[m] -= np.nanmean(a[m]); return a

stVd = demean(stV, stV_ok)
gVd  = demean(gV,  gV_ok)

def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5: return np.nan, int(m.sum())
    x = a[m]-a[m].mean(); y = b[m]-b[m].mean()
    d = np.sqrt((x*x).sum()*(y*y).sum())
    return (float((x*y).sum()/d) if d else np.nan), int(m.sum())

print(f"\n===== {GAL}: ПРЯМАЯ проверка газ vs звёзды (MaNGA DAP) =====")
print(f"valid pix: STELLAR_VEL={int(stV_ok.sum())}, Ha_GVEL={int(gV_ok.sum())}")
print(f"amp: stellar |V| p98={np.nanpercentile(np.abs(stVd),98):.1f}  gas |V| p98={np.nanpercentile(np.abs(gVd),98):.1f} км/с")

print("\nКорреляции (одна и та же сетка, БЕЗ модели/транспонирования):")
r1,n1 = corr(stVd, gVd);            print(f"  STELLAR_VEL  vs  Ha_gas      = {r1:+.3f}   (N={n1})   <== РЕШАЮЩЕЕ")
r2,_  = corr(stVd, demean(nb1,np.isfinite(nb1))); print(f"  STELLAR_VEL  vs  NBursts VEL1 = {r2:+.3f}")
r3,_  = corr(stVd, demean(nb2,np.isfinite(nb2))); print(f"  STELLAR_VEL  vs  NBursts VEL2 = {r3:+.3f}")
r4,_  = corr(gVd,  demean(nb1,np.isfinite(nb1))); print(f"  Ha_gas       vs  NBursts VEL1 = {r4:+.3f}   (sanity: газ=LOSVD1)")
r5,_  = corr(gVd,  demean(nb2,np.isfinite(nb2))); print(f"  Ha_gas       vs  NBursts VEL2 = {r5:+.3f}")

print("\nИТОГ:", "газ СО-вращается с бо́льшей частью звёзд" if r1>0.2 else
      ("газ ПРОТИВО-вращается относительно бо́льшей части звёзд" if r1<-0.2 else
       "знак НЕОПРЕДЕЛЁН (|corr|<0.2) — звёздное поле слабое/неоднозначное"))

# фигура
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
vmax=np.nanpercentile(np.abs(np.concatenate([stVd[np.isfinite(stVd)],gVd[np.isfinite(gVd)]])),98)
fig,ax=plt.subplots(1,4,figsize=(17,4.3))
for a,(t,d) in zip(ax,[("STELLAR_VEL (DAP, по свету)",stVd),("Ha gas (DAP)",gVd),
                       ("NBursts VEL1 (=газ)",demean(nb1,np.isfinite(nb1))),
                       ("NBursts VEL2",demean(nb2,np.isfinite(nb2)))]):
    im=a.imshow(d,origin="lower",cmap="RdBu_r",vmin=-vmax,vmax=vmax); a.set_title(t,fontsize=10)
    a.set_xticks([]);a.set_yticks([]); fig.colorbar(im,ax=a,fraction=0.046)
fig.suptitle(f"{GAL}: газ vs звёзды напрямую (MaNGA DAP, один грид)",fontsize=13)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(f"verify_gas_vs_stars_{GAL}.png",dpi=115)
print(f"saved verify_gas_vs_stars_{GAL}.png")
