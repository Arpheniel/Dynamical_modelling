"""
Газовая металличность по эмиссиям (опциональный дискриминатор §4.3).

ВАЖНО: strong-line калибровки (O3N2, N2) применимы ТОЛЬКО к газу,
ионизованному звёздами (HII / star-forming). У S0 с контрвращением газ часто
ионизован ударно/LINER/AGN — тогда металличность по сильным линиям не имеет
смысла. Поэтому: сначала BPT, и металличность считаем лишь для SF-бинов.

Линии берём из EML_PARAM (chem[2]): FLUX/FLUX_ERR на бин, матч по LINE_ID.
Запускать из рабочей папки галактики.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import chemo_unified as cu

GAL = cu.GALAXY
eml = cu.chem[2].data
LID = np.array(eml['LINE_ID'])     # (Nbin, 66) names
FLUX = np.array(eml['FLUX'], float) # (Nbin, 66)
FERR = np.array(eml['FLUX_ERR'], float)
FLAG = np.array(eml['FLAG'], int)
nbin = FLUX.shape[0]

# имена линий из строки 0 (шаблон фиксирован по бинам)
names0 = [str(x).strip() if not isinstance(x,(bytes,bytearray)) else x.decode('latin1','ignore').strip()
          for x in LID[0]]

def line_idx(substr):
    for j,nm in enumerate(names0):
        if substr.lower() in nm.lower():
            return j
    return None

want = {
    'Ha':   line_idx('H alpha'),
    'Hb':   line_idx('H beta'),
    'O3':   line_idx('[O III] 5008'),
    'N2':   line_idx('[N II] 6585'),
    'S2a':  line_idx('[S II] 6718'),
    'S2b':  line_idx('[S II] 6732'),
    'O2a':  line_idx('[O II] 3727'),
    'O2b':  line_idx('[O II] 3729'),
}
print(f"\n===== ГАЗ-МЕТАЛЛИЧНОСТЬ: {GAL} =====")
print("индексы линий:", want)

def get(key):
    j = want[key]
    f = FLUX[:, j].copy(); e = FERR[:, j].copy(); fl = FLAG[:, j]
    bad = ~np.isfinite(f) | ~np.isfinite(e) | (e <= 0) | (fl < 0)
    f[bad] = np.nan; e[bad] = np.nan
    return f, e

Ha,eHa = get('Ha'); Hb,eHb = get('Hb')
O3,eO3 = get('O3'); N2,eN2 = get('N2')

# S/N на бин
def snr(f,e):
    with np.errstate(invalid='ignore',divide='ignore'):
        return f/e

# Полный поток по галактике (сумма по бинам с положит. потоком и S/N>2)
def total(f,e,thr=2.0):
    m = np.isfinite(f)&np.isfinite(e)&(f>0)&(snr(f,e)>thr)
    return np.nansum(f[m]), int(m.sum())

print("\nИнтегральные потоки (бины с S/N>2):")
for k,(f,e) in [('Ha',(Ha,eHa)),('Hb',(Hb,eHb)),('[OIII]5007',(O3,eO3)),('[NII]6584',(N2,eN2))]:
    t,n = total(f,e)
    print(f"  {k:12s} sum={t:12.3g}  N_bins(S/N>2)={n}")

# ---- BPT по бинам (нужны 4 линии с S/N>3) ----
SN = 3.0
good = (snr(Ha,eHa)>SN)&(snr(Hb,eHb)>SN)&(snr(O3,eO3)>SN)&(snr(N2,eN2)>SN)
good &= np.isfinite(Ha)&np.isfinite(Hb)&np.isfinite(O3)&np.isfinite(N2)
ng = int(np.nansum(good))
print(f"\nБинов с 4 линиями S/N>{SN:.0f}: {ng} из {nbin}")

if ng >= 3:
    x = np.log10(N2[good]/Ha[good])    # log [NII]/Ha
    y = np.log10(O3[good]/Hb[good])    # log [OIII]/Hb
    # демаркации
    def kauff(xx): return 0.61/(xx-0.05)+1.30   # Kauffmann03 SF
    def kewley(xx): return 0.61/(xx-0.47)+1.19  # Kewley01 max starburst
    cls = []
    for xi,yi in zip(x,y):
        kf = kauff(xi) if xi<0.05 else -np.inf
        kw = kewley(xi) if xi<0.47 else -np.inf
        if yi < kf: cls.append('SF')
        elif yi < kw: cls.append('COMP')
        else:
            # Seyfert/LINER split (Schawinski+07 в [NII]-BPT приближённо по [OIII]/Hb)
            cls.append('LINER' if yi < 0.95*xi+0.56 else 'SEY')  # грубо
    cls = np.array(cls)
    import collections
    cnt = collections.Counter(cls)
    print("BPT-классы бинов:", dict(cnt))
    print(f"  медиана log[NII]/Ha = {np.median(x):+.2f}, log[OIII]/Hb = {np.median(y):+.2f}")

    sf = cls=='SF'
    print(f"\nSF-бинов: {int(sf.sum())}")
    if sf.sum() >= 3:
        O3N2 = np.log10((O3[good]/Hb[good])/(N2[good]/Ha[good]))[sf]
        N2r  = np.log10((N2[good]/Ha[good]))[sf]
        OH_O3N2 = 8.533 - 0.214*O3N2     # Marino+2013
        OH_N2   = 8.743 + 0.462*N2r      # Marino+2013
        for nm,oh in [('O3N2',OH_O3N2),('N2',OH_N2)]:
            moh = np.nanmedian(oh)
            print(f"  12+log(O/H) [{nm}] median = {moh:.2f}  -> [O/H] = {moh-8.69:+.2f} (солн.8.69)")
    else:
        print("  SF-бинов мало -> strong-line металличность НЕ считаем (газ не звёздно-ионизован).")

    np.savez(f"bpt_points_{GAL}.npz", x=x, y=y, cls=np.array([str(c) for c in cls]))
    print(f"  [сохранено] bpt_points_{GAL}.npz")
else:
    print("Эмиссии слишком слабы/мало для BPT.")
