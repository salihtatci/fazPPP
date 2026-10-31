"""On isleme: donmus kayit, alici saati sicramalari, GF/MW ile faz kopmasi ve yaylar."""
import math

import numpy as np

from .const import C, F1, F2, LAM1, LAM2, LAM_W, MS


def preprocess(obs, gf_thr=None, mw_thr=4.0, min_arc_sec=600.0):
    t, prn = obs["t"], obs["prn"]
    P1, P2, L1, L2 = (obs[k].copy() for k in ("P1", "P2", "L1", "L2"))
    interval = float(obs["header"]["interval"])
    n0 = len(t)
    order = np.lexsort((t, prn))                       # uydu, sonra zaman
    t, prn, P1, P2, L1, L2 = (a[order] for a in (t, prn, P1, P2, L1, L2))
    C1 = obs["C1"][order]
    same = np.r_[False, prn[1:] == prn[:-1]]

    # 1) donmus L2: bir onceki kayitla birebir ayni deger
    frozen = same & (np.r_[np.nan, np.diff(L2)] == 0.0)
    P2[frozen] = np.nan
    L2[frozen] = np.nan
    n_frozen = int(frozen.sum())

    full = ~(np.isnan(P1) | np.isnan(P2) | np.isnan(L1) | np.isnan(L2))
    n_missing = int((~full).sum())
    t, prn, P1, P2, L1, L2, C1 = (a[full] for a in (t, prn, P1, P2, L1, L2, C1))
    same = np.r_[False, prn[1:] == prn[:-1]]
    dt = np.r_[np.nan, np.diff(t)]

    # 2) alici saati sicramalari: tum uydularda ayni anda d(P1 - lam1 L1) ~ k * 299 792 m
    dcmc = np.r_[np.nan, np.diff(P1) - LAM1 * np.diff(L1)]
    cons = same & (np.abs(dt - interval) < 0.5 * interval)
    ep = np.unique(t)
    jumps = []
    k_ep = np.zeros(len(ep))
    if cons.any():
        tt, vv = t[cons], dcmc[cons]
        o2 = np.argsort(tt, kind="stable")
        tt, vv = tt[o2], vv[o2]
        cut = np.r_[0, np.nonzero(np.diff(tt))[0] + 1, len(tt)]
        for a, b in zip(cut[:-1], cut[1:]):
            if b - a < 2:
                continue
            med = float(np.median(vv[a:b]))
            k = int(round(med / MS))
            if k != 0 and abs(med - k * MS) < 1000.0:
                j = int(np.searchsorted(ep, tt[a]))
                k_ep[j] = k
                jumps.append({"t": float(tt[a]), "median_m": med, "nsat": int(b - a), "k_ms": k})
    K = np.cumsum(k_ep)
    Krec = K[np.searchsorted(ep, t)]
    L1d = L1 + Krec * F1 * 1e-3
    L2d = L2 + Krec * F2 * 1e-3

    # 3) GF ve MW
    GF = LAM1 * L1d - LAM2 * L2d
    MW = ((F1 * LAM1 * L1d - F2 * LAM2 * L2d) / (F1 - F2) - (F1 * P1 + F2 * P2) / (F1 + F2)) / LAM_W
    if gf_thr is None:
        gf_thr = 0.10 if interval <= 30 else min(0.5, 0.10 * interval / 30.0)

    # 4) yaylar
    n = len(t)
    arc = np.full(n, -1, dtype=int)
    outlier = np.zeros(n, dtype=bool)
    arc_no = -1
    prev = None                          # (t, GF)
    mw_sum = 0.0
    mw_n = 0
    gap = 1.5 * interval
    t_last = -1e99
    for i in range(n):
        new = False
        if prev is None or prn[i] != prev[2] or t[i] - t_last > gap:
            new = True
        else:
            ok = abs(GF[i] - prev[1]) <= gf_thr and abs(MW[i] - mw_sum / mw_n) <= mw_thr
            if not ok:
                j = i + 1
                if j < n and prn[j] == prn[i] and t[j] - prev[0] <= 2 * gap and \
                        abs(GF[j] - prev[1]) <= gf_thr and abs(MW[j] - mw_sum / mw_n) <= mw_thr:
                    outlier[i] = True
                    t_last = t[i]
                    continue
                new = True
        if new:
            arc_no += 1
            mw_sum, mw_n = 0.0, 0
        arc[i] = arc_no
        mw_sum += MW[i]
        mw_n += 1
        prev = (t[i], GF[i], prn[i])
        t_last = t[i]

    # yay istatistikleri ve kullanilabilirlik
    min_n = max(3, int(math.ceil(min_arc_sec / interval)))
    arcs = []
    usable_arc = np.zeros(arc_no + 1, dtype=bool)
    for a in range(arc_no + 1):
        m = arc == a
        cnt = int(m.sum())
        if cnt == 0:
            continue
        tm = t[m]
        mwv = MW[m]
        use = cnt >= min_n
        usable_arc[a] = use
        arcs.append({"arc": a, "prn": int(prn[m][0]), "start": float(tm[0]), "end": float(tm[-1]),
                     "n": cnt, "mw_std": float(np.std(mwv)) if cnt > 1 else None, "usable": bool(use)})
    usable = np.zeros(n, dtype=bool)
    good = arc >= 0
    usable[good] = usable_arc[arc[good]]
    usable &= ~outlier

    return {"t": t, "prn": prn, "P1": P1, "P2": P2, "C1": C1, "L1": L1, "L2": L2, "L1d": L1d, "L2d": L2d,
            "GF": GF, "MW": MW, "arc": arc, "outlier": outlier, "usable": usable,
            "jumps": jumps, "K": Krec, "arcs": arcs, "interval": interval,
            "stats": {"n_records": n0, "n_frozen": n_frozen, "n_incomplete": n_missing,
                      "n_complete": n, "n_arcs": len(arcs), "n_usable_arcs": int(usable_arc.sum()),
                      "n_outlier": int(outlier.sum()), "n_usable": int(usable.sum()),
                      "gf_thr": gf_thr, "mw_thr": mw_thr, "min_arc_epochs": min_n}}
