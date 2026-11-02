"""On isleme: donmus kayit, alici saati sicramalari, GF/MW ile faz kopmasi ve yaylar.

GPS (L1/L2) ve Galileo (E1/E5a ya da E5b) kayitlari birlikte islenir; her kaydin
frekanslari uydu sistemine gore secilir.
"""
import math

import numpy as np

from .const import C, F1, F2, GAL, MS


def freqs(prn, f2_gal):
    """Kayit bazinda (f1, f2) [Hz]."""
    isE = np.asarray(prn) >= GAL
    f1 = np.full(len(isE), F1)
    f2 = np.where(isE, f2_gal if f2_gal else np.nan, F2)
    return f1, f2


def preprocess(obs, gf_thr=None, mw_thr=4.0, min_arc_sec=600.0, max_gap=None):
    t, prn = obs["t"], obs["prn"]
    P1, P2, L1, L2 = (obs[k].copy() for k in ("P1", "P2", "L1", "L2"))
    interval = float(obs["header"]["interval"])
    f2_gal = (obs["header"].get("f2") or {}).get("E")
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
    f1, f2 = freqs(prn, f2_gal)
    lam1, lam2 = C / f1, C / f2
    lamw = C / (f1 - f2)

    # 2) alici saati sicramalari: tum uydularda ayni anda d(P1 - lam1 L1) ~ k * 299 792 m
    dcmc = np.r_[np.nan, np.diff(P1) - lam1[1:] * np.diff(L1)]
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
    L1d = L1 + Krec * f1 * 1e-3
    L2d = L2 + Krec * f2 * 1e-3

    # 3) GF ve MW
    GF = lam1 * L1d - lam2 * L2d
    MW = ((f1 * lam1 * L1d - f2 * lam2 * L2d) / (f1 - f2) - (f1 * P1 + f2 * P2) / (f1 + f2)) / lamw
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
    # yuksek hizli veride (ornegin 1 s) birkac epokluk kesinti yeni yay acmaz;
    # kopma varsa GF/MW testi yakalar
    gap = max(1.5 * interval, 30.0) if max_gap is None else max(1.5 * interval, float(max_gap))
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

    # 4b) yuksek hizli veride (<= 5 s) pencere testi: epoktan epoga farkta gurultuye gomulen kucuk
    # kopmalar (ornegin L1 ve L2'de esit devir: GF'de birkac cm, MW'de sifir) icin yayin her
    # noktasinda oncesi ve sonrasi W kayda dogru uydurulur, sinirdaki degerler karsilastirilir.
    # Dogru uydurma, yavas degisimleri (iyonosfer, alcak uydularda yurume sirasinda cok yolluluk)
    # basamak sanmayi onler.
    n_win = 0
    if interval <= 5.0 and arc_no >= 0:
        W = max(5, int(round(20.0 / interval)))
        thr_gf, thr_mw = 0.04, 2.0
        next_arc = arc_no + 1
        xb = (W - 1) / 2.0
        sxx = W * (W * W - 1) / 12.0

        def edge_fits(y, k):
            cy = np.r_[0.0, np.cumsum(y)]
            jj = np.arange(len(y), dtype=float)
            cjy = np.r_[0.0, np.cumsum(jj * y)]
            out = []
            for s0, x0 in ((k - W, W - 0.5), (k, -0.5)):
                sy = cy[s0 + W] - cy[s0]
                sxy = (cjy[s0 + W] - cjy[s0]) - s0 * sy
                ym = sy / W
                b = (sxy - xb * sy) / sxx
                out.append(ym + b * (x0 - xb))
            return out[1] - out[0]

        for a in range(arc_no + 1):
            m = np.where(arc == a)[0]
            if len(m) < 2 * W + 1:
                continue
            tt = t[m]
            k = np.arange(W, len(m) - W + 1)
            dg = edge_fits(GF[m], k)
            dw = edge_fits(MW[m], k)
            span_ok = (tt[np.minimum(k + W - 1, len(m) - 1)] - tt[k - W]) <= 3.0 * W * interval
            score = np.maximum(np.abs(dg) / thr_gf, np.abs(dw) / thr_mw) * span_ok
            cand = np.where(score > 1.0)[0]
            if not len(cand):
                continue
            cuts = []
            grp = [cand[0]]
            for c in cand[1:]:
                if c - grp[-1] <= W:
                    grp.append(c)
                else:
                    cuts.append(grp[int(np.argmax(score[grp]))])
                    grp = [c]
            cuts.append(grp[int(np.argmax(score[grp]))])
            for c in cuts:
                j0 = int(k[c])
                arc[m[j0:]] = next_arc
                next_arc += 1
                n_win += 1
        arc_no = next_arc - 1

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
            "f1": f1, "f2": f2, "isE": prn >= GAL,
            "GF": GF, "MW": MW, "arc": arc, "outlier": outlier, "usable": usable,
            "jumps": jumps, "K": Krec, "arcs": arcs, "interval": interval,
            "stats": {"n_records": n0, "n_frozen": n_frozen, "n_incomplete": n_missing,
                      "n_complete": n, "n_arcs": len(arcs), "n_usable_arcs": int(usable_arc.sum()),
                      "n_outlier": int(outlier.sum()), "n_usable": int(usable.sum()), "n_win_slips": n_win,
                      "gf_thr": gf_thr, "mw_thr": mw_thr, "min_arc_epochs": min_n}}
