"""Kod olculerine Kalman filtresi (fazla suzme, iyonosfer ayrismasiz).

Her uydu ve her faz yayi icin P1 ve P2 ayri ayri, tek durumlu (skaler) bir Kalman
filtresiyle suzulur. Durum: o frekanstaki kod uzakligi x_k [m].

Olcu modeli (frekans i, gamma_i = (f1/fi)^2, gamma_1 = 1):
    P_i   = rho + gamma_i I + e_P          (kod: gurultulu, cok yollu)
    Phi_i = rho - gamma_i I + B_i + e_Phi  (faz: mm duzeyinde, sabit belirsizlik)
    Phi_1 - Phi_2 = (gamma_2 - 1) I + sabit  ->  dI = d(Phi_1 - Phi_2) / (gamma_2 - 1)

Tahmin (zaman guncellemesi): kodun iki epok arasindaki degisimi fazdan bilinir
    x_k^- = x_{k-1}^+ + dPhi_i + 2 gamma_i dI          (iyonosfer ayrismasi yok)
    P_k^- = P_{k-1}^+ + q_k
Guncelleme (olcu guncellemesi): ham kod olcusu z_k = P_i
    v_k = z_k - x_k^-,  S_k = P_k^- + R,  K_k = P_k^- / S_k
    x_k^+ = x_k^- + K_k v_k,  P_k^+ = (1 - K_k) P_k^-
Yenilik testi: |v_k| > thr sqrt(S_k) ise kod aykiri sayilir, guncelleme atlanir; ust uste 3 ret
olursa (filtre aykiri bir degerden baslamis olabilir) filtre o olcuyle yeniden baslatilir.
Yay basinda (ya da faz kopmasinda) filtre ham kodla yeniden baslar: x_0 = z_0, P_0 = R.

q_k, kazancin kararli durumda K = dt / tau olmasini saglayacak sekilde secilir:
    q = R kappa^2 / (1 - kappa),  kappa = dt / tau
tau buyudukce suzme guclenir (tau -> sonsuz: kumulatif ortalama, klasik Hatch filtresi).
"""
import math

import numpy as np

from .const import C


def _filter_arc(z, dphi, dt, tau, R, thr):
    """Tek yay icin skaler Kalman filtresi. z: ham kod, dphi: fazdan kod degisimi (z ile ayni boy)."""
    n = len(z)
    x = np.empty(n)
    Pp = np.empty(n)
    K = np.empty(n)
    rej = np.zeros(n, bool)
    xk, pk = z[0], R
    x[0], Pp[0], K[0] = xk, pk, 1.0
    n_bad = 0
    n_restart = 0
    for k in range(1, n):
        kap = min(dt[k] / tau, 0.999) if tau > 0 else 0.999
        q = R * kap * kap / (1.0 - kap)
        xm = xk + dphi[k]
        pm = pk + q
        v = z[k] - xm
        s = pm + R
        if abs(v) > thr * math.sqrt(s):
            n_bad += 1
            if n_bad >= 3:
                # ust uste 3 ret: filtre buyuk olasilikla aykiri bir ilk degerden basladi; yeniden baslat
                xk, pk = z[k], R
                K[k] = 1.0
                n_bad = 0
                n_restart += 1
            else:
                rej[k] = True
                xk, pk = xm, pm
                K[k] = 0.0
        else:
            n_bad = 0
            g = pm / s
            xk = xm + g * v
            pk = (1.0 - g) * pm
            K[k] = g
        x[k], Pp[k] = xk, pk
    return x, np.sqrt(Pp), K, rej, n_restart


def kalman_code_filter(P, tau=100.0, sig_code=1.0, thr=5.0):
    """On isleme sonucundaki (P) P1 ve P2'yi suzer.

    Doner: dict(P1f, P2f, K1, K2, s1, s2 [m], rej1, rej2 (bool), cmc1_raw, cmc1_f, stats)
    Faz yayi (arc) olmayan kayitlar degismeden birakilir.
    """
    t, arc = P["t"], P["arc"]
    f1, f2 = P["f1"], P["f2"]
    g2 = (f1 / f2) ** 2
    phi1 = C / f1 * P["L1d"]
    phi2 = C / f2 * P["L2d"]
    I1 = (phi1 - phi2) / (g2 - 1.0)              # f1'deki iyonosfer + sabit [m]
    n = len(t)
    out = {k: np.full(n, np.nan) for k in ("K1", "K2", "s1", "s2")}
    out["P1f"] = P["P1"].copy()
    out["P2f"] = P["P2"].copy()
    out["rej1"] = np.zeros(n, bool)
    out["rej2"] = np.zeros(n, bool)
    out["n_restart"] = 0
    # iyonosfer ayrismasiz kod - faz farki (CMC): kod gurultusu + cok yolluluk + yay sabiti
    cmc1 = P["P1"] - phi1 - 2.0 * I1
    cmc2 = P["P2"] - phi2 - 2.0 * g2 * I1
    R = sig_code ** 2
    for a in np.unique(arc[arc >= 0]):
        m = np.where(arc == a)[0]
        m = m[np.argsort(t[m])]
        if len(m) < 2:
            continue
        dt = np.r_[0.0, np.diff(t[m])]
        dph1 = np.r_[0.0, np.diff(phi1[m] + 2.0 * I1[m])]
        dph2 = np.r_[0.0, np.diff(phi2[m] + 2.0 * g2[m] * I1[m])]
        for zk, dk, pk, kk, sk, rk in (("P1", dph1, "P1f", "K1", "s1", "rej1"), ("P2", dph2, "P2f", "K2", "s2", "rej2")):
            x, s, K, rej, nr = _filter_arc(P[zk][m], dk, dt, tau, R, thr)
            out["n_restart"] += nr
            out[pk][m] = x
            out[sk][m] = s
            out[kk][m] = K
            out[rk][m] = rej
    out["cmc1_raw"] = cmc1
    out["cmc2_raw"] = cmc2
    out["cmc1_f"] = out["P1f"] - phi1 - 2.0 * I1
    out["cmc2_f"] = out["P2f"] - phi2 - 2.0 * g2 * I1
    return out


def noise_stats(P, F, usable, skip_sec):
    """Uydu bazinda ham ve suzulmus kod gurultusu: yay ortalamasi cikarilmis CMC'nin RMS'i.

    Her yayin ilk skip_sec saniyesi (filtrenin oturma suresi) disarida birakilir.
    """
    t, arc, prn = P["t"], P["arc"], P["prn"]
    rows = {}
    for a in np.unique(arc[usable & (arc >= 0)]):
        m = np.where((arc == a) & usable)[0]
        if len(m) < 5:
            continue
        m = m[np.argsort(t[m])]
        m = m[t[m] >= t[m[0]] + skip_sec]
        if len(m) < 5:
            continue
        r = rows.setdefault(int(prn[m[0]]), {"n": 0, "r1": 0.0, "f1": 0.0, "r2": 0.0, "f2": 0.0, "rej": 0})
        for key, arr in (("r1", F["cmc1_raw"]), ("f1", F["cmc1_f"]), ("r2", F["cmc2_raw"]), ("f2", F["cmc2_f"])):
            v = arr[m] - np.mean(arr[m])
            r[key] += float(np.sum(v ** 2))
        r["n"] += len(m)
        r["rej"] += int(F["rej1"][m].sum() + F["rej2"][m].sum())
    out = []
    for p_, r in sorted(rows.items()):
        nn = max(r["n"], 1)
        out.append({"prn": p_, "n": r["n"], "raw1": math.sqrt(r["r1"] / nn), "filt1": math.sqrt(r["f1"] / nn),
                    "raw2": math.sqrt(r["r2"] / nn), "filt2": math.sqrt(r["f2"] / nn), "rej": r["rej"]})
    return out
