"""Uydu konum ve saati: SP3 Lagrange ara degerlemesi ve GPS yayin efemerisi."""
import math

import numpy as np

from .const import C, GPS_MU, GPS_OMEGA, F_REL

NPT = 11          # Lagrange nokta sayisi (10. derece)


def _lagrange_weights(tq, tt):
    """tq (m,), tt (m, NPT) -> agirliklar w (m, NPT) ve turev agirliklari dw (m, NPT)."""
    m, n = tt.shape
    w = np.ones((m, n))
    dw = np.zeros((m, n))
    for j in range(n):
        for k in range(n):
            if k == j:
                continue
            w[:, j] *= (tq - tt[:, k]) / (tt[:, j] - tt[:, k])
        # turev: l_j' = sum_k [ prod_{i!=j,k} (t-t_i)/(t_j-t_i) ] / (t_j - t_k)
        s = np.zeros(m)
        for k in range(n):
            if k == j:
                continue
            p = np.ones(m) / (tt[:, j] - tt[:, k])
            for i in range(n):
                if i in (j, k):
                    continue
                p *= (tq - tt[:, i]) / (tt[:, j] - tt[:, i])
            s += p
        dw[:, j] = s
    return w, dw


def sp3_interp(sp3, prn, tq):
    """Uydu konumu [m] ve hizi [m/s], tq (GPS sn) anlarinda. SP3 araligi disi -> NaN."""
    tq = np.atleast_1d(np.asarray(tq, float))
    pos = np.full((len(tq), 3), np.nan)
    vel = np.full((len(tq), 3), np.nan)
    if prn not in sp3["prns"]:
        return pos, vel
    j = sp3["prns"].index(prn)
    ts = sp3["t"]
    xyz = sp3["xyz"][:, j, :]
    n = len(ts)
    if n < NPT:
        return pos, vel
    ok = (tq >= ts[0]) & (tq <= ts[-1])
    if not ok.any():
        return pos, vel
    q = tq[ok]
    s0 = np.clip(np.searchsorted(ts, q) - NPT // 2, 0, n - NPT)
    idx = s0[:, None] + np.arange(NPT)[None, :]
    tt = ts[idx]
    w, dw = _lagrange_weights(q, tt)
    P = xyz[idx]                                  # (m, NPT, 3)
    p = np.einsum("mn,mnk->mk", w, P)
    v = np.einsum("mn,mnk->mk", dw, P)
    bad = np.isnan(P).any(axis=(1, 2))
    p[bad] = np.nan
    v[bad] = np.nan
    pos[ok] = p
    vel[ok] = v
    return pos, vel


def sp3_clock(sp3, prn, tq):
    """Uydu saat hatasi [s], SP3 saatinin dogrusal ara degeri (komsulardan biri yoksa NaN)."""
    tq = np.atleast_1d(np.asarray(tq, float))
    out = np.full(len(tq), np.nan)
    if prn not in sp3["prns"]:
        return out
    j = sp3["prns"].index(prn)
    ts = sp3["t"]
    c = sp3["clk"][:, j]
    k = np.clip(np.searchsorted(ts, tq) - 1, 0, len(ts) - 2)
    t0, t1 = ts[k], ts[k + 1]
    a = (tq - t0) / (t1 - t0)
    v = c[k] * (1 - a) + c[k + 1] * a
    inside = (tq >= ts[0]) & (tq <= ts[-1])
    out[inside] = v[inside]
    return out


# ---------------------------------------------------------------------------
# Yayin efemerisi (IS-GPS-200)
# ---------------------------------------------------------------------------
def select_eph(ephs, prn, t):
    best, bd = None, 1e99
    for e in ephs:
        if e["prn"] != prn or (e.get("health") or 0) != 0:
            continue
        d = abs(t - e["toe_abs"])
        if d < bd:
            best, bd = e, d
    if best is None or bd > 4 * 3600:
        return None
    return best


def bc_clock(e, t):
    dt = t - e["toc"]
    return e["af0"] + e["af1"] * dt + e["af2"] * dt * dt


def bc_position(e, t):
    """Yayin efemerisinden ECEF konum [m] ve rolativistik saat duzeltmesi [s] (t: gonderim ani)."""
    A = e["sqrtA"] ** 2
    n0 = math.sqrt(GPS_MU / A ** 3)
    tk = t - e["toe_abs"]
    n = n0 + e["dn"]
    M = e["M0"] + n * tk
    E = M
    for _ in range(10):
        E = M + e["e"] * math.sin(E)
    nu = math.atan2(math.sqrt(1 - e["e"] ** 2) * math.sin(E), math.cos(E) - e["e"])
    phi = nu + e["omega"]
    du = e["Cus"] * math.sin(2 * phi) + e["Cuc"] * math.cos(2 * phi)
    dr = e["Crs"] * math.sin(2 * phi) + e["Crc"] * math.cos(2 * phi)
    di = e["Cis"] * math.sin(2 * phi) + e["Cic"] * math.cos(2 * phi)
    u = phi + du
    r = A * (1 - e["e"] * math.cos(E)) + dr
    i = e["i0"] + di + e["IDOT"] * tk
    xp, yp = r * math.cos(u), r * math.sin(u)
    OM = e["OMEGA0"] + (e["OMEGADOT"] - GPS_OMEGA) * tk - GPS_OMEGA * e["Toe"]
    x = xp * math.cos(OM) - yp * math.cos(i) * math.sin(OM)
    y = xp * math.sin(OM) + yp * math.cos(i) * math.cos(OM)
    z = yp * math.sin(i)
    drel = F_REL * e["e"] * e["sqrtA"] * math.sin(E)
    return np.array([x, y, z]), drel
