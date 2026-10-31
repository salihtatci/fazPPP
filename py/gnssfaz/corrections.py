"""Troposfer, iyonosfer (Klobuchar), anten, gelgit, Shapiro ve faz sarilmasi duzeltmeleri."""
import math

import numpy as np

from .const import C, GM_E, GM_S, GM_M, RE, LAM_NL


# --------------------------------------------------------------------------- troposfer
def standard_met(h):
    """Olcum yoksa standart atmosfer (Berg): P [mbar], T [C], RH [%]."""
    P = 1013.25 * (1 - 2.2557e-5 * h) ** 5.2568
    T = 15.0 - 6.5e-3 * h
    return P, T, 50.0


def saastamoinen(P, Tc, RH, lat, h):
    """Zenit kuru (ZHD) ve islak (ZWD) gecikme [m]."""
    e = RH / 100.0 * 6.1078 * np.exp(17.27 * Tc / (Tc + 237.3))
    Tk = Tc + 273.15
    zhd = 0.0022768 * P / (1 - 0.00266 * math.cos(2 * lat) - 0.00028 * h / 1000.0)
    zwd = 0.002277 * (1255.0 / Tk + 0.05) * e
    return zhd, zwd


def _cf(E, a, b, c):
    sE = np.sin(np.radians(E))
    return (1 + a / (1 + b / (1 + c))) / (sE + a / (sE + b / (sE + c)))


def herring(E, lat, h, Tc):
    """Herring (1992) kuru ve islak esleme fonksiyonlari."""
    cphi = math.cos(lat)
    hk = h / 1000.0
    dT = np.asarray(Tc) + 273.15 - 283.0
    ad = (1.2320 + 0.0139 * cphi - 0.0209 * hk + 0.00215 * dT) * 1e-3
    bd = (3.1612 - 0.1600 * cphi - 0.0331 * hk + 0.00206 * dT) * 1e-3
    cd = (71.244 - 4.293 * cphi - 0.149 * hk - 0.0021 * dT) * 1e-3
    aw = (0.583 - 0.011 * cphi - 0.052 * hk + 0.0014 * dT) * 1e-3
    bw = (1.402 - 0.102 * cphi - 0.101 * hk + 0.0020 * dT) * 1e-3
    cw = (45.85 - 1.91 * cphi - 1.29 * hk + 0.015 * dT) * 1e-3
    return _cf(E, ad, bd, cd), _cf(E, aw, bw, cw)


# --------------------------------------------------------------------------- Klobuchar
def klobuchar(t_gps, lat, lon, el, az, alpha, beta):
    """Klobuchar L1 iyonosfer gecikmesi [m]. lat/lon [rad], el/az [derece]."""
    phi_u = lat / math.pi
    lam_u = lon / math.pi
    E = np.radians(el) / math.pi
    A = np.radians(az)
    psi = 0.0137 / (E + 0.11) - 0.022
    phi_i = np.clip(phi_u + psi * np.cos(A), -0.416, 0.416)
    lam_i = lam_u + psi * np.sin(A) / np.cos(phi_i * math.pi)
    phi_m = phi_i + 0.064 * np.cos((lam_i - 1.617) * math.pi)
    tl = (4.32e4 * lam_i + np.asarray(t_gps) % 86400.0) % 86400.0
    F = 1.0 + 16.0 * (0.53 - E) ** 3
    amp = sum(alpha[k] * phi_m ** k for k in range(4))
    per = sum(beta[k] * phi_m ** k for k in range(4))
    amp = np.maximum(amp, 0.0)
    per = np.maximum(per, 72000.0)
    x = 2 * math.pi * (tl - 50400.0) / per
    d = np.where(np.abs(x) < 1.57, F * (5e-9 + amp * (1 - x ** 2 / 2 + x ** 4 / 24)), F * 5e-9)
    return C * d


# --------------------------------------------------------------------------- anten
def pcv_interp(table, zen0, dzen, angle):
    zz = zen0 + dzen * np.arange(len(table))
    return np.interp(angle, zz, table)


def sat_body_axes(Xs, sun):
    ez = -Xs / np.linalg.norm(Xs, axis=1)[:, None]
    es = sun - Xs
    es /= np.linalg.norm(es, axis=1)[:, None]
    ey = np.cross(ez, es)
    ey /= np.linalg.norm(ey, axis=1)[:, None]
    ex = np.cross(ey, ez)
    return ex, ey, ez


def receiver_antenna(u, el, rec, arp_hen, Rneu, coef):
    """Alici ARP + PCO + PCV uzaklik duzeltmesi [m] (model uzakliga eklenir).

    coef: {'G01': a, 'G02': b} frekans kombinasyonu katsayilari (L3 icin a3, b3).
    """
    E_, N_, U_ = Rneu[1], Rneu[0], Rneu[2]
    h, e_, n_ = arp_hen
    off = h * U_ + e_ * E_ + n_ * N_
    pcv = np.zeros(len(u))
    if rec is not None:
        pco = sum(coef[f] * rec[f]["pco"] for f in coef) / 1000.0        # (N, E, U)
        off = off + pco[0] * N_ + pco[1] * E_ + pco[2] * U_
        zen = 90.0 - el
        if all("pcv" in rec[f] for f in coef):
            pcv = sum(coef[f] * pcv_interp(rec[f]["pcv"], rec["zen"][0], rec["zen"][1], zen) for f in coef) / 1000.0
    return -(u @ off) + pcv


def satellite_antenna(u, prns, ex, ey, ez, sats, coef):
    """Uydu PCO (govde ekseni) + PCV(nadir) uzaklik duzeltmesi [m]; nadir acisi [derece]."""
    n = len(u)
    d = np.zeros(n)
    nadir = np.degrees(np.arccos(np.clip(-np.einsum("ij,ij->i", u, ez), -1, 1)))
    missing = set()
    for i in range(n):
        s = sats.get(int(prns[i]))
        if s is None:
            missing.add(int(prns[i]))
            continue
        p = sum(coef[f] * s[f]["pco"] for f in coef) / 1000.0
        kay = p[0] * ex[i] + p[1] * ey[i] + p[2] * ez[i]
        pv = 0.0
        if all("pcv" in s[f] for f in coef):
            pv = sum(coef[f] * pcv_interp(s[f]["pcv"], s["zen"][0], s["zen"][1], nadir[i]) for f in coef) / 1000.0
        d[i] = u[i] @ kay + pv
    return d, nadir, missing


# --------------------------------------------------------------------------- gelgit, Shapiro
def solid_tide(X0, sun, moon):
    """IERS 2010 ikinci derece kati yer gelgiti yer degistirmesi [m] (n,3)."""
    rhat = X0 / np.linalg.norm(X0)
    dr = np.zeros_like(sun)
    for cis, GM in ((sun, GM_S), (moon, GM_M)):
        R = np.linalg.norm(cis, axis=1)
        Rh = cis / R[:, None]
        sc = Rh @ rhat
        k = GM / GM_E * RE ** 4 / R ** 3
        dr += k[:, None] * (0.6078 * rhat[None, :] * (1.5 * sc ** 2 - 0.5)[:, None]
                            + 3 * 0.0847 * sc[:, None] * (Rh - sc[:, None] * rhat[None, :]))
    return dr


def shapiro(Xs, X0, rho):
    rs = np.linalg.norm(Xs, axis=1)
    rr = np.linalg.norm(X0)
    return 2 * GM_E / C ** 2 * np.log((rs + rr + rho) / (rs + rr - rho))


# --------------------------------------------------------------------------- faz sarilmasi
def windup_raw(u, ex, ey, Rneu):
    """Wu vd. (1993) faz sarilmasi [devir], [-0.5, 0.5] araliginda (sureklilestirilmemis)."""
    N_, E_ = Rneu[0], Rneu[1]
    k = -u
    exr, eyr = N_, -E_
    ds = ex - k * np.einsum("ij,ij->i", k, ex)[:, None] - np.cross(k, ey)
    drr = exr[None, :] - k * (k @ exr)[:, None] + np.cross(k, np.broadcast_to(eyr, k.shape))
    cosp = np.clip(np.einsum("ij,ij->i", ds, drr) / np.linalg.norm(ds, axis=1) / np.linalg.norm(drr, axis=1), -1, 1)
    ph = np.arccos(cosp) / (2 * math.pi)
    sgn = np.einsum("ij,ij->i", k, np.cross(ds, drr))
    return np.where(sgn < 0, -ph, ph)


def unwrap_by_arc(raw, arc, order):
    """Faz sarilmasini her yay boyunca surekli yap, tam sayi kismini sifirla. order: kayit sirasi."""
    out = np.zeros_like(raw)
    for a in np.unique(arc):
        if a < 0:
            continue
        m = np.where(arc == a)[0]
        m = m[np.argsort(order[m])]
        w = np.unwrap(raw[m] * 2 * math.pi) / (2 * math.pi)
        w -= np.round(w[0])
        out[m] = w
    return out


L3_WINDUP_M = LAM_NL
