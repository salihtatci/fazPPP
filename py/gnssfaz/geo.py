"""Koordinat donusumleri, yukseklik/azimut ve dusuk dogruluklu Gunes/Ay konumu."""
import math

import numpy as np

from .const import WGS_A, WGS_F


def xyz2llh(x, y, z):
    """ECEF -> (enlem [rad], boylam [rad], h [m]) GRS80/WGS84 (Bowring + yineleme)."""
    a, f = WGS_A, WGS_F
    e2 = f * (2 - f)
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - e2))
    h = 0.0
    for _ in range(10):
        N = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return lat, lon, h


def rot_neu(lat, lon):
    """XYZ farkini (N, E, U)'ya ceviren 3x3 matris."""
    sl, cl, sb, cb = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    return np.array([[-sl * cb, -sl * sb, cl],
                     [-sb, cb, 0.0],
                     [cl * cb, cl * sb, sl]])


def elaz(X0, Xs):
    """Alicidan (X0) uydulara (Xs, n x 3) yukseklik ve azimut [derece]."""
    lat, lon, _ = xyz2llh(*X0)
    R = rot_neu(lat, lon)
    d = Xs - X0
    neu = d @ R.T
    rng = np.linalg.norm(d, axis=1)
    el = np.degrees(np.arcsin(np.clip(neu[:, 2] / rng, -1, 1)))
    az = np.degrees(np.arctan2(neu[:, 1], neu[:, 0])) % 360
    return el, az


def dms(rad):
    d = abs(math.degrees(rad))
    g = int(d)
    m = int((d - g) * 60)
    s = (d - g - m / 60) * 3600
    return f"{'-' if rad < 0 else ''}{g}° {m:02d}' {s:08.5f}\""


def sun_moon_ecef(t_gps):
    """Gunes ve Ay'in ECEF konumu [m], Montenbruck & Gill (2000) dusuk dogruluk formulleri.

    t_gps: GPS saniyesi (dizi). Doner: (sun (n,3), moon (n,3))
    """
    t = np.atleast_1d(np.asarray(t_gps, float))
    jd = 2444244.5 + t / 86400.0                          # GPS zaman olcegi
    T = (jd + 51.184 / 86400.0 - 2451545.0) / 36525.0   # TT = GPS + 19 s + 32.184 s
    eps = math.radians(23.43929111)
    asn = math.radians(1 / 3600)
    M = np.radians(357.5256 + 35999.049 * T)
    Ls = math.radians(282.9400) + M + (6892 * np.sin(M) + 72 * np.sin(2 * M)) * asn
    r = 149.619e9 - 2.499e9 * np.cos(M) - 0.021e9 * np.cos(2 * M)
    sun = np.stack([r * np.cos(Ls), r * np.sin(Ls) * math.cos(eps), r * np.sin(Ls) * math.sin(eps)], -1)
    L0 = np.radians(218.31617 + 481267.88088 * T - 1.3972 * T)
    l = np.radians(134.96292 + 477198.86753 * T)
    lp = np.radians(357.52543 + 35999.04944 * T)
    F = np.radians(93.27283 + 483202.01873 * T)
    D = np.radians(297.85027 + 445267.11135 * T)
    lam = L0 + asn * (22640 * np.sin(l) + 769 * np.sin(2 * l) - 4586 * np.sin(l - 2 * D)
                      + 2370 * np.sin(2 * D) - 668 * np.sin(lp) - 412 * np.sin(2 * F)
                      - 212 * np.sin(2 * l - 2 * D) - 206 * np.sin(l + lp - 2 * D)
                      + 192 * np.sin(l + 2 * D) - 165 * np.sin(lp - 2 * D) + 148 * np.sin(l - lp)
                      - 125 * np.sin(D) - 110 * np.sin(l + lp) - 55 * np.sin(2 * F - 2 * D))
    bet = asn * (18520 * np.sin(F + lam - L0 + asn * (412 * np.sin(2 * F) + 541 * np.sin(lp)))
                 - 526 * np.sin(F - 2 * D) + 44 * np.sin(l + F - 2 * D) - 31 * np.sin(-l + F - 2 * D)
                 - 25 * np.sin(-2 * l + F) - 23 * np.sin(lp + F - 2 * D) + 21 * np.sin(-l + F)
                 + 11 * np.sin(-lp + F - 2 * D))
    rm = (385000 - 20905 * np.cos(l) - 3699 * np.cos(2 * D - l) - 2956 * np.cos(2 * D)
          - 570 * np.cos(2 * l) + 246 * np.cos(2 * l - 2 * D) - 205 * np.cos(lp - 2 * D)
          - 171 * np.cos(l + 2 * D) - 152 * np.cos(l + lp - 2 * D)) * 1e3
    ekl = np.stack([rm * np.cos(bet) * np.cos(lam), rm * np.cos(bet) * np.sin(lam), rm * np.sin(bet)], -1)
    moon = np.stack([ekl[:, 0], math.cos(eps) * ekl[:, 1] - math.sin(eps) * ekl[:, 2],
                     math.sin(eps) * ekl[:, 1] + math.cos(eps) * ekl[:, 2]], -1)
    jd_ut = jd - 18.0 / 86400.0                           # UT1 ~ GPS - 18 s (2017+); +-2 s fark onemsiz
    gmst = np.radians((280.46061837 + 360.98564736629 * (jd_ut - 2451545.0)) % 360)
    cg, sg = np.cos(gmst), np.sin(gmst)

    def rot(v):
        return np.stack([cg * v[:, 0] + sg * v[:, 1], -sg * v[:, 0] + cg * v[:, 1], v[:, 2]], -1)
    return rot(sun), rot(moon)
