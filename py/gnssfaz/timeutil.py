"""Zaman donusumleri. Butun ic hesaplar GPS saniyesi (1980-01-06'dan beri) ile yapilir."""
from datetime import datetime, timedelta

GPS0 = datetime(1980, 1, 6)


def to_gps(y, mo, d, h=0, mi=0, s=0.0):
    """Takvim tarihi -> GPS saniyesi (GPS zaman olcegi)."""
    if y < 100:
        y += 2000 if y < 80 else 1900
    dt = datetime(int(y), int(mo), int(d), int(h), int(mi)) - GPS0
    return dt.days * 86400.0 + dt.seconds + float(s)


def from_gps(t):
    return GPS0 + timedelta(seconds=float(t))


def fmt(t, with_date=False):
    d = from_gps(round(t * 1000) / 1000)
    return d.strftime("%Y-%m-%d %H:%M:%S") if with_date else d.strftime("%H:%M:%S")


def day_start(t):
    return t - (t % 86400.0)


def julian_date(t_gps):
    """GPS saniyesi -> Julyen gunu (GPS zaman olcegi)."""
    return 2444244.5 + t_gps / 86400.0
