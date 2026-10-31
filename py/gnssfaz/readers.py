"""RINEX (gozlem, navigasyon, meteoroloji), SP3, ANTEX ve DCB okuyuculari.

Yalnizca GPS kayitlari okunur. RINEX 2.xx ve 3.xx desteklenir; Hatanaka
sikistirmali (.d / .crx) dosyalar desteklenmez.
"""
import math
from datetime import datetime

import numpy as np

from .timeutil import to_gps


class FormatHatasi(ValueError):
    pass


def _f(s):
    s = s.strip().replace("D", "E").replace("d", "e")
    if not s:
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def _lines(text):
    return text.replace("\r\n", "\n").replace("\r", "\n").split("\n")


# ---------------------------------------------------------------------------
# RINEX gozlem dosyasi
# ---------------------------------------------------------------------------
_V3_PREF = {
    "P1": ["C1W", "C1P", "C1Y", "C1C", "C1X"],
    "C1": ["C1C", "C1X"],
    "P2": ["C2W", "C2P", "C2Y", "C2D", "C2X", "C2L", "C2S"],
    "L1": ["L1C", "L1W", "L1P", "L1X"],
    "L2": ["L2W", "L2P", "L2Y", "L2D", "L2X", "L2L", "L2S"],
}


def read_obs(text):
    """RINEX gozlem dosyasini okur.

    Doner: dict(header=..., t, prn, P1, P2, L1, L2, C1) - her kayit bir (epok, uydu).
    """
    L = _lines(text)
    if L and ("COMPACT RINEX" in L[0].upper() or "CRINEX" in L[0].upper()):
        raise FormatHatasi("obs_hatanaka|Hatanaka sikistirmali dosya: once CRX2RNX ile acin.")
    if not L or "RINEX VERSION" not in L[0]:
        raise FormatHatasi("obs_not_rinex|Gozlem dosyasi RINEX degil (ilk satirda 'RINEX VERSION / TYPE' yok).")
    ver = _f(L[0][:9])
    if "COMPACT" in L[0].upper() or "CRINEX" in L[0].upper():
        raise FormatHatasi("obs_hatanaka|Hatanaka sikistirmali dosya: once CRX2RNX ile acin.")
    if L[0][20:21] not in ("O", "o"):
        raise FormatHatasi("obs_not_obs|Bu dosya bir gozlem (O) dosyasi degil.")
    hdr = {"version": ver, "marker": "", "receiver": "", "antenna": "", "approx": [0.0, 0.0, 0.0],
           "delta_hen": [0.0, 0.0, 0.0], "interval": None, "obs_types": [], "leap": None}
    types2 = []
    types3 = {}
    cur_sys = None
    i = 0
    while i < len(L):
        s = L[i]
        lab = s[60:80].strip()
        if lab == "END OF HEADER":
            i += 1
            break
        if lab == "MARKER NAME":
            hdr["marker"] = s[:60].strip()
        elif lab == "REC # / TYPE / VERS":
            hdr["receiver"] = s[20:40].strip()
        elif lab == "ANT # / TYPE":
            hdr["antenna"] = s[20:40].rstrip()
        elif lab == "APPROX POSITION XYZ":
            hdr["approx"] = [_f(s[0:14]), _f(s[14:28]), _f(s[28:42])]
        elif lab == "ANTENNA: DELTA H/E/N":
            hdr["delta_hen"] = [_f(s[0:14]), _f(s[14:28]), _f(s[28:42])]
        elif lab == "INTERVAL":
            hdr["interval"] = _f(s[:10])
        elif lab == "LEAP SECONDS":
            hdr["leap"] = _f(s[:6])
        elif lab == "# / TYPES OF OBSERV":
            types2 += s[6:60].split()
        elif lab == "SYS / # / OBS TYPES":
            if s[0] != " ":
                cur_sys = s[0]
                types3[cur_sys] = []
            types3[cur_sys] += s[7:60].split()
        i += 1
    else:
        raise FormatHatasi("obs_no_header_end|Gozlem dosyasinda 'END OF HEADER' bulunamadi.")

    if ver >= 3:
        types = types3.get("G", [])
        if not types:
            raise FormatHatasi("obs_no_gps_types|RINEX 3 dosyasinda GPS gozlem tipi yok.")
        hdr["obs_types"] = types
        idx = {k: next((types.index(c) for c in v if c in types), None) for k, v in _V3_PREF.items()}
        hdr["used_codes"] = {k: (types[j] if j is not None else None) for k, j in idx.items()}
        recs = _read_body3(L, i, types, idx)
    else:
        types = types2
        hdr["obs_types"] = types
        def pick(*names):
            return next((types.index(n) for n in names if n in types), None)
        idx = {"P1": pick("P1", "C1"), "C1": pick("C1"), "P2": pick("P2", "C2"),
               "L1": pick("L1"), "L2": pick("L2")}
        hdr["used_codes"] = {k: (types[j] if j is not None else None) for k, j in idx.items()}
        recs = _read_body2(L, i, types, idx)

    if idx["L1"] is None or idx["P1"] is None:
        raise FormatHatasi("obs_no_l1|Dosyada L1 fazi veya L1 kodu yok; faz cozumu yapilamaz.")
    if not recs:
        raise FormatHatasi("obs_no_gps|Dosyada GPS gozlemi bulunamadi.")
    arr = np.array(recs, dtype=float)
    out = {"header": hdr, "t": arr[:, 0], "prn": arr[:, 1].astype(int),
           "P1": arr[:, 2], "C1": arr[:, 3], "P2": arr[:, 4], "L1": arr[:, 5], "L2": arr[:, 6]}
    for k in ("P1", "C1", "P2", "L1", "L2"):
        v = out[k]
        v[v == 0.0] = np.nan
    ep = np.unique(out["t"])
    if hdr["interval"] is None or not hdr["interval"] or math.isnan(hdr["interval"]):
        d = np.diff(ep)
        hdr["interval"] = float(np.median(d)) if len(d) else 30.0
    return out


def _val(line, j):
    a = 16 * j
    return _f(line[a:a + 14]) if len(line) > a else math.nan


def _read_body2(L, i, types, idx):
    nt = len(types)
    nl = (nt + 4) // 5
    recs = []
    keys = ("P1", "C1", "P2", "L1", "L2")
    n = len(L)
    while i < n:
        s = L[i]
        if len(s) < 32 or not s[:26].strip():
            i += 1
            continue
        try:
            yy, mo, dd, hh, mi = (int(s[1:3]), int(s[4:6]), int(s[7:9]), int(s[10:12]), int(s[13:15]))
            sec = float(s[15:26])
            flag = int(s[28:29] or 0)
            ns = int(s[29:32])
        except ValueError:
            i += 1
            continue
        if flag > 1:                      # olay kaydi: ns satir atla
            i += 1 + ns
            continue
        sats = []
        row = s[32:68]
        j = i
        while len(sats) < ns:
            for k in range(0, 36, 3):
                c = row[k:k + 3]
                if len(sats) < ns and c.strip():
                    sats.append(c)
            if len(sats) < ns:
                j += 1
                row = L[j][32:68] if j < n else ""
                if j >= n:
                    break
        i = j + 1
        t = to_gps(yy, mo, dd, hh, mi, sec)
        for sv in sats:
            block = []
            for _ in range(nl):
                block.append((L[i] if i < n else "").ljust(80))
                i += 1
            if sv[0] not in ("G", " "):
                continue
            try:
                prn = int(sv[1:3])
            except ValueError:
                continue
            vals = {}
            for k in keys:
                jj = idx[k]
                vals[k] = math.nan if jj is None else _val(block[jj // 5], jj % 5)
            recs.append((t, prn, vals["P1"], vals["C1"], vals["P2"], vals["L1"], vals["L2"]))
    return recs


def _read_body3(L, i, types, idx):
    recs = []
    keys = ("P1", "C1", "P2", "L1", "L2")
    n = len(L)
    t = None
    while i < n:
        s = L[i]
        if s.startswith(">"):
            try:
                y, mo, d, h, mi = int(s[2:6]), int(s[7:9]), int(s[10:12]), int(s[13:15]), int(s[16:18])
                sec = float(s[18:29])
                flag = int(s[31:32] or 0)
                ns = int(s[32:35])
            except ValueError:
                i += 1
                continue
            if flag > 1:
                i += 1 + ns
                t = None
                continue
            t = to_gps(y, mo, d, h, mi, sec)
            i += 1
            continue
        if t is not None and s[:1] == "G" and len(s) > 3:
            try:
                prn = int(s[1:3])
            except ValueError:
                i += 1
                continue
            body = s[3:]
            vals = {}
            for k in keys:
                jj = idx[k]
                vals[k] = math.nan if jj is None else _val(body, jj)
            recs.append((t, prn, vals["P1"], vals["C1"], vals["P2"], vals["L1"], vals["L2"]))
        i += 1
    return recs


# ---------------------------------------------------------------------------
# RINEX navigasyon (GPS yayin efemerisi)
# ---------------------------------------------------------------------------
NAV_KEYS = ["af0", "af1", "af2", "IODE", "Crs", "dn", "M0", "Cuc", "e", "Cus", "sqrtA", "Toe", "Cic",
            "OMEGA0", "Cis", "i0", "Crc", "omega", "OMEGADOT", "IDOT", "codes", "week", "L2P",
            "acc", "health", "TGD", "IODC", "TTM", "fit"]


def read_nav(text):
    L = _lines(text)
    if not L or "RINEX VERSION" not in L[0]:
        raise FormatHatasi("nav_not_rinex|Navigasyon dosyasi RINEX degil.")
    ver = _f(L[0][:9])
    typ = L[0][20:21]
    sysch = L[0][40:41]
    if typ not in ("N", "n"):
        raise FormatHatasi("nav_not_nav|Bu dosya bir navigasyon (N) dosyasi degil.")
    if ver >= 3 and sysch not in ("G", "M", " "):
        raise FormatHatasi("nav_no_gps|Navigasyon dosyasi GPS icermiyor.")
    alpha = beta = None
    i = 0
    while i < len(L):
        s = L[i]
        lab = s[60:80].strip()
        if lab == "ION ALPHA":
            alpha = [_f(s[2 + 12 * k:14 + 12 * k]) for k in range(4)]
        elif lab == "ION BETA":
            beta = [_f(s[2 + 12 * k:14 + 12 * k]) for k in range(4)]
        elif lab == "IONOSPHERIC CORR":
            vals = [_f(s[5 + 12 * k:17 + 12 * k]) for k in range(4)]
            if s[:4] == "GPSA":
                alpha = vals
            elif s[:4] == "GPSB":
                beta = vals
        elif lab == "END OF HEADER":
            i += 1
            break
        i += 1
    eph = []
    n = len(L)
    if ver >= 3:
        while i < n:
            s = L[i]
            if not s.strip():
                i += 1
                continue
            sysc = s[0]
            nlines = {"G": 8, "E": 8, "C": 8, "J": 8, "I": 8, "R": 4, "S": 4}.get(sysc, 8)
            if sysc != "G":
                i += nlines
                continue
            try:
                prn = int(s[1:3])
                y, mo, d, h, mi, sec = int(s[4:8]), int(s[9:11]), int(s[12:14]), int(s[15:17]), int(s[18:20]), float(s[21:23])
            except ValueError:
                i += 1
                continue
            vals = [_f(s[23 + 19 * k:42 + 19 * k]) for k in range(3)]
            for j in range(1, 8):
                ss = L[i + j] if i + j < n else ""
                vals += [_f(ss[4 + 19 * k:23 + 19 * k]) for k in range(4)]
            i += 8
            eph.append(_mk_eph(prn, to_gps(y, mo, d, h, mi, sec), vals))
    else:
        while i < n:
            s = L[i]
            if len(s) < 22 or not s[:2].strip():
                i += 1
                continue
            try:
                prn = int(s[0:2])
                y, mo, d, h, mi, sec = int(s[3:5]), int(s[6:8]), int(s[9:11]), int(s[12:14]), int(s[15:17]), float(s[17:22])
            except ValueError:
                i += 1
                continue
            vals = [_f(s[22 + 19 * k:41 + 19 * k]) for k in range(3)]
            for j in range(1, 8):
                ss = L[i + j] if i + j < n else ""
                vals += [_f(ss[3 + 19 * k:22 + 19 * k]) for k in range(4)]
            i += 8
            eph.append(_mk_eph(prn, to_gps(y, mo, d, h, mi, sec), vals))
    if not eph:
        raise FormatHatasi("nav_no_eph|Navigasyon dosyasinda GPS efemerisi bulunamadi.")
    klob = None
    if alpha and beta and not any(math.isnan(v) for v in alpha + beta):
        klob = {"alpha": alpha, "beta": beta}
    return {"eph": eph, "klobuchar": klob}


def _mk_eph(prn, toc, vals):
    e = {"prn": prn, "toc": toc}
    for k, v in zip(NAV_KEYS, vals):
        e[k] = v
    # Toe: hafta saniyesi -> mutlak GPS saniyesi (toc ile ayni haftaya yakin)
    wk_toc = math.floor(toc / 604800.0)
    toe = wk_toc * 604800.0 + e["Toe"]
    if toe - toc > 302400:
        toe -= 604800
    elif toe - toc < -302400:
        toe += 604800
    e["toe_abs"] = toe
    return e


# ---------------------------------------------------------------------------
# RINEX meteoroloji
# ---------------------------------------------------------------------------
def read_met(text):
    L = _lines(text)
    if not L or "RINEX VERSION" not in L[0]:
        raise FormatHatasi("met_not_rinex|Meteoroloji dosyasi RINEX degil.")
    ver = _f(L[0][:9])
    types = []
    height = math.nan
    i = 0
    while i < len(L):
        s = L[i]
        lab = s[60:80].strip()
        if lab == "# / TYPES OF OBSERV":
            types += s[6:60].split()
        elif lab == "SENSOR POS XYZ/H":
            height = _f(s[42:56])
        elif lab == "END OF HEADER":
            i += 1
            break
        i += 1
    if not types:
        raise FormatHatasi("met_no_types|Meteoroloji dosyasinda gozlem tipi yok.")
    rows = []
    while i < len(L):
        s = L[i]
        tok = s.split()
        if len(tok) >= 6 + len(types):
            try:
                y, mo, d, h, mi = (int(x) for x in tok[:5])
                sec = float(tok[5])
                vals = [float(x) for x in tok[6:6 + len(types)]]
                rows.append([to_gps(y, mo, d, h, mi, sec)] + vals)
            except ValueError:
                pass
        i += 1
    if not rows:
        raise FormatHatasi("met_no_data|Meteoroloji dosyasinda veri bulunamadi.")
    a = np.array(rows)
    out = {"t": a[:, 0], "height": height}
    for j, ty in enumerate(types):
        out[ty] = a[:, 1 + j]
    for need in ("PR", "TD", "HR"):
        if need not in out:
            raise FormatHatasi(f"met_missing_type|Meteoroloji dosyasinda {need} yok (gerekli: PR, TD, HR).")
    return out


# ---------------------------------------------------------------------------
# SP3 kesin yorunge
# ---------------------------------------------------------------------------
def read_sp3(texts):
    """Bir veya birden cok SP3 metnini okur, birlestirir.

    Doner: dict(t=(n,), prns=[...], xyz=(n, nsat, 3) [m], clk=(n, nsat) [s])
    """
    if isinstance(texts, str):
        texts = [texts]
    data = {}
    for text in texts:
        L = _lines(text)
        if not L or L[0][:1] != "#":
            raise FormatHatasi("sp3_bad|SP3 dosyasi '#' ile baslamiyor.")
        t = None
        for s in L:
            if s.startswith("*"):
                tok = s[1:].split()
                t = to_gps(int(tok[0]), int(tok[1]), int(tok[2]), int(tok[3]), int(tok[4]), float(tok[5]))
            elif s.startswith("P") and t is not None:
                sysc = s[1]
                if sysc not in ("G", " "):
                    continue
                try:
                    prn = int(s[2:4])
                except ValueError:
                    continue
                x, y, z = _f(s[4:18]), _f(s[18:32]), _f(s[32:46])
                clk = _f(s[46:60])
                if x == 0 and y == 0 and z == 0:
                    continue
                c = clk * 1e-6 if (not math.isnan(clk) and abs(clk) < 999999) else math.nan
                data.setdefault(t, {})[prn] = (x * 1e3, y * 1e3, z * 1e3, c)
    if not data:
        raise FormatHatasi("sp3_no_gps|SP3 dosyasinda GPS konumu bulunamadi.")
    ts = np.array(sorted(data))
    prns = sorted({p for d in data.values() for p in d})
    xyz = np.full((len(ts), len(prns), 3), np.nan)
    clk = np.full((len(ts), len(prns)), np.nan)
    pi = {p: k for k, p in enumerate(prns)}
    for a, t in enumerate(ts):
        for p, v in data[t].items():
            xyz[a, pi[p]] = v[:3]
            clk[a, pi[p]] = v[3]
    return {"t": ts, "prns": prns, "xyz": xyz, "clk": clk,
            "interval": float(np.median(np.diff(ts))) if len(ts) > 1 else 900.0}


# ---------------------------------------------------------------------------
# ANTEX
# ---------------------------------------------------------------------------
def _atx_date(s):
    tok = s[:43].split()
    if len(tok) < 5:
        return None
    return datetime(int(tok[0]), int(tok[1]), int(tok[2]), int(tok[3]), int(tok[4]))


def read_antex(text, when, rec_type=None):
    """ANTEX'ten verilen tarihte gecerli uydu antenlerini ve (verilmisse) alici antenini okur.

    Doner: {'sat': {prn: ant}, 'rec': ant|None, 'rec_name': str}
    ant = {'block': str, 'zen': (z0, dz), 'G01': {'pco': (3,) mm, 'pcv': (n,) mm}, 'G02': {...}}
    """
    sats, recs = {}, {}
    ant = None
    fr = None
    for s in _lines(text):
        lab = s[60:80].strip()
        if lab == "START OF ANTENNA":
            ant = {"freq": {}, "from": None, "until": None}
        elif ant is None:
            continue
        elif lab == "TYPE / SERIAL NO":
            ant["type"], ant["serial"] = s[0:20], s[20:40].strip()
        elif lab == "ZEN1 / ZEN2 / DZEN":
            tok = s[:20].split()
            ant["zen"] = (float(tok[0]), float(tok[2]))
        elif lab == "VALID FROM":
            ant["from"] = _atx_date(s)
        elif lab == "VALID UNTIL":
            ant["until"] = _atx_date(s)
        elif lab == "START OF FREQUENCY":
            fr = s[3:6]
            ant["freq"][fr] = {}
        elif lab == "NORTH / EAST / UP" and fr:
            ant["freq"][fr]["pco"] = np.array([_f(s[0:10]), _f(s[10:20]), _f(s[20:30])])
        elif fr and s[3:8] == "NOAZI":
            ant["freq"][fr]["pcv"] = np.array([float(x) for x in s[8:].split()])
        elif lab == "END OF FREQUENCY":
            fr = None
        elif lab == "END OF ANTENNA":
            ok = (ant["from"] is None or ant["from"] <= when) and (ant["until"] is None or when <= ant["until"])
            if ok and "G01" in ant["freq"] and "zen" in ant:
                rec = {"block": ant["type"].strip(), "zen": ant["zen"], **ant["freq"]}
                if "G02" not in rec:
                    rec["G02"] = rec["G01"]
                serial = ant["serial"]
                if ant["type"].startswith("BLOCK") or ant["type"].startswith("GPS-"):
                    if serial[:1] == "G" and serial[1:3].isdigit():
                        sats[int(serial[1:3])] = rec
                else:
                    recs[ant["type"].rstrip()] = rec
            ant = None
    out = {"sat": sats, "rec": None, "rec_name": None}
    if rec_type:
        name = rec_type.rstrip()
        cand = [name, name[:16].rstrip().ljust(16) + "NONE", name[:16].rstrip()]
        for c in cand:
            for k in recs:
                if k.rstrip() == c.rstrip() or k.replace(" ", "") == c.replace(" ", ""):
                    out["rec"], out["rec_name"] = recs[k], k
                    break
            if out["rec"] is not None:
                break
    return out


# ---------------------------------------------------------------------------
# CODE DCB (P1-P2), ns
# ---------------------------------------------------------------------------
def read_dcb(text):
    d = {}
    for s in _lines(text):
        if len(s) > 4 and s[0] == "G" and s[1:3].strip().isdigit():
            tok = s[3:].split()
            if tok:
                try:
                    d[int(s[1:3])] = float(tok[0])
                except ValueError:
                    pass
    if not d:
        raise FormatHatasi("dcb_bad|DCB dosyasinda GPS uydu degeri bulunamadi (CODE P1-P2 formati bekleniyor).")
    return d
