"""Tum hesap zinciri: dosyalar -> sonuc sozlugu (JSON'a donusturulebilir)."""
import gzip
import math
import time
from datetime import datetime

import numpy as np

from . import readers, preprocess as pre, orbit, corrections as cr, estimate as est, geo
from .const import (C, F1, F2, LAM1, LAM2, A3, B3, GAMMA, LAM_NL, TECU_L1, OMEGA_E, MS)
from .timeutil import fmt, from_gps, day_start


class Hata(Exception):
    pass


def _r(x, nd=4):
    if x is None:
        return None
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return x
    if math.isnan(xf) or math.isinf(xf):
        return None
    return round(xf, nd)


def _lst(a, nd=4):
    return [_r(v, nd) for v in np.asarray(a).ravel()]


def _dec(n, maxn):
    """Grafik icin seyreltme indeksleri."""
    if n <= maxn:
        return np.arange(n)
    return np.linspace(0, n - 1, maxn).astype(int)


def run(files, params=None, progress=None, default_atx=None):
    params = dict(params or {})
    T0 = time.time()
    warn = []

    def step(msg, frac):
        if progress:
            progress(msg, frac)

    cutoff = float(params.get("cutoff", 10.0))
    sig_phase = float(params.get("sig_phase", 0.01))
    sig_code = float(params.get("sig_code", 1.0))
    zwd_block = float(params.get("zwd_block_h", 2.0)) * 3600.0
    spp_step = float(params.get("spp_step", 300.0))

    # ------------------------------------------------------------------ 1. okuma
    step("read", 0.02)
    obs = readers.read_obs(files["obs"])
    H = obs["header"]
    sp3 = readers.read_sp3(files["sp3"])
    nav = readers.read_nav(files["nav"]) if files.get("nav") else None
    met = readers.read_met(files["met"]) if files.get("met") else None
    dcb = readers.read_dcb(files["dcb"]) if files.get("dcb") else None
    t_first = float(obs["t"].min())
    day0 = day_start(t_first)
    when = from_gps(t_first)
    atx_text = files.get("atx")
    atx_src = "user"
    if not atx_text and default_atx:
        with open(default_atx, "rb") as fh:
            raw = fh.read()
        atx_text = (gzip.decompress(raw) if default_atx.endswith(".gz") else raw).decode("latin-1")
        atx_src = "bundled-igs08"
        if when > datetime(2017, 1, 28):
            warn.append("atx_old")
    atx = readers.read_antex(atx_text, when, H["antenna"]) if atx_text else {"sat": {}, "rec": None, "rec_name": None}
    if atx_text and atx["rec"] is None:
        warn.append("rec_ant_missing")
    if not atx_text:
        warn.append("no_atx")
    if nav is None:
        warn.append("no_nav")
    elif nav["klobuchar"] is None:
        warn.append("no_klobuchar")
    if met is None:
        warn.append("no_met")
    if dcb is None:
        warn.append("no_dcb")

    # ------------------------------------------------------------------ 2. on isleme
    step("preprocess", 0.10)
    P = pre.preprocess(obs)
    t, prn = P["t"], P["prn"]
    n = len(t)
    if n == 0:
        raise Hata("no_dual_freq|Dort olcusu (P1, P2, L1, L2) tam olan kayit yok; cift frekans gerekli.")

    # ------------------------------------------------------------------ 3. yaklasik konum
    X0 = np.array(H["approx"], float)
    x0_src = "header"
    if not np.all(np.isfinite(X0)) or np.linalg.norm(X0) < 6.0e6:
        X0 = _spp_first(P, sp3)
        x0_src = "spp"
        if X0 is None:
            raise Hata("no_approx|Yaklasik konum bulunamadi (RINEX basliginda yok, kod cozumu de basarisiz).")
    lat0, lon0, h0 = geo.xyz2llh(*X0)
    Rn = geo.rot_neu(lat0, lon0)

    # ------------------------------------------------------------------ 4. uydu konum ve saati
    step("orbit", 0.20)
    Xs = np.full((n, 3), np.nan)
    Vs = np.full((n, 3), np.nan)
    dts = np.full(n, np.nan)
    for p_ in np.unique(prn):
        m = prn == p_
        t_sat = t[m] - P["P1"][m] / C
        c0 = orbit.sp3_clock(sp3, int(p_), t_sat)
        ttx = t_sat - c0
        pos, vel = orbit.sp3_interp(sp3, int(p_), ttx)
        Xs[m], Vs[m], dts[m] = pos, vel, c0
    rel = -2 * np.einsum("ij,ij->i", Xs, Vs) / C ** 2
    dts_tot = dts + rel
    valid = np.isfinite(Xs).all(axis=1) & np.isfinite(dts)
    if valid.sum() == 0:
        raise Hata("sp3_no_overlap|SP3 dosyasi gozlem gunuyle ortusmuyor veya uydu saati yok.")
    tau = np.linalg.norm(Xs - X0, axis=1) / C
    for _ in range(2):
        th = OMEGA_E * tau
        Xr = np.column_stack([np.cos(th) * Xs[:, 0] + np.sin(th) * Xs[:, 1],
                              -np.sin(th) * Xs[:, 0] + np.cos(th) * Xs[:, 1], Xs[:, 2]])
        tau = np.linalg.norm(Xr - X0, axis=1) / C
    sagnac = np.linalg.norm(Xr - X0, axis=1) - np.linalg.norm(Xs - X0, axis=1)
    Xs = Xr
    el = np.full(n, np.nan)
    az = np.full(n, np.nan)
    el[valid], az[valid] = geo.elaz(X0, Xs[valid])
    rho0 = np.linalg.norm(Xs - X0, axis=1)
    u = (Xs - X0) / rho0[:, None]

    # ------------------------------------------------------------------ 5. troposfer
    step("atmosphere", 0.30)
    if met is not None:
        Pm = np.interp(t, met["t"], met["PR"])
        Tm_ = np.interp(t, met["t"], met["TD"])
        RH = np.interp(t, met["t"], met["HR"])
        hm = met["height"] if not math.isnan(met.get("height", math.nan)) else h0
        if abs(hm - h0) > 20:                     # sensor yuksekligi farkliysa basinci indirge
            Pm = Pm * ((1 - 2.2557e-5 * h0) / (1 - 2.2557e-5 * hm)) ** 5.2568
    else:
        Ps_, Ts_, RHs_ = cr.standard_met(h0)
        Pm, Tm_, RH = np.full(n, Ps_), np.full(n, Ts_), np.full(n, RHs_)
    ZHD, ZWD = cr.saastamoinen(Pm, Tm_, RH, lat0, h0)
    mh = np.full(n, np.nan)
    mwet = np.full(n, np.nan)
    pos_el = valid & (el > 0)
    mh[pos_el], mwet[pos_el] = cr.herring(el[pos_el], lat0, h0, Tm_[pos_el])
    Ttrop = ZHD * mh + ZWD * mwet

    # ------------------------------------------------------------------ 6. kombinasyonlar ve 4B
    step("corrections", 0.40)
    L1m = LAM1 * P["L1d"]
    L2m = LAM2 * P["L2d"]
    L3 = A3 * L1m + B3 * L2m
    P3 = A3 * P["P1"] + B3 * P["P2"]
    sun, moon = geo.sun_moon_ecef(t)
    ex = ey = ez = None
    d_sat = np.zeros(n)
    nadir = np.full(n, np.nan)
    if valid.any():
        ex, ey, ez = cr.sat_body_axes(np.where(valid[:, None], Xs, 1.0), sun)
        coef3 = {"G01": A3, "G02": B3}
        ds_, nad_, miss = cr.satellite_antenna(u[valid], prn[valid], ex[valid], ey[valid], ez[valid], atx["sat"], coef3)
        d_sat[valid] = ds_
        nadir[valid] = nad_
        if atx_text and miss:
            warn.append("sat_ant_missing:" + ",".join(f"G{p:02d}" for p in sorted(miss)))
    d_rec3 = np.zeros(n)
    d_rec1 = np.zeros(n)
    d_rec2 = np.zeros(n)
    elv = np.where(valid, el, 90.0)
    uv = np.where(valid[:, None], u, 0.0)
    d_rec3 = cr.receiver_antenna(uv, elv, atx["rec"], H["delta_hen"], Rn, {"G01": A3, "G02": B3})
    d_rec1 = cr.receiver_antenna(uv, elv, atx["rec"], H["delta_hen"], Rn, {"G01": 1.0})
    d_rec2 = cr.receiver_antenna(uv, elv, atx["rec"], H["delta_hen"], Rn, {"G02": 1.0})
    dr_tide = cr.solid_tide(X0, sun, moon)
    d_tide = -np.einsum("ij,ij->i", uv, dr_tide)
    tide_up = dr_tide @ Rn[2]
    d_shap = np.where(valid, cr.shapiro(np.where(valid[:, None], Xs, X0 + 2e7), X0, np.where(valid, rho0, 2e7)), 0.0)
    wu = np.zeros(n)
    if ex is not None:
        wraw = np.where(valid, cr.windup_raw(uv + (~valid)[:, None] * np.array([0, 0, 1.0]), ex, ey, Rn), 0.0)
        wu = cr.unwrap_by_arc(wraw, np.where(valid, P["arc"], -1), t)
    d_wu3 = LAM_NL * wu
    common = d_sat + d_tide + d_shap
    corr_P3 = d_rec3 + common
    corr_L3 = corr_P3 + d_wu3

    # ------------------------------------------------------------------ 7. secim ve epok
    sp3_int = sp3["interval"]
    on_sp3 = np.abs(((t - sp3["t"][0] + 1e-3) % sp3_int) - 1e-3) < 1e-3
    base = valid & P["usable"] & (el >= cutoff) & np.isfinite(Ttrop)
    sel = base & on_sp3
    arcs_sel, cnt = np.unique(P["arc"][sel], return_counts=True)
    good_arcs = set(arcs_sel[cnt >= 3].tolist())
    sel &= np.isin(P["arc"], list(good_arcs))
    if sel.sum() < 20:
        raise Hata("too_few_obs|Cozum icin yeterli olcu yok (SP3 epoklarinda, yukseklik sinirinin uzerinde).")
    ep_sel = np.unique(t[sel])
    req = params.get("epoch")
    if req:
        hh, mm, ss = (list(map(float, str(req).split(":"))) + [0, 0])[:3]
        t_req = day0 + hh * 3600 + mm * 60 + ss
    else:
        t_req = ep_sel[len(ep_sel) // 2]
    t_ep = float(ep_sel[np.argmin(np.abs(ep_sel - t_req))])
    if abs(t_ep - t_req) > 0.5:
        warn.append("epoch_snapped")

    dsel = {"t": t[sel], "arc": P["arc"][sel], "Xs": Xs[sel], "el": el[sel], "cdts": C * dts_tot[sel],
            "Tmodel": Ttrop[sel], "mw": mwet[sel]}

    # ------------------------------------------------------------------ 8. L3 cozumu
    step("l3", 0.50)
    sol = est.phase_solution(dsel, L3[sel], P3[sel], corr_L3[sel], corr_P3[sel], X0, t_ep, Rn,
                             sig_phase, sig_code, zwd_block)
    sol_no4b = est.phase_solution(dsel, L3[sel], P3[sel], np.zeros(sel.sum()), np.zeros(sel.sum()), X0, t_ep, Rn,
                                  sig_phase, sig_code, zwd_block)
    Xst = sol["X"]
    latS, lonS, hS = geo.xyz2llh(*Xst)

    ref = params.get("ref_xyz")
    ref_src = "user"
    if ref and len(ref) == 3 and all(r is not None for r in ref):
        Xref = np.array(ref, float)
    else:
        Xref = Xst.copy()
        ref_src = "static"

    def neu(X):
        if X is None:
            return None
        d = Rn @ (np.asarray(X) - Xref)
        return {"N": _r(d[0]), "E": _r(d[1]), "U": _r(d[2]), "H": _r(math.hypot(d[0], d[1])),
                "D3": _r(float(np.linalg.norm(d)))}

    j_ep = int(np.searchsorted(sol["epochs"], t_ep - 1e-3))
    ztd_ep = float(ZHD[sel][np.argmin(np.abs(t[sel] - t_ep))] + ZWD[sel][np.argmin(np.abs(t[sel] - t_ep))]
                   + (sol["zwd_kin"] if sol["has_kin"] else 0.0))
    ztd_model_ep = float(ZHD[sel][np.argmin(np.abs(t[sel] - t_ep))] + ZWD[sel][np.argmin(np.abs(t[sel] - t_ep))])

    # ------------------------------------------------------------------ 9. epok tablolari
    step("tables", 0.60)
    me = (np.abs(t - t_ep) < 1e-3)
    idx_ep = np.where(me)[0]
    idx_ep = idx_ep[np.argsort(-np.nan_to_num(el[idx_ep], nan=-99))]
    used_ep = set(np.where(me & sel)[0].tolist())
    sat_rows, atm_rows, cor_rows, a5_rows = [], [], [], []
    for i in idx_ep:
        ok = valid[i]
        sv = f"G{prn[i]:02d}"
        sat_rows.append({"sv": sv, "X": _r(Xs[i, 0], 3), "Y": _r(Xs[i, 1], 3), "Z": _r(Xs[i, 2], 3),
                         "rho": _r(rho0[i], 3), "tau_ms": _r(tau[i] * 1e3, 3), "dts_us": _r(dts[i] * 1e6, 3),
                         "rel_ns": _r(rel[i] * 1e9, 3), "sagnac": _r(sagnac[i], 3), "el": _r(el[i], 2),
                         "az": _r(az[i], 2), "used": i in used_ep, "valid": bool(ok)})
        atm_rows.append({"sv": sv, "el": _r(el[i], 2), "mh": _r(mh[i], 3), "mw": _r(mwet[i], 3),
                         "T": _r(Ttrop[i], 3), "I1_code": _r((P["P2"][i] - P["P1"][i]) / (GAMMA - 1), 3),
                         "used": i in used_ep})
        cor_rows.append({"sv": sv, "el": _r(el[i], 1), "nadir": _r(nadir[i], 1),
                         "rec": _r(d_rec3[i] * 1e3, 1), "sat": _r(d_sat[i] * 1e3, 1),
                         "tide": _r(d_tide[i] * 1e3, 1), "shapiro": _r(d_shap[i] * 1e3, 1),
                         "windup": _r(d_wu3[i] * 1e3, 1), "total": _r(corr_L3[i] * 1e3, 1),
                         "used": i in used_ep})
        a5_rows.append({"sv": sv, "P1": _r(P["P1"][i], 3), "P2": _r(P["P2"][i], 3), "P3": _r(P3[i], 3),
                        "L1": _r(L1m[i], 3), "L2": _r(L2m[i], 3), "L3": _r(L3[i], 3), "used": i in used_ep})
    k_ep = P["K"][idx_ep[0]] if len(idx_ep) else 0
    check = None
    if len(used_ep):
        ii = np.array(sorted(used_ep))
        k1 = P["P1"][ii] - rho0[ii] + C * dts_tot[ii]
        k3 = P3[ii] - rho0[ii] + C * dts_tot[ii]
        k3t = k3 - Ttrop[ii]
        check = {"spread_P1": _r(np.ptp(k1), 3), "spread_P3": _r(np.ptp(k3), 3), "spread_P3T": _r(np.ptp(k3t), 3),
                 "clock_rough_ms": _r(np.mean(k1) / MS, 5)}

    # ------------------------------------------------------------------ 10. iyonosfer ve TEC
    step("iono", 0.68)
    iono = _ionosphere(P, el, az, valid, cutoff, lat0, lon0, t, t_ep, dcb, nav, sel, idx_ep, used_ep, warn)

    # ------------------------------------------------------------------ 11. tek frekans cozumleri
    step("single", 0.75)
    single = []
    I1_mod = None
    klob = nav["klobuchar"] if nav else None
    if klob is not None:
        I1_mod = np.zeros(n)
        I1_mod[valid] = cr.klobuchar(t[valid], lat0, lon0, el[valid], az[valid], klob["alpha"], klob["beta"])
    dcbm = np.zeros(n)
    if dcb:
        dcbm = np.array([dcb.get(int(p_), 0.0) for p_ in prn]) * 1e-9 * C
    variants = [("L1", "none", None)]
    if I1_mod is not None:
        variants += [("L1", "klobuchar", I1_mod), ("L2", "klobuchar", I1_mod)]
    for fr, model, Iv in variants:
        g = 1.0 if fr == "L1" else GAMMA
        Isel = (Iv[sel] * g) if Iv is not None else np.zeros(sel.sum())
        phase = (L1m if fr == "L1" else L2m)[sel]
        code = (P["P1"] if fr == "L1" else P["P2"])[sel]
        drec = (d_rec1 if fr == "L1" else d_rec2)[sel]
        kd = (B3 * dcbm if fr == "L1" else -A3 * dcbm)[sel]
        wl = (LAM1 if fr == "L1" else LAM2) * wu[sel]
        try:
            s1 = est.phase_solution(dsel, phase, code, common[sel] + drec + wl - Isel,
                                    common[sel] + drec + Isel + kd, X0, t_ep, Rn, sig_phase, sig_code, zwd_block)
        except np.linalg.LinAlgError:
            continue
        jj = int(np.searchsorted(s1["epochs"], t_ep - 1e-3))
        single.append({"freq": fr, "model": model, "X": _lst(s1["X"]), "Xk": _lst(s1["Xk"]) if s1["has_kin"] else None,
                       "neu_static": neu(s1["X"]), "neu_kin": neu(s1["Xk"]), "rms_phase_mm": _r(s1["rms_phase"] * 1e3, 1),
                       "clock_ns": _r(s1["clock"][jj] / C * 1e9, 1), "ztd": _r(ztd_model_ep + (s1["zwd_kin"] if s1["has_kin"] else 0), 4)})
    l3row = {"freq": "L3", "model": "iono-free", "X": _lst(Xst), "Xk": _lst(sol["Xk"]) if sol["has_kin"] else None,
             "neu_static": neu(Xst), "neu_kin": neu(sol["Xk"]), "rms_phase_mm": _r(sol["rms_phase"] * 1e3, 1),
             "clock_ns": _r(sol["clock"][j_ep] / C * 1e9, 1), "ztd": _r(ztd_ep, 4)}
    single.append(l3row)

    # ------------------------------------------------------------------ 12. cok yolluluk
    step("multipath", 0.82)
    mp = _multipath(P, el, az, valid, idx_ep, cutoff)

    # ------------------------------------------------------------------ 13. kod SPP
    step("spp", 0.88)
    spp = _spp_series(P, sp3, nav, X0, lat0, h0, Rn, Xref, t_ep, spp_step, cutoff, met, Pm, Tm_, RH)

    # ------------------------------------------------------------------ 14. grafik verileri
    step("plots", 0.94)
    plots = _plots(P, t, prn, el, az, valid, sel, sol, idx_ep, t_ep, Ttrop, wu, tide_up, iono, mp, spp, day0)

    # ------------------------------------------------------------------ 15. ozet
    solutions = [
        {"name": "l3_static", "X": _lst(Xst), "sX_mm": _lst(sol["sX"] * 1e3, 1), "neu": neu(Xst)},
        {"name": "l3_kin", "X": _lst(sol["Xk"]) if sol["has_kin"] else None,
         "sX_mm": _lst(sol["sXk"] * 1e3, 1) if sol["has_kin"] else None, "neu": neu(sol["Xk"])},
        {"name": "l3_static_no4b", "X": _lst(sol_no4b["X"]), "sX_mm": _lst(sol_no4b["sX"] * 1e3, 1), "neu": neu(sol_no4b["X"])},
    ]
    for s1 in single[:-1]:
        solutions.append({"name": f"{s1['freq'].lower()}_{s1['model']}_kin", "X": s1["Xk"], "sX_mm": None,
                          "neu": s1["neu_kin"]})
    if spp.get("epoch"):
        solutions.append({"name": "spp_epoch", "X": spp["epoch"]["X"], "sX_mm": None, "neu": neu(spp["epoch"]["X"])})
    if spp.get("mean"):
        solutions.append({"name": "spp_mean", "X": spp["mean"], "sX_mm": None, "neu": neu(spp["mean"])})
    solutions.append({"name": "header", "X": _lst(H["approx"]), "sX_mm": None,
                      "neu": neu(H["approx"]) if np.linalg.norm(H["approx"]) > 6e6 else None})

    res = {
        "meta": {"version": "1.0.0", "runtime_s": _r(time.time() - T0, 1), "warnings": warn,
                 "atx_source": atx_src if atx_text else None, "rec_ant_found": atx["rec_name"],
                 "x0_source": x0_src, "ref_source": ref_src, "ref_xyz": _lst(Xref),
                 "params": {"cutoff": cutoff, "sig_phase": sig_phase, "sig_code": sig_code,
                            "zwd_block_h": zwd_block / 3600, "spp_step": spp_step}},
        "station": {"marker": H["marker"], "receiver": H["receiver"], "antenna": H["antenna"].strip(),
                    "approx": _lst(H["approx"]), "delta_hen": _lst(H["delta_hen"]),
                    "obs_types": H["obs_types"], "used_codes": H["used_codes"], "rinex_version": H["version"],
                    "interval": P["interval"], "first": fmt(t_first, True), "last": fmt(float(obs["t"].max()), True),
                    "n_epochs": int(len(np.unique(obs["t"]))), "n_sats": int(len(np.unique(obs["prn"]))),
                    "lat": geo.dms(latS), "lon": geo.dms(lonS), "h": _r(hS, 4),
                    "lat_deg": _r(math.degrees(lat0), 6), "lon_deg": _r(math.degrees(lon0), 6), "h0": _r(h0, 3)},
        "epoch": {"requested": fmt(t_req), "used": fmt(t_ep), "used_full": fmt(t_ep, True), "n_sat": int(len(idx_ep)),
                  "n_used": len(used_ep), "clock_jumps_ms": int(k_ep)},
        "data": {"sp3_interval": sp3_int, "sp3_first": fmt(sp3["t"][0], True), "sp3_last": fmt(sp3["t"][-1], True),
                 "nav_n": len(nav["eph"]) if nav else 0, "klobuchar": klob,
                 "met": {"P": _r(np.interp(t_ep, t, Pm), 1), "T": _r(np.interp(t_ep, t, Tm_), 1),
                         "RH": _r(np.interp(t_ep, t, RH), 1), "source": "met" if met is not None else "standard"},
                 "zhd": _r(np.interp(t_ep, t, ZHD), 4), "zwd": _r(np.interp(t_ep, t, ZWD), 4)},
        "pre": {**P["stats"], "jumps": [{"t": fmt(j["t"]), "median_m": _r(j["median_m"], 3), "nsat": j["nsat"],
                                         "k_ms": j["k_ms"]} for j in P["jumps"]],
                "arcs_epoch": _arcs_epoch(P, idx_ep, prn)},
        "orbit": {"sats": sat_rows, "check": check},
        "atm": {"rows": atm_rows, "tide_up_cm": _r(np.interp(t_ep, t, tide_up) * 100, 1),
                "tide_range_cm": [_r(np.nanmin(tide_up) * 100, 1), _r(np.nanmax(tide_up) * 100, 1)],
                "rec_pco_l3_mm": _lst((A3 * atx["rec"]["G01"]["pco"] + B3 * atx["rec"]["G02"]["pco"]) if atx["rec"] else [0, 0, 0], 2)},
        "corr": {"rows": cor_rows},
        "a5": a5_rows,
        "l3": {"X": _lst(Xst), "sX_mm": _lst(sol["sX"] * 1e3, 1), "Xk": _lst(sol["Xk"]) if sol["has_kin"] else None,
               "sXk_mm": _lst(sol["sXk"] * 1e3, 1) if sol["has_kin"] else None,
               "kin_minus_static_neu_mm": _lst(Rn @ (sol["Xk"] - Xst) * 1e3, 1) if sol["has_kin"] else None,
               "clock_m": _r(sol["clock_kin"] if sol["has_kin"] else sol["clock"][j_ep], 3),
               "clock_ns": _r((sol["clock_kin"] if sol["has_kin"] else sol["clock"][j_ep]) / C * 1e9, 1),
               "sclock_m": _r(sol.get("sclock_kin", 0), 3), "ztd": _r(ztd_ep, 4), "ztd_model": _r(ztd_model_ep, 4),
               "dop": {k: _r(v, 3) for k, v in sol.get("dop", {}).items()}, "n_sat_kin": sol.get("n_sat_kin"),
               "s0": _r(sol["s0"], 3), "dof": sol["dof"], "n_obs": sol["n_obs"], "n_par": sol["n_par"],
               "n_epochs": sol["n_epochs"], "n_arcs": sol["n_arcs"], "n_zwd": sol["n_zwd"],
               "rms_phase_mm": _r(sol["rms_phase"] * 1e3, 1), "rms_code_m": _r(sol["rms_code"], 2),
               "iter_mm": _lst(np.array(sol["iter"]) * 1e3, 1),
               "no4b_static_minus_full_neu_mm": _lst(Rn @ (sol_no4b["X"] - Xst) * 1e3, 1),
               "no4b_rms_phase_mm": _r(sol_no4b["rms_phase"] * 1e3, 1)},
        "iono": iono["summary"],
        "single": single,
        "mp": mp["summary"],
        "spp": {k: v for k, v in spp.items() if k != "series"},
        "solutions": solutions,
        "plots": plots,
    }
    step("done", 1.0)
    return res


# ======================================================================
def _spp_first(P, sp3):
    """Basliksiz dosyalar icin ilk epoklardan kod konumu."""
    t, prn = P["t"], P["prn"]
    for tt in np.unique(t)[:20]:
        m = t == tt
        if m.sum() < 5:
            continue
        P3 = A3 * P["P1"][m] + B3 * P["P2"][m]
        pr = prn[m]
        Xs = np.zeros((m.sum(), 3))
        cd = np.zeros(m.sum())
        ok = np.ones(m.sum(), bool)
        for k, p_ in enumerate(pr):
            tx = tt - P3[k] / C
            c0 = orbit.sp3_clock(sp3, int(p_), [tx])[0]
            pos, _ = orbit.sp3_interp(sp3, int(p_), [tx - c0])
            Xs[k], cd[k] = pos[0], c0 * C
            ok[k] = np.isfinite(pos[0]).all() and np.isfinite(c0)
        if ok.sum() < 5:
            continue
        r = est.spp_epoch(P3[ok], Xs[ok], cd[ok], 0.0)
        if r and np.linalg.norm(r["X"]) > 6e6:
            return r["X"]
    return None


def _arcs_epoch(P, idx_ep, prn):
    rows = []
    for i in idx_ep:
        a = P["arc"][i]
        info = next((x for x in P["arcs"] if x["arc"] == a), None)
        if info is None:
            rows.append({"sv": f"G{prn[i]:02d}", "start": None, "end": None, "n": 0, "min": 0,
                         "mw_std": None, "usable": False, "outlier": bool(P["outlier"][i])})
            continue
        rows.append({"sv": f"G{prn[i]:02d}", "start": fmt(info["start"]), "end": fmt(info["end"]), "n": info["n"],
                     "min": _r((info["end"] - info["start"]) / 60, 1), "mw_std": _r(info["mw_std"], 3),
                     "usable": info["usable"] and not bool(P["outlier"][i]), "outlier": bool(P["outlier"][i])})
    return rows


def _ionosphere(P, el, az, valid, cutoff, lat0, lon0, t, t_ep, dcb, nav, sel, idx_ep, used_ep, warn):
    """Fazin koda seviyelenmesi, ince tabaka, alici DCB'si ve VTEC."""
    R, Hs = 6371000.0, 450000.0
    k4 = GAMMA - 1
    m = valid & P["usable"] & (el >= cutoff)
    P4 = P["P2"] - P["P1"]
    GF = P["GF"]
    w = np.sin(np.radians(np.where(m, el, 10))) ** 2
    lev = np.full(len(t), np.nan)
    arcs = np.unique(P["arc"][m])
    resid = []
    for a in arcs:
        mm = m & (P["arc"] == a)
        c0 = np.sum((P4[mm] - GF[mm]) * w[mm]) / np.sum(w[mm])
        lev[mm] = GF[mm] + c0
        resid.append(P4[mm] - lev[mm])
    code_rms = float(np.sqrt(np.mean(np.concatenate(resid) ** 2))) if resid else None
    ns_tecu = C * 1e-9 / k4 / TECU_L1
    E = np.radians(np.where(m, el, 45))
    A = np.radians(np.where(m, az, 0))
    sz = R / (R + Hs) * np.cos(E)
    M = 1 / np.sqrt(1 - sz ** 2)
    psi = np.pi / 2 - E - np.arcsin(sz)
    fi = np.arcsin(np.sin(lat0) * np.cos(psi) + np.cos(lat0) * np.sin(psi) * np.cos(A))
    la = lon0 + np.arcsin(np.sin(psi) * np.sin(A) / np.cos(fi))
    dcb_src = "file"
    if not dcb:
        dcb = _estimate_sat_dcb(lev, m, P["prn"], t, M, fi, la, lat0, lon0, w, ns_tecu, k4)
        dcb_src = "estimated" if dcb else "none"
    dcbs = np.array([(dcb or {}).get(int(p_), 0.0) for p_ in P["prn"]])
    stec_raw = (lev + C * 1e-9 * dcbs) / k4 / TECU_L1
    win = m & (np.abs(t - t_ep) <= 7200) & np.isfinite(stec_raw)
    out = {"summary": None}
    if win.sum() < 10:
        warn.append("iono_few")
        out["summary"] = {"available": False}
        out["series"] = None
        return out
    dt = (t[win] - t_ep) / 3600
    Am = np.column_stack([M[win], M[win] * dt, M[win] * np.degrees(fi[win] - lat0), M[win] * np.degrees(la[win] - lon0),
                          -ns_tecu * np.ones(win.sum())])
    wp = w[win]
    xh = np.linalg.solve(Am.T @ (Am * wp[:, None]), Am.T @ (wp * stec_raw[win]))
    v = Am @ xh - stec_raw[win]
    dcb_r = float(xh[4])
    stec = stec_raw + ns_tecu * dcb_r
    vtec = stec / M
    I1 = stec * TECU_L1
    rows = []
    klob = nav["klobuchar"] if nav else None
    for i in idx_ep:
        if not m[i]:
            continue
        row = {"sv": f"G{P['prn'][i]:02d}", "el": _r(el[i], 2), "ipp_lat": _r(math.degrees(fi[i]), 3),
               "ipp_lon": _r(math.degrees(la[i]), 3), "stec": _r(stec[i], 2), "vtec": _r(vtec[i], 2),
               "I1": _r(I1[i], 3), "I2": _r(I1[i] * GAMMA, 3), "used": i in used_ep}
        if klob:
            kI = float(cr.klobuchar(np.array([t[i]]), lat0, lon0, np.array([el[i]]), np.array([az[i]]),
                                    klob["alpha"], klob["beta"])[0])
            row["I1_klob"] = _r(kI, 3)
            row["stec_klob"] = _r(kI / TECU_L1, 2)
        rows.append(row)
    summ = {"available": True, "code_rms_m": _r(code_rms, 3), "n_arcs": int(len(arcs)),
            "dcb_rec_ns": _r(dcb_r, 2), "sat_dcb": dcb_src,
            "sat_dcb_values": {f"G{k:02d}": _r(v, 2) for k, v in sorted((dcb or {}).items())
                               if k in set(int(x) for x in P["prn"][m])},
            "vtec_station": _r(xh[0], 2),
            "vtec_rate": _r(xh[1], 2), "model_rms_tecu": _r(float(np.sqrt(np.mean(v ** 2))), 2), "rows": rows}
    if klob and rows:
        diffs = [r["I1"] - r["I1_klob"] for r in rows if r.get("I1") is not None and r.get("I1_klob") is not None]
        if diffs:
            summ["klob_minus"] = {"mean": _r(-np.mean(diffs), 3), "rms": _r(float(np.sqrt(np.mean(np.square(diffs)))), 3)}
    out["summary"] = summ
    ii = np.where(m & np.isfinite(vtec))[0]
    ii = ii[_dec(len(ii), 6000)]
    out["series"] = {"t": t[ii], "prn": P["prn"][ii], "vtec": vtec[ii], "lev": lev[ii], "p4": P4[ii]}
    return out


def _estimate_sat_dcb(lev, m, prn, t, M, fi, la, lat0, lon0, w, ns_tecu, k4):
    """DCB dosyasi yoksa: tum gun ince tabaka modeliyle uydu DCB'lerini kestir (toplam sifir kosulu).

    y = lev/(k4 kappa) = M (sum_k h_k(t) v_k + g_phi dphi + g_lam dlam) - ns_tecu (DCB_s + DCB_r)
    """
    ii = np.where(m & np.isfinite(lev))[0]
    if len(ii) < 200:
        return None
    y = lev[ii] / k4 / TECU_L1
    tt = t[ii]
    nodes = np.arange(tt.min(), tt.max() + 3600.0, 3600.0)
    K = len(nodes)
    sats = np.unique(prn[ii])
    js = np.searchsorted(sats, prn[ii])
    ns = len(sats)
    npar = K + 2 + ns + 1
    A = np.zeros((len(ii) + 1, npar))
    pos = (tt - nodes[0]) / 3600.0
    k0 = np.clip(np.floor(pos).astype(int), 0, K - 2)
    fr = pos - k0
    r = np.arange(len(ii))
    A[r, k0] = M[ii] * (1 - fr)
    A[r, k0 + 1] = M[ii] * fr
    A[r, K] = M[ii] * np.degrees(fi[ii] - lat0)
    A[r, K + 1] = M[ii] * np.degrees(la[ii] - lon0)
    A[r, K + 2 + js] = -ns_tecu
    A[r, -1] = -ns_tecu
    A[-1, K + 2:K + 2 + ns] = 1.0                 # toplam(uydu DCB) = 0
    yy = np.r_[y, 0.0]
    ww = np.r_[w[ii], 1e6]
    try:
        x = np.linalg.solve(A.T @ (A * ww[:, None]), A.T @ (ww * yy))
    except np.linalg.LinAlgError:
        return None
    return {int(s): float(x[K + 2 + k]) for k, s in enumerate(sats)}


def _multipath(P, el, az, valid, idx_ep, cutoff):
    k = 2 / (GAMMA - 1)
    F1m, F2m = LAM1 * P["L1d"], LAM2 * P["L2d"]
    MP1 = P["P1"] - (1 + k) * F1m + k * F2m
    MP2 = P["P2"] - GAMMA * k * F1m + (GAMMA * k - 1) * F2m
    m = P["usable"] & valid
    for a in np.unique(P["arc"][m]):
        mm = m & (P["arc"] == a)
        MP1[mm] -= MP1[mm].mean()
        MP2[mm] -= MP2[mm].mean()
    MP1[~m] = np.nan
    MP2[~m] = np.nan
    rows = []
    for i in idx_ep:
        rows.append({"sv": f"G{P['prn'][i]:02d}", "el": _r(el[i], 2), "az": _r(az[i], 2),
                     "MP1": _r(MP1[i], 3), "MP2": _r(MP2[i], 3)})
    v1 = [r["MP1"] for r in rows if r["MP1"] is not None]
    v2 = [r["MP2"] for r in rows if r["MP2"] is not None]
    bins = [0, 10, 20, 30, 45, 60, 90]
    eb = []
    for a, b in zip(bins[:-1], bins[1:]):
        mm = m & (el > a) & (el <= b)
        eb.append({"bin": f"{a}-{b}", "MP1": _r(np.sqrt(np.nanmean(MP1[mm] ** 2)), 3) if mm.any() else None,
                   "MP2": _r(np.sqrt(np.nanmean(MP2[mm] ** 2)), 3) if mm.any() else None, "n": int(mm.sum())})
    per = []
    for p_ in np.unique(P["prn"][m]):
        mm = m & (P["prn"] == p_)
        per.append({"sv": f"G{p_:02d}", "MP1": _r(np.sqrt(np.nanmean(MP1[mm] ** 2)), 3),
                    "MP2": _r(np.sqrt(np.nanmean(MP2[mm] ** 2)), 3), "n": int(mm.sum())})
    mc = m & (el >= cutoff)
    summ = {"rows": rows, "epoch_rms": [_r(np.sqrt(np.mean(np.square(v1))), 3) if v1 else None,
                                        _r(np.sqrt(np.mean(np.square(v2))), 3) if v2 else None],
            "all_rms": [_r(np.sqrt(np.nanmean(MP1[m] ** 2)), 3), _r(np.sqrt(np.nanmean(MP2[m] ** 2)), 3)],
            "cut_rms": [_r(np.sqrt(np.nanmean(MP1[mc] ** 2)), 3), _r(np.sqrt(np.nanmean(MP2[mc] ** 2)), 3)],
            "bins": eb, "per_sat": per,
            "coef": [_r(1 + k, 4), _r(k, 4), _r(GAMMA * k, 4), _r(GAMMA * k - 1, 4)]}
    ii = np.where(m)[0]
    ii = ii[_dec(len(ii), 8000)]
    return {"summary": summ, "series": {"el": el[ii], "mp1": MP1[ii], "mp2": MP2[ii], "t": P["t"][ii], "prn": P["prn"][ii]}}


def _spp_series(P, sp3, nav, X0, lat0, h0, Rn, Xref, t_ep, step_s, cutoff, met, Pm, Tm_, RH):
    t, prn = P["t"], P["prn"]
    eps = np.unique(t)
    grid = eps[np.abs(((eps - eps[0] + 1e-3) % step_s) - 1e-3) < 1e-3]
    if t_ep not in grid:
        grid = np.sort(np.r_[grid, t_ep])
    src = "broadcast" if nav else "sp3"
    out = []
    epoch_res = None
    for tt in grid:
        m = np.where(np.abs(t - tt) < 1e-3)[0]
        if len(m) < 5:
            continue
        P3 = A3 * P["P1"][m] + B3 * P["P2"][m]
        Xs = np.full((len(m), 3), np.nan)
        cd = np.full(len(m), np.nan)
        for k, i in enumerate(m):
            p_ = int(prn[i])
            tx = tt - P3[k] / C
            if nav:
                e = orbit.select_eph(nav["eph"], p_, tx)
                if e is None:
                    continue
                c0 = orbit.bc_clock(e, tx)
                pos, drel = orbit.bc_position(e, tx - c0)
                Xs[k], cd[k] = pos, (c0 + drel) * C        # P3 icin TGD gerekmez
            else:
                c0 = orbit.sp3_clock(sp3, p_, [tx])[0]
                pos, vel = orbit.sp3_interp(sp3, p_, [tx - c0])
                Xs[k] = pos[0]
                cd[k] = (c0 - 2 * float(pos[0] @ vel[0]) / C ** 2) * C
        ok = np.isfinite(Xs).all(axis=1) & np.isfinite(cd) & np.isfinite(P3)
        if ok.sum() < 5:
            continue
        e_, _ = geo.elaz(X0, Xs[ok])
        keep = e_ >= cutoff
        if keep.sum() < 4:
            continue
        Pe = float(np.interp(tt, t, Pm))
        Te = float(np.interp(tt, t, Tm_))
        He = float(np.interp(tt, t, RH))
        zh, zw = cr.saastamoinen(Pe, Te, He, lat0, h0)
        mh, mw = cr.herring(e_[keep], lat0, h0, Te)
        Tm = zh * mh + zw * mw
        r = est.spp_epoch(P3[ok][keep], Xs[ok][keep], cd[ok][keep], Tm, X0=X0)
        if r is None:
            continue
        d = Rn @ (r["X"] - Xref)
        row = {"t": float(tt), "N": float(d[0]), "E": float(d[1]), "U": float(d[2]), "n": r["n"], "rms": r["rms"],
               "cdt": r["cdt"]}
        out.append(row)
        if abs(tt - t_ep) < 1e-3:
            epoch_res = {"X": _lst(r["X"]), "clock_ns": _r(r["cdt"] / C * 1e9, 1), "n": r["n"], "rms": _r(r["rms"], 2)}
    res = {"source": src, "n": len(out), "epoch": epoch_res}
    if out:
        N = np.array([o["N"] for o in out])
        E = np.array([o["E"] for o in out])
        U = np.array([o["U"] for o in out])
        d3 = np.sqrt(N ** 2 + E ** 2 + U ** 2)
        good = d3 < 50
        res["n_bad"] = int((~good).sum())
        meanNEU = np.array([N[good].mean(), E[good].mean(), U[good].mean()])
        res["mean"] = _lst(Xref + Rn.T @ meanNEU)
        res["rms3d"] = _r(float(np.sqrt(np.mean(d3[good] ** 2))), 2)
        res["series"] = out
    return res


def _plots(P, t, prn, el, az, valid, sel, sol, idx_ep, t_ep, Ttrop, wu, tide_up, iono, mp, spp, day0):
    hrs = lambda x: _lst((np.asarray(x) - day0) / 3600.0, 4)
    eps_prn = sorted({int(prn[i]) for i in idx_ep})
    # gokyuzu izleri
    sky = {}
    for p_ in np.unique(prn):
        m = np.where((prn == p_) & valid & (el > 0))[0]
        m = m[_dec(len(m), 240)]
        sky[f"G{p_:02d}"] = {"az": _lst(az[m], 1), "el": _lst(el[m], 1)}
    sky_ep = [{"sv": f"G{prn[i]:02d}", "az": _r(az[i], 1), "el": _r(el[i], 1), "used": bool(sel[i])} for i in idx_ep]
    # GF / MW (epok uydularinin epoktaki yaylari)
    gfmw = {}
    for i in idx_ep:
        a = P["arc"][i]
        if a < 0:
            continue
        m = np.where(P["arc"] == a)[0]
        m = m[_dec(len(m), 600)]
        gfmw[f"G{prn[i]:02d}"] = {"t": hrs(t[m]), "gf": _lst(P["GF"][m] - P["GF"][m[0]], 3),
                                  "mw": _lst(P["MW"][m] - np.mean(P["MW"][m]), 3)}
    # saat sicramasi: epoktaki uydulardan, yayi en cok sicrama iceren (esitse en uzun) yay
    cj = None
    best = None
    for i in idx_ep:
        a = P["arc"][i]
        if a < 0:
            continue
        m = np.where(P["arc"] == a)[0]
        nj = len(np.unique(P["K"][m]))
        key = (nj, len(m))
        if best is None or key > best[0]:
            best = (key, i, m)
    if best is not None:
        _, i, m = best
        m = m[_dec(len(m), 1500)]
        raw = LAM1 * P["L1"][m] - P["P1"][m]
        cor = LAM1 * P["L1d"][m] - P["P1"][m]
        cj = {"sv": f"G{prn[i]:02d}", "t": hrs(t[m]), "raw_km": _lst((raw - raw[0]) / 1e3, 4),
              "cor_m": _lst(cor - np.median(cor), 3), "n_jumps": int(best[0][0] - 1)}
    # L3 duzeltmeleri, saat, ZTD
    ts = t[sel]
    vv = sol["v_phase"]
    ii = _dec(len(ts), 5000)
    l3 = {"t": hrs(ts[ii]), "v_mm": _lst(vv[ii] * 1e3, 2), "prn": [int(x) for x in prn[sel][ii]],
          "clk_t": hrs(sol["epochs"]), "clk_us": _lst(sol["clock"] / C * 1e6, 3),
          "zwd_t": hrs(sol["zwd_blocks"]), "zwd_mm": _lst(sol["zwd"] * 1e3, 1)}
    # troposfer - yukseklik
    m = np.where(valid & (el > 0) & np.isfinite(Ttrop))[0]
    m = m[_dec(len(m), 4000)]
    trop = {"el": _lst(el[m], 2), "T": _lst(Ttrop[m], 3)}
    # gelgit
    ue = np.unique(t)
    tu = np.interp(ue, t, tide_up)
    k = _dec(len(ue), 600)
    tide = {"t": hrs(ue[k]), "up_cm": _lst(tu[k] * 100, 2)}
    out = {"sky": sky, "sky_epoch": sky_ep, "gfmw": gfmw, "clockjump": cj, "l3": l3, "trop": trop, "tide": tide,
           "epoch_h": _r((t_ep - day0) / 3600.0, 4)}
    if iono.get("series"):
        s = iono["series"]
        out["iono"] = {"t": hrs(s["t"]), "prn": [int(x) for x in s["prn"]], "vtec": _lst(s["vtec"], 2),
                       "lev": _lst(s["lev"], 3), "p4": _lst(s["p4"], 3)}
    s = mp["series"]
    out["mp"] = {"el": _lst(s["el"], 1), "mp1": _lst(s["mp1"], 3), "mp2": _lst(s["mp2"], 3),
                 "t": hrs(s["t"]), "prn": [int(x) for x in s["prn"]]}
    if spp.get("series"):
        S = spp["series"]
        out["spp"] = {"t": hrs([o["t"] for o in S]), "N": _lst([o["N"] for o in S], 3),
                      "E": _lst([o["E"] for o in S], 3), "U": _lst([o["U"] for o in S], 3)}
    return out
