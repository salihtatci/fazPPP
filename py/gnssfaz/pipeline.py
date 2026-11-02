"""Tum hesap zinciri: dosyalar -> sonuc sozlugu (JSON'a donusturulebilir)."""
import gzip
import math
import time
from datetime import datetime

import numpy as np

from . import readers, preprocess as pre, orbit, corrections as cr, estimate as est, geo, codefilter as cf
from .const import (C, F1, F2, F7, LAM1, LAM2, A3, B3, GAMMA, LAM_NL, TECU_L1, OMEGA_E, MS, sv_name)
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


def _load_atx(default_atx, when):
    """Gomulu ANTEX secimi: 27.11.2022 (GPS haftasi 2238) sonrasi igs20, oncesi igs08."""
    if not default_atx:
        return None, None
    if isinstance(default_atx, str):
        default_atx = {"igs08": default_atx}
    key = "igs20" if (when >= datetime(2022, 11, 27) and default_atx.get("igs20")) else "igs08"
    path = default_atx.get(key) or next(iter(default_atx.values()))
    with open(path, "rb") as fh:
        raw = fh.read()
    return (gzip.decompress(raw) if path.endswith(".gz") else raw).decode("latin-1"), "bundled-" + key


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
    _m = str(params.get("mode") or "auto").lower()
    mode_req = "kinematic" if _m.startswith("kin") else ("static" if _m.startswith("sta") else "auto")
    mode = "static" if mode_req == "auto" else mode_req
    systems = "".join(c for c in str(params.get("systems") or "GE").upper() if c in "GE") or "G"
    if "G" not in systems:
        systems = "G" + systems
    rate_req = params.get("rate")
    rate_req = float(rate_req) if rate_req not in (None, "", 0) else None
    min_arc_req = params.get("min_arc_min")
    min_arc_req = float(min_arc_req) if min_arc_req not in (None, "", 0) else None
    code_kf = params.get("code_kf", True) not in (False, "0", 0, "false", "no")
    code_kf_tau = float(params.get("code_kf_tau") or 100.0)
    code_kf_thr = 5.0

    # ------------------------------------------------------------------ 1. okuma
    step("read", 0.02)
    obs = readers.read_obs(files["obs"], systems=systems)
    H = obs["header"]
    sp3 = readers.read_sp3(files["sp3"])
    t_first = float(obs["t"].min())
    t_last = float(obs["t"].max())
    clk = None
    if files.get("clk"):
        clk = readers.read_clk(files["clk"], t_first, t_last)
    nav = readers.read_nav(files["nav"]) if files.get("nav") else None
    met = readers.read_met(files["met"]) if files.get("met") else None
    dcb = readers.read_dcb(files["dcb"]) if files.get("dcb") else None
    day0 = day_start(t_first)
    when = from_gps(t_first)
    atx_text = files.get("atx")
    atx_src = "user"
    if not atx_text and default_atx:
        atx_text, atx_src = _load_atx(default_atx, when)
        if atx_src == "bundled-igs08" and when > datetime(2017, 1, 28):
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
    skipped = [x for x in H.get("systems_skipped", []) if x not in ("E:nodual",)]
    if skipped:
        warn.append("sys_skipped:" + ",".join(skipped))
    if "E:nodual" in H.get("systems_skipped", []):
        warn.append("gal_nodual")
    if clk is not None:
        cc = _clk_sp3_check(clk, sp3, t_first, t_last)
        if cc is not None and cc > 0.3:
            warn.append(f"clk_sp3_mismatch:{cc:.2f} m")
    out_sp3 = float(np.mean((obs["t"] < sp3["t"][0]) | (obs["t"] > sp3["t"][-1])))
    if out_sp3 > 0.02:
        warn.append(f"sp3_partial:{100 * out_sp3:.0f}%")
    if clk is not None:
        out_clk = float(np.mean((obs["t"] < clk["t"][0]) | (obs["t"] > clk["t"][-1])))
        if out_clk > 0.02:
            warn.append(f"clk_partial:{100 * out_clk:.0f}%")

    # ------------------------------------------------------------------ 2. on isleme
    step("preprocess", 0.10)
    span_min = (t_last - t_first) / 60.0
    min_arc = min_arc_req if min_arc_req is not None else float(min(10.0, max(2.0, 0.3 * span_min)))
    if min_arc_req is None and min_arc < 10.0:
        warn.append(f"short_file:{min_arc:.1f}")
    P = pre.preprocess(obs, min_arc_sec=min_arc * 60.0)
    if P["stats"].get("lli_ignored"):
        warn.append("lli_ignored")
    t, prn = P["t"], P["prn"]
    n = len(t)
    if n == 0:
        raise Hata("no_dual_freq|Dort olcusu (P1, P2, L1, L2) tam olan kayit yok; cift frekans gerekli.")
    isE = P["isE"]
    isG = ~isE

    # ------------------------------------------------------------------ 2b. kod olculerine Kalman filtresi
    step("codefilter", 0.15)
    P["P1raw"], P["P2raw"] = P["P1"], P["P2"]
    KF = None
    if code_kf:
        KF = cf.kalman_code_filter(P, tau=code_kf_tau, sig_code=sig_code, thr=code_kf_thr)
        P["P1"], P["P2"] = KF["P1f"], KF["P2f"]
    f1, f2 = P["f1"], P["f2"]
    lam1, lam2 = C / f1, C / f2
    den = f1 ** 2 - f2 ** 2
    a3, b3 = f1 ** 2 / den, -f2 ** 2 / den
    gam = (f1 / f2) ** 2
    lam_nl = C / (f1 + f2)

    # ------------------------------------------------------------------ 3. yaklasik konum
    X0 = np.array(H["approx"], float)
    x0_src = "header"
    Xspp = _spp_first(P, sp3)
    if not np.all(np.isfinite(X0)) or np.linalg.norm(X0) < 6.0e6:
        X0 = Xspp
        x0_src = "spp"
        if X0 is None:
            raise Hata("no_approx|Yaklasik konum bulunamadi (RINEX basliginda yok, kod cozumu de basarisiz).")
    elif Xspp is not None and np.linalg.norm(Xspp - X0) > 1000.0:
        warn.append(f"approx_bad:{np.linalg.norm(Xspp - X0) / 1000:.1f} km")
        X0 = Xspp
        x0_src = "spp"
    lat0, lon0, h0 = geo.xyz2llh(*X0)
    Rn = geo.rot_neu(lat0, lon0)

    # ------------------------------------------------------------------ 4. uydu konum ve saati
    step("orbit", 0.20)
    clk_src = clk if clk is not None else sp3
    Xs = np.full((n, 3), np.nan)
    Vs = np.full((n, 3), np.nan)
    dts = np.full(n, np.nan)
    for p_ in np.unique(prn):
        m = prn == p_
        t_sat = t[m] - P["P1"][m] / C
        c0 = orbit.sp3_clock(clk_src, int(p_), t_sat)
        ttx = t_sat - c0
        pos, vel = orbit.sp3_interp(sp3, int(p_), ttx)
        Xs[m], Vs[m], dts[m] = pos, vel, c0
    rel = -2 * np.einsum("ij,ij->i", Xs, Vs) / C ** 2
    dts_tot = dts + rel
    valid = np.isfinite(Xs).all(axis=1) & np.isfinite(dts)
    if valid.sum() == 0:
        raise Hata("sp3_no_overlap|SP3 dosyasi gozlem gunuyle ortusmuyor veya uydu saati yok.")
    miss_orb = sorted({int(p_) for p_ in np.unique(prn[~valid])} - {int(p_) for p_ in np.unique(prn[valid])})
    if miss_orb:
        warn.append("sat_no_orbit:" + ",".join(sv_name(p_) for p_ in miss_orb))
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
    L1m = lam1 * P["L1d"]
    L2m = lam2 * P["L2d"]
    L3 = a3 * L1m + b3 * L2m
    P3 = a3 * P["P1"] + b3 * P["P2"]
    f2E = (H.get("f2") or {}).get("E")
    e2 = "E07" if (f2E and abs(f2E - F7) < 1) else "E05"
    aE = F1 ** 2 / (F1 ** 2 - f2E ** 2) if f2E else 0.0
    bE = -f2E ** 2 / (F1 ** 2 - f2E ** 2) if f2E else 0.0
    sun, moon = geo.sun_moon_ecef(t)
    ex = ey = ez = None
    d_sat = np.zeros(n)
    nadir = np.full(n, np.nan)
    if valid.any():
        ex, ey, ez = cr.sat_body_axes(np.where(valid[:, None], Xs, 1.0), sun)
        coef3 = {"G": {"G01": A3, "G02": B3}, "E": {"E01": aE, e2: bE}}
        ds_, nad_, miss = cr.satellite_antenna(u[valid], prn[valid], ex[valid], ey[valid], ez[valid], atx["sat"], coef3)
        d_sat[valid] = ds_
        nadir[valid] = nad_
        if atx_text and miss:
            warn.append("sat_ant_missing:" + ",".join(sv_name(p_) for p_ in sorted(miss)))
    elv = np.where(valid, el, 90.0)
    uv = np.where(valid[:, None], u, 0.0)
    d_rec3 = cr.receiver_antenna(uv, elv, atx["rec"], H["delta_hen"], Rn, {"G01": A3, "G02": B3})
    if isE.any():
        recE = None
        if atx["rec"] is not None:
            r_ = atx["rec"]
            recE = {"zen": r_["zen"], "E01": r_.get("E01", r_["G01"]), e2: r_.get(e2, r_.get("E05", r_["G02"]))}
        d_rec3E = cr.receiver_antenna(uv, elv, recE, H["delta_hen"], Rn, {"E01": aE, e2: bE})
        d_rec3 = np.where(isE, d_rec3E, d_rec3)
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
    d_wu3 = lam_nl * wu
    common = d_sat + d_tide + d_shap
    corr_P3 = d_rec3 + common
    corr_L3 = corr_P3 + d_wu3

    # ------------------------------------------------------------------ 7. secim ve epok
    sp3_int = sp3["interval"]
    clk_int = clk["interval"] if clk is not None else None
    if rate_req:
        rate = max(rate_req, float(P["interval"]))
        ref0 = day0
        if clk is None and rate < sp3_int - 1e-6:
            warn.append("clk_from_sp3")
    elif clk is not None:
        rate = max(clk_int, float(P["interval"]))
        ref0 = day0
    else:
        rate = sp3_int
        ref0 = sp3["t"][0]
    on_grid, offgrid = _grid_select(t, rate, ref0, float(P["interval"]), clk is not None or bool(rate_req))
    if offgrid:
        warn.append("epochs_offgrid")
    base = valid & P["usable"] & (el >= cutoff) & np.isfinite(Ttrop)

    def make_sel(min_sat, extra=None):
        s_ = base & on_grid
        if extra is not None:
            s_ &= extra
        for _ in range(3):
            arcs_sel, cnt = np.unique(P["arc"][s_], return_counts=True)
            s_ &= np.isin(P["arc"], arcs_sel[cnt >= 3])
            if min_sat > 1:
                ep_u, ep_c = np.unique(t[s_], return_counts=True)
                s_ &= np.isin(t, ep_u[ep_c >= min_sat])
        return s_

    def subset(mask):
        return {"t": t[mask], "arc": P["arc"][mask], "Xs": Xs[mask], "el": el[mask], "cdts": C * dts_tot[mask],
                "Tmodel": Ttrop[mask], "mw": mwet[mask], "isE": isE[mask]}

    def pick_epoch(s_):
        ep_s = np.unique(t[s_])
        req = params.get("epoch")
        if req:
            hh, mm, ss = (list(map(float, str(req).split(":"))) + [0, 0])[:3]
            tr = day0 + hh * 3600 + mm * 60 + ss
        else:
            tr = ep_s[len(ep_s) // 2]
        return tr, float(ep_s[np.argmin(np.abs(ep_s - tr))])

    def solve(s_, md, phase=None, code=None, cph=None, ccd=None):
        _, te = pick_epoch(s_)
        return est.phase_solution(subset(s_), L3[s_] if phase is None else phase, P3[s_] if code is None else code,
                                  corr_L3[s_] if cph is None else cph, corr_P3[s_] if ccd is None else ccd,
                                  X0, te, Rn, sig_phase, sig_code, zwd_block, mode=md)

    def solve_clean(md):
        """Cozum + artiga dayali ayiklama (kacmis faz kopmasinda yay bolme, aykiri cikarma), en cok 5 tur."""
        ms = 4 if md == "kinematic" else 1
        rem = np.zeros(n, bool)
        cl = {"split": 0, "removed_phase": 0, "removed_code": 0, "iterations": 0}
        s_ = make_sel(ms)
        if s_.sum() < 20:
            return None, s_, cl, rem
        so = None
        for it in range(5):
            so = solve(s_, md)
            cl["iterations"] = it + 1
            if it == 4 or so["s0"] > 20.0:          # kaba model hatasi (ornegin hareketli aliciya statik): ayiklama yok
                break
            if not _clean_residuals(so, s_, P, t, el, sig_phase, sig_code, rem, cl):
                break
            s_new = make_sel(ms, extra=~rem)
            if s_new.sum() < 20:
                break
            s_ = s_new
        return so, s_, cl, rem

    # ------------------------------------------------------------------ 7b/8. statik mi kinematik mi? cozum ve ayiklama
    step("l3", 0.50)
    mode_info = {"requested": mode_req, "chosen": None, "reason": None}
    sol = None
    if mode_req == "auto":
        motion = _motion_test(P, sp3, nav, X0, lat0, h0, Rn, cutoff, met, Pm, Tm_, RH)
        mode_info["spp_spread_m"] = _r(motion, 2) if motion is not None else None
        if motion is not None and motion > 10.0:
            mode, mode_info["reason"] = "kinematic", "spp_motion"
        else:
            mode, mode_info["reason"] = "static", "static_fit"
            st0, s_st, cl_st, rem_st = solve_clean("static")
            if st0 is not None:
                mode_info["s0_static"] = _r(st0["s0"], 3)
                sol, sel, clean, removed = st0, s_st, cl_st, rem_st
                if st0["s0"] > 2.5 or st0["rms_phase"] > 0.04:
                    k0, s_k, cl_k, rem_k = solve_clean("kinematic")
                    if k0 is not None:
                        mode_info["s0_kinematic"] = _r(k0["s0"], 3)
                        if k0["s0"] < 0.5 * st0["s0"]:
                            mode, mode_info["reason"] = "kinematic", "residuals"
                            sol, sel, clean, removed = k0, s_k, cl_k, rem_k
        if mode == "kinematic":
            warn.append("auto_kin")
    mode_info["chosen"] = mode
    if clk is None and mode == "kinematic":
        warn.append("kin_no_clk")
    min_sat = 4 if mode == "kinematic" else 1
    if sol is None or (sol.get("mode") != mode):
        sol, sel, clean, removed = solve_clean(mode)
    if sol is None:
        raise Hata("too_few_obs|Cozum icin yeterli olcu yok (SP3 epoklarinda, yukseklik sinirinin uzerinde).")
    if clean["split"] or clean["removed_phase"] or clean["removed_code"]:
        warn.append("cleaned:" + f"{clean['split']}/{clean['removed_phase']}/{clean['removed_code']}")
    if mode == "static" and sol["s0"] > 3.0:
        warn.append("static_bad_fit")
    elif sol["s0"] > 3.0:
        warn.append("bad_fit")
    t_req, t_ep = pick_epoch(sel)
    if abs(t_ep - t_req) > 0.5:
        warn.append("epoch_snapped")
    dsel = subset(sel)
    sol_no4b = solve(sel, mode, cph=np.zeros(sel.sum()), ccd=np.zeros(sel.sum()))
    sol_raw = None
    if KF is not None:
        P3raw = a3 * P["P1raw"] + b3 * P["P2raw"]
        sol_raw = solve(sel, mode, code=P3raw[sel])
    Xst = sol["X"]
    latS, lonS, hS = geo.xyz2llh(*Xst)

    ref = params.get("ref_xyz")
    ref_src = "user"
    if ref and len(ref) == 3 and all(r is not None for r in ref):
        Xref = np.array(ref, float)
    else:
        Xref = Xst.copy()
        ref_src = "static" if mode == "static" else "epoch"

    def neu(X):
        if X is None:
            return None
        d = Rn @ (np.asarray(X) - Xref)
        return {"N": _r(d[0]), "E": _r(d[1]), "U": _r(d[2]), "H": _r(math.hypot(d[0], d[1])),
                "D3": _r(float(np.linalg.norm(d)))}

    j_ep = int(np.searchsorted(sol["epochs"], t_ep - 1e-3))
    k_sel = np.argmin(np.abs(t[sel] - t_ep))
    zwd_add = sol["zwd_kin"] if (sol["has_kin"] or mode == "kinematic") else 0.0
    ztd_ep = float(ZHD[sel][k_sel] + ZWD[sel][k_sel] + zwd_add)
    ztd_model_ep = float(ZHD[sel][k_sel] + ZWD[sel][k_sel])
    kin_like = sol["has_kin"] or mode == "kinematic"

    # ------------------------------------------------------------------ 9. epok tablolari
    step("tables", 0.60)
    me = (np.abs(t - t_ep) < 1e-3)
    idx_ep = np.where(me)[0]
    idx_ep = idx_ep[np.argsort(-np.nan_to_num(el[idx_ep], nan=-99))]
    used_ep = set(np.where(me & sel)[0].tolist())
    sat_rows, atm_rows, cor_rows, a5_rows = [], [], [], []
    for i in idx_ep:
        ok = valid[i]
        sv = sv_name(prn[i])
        sat_rows.append({"sv": sv, "X": _r(Xs[i, 0], 3), "Y": _r(Xs[i, 1], 3), "Z": _r(Xs[i, 2], 3),
                         "rho": _r(rho0[i], 3), "tau_ms": _r(tau[i] * 1e3, 3), "dts_us": _r(dts[i] * 1e6, 3),
                         "rel_ns": _r(rel[i] * 1e9, 3), "sagnac": _r(sagnac[i], 3), "el": _r(el[i], 2),
                         "az": _r(az[i], 2), "used": i in used_ep, "valid": bool(ok)})
        atm_rows.append({"sv": sv, "el": _r(el[i], 2), "mh": _r(mh[i], 3), "mw": _r(mwet[i], 3),
                         "T": _r(Ttrop[i], 3), "I1_code": _r((P["P2"][i] - P["P1"][i]) / (gam[i] - 1), 3),
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
    used_g = [i for i in sorted(used_ep) if isG[i]]
    if len(used_g):
        ii = np.array(used_g)
        k1 = P["P1"][ii] - rho0[ii] + C * dts_tot[ii]
        k3 = P3[ii] - rho0[ii] + C * dts_tot[ii]
        k3t = k3 - Ttrop[ii]
        check = {"spread_P1": _r(np.ptp(k1), 3), "spread_P3": _r(np.ptp(k3), 3), "spread_P3T": _r(np.ptp(k3t), 3),
                 "clock_rough_ms": _r(np.mean(k1) / MS, 5)}

    # ------------------------------------------------------------------ 10. iyonosfer ve TEC (GPS)
    step("iono", 0.68)
    iono = _ionosphere(P, el, az, valid & isG, cutoff, lat0, lon0, t, t_ep, dcb, nav, sel, idx_ep, used_ep, warn)

    # ------------------------------------------------------------------ 11. tek frekans cozumleri (GPS)
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
    sel_g = sel & isG
    if min_sat > 1:
        for _ in range(3):
            ep_u, ep_c = np.unique(t[sel_g], return_counts=True)
            sel_g &= np.isin(t, ep_u[ep_c >= 4])
            arcs_sel, cnt = np.unique(P["arc"][sel_g], return_counts=True)
            sel_g &= np.isin(P["arc"], arcs_sel[cnt >= 3])
    dsel_g = subset(sel_g)
    variants = [("L1", "none", None)]
    if I1_mod is not None:
        variants += [("L1", "klobuchar", I1_mod), ("L2", "klobuchar", I1_mod)]
    if sel_g.sum() < 20:
        variants = []
    for fr, model, Iv in variants:
        g = 1.0 if fr == "L1" else GAMMA
        Isel = (Iv[sel_g] * g) if Iv is not None else np.zeros(sel_g.sum())
        phase = (L1m if fr == "L1" else L2m)[sel_g]
        code = (P["P1"] if fr == "L1" else P["P2"])[sel_g]
        drec = (d_rec1 if fr == "L1" else d_rec2)[sel_g]
        kd = (B3 * dcbm if fr == "L1" else -A3 * dcbm)[sel_g]
        wl = (LAM1 if fr == "L1" else LAM2) * wu[sel_g]
        try:
            s1 = est.phase_solution(dsel_g, phase, code, common[sel_g] + drec + wl - Isel,
                                    common[sel_g] + drec + Isel + kd, X0, t_ep, Rn, sig_phase, sig_code, zwd_block,
                                    mode=mode)
        except np.linalg.LinAlgError:
            continue
        jj = int(np.searchsorted(s1["epochs"], t_ep - 1e-3))
        jj = min(jj, len(s1["epochs"]) - 1)
        if mode == "kinematic":
            Xk1, neu_st = s1["X"], None
            ztd1 = ztd_model_ep + s1.get("zwd_kin", 0.0)
        else:
            Xk1, neu_st = (s1["Xk"] if s1["has_kin"] else None), neu(s1["X"])
            ztd1 = ztd_model_ep + (s1["zwd_kin"] if s1["has_kin"] else 0)
        single.append({"freq": fr, "model": model, "X": _lst(s1["X"]), "Xk": _lst(Xk1) if Xk1 is not None else None,
                       "neu_static": neu_st, "neu_kin": neu(Xk1), "rms_phase_mm": _r(s1["rms_phase"] * 1e3, 1),
                       "clock_ns": _r(s1["clock"][jj] / C * 1e9, 1), "ztd": _r(ztd1, 4)})
    if mode == "kinematic":
        l3row = {"freq": "L3", "model": "iono-free", "X": _lst(Xst), "Xk": _lst(Xst), "neu_static": None,
                 "neu_kin": neu(Xst), "rms_phase_mm": _r(sol["rms_phase"] * 1e3, 1),
                 "clock_ns": _r(sol["clock"][j_ep] / C * 1e9, 1), "ztd": _r(ztd_ep, 4)}
    else:
        l3row = {"freq": "L3", "model": "iono-free", "X": _lst(Xst), "Xk": _lst(sol["Xk"]) if sol["has_kin"] else None,
                 "neu_static": neu(Xst), "neu_kin": neu(sol["Xk"]), "rms_phase_mm": _r(sol["rms_phase"] * 1e3, 1),
                 "clock_ns": _r(sol["clock"][j_ep] / C * 1e9, 1), "ztd": _r(ztd_ep, 4)}
    single.append(l3row)

    # ------------------------------------------------------------------ 12. cok yolluluk (GPS)
    step("multipath", 0.82)
    mp = _multipath(P, el, az, valid & isG, idx_ep, cutoff)

    # ------------------------------------------------------------------ 13. kod SPP (GPS)
    step("spp", 0.88)
    spp = _spp_series(P, sp3, nav, X0, lat0, h0, Rn, Xref, t_ep, spp_step, cutoff, met, Pm, Tm_, RH)
    kfres = {"enabled": False}
    if KF is not None:
        kfres = _codefilter_report(P, KF, sel, idx_ep, sol, sol_raw, mode, Rn, Xref, Xst, sp3, nav, X0, lat0, h0,
                                   t_ep, spp, spp_step, cutoff, met, Pm, Tm_, RH, rate, day0, code_kf_tau, code_kf_thr)
    if mode == "static" and spp.get("series") and len(spp["series"]) >= 5:
        Ns = np.array([o["N"] for o in spp["series"]])
        Es = np.array([o["E"] for o in spp["series"]])
        hor = np.hypot(Ns - np.median(Ns), Es - np.median(Es))
        if np.percentile(hor, 90) > 15.0:
            warn.append("rx_moving")

    # ------------------------------------------------------------------ 14. grafik verileri
    step("plots", 0.94)
    plots = _plots(P, t, prn, el, az, valid, sel, sol, idx_ep, t_ep, Ttrop, wu, tide_up, iono, mp, spp, day0)
    kin = _kin_track(sol, Rn, Xref, day0, rate) if mode == "kinematic" else None

    # ------------------------------------------------------------------ 15. ozet
    if mode == "kinematic":
        solutions = [
            {"name": "l3_kin_epoch", "X": _lst(Xst), "sX_mm": _lst(sol["sX"] * 1e3, 1), "neu": neu(Xst)},
            {"name": "l3_kin_no4b", "X": _lst(sol_no4b["X"]), "sX_mm": _lst(sol_no4b["sX"] * 1e3, 1),
             "neu": neu(sol_no4b["X"])},
        ]
    else:
        solutions = [
            {"name": "l3_static", "X": _lst(Xst), "sX_mm": _lst(sol["sX"] * 1e3, 1), "neu": neu(Xst)},
            {"name": "l3_kin", "X": _lst(sol["Xk"]) if sol["has_kin"] else None,
             "sX_mm": _lst(sol["sXk"] * 1e3, 1) if sol["has_kin"] else None, "neu": neu(sol["Xk"])},
            {"name": "l3_static_no4b", "X": _lst(sol_no4b["X"]), "sX_mm": _lst(sol_no4b["sX"] * 1e3, 1),
             "neu": neu(sol_no4b["X"])},
        ]
    for s1 in single[:-1]:
        solutions.append({"name": f"{s1['freq'].lower()}_{s1['model']}_kin", "X": s1["Xk"], "sX_mm": None,
                          "neu": s1["neu_kin"]})
    if spp.get("epoch"):
        solutions.append({"name": "spp_epoch", "X": spp["epoch"]["X"], "sX_mm": None, "neu": neu(spp["epoch"]["X"])})
    if spp.get("mean") and mode == "static":
        solutions.append({"name": "spp_mean", "X": spp["mean"], "sX_mm": None, "neu": neu(spp["mean"])})
    solutions.append({"name": "header", "X": _lst(H["approx"]), "sX_mm": None,
                      "neu": neu(H["approx"]) if np.linalg.norm(H["approx"]) > 6e6 else None})

    n_gal = int(len(np.unique(prn[isE])))
    res = {
        "meta": {"version": "1.3.0", "runtime_s": _r(time.time() - T0, 1), "warnings": warn, "mode_auto": mode_info,
                 "clean": clean,
                 "atx_source": atx_src if atx_text else None, "rec_ant_found": atx["rec_name"],
                 "x0_source": x0_src, "ref_source": ref_src, "ref_xyz": _lst(Xref),
                 "params": {"cutoff": cutoff, "sig_phase": sig_phase, "sig_code": sig_code,
                            "zwd_block_h": zwd_block / 3600, "spp_step": spp_step, "mode": mode,
                            "systems": systems, "rate": rate, "min_arc_min": min_arc, "mode_req": mode_req,
                            "code_kf": bool(code_kf), "code_kf_tau": code_kf_tau}},
        "station": {"marker": H["marker"], "receiver": H["receiver"], "antenna": H["antenna"].strip(),
                    "approx": _lst(H["approx"]), "delta_hen": _lst(H["delta_hen"]),
                    "obs_types": H["obs_types"], "used_codes": H["used_codes"], "rinex_version": H["version"],
                    "obs_types_E": H.get("obs_types_E"), "used_codes_E": H.get("used_codes_E"),
                    "interval": P["interval"], "first": fmt(t_first, True), "last": fmt(t_last, True),
                    "n_epochs": int(len(np.unique(obs["t"]))), "n_sats": int(len(np.unique(obs["prn"]))),
                    "n_gps": int(len(np.unique(prn[isG]))), "n_gal": n_gal,
                    "lat": geo.dms(latS), "lon": geo.dms(lonS), "h": _r(hS, 4),
                    "lat_deg": _r(math.degrees(lat0), 6), "lon_deg": _r(math.degrees(lon0), 6), "h0": _r(h0, 3)},
        "epoch": {"requested": fmt(t_req), "used": fmt(t_ep), "used_full": fmt(t_ep, True), "n_sat": int(len(idx_ep)),
                  "n_used": len(used_ep), "clock_jumps_ms": int(k_ep)},
        "data": {"sp3_interval": sp3_int, "sp3_first": fmt(sp3["t"][0], True), "sp3_last": fmt(sp3["t"][-1], True),
                 "clk": ({"interval": clk_int, "first": fmt(clk["t"][0], True), "last": fmt(clk["t"][-1], True),
                          "n_sat": len(clk["prns"]), "ac": clk.get("ac")} if clk is not None else None),
                 "sat_clock": "clk" if clk is not None else "sp3", "rate": rate,
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
        "l3": {"mode": mode, "X": _lst(Xst), "sX_mm": _lst(sol["sX"] * 1e3, 1),
               "Xk": _lst(sol["Xk"]) if sol["has_kin"] else None,
               "sXk_mm": _lst(sol["sXk"] * 1e3, 1) if sol["has_kin"] else None,
               "kin_minus_static_neu_mm": _lst(Rn @ (sol["Xk"] - Xst) * 1e3, 1) if sol["has_kin"] else None,
               "clock_m": _r(sol["clock_kin"] if kin_like else sol["clock"][j_ep], 3),
               "clock_ns": _r((sol["clock_kin"] if kin_like else sol["clock"][j_ep]) / C * 1e9, 1),
               "sclock_m": _r(sol.get("sclock_kin", 0), 3), "ztd": _r(ztd_ep, 4), "ztd_model": _r(ztd_model_ep, 4),
               "isb_ns": _r(sol["isb"] / C * 1e9, 2) if sol.get("isb") is not None else None,
               "sisb_ns": _r(sol["sisb"] / C * 1e9, 2) if sol.get("sisb") is not None else None,
               "n_gal_used": int(len(np.unique(prn[sel & isE]))), "n_gps_used": int(len(np.unique(prn[sel & isG]))),
               "dop": {k: _r(v, 3) for k, v in sol.get("dop", {}).items()}, "n_sat_kin": sol.get("n_sat_kin"),
               "s0": _r(sol["s0"], 3), "dof": sol["dof"], "n_obs": sol["n_obs"], "n_par": sol["n_par"],
               "n_epochs": sol["n_epochs"], "n_arcs": sol["n_arcs"], "n_zwd": sol["n_zwd"],
               "rms_phase_mm": _r(sol["rms_phase"] * 1e3, 1), "rms_code_m": _r(sol["rms_code"], 2),
               "iter_mm": _lst(np.array(sol["iter"]) * 1e3, 1),
               "no4b_static_minus_full_neu_mm": _lst(Rn @ (sol_no4b["X"] - Xst) * 1e3, 1),
               "no4b_rms_phase_mm": _r(sol_no4b["rms_phase"] * 1e3, 1)},
        "kin": kin,
        "codefilter": kfres,
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
        m = (t == tt) & ~P["isE"]
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
        if info is None and a >= 0:                    # artiga dayali ayiklamada bolunmus yay
            mm = np.where(P["arc"] == a)[0]
            if len(mm):
                tm = P["t"][mm]
                info = {"arc": a, "start": float(tm.min()), "end": float(tm.max()), "n": int(len(mm)),
                        "mw_std": float(np.std(P["MW"][mm])) if len(mm) > 1 else None, "usable": True}
        if info is None:
            rows.append({"sv": sv_name(prn[i]), "start": None, "end": None, "n": 0, "min": 0,
                         "mw_std": None, "usable": False, "outlier": bool(P["outlier"][i])})
            continue
        rows.append({"sv": sv_name(prn[i]), "start": fmt(info["start"]), "end": fmt(info["end"]), "n": info["n"],
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
    P1r, P2r = P.get("P1raw", P["P1"]), P.get("P2raw", P["P2"])        # cok yolluluk ham kodun ozelligi
    MP1 = P1r - (1 + k) * F1m + k * F2m
    MP2 = P2r - GAMMA * k * F1m + (GAMMA * k - 1) * F2m
    m = P["usable"] & valid
    for a in np.unique(P["arc"][m]):
        mm = m & (P["arc"] == a)
        MP1[mm] -= MP1[mm].mean()
        MP2[mm] -= MP2[mm].mean()
    MP1[~m] = np.nan
    MP2[~m] = np.nan
    rows = []
    for i in idx_ep:
        if P["isE"][i]:
            continue
        rows.append({"sv": sv_name(P["prn"][i]), "el": _r(el[i], 2), "az": _r(az[i], 2),
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


def _spp_series(P, sp3, nav, X0, lat0, h0, Rn, Xref, t_ep, step_s, cutoff, met, Pm, Tm_, RH, codes=None, grid=None):
    t, prn = P["t"], P["prn"]
    C1, C2 = codes if codes is not None else (P["P1"], P["P2"])
    eps = np.unique(t)
    if grid is None:
        grid = eps[np.abs(((eps - eps[0] + 1e-3) % step_s) - 1e-3) < 1e-3]
    if t_ep not in grid:
        grid = np.sort(np.r_[grid, t_ep])
    src = "broadcast" if nav else "sp3"
    out = []
    epoch_res = None
    # uydu konum ve saatleri: izgara epoklarindaki butun GPS kayitlari icin uydu basina toplu
    gk = np.round(t * 1000.0).astype(np.int64)
    gm = np.where(np.isin(gk, np.round(np.asarray(grid) * 1000.0).astype(np.int64)) & ~P["isE"])[0]
    P3g = A3 * C1[gm] + B3 * C2[gm]
    txg = t[gm] - P3g / C
    Xg = np.full((len(gm), 3), np.nan)
    cdg = np.full(len(gm), np.nan)
    for p_ in np.unique(prn[gm]):
        jj = np.where(prn[gm] == p_)[0]
        if nav:
            for j in jj:
                e = orbit.select_eph(nav["eph"], int(p_), txg[j])
                if e is None:
                    continue
                c0 = orbit.bc_clock(e, txg[j])
                pos, drel = orbit.bc_position(e, txg[j] - c0)
                Xg[j], cdg[j] = pos, (c0 + drel) * C        # P3 icin TGD gerekmez
        else:
            c0 = orbit.sp3_clock(sp3, int(p_), txg[jj])
            pos, vel = orbit.sp3_interp(sp3, int(p_), txg[jj] - c0)
            Xg[jj] = pos
            cdg[jj] = (c0 - 2 * np.einsum("ij,ij->i", pos, vel) / C ** 2) * C
    og = np.argsort(t[gm], kind="stable")
    tg_sorted = t[gm][og]
    for tt in grid:
        a0, a1 = np.searchsorted(tg_sorted, tt - 1e-3), np.searchsorted(tg_sorted, tt + 1e-3)
        sl = og[a0:a1]
        if len(sl) < 5:
            continue
        P3 = P3g[sl]
        Xs = Xg[sl]
        cd = cdg[sl]
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
        sky[sv_name(p_)] = {"az": _lst(az[m], 1), "el": _lst(el[m], 1)}
    sky_ep = [{"sv": sv_name(prn[i]), "az": _r(az[i], 1), "el": _r(el[i], 1), "used": bool(sel[i])} for i in idx_ep]
    # GF / MW (epok uydularinin epoktaki yaylari)
    gfmw = {}
    for i in idx_ep:
        a = P["arc"][i]
        if a < 0:
            continue
        m = np.where(P["arc"] == a)[0]
        m = m[_dec(len(m), 600)]
        gfmw[sv_name(prn[i])] = {"t": hrs(t[m]), "gf": _lst(P["GF"][m] - P["GF"][m[0]], 3),
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
        P1r = P.get("P1raw", P["P1"])
        raw = LAM1 * P["L1"][m] - P1r[m]
        cor = LAM1 * P["L1d"][m] - P1r[m]
        cj = {"sv": sv_name(prn[i]), "t": hrs(t[m]), "raw_km": _lst((raw - raw[0]) / 1e3, 4),
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


def _kin_track(sol, Rn, Xref, day0, rate):
    """Kinematik cozumun epok epok rotasi (NEU, referansa gore), hiz ve ozet."""
    tr = sol["track"]
    tt = tr["t"]
    X = tr["X"]
    d = (X - Xref) @ Rn.T                                   # (ne, 3): N, E, U
    dt = np.diff(tt)
    sp = np.full(len(tt), np.nan)
    if len(tt) > 1:
        v = np.hypot(np.diff(d[:, 0]), np.diff(d[:, 1])) / np.where(dt > 0, dt, np.nan)
        sp[1:] = v
    step = np.hypot(np.diff(d[:, 0]), np.diff(d[:, 1])) if len(tt) > 1 else np.zeros(0)
    moving = sp > 0.2
    sN = tr["sNEU"] * 1e3
    rows = []
    for k in range(len(tt)):
        rows.append({"t": fmt(tt[k]), "X": _r(X[k, 0], 4), "Y": _r(X[k, 1], 4), "Z": _r(X[k, 2], 4),
                     "N": _r(d[k, 0], 4), "E": _r(d[k, 1], 4), "U": _r(d[k, 2], 4),
                     "sN": _r(sN[k, 0], 1), "sE": _r(sN[k, 1], 1), "sU": _r(sN[k, 2], 1),
                     "n": int(tr["nsat"][k]), "pdop": _r(tr["pdop"][k], 2), "v": _r(sp[k], 3)})
    summ = {"n_epochs": int(len(tt)), "rate": rate, "first": fmt(tt[0]), "last": fmt(tt[-1]),
            "span_N": _r(np.ptp(d[:, 0]), 3), "span_E": _r(np.ptp(d[:, 1]), 3), "span_U": _r(np.ptp(d[:, 2]), 3),
            "path_m": _r(float(step.sum()), 2), "v_max": _r(np.nanmax(sp), 3) if np.isfinite(sp).any() else None,
            "v_mean_moving": _r(float(np.nanmean(sp[moving])), 3) if moving.any() else None,
            "moving_frac": _r(float(moving.sum()) / max(len(tt) - 1, 1), 3),
            "sig_med_mm": [_r(np.median(sN[:, 0]), 1), _r(np.median(sN[:, 1]), 1), _r(np.median(sN[:, 2]), 1)],
            "j_epoch": int(tr["j"])}
    return {"summary": summ, "rows": rows,
            "plot": {"t": _lst((tt - day0) / 3600.0, 5), "N": _lst(d[:, 0], 4), "E": _lst(d[:, 1], 4),
                     "U": _lst(d[:, 2], 4), "sN": _lst(sN[:, 0], 1), "sE": _lst(sN[:, 1], 1), "sU": _lst(sN[:, 2], 1),
                     "n": [int(x) for x in tr["nsat"]], "pdop": _lst(tr["pdop"], 2), "v": _lst(sp, 3)}}


def _codefilter_report(P, KF, sel, idx_ep, sol, sol_raw, mode, Rn, Xref, Xst, sp3, nav, X0, lat0, h0,
                       t_ep, spp, spp_step, cutoff, met, Pm, Tm_, RH, rate, day0, tau, thr):
    """Kod Kalman filtresinin raporu: gurultu azalmasi, kod cozumune ve faz cozumune etkisi, grafikler."""
    t, prn = P["t"], P["prn"]
    st = cf.noise_stats(P, KF, P["usable"], skip_sec=3.0 * tau)
    rows = []
    for r in st:
        rows.append({"sv": sv_name(r["prn"]), "n": r["n"], "raw1": _r(r["raw1"], 3), "filt1": _r(r["filt1"], 3),
                     "red1": _r(100 * (1 - r["filt1"] / r["raw1"]), 0) if r["raw1"] > 0 else None,
                     "raw2": _r(r["raw2"], 3), "filt2": _r(r["filt2"], 3),
                     "red2": _r(100 * (1 - r["filt2"] / r["raw2"]), 0) if r["raw2"] > 0 else None, "rej": r["rej"]})
    tot = None
    if st:
        w = np.array([r["n"] for r in st], float)
        av = lambda k: float(np.sqrt(np.average([r[k] ** 2 for r in st], weights=w)))
        tot = {"raw1": _r(av("raw1"), 3), "filt1": _r(av("filt1"), 3), "raw2": _r(av("raw2"), 3), "filt2": _r(av("filt2"), 3)}

    # kod cozumu (P3, GPS) ham ve suzulmus: faz cozumune gore fark
    eps = np.unique(t[sel])
    if mode == "kinematic":
        tr = sol["track"]
        ref_pos = {round(float(a), 3): tr["X"][k] for k, a in enumerate(tr["t"])}
        stepn = max(1, int(math.ceil(len(eps) / 300.0)))
        grid = eps[::stepn]
    else:
        ref_pos = None
        grid = None
    common = dict(sp3=sp3, nav=nav, X0=X0, lat0=lat0, h0=h0, Rn=Rn, Xref=Xref, t_ep=t_ep, step_s=spp_step,
                  cutoff=cutoff, met=met, Pm=Pm, Tm_=Tm_, RH=RH)
    s_raw = _spp_series(P, codes=(P["P1raw"], P["P2raw"]), grid=grid, **common)
    s_flt = spp if grid is None else _spp_series(P, grid=grid, **common)

    def diffs(series):
        out = {}
        for o in series.get("series") or []:
            if ref_pos is None:
                X = Xst
            else:
                X = ref_pos.get(round(o["t"], 3))
                if X is None:
                    continue
            d = Rn @ (np.asarray(X) - Xref)
            out[round(o["t"], 3)] = np.array([o["N"] - d[0], o["E"] - d[1], o["U"] - d[2]])
        return out
    dr, dfl = diffs(s_raw), diffs(s_flt)
    ks = sorted(set(dr) & set(dfl))
    spp_cmp = None
    n_bad_cmp = 0
    if ks:
        A = np.array([dr[k] for k in ks])
        B = np.array([dfl[k] for k in ks])
        good = (np.linalg.norm(A, axis=1) < 50.0) & (np.linalg.norm(B, axis=1) < 50.0)   # bozuk epoklar (SPP ile ayni olcut)
        n_bad_cmp = int((~good).sum())
        A, B = A[good], B[good]
        ks = [k for k, g in zip(ks, good) if g]
    if ks:
        rms = lambda X: {"N": _r(np.sqrt(np.mean(X[:, 0] ** 2)), 3), "E": _r(np.sqrt(np.mean(X[:, 1] ** 2)), 3),
                         "U": _r(np.sqrt(np.mean(X[:, 2] ** 2)), 3),
                         "H": _r(np.sqrt(np.mean(X[:, 0] ** 2 + X[:, 1] ** 2)), 3),
                         "D3": _r(np.sqrt(np.mean(np.sum(X ** 2, 1))), 3)}
        spp_cmp = {"n": len(ks), "n_bad": n_bad_cmp, "raw": rms(A), "filt": rms(B),
                   "plot": {"t": _lst((np.array(ks) - day0) / 3600.0, 5),
                            "H_raw": _lst(np.hypot(A[:, 0], A[:, 1]), 3), "H_filt": _lst(np.hypot(B[:, 0], B[:, 1]), 3),
                            "U_raw": _lst(A[:, 2], 3), "U_filt": _lst(B[:, 2], 3)}}

    # faz cozumune etkisi
    phase_eff = None
    if sol_raw is not None:
        if mode == "kinematic":
            d = (sol["track"]["X"] - sol_raw["track"]["X"]) @ Rn.T * 1e3
            dd = {"rms_neu_mm": _lst(np.sqrt(np.mean(d ** 2, 0)), 2), "max_neu_mm": _lst(np.max(np.abs(d), 0), 1)}
        else:
            d = Rn @ (sol["X"] - sol_raw["X"]) * 1e3
            dd = {"neu_mm": _lst(d, 2)}
        phase_eff = {**dd, "rms_code_raw": _r(sol_raw["rms_code"], 3), "rms_code_filt": _r(sol["rms_code"], 3),
                     "rms_phase_raw_mm": _r(sol_raw["rms_phase"] * 1e3, 2), "rms_phase_filt_mm": _r(sol["rms_phase"] * 1e3, 2),
                     "s0_raw": _r(sol_raw["s0"], 3), "s0_filt": _r(sol["s0"], 3)}

    # ornek uydu: epoktaki uydulardan en uzun yay
    plot = None
    best = None
    for i in idx_ep:
        a = P["arc"][i]
        if a < 0:
            continue
        m = np.where(P["arc"] == a)[0]
        if best is None or len(m) > len(best[1]):
            best = (i, m)
    if best is not None:
        i, m = best
        m = m[np.argsort(t[m])]
        mm = m[_dec(len(m), 2500)]
        c1r = KF["cmc1_raw"][m]
        c1f = KF["cmc1_f"][m]
        off = np.nanmean(c1r)
        plot = {"sv": sv_name(prn[i]), "t": _lst((t[mm] - day0) / 3600.0, 5),
                "cmc1_raw": _lst(KF["cmc1_raw"][mm] - off, 3), "cmc1_f": _lst(KF["cmc1_f"][mm] - off, 3),
                "cmc2_raw": _lst(KF["cmc2_raw"][mm] - np.nanmean(KF["cmc2_raw"][m]), 3),
                "cmc2_f": _lst(KF["cmc2_f"][mm] - np.nanmean(KF["cmc2_raw"][m]), 3),
                "K1": _lst(KF["K1"][mm], 4), "s1": _lst(KF["s1"][mm], 3),
                "rej_t": _lst((t[m][KF["rej1"][m] | KF["rej2"][m]] - day0) / 3600.0, 5)}
    used = P["arc"] >= 0
    return {"enabled": True, "tau": tau, "thr": thr, "rows": rows, "total": tot,
            "n_rej": int(KF["rej1"][used].sum() + KF["rej2"][used].sum()), "n_restart": int(KF["n_restart"]),
            "n_filtered": int(used.sum()), "k_ss": _r(min(float(P["interval"]) / tau, 1.0), 3),
            "spp_cmp": spp_cmp, "phase_effect": phase_eff, "plot": plot}


def _grid_select(t, rate, ref0, interval, interp_ok):
    """Islem epoklari: ref0 + k*rate dugumleri. Gozlem epoklari dugumlere oturmuyorsa (kaydirilmis ya da
    ondalikli zaman etiketleri) her dugume en yakin gozlem epogu secilir. Doner: (maske, kaydirildi_mi)."""
    ep = np.unique(t)
    exact = np.abs(((ep - ref0 + 1e-3) % rate) - 1e-3) < 1e-3
    k0 = math.ceil((ep[0] - ref0) / rate - 1e-9)
    k1 = math.floor((ep[-1] - ref0) / rate + 1e-9)
    n_nodes = max(k1 - k0 + 1, 1)
    if exact.sum() >= 0.5 * min(n_nodes, len(ep)):
        return np.isin(t, ep[exact]), False
    tol = 0.5 * rate if interp_ok else max(0.5 * interval, 1.0)
    nodes = ref0 + rate * np.arange(k0 - 1, k1 + 2)
    j = np.clip(np.searchsorted(ep, nodes), 1, len(ep) - 1)
    left, right = ep[j - 1], ep[j]
    best = np.where(np.abs(nodes - left) <= np.abs(right - nodes), left, right)
    ok = np.abs(best - nodes) <= tol + 1e-6
    chosen = np.unique(best[ok])
    return np.isin(t, chosen), True


def _clk_sp3_check(clk, sp3, t0, t1):
    """CLK ve SP3 uydu saatlerinin SP3 epoklarinda tutarliligi [m]: epok ortalamasi cikarilmis uydu farklarinin
    medyan RMS'i. Ayni analiz merkezinin urunlerinde mm, farkli merkezlerde dm-m duzeyindedir."""
    ts = sp3["t"][(sp3["t"] >= t0 - 1) & (sp3["t"] <= t1 + 1)]
    if len(ts) < 2:
        return None
    rows = []
    for p in sp3["prns"]:
        a = orbit.sp3_clock(clk, p, ts)
        b = sp3["clk"][np.isin(sp3["t"], ts), sp3["prns"].index(p)]
        rows.append((a - b) * C)
    D = np.array(rows)                       # (uydu, epok)
    if not np.isfinite(D).any():
        return None
    D = D - np.nanmedian(D, axis=0)[None, :]
    rms = np.sqrt(np.nanmean(D ** 2, axis=1))
    rms = rms[np.isfinite(rms)]
    return float(np.median(rms)) if len(rms) else None


def _motion_test(P, sp3, nav, X0, lat0, h0, Rn, cutoff, met, Pm, Tm_, RH):
    """Kod cozumune (GPS, en cok ~120 epok) gore yatay dagilim: medyandan uzakliklarin %90'i [m]."""
    ep = np.unique(P["t"][~P["isE"]])
    if len(ep) < 5:
        return None
    grid = ep[np.linspace(0, len(ep) - 1, min(120, len(ep))).astype(int)]
    s = _spp_series(P, sp3, nav, X0, lat0, h0, Rn, X0, grid[0], 300.0, cutoff, met, Pm, Tm_, RH, grid=grid)
    ser = s.get("series") or []
    if len(ser) < 5:
        return None
    N = np.array([o["N"] for o in ser])
    E = np.array([o["E"] for o in ser])
    ok = np.array([o["rms"] for o in ser]) < 20.0
    if ok.sum() < 5:
        return None
    hor = np.hypot(N[ok] - np.median(N[ok]), E[ok] - np.median(E[ok]))
    return float(np.percentile(hor, 90))


def _clean_residuals(sol, sel, P, t, el, sig_phase, sig_code, removed, clean):
    """Artiga dayali ayiklama. Faz artiklarinda basamak bulunan yaylar bolunur (kacmis faz kopmasi),
    tek tek buyuk faz ya da kod artiklari cikarilir. Degisiklik olduysa True."""
    idx = np.where(sel)[0]
    se = np.sin(np.radians(el[idx]))
    vph = sol["v_phase"] * se                     # yukseklige gore olceklenmis
    vcd = sol["v_code"] * se
    sph = 1.4826 * float(np.median(np.abs(vph))) + 1e-4
    scd = 1.4826 * float(np.median(np.abs(vcd))) + 1e-3
    changed = False
    # 1) basamak: her yayda pencere ortalamalari farki
    arc = P["arc"]
    next_id = int(arc.max()) + 1
    cands = []
    for a in np.unique(arc[idx]):
        m = np.where(arc[idx] == a)[0]
        if len(m) < 8:
            continue
        m = m[np.argsort(t[idx][m])]
        r = vph[m]
        tt = t[idx][m]
        dtm = float(np.median(np.diff(tt))) if len(tt) > 1 else 30.0
        W = int(np.clip(round(120.0 / max(dtm, 1e-3)), 3, max(3, len(m) // 4)))   # ~2 dk pencere
        if len(m) < 2 * W:
            continue
        c = np.r_[0.0, np.cumsum(r)]
        k = np.arange(W, len(m) - W + 1)
        stepv = (c[k + W] - c[k]) / W - (c[k] - c[k - W]) / W
        thr = max(0.025, 5.0 * sph * math.sqrt(2.0 / W))
        j = int(np.argmax(np.abs(stepv)))
        if abs(stepv[j]) > thr:
            cands.append((abs(stepv[j]), a, t[idx][m][k[j]]))
    if cands:
        # yalniz en belirgin basamaklar: bozuk yaylar oteki uydularin artiklarini da kaydirir
        big = max(c[0] for c in cands)
        for sz, a, t_cut in cands:
            if sz >= 0.5 * big:
                allm = np.where((arc == a) & (t >= t_cut - 1e-6))[0]
                arc[allm] = next_id
                next_id += 1
                clean["split"] += 1
                changed = True
        return True
    # 2) tek tek aykirilar (en buyukler, en cok kayitlarin %2'si)
    lim = max(1, int(0.02 * len(idx)))
    bad_ph = np.where(np.abs(vph) > max(0.03, 5.0 * sph))[0]
    bad_cd = np.where(np.abs(vcd) > max(5.0, 5.0 * scd))[0]
    if len(bad_ph):
        bad_ph = bad_ph[np.argsort(-np.abs(vph[bad_ph]))][:lim]
        removed[idx[bad_ph]] = True
        clean["removed_phase"] += len(bad_ph)
        changed = True
    if len(bad_cd):
        bad_cd = bad_cd[np.argsort(-np.abs(vcd[bad_cd]))][:lim]
        removed[idx[bad_cd]] = True
        clean["removed_code"] += len(bad_cd)
        changed = True
    return changed
