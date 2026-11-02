"""En kucuk kareler: faz + kod cozumu (statik, secilen epokta kinematik ya da tum epoklarda kinematik), DOP, kod SPP."""
import math

import numpy as np

from .const import C, OMEGA_E


def dop(u, Rneu):
    G = np.column_stack([-u, np.ones(len(u))])
    Q = np.linalg.inv(G.T @ G)
    Qn = Rneu @ Q[:3, :3] @ Rneu.T
    return {"GDOP": float(math.sqrt(np.trace(Q))), "PDOP": float(math.sqrt(np.trace(Q[:3, :3]))),
            "HDOP": float(math.sqrt(Qn[0, 0] + Qn[1, 1])), "VDOP": float(math.sqrt(Qn[2, 2])),
            "TDOP": float(math.sqrt(Q[3, 3]))}


def _sparse_global(n, ie, iy, iz, mw, u, isE, static_rows, gX, gz, gi, gb):
    """Global bilinmeyenlerin (statik konum, ZWD, ISB, belirsizlik) seyrek tasarim ogeleri.

    Satirlar: 0..n-1 faz, n..2n-1 kod. Doner: (satir, sutun, deger) dizileri.
    """
    r_all = np.arange(2 * n)
    rec = np.r_[np.arange(n), np.arange(n)]
    R, Cc, V = [], [], []
    if gX is not None:
        rr = r_all[np.r_[static_rows, static_rows]]
        kk = rec[rr]
        for a in range(3):
            R.append(rr); Cc.append(np.full(len(rr), gX + a)); V.append(-u[kk, a])
    R.append(r_all); Cc.append(gz + iz[rec]); V.append(mw[rec])
    if gi is not None:
        rr = r_all[np.r_[isE, isE]]
        R.append(rr); Cc.append(np.full(len(rr), gi)); V.append(np.ones(len(rr)))
    R.append(np.arange(n)); Cc.append(gb + iy); V.append(np.ones(n))
    return np.concatenate(R), np.concatenate(Cc).astype(int), np.concatenate(V)


def phase_solution(d, phase, code, corr_phase, corr_code, X0, kin_epoch, Rneu,
                   sig_phase=0.01, sig_code=1.0, zwd_block=7200.0, iters=5, mode="static"):
    """Faz + kod EKK cozumu (epok parametreleri epok epok on elimine edilir).

    d: dict, kayit dizileri (t, arc, Xs(n,3), el, cdts, Tmodel, mw, istege bagli isE)
    phase, code: olculer [m]; corr_*: model uzakliga eklenecek duzeltmeler [m]
    mode="static": bilinmeyenler statik konum (3), kin_epoch'ta ayri konum (3), her epokta
        saat, zwd_block saniyelik islak troposfer, Galileo varsa GPS-Galileo sistemler arasi
        saat farki (ISB), her yay icin ondalikli belirsizlik.
    mode="kinematic": statik konum yerine her epokta ayri konum (3).
    """
    t = d["t"]
    n = len(t)
    kinematic = mode == "kinematic"
    ep, ie = np.unique(t, return_inverse=True)
    yy, iy = np.unique(d["arc"], return_inverse=True)
    zb = np.floor((t - t.min()) / zwd_block).astype(int)
    zz, iz = np.unique(zb, return_inverse=True)
    ne, ny, nz = len(ep), len(yy), len(zz)
    isE = np.asarray(d.get("isE", np.zeros(n, bool)), bool)
    has_isb = bool(isE.any() and (~isE).any())

    # epok konumu tasiyan epoklar
    haspos = np.full(ne, kinematic)
    j_kin = -1
    has_kin = False
    if not kinematic and kin_epoch is not None:
        km = np.abs(t - kin_epoch) < 1e-3
        if km.sum() >= 4:
            has_kin = True
            j_kin = int(ie[np.argmax(km)])
            haspos[j_kin] = True
    j_sel = int(np.argmin(np.abs(ep - kin_epoch))) if kin_epoch is not None else ne // 2
    static_rows = ~haspos[ie]

    gX = None if kinematic else 0
    gz = 0 if kinematic else 3
    gi = gz + nz if has_isb else None
    gb = gz + nz + (1 if has_isb else 0)
    ng = gb + ny

    Xs = d["Xs"]
    s2 = np.sin(np.radians(d["el"])) ** 2
    cdts = d["cdts"]
    Tm = d["Tmodel"]
    mw = d["mw"]
    X = np.array(X0, float).copy()
    Xe = np.tile(X, (ne, 1))
    rho0 = np.linalg.norm(Xs - X, axis=1)
    k0 = code - rho0 + cdts - Tm - corr_code
    o = np.argsort(ie, kind="stable")
    cut = np.r_[0, np.nonzero(np.diff(ie[o]))[0] + 1, n]
    clk = np.array([np.median(k0[o[a:b]]) for a, b in zip(cut[:-1], cut[1:])])
    zwd = np.zeros(nz)
    isb = float(np.median(k0[isE] - clk[ie][isE])) if has_isb else 0.0
    B = np.zeros(ny)
    p = np.concatenate([s2 / sig_phase ** 2, s2 / sig_code ** 2])
    erow = np.r_[ie, ie]
    hist = []
    for it in range(iters):
        Xa = np.where(haspos[ie][:, None], Xe[ie], X)
        D = Xs - Xa
        rho = np.linalg.norm(D, axis=1)
        u = D / rho[:, None]
        com = rho + clk[ie] - cdts + Tm + zwd[iz] * mw + isb * isE
        l = np.concatenate([phase - (com + B[iy] + corr_phase), code - (com + corr_code)])

        # epok tasarimi: [-u * haspos, 1]
        a = np.ones((2 * n, 4))
        a[:, :3] = np.r_[-u, -u] * haspos[erow][:, None]
        pa = a * p[:, None]
        Nee = np.zeros((ne, 4, 4))
        np.add.at(Nee, erow, pa[:, :, None] * a[:, None, :])
        Nee[~haspos, 0, 0] = Nee[~haspos, 1, 1] = Nee[~haspos, 2, 2] = 1.0
        be = np.zeros((ne, 4))
        np.add.at(be, erow, pa * l[:, None])

        # global tasarim (seyrek)
        gr, gc, gv = _sparse_global(n, ie, iy, iz, mw, u, isE, static_rows, gX, gz, gi, gb)
        Neg = np.zeros((ne, ng, 4))
        np.add.at(Neg, (erow[gr], gc), pa[gr] * gv[:, None])
        bg = np.zeros(ng)
        np.add.at(bg, gc, p[gr] * gv * l[gr])
        # Ngg: ayni satirdaki oge ciftleri
        Ngg = np.zeros((ng, ng))
        order = np.argsort(gr, kind="stable")
        gr_s, gc_s, gv_s = gr[order], gc[order], gv[order]
        starts = np.r_[0, np.nonzero(np.diff(gr_s))[0] + 1]
        counts = np.diff(np.r_[starts, len(gr_s)])
        kmax = int(counts.max())
        for k1 in range(kmax):
            m1 = counts > k1
            i1 = starts[m1] + k1
            for k2 in range(kmax):
                m2 = counts[m1] > k2
                i1b = i1[m2]
                i2 = starts[m1][m2] + k2
                np.add.at(Ngg, (gc_s[i1b], gc_s[i2]), p[gr_s[i1b]] * gv_s[i1b] * gv_s[i2])

        Ninv = np.linalg.inv(Nee)
        NegT = np.transpose(Neg, (0, 2, 1))                      # (ne, 4, ng)
        H = np.einsum("eab,ebj->eaj", Ninv, NegT, optimize=True)
        Nred = Ngg - np.einsum("eai,eaj->ij", NegT, H, optimize=True)
        bred = bg - np.einsum("eai,ea->i", H, be, optimize=True)
        try:
            dg = np.linalg.solve(Nred, bred)
        except np.linalg.LinAlgError:
            dg = np.linalg.lstsq(Nred, bred, rcond=None)[0]
        de = np.einsum("eab,eb->ea", Ninv, be - np.einsum("eaj,j->ea", NegT, dg, optimize=True))
        de[~haspos, :3] = 0.0

        if gX is not None:
            X += dg[gX:gX + 3]
        Xe[haspos] += de[haspos, :3]
        clk += de[:, 3]
        zwd += dg[gz:gz + nz]
        if has_isb:
            isb += dg[gi]
        B += dg[gb:]
        if kinematic:
            step = float(np.max(np.linalg.norm(de[:, :3], axis=1)))
        else:
            step = float(np.linalg.norm(dg[0:3]))
        hist.append(step)
        if it >= 1 and step < 1e-5:
            break

    # artiklar (son iterasyonun dogrusal modeli)
    vg = np.zeros(2 * n)
    np.add.at(vg, gr, gv * dg[gc])
    v = np.einsum("ra,ra->r", a, de[erow]) + vg - l
    n_ep_par = int(ne + 3 * haspos.sum())
    dof = 2 * n - (n_ep_par + ng)
    s0 = float(math.sqrt((v ** 2 * p).sum() / max(dof, 1)))
    try:
        Qgg = np.linalg.inv(Nred)
    except np.linalg.LinAlgError:
        Qgg = np.linalg.pinv(Nred)
    HQ = np.einsum("eai,ij->eaj", H, Qgg, optimize=True)
    Qee = Ninv + np.einsum("eaj,ebj->eab", HQ, H, optimize=True)
    sg = s0 * np.sqrt(np.abs(np.diag(Qgg)))
    se = s0 * np.sqrt(np.abs(np.einsum("eaa->ea", Qee)))

    out = {"epochs": ep, "clock": clk, "sclock": se[:, 3], "zwd_blocks": zz * zwd_block + t.min(),
           "zwd": zwd, "szwd": sg[gz:gz + nz], "v_phase": v[:n], "v_code": v[n:], "s0": s0, "dof": dof,
           "n_obs": 2 * n, "n_par": n_ep_par + ng, "n_epochs": ne, "n_arcs": ny, "n_zwd": nz, "iter": hist,
           "rms_phase": float(np.sqrt((v[:n] ** 2).mean())), "rms_code": float(np.sqrt((v[n:] ** 2).mean())),
           "has_kin": has_kin, "iz": iz, "ie": ie, "mode": mode,
           "isb": isb if has_isb else None, "sisb": float(sg[gi]) if has_isb else None}
    if kinematic:
        j = j_sel
        out["X"] = Xe[j].copy()
        out["sX"] = se[j, :3]
        out["Xk"] = None
        out["sXk"] = None
        Qp = s0 ** 2 * Qee[:, :3, :3]
        Qn = np.einsum("ab,ebc,dc->ead", Rneu, Qp, Rneu)
        nsat = np.bincount(ie, minlength=ne)
        pdop = np.full(ne, np.nan)
        cnt = np.bincount(ie, weights=None, minlength=ne)
        for e in range(ne):
            m = o[cut[e]:cut[e + 1]]
            if len(m) >= 4:
                try:
                    pdop[e] = dop(u[m], Rneu)["PDOP"]
                except np.linalg.LinAlgError:
                    pass
        out["track"] = {"t": ep, "X": Xe, "sNEU": np.sqrt(np.abs(np.einsum("eaa->ea", Qn))), "nsat": nsat,
                        "pdop": pdop, "j": j}
        m = o[cut[j]:cut[j + 1]]
        out["dop"] = dop(u[m], Rneu) if len(m) >= 4 else {}
        out["n_sat_kin"] = int(len(m))
        out["clock_kin"] = float(clk[j])
        out["sclock_kin"] = float(se[j, 3])
        out["zwd_kin"] = float(zwd[iz[m[0]]]) if len(m) else 0.0
        out["j_sel"] = j
    else:
        out["X"] = X
        out["sX"] = sg[0:3]
        out["Xk"] = Xe[j_kin].copy() if has_kin else None
        out["sXk"] = se[j_kin, :3] if has_kin else None
        if has_kin:
            m = o[cut[j_kin]:cut[j_kin + 1]]
            Dk = Xs[m] - Xe[j_kin]
            out["dop"] = dop(Dk / np.linalg.norm(Dk, axis=1)[:, None], Rneu)
            out["n_sat_kin"] = int(len(m))
            out["clock_kin"] = float(clk[j_kin])
            out["sclock_kin"] = float(se[j_kin, 3])
            out["zwd_kin"] = float(zwd[iz[m[0]]])
    return out


def spp_epoch(Ps, Xs, cdts, Tm, X0=None, iters=8, el=None, w_el=True):
    """Tek epok kod konumu: Ps (m,), Xs (m,3) Dunya donusu uygulanmamis, cdts [m], Tm [m]."""
    m = len(Ps)
    if m < 4:
        return None
    X = np.zeros(3) if X0 is None else np.array(X0, float)
    cdt = 0.0
    for _ in range(iters):
        rr = np.linalg.norm(Xs - X, axis=1)
        th = OMEGA_E * rr / C
        Xr = np.column_stack([np.cos(th) * Xs[:, 0] + np.sin(th) * Xs[:, 1],
                              -np.sin(th) * Xs[:, 0] + np.cos(th) * Xs[:, 1], Xs[:, 2]])
        D = Xr - X
        rho = np.linalg.norm(D, axis=1)
        u = D / rho[:, None]
        l = Ps - (rho + cdt - cdts + (Tm if np.linalg.norm(X) > 6e6 else 0.0))
        A = np.column_stack([-u, np.ones(m)])
        dx = np.linalg.lstsq(A, l, rcond=None)[0]
        X += dx[:3]
        cdt += dx[3]
        if np.linalg.norm(dx[:3]) < 1e-4:
            break
    res = l - A @ dx
    return {"X": X, "cdt": cdt, "rms": float(np.sqrt((res ** 2).mean())), "n": m}


# ======================================================================
# Kalman filtresi (genisletilmis, EKF) ve RTS yumusatici
# ======================================================================
def _q_neu(Rneu, qh, qv):
    """NEU'da diag(qh, qh, qv) -> ECEF (3x3)."""
    return Rneu.T @ np.diag([qh, qh, qv]) @ Rneu


def kalman_solution(d, phase, code, corr_phase, corr_code, X0, kin_epoch, Rneu,
                    sig_phase=0.01, sig_code=1.0, dyn="pv", q_h=0.5, q_v=0.01, q_zwd=1e-8,
                    smooth=True, iters=2, thr=5.0, sig_amb0=10.0, max_mem_mb=300.0):
    """Faz + kod olculerini epok epok isleyen genisletilmis Kalman filtresi ve RTS yumusatici.

    Durum vektoru x = [X (3), V (3, yalniz dyn='pv'), c*dt_r, ZWD, ISB (Galileo varsa), B_1..B_ny]
    Dinamik model (konum):
        'static': X sabit (Q = 0)
        'rw'    : rastgele yuruyus, Q = diag(q_h, q_h, q_v) * dt               [q: m^2/s]
        'pv'    : sabit hiz, beyaz ivme gurultusu; Q_XX = S dt^3/3, Q_XV = S dt^2/2,
                  Q_VV = S dt, S = diag(q_h, q_h, q_v)                         [q: m^2/s^3]
    Saat her epokta yeniden (beyaz gurultu), ZWD rastgele yuruyus (q_zwd m^2/s), ISB ve yay
    belirsizlikleri sabit. Yay basladiginda (ya da fazda kopma bulundugunda) belirsizlik
    L3 - P3 farkindan ve buyuk varyansla yeniden baslatilir.
    Olcum guncellemesi: yinelemeli EKF (iters), agirlik sin^2(E)/sigma^2.
    Normalize yenilik |v|/sqrt(S) > thr: kod satiri atilir, faz satirinda yay sifirlanir.
    """
    t = d["t"]
    n = len(t)
    ep, ie = np.unique(t, return_inverse=True)
    yy, iy = np.unique(d["arc"], return_inverse=True)
    ne, ny = len(ep), len(yy)
    isE = np.asarray(d.get("isE", np.zeros(n, bool)), bool)
    has_isb = bool(isE.any() and (~isE).any())
    pv = dyn == "pv"
    nv = 3 if pv else 0
    iX, iV = 0, 3
    ic = 3 + nv
    iz = ic + 1
    ii = iz + 1 if has_isb else None
    ib = iz + 1 + (1 if has_isb else 0)
    ns = ib + ny
    Xs = d["Xs"]
    s2 = np.sin(np.radians(d["el"])) ** 2
    cdts, Tm, mw = d["cdts"], d["Tmodel"], d["mw"]
    o = np.argsort(ie, kind="stable")
    cut = np.r_[0, np.nonzero(np.diff(ie[o]))[0] + 1, n]
    first_ep = np.full(ny, ne, int)
    np.minimum.at(first_ep, iy, ie)
    do_smooth = bool(smooth) and (2.0 * ne * ns * ns * 8 / 1e6 <= max_mem_mb)

    # baslangic durumu
    x = np.zeros(ns)
    P = np.zeros((ns, ns))
    x[iX:iX + 3] = X0
    P[iX:iX + 3, iX:iX + 3] = np.eye(3) * 100.0 ** 2
    if pv:
        P[iV:iV + 3, iV:iV + 3] = np.eye(3) * 2.0 ** 2
    P[ic, ic] = 1e3 ** 2
    P[iz, iz] = 0.5 ** 2
    if has_isb:
        P[ii, ii] = 1e3 ** 2
    P[ib:, ib:] = np.eye(ny) * 1e4

    Qs = _q_neu(Rneu, q_h, q_v)
    keep = (do_smooth, ne, ns)
    if do_smooth:
        Xp = np.zeros((ne, ns)); Pp = np.zeros((ne, ns, ns))
        Xf = np.zeros((ne, ns)); Pf = np.zeros((ne, ns, ns))
        zero_mask = np.zeros((ne, ns), bool)
        dts_ep = np.zeros(ne)
    Xfw = np.zeros((ne, 3)); Pfw = np.zeros((ne, 3, 3))
    clk_f = np.zeros(ne); sclk_f = np.zeros(ne); zwd_f = np.zeros(ne)
    n_rej = 0
    n_reset = 0
    rej = np.zeros(2 * n, bool)
    reset_log = []
    t_prev = None

    def amb_init(m):
        # L3 - P3 (duzeltmelerle) ~ B + kod gurultusu
        return (phase[m] - corr_phase[m]) - (code[m] - corr_code[m])

    for e in range(ne):
        m = o[cut[e]:cut[e + 1]]
        dt = 0.0 if t_prev is None else float(ep[e] - t_prev)
        t_prev = ep[e]
        # ---------------- zaman guncellemesi (tahmin)
        zmask = np.zeros(ns, bool)
        if e > 0:
            Phi = np.eye(ns)
            Q = np.zeros((ns, ns))
            if pv:
                Phi[iX:iX + 3, iV:iV + 3] = np.eye(3) * dt
                Q[iX:iX + 3, iX:iX + 3] = Qs * dt ** 3 / 3.0
                Q[iX:iX + 3, iV:iV + 3] = Qs * dt ** 2 / 2.0
                Q[iV:iV + 3, iX:iX + 3] = Qs * dt ** 2 / 2.0
                Q[iV:iV + 3, iV:iV + 3] = Qs * dt
            elif dyn == "rw":
                Q[iX:iX + 3, iX:iX + 3] = Qs * dt
            Q[iz, iz] = q_zwd * dt
            x = Phi @ x
            P = Phi @ P @ Phi.T + Q
        # saat: beyaz gurultu (her epokta yeniden), yeni yaylar: yeniden baslat
        zmask[ic] = True
        new_arcs = np.unique(iy[m][first_ep[iy[m]] == e])
        zmask[ib + new_arcs] = True
        # saat icin kaba deger: kod - geometri
        rho = np.linalg.norm(Xs[m] - x[iX:iX + 3], axis=1)
        isb_m = x[ii] * isE[m] if has_isb else 0.0
        kc = code[m] - (rho - cdts[m] + Tm[m] + x[iz] * mw[m] + isb_m + corr_code[m])
        if e == 0 and has_isb and isE[m].any() and (~isE[m]).any():
            x[ii] = float(np.median(kc[isE[m]]) - np.median(kc[~isE[m]]))
            kc = kc - x[ii] * isE[m]
        x[ic] = float(np.median(kc))
        for a_ in new_arcs:
            k_ = m[iy[m] == a_]
            x[ib + a_] = float(amb_init(k_)[0])
        for j in np.where(zmask)[0]:
            P[j, :] = 0.0
            P[:, j] = 0.0
        P[ic, ic] = 1e3 ** 2 if e == 0 else 300.0 ** 2
        P[ib + new_arcs, ib + new_arcs] = sig_amb0 ** 2

        # ---------------- olcum guncellemesi (yinelemeli EKF), yenilik testi
        for attempt in range(3):
            if do_smooth:
                Xp[e] = x; Pp[e] = P; zero_mask[e] = zmask; dts_ep[e] = dt
            xm = x.copy()
            Pm = P.copy()
            xi = xm.copy()
            mm = len(m)
            ok = ~rej[np.r_[m, n + m]]
            for it in range(iters):
                D = Xs[m] - xi[iX:iX + 3]
                rho = np.linalg.norm(D, axis=1)
                u = D / rho[:, None]
                com = rho + xi[ic] - cdts[m] + Tm[m] + xi[iz] * mw[m] + (xi[ii] * isE[m] if has_isb else 0.0)
                h = np.r_[com + xi[ib + iy[m]] + corr_phase[m], com + corr_code[m]]
                z = np.r_[phase[m], code[m]]
                H = np.zeros((2 * mm, ns))
                H[:, iX:iX + 3] = np.r_[-u, -u]
                H[:, ic] = 1.0
                H[:, iz] = np.r_[mw[m], mw[m]]
                if has_isb:
                    H[:, ii] = np.r_[isE[m], isE[m]].astype(float)
                H[np.arange(mm), ib + iy[m]] = 1.0
                Rv = np.r_[sig_phase ** 2 / s2[m], sig_code ** 2 / s2[m]]
                v = z - h - H @ (xm - xi)
                Hk, vk, Rk = H[ok], v[ok], Rv[ok]
                PHt = Pm @ Hk.T
                Sm = Hk @ PHt + np.diag(Rk)
                K = np.linalg.solve(Sm, PHt.T).T
                xi = xm + K @ vk
            # yenilik testi (tahmin edilen duruma gore)
            wn = np.abs(vk) / np.sqrt(np.diag(Sm))
            bad = np.where(ok)[0][wn > thr]
            if attempt < 2 and len(bad):
                bad_code = bad[bad >= mm]
                bad_ph = bad[bad < mm]
                if len(bad_code):
                    rej[n + m[bad_code - mm]] = True
                    n_rej += len(bad_code)
                for r in bad_ph:
                    a_ = iy[m[r]]
                    j = ib + a_
                    x[j] = float(amb_init(np.array([m[r]]))[0])
                    P[j, :] = 0.0
                    P[:, j] = 0.0
                    P[j, j] = sig_amb0 ** 2
                    zmask[j] = True
                    n_reset += 1
                    reset_log.append((float(ep[e]), int(a_)))
                continue
            break
        IKH = np.eye(ns) - K @ Hk
        x = xi
        P = IKH @ Pm @ IKH.T + K @ np.diag(Rk) @ K.T
        P = 0.5 * (P + P.T)
        if do_smooth:
            Xf[e] = x; Pf[e] = P
        Xfw[e] = x[iX:iX + 3]
        Pfw[e] = P[iX:iX + 3, iX:iX + 3]
        clk_f[e] = x[ic]; sclk_f[e] = math.sqrt(max(P[ic, ic], 0)); zwd_f[e] = x[iz]

    # ---------------- RTS yumusatici
    if do_smooth:
        Xsm = Xf.copy(); Psm = Pf.copy()
        for e in range(ne - 2, -1, -1):
            dt = dts_ep[e + 1]
            Phi = np.eye(ns)
            if pv:
                Phi[iX:iX + 3, iV:iV + 3] = np.eye(3) * dt
            Phi[zero_mask[e + 1], :] = 0.0
            A = Phi @ Pf[e]
            try:
                Cm = np.linalg.solve(Pp[e + 1], A).T
            except np.linalg.LinAlgError:
                Cm = (np.linalg.pinv(Pp[e + 1]) @ A).T
            Xsm[e] = Xf[e] + Cm @ (Xsm[e + 1] - Xp[e + 1])
            Psm[e] = Pf[e] + Cm @ (Psm[e + 1] - Pp[e + 1]) @ Cm.T
        Xe = Xsm[:, iX:iX + 3]
        Pe = Psm[:, iX:iX + 3, iX:iX + 3]
        clk = Xsm[:, ic]
        sclk = np.sqrt(np.abs(Psm[:, ic, ic]))
        zwd = Xsm[:, iz]
        xfin = Xsm
    else:
        Xe, Pe, clk, sclk, zwd = Xfw, Pfw, clk_f, sclk_f, zwd_f
        xfin = None

    # ---------------- artiklar (yumusatilmis ya da ileri durumla)
    v_ph = np.zeros(n)
    v_cd = np.zeros(n)
    for e in range(ne):
        m = o[cut[e]:cut[e + 1]]
        D = Xs[m] - Xe[e]
        rho = np.linalg.norm(D, axis=1)
        isb_e = (xfin[e, ii] if xfin is not None else x[ii]) if has_isb else 0.0
        com = rho + clk[e] - cdts[m] + Tm[m] + zwd[e] * mw[m] + isb_e * isE[m]
        B_e = xfin[e, ib + iy[m]] if xfin is not None else x[ib + iy[m]]
        v_ph[m] = (com + B_e + corr_phase[m]) - phase[m]
        v_cd[m] = (com + corr_code[m]) - code[m]
    w_ph = v_ph ** 2 * s2 / sig_phase ** 2
    w_cd = (v_cd ** 2 * s2 / sig_code ** 2)[~rej[n:]]
    nobs = 2 * n - int(rej.sum())
    s0 = float(math.sqrt((w_ph.sum() + w_cd.sum()) / max(nobs, 1)))
    isb = (float(xfin[0, ii]) if xfin is not None else float(x[ii])) if has_isb else None
    sisb = (float(math.sqrt(abs(Psm[0, ii, ii]))) if xfin is not None else float(math.sqrt(abs(P[ii, ii])))) if has_isb else None

    Qn = np.einsum("ab,ebc,dc->ead", Rneu, Pe, Rneu)
    Qf = np.einsum("ab,ebc,dc->ead", Rneu, Pfw, Rneu)
    nsat = np.bincount(ie, minlength=ne)
    pdop = np.full(ne, np.nan)
    for e in range(ne):
        m = o[cut[e]:cut[e + 1]]
        if len(m) >= 4:
            D = Xs[m] - Xe[e]
            try:
                pdop[e] = dop(D / np.linalg.norm(D, axis=1)[:, None], Rneu)["PDOP"]
            except np.linalg.LinAlgError:
                pass
    j = int(np.argmin(np.abs(ep - kin_epoch))) if kin_epoch is not None else ne // 2
    m = o[cut[j]:cut[j + 1]]
    D = Xs[m] - Xe[j]
    out = {"X": Xe[j].copy(), "sX": np.sqrt(np.abs(np.diag(Pe[j]))), "Xk": None, "sXk": None, "has_kin": False,
           "epochs": ep, "clock": clk, "sclock": sclk, "zwd_blocks": ep, "zwd": zwd, "szwd": None,
           "v_phase": v_ph, "v_code": v_cd, "s0": s0, "dof": nobs, "n_obs": 2 * n, "n_par": ns, "n_epochs": ne,
           "n_arcs": ny, "n_zwd": ne, "iter": [], "rms_phase": float(np.sqrt((v_ph ** 2).mean())),
           "rms_code": float(np.sqrt((v_cd[~rej[n:]] ** 2).mean())), "iz": np.arange(n), "ie": ie, "mode": "kalman",
           "isb": isb, "sisb": sisb,
           "track": {"t": ep, "X": Xe, "sNEU": np.sqrt(np.abs(np.einsum("eaa->ea", Qn))), "nsat": nsat, "pdop": pdop, "j": j},
           "dop": dop(D / np.linalg.norm(D, axis=1)[:, None], Rneu) if len(m) >= 4 else {},
           "n_sat_kin": int(len(m)), "clock_kin": float(clk[j]), "sclock_kin": float(sclk[j]),
           "zwd_kin": float(zwd[j]), "j_sel": j,
           "kf": {"dyn": dyn, "q_h": q_h, "q_v": q_v, "q_zwd": q_zwd, "smoothed": do_smooth, "smooth_req": bool(smooth),
                  "n_state": ns, "n_reject": int(n_rej), "n_reset": int(n_reset), "thr": thr,
                  "X_fwd": Xfw, "sNEU_fwd": np.sqrt(np.abs(np.einsum("eaa->ea", Qf))),
                  "resets": reset_log}}
    return out
