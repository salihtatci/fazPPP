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
        # her epokta etkin global sutunlar (o epokta gorunen yaylar, ZWD, ISB, statik konum): (ne, K)
        if it == 0:
            key = erow[gr].astype(np.int64) * (ng + 1) + gc
            ukey = np.unique(key)
            ue = (ukey // (ng + 1)).astype(np.intp)          # Pyodide (wasm32): indisler intp olmali
            starts_e = np.searchsorted(ue, np.arange(ne))
            cnt_e = np.bincount(ue, minlength=ne)
            K = int(cnt_e.max())
            Cidx = np.full((ne, K), ng, dtype=np.intp)          # ng: bos (sifir) sutun
            posk = (np.arange(len(ukey)) - starts_e[ue]).astype(np.intp)
            Cidx[ue, posk] = (ukey % (ng + 1)).astype(np.intp)
            kk = posk[np.searchsorted(ukey, key)]            # her ogenin epok icindeki sirasi
            ek = erow[gr]
        NegS = np.zeros((ne, K, 4))
        np.add.at(NegS, (ek, kk), pa[gr] * gv[:, None])
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
        NegT = np.ascontiguousarray(np.transpose(NegS, (0, 2, 1)))  # (ne, 4, K)
        H = np.matmul(Ninv, NegT)                                   # (ne, 4, K)
        Mred = np.matmul(np.transpose(NegT, (0, 2, 1)), H)          # (ne, K, K)
        Nx = np.zeros((ng + 1, ng + 1))
        np.add.at(Nx, (np.broadcast_to(Cidx[:, :, None], Mred.shape), np.broadcast_to(Cidx[:, None, :], Mred.shape)), Mred)
        Nred = Ngg - Nx[:ng, :ng]
        bx = np.zeros(ng + 1)
        np.add.at(bx, Cidx, np.matmul(np.transpose(H, (0, 2, 1)), be[:, :, None])[:, :, 0])
        bred = bg - bx[:ng]
        try:
            dg = np.linalg.solve(Nred, bred)
        except np.linalg.LinAlgError:
            dg = np.linalg.lstsq(Nred, bred, rcond=None)[0]
        dgx = np.r_[dg, 0.0]
        de = np.matmul(Ninv, (be - np.matmul(NegT, dgx[Cidx][:, :, None])[:, :, 0])[:, :, None])[:, :, 0]
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
    Qx = np.zeros((ng + 1, ng + 1))
    Qx[:ng, :ng] = Qgg
    Qsub = Qx[Cidx[:, :, None], Cidx[:, None, :]]                  # (ne, K, K)
    HQ = np.matmul(H, Qsub)
    Qee = Ninv + np.matmul(HQ, np.transpose(H, (0, 2, 1)))
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
        Qn = np.matmul(np.matmul(Rneu, Qp), Rneu.T)
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
