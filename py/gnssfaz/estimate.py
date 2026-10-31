"""En kucuk kareler: faz + kod cozumu (statik ve secilen epokta kinematik), DOP, kod SPP."""
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


def phase_solution(d, phase, code, corr_phase, corr_code, X0, kin_epoch, Rneu,
                   sig_phase=0.01, sig_code=1.0, zwd_block=7200.0, iters=5):
    """Faz + kod EKK cozumu.

    d: dict, kayit dizileri (t, arc, Xs(n,3), el, cdts, Tmodel, mw)
    phase, code: olculer [m]; corr_*: model uzakliga eklenecek duzeltmeler [m]
    Bilinmeyenler: statik konum (3), kin_epoch'ta ayri konum (3), her epokta saat,
    zwd_block saniyelik islak troposfer duzeltmeleri, her yay icin ondalikli belirsizlik.
    """
    t = d["t"]
    n = len(t)
    ep, ie = np.unique(t, return_inverse=True)
    yy, iy = np.unique(d["arc"], return_inverse=True)
    zb = np.floor((t - t.min()) / zwd_block).astype(int)
    zz, iz = np.unique(zb, return_inverse=True)
    ne, ny, nz = len(ep), len(yy), len(zz)
    km = np.abs(t - kin_epoch) < 1e-3 if kin_epoch is not None else np.zeros(n, bool)
    has_kin = bool(km.sum() >= 4)
    if not has_kin:
        km[:] = False
    i_s, i_z = 6, 6 + ne
    i_b = i_z + nz
    npar = i_b + ny
    Xs = d["Xs"]
    s2 = np.sin(np.radians(d["el"])) ** 2
    cdts = d["cdts"]
    Tm = d["Tmodel"]
    mw = d["mw"]
    X = np.array(X0, float).copy()
    Xk = X.copy()
    rho0 = np.linalg.norm(Xs - X, axis=1)
    k0 = code - rho0 + cdts - Tm - corr_code
    clk = np.array([np.median(k0[ie == j]) for j in range(ne)])
    zwd = np.zeros(nz)
    B = np.zeros(ny)
    p = np.concatenate([s2 / sig_phase ** 2, s2 / sig_code ** 2])
    hist = []
    for it in range(iters):
        Xa = np.where(km[:, None], Xk, X)
        D = Xs - Xa
        rho = np.linalg.norm(D, axis=1)
        u = D / rho[:, None]
        com = rho + clk[ie] - cdts + Tm + zwd[iz] * mw
        l = np.concatenate([phase - (com + B[iy] + corr_phase), code - (com + corr_code)])
        A = np.zeros((2 * n, npar))
        for blk, r in enumerate((np.arange(n), np.arange(n) + n)):
            A[r[~km], 0:3] = -u[~km]
            if has_kin:
                A[r[km], 3:6] = -u[km]
            A[r, i_s + ie] = 1.0
            A[r, i_z + iz] = mw
            if blk == 0:
                A[r, i_b + iy] = 1.0
        cols = np.ones(npar, bool)
        if not has_kin:
            cols[3:6] = False
        Ac = A[:, cols]
        N = Ac.T @ (Ac * p[:, None])
        dxc = np.linalg.solve(N, Ac.T @ (p * l))
        dx = np.zeros(npar)
        dx[cols] = dxc
        X += dx[0:3]
        Xk += dx[3:6]
        clk += dx[i_s:i_z]
        zwd += dx[i_z:i_b]
        B += dx[i_b:]
        hist.append(float(np.linalg.norm(dx[0:3])))
        if it >= 1 and hist[-1] < 1e-5:
            break
    v = Ac @ dxc - l
    dof = 2 * n - int(cols.sum())
    s0 = float(math.sqrt((v ** 2 * p).sum() / dof))
    Q = np.zeros((npar, npar))
    Q[np.ix_(cols, cols)] = np.linalg.inv(N)
    sig = s0 * np.sqrt(np.abs(np.diag(Q)))
    out = {"X": X, "sX": sig[0:3], "Xk": Xk if has_kin else None, "sXk": sig[3:6] if has_kin else None,
           "epochs": ep, "clock": clk, "sclock": sig[i_s:i_z], "zwd_blocks": zz * zwd_block + t.min(),
           "zwd": zwd, "v_phase": v[:n], "v_code": v[n:], "s0": s0, "dof": dof, "n_obs": 2 * n,
           "n_par": int(cols.sum()), "n_epochs": ne, "n_arcs": ny, "n_zwd": nz, "iter": hist,
           "rms_phase": float(np.sqrt((v[:n] ** 2).mean())), "rms_code": float(np.sqrt((v[n:] ** 2).mean())),
           "has_kin": has_kin, "iz": iz, "ie": ie}
    if has_kin:
        D = Xs[km] - Xk
        out["dop"] = dop(D / np.linalg.norm(D, axis=1)[:, None], Rneu)
        out["n_sat_kin"] = int(km.sum())
        j = int(np.searchsorted(ep, kin_epoch - 1e-3))
        out["clock_kin"] = float(clk[j])
        out["sclock_kin"] = float(sig[i_s + j])
        out["zwd_kin"] = float(zwd[iz[np.argmax(km)]])
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
