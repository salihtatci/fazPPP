"""Komut satirindan calistirma: python py/cli.py gozlem.o yorunge.sp3 [--clk x.clk] [--nav x.n] [--met x.m] ...

Sonuc JSON olarak yazilir. Ornekler:
    python py/cli.py examples/ista_2014_240/ISTA2400.14o examples/ista_2014_240/igs18074.sp3 \
        --nav examples/ista_2014_240/ISTA2400.14n --met examples/ista_2014_240/ISTA2400.14m \
        --epoch 13:15:00 --out sonuc.json

    # hareketli alici, 30 s saat dosyasiyla kinematik, GPS + Galileo
    python py/cli.py gozlem.23O COD0MGXFIN_..._ORB.SP3.gz --clk COD0MGXFIN_..._CLK.CLK.gz --mode kinematic
"""
import argparse
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gnssfaz import pipeline  # noqa: E402


def _read(path):
    if path is None:
        return None
    raw = open(path, "rb").read()
    if path.endswith(".gz"):
        raw = gzip.decompress(raw)
    return raw.decode("latin-1")


def main():
    ap = argparse.ArgumentParser(description="GPS tasiyici faz ile istasyon konumu (FazPPP)")
    ap.add_argument("obs")
    ap.add_argument("sp3", nargs="+")
    ap.add_argument("--clk", nargs="+", help="RINEX saat dosyalari (ornegin 30 s)")
    ap.add_argument("--mode", choices=["static", "kinematic", "kalman"], default="static",
                    help="static: statik EKK; kinematic: toplu kinematik EKK; kalman: Kalman filtresi + RTS")
    ap.add_argument("--kf-dyn", choices=["pv", "rw", "static"], default="pv",
                    help="Kalman dinamik modeli: pv sabit hiz, rw rastgele yuruyus, static sabit konum")
    ap.add_argument("--kf-qh", type=float, help="yatay surec gurultusu (pv: m^2/s^3, rw: m^2/s)")
    ap.add_argument("--kf-qv", type=float, help="dusey surec gurultusu")
    ap.add_argument("--no-smooth", action="store_true", help="RTS yumusatmayi kapat (yalniz ileri filtre)")
    ap.add_argument("--rate", type=float, help="islem araligi [s]; verilmezse CLK (yoksa SP3) araligi")
    ap.add_argument("--systems", default="GE", help="G veya GE")
    ap.add_argument("--min-arc", type=float, default=10.0, help="en kisa kullanilir yay [dakika]")
    ap.add_argument("--nav")
    ap.add_argument("--met")
    ap.add_argument("--atx")
    ap.add_argument("--dcb")
    ap.add_argument("--epoch", help="SS:DD:ss (GPS zamani)")
    ap.add_argument("--cutoff", type=float, default=10.0)
    ap.add_argument("--ref", nargs=3, type=float, metavar=("X", "Y", "Z"))
    ap.add_argument("--out", default="fazppp_sonuc.json")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    files = {"obs": _read(a.obs), "sp3": [_read(p) for p in a.sp3], "nav": _read(a.nav), "met": _read(a.met),
             "atx": _read(a.atx), "dcb": _read(a.dcb), "clk": [_read(p) for p in a.clk] if a.clk else None}
    res = pipeline.run(files, {"epoch": a.epoch, "cutoff": a.cutoff, "ref_xyz": a.ref, "mode": a.mode,
                               "rate": a.rate, "systems": a.systems, "min_arc_min": a.min_arc,
                               "kf_dyn": a.kf_dyn, "kf_qh": a.kf_qh, "kf_qv": a.kf_qv, "kf_smooth": not a.no_smooth},
                       progress=lambda m, f: print(f"[{f:4.0%}] {m}", file=sys.stderr),
                       default_atx={"igs08": os.path.join(here, "..", "data", "igs08.atx.gz"),
                                    "igs20": os.path.join(here, "..", "data", "igs20.atx.gz")})
    json.dump(res, open(a.out, "w"), indent=1)
    L = res["l3"]
    if res.get("kin"):
        K = res["kin"]["summary"]
        print(f"Kinematik L3: {K['n_epochs']} epok ({K['first']}-{K['last']}, {K['rate']:g} s), "
              f"yol {K['path_m']} m, medyan sigma NEU {K['sig_med_mm']} mm")
        print(f"  {res['epoch']['used']}: X={L['X'][0]:.4f} Y={L['X'][1]:.4f} Z={L['X'][2]:.4f}")
        if K.get("kf"):
            F = K["kf"]
            print(f"  Kalman ({F['dyn']}, q={F['q_h']}/{F['q_v']}, yumusatma={F['smoothed']}): "
                  f"atilan kod {F['n_reject']}, sifirlanan yay {F['n_reset']}, yakinsama {F['conv_min']} dk, "
                  f"ileri-yumusatilmis RMS NEU {F['fwd_minus_smooth_rms_mm']} mm (son yari {F['fwd_minus_smooth_rms_late_mm']})")
    else:
        print(f"Statik L3 : X={L['X'][0]:.4f} Y={L['X'][1]:.4f} Z={L['X'][2]:.4f}")
    if L.get("isb_ns") is not None:
        print(f"GPS-Galileo ISB: {L['isb_ns']} +- {L['sisb_ns']} ns")
    print(f"s0={L['s0']}  faz RMS={L['rms_phase_mm']} mm  kod RMS={L['rms_code_m']} m")
    if L["Xk"]:
        print(f"Kinematik {res['epoch']['used']}: X={L['Xk'][0]:.4f} Y={L['Xk'][1]:.4f} Z={L['Xk'][2]:.4f}")
    print(f"Sonuc: {a.out}")


if __name__ == "__main__":
    main()
