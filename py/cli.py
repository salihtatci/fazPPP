"""Komut satirindan calistirma: python py/cli.py gozlem.o yorunge.sp3 [--nav x.n] [--met x.m] ...

Sonuc JSON olarak yazilir. Ornek:
    python py/cli.py examples/ista_2014_240/ISTA2400.14o examples/ista_2014_240/igs18074.sp3 \
        --nav examples/ista_2014_240/ISTA2400.14n --met examples/ista_2014_240/ISTA2400.14m \
        --epoch 13:15:00 --out sonuc.json
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
             "atx": _read(a.atx), "dcb": _read(a.dcb)}
    res = pipeline.run(files, {"epoch": a.epoch, "cutoff": a.cutoff, "ref_xyz": a.ref},
                       progress=lambda m, f: print(f"[{f:4.0%}] {m}", file=sys.stderr),
                       default_atx=os.path.join(here, "..", "data", "igs08.atx.gz"))
    json.dump(res, open(a.out, "w"), indent=1)
    L = res["l3"]
    print(f"Statik L3 : X={L['X'][0]:.4f} Y={L['X'][1]:.4f} Z={L['X'][2]:.4f}")
    if L["Xk"]:
        print(f"Kinematik {res['epoch']['used']}: X={L['Xk'][0]:.4f} Y={L['Xk'][1]:.4f} Z={L['Xk'][2]:.4f}")
    print(f"Sonuc: {a.out}")


if __name__ == "__main__":
    main()
