"""Fiziksel sabitler ve GPS frekanslari."""
import math

C = 299792458.0                      # isik hizi [m/s]
F0 = 10.23e6
F1 = 154 * F0                        # 1575.42 MHz
F2 = 120 * F0                        # 1227.60 MHz
LAM1 = C / F1
LAM2 = C / F2
GAMMA = (F1 / F2) ** 2               # 1.6469
A3 = F1 ** 2 / (F1 ** 2 - F2 ** 2)   # 2.5457  (L3 katsayisi)
B3 = -F2 ** 2 / (F1 ** 2 - F2 ** 2)  # -1.5457
LAM_W = C / (F1 - F2)                # genis serit 0.862 m
LAM_NL = C / (F1 + F2)               # dar serit 0.107 m
MS = C * 1e-3                        # 1 ms -> 299 792.458 m
TECU_L1 = 40.3e16 / F1 ** 2          # 1 TECU -> L1'de 0.162 m

OMEGA_E = 7.2921151467e-5            # Dunya donus hizi [rad/s]
GM_E = 3.986004418e14
GM_S = 1.32712442076e20
GM_M = 4.9028e12
RE = 6378136.6
WGS_A = 6378137.0
WGS_F = 1 / 298.257223563
GPS_MU = 3.986005e14                 # IS-GPS-200
GPS_OMEGA = 7.2921151467e-5
F_REL = -4.442807633e-10
DEG = math.pi / 180
