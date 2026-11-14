#!/usr/bin/env python3
"""Seawater: temperature of maximum density and freezing point vs Absolute Salinity (TEOS-10, p = 0).
Run: uv run --no-project --with gsw --with scipy python analysis/seawater_tmd.py
Result (2026-11-13): slope ~0.22 K per g/kg; TMD = freezing point at SA = 23.85 g/kg."""
import numpy as np, gsw
# Temperature of maximum density of seawater at the sea surface (p = 0 dbar) vs Absolute Salinity, TEOS-10
T = np.linspace(-3, 5, 8001)
print("SA_g/kg  TMD_C   T_freeze_C")
for SA in [0, 5, 10, 15, 20, 24, 24.5, 25, 26, 30, 35]:
    CT = gsw.CT_from_t(SA, T, 0)
    rho = gsw.rho(SA, CT, 0)
    i = np.argmax(rho)
    tf = gsw.t_freezing(SA, 0, 0)   # air-free seawater
    print(f"{SA:6.1f}  {T[i]:6.2f}  {tf:6.2f}", "(max at edge)" if i in (0, len(T)-1) else "")
# crossing salinity
from scipy.optimize import brentq
def f(SA):
    CT = gsw.CT_from_t(SA, T, 0); return T[np.argmax(gsw.rho(SA, CT, 0))] - gsw.t_freezing(SA, 0, 0)
print("crossing SA =", round(brentq(f, 20, 28), 2), "g/kg")
