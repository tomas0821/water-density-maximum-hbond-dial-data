"""Real-water density at 0.101325 MPa for the comparison curve.

238-300 K: IAPWS G12-15 guideline for supercooled water (Holten, Sengers & Anisimov,
J. Phys. Chem. Ref. Data 43, 043101, 2014), as implemented in the `iapws` package.
300-330 K: IAPWS-95 (Wagner & Pruss, J. Phys. Chem. Ref. Data 31, 387, 2002).
Run with:  uv run --no-project --with iapws python analysis/experiment.py
"""
import csv
import numpy as np
from iapws import IAPWS95
from iapws._iapws import _Supercooled

rows = []
for T in np.arange(238.0, 330.01, 0.5):
    if T <= 300.0:
        rho, src = _Supercooled(T, 0.101325)["rho"], "IAPWS G12-15"
    else:
        rho, src = IAPWS95(T=T, P=0.101325).rho, "IAPWS-95"
    rows.append((round(T, 2), round(rho / 1000.0, 6), src))
with open("data/real_water_density_1atm.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["T_K", "density_g_cm3", "source"])
    w.writerows(rows)
Ts = np.array([r[0] for r in rows]); rh = np.array([r[1] for r in rows])
print(len(rows), "rows; maximum at", Ts[rh.argmax()], "K, rho_max =", rh.max())
