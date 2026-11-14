"""Classroom activity: find the temperature at which (simulated) water is densest.

Needs Python 3 with numpy and matplotlib (both included in Anaconda, or use Google Colab).
The data file density_vs_temperature.csv has one row per simulation:
    model, lambda, T_K, density_g_cm3, sem_g_cm3, ...
"lambda" is the hydrogen-bond dial: 1 = the TIP4P/2005 model, smaller = weaker hydrogen bonds.
The script skips the few runs marked "no" in the last column: they were too cold and sluggish to settle.

Tasks
1. Run the script as it is. It plots the density of the TIP4P/2005 model and of real water.
2. Read off the temperature of the highest density. Compare with real water (3.98 C).
3. Fit a parabola rho = a + b*T + c*T^2 to the points near the maximum (done below) and
   find the maximum from the vertex formula T = -b/(2c). Is it close to your reading?
4. Change MODEL below to another dial setting and repeat. How does the maximum move?
   (If the densest point is the first or last one, the script says so instead of fitting:
   the maximum then lies outside the simulated range.)
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

MODEL = "t4p2005_l1.00"      # try t4p2005_l1.05, t4p2005_l0.95, spce, tip3p ...

T, rho, err = [], [], []
with open("density_vs_temperature.csv") as fh:
    for row in csv.DictReader(fh):
        if row["model"] == MODEL and row.get("used_in_fits_and_figures", "yes") == "yes":
            # (runs flagged "no" were too sluggish to settle: the molecules barely moved)
            T.append(float(row["T_K"]) - 273.15)      # in degrees Celsius
            rho.append(float(row["density_g_cm3"]))
            err.append(float(row["sem_g_cm3"]))
T, rho, err = np.array(T), np.array(rho), np.array(err)

real = np.loadtxt("real_water_density_1atm.csv", delimiter=",", skiprows=1, usecols=(0, 1))

plt.errorbar(T, rho, yerr=err, fmt="o", label=f"simulation: {MODEL}")
plt.plot(real[:, 0] - 273.15, real[:, 1], "k-", label="real water")

k = int(np.argmax(rho))
if k == 0 or k == len(T) - 1:
    # the densest point is at the end of the data: a parabola would extrapolate nonsense
    side = "below the coldest" if k == 0 else "above the warmest"
    message = f"{MODEL}: no maximum inside the data; it lies {side} temperature simulated ({T[k]:.0f} °C)"
else:
    # parabola through the (up to) five points around the highest density
    near = slice(max(k - 2, 0), min(k + 3, len(T)))
    c, b, a = np.polyfit(T[near], rho[near], 2)
    T_max = -b / (2 * c)
    Tfit = np.linspace(T[near].min(), T[near].max(), 100)
    plt.plot(Tfit, a + b * Tfit + c * Tfit**2, "--", label=f"parabola: maximum at {T_max:.1f} °C")
    message = f"{MODEL}: the parabola puts the density maximum at {T_max:.1f} °C"

plt.xlabel("temperature (°C)")
plt.ylabel("density (g/cm³)")
plt.legend()
plt.show()
print(message)
