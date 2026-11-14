#!/usr/bin/env python3
"""Fit temperatures of maximum density (TMD) and make the paper's figures.

Input : data/summary.json   (from analyze_runs.py; keys "<model>/T<K>")
        data/real_water_density_1atm.csv  (IAPWS reference)
Output: data/tmd_table.csv, data/density_table.csv, figures/fig2..fig4 (pdf + png)

TMD: a weighted cubic is fitted to the points within +-35 K of the highest sampled
density (at least 5 points; otherwise a quadratic), and the maximum of the fitted curve
inside the sampled window is taken. Uncertainty: 2000 Monte Carlo refits with each
density redrawn from a normal distribution of width equal to its blocking error.
TMD_err_K is the standard deviation of the refits whose maximum lies inside the window;
TMD_p16_K/TMD_p84_K are the 16th/84th percentiles of ALL refits, and frac_fits_inside is
the fraction inside. If frac_fits_inside < 0.9 the maximum may lie at or beyond the
window edge and TMD_err_K understates the uncertainty: quote the percentiles instead.
If the highest density sits at the lowest sampled temperature, the TMD is reported
as an upper bound ("below T_min").
"""
import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
plt.rcParams.update({"font.size": 9.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.8})

DIAL = ["t4p2005_l1.05", "t4p2005_l1.00", "t4p2005_l0.97", "t4p2005_l0.95", "t4p2005_l0.93", "t4p2005_l0.90", "t4p2005_l0.85"]
DIAL_SUB = ["t4p2005_l1.05", "t4p2005_l1.00", "t4p2005_l0.95", "t4p2005_l0.90", "t4p2005_l0.85"]   # g(r) panel
STANDARD = ["t4p2005_l1.00", "spce", "tip3p"]
LABEL = {"spce": "SPC/E", "tip3p": "TIP3P", "t4p2005_l1.00": "TIP4P/2005"}
# colour-blind-safe sequential-ish palette for the dial (strong -> weak)
DIAL_COL = {"t4p2005_l1.05": "#7f2704", "t4p2005_l1.00": "#d94801", "t4p2005_l0.97": "#fd8d3c",
            "t4p2005_l0.95": "#08519c", "t4p2005_l0.93": "#3182bd", "t4p2005_l0.90": "#6baed6",
            "t4p2005_l0.85": "#9e9ac8"}
STD_COL = {"t4p2005_l1.00": "#d94801", "spce": "#238b45", "tip3p": "#6a51a3"}


LAM_FIX = {"t4p2005_l0.97": 0.975, "t4p2005_l0.93": 0.925}   # folder names were rounded


def lam_of(model):
    model = model.split("_dt")[0]          # t4p2005_l1.00_dt05 = time-step check at lambda = 1
    if model in LAM_FIX:
        return LAM_FIX[model]
    return float(model.split("_l")[1]) if model.startswith("t4p2005_l") else None


def label_of(model):
    if model in LABEL and model != "t4p2005_l1.00":
        return LABEL[model]
    lam = lam_of(model)
    txt = f"{lam:.3f}"
    txt = txt[:-1] if txt.endswith("0") else txt          # 0.90, 1.00, 1.05, 0.975
    return rf"$\lambda$ = {txt}" + (" (TIP4P/2005)" if abs(lam - 1) < 1e-9 else "")


D_MIN = 0.1   # 1e-5 cm^2/s. Below this the molecules move < ~2.5 nm (RMS) in the 10 ns run: too sluggish to trust


def load():
    with open(os.path.join(BASE, "data", "summary.json")) as fh:
        raw = json.load(fh)
    models = {}
    for key, res in raw.items():
        if "rho" not in res:
            continue
        res["sluggish"] = res.get("D_1e-5_cm2_s", 1.0) < D_MIN
        model, t = key.split("/T")
        models.setdefault(model, []).append((float(t), res))
    for m in models:
        models[m].sort(key=lambda x: x[0])
    exp = np.loadtxt(os.path.join(BASE, "data", "real_water_density_1atm.csv"), delimiter=",",
                     skiprows=1, usecols=(0, 1))
    return models, exp


def fit_tmd(T, rho, sem, rng, nmc=2000):
    T, rho, sem = map(np.asarray, (T, rho, sem))
    imax = int(np.argmax(rho))
    if imax == 0:
        return dict(status="below", bound=float(T[0]))
    if imax == len(T) - 1:
        return dict(status="above", bound=float(T[-1]))
    if T[imax] - T[imax - 1] > 15.0 or T[imax + 1] - T[imax] > 15.0:
        # neighbouring temperatures still missing: the maximum is not bracketed yet
        return dict(status="gap", bound=float(T[imax]))
    sel = np.abs(T - T[imax]) <= 35.0
    Ts, rs, ss = T[sel], rho[sel], np.maximum(sem[sel], 1e-5)
    deg = 3 if len(Ts) >= 5 else 2
    grid = np.linspace(Ts.min(), Ts.max(), 2001)

    def one(r):
        c = np.polyfit(Ts, r, deg, w=1.0 / ss)
        y = np.polyval(c, grid)
        k = int(np.argmax(y))
        return grid[k], y[k]

    t0, r0 = one(rs)
    draws = np.array([one(rs + rng.normal(0.0, ss)) for _ in range(nmc)])
    inside = (draws[:, 0] > Ts.min() + 1e-6) & (draws[:, 0] < Ts.max() - 1e-6)
    p16, p84 = np.percentile(draws[:, 0], [16, 84])        # all refits, including those stuck at the window edge
    return dict(status="ok", tmd=float(t0), tmd_err=float(np.std(draws[inside, 0])) if inside.any() else float("nan"),
                tmd_p16=float(p16), tmd_p84=float(p84),
                rho_max=float(r0), rho_max_err=float(np.std(draws[inside, 1])) if inside.any() else float("nan"),
                deg=deg, npts=int(sel.sum()), frac_inside=float(inside.mean()))


def main():
    rng = np.random.default_rng(20261113)
    models, exp = load()
    # ---------------------------------------------------------- tables
    tmd = {}
    with open(os.path.join(BASE, "data", "density_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model", "lambda", "T_K", "density_g_cm3", "sem_g_cm3", "first_half", "second_half",
                    "hbonds_per_molecule", "tetrahedral_q", "frac_4_hbonds", "E_elec_kJ_per_mol",
                    "E_pot_kJ_per_mol", "D_1e-5_cm2_s", "used_in_fits_and_figures"])
        for m, pts in sorted(models.items()):
            for t, r in pts:
                w.writerow([m, lam_of(m) or "", int(t), f"{r['rho']:.5f}", f"{r['rho_sem']:.5f}",
                            f"{r['rho_first_half']:.5f}", f"{r['rho_second_half']:.5f}",
                            f"{r.get('hb_mean', float('nan')):.3f}", f"{r.get('q_mean', float('nan')):.4f}",
                            f"{r.get('frac4_mean', float('nan')):.3f}", f"{r.get('E_elec_per_mol_kJ', float('nan')):.3f}",
                            f"{r.get('E_pot_per_mol_kJ', float('nan')):.3f}", f"{r.get('D_1e-5_cm2_s', float('nan')):.4f}",
                            "no (too sluggish)" if r["sluggish"] else "yes"])
    # fits and figures use only the points where the liquid is mobile enough to have settled
    models = {m: [p for p in pts if not p[1]["sluggish"]] for m, pts in models.items()}
    models = {m: pts for m, pts in models.items() if pts}
    for m, pts in models.items():
        T = [p[0] for p in pts]; rho = [p[1]["rho"] for p in pts]; sem = [p[1]["rho_sem"] for p in pts]
        if len(T) >= 4:
            tmd[m] = fit_tmd(T, rho, sem, rng)
    # real water: vertex of the parabola through the highest grid point and its two neighbours (grid is 0.5 K)
    iexp = int(np.argmax(exp[:, 1]))
    cexp = np.polyfit(exp[iexp - 1:iexp + 2, 0] - exp[iexp, 0], exp[iexp - 1:iexp + 2, 1], 2)
    t_real = exp[iexp, 0] - cexp[1] / (2 * cexp[0])
    rho_real = float(np.polyval(cexp, t_real - exp[iexp, 0]))
    with open(os.path.join(BASE, "data", "tmd_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model", "lambda", "status", "TMD_K", "TMD_err_K", "rho_max", "rho_max_err", "fit_deg", "npts",
                    "frac4_at_TMD", "hbonds_at_TMD", "q_at_TMD", "TMD_p16_K", "TMD_p84_K", "frac_fits_inside"])
        w.writerow(["real water (IAPWS)", "", "ok", f"{t_real:.2f}", "", f"{rho_real:.6f}", "", "", "", "", "", "", "", "", ""])
        for m in sorted(tmd):
            d = tmd[m]
            if d["status"] == "ok":
                # structure at the maximum: linear interpolation between the sampled temperatures
                st = [(p[0], p[1]) for p in models[m] if "frac4_mean" in p[1]]
                at = {}
                if len(st) >= 2 and st[0][0] <= d["tmd"] <= st[-1][0]:      # interpolate, never extrapolate
                    Ts_ = np.array([p[0] for p in st])
                    for key in ("frac4_mean", "hb_mean", "q_mean"):
                        at[key] = float(np.interp(d["tmd"], Ts_, [p[1][key] for p in st]))
                d.update({k + "_at_tmd": v for k, v in at.items()})
                w.writerow([m, lam_of(m) or "", "ok", f"{d['tmd']:.1f}", f"{d['tmd_err']:.1f}", f"{d['rho_max']:.5f}",
                            f"{d['rho_max_err']:.5f}", d["deg"], d["npts"],
                            f"{at.get('frac4_mean', float('nan')):.3f}", f"{at.get('hb_mean', float('nan')):.3f}",
                            f"{at.get('q_mean', float('nan')):.4f}", f"{d['tmd_p16']:.1f}", f"{d['tmd_p84']:.1f}",
                            f"{d['frac_inside']:.3f}"])
            else:
                w.writerow([m, lam_of(m) or "", d["status"], f"{d['bound']:.0f}", "", "", "", "", "", "", "", "", "", "", ""])
    for m in sorted(tmd):
        print(m, tmd[m])

    # ---------------------------------------------------------- fig 2: standard models
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.plot(exp[:, 0] - 273.15, exp[:, 1], color="k", lw=1.6, label="real water")
    for m in STANDARD:
        if m not in models:
            continue
        T = np.array([p[0] for p in models[m]]); r = np.array([p[1]["rho"] for p in models[m]])
        s = np.array([p[1]["rho_sem"] for p in models[m]])
        ax.errorbar(T - 273.15, r, yerr=s, fmt="o-", ms=3.5, lw=1, capsize=1.5, color=STD_COL[m], label=LABEL[m])
        d = tmd.get(m, {})
        if d.get("status") == "ok":
            ax.axvline(d["tmd"] - 273.15, color=STD_COL[m], lw=0.7, ls=":")
    ax.set_xlabel("temperature (°C)"); ax.set_ylabel(r"density (g cm$^{-3}$)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(BASE, "figures", "fig2_models.pdf")); fig.savefig(os.path.join(BASE, "figures", "fig2_models.png"), dpi=200)
    plt.close(fig)

    # ---------------------------------------------------------- fig 3: the dial
    # (a) absolute density for five settings; (b) T_MD of all seven settings against lambda^2.
    # One shared legend under both panels.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.7), gridspec_kw={"width_ratios": [1.35, 1]})
    handles = {}
    h, = a1.plot(exp[:, 0] - 273.15, exp[:, 1], color="k", lw=1.8)
    handles["real water"] = h
    ylo, yhi = 1.0, 1.0
    for m in DIAL_SUB:
        if m not in models:
            continue
        T = np.array([p[0] for p in models[m]]); r = np.array([p[1]["rho"] for p in models[m]])
        s_ = np.array([p[1]["rho_sem"] for p in models[m]])
        h = a1.errorbar(T - 273.15, r, yerr=s_, fmt="o-", ms=3.2, lw=1, capsize=1.5, color=DIAL_COL[m])
        handles[label_of(m)] = h
        ylo, yhi = min(ylo, r.min()), max(yhi, r.max())
        d = tmd.get(m, {})
        if d.get("status") == "ok":
            a1.plot(d["tmd"] - 273.15, d["rho_max"] + 0.003, marker="v", ms=6.5, color=DIAL_COL[m], mec="k", mew=0.6, zorder=5)
    a1.set_ylim(ylo - 0.003, yhi + 0.007)
    a1.set_xlabel("temperature (°C)"); a1.set_ylabel(r"density (g cm$^{-3}$)")
    a1.set_title("(a)", loc="left", fontsize=9.5)
    t1 = tmd.get("t4p2005_l1.00", {}).get("tmd", 277.0)      # lambda^2 line through the standard model
    x = np.linspace(0.68, 1.14, 10)
    a2.plot(x, t1 * x - 273.15, ls="--", color="0.4", lw=1)
    a2.text(0.80, t1 * 0.80 - 273.15 + 9, r"$T_\mathrm{MD}\propto\lambda^2$", fontsize=8, color="0.3", ha="right")
    # straight line through the bracketed maxima of the dial (the simulations' own trend)
    fx = [lam_of(m) ** 2 for m in DIAL if tmd.get(m, {}).get("status") == "ok"]
    fy = [tmd[m]["tmd"] for m in DIAL if tmd.get(m, {}).get("status") == "ok"]
    if len(fx) >= 3:
        cfit = np.polyfit(fx, fy, 1)
        xx = np.linspace(min(fx) - 0.04, max(fx) + 0.03, 10)
        a2.plot(xx, np.polyval(cfit, xx) - 273.15, color="#d94801", lw=0.9, alpha=0.6, zorder=1)
        a2.text(0.69, 30, f"simulations: {cfit[0] / 100:.1f} K\nper 1% of electrostatic energy", fontsize=7.5,
                color="#d94801", ha="left")
    a2.axhline(t_real - 273.15, color="k", lw=0.8, ls=":")
    a2.text(0.69, t_real - 273.15 + 3, "real water, 4 °C", fontsize=7.5)
    for m in DIAL:
        d = tmd.get(m, {})
        if not d:
            continue
        x2 = lam_of(m) ** 2
        if d.get("status") == "ok":
            # poorly bracketed maximum (many refits at the window edge): show the 16-84 % range of all refits
            yerr = d["tmd_err"] if d["frac_inside"] >= 0.9 else [[d["tmd"] - d["tmd_p16"]], [d["tmd_p84"] - d["tmd"]]]
            h = a2.errorbar(x2, d["tmd"] - 273.15, yerr=yerr, fmt="o", ms=5.5, capsize=2, color=DIAL_COL[m],
                            mec="k", mew=0.5, zorder=4)
        elif d.get("status") in ("below", "above"):
            sgn = -1 if d["status"] == "below" else 1
            y0 = d["bound"] - 273.15
            a2.annotate("", xy=(x2, y0 + sgn * 14), xytext=(x2, y0),
                        arrowprops=dict(arrowstyle="-|>", color=DIAL_COL[m], lw=1.4))
            h, = a2.plot(x2, y0, marker="_", ms=9, mew=2, color=DIAL_COL[m], ls="none")
        else:
            continue
        handles.setdefault(label_of(m), h)
    lows = [tmd[m]["bound"] - 22 for m in DIAL if tmd.get(m, {}).get("status") == "below"]
    lows += [tmd[m]["tmd_p16"] - 8 for m in DIAL if tmd.get(m, {}).get("status") == "ok"]   # keep error bars visible
    if lows:
        a2.set_ylim(bottom=min(lows) - 273.15)
    a2.set_xlabel(r"$\lambda^2$ (electrostatic energy, relative)")
    a2.set_ylabel("temperature of maximum density (°C)")
    a2.set_title("(b)", loc="left", fontsize=9.5)
    order = ["real water"] + [label_of(m) for m in DIAL]
    hs = [handles[k] for k in order if k in handles]
    fig.legend(hs, [k for k in order if k in handles], loc="lower center", ncol=4, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(os.path.join(BASE, "figures", "fig3_dial.pdf")); fig.savefig(os.path.join(BASE, "figures", "fig3_dial.png"), dpi=200)
    plt.close(fig)

    # ---------------------------------------------------------- fig 4: structure
    have_struct = any("hb_mean" in p[1] for m in DIAL if m in models for p in models[m])
    if have_struct:
        fig, (b1, b2) = plt.subplots(1, 2, figsize=(7.2, 3.1))
        # (a) hydrogen bonds per molecule; triangles: value at each model's density maximum
        at_max = []
        for m in DIAL_SUB + ["spce", "tip3p"]:
            if m not in models:
                continue
            pts = [(p[0], p[1]) for p in models[m] if "hb_mean" in p[1]]
            if not pts:
                continue
            col = DIAL_COL.get(m, STD_COL.get(m))
            T = np.array([p[0] for p in pts]); hb = np.array([p[1]["hb_mean"] for p in pts])
            b1.plot(T - 273.15, hb, "o", ms=3.0, lw=1, color=col, label=label_of(m),
                    ls="-" if m.startswith("t4p") else "--")
            d = tmd.get(m, {})
            if d.get("status") == "ok" and T.min() <= d["tmd"] <= T.max():      # interpolate, never extrapolate
                y = float(np.interp(d["tmd"], T, hb)); at_max.append(y)
                b1.plot(d["tmd"] - 273.15, y, marker="v", ms=7, color=col, mec="k", mew=0.6, zorder=5)
        if at_max:
            b1.axhline(np.mean(at_max), color="0.5", lw=0.8, ls=":")
        b1.set_xlabel("temperature (°C)"); b1.set_ylabel("hydrogen bonds per molecule")
        b1.legend(frameon=False, fontsize=6.8, ncol=2, loc="lower left"); b1.set_title("(a)", loc="left", fontsize=9.5)
        # (b) at corresponding temperatures T = lambda^2 * Tref (Eq. 6 of the paper): the settings then differ only in
        # the depth of the Lennard-Jones well (eps/lambda^2), so the panel isolates its effect. Nearest simulated T.
        Tref = float(sys.argv[1]) if len(sys.argv) > 1 else 260.0
        for m in DIAL_SUB:
            if m not in models:
                continue
            cand = [p for p in models[m] if "g_oo" in p[1]]
            if not cand:
                continue
            target = lam_of(m) ** 2 * Tref
            t_, res = min(cand, key=lambda p: abs(p[0] - target))
            if abs(t_ - target) > 6.0:
                continue
            r = np.array(res["r_nm"]) * 10; g = np.array(res["g_oo"])
            lab = label_of(m).replace(" (TIP4P/2005)", "")
            b2.plot(r, g, lw=1.2, color=DIAL_COL[m], label=f"{lab}, {t_:.0f} K")
        # minimum of the oxygen-oxygen Lennard-Jones potential, 2^(1/6) sigma = 3.55 A (not changed by the dial)
        b2.axvline(3.546, color="0.45", ls=":", lw=0.9)
        b2.text(3.62, 0.97, "bottom of the O–O\nLennard-Jones well", transform=b2.get_xaxis_transform(), fontsize=7, va="top", color="0.3")
        b2.set_xlim(2.2, 8.0); b2.set_xlabel("oxygen–oxygen distance (Å)"); b2.set_ylabel(r"$g_\mathrm{OO}(r)$")
        b2.set_title(rf"(b) at $T = \lambda^2 \times {Tref:.0f}$ K", loc="left", fontsize=9.5)
        b2.legend(frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 0.80))
        fig.tight_layout(); fig.savefig(os.path.join(BASE, "figures", "fig4_structure.pdf")); fig.savefig(os.path.join(BASE, "figures", "fig4_structure.png"), dpi=200)
        plt.close(fig)


if __name__ == "__main__":
    main()
