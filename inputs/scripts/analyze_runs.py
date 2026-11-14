#!/usr/bin/env python3
"""Per-state-point analysis for the water hydrogen-bond-dial study (numpy only).

Run on the cluster from /path/to/water-tmd:
    python3 analyze_runs.py runs <runlist> out.json
For every finished state point (DONE marker) it reports
  * mean density and its standard error from Flyvbjerg-Petersen blocking of density.xvg,
    plus first-half / second-half means as a drift check;
  * if frames.gro exists (written by gmx trjconv, O and H positions every 100 ps):
      - hydrogen bonds per molecule, geometric criterion of Luzar & Chandler (1996):
        R_OO < 0.35 nm and angle(O_d-H, O_d...O_a) < 30 degrees;
      - Errington-Debenedetti tetrahedral order q of each oxygen's 4 nearest oxygens;
      - fraction of molecules with exactly four hydrogen bonds;
      - g_OO(r) up to 0.9 nm.
"""
import json
import math
import os
import sys

import numpy as np


def read_xvg(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line[:1] in "#@":
                continue
            parts = line.split()
            if len(parts) >= 2:
                rows.append((float(parts[0]), float(parts[1])))
    return np.array(rows)


def read_xvg_all(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line[:1] in "#@":
                continue
            parts = line.split()
            if parts:
                rows.append([float(v) for v in parts])
    return np.array(rows)


def blocking_sem(x):
    """Flyvbjerg-Petersen blocking: return the plateau estimate of the SEM."""
    x = np.asarray(x, float)
    sems = []
    while len(x) >= 32:
        n = len(x)
        sems.append(np.std(x, ddof=1) / math.sqrt(n))
        if n % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
    if not sems:
        return float("nan"), []
    # plateau: take the maximum over the last levels that still have >= 32 blocks
    return float(max(sems[-4:])), [float(s) for s in sems]


def read_gro_frames(path):
    """Yield (names, xyz[nm], box[3]) for each frame of a multi-frame .gro file."""
    with open(path) as fh:
        while True:
            title = fh.readline()
            if not title:
                return
            n = int(fh.readline())
            names, xyz = [], np.empty((n, 3))
            for i in range(n):
                line = fh.readline()
                names.append(line[10:15].strip())
                xyz[i] = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
            box = np.array([float(v) for v in fh.readline().split()[:3]])
            yield names, xyz, box


def frame_structure(names, xyz, box, rbins):
    names = np.array(names)
    O = xyz[names == "OW"]
    H1 = xyz[names == "HW1"]
    H2 = xyz[names == "HW2"]
    n = len(O)
    d = O[None, :, :] - O[:, None, :]
    d -= box * np.round(d / box)
    r2 = np.einsum("ijk,ijk->ij", d, d)
    np.fill_diagonal(r2, np.inf)
    r = np.sqrt(r2)
    # g_OO(r)
    hist, _ = np.histogram(r[np.triu_indices(n, 1)], bins=rbins)
    # tetrahedral order from the 4 nearest O neighbours
    nn = np.argsort(r2, axis=1)[:, :4]
    vec = d[np.arange(n)[:, None], nn]                   # (n,4,3) vectors i -> j
    vec /= np.linalg.norm(vec, axis=2)[:, :, None]
    q = np.ones(n)
    for a in range(3):
        for b in range(a + 1, 4):
            c = np.einsum("ij,ij->i", vec[:, a], vec[:, b])
            q -= 0.375 * (c + 1.0 / 3.0) ** 2
    # hydrogen bonds: donor i (via H1 or H2) to acceptor j
    close = r < 0.35
    nhb = np.zeros(n, int)
    for H in (H1, H2):
        oh = H - O
        oh -= box * np.round(oh / box)
        oh /= np.linalg.norm(oh, axis=1)[:, None]
        cosang = np.einsum("ijk,ik->ij", d, oh) / np.where(r > 0, r, 1.0)
        hb = close & (cosang > math.cos(math.radians(30.0)))
        nhb += hb.sum(axis=1)      # bonds donated by i
        nhb += hb.sum(axis=0)      # bonds accepted by j
    vol = float(np.prod(box))
    return hist, q, nhb, vol, n


def main():
    root, runlist, out = sys.argv[1], sys.argv[2], sys.argv[3]
    rbins = np.linspace(0.0, 0.9, 181)
    results = {}
    for rel in [l.strip() for l in open(runlist) if l.strip()]:
        d = os.path.join(root, rel)
        if not os.path.exists(os.path.join(d, "DONE")):
            continue
        cache = os.path.join(d, "analysis.json")
        post = os.path.join(d, "POST_DONE")
        if os.path.exists(cache) and os.path.exists(post) and os.path.getmtime(cache) > os.path.getmtime(post):
            try:
                with open(cache) as fh:
                    results[rel] = json.load(fh)    # analysed before: reuse
                continue
            except (ValueError, OSError):
                pass                                # damaged cache: recompute
        res = {}
        top = os.path.join(d, "topol.top")
        if os.path.exists(top):
            first = open(top).readline()
            if "lam =" in first:
                res["lam"] = float(first.split("lam =")[1].split()[0])
        dx = os.path.join(d, "density.xvg")
        if os.path.exists(dx):
            a = read_xvg(dx)
            rho = a[:, 1] / 1000.0                       # kg/m^3 -> g/cm^3
            sem, ladder = blocking_sem(rho)
            half = len(rho) // 2
            res.update(rho=float(rho.mean()), rho_sem=sem, rho_first_half=float(rho[:half].mean()),
                       rho_second_half=float(rho[half:].mean()), n_samples=int(len(rho)),
                       t_ps=[float(a[0, 0]), float(a[-1, 0])], blocking=ladder)
        if "rho" in res:
            # 100-ps block means of the density for convergence plots
            nb = int(round((a[-1, 0] - a[0, 0]) / 100.0))
            if nb >= 10:
                res["rho_blocks_100ps"] = [round(float(b.mean()), 6) for b in np.array_split(rho, nb)]
        et = os.path.join(d, "energy_terms.xvg")
        if os.path.exists(et):
            legends = []
            with open(et) as fh:
                for line in fh:
                    if line.startswith("@ s") and "legend" in line:
                        legends.append(line.split('"')[1])
            e = read_xvg_all(et)
            nmol = 1000.0
            for k, name in enumerate(legends):
                col = e[:, k + 1]
                key = name.replace(" ", "_").replace(".", "").replace("(", "").replace(")", "")
                res["E_" + key] = float(col.mean())
            if "E_Coulomb_SR" in res and "E_Coul_recip" in res:
                res["E_elec_per_mol_kJ"] = (res["E_Coulomb_SR"] + res["E_Coul_recip"]) / nmol
            if "E_Potential" in res:
                res["E_pot_per_mol_kJ"] = res["E_Potential"] / nmol
        # GROMACS 2025 gmx msd puts D in the xvg legend:
        #   @ s0 legend "D[   name OW] = 3.9136 (+/- 0.02) (1e-5 cm^2/s)"
        for ml in (os.path.join(d, "msd.xvg"), os.path.join(d, "msd.log")):
            if not os.path.exists(ml) or "D_1e-5_cm2_s" in res:
                continue
            for line in open(ml):
                if "D[" in line and "(+/-" in line:
                    try:
                        res["D_1e-5_cm2_s"] = float(line.split("] =")[1].split("(")[0].strip())
                        res["D_err"] = float(line.split("(+/-")[1].split(")")[0])
                    except (IndexError, ValueError):
                        pass
        fg = os.path.join(d, "frames.gro")
        if os.path.exists(fg):
            G = np.zeros(len(rbins) - 1)
            qs, hbs, n4, vols, nf = [], [], [], [], 0
            for names, xyz, box in read_gro_frames(fg):
                hist, q, nhb, vol, n = frame_structure(names, xyz, box, rbins)
                G += hist
                qs.append(q.mean()); hbs.append(nhb.mean()); n4.append(np.mean(nhb == 4)); vols.append(vol)
                nf += 1
            rmid = 0.5 * (rbins[1:] + rbins[:-1])
            shell = 4.0 / 3.0 * math.pi * (rbins[1:] ** 3 - rbins[:-1] ** 3)
            dens_pairs = n * (n - 1) / 2.0 / np.mean(vols)
            g = G / nf / (shell * dens_pairs)
            res.update(q_mean=float(np.mean(qs)), q_sem=float(np.std(qs, ddof=1) / math.sqrt(nf)),
                       hb_mean=float(np.mean(hbs)), hb_sem=float(np.std(hbs, ddof=1) / math.sqrt(nf)),
                       frac4_mean=float(np.mean(n4)), nframes=nf,
                       r_nm=[round(float(x), 4) for x in rmid], g_oo=[round(float(x), 4) for x in g])
        results[rel] = res
        if os.path.exists(post) and "hb_mean" in res and "D_1e-5_cm2_s" in res:
            with open(cache, "w") as fh:
                json.dump(res, fh)
        print(rel, {k: v for k, v in res.items() if k in ("rho", "rho_sem", "hb_mean", "q_mean", "frac4_mean", "D_1e-5_cm2_s")}, flush=True)
    with open(out, "w") as fh:
        json.dump(results, fh)


if __name__ == "__main__":
    main()
