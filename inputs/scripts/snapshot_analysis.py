"""Shared definitions for the snapshot figures: read a GROMACS .gro file, and the distance/angle cut-offs.

Hydrogen bonds use the same rule as the paper (Luzar-Chandler: R_OO < 0.35 nm, angle O-H...O_acceptor < 30 deg).
An "interstitial" pair is two oxygens closer than R_INT that are NOT hydrogen bonded to each other: the molecules
that fill the gap between the first and second shells of the network (Sec. V C of the paper).
"""
import numpy as np

R_HB, ANG = 0.35, 30.0      # nm, deg
R_INT = 0.33                # nm: non-bonded neighbor inside the first minimum of g_OO (3.3 A)


def read_gro(path):
    """Atom names, coordinates (nm) and box lengths (nm) of a rectangular-box .gro file."""
    lines = open(path).read().splitlines()
    n = int(lines[1]); names = []; xyz = np.empty((n, 3))
    for i in range(n):
        l = lines[2 + i]; names.append(l[10:15].strip()); xyz[i] = (float(l[20:28]), float(l[28:36]), float(l[36:44]))
    box = np.array([float(v) for v in lines[2 + n].split()[:3]])
    return np.array(names), xyz, box
