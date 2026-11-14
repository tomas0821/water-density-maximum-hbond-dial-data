"""Render slabs of the simulated liquid with OVITO (Fig. 4(c,d) of the paper; run from the project root).

For each configuration (data/snapshots/<run>.gro, the final frame of that production run): a 6 A slab seen
along z, 20 A across, molecules as balls and sticks, hydrogen bonds (Luzar-Chandler, as in the paper) as thin blue lines from
H to the acceptor O, and "interstitial" pairs -- two oxygens closer than R_INT that are NOT hydrogen bonded -- as
orange sticks. Writes figures/ovito/snap_<run>.png and data/snapshots/snapshot_stats.json (whole-box statistics).

Usage (from the project root, OVITO 3 Python module): python figures/ovito/render_snapshots.py [--preview]
(headless: bash figures/ovito/render.sh figures/ovito/render_snapshots.py)
"""
import json, math, os, sys
import numpy as np
from ovito.data import DataCollection, ParticleType
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import TachyonRenderer, Viewport

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from snapshot_analysis import read_gro, R_HB, ANG, R_INT   # same hydrogen-bond rule as the paper

SLAB = 6.0              # A
WINDOW = 20.0           # A: side of the square region shown (centered in the box)
RUNS = ["t4p2005_l1.00_T260", "t4p2005_l0.85_T190"]          # corresponding temperatures: T/lambda^2 = 260, 263 K
COL = dict(O=(0.80, 0.16, 0.13), H=(0.93, 0.93, 0.93), cov=(0.55, 0.55, 0.57),
           hb=(0.19, 0.51, 0.74), inter=(0.95, 0.50, 0.08))
PREVIEW = "--preview" in sys.argv


def build(run):
    names, xyz, box = read_gro(f"data/snapshots/{run}.gro")
    O, H1, H2 = xyz[names == "OW"], xyz[names == "HW1"], xyz[names == "HW2"]
    n = len(O)
    d = O[None, :, :] - O[:, None, :]; d -= box * np.round(d / box)
    r = np.sqrt((d ** 2).sum(-1)); np.fill_diagonal(r, np.inf)
    Hs, hb, hb_pairs = [], np.zeros((n, n), bool), []
    for k, H in enumerate((H1, H2)):
        oh = H - O; oh -= box * np.round(oh / box); Hs.append(O + oh)          # whole molecules
        u = oh / np.linalg.norm(oh, axis=1)[:, None]
        cos = np.einsum("ijk,ik->ij", d, u) / np.where(np.isfinite(r), r, 1.0)
        m = (r < R_HB) & (cos > math.cos(math.radians(ANG)))
        hb |= m
        hb_pairs += [(i, k, j) for i, j in zip(*np.nonzero(m))]
    bonded = hb | hb.T
    inter = (r < R_INT) & ~bonded
    nhb = hb.sum(1) + hb.sum(0)
    stats = dict(hb_mean=float(nhb.mean()), frac4=float((nhb == 4).mean()),
                 interstitial_pairs_per_molecule=float(inter.sum() / 2 / n),
                 frac_with_interstitial=float(inter.any(1).mean()), R_INT_nm=R_INT)
    # slab: molecules whose oxygen lies within SLAB/2 of the box middle along z
    L = box * 10.0
    Oa, H1a, H2a = O * 10.0, Hs[0] * 10.0, Hs[1] * 10.0
    zc = L[2] / 2
    dz = Oa[:, 2] - zc; dz -= L[2] * np.round(dz / L[2])
    keep = np.nonzero(np.abs(dz) < SLAB / 2)[0]
    idx = {int(i): q for q, i in enumerate(keep)}
    pos = np.vstack([Oa[keep], H1a[keep], H2a[keep]])
    shift = np.array([[0.0, 0.0, zc - Oa[i, 2] + dz[i]] for i in keep])        # bring the slab together in z
    pos += np.vstack([shift, shift, shift])
    m = len(keep)
    ptype = np.r_[np.ones(m, int), 2 * np.ones(2 * m, int)]
    top, pimg, col, wid = [], [], [], []

    def add(a, b, c, w):
        v = pos[b] - pos[a]
        s = -np.round(v[:2] / L[:2])
        top.append((a, b)); pimg.append((int(s[0]), int(s[1]), 0)); col.append(c); wid.append(w)

    for q in range(m):                                                      # covalent O-H
        add(q, m + q, COL["cov"], 0.30); add(q, 2 * m + q, COL["cov"], 0.30)
    for i, k, j in hb_pairs:                                                # H ... acceptor O
        if i in idx and j in idx:
            add(m * (k + 1) + idx[i], idx[j], COL["hb"], 0.16)
    for i, j in zip(*np.nonzero(np.triu(inter))):                           # interstitial O ... O
        if i in idx and j in idx:
            add(idx[i], idx[j], COL["inter"], 0.40)
    dc = DataCollection()
    cell = dc.create_cell(np.array([[L[0], 0, 0, 0], [0, L[1], 0, 0], [0, 0, L[2], 0]]), pbc=(True, True, False))
    cell.vis.enabled = False
    parts = dc.create_particles(count=len(pos))
    parts.create_property("Position", data=pos)
    tp = parts.create_property("Particle Type", data=ptype)
    tp.types.append(ParticleType(id=1, name="O", color=COL["O"], radius=0.62))
    tp.types.append(ParticleType(id=2, name="H", color=COL["H"], radius=0.34))
    bonds = parts.create_bonds(count=len(top))
    bonds.create_property("Topology", data=np.array(top))
    bonds.create_property("Periodic Image", data=np.array(pimg))
    bonds.create_property("Color", data=np.array(col))
    bonds.create_property("Width", data=np.array(wid))
    return dc, L, stats


def render(run):
    dc, L, stats = build(run)
    pipe = Pipeline(source=StaticSource(data=dc)); pipe.add_to_scene()
    c = L / 2
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0, 0, -1), camera_pos=(c[0], c[1], 80.0), fov=WINDOW / 2)
    vp.camera_up = (0, 1, 0)
    rend = (TachyonRenderer(ambient_occlusion=False, antialiasing_samples=4, shadows=False) if PREVIEW else
            TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.85, ambient_occlusion_samples=16,
                            direct_light_intensity=0.9, shadows=False, antialiasing_samples=16))
    out = f"figures/ovito/snap_{run}.png"
    vp.render_image(size=(1000, 1000), filename=out, background=(1, 1, 1), renderer=rend)
    pipe.remove_from_scene()
    print(out, {k: round(v, 3) for k, v in stats.items()})
    return stats


if __name__ == "__main__":
    allstats = {run: render(run) for run in RUNS}
    json.dump(allstats, open("data/snapshots/snapshot_stats.json", "w"), indent=1)
