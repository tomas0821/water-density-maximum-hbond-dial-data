#!/usr/bin/env python3
"""Generate GROMACS inputs for the water 'hydrogen-bond dial' study.

Models
------
* TIP4P/2005 (Abascal & Vega, JCP 123, 234505, 2005) with both partial charges
  multiplied by a factor lam ("hydrogen-bond dial"); lam = 1 is the original model.
  Geometry and the oxygen Lennard-Jones site are left unchanged, as in
  Horstmann & Vogel, JCP 154, 054502 (2021).
* SPC/E (Berendsen, Grigera & Straatsma, JPC 91, 6269, 1987).
* TIP3P (Jorgensen et al., JCP 79, 926, 1983), original form (no LJ on H).

Every state point gets its own directory with topol.top, em.mdp, eq.mdp and
prod.mdp. Run lists for the worker job are written to runlists/.

Usage: python3 make_inputs.py [--outdir runs] [--nmol 1000] [--dt 0.001] [--eq-ns 2] [--prod-ns 10]
(defaults = the settings used in the paper)
"""
import argparse
import os

KB = 0.0083144626  # kJ mol^-1 K^-1

# ---------------------------------------------------------------- topologies
TIP4P2005_TOP = """; TIP4P/2005 with partial charges scaled by lam = {lam:.3f}
; q_H = {qh:+.6f} e, q_M = {qm:+.6f} e (original: +0.5564, -1.1128)
[ defaults ]
; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
  1        2          no         1.0      1.0

[ atomtypes ]
; name     at.num  mass      charge  ptype  sigma(nm)  epsilon(kJ/mol)
  OW_t4p   8       15.9994   0.0     A      0.315890   0.774908
  HW_t4p   1        1.0080   0.0     A      0.000000   0.000000
  MW       0        0.0000   0.0     D      0.000000   0.000000

[ moleculetype ]
; molname  nrexcl
  SOL      2

[ atoms ]
; nr  type    resnr  residue  atom  cgnr  charge       mass
  1   OW_t4p  1      SOL      OW    1      0.000000    15.9994
  2   HW_t4p  1      SOL      HW1   1     {qh:+.6f}     1.0080
  3   HW_t4p  1      SOL      HW2   1     {qh:+.6f}     1.0080
  4   MW      1      SOL      MW    1     {qm:+.6f}     0.0000

[ settles ]
; OW  funct  doh      dhh
  1   1      0.09572  0.15139007

[ virtual_sites3 ]
; site  from      funct  a           b
  4     1  2  3   1      0.13193777  0.13193777

[ exclusions ]
  1 2 3 4
  2 1 3 4
  3 1 2 4
  4 1 2 3

[ system ]
TIP4P/2005 lam={lam:.3f}

[ molecules ]
SOL {nmol}
"""

THREE_SITE_TOP = """; {name}
[ defaults ]
; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
  1        2          no         1.0      1.0

[ atomtypes ]
; name   at.num  mass      charge  ptype  sigma(nm)  epsilon(kJ/mol)
  OW_x   8       15.9994   0.0     A      {sig:.6f}   {eps:.6f}
  HW_x   1        1.0080   0.0     A      0.000000   0.000000

[ moleculetype ]
; molname  nrexcl
  SOL      2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge      mass
  1   OW_x  1      SOL      OW    1     {qo:+.6f}   15.9994
  2   HW_x  1      SOL      HW1   1     {qh:+.6f}    1.0080
  3   HW_x  1      SOL      HW2   1     {qh:+.6f}    1.0080

[ settles ]
; OW  funct  doh     dhh
  1   1      {doh:.5f}  {dhh:.8f}

[ exclusions ]
  1 2 3
  2 1 3
  3 1 2

[ system ]
{name}

[ molecules ]
SOL {nmol}
"""

THREE_SITE = {
    # name: sigma(nm), epsilon(kJ/mol), qO, qH, dOH(nm), dHH(nm)
    "spce": dict(name="SPC/E", sig=0.316557, eps=0.650194, qo=-0.8476, qh=0.4238,
                 doh=0.1, dhh=0.16329809),
    "tip3p": dict(name="TIP3P (original, no LJ on H)", sig=0.315061, eps=0.636386,
                  qo=-0.834, qh=0.417, doh=0.09572, dhh=0.15139007),
}

# ---------------------------------------------------------------- mdp files
MDP_COMMON = """cutoff-scheme     = Verlet
nstcalcenergy     = 100
nstenergy         = 100       ; 0.1 ps at dt = 1 fs
nstlog            = 100000
coulombtype       = PME
rcoulomb          = 0.9
fourierspacing    = 0.10
pme-order         = 4
ewald-rtol        = 1e-5
vdwtype           = Cut-off
vdw-modifier      = Potential-shift
rvdw              = 0.9
DispCorr          = EnerPres
tcoupl            = V-rescale
tc-grps           = System
tau-t             = 1.0
ref-t             = {T}
pcoupl            = C-rescale
pcoupltype        = isotropic
tau-p             = 2.0
ref-p             = 1.0
compressibility   = 4.5e-5
"""

EM_MDP = """integrator        = steep
emtol             = 500
nsteps            = 5000
cutoff-scheme     = Verlet
coulombtype       = PME
rcoulomb          = 0.9
vdwtype           = Cut-off
rvdw              = 0.9
DispCorr          = EnerPres
"""

EQ_MDP = """integrator        = md
dt                = {dt}
nsteps            = {nsteps}
gen-vel           = yes
gen-temp          = {T}
gen-seed          = -1
continuation      = no
nstxout-compressed = 0
""" + MDP_COMMON

PROD_MDP = """integrator        = md
dt                = {dt}
nsteps            = {nsteps}
gen-vel           = no
continuation      = yes
nstxout-compressed = {nxtc}    ; every 10 ps
compressed-x-precision = 1000
""" + MDP_COMMON

# ---------------------------------------------------------------- state points
LAMBDAS = {
    1.05: range(250, 331, 10),
    1.00: range(240, 321, 10),
    0.95: range(220, 301, 10),
    0.90: range(200, 291, 10),
    0.85: range(180, 271, 10),
}
THREE_SITE_T = {
    "spce": range(210, 301, 10),
    "tip3p": range(160, 251, 10),
}


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="runs")
    ap.add_argument("--nmol", type=int, default=1000)
    ap.add_argument("--dt", type=float, default=0.001)
    ap.add_argument("--eq-ns", type=float, default=2.0)
    ap.add_argument("--prod-ns", type=float, default=10.0)
    ap.add_argument("--points", default="",
                    help="only these TIP4P/2005-lambda points, e.g. '0.975:230-310,0.90:170-190' (10 K steps)")
    args = ap.parse_args()
    if args.points:
        LAMBDAS.clear(); THREE_SITE_T.clear(); THREE_SITE.clear()
        for item in args.points.split(","):
            lam, rng = item.split(":")
            lo, hi = (int(x) for x in rng.split("-"))
            LAMBDAS[float(lam)] = range(lo, hi + 1, 10)

    nsteps_eq = int(round(args.eq_ns * 1000 / args.dt))
    nsteps_prod = int(round(args.prod_ns * 1000 / args.dt))
    nxtc = int(round(10.0 / args.dt))  # 10 ps

    entries = []  # (model_key, init_box, T, dir)
    for lam, temps in LAMBDAS.items():
        key = f"t4p2005_l{lam:.2f}"
        top = TIP4P2005_TOP.format(lam=lam, qh=0.5564 * lam, qm=-1.1128 * lam, nmol=args.nmol)
        for T in temps:
            entries.append((key, "tip4p", T, top))
    for key, p in THREE_SITE.items():
        top = THREE_SITE_TOP.format(nmol=args.nmol, **p)
        for T in THREE_SITE_T[key]:
            entries.append((key, "spc", T, top))

    os.makedirs(os.path.join(args.outdir, "runlists"), exist_ok=True)
    allruns = []
    for key, init, T, top in entries:
        d = os.path.join(args.outdir, key, f"T{T}")
        os.makedirs(d, exist_ok=True)
        write(os.path.join(d, "topol.top"), top)
        write(os.path.join(d, "em.mdp"), EM_MDP)
        write(os.path.join(d, "eq.mdp"), EQ_MDP.format(dt=args.dt, nsteps=nsteps_eq, T=T))
        write(os.path.join(d, "prod.mdp"), PROD_MDP.format(dt=args.dt, nsteps=nsteps_prod, nxtc=nxtc, T=T))
        write(os.path.join(d, "init_type"), init + "\n")
        allruns.append(os.path.join(key, f"T{T}"))
    listname = "extra.txt" if args.points else "all.txt"
    write(os.path.join(args.outdir, "runlists", listname), "\n".join(allruns) + "\n")
    print(f"{len(allruns)} state points written under {args.outdir}/ "
          f"(dt={args.dt} ps, eq {args.eq_ns} ns, prod {args.prod_ns} ns, N={args.nmol})")


if __name__ == "__main__":
    main()
