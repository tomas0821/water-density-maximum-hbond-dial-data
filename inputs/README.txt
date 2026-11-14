One set of GROMACS input files per model, written by scripts/make_inputs.py.
The files shown are for 260 K (TIP3P: 250 K); for another temperature change ref-t and gen-temp
in eq.mdp and ref-t in prod.mdp. init_type says which water box gmx solvate starts from
(tip4p.gro or spc216.gro, both shipped with GROMACS).

Order of commands for one temperature (see scripts/worker.sh for the exact options):
  gmx solvate -cs tip4p.gro -box 3.2 3.2 3.2 -maxsol 1000 -o conf.gro      (spc216.gro for SPC/E, TIP3P)
  gmx grompp -f em.mdp   -c conf.gro -p topol.top -o em.tpr   && gmx mdrun -deffnm em
  gmx grompp -f eq.mdp   -c em.gro   -p topol.top -o eq.tpr   && gmx mdrun -deffnm eq
  gmx grompp -f prod.mdp -c eq.gro -t eq.cpt -p topol.top -o prod.tpr && gmx mdrun -deffnm prod
  echo Density | gmx energy -f prod.edr -o density.xvg

t4p2005_l1.00_dt05: the time-step check (0.5 fs instead of 1 fs; same topology as t4p2005_l1.00).
t4p2005_eps1.108_mapcheck: Problem 6 (test of the scaling argument): standard TIP4P/2005 charges with the
oxygen Lennard-Jones epsilon multiplied by 1/0.95^2 = 1.108, at 299.169 K and 1.108033 bar. It should reproduce
the density of lambda = 0.95 at 270 K and 1 bar (0.99016 +- 0.00021 g/cm3); our run gave 0.99014 +- 0.00021.
For the second pair use ref-t = gen-temp = 254.848 K (compare lambda = 0.95 at 230 K).
Folder names l0.97 and l0.93 are the dial settings 0.975 and 0.925.
scripts/: input generator, cluster worker, post-processing and analysis scripts used for the paper.
