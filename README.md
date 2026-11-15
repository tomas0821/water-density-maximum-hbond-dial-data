# Water's density maximum under a hydrogen-bond dial — simulation data

Molecular dynamics data, inputs and analysis scripts for a study of how the temperature of maximum density of
liquid water moves when the hydrogen bonds of a water model are weakened or strengthened. In the TIP4P/2005
model every partial charge is multiplied by a factor λ (the "hydrogen-bond dial"; λ = 1 is the standard
model), which scales the electrostatic energy by λ² while the geometry and the oxygen–oxygen Lennard-Jones
interaction stay fixed. SPC/E and TIP3P are included for comparison.

- 95 production simulations at 1 bar (λ = 0.85, 0.90, 0.925, 0.95, 0.975, 1.00, 1.05; SPC/E; TIP3P), 10 K apart,
  each 2 ns of equilibration and 10 ns of production with 1000 molecules, GROMACS 2025.1, 1 fs time step;
- 3 time-step checks (λ = 1 at 0.5 fs) and 2 runs that test the scaling relation between the dial and the
  Lennard-Jones depth (`mapping_check.json`);
- tables of densities, hydrogen-bond statistics, energies and diffusion coefficients, fitted temperatures of
  maximum density, and real-water reference densities (IAPWS);
- final configurations of four runs at corresponding temperatures (drawn in the article's snapshot figure);
- GROMACS inputs for every model and the scripts that produced every table.

## Files

| File | What it is |
|---|---|
| `density_vs_temperature.csv` | One row per simulation: model, dial setting λ, temperature (K), mean density and its standard error (g cm⁻³), the means of the first and second halves of the run (a drift check), hydrogen bonds per molecule, tetrahedral order *q*, fraction of molecules with four hydrogen bonds, electrostatic and total potential energy per molecule (kJ mol⁻¹), self-diffusion coefficient (10⁻⁵ cm² s⁻¹), and whether the run was used in the fits ("no" = too sluggish to settle: D < 0.1 × 10⁻⁵ cm² s⁻¹). The `t4p2005_l1.00_dt05` rows are the time-step check (0.5 fs). |
| `tmd_table.csv` | Fitted temperature of maximum density and maximum density for every model, with uncertainties (`TMD_p16_K`/`TMD_p84_K`: 16–84 % range of all Monte Carlo refits; `frac_fits_inside`: fraction of refits with an interior maximum), and the hydrogen bonds per molecule and four-bonded fraction at the maximum. `below`/`above` rows give the bound when the maximum lies outside the simulated range. |
| `per_run_analysis.json` | Full analysis of every run, keyed `"<model>/T<K>"`: mean density `rho` and its blocking error `rho_sem` (g cm⁻³), the blocking ladder `blocking`, 100 block means of 100 ps each `rho_blocks_100ps`, hydrogen bonds `hb_mean`, tetrahedral order `q_mean`, four-bonded fraction `frac4_mean`, energy terms `E_*` (kJ mol⁻¹, whole box), diffusion coefficient `D_1e-5_cm2_s`, and the O–O pair distribution function (`r_nm`, `g_oo`). |
| `mapping_check.json` | Same analysis for the two scaling-relation runs: standard TIP4P/2005 charges with the oxygen Lennard-Jones ε multiplied by 1/0.95² at 299.169 K and 254.848 K and 1.108033 bar, which should reproduce λ = 0.95 at 270 K and 230 K and 1 bar. |
| `snapshots/` | Final configurations (GROMACS `.gro`, 1000 molecules) of λ = 1.00, 0.95, 0.90 and 0.85 at 260, 230, 210 and 190 K, i.e. at corresponding temperatures T = λ² × 260 K; the article's snapshot figure shows the first and the last. `snapshot_stats.json` gives, for each pictured configuration, the hydrogen bonds per molecule, the four-bonded fraction, and the non-bonded oxygen pairs closer than 3.3 Å ("interstitial" pairs) per molecule. Drawn and counted by `inputs/scripts/render_snapshots.py` (needs the OVITO 3 Python module; expects the files in `data/snapshots/`). |
| `real_water_density_1atm.csv` | Density of real water at 1 atm, 238–330 K (IAPWS reference formulations). |
| `activity_find_the_maximum.py` | Classroom activity: plot the data, fit a parabola, find the maximum (runs in Google Colab). The tabulated maxima use a weighted cubic over ±35 K, so the parabola's maxima differ by up to a few kelvin. |
| `inputs/` | GROMACS topology and run files for every model (see `inputs/README.txt`) and the scripts used for the analysis (`inputs/scripts/`). |

**Model names.** `t4p2005_l1.00` is TIP4P/2005; `t4p2005_lX` is TIP4P/2005 with all charges multiplied by X
(the folders `l0.97` and `l0.93` are λ = 0.975 and 0.925); `spce` is SPC/E; `tip3p` is TIP3P;
`t4p2005_eps1.108_mapcheck` is the scaling-relation check.

**Running your own simulation.** Install GROMACS (free), copy one folder from `inputs/`, and run
`gmx solvate` → `gmx grompp`/`gmx mdrun` for `em.mdp`, `eq.mdp`, `prod.mdp` in turn (see `inputs/scripts/worker.sh`).
One temperature takes about 1.5 hours on a gaming laptop (NVIDIA RTX 3050, 189 ns/day at these settings),
so the six temperatures around a maximum run overnight.

## Citation

Archived on Zenodo: version 1.0.0, https://doi.org/10.5281/zenodo.24248781 (all versions:
https://doi.org/10.5281/zenodo.24248780). See `CITATION.cff`. The accompanying article is in preparation.

## License

Data and documentation: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (`LICENSE-data.md`).
Code (`*.py`, `*.sh`): MIT (`LICENSE`).

## Acknowledgment

This research was partially supported by a machine allocation on Kabré supercomputer at the Costa Rica National
High Technology Center.
