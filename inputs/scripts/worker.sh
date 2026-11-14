#!/bin/bash
# Worker job for the water hydrogen-bond-dial study (GPU cluster nodes).
# Runs the state points listed in RUNLIST, NSLOT at a time on one GPU, each with
# NTOMP OpenMP threads pinned to its own cores. Each state point: solvate -> EM ->
# NPT equilibration -> NPT production. Restartable: finished points are skipped,
# an interrupted production run resumes from its checkpoint.
#
# Usage (inside sbatch): worker.sh RUNLIST NSLOT NTOMP
#SBATCH --job-name=wtmd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x-%j.out

set -u
RUNLIST=${1:?runlist}
NSLOT=${2:-5}
NTOMP=${3:-4}
ROOT=${ROOT:-/path/to/water-tmd/runs}
GMXROOT=/opt/tools/gromacs-2025.1
export PATH=$GMXROOT/bin:$PATH
export LD_LIBRARY_PATH=/opt/tools/cuda-12.4.0/lib64:$GMXROOT/lib64:${LD_LIBRARY_PATH:-}
export GMX_MAXBACKUP=-1
TOPDIR=$GMXROOT/share/gromacs/top

echo "host $(hostname)  start $(date)  runlist $RUNLIST  slots $NSLOT x $NTOMP threads"
nvidia-smi -L 2>/dev/null

run_point() {   # $1 = relative dir, $2 = pin offset
  local rel=$1 off=$2 d=$ROOT/$1
  cd "$d" || return 1
  [ -f DONE ] && { echo "[$rel] already done"; return 0; }
  local mdflags="-ntmpi 1 -ntomp $NTOMP -pin on -pinoffset $off -pinstride 1 -nb gpu -pme gpu -bonded cpu -update auto"
  if [ ! -f eq.gro ]; then
    if [ ! -f em.gro ]; then
      local init; init=$(cat init_type)
      local src=$TOPDIR/tip4p.gro; [ "$init" = spc ] && src=$TOPDIR/spc216.gro
      gmx solvate -cs "$src" -box 3.2 3.2 3.2 -maxsol 1000 -o conf.gro > solvate.log 2>&1 || { echo FAILED-solvate > FAILED; return 1; }
      gmx grompp -f em.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 1 > grompp_em.log 2>&1 || { echo FAILED-grompp-em > FAILED; return 1; }
      # steepest descent: PME cannot run on the GPU, so only the short-range part goes there
      gmx mdrun -deffnm em -ntmpi 1 -ntomp $NTOMP -pin on -pinoffset $off -pinstride 1 -nb gpu -pme cpu > mdrun_em.log 2>&1 || { echo FAILED-em > FAILED; return 1; }
    fi
    gmx grompp -f eq.mdp -c em.gro -p topol.top -o eq.tpr -maxwarn 1 > grompp_eq.log 2>&1 || { echo FAILED-grompp-eq > FAILED; return 1; }
    gmx mdrun -deffnm eq $mdflags > mdrun_eq.log 2>&1 || { echo FAILED-eq > FAILED; return 1; }
  fi
  if [ ! -f prod.tpr ]; then
    gmx grompp -f prod.mdp -c eq.gro -t eq.cpt -p topol.top -o prod.tpr -maxwarn 1 > grompp_prod.log 2>&1 || { echo FAILED-grompp-prod > FAILED; return 1; }
  fi
  gmx mdrun -deffnm prod -cpi prod.cpt $mdflags > mdrun_prod.log 2>&1 || { echo FAILED-prod > FAILED; return 1; }
  [ -f prod.gro ] || { echo FAILED-noprodgro > FAILED; return 1; }
  printf "Density\n" | gmx energy -f prod.edr -o density.xvg > energy.log 2>&1
  rm -f FAILED; date > DONE
  echo "[$rel] done $(date)"
}

mapfile -t POINTS < <(grep -v '^\s*$' "$RUNLIST")
echo "${#POINTS[@]} state points (shared queue: a slot claims a point with an atomic mkdir)"
for ((s = 0; s < NSLOT; s++)); do
  (
    for p in "${POINTS[@]}"; do
      [ -f "$ROOT/$p/DONE" ] && continue
      mkdir "$ROOT/$p/CLAIM" 2>/dev/null || continue
      echo "$(hostname) slot $s job ${SLURM_JOB_ID:-?} $(date)" > "$ROOT/$p/CLAIM/who"
      run_point "$p" $((s * NTOMP))
    done
  ) &
done
wait
echo "all slots finished $(date)"
