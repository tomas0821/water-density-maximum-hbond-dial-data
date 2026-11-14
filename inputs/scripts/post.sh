#!/bin/bash
# Post-process finished state points on the cluster (light, niced).
# Usage: bash post.sh            (all DONE points without POST_DONE)
#        ONLY=t4p2005_l1.00/T280 bash post.sh   (one point, for testing)
ROOT=/path/to/water-tmd/runs
GMXROOT=/opt/tools/gromacs-2025.1
export PATH=$GMXROOT/bin:$PATH
export LD_LIBRARY_PATH=/opt/tools/cuda-12.4.0/lib64:$GMXROOT/lib64:${LD_LIBRARY_PATH:-}
export GMX_MAXBACKUP=-1
cd $ROOT || exit 1
for rel in $(ls -d */T*/DONE 2>/dev/null | sed 's|/DONE$||'); do
  [ -n "$ONLY" ] && [ "$rel" != "$ONLY" ] && continue
  d=$ROOT/$rel
  [ -f $d/POST_DONE ] && continue
  cd $d || continue
  printf "Potential\nCoulomb-(SR)\nCoul.-recip.\nLJ-(SR)\nDisper.-corr.\nVolume\nDensity\n\n" | nice -n 19 gmx energy -f prod.edr -o energy_terms.xvg > energy_terms.log 2>&1
  echo System | nice -n 19 gmx trjconv -f prod.xtc -s prod.tpr -o frames.gro -dt 100 > trjconv.log 2>&1
  nice -n 19 gmx msd -f prod.xtc -s prod.tpr -sel "name OW" -o msd.xvg -trestart 10 -maxtau 3000 -beginfit 300 -endfit 2500 > msd.log 2>&1
  if [ -s energy_terms.xvg ] && [ -s frames.gro ] && grep -q "legend \"D\[" msd.xvg; then
    touch POST_DONE; echo "post-processed $rel"
  else
    echo "POST FAILED $rel"
  fi
  cd $ROOT
done
