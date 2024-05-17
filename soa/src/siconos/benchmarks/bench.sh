#!/bin/sh
cd `dirname $0`
time=`which time`
for n in `seq $1 $2 $3`; do
  echo -n "$n "; $time -f 'RESULTS : -- %E %e %M %R --' python3 bench-disks.py bullet $n 2>&1 |\
    tee bench-disks-$n-bullet.log |\
    sed -n "s/RESULTS : -- \(.*\) --/\1/gp" |\
    tr '\n' ' '
  grep -q "step 2000 of 2000" bench-disks-$n-bullet.log ||\
    { echo "bad bench-disks-$n-bullet.log" ; exit 1; }

  $time -f 'RESULTS : -- %E %e %M %R --' python3 bench-disks.py vnative $n 2>&1 |\
    tee bench-disks-$n-vkernel.log | sed -n "s/RESULTS : -- \(.*\) --/\1/gp" |\
    tr '\n' ' '
  grep -q "step 2000 of 2000" bench-disks-$n-vkernel.log ||\
    { echo "bad bench-disks-$n-vkernel.log" ; exit 1; }
  echo
done
