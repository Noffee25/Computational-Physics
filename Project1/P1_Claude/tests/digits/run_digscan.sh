#!/bin/bash
# digit scans: 9 generators x {double, single}, N = 1e8, seed 777 (MIDSQ 12345678)
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits
N=100000000
jobs=()
for prec in 0 1; do
  for g in 1 2 3 4 5 6 7 8 9; do
    seed=777; [ $g -eq 8 ] && seed=12345678
    echo "./digscan $g $prec $N $seed out/dig_g${g}_p${prec}.txt"
  done
done > joblist.txt
echo "./digscan 9 0 10000000 777 out/dig_g9_p0_N1e7.txt" >> joblist.txt
echo "./digscan 9 1 10000000 777 out/dig_g9_p1_N1e7.txt" >> joblist.txt
cat joblist.txt | xargs -P 3 -I{} bash -c "{}"
