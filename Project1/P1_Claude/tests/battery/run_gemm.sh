#!/bin/bash
# stream 1: AMX/GEMM-heavy jobs (lag correlations, three-point maps)
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery
for g in 2 1 7 6 4 3 5 9; do ./lagcorr $g 777 29 300 out/lag/lag_$g.txt; done
./lagcorr 8 12345678 24 300 out/lag/lag_8.txt
# extra seeds for convergence plot (lags 1-3 only matter)
for s in 1 314159 65539; do ./lagcorr 2 $s 29 8 out/lag/lagconv_2_$s.txt; ./lagcorr 1 $s 29 8 out/lag/lagconv_1_$s.txt; done
for s in 1 314159; do ./lagcorr 7 $s 29 8 out/lag/lagconv_7_$s.txt; done
for g in 5 1 7 2 9 4 6 3; do ./c3map $g 777 24 256 out/c3/c3_$g.bin; done
echo GEMM_DONE
