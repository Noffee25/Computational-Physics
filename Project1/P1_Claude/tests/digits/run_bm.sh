#!/bin/bash
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits
L=50000
{
for j in 1 32; do echo "./bm 5 777 $j $L out/bm_g5_j$j.txt"; done
for j in 1 32 53; do echo "./bm 7 777 $j $L out/bm_g7_j$j.txt"; done
for j in 1 53; do echo "./bm 6 777 $j $L out/bm_g6_j$j.txt"; done
for j in 1 31 53; do echo "./bm 4 777 $j $L out/bm_g4_j$j.txt"; done
for j in 1 22; do echo "./bm 3 777 $j $L out/bm_g3_j$j.txt"; done
for j in 1 14 16 18 20 24 28; do echo "./bm 1 777 $j $L out/bm_g1_j$j.txt"; done
for j in 1 16 20 24; do echo "./bm 2 777 $j $L out/bm_g2_j$j.txt"; done
for j in 1 40 53; do echo "./bm 9 777 $j $L out/bm_g9_j$j.txt"; done
echo "./bm 8 12345678 1 $L out/bm_g8_j1.txt"
} > joblist_bm.txt
cat joblist_bm.txt | xargs -P 2 -I{} bash -c "{}"
