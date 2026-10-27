#!/bin/bash
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/digits
N=50000000
for prec in 0 1; do
  for g in 1 2 3 4 5 6 7 8 9; do
    seed=777; [ $g -eq 8 ] && seed=12345678
    echo "./bitscan $g $prec $N $seed out/bit_g${g}_p${prec}.txt out/dump_g${g}_p${prec}.txt 65536"
  done
done > joblist_bits.txt
cat joblist_bits.txt | xargs -P 2 -I{} bash -c "{}"
