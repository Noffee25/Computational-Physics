#!/bin/bash
# stream 2: cell tests, Knuth battery, triples, second level, deep bins
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/battery
seedof() { if [ $1 -eq 8 ]; then echo 12345678; else echo 777; fi; }
for g in 1 2 3 4 5 6 7 8 9; do
  s=$(seedof $g)
  ./cells $g $s 2 0 10000000 9 16 32 64 128 256 512 1000 2048 4096 out/cells/c_${g}_d2.txt
  ./cells $g $s 3 0 10000000 8 4 8 16 32 64 100 128 256 out/cells/c_${g}_d3.txt
  ./cells $g $s 4 0 10000000 5 4 8 16 32 64 out/cells/c_${g}_d4.txt
  ./cells $g $s 3 1 10000000 3 32 64 100 out/cells/c_${g}_lag.txt
done
echo CELLS_DONE
for g in 1 2 3 4 5 6 7 8 9; do
  s=$(seedof $g)
  ./knuth $g $s 0 24 out/kn/kn_${g}_v0.txt
  ./knuth $g $s 1 24 out/kn/kn_${g}_v1.txt
done
echo KNUTH_DONE
for g in 1 2 4 7 6 9; do ./triples $g 777 29 1 2 out/tri/tri_${g}_1_2.txt; done
./triples 5 777 29 103 250 out/tri/tri_5_103_250.txt
./triples 7 777 29 103 250 out/tri/tri_7_103_250.txt
./triples 1 777 29 2 4 out/tri/tri_1_2_4.txt
./triples 7 777 29 2 4 out/tri/tri_7_2_4.txt
echo TRIPLES_DONE
for g in 1 2 3 4 5 6 7 8 9; do s=$(seedof $g); ./secondlevel $g $s 29 20 23 out/sl/sl_$g.txt; done
echo SL_DONE
./deepbins 1 777 1000000000 10000000 out/deep/deep_1.txt
./deepbins 7 777 1000000000 10000000 out/deep/deep_7.txt
./deepbins 2 777 1000000000 10000000 out/deep/deep_2.txt
echo DEEP_DONE
