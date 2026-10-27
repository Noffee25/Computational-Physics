for s in $(seq 1001 1600); do ./hsim logis $s 100 10000 1 0; done > logis_seeds_mine.out
for s in $(seq 1001 1600); do ./hsim sm64 $s 100 10000 1 0; done > sm64_seeds_mine.out
