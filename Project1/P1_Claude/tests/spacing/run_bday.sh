#!/bin/bash
cd /Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude/tests/spacing
for spec in "RANDU 1 777 0" "RANDM 2 777 0" "RAND 3 777 0" "MINSTD 4 777 0" "R250 5 777 0" "RANNUM 6 777 0" "MT 7 777 0" "LOGISTIC 9 777 0" "MIDSQ 8 12345678 0" "MTsingle 7 777 1"; do
  set -- $spec
  ./birthday $2 $3 $4 plan_bday.txt bday/$1.txt
done
